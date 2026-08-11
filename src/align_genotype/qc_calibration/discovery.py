"""Metamist cohort discovery - which datasets and MultiQC reports go into a manifest.

The query function is injected as an argument (``QueryFn``), so the selection logic -
eligibility filtering, latest-analysis-wins - is fully unit-testable with a fake,
leaving only a ~4-line GraphQL adapter (``default_query``) untested. Discovery is
deliberately a separate step from parsing, not something later commands do implicitly:
the resulting manifest is meant to be reviewed and edited by the operator - dropping a
bad cohort, pinning an older analysis, substituting a local file - before anything
expensive is parsed.
"""

import logging
from collections.abc import Callable
from datetime import datetime as dt
from datetime import timezone
from typing import Any

from align_genotype.qc_calibration import tomlio
from align_genotype.qc_calibration.manifest import Cohort, Manifest

QueryFn = Callable[[str, dict[str, Any] | None], dict[str, Any]]

DATASETS_QUERY = """
    query Datasets {
        myProjects {
            name
            dataset
            meta
        }
    }
"""

ANALYSES_QUERY = """
    query DatasetData($datasetName: String!) {
        project(name: $datasetName) {
            analyses(status: {eq: COMPLETED}, type: {eq: "qc"}) {
                id
                meta
                output
                timestampCompleted
            }
        }
    }
"""

# Staging/derivative projects that must never be treated as source cohorts, even though
# they carry `is_seqr: true` like real source datasets do. `is_seqr` marks a dataset as
# part of the seqr pipeline family; these name tokens instead flag seqr-loader *staging*
# projects (test/training copies, or the "-seqr" loader project itself) that shadow a
# real dataset. Excluding by name on top of requiring `is_seqr` is not a contradiction -
# it is filtering two different things - so do not "simplify" this to just `is_seqr`.
EXCLUDED_NAME_TOKENS = ('test', 'training', 'seqr')


class DiscoveryError(RuntimeError):
    """Discovery found no cohorts to put in the manifest."""


def default_query(query_text: str, variables: dict[str, Any] | None = None) -> dict[str, Any]:
    """Run a GraphQL query against Metamist.

    `metamist` is imported lazily so that nothing else in this package - and no test -
    needs it installed; only this ~4-line adapter ever touches it.
    """
    from metamist.graphql import gql, query  # noqa: PLC0415 - lazy import, see docstring

    return query(gql(query_text), variables=variables or {})


def is_eligible(project: dict[str, Any]) -> bool:
    """Whether `project` is a real source dataset, not a staging/derivative project."""
    meta = project.get('meta') or {}
    if not meta.get('is_seqr', False):
        return False
    name = project.get('name', '')
    return not any(token in name for token in EXCLUDED_NAME_TOKENS)


def latest_analysis(
    analyses: list[dict[str, Any]],
    seq_type: str,
    dataset_label: str | None = None,
) -> dict[str, Any] | None:
    """The most recently completed analysis of `seq_type` that produced an output.

    Selection is by `timestampCompleted`, not by id or list order: a higher id with an
    older timestamp must lose. This is talking to a live external service, so malformed
    rows are expected: an analysis with a missing or non-string `timestampCompleted`
    can't be ordered against the others, so it is excluded rather than compared - it's
    better to silently drop a candidate we can't rank than to silently pick the wrong
    one, or crash with a `TypeError`/`KeyError` instead of a clean `DiscoveryError`.
    `dataset_label` is only used to name the dataset in that warning.
    """
    candidates = []
    for analysis in analyses:
        if (analysis.get('meta') or {}).get('sequencing_type') != seq_type or not analysis.get('output'):
            continue
        timestamp = analysis.get('timestampCompleted')
        if not isinstance(timestamp, str):
            logging.warning(
                f'discovery: dataset {dataset_label!r} analysis {analysis.get("id")!r} has no usable '
                f'timestampCompleted ({timestamp!r}); excluding it from latest-analysis selection',
            )
            continue
        candidates.append(analysis)
    if not candidates:
        return None
    return max(candidates, key=lambda a: a['timestampCompleted'])


def build_manifest(seq_type: str, query_fn: QueryFn = default_query, generated: str | None = None) -> Manifest:
    """Discover eligible cohorts for `seq_type` and assemble a `Manifest`."""
    if generated is None:
        generated = dt.now(tz=timezone.utc).isoformat(timespec='seconds')

    projects = query_fn(DATASETS_QUERY, None).get('myProjects', [])
    eligible = [p for p in projects if is_eligible(p)]

    cohorts: list[Cohort] = []
    for project in eligible:
        label = project.get('dataset') or project.get('name')
        try:
            tomlio.require_bare_key(label, 'cohort label')
        except ValueError:
            logging.warning(f'Skipping dataset {label!r}: not a bare TOML key, cannot be used as a cohort label')
            continue

        result = query_fn(ANALYSES_QUERY, {'datasetName': label})
        analyses = result.get('project', {}).get('analyses', [])
        analysis = latest_analysis(analyses, seq_type, label)
        if analysis is None:
            continue

        cohorts.append(
            Cohort(
                label=label,
                uri=analysis['output'],
                analysis_id=int(analysis['id']),
                timestamp=analysis.get('timestampCompleted'),
            ),
        )

    if not cohorts:
        raise DiscoveryError(
            f'discovery found no cohorts for seq_type={seq_type!r}; check the is_seqr eligibility filter '
            'and whether any dataset has a completed qc analysis of this sequencing type',
        )

    return Manifest(seq_type=seq_type, generated=generated, cohorts=tuple(cohorts))
