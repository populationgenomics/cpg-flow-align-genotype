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

# Ported verbatim from fetch_multiqc_json_paths.py: 'test' not in name and 'training'
# not in name and 'seqr' not in name and meta.get('is_seqr', False). That script carried
# no comment, so only the behaviour below is load-bearing and deliberately preserved;
# the *probable* reason - a guess, not documented original intent - is that these
# substrings flag seqr-loader staging/test/training projects that also carry
# `is_seqr: true`, distinct from `is_seqr` marking a real source dataset. Excluding by
# name on top of requiring `is_seqr` therefore is not a contradiction, so do not
# "simplify" this to just `is_seqr`.
#
# Despite the name, this is substring matching, not whole-token matching: any dataset
# whose name contains one of these anywhere - not just as a hyphen-delimited token - is
# excluded.
EXCLUDED_NAME_SUBSTRINGS = ('test', 'training', 'seqr')


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
    """Whether `project` is a real source dataset, not a staging/derivative project.

    Tolerates malformed rows - `meta` that isn't a dict, a missing/`None` `name` - by
    treating them as ineligible rather than raising: this talks to a live external
    service across dozens of datasets, and one bad row shouldn't take the run down.
    """
    meta = project.get('meta')
    if not isinstance(meta, dict) or not meta.get('is_seqr', False):
        return False
    name = project.get('name') or ''
    if any(substring in name for substring in EXCLUDED_NAME_SUBSTRINGS):
        logging.info(
            f'discovery: excluding dataset {name!r}: name matches excluded substring(s) {EXCLUDED_NAME_SUBSTRINGS!r}'
        )
        return False
    return True


def latest_analysis(
    analyses: list[dict[str, Any]],
    seq_type: str,
    dataset_label: str,
    skipped_timestamps: list[tuple[str, Any]] | None = None,
) -> dict[str, Any] | None:
    """The most recently completed analysis of `seq_type` that produced an output.

    Selection is by `timestampCompleted`, not by id or list order: a higher id with an
    older timestamp must lose. This talks to a live external service, so malformed rows
    are expected: an analysis with a missing or non-string `timestampCompleted` can't be
    ordered against the others, so it is excluded rather than compared - better to
    silently drop a candidate we can't rank than to silently pick the wrong one, or crash
    with a `TypeError`/`KeyError` instead of a clean `DiscoveryError`. The resulting
    failure mode is "a slightly older report of the same cohort" - a data-currency
    problem the manifest's recorded `timestamp` makes inspectable, and one an operator
    can fix by pinning a different analysis - not a wrong-answer problem, so skipping
    rather than failing loudly is not an over-correction.

    `dataset_label` names the dataset in the warning logged for each skip; it is
    required because a defaulted parameter whose only effect is a log message is a
    smell, and the single production caller always has a real label to pass. If
    `skipped_timestamps` is given, `(dataset_label, analysis id)` is appended to it for
    each skip, so a caller iterating many datasets can emit one combined summary instead
    of relying on individual warnings not to scroll past unnoticed.
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
            if skipped_timestamps is not None:
                skipped_timestamps.append((dataset_label, analysis.get('id')))
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
    skipped_timestamps: list[tuple[str, Any]] = []
    for project in eligible:
        label = project.get('dataset') or project.get('name')
        if not label:
            logging.warning('discovery: skipping a project with no dataset/name to use as a cohort label')
            continue
        try:
            tomlio.require_bare_key(label, 'cohort label')
        except ValueError:
            logging.warning(f'Skipping dataset {label!r}: not a bare TOML key, cannot be used as a cohort label')
            continue

        result = query_fn(ANALYSES_QUERY, {'datasetName': label})
        # `result['project']` can be present-but-null (e.g. no read access to that
        # project), so `.get('project', {})` is not enough - a `None` value wins over
        # the default and `.get('analyses', ...)` would then crash on `None`.
        analyses = (result.get('project') or {}).get('analyses') or []
        analysis = latest_analysis(analyses, seq_type, label, skipped_timestamps)
        if analysis is None:
            logging.info(f'{label}: no completed {seq_type} QC analysis with an output; skipping')
            continue

        cohorts.append(
            Cohort(
                label=label,
                uri=analysis['output'],
                analysis_id=int(analysis['id']),
                # `latest_analysis` only ever returns candidates with a string
                # `timestampCompleted`, so indexing (not `.get`) documents that invariant.
                timestamp=analysis['timestampCompleted'],
            ),
        )

    if skipped_timestamps:
        distinct_datasets = {label for label, _analysis_id in skipped_timestamps}
        logging.warning(
            f'discovery: {len(skipped_timestamps)} analyses across {len(distinct_datasets)} datasets had '
            'unusable timestamps and were not considered',
        )

    if not cohorts:
        raise DiscoveryError(
            f'discovery found no cohorts for seq_type={seq_type!r}; check the is_seqr eligibility filter, '
            f'the {EXCLUDED_NAME_SUBSTRINGS!r} name-substring exclusions, and whether any dataset has a '
            'completed qc analysis of this sequencing type',
        )

    return Manifest(seq_type=seq_type, generated=generated, cohorts=tuple(cohorts))
