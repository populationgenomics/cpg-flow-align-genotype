"""Find each dataset's latest CramMultiQC report in Metamist.

`CramMultiQC` and `GvcfMultiQC` both register `analysis_type='qc'`, and each registers
one analysis per entry in `analysis_keys` - so a plain type-and-sequencing-type query
returns GVCF reports and HTML files alongside the CRAM JSON we want. CPG Flow merges the
stage name into analysis meta, so the stage filter runs server-side and the output suffix
picks the JSON of the two.

The query function is injected so the selection logic is testable with a fake, leaving
only the small adapter untested.

An analysis with an unrankable `timestampCompleted` is skipped with a per-row warning,
not aggregated into a summary the way the old cohort-wide discovery did. That aggregation
existed because the old module swept every project in Metamist, where dozens of scattered
warnings genuinely needed collecting into one line. This module is called once per
dataset in the operator's own `input_cohorts`, so the worst case is one warning per
dataset in a run the operator deliberately assembled - legible on its own in a driver
log - and reintroducing aggregation would mean either module-level mutable state or an
accumulator threaded through the stage layer, neither of which is worth it for a rare
condition at this scale.
"""

import functools
import logging
from collections.abc import Callable
from dataclasses import dataclass
from typing import Any

from metamist.graphql import gql, query

# The CramMultiQC output we want; the stage also registers its HTML under the same type.
JSON_SUFFIX = 'multiqc_data.json'
CRAM_MULTIQC_STAGE = 'CramMultiQC'

QueryFn = Callable[[str, dict[str, Any] | None], dict[str, Any]]

# `meta` is an opaque `JSON` scalar in Metamist's schema - there is no typed filter
# object for it (unlike `status`/`type`, which have `{eq: ...}`-style filter inputs), so
# neither the flat key-value form used here nor a nested `{'eq': ...}` form is
# schema-enforced; the server just receives raw JSON. This flat form matches the one
# other place in this repo that filters analyses by stage,
# `scripts/build_vntyper_index.py`'s `REPORT_QUERY`, which is shipped and working -
# don't "modernise" this to the nested form without re-confirming against a live
# Metamist first.
ANALYSES_QUERY = gql(
    """
    query CramMultiqc($dataset: String!, $analysisType: String!, $metaFilter: JSON) {
        project(name: $dataset) {
            analyses(status: {eq: COMPLETED}, type: {eq: $analysisType}, meta: $metaFilter) {
                id
                output
                timestampCompleted
            }
        }
    }
    """,
)


@dataclass(frozen=True)
class MultiqcReport:
    """One dataset's MultiQC report, and the analysis it was registered under."""

    dataset: str
    uri: str
    analysis_id: int
    timestamp: str


def default_query(query_text: str, variables: dict[str, Any] | None = None) -> dict[str, Any]:
    return query(query_text, variables=variables or {})


def latest_cram_multiqc(
    dataset: str,
    seq_type: str,
    query_fn: QueryFn = default_query,
) -> MultiqcReport | None:
    """The newest completed CramMultiQC JSON for `dataset`, or None if there is none."""
    result = query_fn(
        ANALYSES_QUERY,
        {
            'dataset': dataset,
            'analysisType': 'qc',
            'metaFilter': {'stage': CRAM_MULTIQC_STAGE, 'sequencing_type': seq_type},
        },
    )
    # `project` can be present-but-null - no read access, say - so `.get('project', {})`
    # is not enough: the null wins over the default and the next `.get` crashes on None.
    analyses = (result.get('project') or {}).get('analyses') or []

    candidates = []
    for analysis in analyses:
        output = analysis.get('output')
        if not output or not str(output).endswith(JSON_SUFFIX):
            continue
        timestamp = analysis.get('timestampCompleted')
        if not isinstance(timestamp, str):
            # Cannot be ordered against the others. Dropping a candidate we cannot rank
            # beats picking the wrong one, and the cost is at worst a slightly older
            # report - inspectable, because the report records what it used.
            logging.warning(
                f'{dataset}: analysis {analysis.get("id")!r} has no usable timestampCompleted '
                f'({timestamp!r}); excluding it from latest-report selection',
            )
            continue
        candidates.append(analysis)

    if not candidates:
        return None
    newest = max(candidates, key=lambda a: a['timestampCompleted'])
    return MultiqcReport(
        dataset=dataset,
        uri=str(newest['output']),
        analysis_id=int(newest['id']),
        timestamp=newest['timestampCompleted'],
    )


@functools.cache
def cached_latest_cram_multiqc(dataset: str, seq_type: str) -> MultiqcReport | None:
    """`latest_cram_multiqc`, memoised for the driver.

    `expected_outputs` keys the per-dataset output path on the analysis ID and is called
    repeatedly during DAG assembly, so without this each call would be a GraphQL round
    trip.

    A `None` result (no report yet for this dataset/seq_type) is cached too, for the life
    of the process - safe for one-shot DAG assembly, where a dataset with no report at
    first lookup has no report for the rest of that assembly. Anything longer-lived that
    reuses this function (a service, a long-running loop) would need an explicit
    invalidation strategy, since a report completing after the first lookup would
    otherwise never be seen.
    """
    return latest_cram_multiqc(dataset, seq_type)
