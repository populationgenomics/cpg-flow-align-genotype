"""Distil one MultiQC report into one dataset's values file.

Extraction goes through `check_multiqc.normalise_sections` and
`check_multiqc.gather_metric_values`, so calibration sees precisely what enforcement
sees. A metric configured but absent from the report is recorded present-and-empty
rather than omitted: that distinction is what the report's presence matrix is built
from, and a silently absent metric is how a gate goes inert.
"""

import math
from typing import Any

from align_genotype.qc_calibration.settings import CalibrationSettings
from align_genotype.qc_calibration.values import DatasetValues, MetricValues
from align_genotype.scripts import check_multiqc


class ExtractError(RuntimeError):
    """A MultiQC report could not be read, or holds nothing to extract."""


def extract(
    document: Any,
    settings: CalibrationSettings,
    *,
    dataset: str,
    analysis_id: int,
    timestamp: str,
    uri: str,
    generated: str,
) -> DatasetValues:
    """Extract every configured metric from one parsed MultiQC document."""
    # `[]`, `null`, `42` and `"str"` are all valid JSON, so a truncated or wrong-file URI
    # can parse cleanly and then fail on `.get` with a bare AttributeError naming neither
    # the dataset nor the path.
    if not isinstance(document, dict):
        raise ExtractError(f'{dataset}: report is a {type(document).__name__}, not an object, in {uri}')

    version = str(document.get('config_version', 'unknown'))
    sections = check_multiqc.normalise_sections(document.get('report_general_stats_data'))
    if not sections:
        raise ExtractError(f'{dataset}: no usable report_general_stats_data (multiqc {version}) in {uri}')

    metrics: dict[str, MetricValues] = {}
    for metric in settings.metrics:
        entries, n_non_numeric = check_multiqc.gather_metric_values(sections, metric.key)
        finite = tuple(
            # Production derives a sequencing group ID the same way: MultiQC runs with
            # --replace-names against the dataset's rich ID map, so a sample key can read
            # `CPG1|EXTID`.
            (section, sample.split('|', 1)[0], value)
            for section, sample, value in entries
            if math.isfinite(value)
        )
        # Two disjoint kinds of loss - values float() refused, and values it accepted that
        # came back nan/inf - so summing them cannot double-count.
        metrics[metric.key] = MetricValues(
            entries=finite,
            n_dropped=n_non_numeric + (len(entries) - len(finite)),
        )

    return DatasetValues(
        dataset=dataset,
        seq_type=settings.seq_type,
        analysis_id=analysis_id,
        timestamp=timestamp,
        uri=uri,
        multiqc_version=version,
        generated=generated,
        n_sequencing_groups=len(
            {sample.split('|', 1)[0] for section in sections.values() for sample in section},
        ),
        section_sizes={name: len(section) for name, section in sections.items()},
        metrics=metrics,
    )
