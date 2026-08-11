"""Survey and extract in a single parse of each MultiQC report.

The survey answers "is every metric I intend to gate actually present in general
stats, in every cohort?". Getting that wrong is how a gate goes silently inert - the
genome reads-mapped gate was configured on ``PCT_PF_READS_ALIGNED``, which MultiQC only
writes to ``report_saved_raw_data``, so it checked nothing at all. A gated metric
missing from any cohort is therefore fatal here, and the cache it writes is marked
incomplete so nothing downstream will run on it.

Surveying and extracting are one pass because both need the same expensive read: these
reports run to hundreds of megabytes, and re-reading them to answer a second question
would double the only slow step in the workflow.

Extraction goes through ``check_multiqc.normalise_sections`` and
``check_multiqc.gather_metric_values`` so calibration sees precisely what enforcement
will see.
"""

import gc
import json
import math
from dataclasses import dataclass
from datetime import datetime, timezone
from typing import Any

from align_genotype.qc_calibration.cache import CohortValues, ValueCache
from align_genotype.qc_calibration.manifest import Cohort, Manifest
from align_genotype.qc_calibration.spec import CalibrationSpec
from align_genotype.scripts import check_multiqc


class CollectError(RuntimeError):
    """A report could not be read, or is unusable for calibration."""


@dataclass(frozen=True)
class SurveyRow:
    """What one report looked like: provenance, shape, and where each metric lives."""

    label: str
    uri: str
    multiqc_version: str
    shape: str
    n_samples: int
    section_sizes: dict[str, int]
    where: dict[str, tuple[str, ...]]
    section_keys: dict[str, tuple[str, ...]]
    n_dropped: int
    n_dropped_by_metric: dict[str, int]


@dataclass(frozen=True)
class CollectResult:
    cache: ValueCache
    rows: tuple[SurveyRow, ...]
    missing_gated: dict[str, tuple[str, ...]]
    failures: tuple[tuple[str, str], ...]

    @property
    def ok(self) -> bool:
        return not self.missing_gated and not self.failures


def sections_carrying(sections: dict[str, Any], metric: str) -> tuple[str, ...]:
    """Section names in which at least one sample carries `metric`.

    Presence is judged on the key, not on whether its value parses as a number, so a
    metric that is present but entirely non-numeric still reads as present. Judging
    presence on numeric-ness would report a renamed key and an all-'?' column
    identically, and those need different fixes - rename the metric in the spec versus
    investigate why the pipeline emitted no numbers. The two stay distinguishable by
    reading this against the extracted values: an empty entry here means the key is
    absent or renamed, while a non-empty entry whose value list is empty means the key
    is present but every value was unusable.
    """
    return tuple(
        name
        for name, section in sections.items()
        if any(metric in values for values in section.values() if isinstance(values, dict))
    )


def _all_keys(section: dict[str, Any]) -> tuple[str, ...]:
    keys: set[str] = set()
    for values in section.values():
        if isinstance(values, dict):
            keys |= set(values)
    return tuple(sorted(keys))


def _extract(sections: dict[str, Any], spec: CalibrationSpec) -> tuple[dict[str, list[float]], dict[str, int]]:
    """Each metric's finite values, and how many values it lost, per metric.

    Drops are attributed per metric rather than summed, because a single total can't be
    unpicked: two metrics each losing one value and one metric losing two are both
    reported as 2. A metric where most samples are Picard's ``'?'`` placeholder is a
    real signal - production logs it per metric in ``check_multiqc`` - and it needs to
    survive into the survey report rather than being averaged away.

    A sample carrying `metric` in two sections contributes its value twice. That is not
    hypothetical: MultiQC 1.33 can split one tool across sections (``picard_1`` and
    ``picard_4`` both sit in the Picard namespace), and `normalise_sections` keeps them
    separate on purpose. The duplication is kept rather than de-duplicated here because
    it is what production does - ``_relative_flags_for_metric`` feeds the same doubled
    list to ``robust_threshold``, so de-duplicating would make calibration disagree with
    enforcement, which is the one thing this module exists to prevent. The consequence
    for callers: ``len(values[metric])`` counts *values*, not samples, and may exceed
    ``n_samples``. Anything treating it as a cohort size, or taking a percentile or MAD
    from it, is weighting duplicated samples twice.
    """
    values: dict[str, list[float]] = {}
    n_dropped_by_metric: dict[str, int] = {}
    for metric in spec.metrics:
        entries, n_non_numeric = check_multiqc.gather_metric_values(sections, metric.key)
        finite = [value for _, _, value in entries if math.isfinite(value)]
        values[metric.key] = finite
        # Two disjoint kinds of loss: values float() refused, and values it accepted
        # that came back nan/inf. Summing them can't double-count.
        n_dropped_by_metric[metric.key] = n_non_numeric + (len(entries) - len(finite))
    return values, n_dropped_by_metric


def collect_cohort(cohort: Cohort, spec: CalibrationSpec) -> tuple[CohortValues, SurveyRow]:
    """Parse one report and distil it to values plus a survey row."""
    from cpg_utils import to_path  # noqa: PLC0415

    with to_path(cohort.uri).open() as f:
        document = json.load(f)

    # `[]`, `null`, `42` and `"str"` are all valid JSON, so a truncated or wrong-file
    # URI can parse cleanly and then fail on `.get` with a bare AttributeError naming
    # neither the cohort nor the path. Check the shape and say which report it was.
    if not isinstance(document, dict):
        raise CollectError(f'{cohort.label}: report is a {type(document).__name__}, not an object, in {cohort.uri}')

    version = str(document.get('config_version', 'unknown'))
    raw = document.get('report_general_stats_data')
    shape = 'dict' if isinstance(raw, dict) else 'list' if isinstance(raw, list) else type(raw).__name__
    sections = check_multiqc.normalise_sections(raw)
    if not sections:
        raise CollectError(
            f'{cohort.label}: no usable report_general_stats_data (multiqc {version}, shape {shape}) in {cohort.uri}',
        )

    n_samples = len({sample for section in sections.values() for sample in section})
    values, n_dropped_by_metric = _extract(sections, spec)
    n_dropped = sum(n_dropped_by_metric.values())
    row = SurveyRow(
        label=cohort.label,
        uri=cohort.uri,
        multiqc_version=version,
        shape=shape,
        n_samples=n_samples,
        section_sizes={name: len(section) for name, section in sections.items()},
        where={metric.key: sections_carrying(sections, metric.key) for metric in spec.metrics},
        section_keys={name: _all_keys(section) for name, section in sections.items()},
        n_dropped=n_dropped,
        n_dropped_by_metric=n_dropped_by_metric,
    )
    cohort_values = CohortValues(
        label=cohort.label,
        n_samples=n_samples,
        multiqc_version=version,
        shape=shape,
        n_dropped=n_dropped,
        values=values,
    )
    return cohort_values, row


def collect_all(manifest: Manifest, spec: CalibrationSpec, generated: str | None = None) -> CollectResult:
    """Collect every cohort in `manifest`, one report at a time.

    A cohort that can't be read is recorded and the rest continue, so one bad URI
    doesn't waste a long run - but the result is not `ok`, and the caller must exit
    non-zero.

    That containment is why the per-cohort catch is deliberately broad. The failures
    worth surviving here are mostly not `OSError` or `ValueError`: a 403 on one
    dataset's bucket raises `google.api_core.exceptions.Forbidden`, an expired
    credential raises `google.auth.exceptions.RefreshError`, and cloudpathlib adds its
    own `MissingCredentialsError` and `NoStatError` - none of which share an ancestor
    with the two obvious ones. JSON also permits unbounded integer literals, so a
    400-digit number reaches `float()` and raises `OverflowError`, an `ArithmeticError`.
    The asymmetry settles it: over-catching costs one cohort recorded as a failure that
    an operator then reads, while under-catching throws away every report already parsed
    - a 403 at cohort 9 of 10 would discard the eight successes behind it. `Exception`
    is the boundary rather than `BaseException`, so Ctrl-C and `SystemExit` still abort
    the run immediately instead of being logged as ten cohort failures.
    """
    if manifest.seq_type != spec.seq_type:
        raise CollectError(f'manifest is for {manifest.seq_type!r} but the spec is for {spec.seq_type!r}')

    collected: list[CohortValues] = []
    rows: list[SurveyRow] = []
    missing_gated: dict[str, tuple[str, ...]] = {}
    failures: list[tuple[str, str]] = []

    for cohort in manifest.cohorts:
        try:
            values, row = collect_cohort(cohort, spec)
        except Exception as exc:  # noqa: BLE001 - one bad cohort must not cost the whole run
            failures.append((cohort.label, str(exc)))
            continue
        finally:
            gc.collect()  # release the parsed report before the next one is read

        collected.append(values)
        rows.append(row)
        if gaps := tuple(metric.key for metric in spec.gated if not row.where.get(metric.key)):
            missing_gated[cohort.label] = gaps

    cache = ValueCache(
        seq_type=spec.seq_type,
        generated=generated or datetime.now(tz=timezone.utc).isoformat(timespec='seconds'),
        complete=not missing_gated and not failures,
        metrics=spec.metric_keys,
        cohorts=tuple(collected),
    )
    return CollectResult(cache=cache, rows=tuple(rows), missing_gated=missing_gated, failures=tuple(failures))
