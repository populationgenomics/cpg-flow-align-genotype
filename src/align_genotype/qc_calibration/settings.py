"""Calibration settings, read from the `[qc_calibration]` config block.

Metric lists are per sequencing type because the Picard module differs: genome uses
`CollectWgsMetrics`, exome `CollectHsMetrics`. That is why the candidate metric set is a
config parameter and not a Python constant - an exome spec copied onto a genome run would
gate keys that do not exist.
"""

from dataclasses import dataclass, field
from typing import Any

from cpg_utils import config

from align_genotype.scripts import check_multiqc

DIRECTIONS = ('min', 'max')
UNITS = ('x', 'frac', '%')

_METRIC_KEYS = {'direction', 'unit', 'relative'}


class SettingsError(ValueError):
    """The `[qc_calibration]` config block is missing something or is inconsistent."""


@dataclass(frozen=True)
class MetricSpec:
    """One metric's direction, display unit and whether to evaluate a relative tier.

    `direction`: 'min' = higher is better, so a value below the threshold is flagged;
    'max' is the mirror. `unit` affects rounding and display only.
    """

    key: str
    direction: str
    unit: str = 'frac'
    relative: bool = False


@dataclass(frozen=True)
class Bars:
    """The advisory bars a relative tier's verdict is judged against.

    Growth and merge churn have separate bars because they model different things.
    Growth - a before-slice re-scored against the threshold the whole dataset produces -
    is a forecast: datasets accrete sequencing groups over time, which is what a shipped
    tier actually faces. Merge - two whole projects pooled into one run - is a stress
    test, not a prediction of anything scheduled. Hence the looser merge bar.
    """

    max_warn_rate: float = 0.10
    max_growth_churn: float = 0.02
    max_merge_churn: float = 0.05


@dataclass(frozen=True)
class CalibrationSettings:
    seq_type: str
    metrics: tuple[MetricSpec, ...]
    k: float = 3.5
    min_samples: int = 50
    bars: Bars = field(default_factory=Bars)

    @property
    def metric_keys(self) -> tuple[str, ...]:
        return tuple(m.key for m in self.metrics)

    @property
    def relative_metrics(self) -> tuple[MetricSpec, ...]:
        return tuple(m for m in self.metrics if m.relative)

    def metric(self, key: str) -> MetricSpec:
        for m in self.metrics:
            if m.key == key:
                return m
        raise KeyError(f'{key!r} is not a configured calibration metric; known: {list(self.metric_keys)}')


def parse_metric(key: str, raw: Any) -> MetricSpec:
    """Validate one `[qc_calibration.<seq_type>.metrics.<KEY>]` table."""
    if not isinstance(raw, dict):
        raise SettingsError(f'metric {key!r}: must be a table of direction/unit/relative keys')
    if unknown := sorted(set(raw) - _METRIC_KEYS):
        raise SettingsError(f'metric {key!r}: unknown key(s) {unknown}; expected {sorted(_METRIC_KEYS)}')
    if 'direction' not in raw:
        raise SettingsError(f'metric {key!r}: missing required key: direction')
    direction = raw['direction']
    if direction not in DIRECTIONS:
        raise SettingsError(f"metric {key!r}: direction must be 'min' or 'max', got {direction!r}")
    unit = raw.get('unit', 'frac')
    if unit not in UNITS:
        raise SettingsError(f'metric {key!r}: unit must be one of {list(UNITS)}, got {unit!r}')
    relative = raw.get('relative', False)
    # `bool('false')` is True, so a quoted boolean must be rejected rather than coerced.
    if not isinstance(relative, bool):
        raise SettingsError(f'metric {key!r}: relative must be true or false, got {relative!r}')
    return MetricSpec(key=key, direction=direction, unit=unit, relative=relative)


def enabled() -> bool:
    """Whether the calibration stages should queue any jobs at all.

    Defaults to False so that the stages sit inert in the DAG of an ordinary production
    run. They carry no `required_stages` dependency, so without this they would queue a
    job per dataset and register Metamist analyses on every invocation.
    """
    return bool(config.config_retrieve(['qc_calibration', 'enabled'], False))


def load() -> CalibrationSettings:
    """Read the settings for the run's sequencing type."""
    seq_type = config.config_retrieve(['workflow', 'sequencing_type'])
    raw_metrics = config.config_retrieve(['qc_calibration', seq_type, 'metrics'], {})
    if not raw_metrics:
        raise SettingsError(
            f'no calibration metrics configured for sequencing type {seq_type!r}; '
            f'add [qc_calibration.{seq_type}.metrics.<KEY>] tables',
        )
    return CalibrationSettings(
        seq_type=seq_type,
        metrics=tuple(parse_metric(key, value) for key, value in raw_metrics.items()),
        k=float(config.config_retrieve(['qc_calibration', 'k'], 3.5)),
        min_samples=int(config.config_retrieve(['qc_calibration', 'min_samples'], 50)),
        bars=Bars(
            max_warn_rate=float(config.config_retrieve(['qc_calibration', 'max_warn_rate'], 0.10)),
            max_growth_churn=float(config.config_retrieve(['qc_calibration', 'max_growth_churn'], 0.02)),
            max_merge_churn=float(config.config_retrieve(['qc_calibration', 'max_merge_churn'], 0.05)),
        ),
    )


def current_thresholds(seq_type: str) -> dict[str, dict[str, float]]:
    """The shipped absolute thresholds, as `{metric: {severity: value}}`.

    Read through `check_multiqc.load_thresholds` rather than a second config reader, so
    the "current" column in the report is exactly what production enforces today.
    """
    by_direction = check_multiqc.load_thresholds(seq_type)
    return {metric: dict(tiers) for by_metric in by_direction.values() for metric, tiers in by_metric.items()}
