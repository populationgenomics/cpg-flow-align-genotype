"""The calibration spec - what to gate, in which direction, at what value.

A spec plus a manifest plus a value cache is everything needed to reproduce a
calibration, which is why thresholds live here rather than in editable Python
constants (the pain point in the original testing_scripts/ workflow).

Validation encodes the threshold-design rules structurally: a cohort-relative warn
tier always keeps an absolute fail gate behind it, and never coexists with an
absolute warn (the relative tier *is* the warn tier).
"""

from dataclasses import dataclass, replace
from pathlib import Path
from typing import Any

from align_genotype.qc_calibration import tomlio

DIRECTIONS = ('min', 'max')
UNITS = ('x', 'frac', '%')

_METRIC_KEYS = {'direction', 'unit', 'gated', 'fail', 'warn', 'relative', 'reviewed', 'rationale'}
_RELATIVE_KEYS = {'k', 'min_cohort'}


class SpecError(ValueError):
    """The calibration spec is missing something or is internally inconsistent."""


@dataclass(frozen=True)
class RelativeSpec:
    """Cohort-relative (MAD) warn tier settings; mirrors [qc_thresholds.*.relative.*]."""

    k: float = 3.5
    min_cohort: int = 50


@dataclass(frozen=True)
class MetricSpec:
    """One metric's direction, candidate thresholds and review state.

    `direction`: 'min' = higher is better (flag below the threshold); 'max' = lower is
    better (flag above it). `gated` metrics are enforced and are fatal if a cohort is
    missing them; un-gated metrics are surveyed and profiled for the record only.
    `reviewed` is False on anything `suggest` seeded and not yet signed off; emit-config
    refuses to emit an unreviewed gated metric.
    """

    key: str
    direction: str
    unit: str = 'frac'
    gated: bool = True
    fail: float | None = None
    warn: float | None = None
    relative: RelativeSpec | None = None
    reviewed: bool = False
    rationale: str = ''


@dataclass(frozen=True)
class CalibrationSpec:
    seq_type: str
    cache: str
    metrics: tuple[MetricSpec, ...]

    @property
    def gated(self) -> tuple[MetricSpec, ...]:
        return tuple(m for m in self.metrics if m.gated)

    @property
    def metric_keys(self) -> tuple[str, ...]:
        return tuple(m.key for m in self.metrics)

    def metric(self, key: str) -> MetricSpec:
        for m in self.metrics:
            if m.key == key:
                return m
        raise KeyError(f'{key!r} is not in this calibration spec')

    def with_metric(self, replace_key: str, **changes: Any) -> 'CalibrationSpec':
        """Return a copy with one metric's fields replaced; order preserved."""
        updated = tuple(replace(m, **changes) if m.key == replace_key else m for m in self.metrics)
        return replace(self, metrics=updated)


def _parse_relative(key: str, raw: Any) -> RelativeSpec:
    if not isinstance(raw, dict):
        raise SpecError(f'metric {key!r}: [metrics.{key}.relative] must be a table')
    if unknown := sorted(set(raw) - _RELATIVE_KEYS):
        raise SpecError(f'metric {key!r}: unknown relative key(s) {unknown}; expected {sorted(_RELATIVE_KEYS)}')
    return RelativeSpec(k=float(raw.get('k', 3.5)), min_cohort=int(raw.get('min_cohort', 50)))


def _parse_metric(key: str, raw: Any) -> MetricSpec:
    try:
        tomlio.require_bare_key(key, 'metric key')
    except ValueError as exc:
        raise SpecError(str(exc)) from exc
    if not isinstance(raw, dict):
        raise SpecError(f'metric {key!r}: [metrics.{key}] must be a table')
    if unknown := sorted(set(raw) - _METRIC_KEYS):
        raise SpecError(f'metric {key!r}: unknown key(s) {unknown}; expected {sorted(_METRIC_KEYS)}')
    if 'direction' not in raw:
        raise SpecError(f'metric {key!r}: missing required key: direction')
    direction = raw['direction']
    if direction not in DIRECTIONS:
        raise SpecError(f"metric {key!r}: direction must be 'min' or 'max', got {direction!r}")
    unit = raw.get('unit', 'frac')
    if unit not in UNITS:
        raise SpecError(f'metric {key!r}: unit must be one of {list(UNITS)}, got {unit!r}')

    metric = MetricSpec(
        key=key,
        direction=direction,
        unit=unit,
        gated=bool(raw.get('gated', True)),
        fail=raw.get('fail'),
        warn=raw.get('warn'),
        relative=_parse_relative(key, raw['relative']) if 'relative' in raw else None,
        reviewed=bool(raw.get('reviewed', False)),
        rationale=str(raw.get('rationale', '')),
    )
    _validate_metric(metric)
    return metric


def _validate_metric(metric: MetricSpec) -> None:
    key = metric.key
    has_tier = metric.fail is not None or metric.warn is not None or metric.relative is not None
    if not metric.gated:
        if has_tier:
            raise SpecError(
                f'metric {key!r}: non-gated metric must not define fail, warn or relative - '
                f'set gated = true to enforce it, or remove the thresholds.',
            )
        return
    if not has_tier:
        raise SpecError(f'metric {key!r}: gated metric defines no fail, warn or relative tier')
    if metric.relative is not None:
        if metric.fail is None:
            raise SpecError(
                f'metric {key!r}: a relative tier requires an absolute fail gate behind it - '
                f'cohort-relative flagging is warn-only.',
            )
        if metric.warn is not None:
            raise SpecError(
                f'metric {key!r}: the relative tier replaces the absolute warn; remove `warn` '
                f'or remove the [metrics.{key}.relative] block.',
            )


def loads(text: str) -> CalibrationSpec:
    """Parse a calibration spec from a TOML string."""
    return _from_dict(tomlio.loads(text))


def load(path: str | Path) -> CalibrationSpec:
    """Load a calibration spec from a local or cloud path."""
    return _from_dict(tomlio.load_path(path))


def _from_dict(raw: dict[str, Any]) -> CalibrationSpec:
    for required in ('seq_type', 'cache'):
        if required not in raw:
            raise SpecError(f'calibration spec: missing required key: {required}')
    metrics_raw = raw.get('metrics', {})
    if not metrics_raw:
        raise SpecError('calibration spec defines no metrics; add at least one [metrics.<KEY>] table')
    return CalibrationSpec(
        seq_type=str(raw['seq_type']),
        cache=str(raw['cache']),
        metrics=tuple(_parse_metric(key, value) for key, value in metrics_raw.items()),
    )


def dumps(spec: CalibrationSpec) -> str:
    """Render a spec back to TOML, preserving metric order."""
    lines = [tomlio.fmt_kv('seq_type', spec.seq_type), tomlio.fmt_kv('cache', spec.cache)]
    for metric in spec.metrics:
        lines += ['', f'[metrics.{metric.key}]', tomlio.fmt_kv('direction', metric.direction)]
        lines.append(tomlio.fmt_kv('unit', metric.unit))
        lines.append(tomlio.fmt_kv('gated', metric.gated))
        if metric.fail is not None:
            lines.append(tomlio.fmt_kv('fail', metric.fail))
        if metric.warn is not None:
            lines.append(tomlio.fmt_kv('warn', metric.warn))
        lines.append(tomlio.fmt_kv('reviewed', metric.reviewed))
        lines.append(tomlio.fmt_kv('rationale', metric.rationale))
        # The relative sub-table must follow all of its parent's scalar keys.
        if metric.relative is not None:
            lines += [
                f'[metrics.{metric.key}.relative]',
                tomlio.fmt_kv('k', metric.relative.k),
                tomlio.fmt_kv('min_cohort', metric.relative.min_cohort),
            ]
    return '\n'.join(lines) + '\n'


def save(spec: CalibrationSpec, path: str | Path) -> None:
    """Write a spec to disk."""
    Path(path).write_text(dumps(spec))
