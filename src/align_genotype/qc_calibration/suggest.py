"""Seed a first draft of thresholds from the cohort distributions.

This module never writes a decision - it writes a starting point. Every seed lands
with `reviewed = false`, and `emit-config` refuses to emit an unreviewed gated metric,
so a percentile can propose but cannot itself reach production config. That interlock
exists because the judgement a calibration actually turns on is not something a
distribution can supply: preserving a lab's intent for a hard gate unless the data
clearly contradicts it, refusing to hard-fail on metrics that track ancestry, biology
or chemistry rather than sample quality, and sometimes overriding the lab outright -
genome duplication `fail` was relaxed from 25 to 40 once the data showed 25 would have
hard-failed a quarter to a third of legitimately higher-duplication preps. A seeded
number can't know any of that; a human has to look at it.

Which tail seeds which tier: for `min` metrics (higher is better) the bad samples sit
in the low tail, so `fail` comes from the lowest per-cohort p1 and `warn` from the
lowest per-cohort p5. For `max` metrics (lower is better) it's the high tail - the
highest p99 for `fail`, the highest p95 for `warn`. Taking the worst per-cohort value
rather than pooling or averaging means a seeded threshold doesn't already flag a large
slice of the most marginal cohort in the set.

A metric with a `relative` block gets a `fail` seed only: its warn tier is cohort-
derived by construction, and an absolute warn alongside it would double-flag the same
samples. `spec.with_metric` enforces this structurally - seeding a `warn` next to a
`relative` block raises `SpecError` rather than producing an inconsistent spec.
"""

from dataclasses import dataclass
from typing import Any

import numpy as np

from align_genotype.qc_calibration.cache import ValueCache
from align_genotype.qc_calibration.spec import CalibrationSpec, MetricSpec

_SEED_PERCENTILES: dict[str, dict[str, int]] = {
    'min': {'fail': 1, 'warn': 5},
    'max': {'fail': 99, 'warn': 95},
}

# 'x' (coverage) and '%' round to whole numbers for readability; 'frac' keeps 2 dp
# since a fraction rounded to an integer would collapse to 0 or 1.
_INTEGER_UNITS = {'x', '%'}
_FRACTION_DP = 2


@dataclass(frozen=True)
class Seeded:
    """One metric's newly-seeded thresholds, and the evidence behind them."""

    key: str
    fail: float | None
    warn: float | None
    evidence: str


def _round_for_unit(value: float, unit: str) -> float:
    """Round a raw percentile to a sensible reading for `unit`, as a plain Python scalar.

    `value` may already be a numpy float64 (every caller's is, since it comes straight
    out of `np.percentile`). `round(np.float64(...), 2)` returns another `np.float64`,
    which is why the `frac` branch goes through `float()` first rather than relying on
    `round` alone. The integer branch calls `round` with no `ndigits` on a plain
    `float`, which returns a genuine Python `int` - no separate `int(...)` cast needed.
    Both branches must hand back a Python builtin: these values are written to TOML,
    and the spec's own validation (`_require_number`) expects real numbers.
    """
    if unit in _INTEGER_UNITS:
        return round(float(value))
    return round(float(value), _FRACTION_DP)


def _tail(cache: ValueCache, metric: MetricSpec, percentile: int) -> tuple[float, list[tuple[str, float]]] | None:
    """The worst per-cohort value of `percentile` across every cohort with data.

    'Worst' is direction-dependent: for a `min` metric (higher is better) the worst
    case is the *lowest* per-cohort percentile, so a seeded threshold doesn't already
    flag a slice of the most marginal cohort's good samples. For a `max` metric (lower
    is better) the worst case is the *highest* per-cohort percentile, for the mirror
    reason. Cohorts with no data for this metric are skipped; returns None if none has
    any. Also returns the raw per-cohort readings, for the evidence string.
    """
    readings = []
    for label in cache.labels:
        values = cache.series(label, metric.key)
        if values.size == 0:
            continue
        readings.append((label, float(np.percentile(values, percentile))))
    if not readings:
        return None
    worst = min(readings, key=lambda r: r[1]) if metric.direction == 'min' else max(readings, key=lambda r: r[1])
    return worst[1], readings


def _seed_metric(cache: ValueCache, metric: MetricSpec) -> Seeded | None:
    """The candidate fail/warn seed for one gated metric, or None if it has no data."""
    pcts = _SEED_PERCENTILES[metric.direction]
    fail_tail = _tail(cache, metric, pcts['fail'])
    if fail_tail is None:
        return None
    fail_value, fail_readings = fail_tail
    fail = _round_for_unit(fail_value, metric.unit)

    n_cohorts = len(fail_readings)
    evidence_parts = [
        f'p{pcts["fail"]}={fail_value:.4g} (fail) across {n_cohorts} cohorts '
        f'[{", ".join(f"{label}={v:.4g}" for label, v in fail_readings)}]',
    ]

    warn = None
    if metric.relative is None:
        warn_value, warn_readings = _tail(cache, metric, pcts['warn'])  # same cohorts have data for both tails
        warn = _round_for_unit(warn_value, metric.unit)
        evidence_parts.append(
            f'p{pcts["warn"]}={warn_value:.4g} (warn) across {len(warn_readings)} cohorts '
            f'[{", ".join(f"{label}={v:.4g}" for label, v in warn_readings)}]',
        )

    return Seeded(key=metric.key, fail=fail, warn=warn, evidence='; '.join(evidence_parts))


def seed(cache: ValueCache, spec: CalibrationSpec) -> tuple[CalibrationSpec, tuple[Seeded, ...]]:
    """Seed fail/warn thresholds for every gated, unreviewed metric with data.

    Already-`reviewed` metrics and un-gated metrics are left untouched; a metric with
    no data in any cohort is skipped rather than seeded with a fabricated value. Every
    seed is applied with `reviewed = False` and a `rationale` recording the evidence,
    via `spec.with_metric`, which re-validates the changed metric - the same guard that
    stops a `relative` metric from ever being seeded with a `warn`. Returns a new spec
    (the input is never mutated) plus what changed, in spec order.
    """
    updated = spec
    seeded: list[Seeded] = []
    for metric in spec.gated:
        if metric.reviewed:
            continue
        result = _seed_metric(cache, metric)
        if result is None:
            continue
        changes: dict[str, Any] = {'fail': result.fail, 'reviewed': False, 'rationale': result.evidence}
        if metric.relative is None:
            changes['warn'] = result.warn
        updated = updated.with_metric(metric.key, **changes)
        seeded.append(result)
    return updated, tuple(seeded)
