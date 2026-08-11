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

The evidence text this module writes into `rationale` is aggregated across cohorts -
a cohort count and the percentile value, never a per-cohort reading or label. `emit`
pastes `rationale` verbatim into a file that is committed and pushed to a public repo,
so a label here would leak a real CPG dataset name into that file the moment an
operator reviews and signs off - exactly the point the review gate cannot catch, since
signing off is the expected next step, not a red flag. `emit.render` also checks this
independently before printing anything, in case a hand-edited rationale reintroduces
a label; this module's job is to make sure it never has to.
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

# Marks a rationale as machine-written so a re-run can safely refresh it. An operator's
# own prose - anything not starting with this marker, and any non-empty rationale
# `seed` didn't itself write - is left alone; only a metric with no rationale, or one
# still carrying a previous seed's text, gets overwritten on the next run.
#
# `_seed_metric`'s evidence text is written to read as the rest of this sentence (it
# opens with 'from cohort percentiles...'), not as a second, redundant "Seeded" of its
# own - the rendered rationale is `RATIONALE_MARKER + evidence`, read verbatim by an
# operator and later pasted into config_template.toml.
RATIONALE_MARKER = 'Seeded: '


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


def _tail(cache: ValueCache, metric: MetricSpec, percentile: int) -> float | None:
    """The worst per-cohort value of `percentile` across every cohort with data.

    'Worst' is direction-dependent: for a `min` metric (higher is better) the worst
    case is the *lowest* per-cohort percentile, so a seeded threshold doesn't already
    flag a slice of the most marginal cohort's good samples. For a `max` metric (lower
    is better) the worst case is the *highest* per-cohort percentile, for the mirror
    reason. Cohorts with no data for this metric are skipped; returns None if none has
    any. Returns a bare value, never a per-cohort breakdown - that breakdown would be
    keyed by cohort label, and this module must never carry a label into a `rationale`.
    """
    readings = [
        float(np.percentile(values, percentile))
        for label in cache.labels
        if (values := cache.series(label, metric.key)).size
    ]
    if not readings:
        return None
    return min(readings) if metric.direction == 'min' else max(readings)


def _n_with_data(cache: ValueCache, metric_key: str) -> int:
    """How many cohorts have at least one finite value for `metric_key`."""
    return sum(1 for label in cache.labels if cache.series(label, metric_key).size)


def _seed_metric(cache: ValueCache, metric: MetricSpec) -> Seeded | None:
    """The candidate fail/warn seed for one gated metric, or None if it has no data."""
    pcts = _SEED_PERCENTILES[metric.direction]
    fail_raw = _tail(cache, metric, pcts['fail'])
    if fail_raw is None:
        return None
    fail = _round_for_unit(fail_raw, metric.unit)
    n_cohorts = _n_with_data(cache, metric.key)
    evidence = f'from cohort percentiles across {n_cohorts} cohorts (fail p{pcts["fail"]}={fail_raw:.4g}'

    warn = None
    if metric.relative is None:
        warn_raw = _tail(cache, metric, pcts['warn'])
        if warn_raw is None:
            # Unreachable: `_tail` walks the same per-cohort series as the fail lookup
            # above, just at a different percentile, so a cohort with fail-percentile
            # data always has warn-percentile data too. Raised explicitly (not a plain
            # `assert`, which `ruff` bans outside `test/**` since it's stripped under
            # `-O`) so a violated invariant surfaces as a real exception, and so mypy
            # narrows `warn_raw` to `float` from here on.
            raise AssertionError('a cohort with fail-percentile data must also have warn-percentile data')
        warn = _round_for_unit(warn_raw, metric.unit)
        evidence += f', warn p{pcts["warn"]}={warn_raw:.4g}'
    evidence += ').'

    return Seeded(key=metric.key, fail=fail, warn=warn, evidence=evidence)


def seed(cache: ValueCache, spec: CalibrationSpec) -> tuple[CalibrationSpec, tuple[Seeded, ...]]:
    """Seed fail/warn thresholds for every gated, unreviewed metric with data.

    Already-`reviewed` metrics and un-gated metrics are left untouched; a metric with
    no data in any cohort is skipped rather than seeded with a fabricated value. Every
    seed refreshes `fail`/`warn` and `reviewed = False` via `spec.with_metric`, which
    re-validates the changed metric - the same guard that stops a `relative` metric from
    ever being seeded with a `warn`. `rationale` is only overwritten when it is empty or
    still carries a previous seed's `RATIONALE_MARKER` prefix, so an operator's own
    prose on a metric they've annotated but not yet reviewed survives a re-run. Returns
    a new spec (the input is never mutated) plus what changed, in spec order.
    """
    updated = spec
    seeded: list[Seeded] = []
    for metric in spec.gated:
        if metric.reviewed:
            continue
        result = _seed_metric(cache, metric)
        if result is None:
            continue
        changes: dict[str, Any] = {'fail': result.fail, 'reviewed': False}
        if metric.relative is None:
            changes['warn'] = result.warn
        if not metric.rationale or metric.rationale.startswith(RATIONALE_MARKER):
            changes['rationale'] = RATIONALE_MARKER + result.evidence
        updated = updated.with_metric(metric.key, **changes)
        seeded.append(result)
    return updated, tuple(seeded)
