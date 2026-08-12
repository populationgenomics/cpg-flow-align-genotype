"""Candidate fixed thresholds, seeded from the per-dataset distribution tails.

These are starting points, not answers, and the report says so. What a percentile cannot
supply is the judgement a threshold actually turns on: preserving the lab's intent for a
hard gate unless the data clearly contradicts it, refusing to hard-fail on metrics that
track ancestry, biology or chemistry rather than sample quality, and overriding the lab
only with a measured count of good sequencing groups their line would have discarded.

Which tail seeds which tier: for `min` metrics (higher is better) the bad values sit in
the low tail, so `fail` comes from the lowest per-dataset p1 and `warn` from the lowest
p5. For `max` metrics it is the high tail - highest p99 for `fail`, highest p95 for
`warn`. Taking the worst per-dataset value rather than pooling means a candidate does not
already flag a large slice of the most marginal dataset in the set.
"""

from dataclasses import dataclass

import numpy as np

from align_genotype.qc_calibration.settings import MetricSpec
from align_genotype.qc_calibration.values import MetricValues

_SEED_PERCENTILES: dict[str, dict[str, int]] = {
    'min': {'fail': 1, 'warn': 5},
    'max': {'fail': 99, 'warn': 95},
}

# 'x' (coverage) and '%' round to whole numbers for readability; 'frac' keeps 2 dp, since
# a fraction rounded to an integer would collapse to 0 or 1.
_INTEGER_UNITS = frozenset({'x', '%'})
_FRACTION_DP = 2


@dataclass(frozen=True)
class Candidate:
    """One metric's candidate thresholds and the evidence behind them."""

    metric: str
    fail: float | None
    warn: float | None
    basis: str


def _round_for_unit(value: float, unit: str) -> float | int:
    """Round a raw percentile to a sensible reading for `unit`, as a Python builtin.

    `round(np.float64(...), 2)` returns another `np.float64`, which serialises to JSON
    only by accident, so the fraction branch casts first. The integer branch calls
    `round` with no `ndigits` on a plain float, which returns a genuine `int`.
    """
    if unit in _INTEGER_UNITS:
        return round(float(value))
    return round(float(value), _FRACTION_DP)


def _tail(by_dataset: dict[str, MetricValues], metric: MetricSpec, percentile: int) -> float | None:
    """The worst per-dataset value of `percentile`, or None if no dataset has data.

    'Worst' is direction-dependent: the lowest per-dataset percentile for a `min` metric,
    the highest for a `max` one.
    """
    readings = [
        float(np.percentile(values, percentile))
        for metric_values in by_dataset.values()
        if (values := metric_values.array).size
    ]
    if not readings:
        return None
    return min(readings) if metric.direction == 'min' else max(readings)


def candidate(by_dataset: dict[str, MetricValues], metric: MetricSpec) -> Candidate | None:
    """Candidate thresholds for one metric, or None when no dataset has any values.

    A metric with a relative tier gets a `fail` candidate only: its warn tier is derived
    per run from the dataset's own spread, and an absolute warn alongside it would
    double-flag the same values and disagree about which threshold was breached.
    """
    pcts = _SEED_PERCENTILES[metric.direction]
    fail_raw = _tail(by_dataset, metric, pcts['fail'])
    if fail_raw is None:
        return None

    n_with_data = sum(1 for v in by_dataset.values() if v.array.size)
    plural = '' if n_with_data == 1 else 's'
    basis = f'worst per-dataset p{pcts["fail"]}={fail_raw:.4g} across {n_with_data} dataset{plural}'

    warn: float | None = None
    if metric.relative:
        basis += '; warn tier is dataset-relative, so no absolute warn is proposed'
    else:
        warn_raw = _tail(by_dataset, metric, pcts['warn'])
        # Unreachable: `_tail` walks the same per-dataset arrays at a different
        # percentile, so a dataset with fail-percentile data always has warn-percentile
        # data. Raised rather than asserted, since ruff bans `assert` outside test/**.
        if warn_raw is None:
            raise AssertionError('a dataset with fail-percentile data must also have warn-percentile data')
        warn = _round_for_unit(warn_raw, metric.unit)
        basis += f', warn p{pcts["warn"]}={warn_raw:.4g}'

    return Candidate(
        metric=metric.key,
        fail=_round_for_unit(fail_raw, metric.unit),
        warn=warn,
        basis=basis,
    )
