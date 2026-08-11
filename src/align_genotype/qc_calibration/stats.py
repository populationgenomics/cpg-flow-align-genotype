"""Numeric analysis over cached values - percentiles, flag rates, cohort-growth churn.

No I/O and no config: everything here is a pure function of arrays already in memory,
which is what lets the flagrates/mad tuning loop run in under a second. An operator
iterating on candidate thresholds is the primary use, and that loop has to be fast
enough to stay interactive.
"""

from dataclasses import dataclass

import numpy as np

from align_genotype.scripts import check_multiqc

PERCENTILES: tuple[int, ...] = (1, 5, 10, 25, 50, 75, 90, 95, 99)

# A healthy cohort should sit near 0% fail and single-digit % warn. Beyond these, the
# candidate threshold gets flagged for a second look - it is advice, not a rejection:
# "healthy cohort" is the operator's judgement, not a computable property.
FAIL_RATE_LIMIT = 0.02
WARN_RATE_LIMIT = 0.10


def percentiles(values: np.ndarray, pcts: tuple[int, ...] = PERCENTILES) -> tuple[float, ...]:
    """Percentiles of `values`, in the order given. Empty input yields an empty tuple."""
    if values.size == 0:
        return ()
    return tuple(float(p) for p in np.percentile(values, pcts))


def breach(values: np.ndarray, threshold: float, direction: str) -> np.ndarray:
    """Boolean mask of samples on the bad side of `threshold`.

    'min' = higher is better, so low values breach; 'max' = lower is better. Strict
    comparisons, matching the `<` and `>` production uses - a value exactly on the
    threshold is not flagged.
    """
    return values < threshold if direction == 'min' else values > threshold


def flag_rates(
    values: np.ndarray,
    direction: str,
    fail: float | None = None,
    warn: float | None = None,
) -> tuple[float, float]:
    """``(fail_rate, warn_rate)`` as fractions, warn excluding samples already failing.

    Mirrors production, which evaluates fail before warn and records one flag per
    metric at the worst tier. Counting failing samples in the warn rate would
    overstate what an operator actually sees in the report. Empty input yields
    ``(nan, nan)`` so callers can distinguish "no samples" from "nothing flagged".
    """
    if values.size == 0:
        return (float('nan'), float('nan'))
    is_fail = breach(values, fail, direction) if fail is not None else np.zeros(values.size, dtype=bool)
    is_warn = breach(values, warn, direction) & ~is_fail if warn is not None else np.zeros(values.size, dtype=bool)
    return float(is_fail.mean()), float(is_warn.mean())


def needs_review(fail_rate: float, warn_rate: float) -> bool:
    """Whether a candidate threshold's flag rate is outside the healthy-cohort guardrail.

    NaN (an absent metric) is not "needs review" - there is nothing to look at.
    """
    if np.isnan(fail_rate) or np.isnan(warn_rate):
        return False
    return bool(fail_rate > FAIL_RATE_LIMIT or warn_rate > WARN_RATE_LIMIT)


@dataclass(frozen=True)
class ChurnResult:
    """How a cohort-relative threshold moved, and who changed status because of it."""

    threshold_before: float
    threshold_after: float
    n_initial: int
    flagged_before: int
    flagged_after: int
    flips: int

    @property
    def flip_rate(self) -> float:
        return self.flips / self.n_initial if self.n_initial else 0.0


def churn(initial: np.ndarray, grown: np.ndarray, direction: str, k: float) -> ChurnResult | None:
    """Re-score the *initial* samples against the *grown* cohort's threshold.

    `flips` counts samples whose flag status changes purely because the cohort grew -
    each one would be a spurious "updated" flag in the database, which is the cost that
    decides whether a cohort-relative tier is safe to adopt. Returns None when either
    cohort has a degenerate (zero) MAD, since no threshold exists to compare.

    Thresholds are rounded to 4 dp exactly as production does, so the simulation
    measures the churn operators would actually see rather than sub-0.0001 jitter.
    """
    before_threshold = check_multiqc.robust_threshold(list(initial), direction, k)
    after_threshold = check_multiqc.robust_threshold(list(grown), direction, k)
    if before_threshold is None or after_threshold is None:
        return None
    before_threshold = round(before_threshold, 4)
    after_threshold = round(after_threshold, 4)
    before = breach(initial, before_threshold, direction)
    after = breach(initial, after_threshold, direction)
    return ChurnResult(
        threshold_before=before_threshold,
        threshold_after=after_threshold,
        n_initial=int(initial.size),
        flagged_before=int(before.sum()),
        flagged_after=int(after.sum()),
        flips=int((before != after).sum()),
    )
