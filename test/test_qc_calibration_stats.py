"""Unit tests for the calibration statistics."""

import numpy as np
import pytest

from align_genotype.qc_calibration import stats


def test_percentiles_of_a_known_range():
    values = np.arange(1, 101, dtype=float)
    result = dict(zip(stats.PERCENTILES, stats.percentiles(values), strict=True))
    assert result[1] == pytest.approx(1.99)
    assert result[25] == pytest.approx(25.75)
    assert result[50] == pytest.approx(50.5)
    assert result[99] == pytest.approx(99.01)


def test_percentiles_of_empty_is_empty():
    assert stats.percentiles(np.array([])) == ()


def test_percentiles_of_single_value():
    assert stats.percentiles(np.array([7.0])) == tuple(7.0 for _ in stats.PERCENTILES)


def test_breach_min_flags_below():
    np.testing.assert_array_equal(
        stats.breach(np.array([1.0, 5.0, 10.0]), 5.0, 'min'),
        np.array([True, False, False]),
    )


def test_breach_max_flags_above():
    np.testing.assert_array_equal(
        stats.breach(np.array([1.0, 5.0, 10.0]), 5.0, 'max'),
        np.array([False, False, True]),
    )


def test_breach_is_strict_at_the_threshold():
    """A value exactly on the threshold is not a breach, matching production's < and >."""
    np.testing.assert_array_equal(stats.breach(np.array([5.0]), 5.0, 'min'), np.array([False]))
    np.testing.assert_array_equal(stats.breach(np.array([5.0]), 5.0, 'max'), np.array([False]))


def test_flag_rates_warn_excludes_fail():
    """Production records one flag at the worst tier, so warn must not double-count fails."""
    values = np.array([1.0, 2.0, 3.0, 100.0])  # max metric: fail>50, warn>2
    fail_rate, warn_rate = stats.flag_rates(values, 'max', fail=50, warn=2)
    assert fail_rate == pytest.approx(0.25)  # just the 100
    assert warn_rate == pytest.approx(0.25)  # just the 3; the 100 is already a fail


def test_flag_rates_min_direction_warn_excludes_fail():
    values = np.array([1.0, 20.0, 30.0, 100.0])  # min metric: fail<15, warn<25
    fail_rate, warn_rate = stats.flag_rates(values, 'min', fail=15, warn=25)
    assert fail_rate == pytest.approx(0.25)  # just the 1.0
    assert warn_rate == pytest.approx(0.25)  # just the 20.0


def test_flag_rates_with_only_one_tier():
    values = np.array([1.0, 2.0, 3.0, 100.0])
    assert stats.flag_rates(values, 'max', fail=50, warn=None) == (pytest.approx(0.25), 0.0)
    assert stats.flag_rates(values, 'max', fail=None, warn=2) == (0.0, pytest.approx(0.5))


def test_flag_rates_with_no_tiers_is_zero():
    assert stats.flag_rates(np.array([1.0, 2.0]), 'max', fail=None, warn=None) == (0.0, 0.0)


def test_flag_rates_of_empty_is_nan():
    fail_rate, warn_rate = stats.flag_rates(np.array([]), 'max', fail=1, warn=2)
    assert np.isnan(fail_rate) and np.isnan(warn_rate)


def test_needs_review_uses_the_healthy_cohort_guardrail():
    assert stats.needs_review(0.0, 0.04) is False
    assert stats.needs_review(0.05, 0.04) is True  # too many fails
    assert stats.needs_review(0.0, 0.30) is True  # too many warns


def test_needs_review_is_false_for_nan():
    """An absent metric shouldn't be reported as needing review."""
    assert stats.needs_review(float('nan'), float('nan')) is False


# --- churn ---------------------------------------------------------------------

INITIAL = np.array([1.0, 2.0, 3.0, 4.0, 20.0])
GROWN = np.concatenate([INITIAL, np.array([18.0, 19.0, 20.0, 21.0, 22.0])])


def test_churn_measures_flips_caused_purely_by_cohort_growth():
    result = stats.churn(INITIAL, GROWN, 'max', 3.5)
    assert result.threshold_before == pytest.approx(8.189, abs=0.001)
    assert result.threshold_after == pytest.approx(34.067, abs=0.001)
    assert result.flagged_before == 1  # 20 is an outlier in the small cohort
    assert result.flagged_after == 0  # once the cohort includes similar values, it isn't
    assert result.flips == 1
    assert result.flip_rate == pytest.approx(0.2)


def test_churn_rounds_to_four_dp_like_production():
    result = stats.churn(INITIAL, GROWN, 'max', 3.5)
    assert result.threshold_before == round(result.threshold_before, 4)
    assert result.threshold_after == round(result.threshold_after, 4)


def test_churn_returns_none_on_degenerate_mad():
    identical = np.array([5.0] * 10)
    assert stats.churn(identical, identical, 'max', 3.5) is None


def test_churn_returns_none_when_either_cohort_is_degenerate():
    identical = np.array([5.0] * 10)
    varied = np.array([1.0, 2.0, 3.0, 4.0, 20.0])
    assert stats.churn(identical, varied, 'max', 3.5) is None
    assert stats.churn(varied, identical, 'max', 3.5) is None


def test_churn_of_empty_initial_is_none():
    assert stats.churn(np.array([]), GROWN, 'max', 3.5) is None


def test_churn_zero_flips_when_cohort_is_unchanged():
    result = stats.churn(INITIAL, INITIAL, 'max', 3.5)
    assert result.flips == 0
    assert result.flip_rate == 0.0
    assert result.threshold_before == result.threshold_after
