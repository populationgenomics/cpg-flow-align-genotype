"""Unit tests for cohort-relative (MAD) evaluation."""

import os
from collections.abc import Iterator
from typing import Any, NoReturn

import numpy as np
import pytest

from cpg_utils import config

from align_genotype.qc_calibration import relative as relative_mod
from align_genotype.qc_calibration import spec as spec_mod
from align_genotype.qc_calibration.cache import CohortValues, ValueCache
from align_genotype.qc_calibration.relative import CohortMad, HomogeneousChurn, MadEvaluation
from align_genotype.qc_calibration.stats import ChurnResult


def _metric(k: float = 3.5, min_cohort: int = 5) -> spec_mod.MetricSpec:
    return spec_mod.loads(
        'seq_type = "genome"\ncache = "c.json"\n'
        '[metrics.DUP]\ndirection = "max"\nunit = "%"\nfail = 1000.0\n'
        f'[metrics.DUP.relative]\nk = {k}\nmin_cohort = {min_cohort}\n',
    ).metric('DUP')


METRIC = _metric()

# 49 evenly spread values plus one extreme: median 25.5, raw MAD 12.5,
# so the modified-z threshold lands at 25.5 + 3.5*12.5/0.6745 = 90.3629.
OUTLIER_COHORT = [*range(1, 50), 500.0]

# A cohort whose churn depends on *which* 60% is taken as the before-slice. Ordered by
# value, so the leading 30 are the tight 10.0/10.5 body plus the five 20.0s: raw MAD 0.5,
# threshold 13.0945, and those five breach it. A seeded-shuffled 30 mixes in the 20.5
# tail, and the full cohort's threshold is 42.4924, which none of them breach.
ORDERING_SENSITIVE_COHORT = [10.0] * 13 + [10.5] * 12 + [20.0] * 5 + [20.5] * 20


@pytest.fixture(autouse=True)
def _restore_config_paths() -> Iterator[None]:
    """Keep the global cpg-utils config state from leaking between tests."""
    previous = os.environ.get('CPG_CONFIG_PATH', '')
    yield
    config.set_config_paths([p for p in previous.split(',') if p])


def _cache(**cohorts: list[float]) -> ValueCache:
    """A cache where every cohort carries one value per sample (no duplication)."""
    return ValueCache(
        seq_type='genome',
        generated='x',
        complete=True,
        metrics=('DUP',),
        cohorts=tuple(
            CohortValues(label, len(values), '1.33', 'dict', 0, {'DUP': [float(v) for v in values]})
            for label, values in cohorts.items()
        ),
    )


def _cache_with_sample_count(label: str, values: list[float], n_samples: int) -> ValueCache:
    """A cache whose value count and sample count deliberately disagree."""
    return ValueCache(
        seq_type='genome',
        generated='x',
        complete=True,
        metrics=('DUP',),
        cohorts=(CohortValues(label, n_samples, '1.33', 'dict', 0, {'DUP': [float(v) for v in values]}),),
    )


# --- per-cohort numbers -----------------------------------------------------------


def test_evaluate_reports_median_mad_and_threshold():
    evaluation = relative_mod.evaluate(_cache(**{'dataset-a': OUTLIER_COHORT}), METRIC, 'genome')
    cohort = evaluation.cohorts[0]
    assert cohort.label == 'dataset-a'
    assert cohort.n_values == 50
    assert cohort.median == pytest.approx(25.5)
    assert cohort.mad_raw == pytest.approx(12.5)
    assert cohort.threshold == pytest.approx(90.363, abs=0.001)


def test_warn_count_comes_from_the_production_relative_flags_path():
    evaluation = relative_mod.evaluate(_cache(**{'dataset-a': OUTLIER_COHORT}), METRIC, 'genome')
    assert evaluation.cohorts[0].n_warn == 1
    assert evaluation.cohorts[0].warn_rate == pytest.approx(0.02)


def test_cohort_below_min_cohort_is_skipped():
    evaluation = relative_mod.evaluate(_cache(**{'dataset-a': [1.0, 2.0, 3.0]}), METRIC, 'genome')
    cohort = evaluation.cohorts[0]
    assert cohort.skipped is not None
    assert 'min_cohort' in cohort.skipped
    assert cohort.n_warn == 0
    assert cohort.threshold is None


def test_degenerate_mad_cohort_is_skipped_and_points_at_the_absolute_gate():
    evaluation = relative_mod.evaluate(_cache(**{'dataset-a': [5.0] * 10}), METRIC, 'genome')
    cohort = evaluation.cohorts[0]
    assert cohort.threshold is None
    assert 'zero MAD' in cohort.skipped
    # The operator needs to be told what to do instead, not just what failed.
    assert 'absolute gate' in cohort.skipped
    assert cohort.n_warn == 0


def test_a_metric_absent_from_the_cohort_is_reported_as_absent_not_as_too_small():
    """`cohort 0 < min_cohort 5` would send an operator to the wrong problem."""
    evaluation = relative_mod.evaluate(_cache(**{'dataset-a': []}), METRIC, 'genome')
    cohort = evaluation.cohorts[0]
    assert cohort.n_values == 0
    assert 'no values' in cohort.skipped
    assert 'min_cohort' not in cohort.skipped
    assert cohort.threshold is None
    assert np.isnan(cohort.median)
    assert np.isnan(cohort.mad_raw)


def test_value_count_and_sample_count_are_reported_separately():
    """`n_values > n_samples` is the visible signal that a metric spans two sections."""
    cache = _cache_with_sample_count('dataset-a', [float(v) for v in OUTLIER_COHORT], n_samples=26)
    cohort = relative_mod.evaluate(cache, METRIC, 'genome').cohorts[0]
    assert cohort.n_values == 50
    assert cohort.n_samples == 26
    assert cohort.duplicated is True
    # The rate stays per-value, consistent with the threshold reported beside it.
    assert cohort.warn_rate == pytest.approx(1 / 50)


def test_a_cohort_with_one_value_per_sample_is_not_reported_as_duplicated():
    cohort = relative_mod.evaluate(_cache(**{'dataset-a': OUTLIER_COHORT}), METRIC, 'genome').cohorts[0]
    assert cohort.n_values == cohort.n_samples == 50
    assert cohort.duplicated is False


# --- config plumbing --------------------------------------------------------------


def test_evaluate_restores_previous_config_paths(tmp_path):
    existing = tmp_path / 'existing.toml'
    existing.write_text('[workflow]\nsequencing_type = "genome"\n')
    config.set_config_paths([str(existing)])
    relative_mod.evaluate(_cache(**{'dataset-a': OUTLIER_COHORT}), METRIC, 'genome')
    assert config.get_config_paths() == [str(existing)]


def test_evaluate_restores_config_paths_even_on_failure(tmp_path, monkeypatch):
    """The context manager must not leak global config state when the body raises."""
    existing = tmp_path / 'existing.toml'
    existing.write_text('[workflow]\nsequencing_type = "genome"\n')
    config.set_config_paths([str(existing)])
    monkeypatch.setattr(relative_mod.check_multiqc, 'relative_flags', lambda *_a, **_k: 1 / 0)
    with pytest.raises(ZeroDivisionError):
        relative_mod.evaluate(_cache(**{'dataset-a': OUTLIER_COHORT}), METRIC, 'genome')
    assert config.get_config_paths() == [str(existing)]


def test_an_unreadable_previous_config_does_not_mask_the_bodys_exception(tmp_path, monkeypatch):
    """Restoring re-reads every previous path; that failure must not win the traceback.

    Under analysis-runner the previous paths are `gs://` URLs, so a transient read
    failure here is real - and it would otherwise surface as a `ValueError` about
    missing config files, hiding what actually went wrong.
    """
    existing = tmp_path / 'existing.toml'
    existing.write_text('[workflow]\nsequencing_type = "genome"\n')
    config.set_config_paths([str(existing)])

    def _unlink_then_fail(*_a: Any, **_k: Any) -> NoReturn:
        existing.unlink()
        raise RuntimeError('the real problem')

    monkeypatch.setattr(relative_mod.check_multiqc, 'relative_flags', _unlink_then_fail)
    with pytest.raises(RuntimeError, match='the real problem'):
        relative_mod.evaluate(_cache(**{'dataset-a': OUTLIER_COHORT}), METRIC, 'genome')


def test_an_unreadable_previous_config_does_not_fail_a_successful_evaluation(tmp_path, monkeypatch):
    existing = tmp_path / 'existing.toml'
    existing.write_text('[workflow]\nsequencing_type = "genome"\n')
    config.set_config_paths([str(existing)])
    real_relative_flags = relative_mod.check_multiqc.relative_flags

    def _unlink_then_delegate(*args: Any, **kwargs: Any) -> list:
        existing.unlink(missing_ok=True)
        return real_relative_flags(*args, **kwargs)

    monkeypatch.setattr(relative_mod.check_multiqc, 'relative_flags', _unlink_then_delegate)
    evaluation = relative_mod.evaluate(_cache(**{'dataset-a': OUTLIER_COHORT}), METRIC, 'genome')
    assert evaluation.cohorts[0].n_warn == 1


def test_metric_without_a_relative_block_raises():
    plain = spec_mod.loads(
        'seq_type = "genome"\ncache = "c.json"\n[metrics.M]\ndirection = "max"\nfail = 1\n',
    ).metric('M')
    with pytest.raises(ValueError, match=r'no .* relative'):
        relative_mod.evaluate(_cache(**{'dataset-a': OUTLIER_COHORT}), plain, 'genome')


def test_warn_count_tracks_the_k_written_into_the_production_config():
    """Proof the warn count really comes from the temp config, not a silent no-op.

    With k=3.5 only the 500 outlier clears the line; a k that small pulls the
    threshold down below most of the cohort. If `relative_flags` were reading an
    empty spec, both counts would be 0.
    """
    loose = relative_mod.evaluate(_cache(**{'dataset-a': OUTLIER_COHORT}), METRIC, 'genome')
    tight = relative_mod.evaluate(_cache(**{'dataset-a': OUTLIER_COHORT}), _metric(k=0.1), 'genome')
    assert loose.cohorts[0].n_warn == 1
    assert tight.cohorts[0].n_warn > 1


def test_evaluate_works_from_a_cold_start_and_leaves_the_env_key_absent(monkeypatch):
    """`get_config_paths()` raises when nothing was ever set - the restore must cope."""
    monkeypatch.delenv('CPG_CONFIG_PATH', raising=False)
    config.set_config_paths([])
    evaluation = relative_mod.evaluate(_cache(**{'dataset-a': OUTLIER_COHORT}), METRIC, 'genome')
    assert evaluation.cohorts[0].n_warn == 1
    # Absent, not present-but-empty: a subprocess testing for the key must see what it
    # would have seen had `evaluate` never run.
    assert 'CPG_CONFIG_PATH' not in os.environ
    with pytest.raises(config.ConfigError):
        config.get_config_paths()


# --- churn --------------------------------------------------------------------


def test_evaluate_simulates_homogeneous_and_heterogeneous_growth():
    evaluation = relative_mod.evaluate(
        _cache(**{'dataset-a': OUTLIER_COHORT, 'dataset-b': [float(v) for v in range(200, 250)]}),
        METRIC,
        'genome',
    )
    assert {h.label for h in evaluation.homogeneous} == {'dataset-a', 'dataset-b'}
    assert {(a, b) for a, b, _ in evaluation.heterogeneous} == {
        ('dataset-a', 'dataset-b'),
        ('dataset-b', 'dataset-a'),
    }


def test_cohorts_too_small_are_excluded_from_churn():
    evaluation = relative_mod.evaluate(_cache(**{'dataset-a': [1.0, 2.0, 3.0]}), METRIC, 'genome')
    assert evaluation.homogeneous == ()
    assert evaluation.heterogeneous == ()


def test_both_merge_directions_are_simulated_and_no_cohort_is_self_paired():
    """A merge churns both projects' flag sets, so `n*(n-1)` entries, not `n*(n-1)/2`.

    One direction per pair would leave the headline merge figure depending on cache
    insertion order.
    """
    labels = ('dataset-a', 'dataset-b', 'dataset-c')
    evaluation = relative_mod.evaluate(
        _cache(
            **{
                'dataset-a': OUTLIER_COHORT,
                'dataset-b': [float(v) for v in range(200, 250)],
                'dataset-c': [float(v) for v in range(300, 350)],
            }
        ),
        METRIC,
        'genome',
    )
    pairs = {(a, b) for a, b, _ in evaluation.heterogeneous}
    assert pairs == {(a, b) for a in labels for b in labels if a != b}
    assert len(evaluation.heterogeneous) == len(labels) * (len(labels) - 1)


def test_the_two_merge_directions_of_one_pair_measure_different_churn():
    """Which cohort's flag set is being re-scored changes the answer.

    `dataset-a` is a tight 1..49 body plus one extreme; `dataset-b` sits far above it.
    Pooling them disturbs the two flag sets by different amounts, which is why simulating
    only one direction under-reports - measured at 59.09% vs 64.71% on the real WGS set.
    """
    evaluation = relative_mod.evaluate(
        _cache(
            **{
                'dataset-a': OUTLIER_COHORT,
                'dataset-b': [float(v) for v in range(200, 250)],
            }
        ),
        METRIC,
        'genome',
    )
    by_pair = {(a, b): result.flip_rate for a, b, result in evaluation.heterogeneous}
    assert by_pair[('dataset-a', 'dataset-b')] != by_pair[('dataset-b', 'dataset-a')]
    assert evaluation.max_merge_churn == pytest.approx(max(by_pair.values()))


def test_both_before_slice_orderings_are_simulated_and_the_worst_one_counts():
    """Churn must not swing on MultiQC's JSON key order alone.

    Exact rates pin the seeded shuffle: an unseeded or differently-seeded RNG would not
    reproduce 13/30. (numpy guarantees `default_rng` stream stability across versions.)
    """
    evaluation = relative_mod.evaluate(
        _cache(**{'dataset-a': ORDERING_SENSITIVE_COHORT}),
        METRIC,
        'genome',
    )
    growth = evaluation.homogeneous[0]
    assert growth.label == 'dataset-a'
    assert growth.ordered.flip_rate == pytest.approx(5 / 30)
    assert growth.shuffled.flip_rate == pytest.approx(13 / 30)
    assert growth.flip_rate == pytest.approx(13 / 30)
    assert evaluation.max_growth_churn == pytest.approx(13 / 30)


def test_a_before_slice_below_min_cohort_is_not_simulated():
    """Production skips a cohort that size outright, so its churn cannot occur.

    With the shipped `min_cohort = 50`, a 60-sample cohort's 60% slice is 36 - a state
    production never reaches. A 100-sample cohort's slice is 60, which it does.
    """
    metric = _metric(min_cohort=50)
    small = relative_mod.evaluate(_cache(**{'dataset-a': [float(v) for v in range(60)]}), metric, 'genome')
    assert small.cohorts[0].skipped is None  # the cohort itself is big enough
    assert small.cohorts[0].threshold is not None
    assert small.homogeneous == ()

    large = relative_mod.evaluate(_cache(**{'dataset-a': [float(v) for v in range(100)]}), metric, 'genome')
    assert {h.label for h in large.homogeneous} == {'dataset-a'}


# --- HomogeneousChurn -------------------------------------------------------------


def _churn(flip_rate: float) -> ChurnResult:
    return ChurnResult(1.0, 1.0, 10000, 0, 0, round(flip_rate * 10000))


def test_homogeneous_flip_rate_is_the_worst_of_the_two_orderings():
    growth = HomogeneousChurn('dataset-a', ordered=_churn(0.01), shuffled=_churn(0.10))
    assert growth.flip_rate == pytest.approx(0.10)
    assert HomogeneousChurn('dataset-a', ordered=_churn(0.10), shuffled=_churn(0.01)).flip_rate == pytest.approx(0.10)


def test_homogeneous_flip_rate_tolerates_a_degenerate_ordering():
    assert HomogeneousChurn('dataset-a', ordered=None, shuffled=_churn(0.05)).flip_rate == pytest.approx(0.05)
    assert HomogeneousChurn('dataset-a', ordered=_churn(0.05), shuffled=None).flip_rate == pytest.approx(0.05)
    assert HomogeneousChurn('dataset-a', ordered=None, shuffled=None).flip_rate == 0.0


def test_ordering_sensitive_only_when_the_orderings_straddle_the_bar():
    assert HomogeneousChurn('dataset-a', _churn(0.01), _churn(0.10)).ordering_sensitive is True
    assert HomogeneousChurn('dataset-a', _churn(0.01), _churn(0.015)).ordering_sensitive is False
    assert HomogeneousChurn('dataset-a', _churn(0.10), _churn(0.245)).ordering_sensitive is False
    # One ordering degenerate: nothing to disagree with.
    assert HomogeneousChurn('dataset-a', None, _churn(0.10)).ordering_sensitive is False


# --- verdict ------------------------------------------------------------------


def _evaluation(
    warn_rates: list[float],
    flip_rates: list[float],
    merge_rates: tuple[float, ...] = (),
) -> MadEvaluation:
    cohorts = tuple(
        CohortMad(
            f'cohort-{i}',
            n_values=10000,
            n_samples=10000,
            median=1.0,
            mad_raw=1.0,
            threshold=5.0,
            n_warn=round(rate * 10000),
            skipped=None,
        )
        for i, rate in enumerate(warn_rates)
    )
    churns = tuple(
        HomogeneousChurn(f'cohort-{i}', ordered=_churn(rate), shuffled=None) for i, rate in enumerate(flip_rates)
    )
    merges = tuple((f'cohort-{i}', f'cohort-{i + 1}', _churn(rate)) for i, rate in enumerate(merge_rates))
    return MadEvaluation('DUP', 'max', cohorts, churns, merges)


def test_verdict_recommends_when_warn_and_churn_are_low():
    evaluation = _evaluation(warn_rates=[0.0, 0.042], flip_rates=[0.01, 0.02])
    assert evaluation.max_warn_rate == pytest.approx(0.042)
    assert evaluation.max_churn == pytest.approx(0.02)
    assert evaluation.verdict == 'RECOMMEND'
    assert evaluation.verdict_reason == ''


def test_verdict_rejects_on_high_churn():
    evaluation = _evaluation(warn_rates=[0.02], flip_rates=[0.245])
    assert evaluation.verdict == 'REJECT'
    assert 'churn' in evaluation.verdict_reason


def test_the_shipped_genome_duplication_tier_is_rejected_on_the_full_cohort_set():
    """The shipped tier does not clear these bars, and the tool must say so.

    Genome `reads_duplicated_percent` is in production, adopted on a hand-picked subset -
    three cohorts for growth, two ordered pairs for merge - which measured ~1% growth and
    ~2.1% merge. This module reproduces those numbers exactly on that same subset. Across
    all 10 WGS cohorts and every ordered pair it measures 8.61% growth and 64.71% merge,
    so both bars are breached.

    Pinned deliberately, and not as a RECOMMEND: shipping the tool honest is the decision,
    and the gap between this verdict and the shipped tier is a live QC question about the
    tier rather than a sign these constants need loosening.
    """
    evaluation = _evaluation(warn_rates=[0.03], flip_rates=[0.0861], merge_rates=(0.6471,))
    assert evaluation.max_growth_churn == pytest.approx(0.0861)
    assert evaluation.max_merge_churn == pytest.approx(0.6471)
    assert evaluation.verdict == 'REJECT'
    assert 'growth churn' in evaluation.verdict_reason
    assert 'merge churn' in evaluation.verdict_reason


def test_the_shipped_exome_zero_coverage_tier_is_rejected_on_merge_alone():
    """Exome ZERO_CVG_TARGETS_PCT: 1.0% growth clears, 51.11% merge does not."""
    evaluation = _evaluation(warn_rates=[0.03], flip_rates=[0.010], merge_rates=(0.5111,))
    assert evaluation.verdict == 'REJECT'
    assert 'growth' not in evaluation.verdict_reason
    assert 'merge churn' in evaluation.verdict_reason


def test_verdict_rejects_just_above_the_growth_bar_and_names_growth():
    """The growth bar is strict: 0.02 clears, 0.0201 does not."""
    assert _evaluation(warn_rates=[0.0], flip_rates=[0.02]).verdict == 'RECOMMEND'
    over = _evaluation(warn_rates=[0.0], flip_rates=[0.0201])
    assert over.verdict == 'REJECT'
    assert 'growth churn' in over.verdict_reason
    assert 'merge' not in over.verdict_reason


def test_verdict_rejects_well_over_the_merge_bar_and_names_merge():
    """24.5% is what rejected the exome PCT_SELECTED_BASES / PCT_OFF_BAIT candidates.

    Merge churn is looser than growth churn, not advisory - a metric this unstable must
    still be refused.
    """
    over = _evaluation(warn_rates=[0.0], flip_rates=[0.0], merge_rates=(0.245,))
    assert over.verdict == 'REJECT'
    assert 'merge churn' in over.verdict_reason
    assert 'growth' not in over.verdict_reason


def test_verdict_is_strict_at_the_merge_bar():
    assert _evaluation(warn_rates=[0.0], flip_rates=[0.0], merge_rates=(0.05,)).verdict == 'RECOMMEND'
    assert _evaluation(warn_rates=[0.0], flip_rates=[0.0], merge_rates=(0.0501,)).verdict == 'REJECT'


def test_verdict_reason_names_each_failing_simulation_separately():
    evaluation = _evaluation(warn_rates=[0.4], flip_rates=[0.0201], merge_rates=(0.245,))
    assert evaluation.verdict == 'REJECT'
    assert 'warn rate' in evaluation.verdict_reason
    assert 'growth churn' in evaluation.verdict_reason
    assert 'merge churn' in evaluation.verdict_reason


def test_max_churn_is_the_headline_across_both_simulations():
    """Kept as a display figure; the verdict judges the two separately."""
    evaluation = _evaluation(warn_rates=[0.0], flip_rates=[0.01], merge_rates=(0.04,))
    assert evaluation.max_churn == pytest.approx(0.04)
    assert evaluation.verdict == 'RECOMMEND'


def test_verdict_rejects_just_above_the_warn_bar():
    assert _evaluation(warn_rates=[0.10], flip_rates=[0.0]).verdict == 'RECOMMEND'
    assert _evaluation(warn_rates=[0.1001], flip_rates=[0.0]).verdict == 'REJECT'


def test_verdict_rejects_on_high_warn_rate():
    evaluation = _evaluation(warn_rates=[0.4], flip_rates=[0.0])
    assert evaluation.verdict == 'REJECT'
    assert 'warn' in evaluation.verdict_reason


def test_verdict_reports_both_reasons_when_warn_and_churn_both_fail():
    evaluation = _evaluation(warn_rates=[0.4], flip_rates=[0.245])
    assert evaluation.verdict == 'REJECT'
    assert 'warn' in evaluation.verdict_reason
    assert 'churn' in evaluation.verdict_reason


def test_max_warn_rate_ignores_skipped_cohorts():
    """Documents the invariant the skip filter rests on, rather than the filter itself.

    A skipped cohort has `n_warn == 0` and so `warn_rate == 0.0`; adding zeros to `max()`
    cannot change the result, which means this cannot be written as a test that fails
    when the filter is removed. What it does pin is the invariant - if a future change
    ever recorded a non-zero count for a cohort production skips, the filter starts
    mattering and this test documents why it is there.
    """
    skipped = CohortMad(
        'dataset-b',
        n_values=3,
        n_samples=3,
        median=1.0,
        mad_raw=0.0,
        threshold=None,
        n_warn=0,
        skipped='cohort 3 < min_cohort 5',
    )
    warned = CohortMad(
        'dataset-a',
        n_values=100,
        n_samples=100,
        median=1.0,
        mad_raw=1.0,
        threshold=5.0,
        n_warn=40,
        skipped=None,
    )
    assert skipped.warn_rate == 0.0  # the invariant
    assert MadEvaluation('DUP', 'max', (warned, skipped), (), ()).max_warn_rate == pytest.approx(0.4)


def test_max_churn_of_no_churn_data_is_zero():
    assert MadEvaluation('DUP', 'max', (), (), ()).max_churn == 0.0


def test_evaluation_lists_the_ordering_sensitive_cohorts():
    evaluation = MadEvaluation(
        'DUP',
        'max',
        (),
        (
            HomogeneousChurn('dataset-a', _churn(0.01), _churn(0.10)),
            HomogeneousChurn('dataset-b', _churn(0.01), _churn(0.015)),
        ),
        (),
    )
    assert evaluation.ordering_sensitive == ('dataset-a',)
