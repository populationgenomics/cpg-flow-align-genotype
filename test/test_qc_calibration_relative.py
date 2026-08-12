"""Unit tests for dataset-relative (MAD) tier evaluation."""

import numpy as np
import pytest

from align_genotype.qc_calibration import relative as relative_mod
from align_genotype.qc_calibration import settings as settings_mod
from align_genotype.qc_calibration import values as values_mod

METRIC = settings_mod.MetricSpec(key='dup_pct', direction='max', unit='%', relative=True)


def settings(min_samples=4, **bars: float) -> settings_mod.CalibrationSettings:
    return settings_mod.CalibrationSettings(
        seq_type='genome',
        metrics=(METRIC,),
        k=3.5,
        min_samples=min_samples,
        bars=settings_mod.Bars(**bars),
    )


def metric_values(values, section='samtools') -> values_mod.MetricValues:
    return values_mod.MetricValues(
        entries=tuple((section, f'CPG{i}', float(v)) for i, v in enumerate(values)),
        n_dropped=0,
    )


def test_derives_median_mad_and_threshold_per_dataset():
    by_dataset = {'ds-a': metric_values([10.0, 10.0, 12.0, 12.0, 40.0])}
    result = relative_mod.evaluate(by_dataset, METRIC, settings())
    (dataset,) = result.datasets
    assert dataset.dataset == 'ds-a'
    assert dataset.median == pytest.approx(12.0)
    assert dataset.mad_raw == pytest.approx(2.0)
    # median + k*MAD/0.6745 = 12 + 3.5*2/0.6745 = 22.3781...
    assert dataset.threshold == pytest.approx(22.3781, abs=1e-4)
    assert dataset.n_warn == 1  # only the 40.0
    assert dataset.skipped is None


def test_threshold_matches_production_rounded_to_four_dp():
    """The displayed threshold must be the number production compares against."""
    from align_genotype.scripts import check_multiqc  # noqa: PLC0415

    values = [10.0, 10.0, 12.0, 12.0, 40.0]
    expected = round(check_multiqc.robust_threshold(values, 'max', 3.5), 4)
    result = relative_mod.evaluate({'ds-a': metric_values(values)}, METRIC, settings())
    assert result.datasets[0].threshold == expected


def test_dataset_below_min_samples_is_skipped_with_a_size_reason():
    result = relative_mod.evaluate({'ds-a': metric_values([1.0, 2.0])}, METRIC, settings(min_samples=50))
    (dataset,) = result.datasets
    assert dataset.threshold is None
    assert dataset.n_warn == 0
    assert '2 values < min_samples 50' in dataset.skipped


def test_dataset_with_no_values_says_so_rather_than_too_small():
    """A collection problem and a size problem send an operator to different places."""
    result = relative_mod.evaluate({'ds-a': metric_values([])}, METRIC, settings())
    (dataset,) = result.datasets
    assert 'no values' in dataset.skipped
    # Distinguished from the below-min_samples branch, which names a size instead.
    assert 'min_samples' not in dataset.skipped
    # n_values == 0 must not raise ZeroDivisionError in the warn-rate property.
    assert dataset.warn_rate == 0.0


def test_zero_mad_dataset_is_skipped_as_degenerate():
    result = relative_mod.evaluate({'ds-a': metric_values([5.0] * 10)}, METRIC, settings())
    (dataset,) = result.datasets
    assert dataset.threshold is None
    assert 'zero MAD' in dataset.skipped


def test_warn_rate_is_per_value_and_duplication_is_flagged():
    by_dataset = {
        'ds-a': values_mod.MetricValues(
            entries=(
                ('picard_1', 'CPG0', 10.0),
                ('picard_4', 'CPG0', 10.0),
                ('picard_1', 'CPG1', 12.0),
                ('picard_4', 'CPG1', 12.0),
                ('picard_1', 'CPG2', 12.0),
                ('picard_4', 'CPG2', 12.0),
            ),
            n_dropped=0,
        ),
    }
    (dataset,) = relative_mod.evaluate(by_dataset, METRIC, settings()).datasets
    assert dataset.n_values == 6
    assert dataset.n_groups_with_values == 3
    assert dataset.duplicated is True


def test_duplicated_is_false_when_every_value_is_a_distinct_sequencing_group():
    """The mirror of the case above: no section carries two values per group."""
    by_dataset = {'ds-a': metric_values([10.0, 10.0, 12.0, 12.0, 40.0])}
    (dataset,) = relative_mod.evaluate(by_dataset, METRIC, settings()).datasets
    assert dataset.n_values == dataset.n_groups_with_values == 5
    assert dataset.duplicated is False


def test_warn_rate_doubles_when_the_flagged_value_is_duplicated_across_sections():
    """The per-value warn rate, not a per-group one - `duplicated` is the signal it's in play."""
    singleton = {
        'ds-a': values_mod.MetricValues(
            entries=(
                ('samtools', 'CPG0', 10.0),
                ('samtools', 'CPG1', 10.5),
                ('samtools', 'CPG2', 11.0),
                ('samtools', 'CPG3', 11.5),
                ('samtools', 'CPG4', 12.0),
                ('samtools', 'CPG5', 40.0),
            ),
            n_dropped=0,
        ),
    }
    doubled = {
        'ds-a': values_mod.MetricValues(
            entries=singleton['ds-a'].entries + tuple(
                ('picard', sg, v) for _, sg, v in singleton['ds-a'].entries
            ),
            n_dropped=0,
        ),
    }
    (single,) = relative_mod.evaluate(singleton, METRIC, settings()).datasets
    (double,) = relative_mod.evaluate(doubled, METRIC, settings()).datasets
    assert single.duplicated is False
    assert double.duplicated is True
    assert double.n_warn == 2 * single.n_warn
    assert double.warn_rate == pytest.approx(single.warn_rate)


def test_growth_churn_reports_both_orderings_and_uses_the_worse():
    rng = np.random.default_rng(1)
    # Leading 60% is tight, trailing 40% is high - the batch-ordering effect.
    ordered = [*list(rng.normal(10, 0.5, 30)), *list(rng.normal(20, 0.5, 20))]
    result = relative_mod.evaluate({'ds-a': metric_values(ordered)}, METRIC, settings(min_samples=10))
    (growth,) = result.growth
    assert growth.ordered is not None
    assert growth.shuffled is not None
    assert growth.flip_rate == max(growth.ordered.flip_rate, growth.shuffled.flip_rate)


def test_growth_is_not_simulated_when_the_before_slice_is_below_min_samples():
    """Production would have skipped a dataset that size, so churn against it is fiction."""
    result = relative_mod.evaluate({'ds-a': metric_values(range(1, 21))}, METRIC, settings(min_samples=20))
    assert result.growth == ()


def test_growth_is_reproducible_across_repeated_evaluations():
    """A verdict that changes on re-run, with no input change, is worthless."""
    values = [9.84, 8.96, 8.33, 9.51, 9.95, 11.77, 10.13, 10.98, 9.50, 8.82, 28.07, 28.55, 34.26]
    by_dataset = {'ds-a': metric_values(values)}
    first = relative_mod.evaluate(by_dataset, METRIC, settings(min_samples=5))
    second = relative_mod.evaluate(by_dataset, METRIC, settings(min_samples=5))
    assert first.growth == second.growth
    assert first.max_growth_churn == second.max_growth_churn


def test_growth_excludes_a_dataset_whose_whole_array_is_also_degenerate():
    """Both before and after thresholds are None on an all-identical dataset; no entry, not a zero."""
    result = relative_mod.evaluate({'ds-flat': metric_values([5.0] * 10)}, METRIC, settings(min_samples=4))
    assert result.growth == ()


def test_ordering_sensitivity_is_flagged_when_the_two_slices_disagree():
    """A metric where growth churn depends on which 60% you take needs a call-out."""
    values = [9.84, 8.96, 8.33, 9.51, 9.95, 11.77, 10.13, 10.98, 9.50, 8.82, 28.07, 28.55, 34.26]
    result = relative_mod.evaluate({'ds-a': metric_values(values)}, METRIC, settings(min_samples=5))
    (growth,) = result.growth
    assert growth.ordered.flip_rate > 0
    assert growth.shuffled.flip_rate == 0
    assert growth.ordering_sensitive(0.02) is True
    assert growth.ordering_sensitive(0.5) is False
    assert result.ordering_sensitive == ('ds-a',)


def test_merge_churn_simulates_both_directions_of_every_pair():
    by_dataset = {
        'ds-a': metric_values([10.0, 10.5, 11.0, 11.5, 12.0, 40.0]),
        'ds-b': metric_values([30.0, 30.5, 31.0, 31.5, 32.0, 60.0]),
    }
    result = relative_mod.evaluate(by_dataset, METRIC, settings())
    assert {(a, b) for a, b, _ in result.merge} == {('ds-a', 'ds-b'), ('ds-b', 'ds-a')}


def test_merge_churn_never_self_pairs():
    by_dataset = {'ds-a': metric_values([10.0, 11.0, 12.0, 13.0, 40.0])}
    assert relative_mod.evaluate(by_dataset, METRIC, settings()).merge == ()


def test_merge_churn_omits_a_direction_whose_before_threshold_is_degenerate():
    """A flat dataset has no threshold to merge *from*; only the other direction is measurable."""
    by_dataset = {
        'ds-flat': metric_values([5.0] * 6),
        'ds-normal': metric_values([10.0, 10.5, 11.0, 11.5, 12.0, 40.0]),
    }
    result = relative_mod.evaluate(by_dataset, METRIC, settings())
    pairs = {(a, b) for a, b, _ in result.merge}
    assert pairs == {('ds-normal', 'ds-flat')}


def test_verdict_recommends_when_every_bar_is_cleared():
    by_dataset = {'ds-a': metric_values([10.0, 10.5, 11.0, 11.5, 12.0, 11.2, 10.8, 11.4])}
    result = relative_mod.evaluate(by_dataset, METRIC, settings())
    assert result.verdict == 'RECOMMEND'
    assert result.verdict_reason == ''


def test_verdict_names_the_bar_that_was_missed():
    by_dataset = {'ds-a': metric_values([10.0, 10.5, 11.0, 11.5, 12.0, 40.0])}
    result = relative_mod.evaluate(by_dataset, METRIC, settings(max_warn_rate=0.0))
    assert result.verdict == 'REJECT'
    assert 'warn rate' in result.verdict_reason


def test_verdict_names_growth_churn_when_only_that_bar_is_missed():
    """The growth-specific reason string, isolated from warn-rate and merge-churn wording."""
    values = [9.84, 8.96, 8.33, 9.51, 9.95, 11.77, 10.13, 10.98, 9.50, 8.82, 28.07, 28.55, 34.26]
    result = relative_mod.evaluate(
        {'ds-a': metric_values(values)},
        METRIC,
        settings(min_samples=5, max_warn_rate=1.0),
    )
    assert result.verdict == 'REJECT'
    assert 'growth churn' in result.verdict_reason
    assert 'merge churn' not in result.verdict_reason
    assert 'warn rate' not in result.verdict_reason


def test_verdict_distinguishes_growth_from_merge_churn():
    """Conflating them would send an operator to the wrong table."""
    by_dataset = {
        'ds-a': metric_values([10.0, 10.5, 11.0, 11.5, 12.0, 40.0]),
        'ds-b': metric_values([30.0, 30.5, 31.0, 31.5, 32.0, 60.0]),
    }
    result = relative_mod.evaluate(by_dataset, METRIC, settings(max_merge_churn=0.0))
    assert 'merge churn' in result.verdict_reason
    assert 'growth churn' not in result.verdict_reason


def test_verdict_joins_every_missed_bar_when_more_than_one_is_missed():
    by_dataset = {
        'ds-a': metric_values([10.0, 10.5, 11.0, 11.5, 12.0, 40.0]),
        'ds-b': metric_values([30.0, 30.5, 31.0, 31.5, 32.0, 60.0]),
    }
    result = relative_mod.evaluate(by_dataset, METRIC, settings(max_warn_rate=0.0, max_merge_churn=0.0))
    assert 'warn rate' in result.verdict_reason
    assert 'merge churn' in result.verdict_reason
    assert '; ' in result.verdict_reason


def test_evaluate_refuses_a_metric_with_no_relative_tier():
    plain = settings_mod.MetricSpec(key='dup_pct', direction='max', unit='%', relative=False)
    with pytest.raises(ValueError, match='no configured relative tier'):
        relative_mod.evaluate({'ds-a': metric_values([1.0, 2.0, 3.0])}, plain, settings())


def test_max_warn_rate_ignores_skipped_datasets():
    by_dataset = {
        'ds-a': metric_values([10.0, 10.5, 11.0, 11.5, 12.0, 40.0]),
        'ds-tiny': metric_values([1.0]),
    }
    result = relative_mod.evaluate(by_dataset, METRIC, settings())
    skipped = next(d for d in result.datasets if d.dataset == 'ds-tiny')
    assert skipped.skipped is not None
    assert result.max_warn_rate == max(d.warn_rate for d in result.datasets if d.skipped is None)


def test_global_config_detour_is_not_reintroduced():
    """The hazard this rewrite removes: mutating cpg_utils.config's installed paths.

    This is a source-substring guard, not a semantic one: it catches the detour coming
    back the way it would plausibly come back - someone re-adding a call or import who
    doesn't know why it was removed. It does not catch a determined reintroduction via a
    fully-qualified import (``import align_genotype.qc_calibration.tomlio`` contains no
    literal ``import tomlio``) or dynamic dispatch (``getattr(check_multiqc, 'relative' +
    '_flags')``). Treat a pass here as "nobody did this by accident", not as a hard
    barrier.

    The checks look for call/import sites rather than bare substrings for a narrower
    reason: this module's own docstring names `_relative_flags_for_metric` (production's
    function) as the thing whose rounding it mirrors, and that name contains
    `relative_flags` as a substring - a bare substring check would fail on the docstring
    itself.
    """
    import inspect  # noqa: PLC0415

    source = inspect.getsource(relative_mod)
    for banned in ('set_config_paths', 'import tomlio', 'from align_genotype.qc_calibration import tomlio'):
        assert banned not in source
    assert 'check_multiqc.relative_flags(' not in source
    assert '_production_config' not in source
    assert '_restore_config_paths' not in source
    assert '_warn_count' not in source
