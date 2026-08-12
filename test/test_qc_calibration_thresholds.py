"""Unit tests for candidate fixed thresholds."""

import pytest

from align_genotype.qc_calibration import settings as settings_mod
from align_genotype.qc_calibration import thresholds as thresholds_mod
from align_genotype.qc_calibration import values as values_mod


def metric_values(values) -> values_mod.MetricValues:
    return values_mod.MetricValues(
        entries=tuple(('s', f'CPG{i}', float(v)) for i, v in enumerate(values)),
        n_dropped=0,
    )


def test_min_metric_seeds_from_the_lowest_per_dataset_tails():
    """Worst-case per dataset, so a candidate does not already flag the marginal one."""
    metric = settings_mod.MetricSpec(key='COV', direction='min', unit='x')
    by_dataset = {'ds-a': metric_values(range(30, 131)), 'ds-b': metric_values(range(10, 111))}
    candidate = thresholds_mod.candidate(by_dataset, metric)
    # ds-b is the marginal dataset: its p1 is ~11 and its p5 is ~15.
    assert candidate.fail == 11
    assert candidate.warn == 15


def test_max_metric_seeds_from_the_highest_per_dataset_tails():
    metric = settings_mod.MetricSpec(key='DUP', direction='max', unit='%')
    by_dataset = {'ds-a': metric_values(range(1, 102)), 'ds-b': metric_values(range(20, 121))}
    candidate = thresholds_mod.candidate(by_dataset, metric)
    assert candidate.fail == 119
    assert candidate.warn == 115


def test_x_and_percent_units_round_to_whole_numbers():
    metric = settings_mod.MetricSpec(key='COV', direction='min', unit='x')
    candidate = thresholds_mod.candidate({'ds': metric_values([10.4, 20.0, 30.0])}, metric)
    assert isinstance(candidate.fail, int)


def test_frac_unit_keeps_two_decimal_places():
    """A fraction rounded to an integer would collapse to 0 or 1."""
    metric = settings_mod.MetricSpec(key='PCT_20X', direction='min', unit='frac')
    candidate = thresholds_mod.candidate({'ds': metric_values([0.9012, 0.95, 0.97])}, metric)
    assert candidate.fail == pytest.approx(0.9, abs=0.01)
    assert isinstance(candidate.fail, float)


def test_a_relative_metric_gets_no_absolute_warn():
    """The relative tier is the warn tier; both would double-flag the same values."""
    metric = settings_mod.MetricSpec(key='DUP', direction='max', unit='%', relative=True)
    candidate = thresholds_mod.candidate({'ds': metric_values(range(1, 102))}, metric)
    assert candidate.fail is not None
    assert candidate.warn is None
    assert 'relative' in candidate.basis


def test_a_metric_with_no_values_anywhere_yields_no_candidate():
    metric = settings_mod.MetricSpec(key='ABSENT', direction='min', unit='x')
    assert thresholds_mod.candidate({'ds-a': metric_values([]), 'ds-b': metric_values([])}, metric) is None


def test_no_datasets_at_all_yields_no_candidate():
    metric = settings_mod.MetricSpec(key='COV', direction='min', unit='x')
    assert thresholds_mod.candidate({}, metric) is None


def test_datasets_without_values_are_skipped_not_counted():
    metric = settings_mod.MetricSpec(key='COV', direction='min', unit='x')
    by_dataset = {'ds-a': metric_values(range(30, 131)), 'ds-empty': metric_values([])}
    candidate = thresholds_mod.candidate(by_dataset, metric)
    assert '1 dataset' in candidate.basis


def test_basis_names_the_percentiles_used():
    metric = settings_mod.MetricSpec(key='COV', direction='min', unit='x')
    candidate = thresholds_mod.candidate({'ds': metric_values(range(30, 131))}, metric)
    assert 'p1' in candidate.basis
    assert 'p5' in candidate.basis


def test_basis_pluralises_dataset_count():
    metric = settings_mod.MetricSpec(key='COV', direction='min', unit='x')
    by_dataset = {'ds-a': metric_values(range(30, 131)), 'ds-b': metric_values(range(20, 121))}
    assert '2 datasets' in thresholds_mod.candidate(by_dataset, metric).basis


def test_candidate_reports_the_metric_key():
    metric = settings_mod.MetricSpec(key='COV', direction='min', unit='x')
    assert thresholds_mod.candidate({'ds': metric_values([1.0, 2.0, 3.0])}, metric).metric == 'COV'
