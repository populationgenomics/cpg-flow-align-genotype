"""Unit tests for the calibration value cache."""

import json

import numpy as np
import pytest

from align_genotype.qc_calibration import cache as cache_mod
from align_genotype.qc_calibration import spec as spec_mod
from align_genotype.qc_calibration.cache import CacheError, CohortValues, ValueCache

SPEC = spec_mod.loads(
    'seq_type = "genome"\ncache = "c.json"\n'
    '\n[metrics.MEDIAN_COVERAGE]\ndirection = "min"\nunit = "x"\nfail = 15\n'
    '\n[metrics.FREEMIX]\ndirection = "max"\nfail = 0.04\n',
)


def _cache(complete: bool = True, metrics: tuple[str, ...] = ('MEDIAN_COVERAGE', 'FREEMIX')) -> ValueCache:
    return ValueCache(
        seq_type='genome',
        generated='2026-08-11T13:52:00',
        complete=complete,
        metrics=metrics,
        cohorts=(
            CohortValues(
                label='dataset-a',
                n_samples=3,
                multiqc_version='1.33',
                shape='dict',
                n_dropped=1,
                values={'MEDIAN_COVERAGE': [30.0, 28.0, 12.0], 'FREEMIX': [0.001, 0.002, 0.05]},
            ),
        ),
    )


def test_save_and_load_round_trip(tmp_path):
    path = tmp_path / 'values.json'
    original = _cache()
    cache_mod.save(original, path)
    assert cache_mod.load(path) == original


def test_saved_json_is_readable_and_shaped_as_documented(tmp_path):
    path = tmp_path / 'values.json'
    cache_mod.save(_cache(), path)
    raw = json.loads(path.read_text())
    assert raw['seq_type'] == 'genome'
    assert raw['complete'] is True
    assert raw['metrics'] == ['MEDIAN_COVERAGE', 'FREEMIX']
    assert raw['cohorts']['dataset-a']['n_samples'] == 3
    assert raw['cohorts']['dataset-a']['values']['MEDIAN_COVERAGE'] == [30.0, 28.0, 12.0]


def test_save_refuses_to_write_nan(tmp_path):
    """NaN is filtered during collect; if one reaches here the cache must not go silently invalid."""
    bad = ValueCache(
        seq_type='genome',
        generated='x',
        complete=True,
        metrics=('M',),
        cohorts=(CohortValues('dataset-a', 1, '1.33', 'dict', 0, {'M': [float('nan')]}),),
    )
    with pytest.raises(ValueError, match=r'Out of range float values|NaN'):
        cache_mod.save(bad, tmp_path / 'values.json')


def test_save_handles_numpy_floats(tmp_path):
    """collect casts to float, but the sink shouldn't corrupt if one slips through."""
    c = ValueCache(
        seq_type='genome',
        generated='x',
        complete=True,
        metrics=('M',),
        cohorts=(CohortValues('dataset-a', 2, '1.33', 'dict', 0, {'M': [np.float64(1.5), np.float64(2.5)]}),),
    )
    path = tmp_path / 'values.json'
    cache_mod.save(c, path)
    assert json.loads(path.read_text())['cohorts']['dataset-a']['values']['M'] == [1.5, 2.5]


def test_labels_and_cohort_lookup():
    c = _cache()
    assert c.labels == ('dataset-a',)
    assert c.cohort('dataset-a').n_samples == 3
    with pytest.raises(KeyError, match='dataset-z'):
        c.cohort('dataset-z')


def test_series_returns_float_array():
    series = _cache().series('dataset-a', 'MEDIAN_COVERAGE')
    assert series.dtype == np.float64
    np.testing.assert_array_equal(series, np.array([30.0, 28.0, 12.0]))


def test_series_filters_non_finite_defensively():
    c = ValueCache(
        'genome',
        'x',
        True,
        ('M',),
        (CohortValues('dataset-a', 4, '1.33', 'dict', 0, {'M': [1.0, None, float('inf'), float('nan')]}),),
    )
    np.testing.assert_array_equal(c.series('dataset-a', 'M'), np.array([1.0]))


def test_series_for_absent_metric_is_empty():
    assert _cache().series('dataset-a', 'NOT_COLLECTED').size == 0


def test_series_for_empty_list_is_empty():
    c = ValueCache('genome', 'x', True, ('M',), (CohortValues('dataset-a', 0, '1.33', 'dict', 0, {'M': []}),))
    assert c.series('dataset-a', 'M').size == 0


def test_require_usable_accepts_a_superset_cache():
    cache_mod.require_usable(_cache(metrics=('MEDIAN_COVERAGE', 'FREEMIX', 'PCT_20X')), SPEC)  # does not raise


def test_require_usable_rejects_missing_metrics():
    with pytest.raises(CacheError, match=r"missing.*\['FREEMIX'\].*qc_calibrate collect"):
        cache_mod.require_usable(_cache(metrics=('MEDIAN_COVERAGE',)), SPEC)


def test_require_usable_rejects_incomplete_cache():
    with pytest.raises(CacheError, match='incomplete'):
        cache_mod.require_usable(_cache(complete=False), SPEC)


def test_require_usable_reports_incompleteness_before_missing_metrics():
    """An incomplete cache is the more actionable problem; report it first."""
    with pytest.raises(CacheError, match='incomplete'):
        cache_mod.require_usable(_cache(complete=False, metrics=('MEDIAN_COVERAGE',)), SPEC)
