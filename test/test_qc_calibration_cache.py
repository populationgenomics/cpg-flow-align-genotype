"""Unit tests for the calibration value cache."""

import io
import json
from pathlib import Path

import numpy as np
import pytest

import cpg_utils

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


class _FakeCloudFile:
    """File-like double returned by `_FakeCloudPath.open`; records what got written."""

    def __init__(self, owner: '_FakeCloudPath') -> None:
        self._owner = owner
        self._buffer = io.StringIO()

    def __enter__(self) -> io.StringIO:
        return self._buffer

    def __exit__(self, *exc_info: object) -> bool:
        self._owner.contents = self._buffer.getvalue()
        return False


class _FakeCloudPath:
    """Minimal double for `cloudpathlib.CloudPath`: has `.open` but is not a `pathlib.Path`.

    `save` branches on `isinstance(target, Path)` to decide whether the tmp-then-replace
    dance is safe; this stands in for the "no" case without touching the network.
    """

    def __init__(self) -> None:
        self.contents: str | None = None

    def open(self, mode: str = 'r') -> _FakeCloudFile:
        assert mode == 'w'
        return _FakeCloudFile(self)


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


def test_save_over_existing_cache_is_atomic_on_nan_failure(tmp_path):
    """A failed save must not destroy a previously-saved valid cache at the same path.

    `save` used to stream JSON straight into the open target file; hitting a NaN
    partway through left a truncated file in place of whatever was there before -
    destroying up to ten minutes of report-parsing work for no benefit.
    """
    path = tmp_path / 'values.json'
    good = _cache()
    cache_mod.save(good, path)

    bad = ValueCache(
        seq_type='genome',
        generated='x',
        complete=True,
        metrics=('M',),
        cohorts=(CohortValues('dataset-a', 1, '1.33', 'dict', 0, {'M': [float('nan')]}),),
    )
    with pytest.raises(ValueError, match=r'Out of range float values|NaN'):
        cache_mod.save(bad, path)

    assert cache_mod.load(path) == good
    assert list(tmp_path.glob('*.tmp')) == []


def test_save_leaves_no_tmp_file_on_success(tmp_path):
    path = tmp_path / 'values.json'
    cache_mod.save(_cache(), path)
    assert list(tmp_path.glob('*.tmp')) == []


def test_save_over_existing_cache_survives_replace_failure(tmp_path, monkeypatch):
    """Pins the actual atomic-write mechanism, not just the pre-write NaN guard.

    `test_save_over_existing_cache_is_atomic_on_nan_failure` above passes purely
    because `json.dumps` raises before any file is touched - it never exercises the
    temp-file-then-replace machinery at all. This forces `Path.replace` itself to fail
    *after* the temp file has been fully written, which is the scenario that
    machinery exists to survive (e.g. a crash or a full disk during the rename).
    """
    path = tmp_path / 'values.json'
    good = _cache()
    cache_mod.save(good, path)

    def boom(_self: Path, _target: Path) -> Path:
        raise OSError('simulated failure during rename')

    monkeypatch.setattr(Path, 'replace', boom)

    with pytest.raises(OSError, match='simulated failure during rename'):
        cache_mod.save(_cache(metrics=('MEDIAN_COVERAGE',)), path)

    assert cache_mod.load(path) == good
    assert list(tmp_path.glob('*.tmp')) == []


def test_save_writes_directly_for_non_local_targets(monkeypatch):
    """A CloudPath-like target isn't a `pathlib.Path`, so `save` must skip the
    tmp-then-replace dance entirely and write straight through `.open('w')` - see the
    docstring on `save` for why that dance is actively harmful on a real CloudPath."""
    fake = _FakeCloudPath()
    monkeypatch.setattr(cpg_utils, 'to_path', lambda _path: fake)

    cache_mod.save(_cache(), 'gs://some-bucket/values.json')

    assert fake.contents is not None
    written = json.loads(fake.contents)
    assert written['seq_type'] == 'genome'
    assert written['cohorts']['dataset-a']['values']['MEDIAN_COVERAGE'] == [30.0, 28.0, 12.0]


def test_save_raises_cache_error_on_non_numeric_value(tmp_path):
    """`series` tolerates a bad hand-edited value on read; `save` must not write one
    back out silently - see `_coerce_metric_values`."""
    c = ValueCache(
        seq_type='genome',
        generated='x',
        complete=True,
        metrics=('M',),
        cohorts=(CohortValues('dataset-a', 2, '1.33', 'dict', 0, {'M': [1.0, None]}),),
    )
    with pytest.raises(CacheError, match=r"dataset-a.*'M'.*non-numeric"):
        cache_mod.save(c, tmp_path / 'values.json')


def test_save_handles_numpy_floats(tmp_path):
    """collect casts to float, but the sink shouldn't corrupt if one slips through.

    `np.float64` subclasses `float`, so it round-trips through `json.dumps` even
    without a cast - it's `np.int64` that actually pins the `float()` cast, since
    `json.dumps` raises `TypeError` on a bare numpy integer.
    """
    c = ValueCache(
        seq_type='genome',
        generated='x',
        complete=True,
        metrics=('M',),
        cohorts=(CohortValues('dataset-a', 2, '1.33', 'dict', 0, {'M': [np.float64(1.5), np.int64(2)]}),),
    )
    path = tmp_path / 'values.json'
    cache_mod.save(c, path)
    assert json.loads(path.read_text())['cohorts']['dataset-a']['values']['M'] == [1.5, 2.0]


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


def test_series_raises_cache_error_on_non_numeric_junk():
    """A hand-edited cache can contain outright junk, not just `None` - `series` should
    name the cohort and metric rather than let a bare numpy exception through."""
    c = ValueCache('genome', 'x', True, ('M',), (CohortValues('dataset-a', 1, '1.33', 'dict', 0, {'M': ['abc']}),))
    with pytest.raises(CacheError, match=r"dataset-a.*'M'.*non-numeric"):
        c.series('dataset-a', 'M')


def test_load_raises_cache_error_on_missing_key(tmp_path):
    path = tmp_path / 'values.json'
    path.write_text(json.dumps({'seq_type': 'genome'}))
    with pytest.raises(CacheError, match='not a usable value cache'):
        cache_mod.load(path)


def test_load_raises_cache_error_on_wrong_shaped_value(tmp_path):
    path = tmp_path / 'values.json'
    path.write_text(
        json.dumps({'seq_type': 'genome', 'generated': 'x', 'complete': True, 'metrics': [], 'cohorts': 'nope'})
    )
    with pytest.raises(CacheError, match='not a usable value cache'):
        cache_mod.load(path)


def test_load_raises_cache_error_on_malformed_json(tmp_path):
    path = tmp_path / 'values.json'
    path.write_text('{not valid json')
    with pytest.raises(CacheError, match='not a usable value cache'):
        cache_mod.load(path)


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
