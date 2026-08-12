"""Unit tests for the per-dataset values file."""

import json
from typing import Any

import numpy as np
import pytest

from align_genotype.qc_calibration import values as values_mod


def make_values(**overrides: Any) -> values_mod.DatasetValues:
    defaults = {
        'dataset': 'dataset-a',
        'seq_type': 'genome',
        'analysis_id': 42,
        'timestamp': '2026-06-01T00:00:00',
        'uri': 'gs://bucket/multiqc_data.json',
        'multiqc_version': '1.33',
        'generated': '2026-08-12T00:00:00',
        'n_sequencing_groups': 3,
        'section_sizes': {'picard_1': 3},
        'metrics': {
            'MEDIAN_COVERAGE': values_mod.MetricValues(
                entries=(('picard_1', 'CPG1', 30.0), ('picard_1', 'CPG2', 34.0), ('picard_1', 'CPG3', 38.0)),
                n_dropped=1,
            ),
        },
    }
    return values_mod.DatasetValues(**{**defaults, **overrides})


def test_metric_values_exposes_an_array_in_entry_order():
    metric = make_values().metrics['MEDIAN_COVERAGE']
    np.testing.assert_array_equal(metric.array, np.array([30.0, 34.0, 38.0]))


def test_metric_values_counts_values_and_sequencing_groups_separately():
    """One sequencing group in two sections yields two values but one group."""
    metric = values_mod.MetricValues(
        entries=(('picard_1', 'CPG1', 1.0), ('picard_4', 'CPG1', 1.0), ('picard_1', 'CPG2', 2.0)),
        n_dropped=0,
    )
    assert metric.n_values == 3
    assert metric.n_sequencing_groups == 2
    assert metric.duplicated is True
    assert metric.sections == ('picard_1', 'picard_4')


def test_metric_values_not_duplicated_when_one_section():
    assert make_values().metrics['MEDIAN_COVERAGE'].duplicated is False


def test_empty_metric_values_gives_an_empty_array():
    metric = values_mod.MetricValues(entries=(), n_dropped=0)
    assert metric.array.size == 0
    assert metric.n_values == 0
    assert metric.n_sequencing_groups == 0
    assert metric.sections == ()
    assert metric.duplicated is False


def test_array_filters_non_finite_defensively():
    """extract already drops these; the file is JSON and could be hand-edited."""
    metric = values_mod.MetricValues(
        entries=(('s', 'CPG1', 1.0), ('s', 'CPG2', float('nan')), ('s', 'CPG3', float('inf'))),
        n_dropped=0,
    )
    np.testing.assert_array_equal(metric.array, np.array([1.0]))


def test_dataset_values_metric_returns_empty_for_an_absent_metric():
    assert make_values().metric('NOPE').entries == ()


def test_dataset_values_metric_returns_the_stored_metric_when_present():
    assert make_values().metric('MEDIAN_COVERAGE').n_dropped == 1


def test_round_trips_through_json(tmp_path):
    path = tmp_path / 'values.json'
    original = make_values()
    values_mod.save(original, path)
    assert values_mod.load(path) == original


def test_saved_json_is_readable_and_shaped_as_documented(tmp_path):
    path = tmp_path / 'values.json'
    values_mod.save(make_values(), path)
    payload = json.loads(path.read_text())
    assert payload['dataset'] == 'dataset-a'
    assert payload['analysis_id'] == 42
    assert payload['metrics']['MEDIAN_COVERAGE']['n_dropped'] == 1
    assert payload['metrics']['MEDIAN_COVERAGE']['entries'][0] == ['picard_1', 'CPG1', 30.0]


def test_save_refuses_a_non_finite_value(tmp_path):
    """A NaN reaching the file would poison every percentile downstream."""
    broken = make_values(
        metrics={'X': values_mod.MetricValues(entries=(('s', 'CPG1', float('nan')),), n_dropped=0)},
    )
    with pytest.raises(values_mod.ValuesError, match='non-finite'):
        values_mod.save(broken, tmp_path / 'values.json')


def test_save_writes_nothing_when_it_refuses(tmp_path):
    """A rejected save must not leave a partial file behind."""
    path = tmp_path / 'values.json'
    broken = make_values(
        metrics={'X': values_mod.MetricValues(entries=(('s', 'CPG1', float('nan')),), n_dropped=0)},
    )
    with pytest.raises(values_mod.ValuesError):
        values_mod.save(broken, path)
    assert not path.exists()


def test_load_names_the_file_when_the_shape_is_wrong(tmp_path):
    path = tmp_path / 'values.json'
    path.write_text('{"dataset": "a"}')
    with pytest.raises(values_mod.ValuesError, match=str(path)):
        values_mod.load(path)


def test_load_names_the_file_when_the_json_is_malformed(tmp_path):
    path = tmp_path / 'values.json'
    path.write_text('{not json')
    with pytest.raises(values_mod.ValuesError, match=str(path)):
        values_mod.load(path)


def test_load_rejects_a_boolean_count(tmp_path):
    """`bool` subclasses `int`, so a bare `int()` would silently turn `true` into `1`."""
    payload = {
        'dataset': 'a',
        'seq_type': 'genome',
        'analysis_id': True,
        'timestamp': 't',
        'uri': 'u',
        'multiqc_version': 'v',
        'generated': 'g',
        'n_sequencing_groups': 1,
        'section_sizes': {},
        'metrics': {},
    }
    path = tmp_path / 'values.json'
    path.write_text(json.dumps(payload))
    with pytest.raises(values_mod.ValuesError, match='integer'):
        values_mod.load(path)


def test_load_rejects_a_boolean_metric_value(tmp_path):
    """Same `bool`-is-an-`int` trap, but for a value that feeds every percentile."""
    payload = {
        'dataset': 'a',
        'seq_type': 'genome',
        'analysis_id': 1,
        'timestamp': 't',
        'uri': 'u',
        'multiqc_version': 'v',
        'generated': 'g',
        'n_sequencing_groups': 1,
        'section_sizes': {},
        'metrics': {'X': {'n_dropped': 0, 'entries': [['s', 'CPG1', True]]}},
    }
    path = tmp_path / 'values.json'
    path.write_text(json.dumps(payload))
    with pytest.raises(values_mod.ValuesError, match='number'):
        values_mod.load(path)
