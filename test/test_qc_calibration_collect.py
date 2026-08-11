"""Unit tests for the collect step: survey + extraction in a single parse."""

import json
from typing import Any

import pytest

from align_genotype.qc_calibration import collect as collect_mod
from align_genotype.qc_calibration import manifest as manifest_mod
from align_genotype.qc_calibration import spec as spec_mod
from align_genotype.qc_calibration.collect import CollectError

SPEC = spec_mod.loads(
    'seq_type = "genome"\ncache = "c.json"\n'
    '\n[metrics.MEDIAN_COVERAGE]\ndirection = "min"\nunit = "x"\nfail = 15\n'
    '\n[metrics.FREEMIX]\ndirection = "max"\nfail = 0.04\n'
    '\n[metrics.error_rate]\ndirection = "max"\ngated = false\n',
)

DICT_SECTIONS = {
    'picard': {'CPG1|S1': {'MEDIAN_COVERAGE': 30}, 'CPG2|S2': {'MEDIAN_COVERAGE': '?'}},
    'verifybamid': {'CPG1|S1': {'FREEMIX': 0.001}, 'CPG2|S2': {'FREEMIX': 0.05}},
}


def _write_report(tmp_path, name: str, sections: Any, version: str = '1.33') -> str:
    path = tmp_path / name
    path.write_text(json.dumps({'config_version': version, 'report_general_stats_data': sections}))
    return str(path)


def _manifest(**cohorts: str) -> manifest_mod.Manifest:
    body = ''.join(f'[cohorts.{label}]\nuri = "{uri}"\n\n' for label, uri in cohorts.items())
    return manifest_mod.loads(f'seq_type = "genome"\ngenerated = "x"\n\n{body}')


# --- single cohort ---------------------------------------------------------------


def test_collect_cohort_extracts_values(tmp_path):
    uri = _write_report(tmp_path, 'a.json', DICT_SECTIONS)
    values, row = collect_mod.collect_cohort(manifest_mod.Cohort('dataset-a', uri), SPEC)
    assert values.values['MEDIAN_COVERAGE'] == [30.0]
    assert sorted(values.values['FREEMIX']) == [0.001, 0.05]
    assert row.multiqc_version == '1.33'
    assert row.shape == 'dict'


def test_collect_cohort_counts_non_numeric_drops(tmp_path):
    uri = _write_report(tmp_path, 'a.json', DICT_SECTIONS)
    values, row = collect_mod.collect_cohort(manifest_mod.Cohort('dataset-a', uri), SPEC)
    assert values.n_dropped == 1  # the Picard '?' placeholder
    assert row.n_dropped == 1


def test_collect_cohort_filters_non_finite(tmp_path):
    """float('nan') and 'inf' coerce fine but must not reach the cache."""
    sections = {
        'picard': {'S1': {'MEDIAN_COVERAGE': 30}, 'S2': {'MEDIAN_COVERAGE': 'nan'}, 'S3': {'MEDIAN_COVERAGE': 'inf'}},
    }
    uri = _write_report(tmp_path, 'a.json', sections)
    spec = spec_mod.loads(
        'seq_type = "genome"\ncache = "c.json"\n[metrics.MEDIAN_COVERAGE]\ndirection = "min"\nfail = 1\n',
    )
    values, _ = collect_mod.collect_cohort(manifest_mod.Cohort('dataset-a', uri), spec)
    assert values.values['MEDIAN_COVERAGE'] == [30.0]
    assert values.n_dropped == 2


def test_collect_cohort_dedupes_samples_across_sections(tmp_path):
    uri = _write_report(tmp_path, 'a.json', DICT_SECTIONS)
    values, _ = collect_mod.collect_cohort(manifest_mod.Cohort('dataset-a', uri), SPEC)
    assert values.n_samples == 2  # CPG1|S1 and CPG2|S2 appear in both sections


def test_collect_cohort_records_which_sections_carry_each_metric(tmp_path):
    uri = _write_report(tmp_path, 'a.json', DICT_SECTIONS)
    _, row = collect_mod.collect_cohort(manifest_mod.Cohort('dataset-a', uri), SPEC)
    assert row.where['MEDIAN_COVERAGE'] == ('picard',)
    assert row.where['FREEMIX'] == ('verifybamid',)
    assert row.where['error_rate'] == ()


def test_collect_cohort_handles_list_shaped_report(tmp_path):
    """MultiQC v1.14 stores general stats as a positional list."""
    uri = _write_report(
        tmp_path,
        'a.json',
        [{'S1': {'MEDIAN_COVERAGE': 30, 'FREEMIX': 0.01}}],
        version='1.14',
    )
    values, row = collect_mod.collect_cohort(manifest_mod.Cohort('dataset-a', uri), SPEC)
    assert row.shape == 'list'
    assert row.where['MEDIAN_COVERAGE'] == ('section_0',)
    assert values.values['MEDIAN_COVERAGE'] == [30.0]


def test_collect_cohort_presence_is_by_key_not_by_numeric_value(tmp_path):
    """A metric present but entirely non-numeric still counts as present; the drop count tells the story."""
    sections = {'picard': {'S1': {'MEDIAN_COVERAGE': '?'}}}
    uri = _write_report(tmp_path, 'a.json', sections)
    spec = spec_mod.loads(
        'seq_type = "genome"\ncache = "c.json"\n[metrics.MEDIAN_COVERAGE]\ndirection = "min"\nfail = 1\n',
    )
    values, row = collect_mod.collect_cohort(manifest_mod.Cohort('dataset-a', uri), spec)
    assert row.where['MEDIAN_COVERAGE'] == ('picard',)
    assert values.values['MEDIAN_COVERAGE'] == []
    assert values.n_dropped == 1


def test_collect_cohort_records_section_keys_for_rename_mapping(tmp_path):
    uri = _write_report(tmp_path, 'a.json', DICT_SECTIONS)
    _, row = collect_mod.collect_cohort(manifest_mod.Cohort('dataset-a', uri), SPEC)
    assert row.section_keys['picard'] == ('MEDIAN_COVERAGE',)


def test_collect_cohort_raises_without_general_stats(tmp_path):
    path = tmp_path / 'a.json'
    path.write_text(json.dumps({'config_version': '1.33', 'report_saved_raw_data': {}}))
    with pytest.raises(CollectError, match='no usable report_general_stats_data'):
        collect_mod.collect_cohort(manifest_mod.Cohort('dataset-a', str(path)), SPEC)


# --- all cohorts -----------------------------------------------------------------


def test_collect_all_builds_a_complete_cache(tmp_path):
    a = _write_report(tmp_path, 'a.json', DICT_SECTIONS)
    b = _write_report(tmp_path, 'b.json', DICT_SECTIONS)
    result = collect_mod.collect_all(
        _manifest(**{'dataset-a': a, 'dataset-b': b}),
        SPEC,
        generated='2026-08-11T00:00:00',
    )
    assert result.ok
    assert result.cache.complete is True
    assert result.cache.labels == ('dataset-a', 'dataset-b')
    assert result.cache.metrics == ('MEDIAN_COVERAGE', 'FREEMIX', 'error_rate')


def test_collect_all_flags_missing_gated_metric_and_marks_cache_incomplete(tmp_path):
    good = _write_report(tmp_path, 'a.json', DICT_SECTIONS)
    without_freemix = _write_report(tmp_path, 'b.json', {'picard': {'S1': {'MEDIAN_COVERAGE': 30}}})
    result = collect_mod.collect_all(_manifest(**{'dataset-a': good, 'dataset-b': without_freemix}), SPEC)
    assert not result.ok
    assert result.missing_gated == {'dataset-b': ('FREEMIX',)}
    assert result.cache.complete is False


def test_collect_all_ignores_missing_ungated_metric(tmp_path):
    """error_rate is absent from every fixture but is gated = false, so it isn't fatal."""
    a = _write_report(tmp_path, 'a.json', DICT_SECTIONS)
    result = collect_mod.collect_all(_manifest(**{'dataset-a': a}), SPEC)
    assert result.ok


def test_collect_all_records_unreadable_cohort_but_continues(tmp_path):
    good = _write_report(tmp_path, 'a.json', DICT_SECTIONS)
    result = collect_mod.collect_all(
        _manifest(**{'dataset-a': good, 'dataset-b': str(tmp_path / 'nope.json')}),
        SPEC,
    )
    assert not result.ok
    assert result.cache.labels == ('dataset-a',)
    assert [label for label, _ in result.failures] == ['dataset-b']


def test_collect_all_records_malformed_json_as_a_failure(tmp_path):
    good = _write_report(tmp_path, 'a.json', DICT_SECTIONS)
    bad = tmp_path / 'bad.json'
    bad.write_text('{not json')
    result = collect_mod.collect_all(_manifest(**{'dataset-a': good, 'dataset-b': str(bad)}), SPEC)
    assert not result.ok
    assert [label for label, _ in result.failures] == ['dataset-b']


def test_collect_all_rejects_seq_type_mismatch(tmp_path):
    a = _write_report(tmp_path, 'a.json', DICT_SECTIONS)
    exome_manifest = manifest_mod.loads(f'seq_type = "exome"\ngenerated = "x"\n\n[cohorts.dataset-a]\nuri = "{a}"\n')
    with pytest.raises(CollectError, match="manifest is for 'exome' but the spec is for 'genome'"):
        collect_mod.collect_all(exome_manifest, SPEC)
