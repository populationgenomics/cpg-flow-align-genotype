"""Unit tests for extracting one MultiQC report into a values file."""

import numpy as np
import pytest

from align_genotype.qc_calibration import extract as extract_mod
from align_genotype.qc_calibration import settings as settings_mod

SETTINGS = settings_mod.CalibrationSettings(
    seq_type='genome',
    metrics=(
        settings_mod.MetricSpec(key='MEDIAN_COVERAGE', direction='min', unit='x'),
        settings_mod.MetricSpec(key='FREEMIX', direction='max', unit='frac'),
    ),
)

PROVENANCE = {
    'dataset': 'dataset-a',
    'analysis_id': 42,
    'timestamp': '2026-06-01T00:00:00',
    'uri': 'gs://bucket/multiqc_data.json',
    'generated': '2026-08-12T00:00:00',
}


def document(general_stats, version='1.33') -> dict:
    return {'config_version': version, 'report_general_stats_data': general_stats}


def test_extracts_values_from_the_v133_dict_shape():
    doc = document(
        {
            'picard_1': {'CPG1': {'MEDIAN_COVERAGE': 30.0}, 'CPG2': {'MEDIAN_COVERAGE': 34.0}},
            'verifybamid': {'CPG1': {'FREEMIX': 0.001}},
        },
    )
    result = extract_mod.extract(doc, SETTINGS, **PROVENANCE)
    assert result.dataset == 'dataset-a'
    assert result.seq_type == 'genome'
    assert result.analysis_id == 42
    assert result.multiqc_version == '1.33'
    assert result.n_sequencing_groups == 2
    assert result.section_sizes == {'picard_1': 2, 'verifybamid': 1}
    np.testing.assert_array_equal(result.metric('MEDIAN_COVERAGE').array, np.array([30.0, 34.0]))
    np.testing.assert_array_equal(result.metric('FREEMIX').array, np.array([0.001]))


def test_extracts_values_from_the_v114_list_shape():
    """v1.14 stores general stats as a positional list; both shapes exist in real reports."""
    doc = document([{'CPG1': {'MEDIAN_COVERAGE': 30.0}}, {'CPG1': {'FREEMIX': 0.001}}], version='1.14')
    result = extract_mod.extract(doc, SETTINGS, **PROVENANCE)
    assert sorted(result.section_sizes) == ['section_0', 'section_1']
    np.testing.assert_array_equal(result.metric('MEDIAN_COVERAGE').array, np.array([30.0]))


def test_strips_the_rich_id_suffix_to_get_the_sequencing_group():
    """MultiQC runs with --replace-names, so sample keys can read CPG1|EXTID."""
    doc = document({'picard_1': {'CPG1|EXT1': {'MEDIAN_COVERAGE': 30.0}}})
    result = extract_mod.extract(doc, SETTINGS, **PROVENANCE)
    assert result.metric('MEDIAN_COVERAGE').entries == (('picard_1', 'CPG1', 30.0),)
    assert result.n_sequencing_groups == 1


def test_one_group_in_two_sections_yields_two_values_and_one_group():
    doc = document(
        {
            'picard_1': {'CPG1': {'MEDIAN_COVERAGE': 30.0}},
            'picard_4': {'CPG1': {'MEDIAN_COVERAGE': 30.0}},
        },
    )
    metric = extract_mod.extract(doc, SETTINGS, **PROVENANCE).metric('MEDIAN_COVERAGE')
    assert (metric.n_values, metric.n_groups_with_values) == (2, 1)


def test_picard_question_mark_placeholder_is_dropped_and_counted():
    doc = document({'picard_1': {'CPG1': {'MEDIAN_COVERAGE': '?'}, 'CPG2': {'MEDIAN_COVERAGE': 34.0}}})
    metric = extract_mod.extract(doc, SETTINGS, **PROVENANCE).metric('MEDIAN_COVERAGE')
    np.testing.assert_array_equal(metric.array, np.array([34.0]))
    assert metric.n_dropped == 1


def test_non_finite_value_is_dropped_and_counted():
    doc = document({'picard_1': {'CPG1': {'MEDIAN_COVERAGE': float('nan')}, 'CPG2': {'MEDIAN_COVERAGE': 34.0}}})
    metric = extract_mod.extract(doc, SETTINGS, **PROVENANCE).metric('MEDIAN_COVERAGE')
    np.testing.assert_array_equal(metric.array, np.array([34.0]))
    assert metric.n_dropped == 1


def test_numeric_strings_are_coerced():
    doc = document({'picard_1': {'CPG1': {'MEDIAN_COVERAGE': '30.5'}}})
    metric = extract_mod.extract(doc, SETTINGS, **PROVENANCE).metric('MEDIAN_COVERAGE')
    np.testing.assert_array_equal(metric.array, np.array([30.5]))
    assert metric.n_dropped == 0


def test_boolean_value_is_coerced_to_a_float_not_dropped():
    """Pinning inherited behaviour, not endorsing it.

    `check_multiqc.gather_metric_values` coerces with a bare `float()`, and
    `float(True) == 1.0` raises nothing, so a JSON `true` is silently accepted as a real
    reading rather than dropped as non-numeric. `math.isfinite(True)` is also `True`, so
    this module's own finiteness filter does not catch it either. This module goes
    through `check_multiqc`'s functions deliberately so calibration matches enforcement
    exactly, so the fix (if any) belongs there, not here - this test exists so the
    behaviour is documented rather than rediscovered as a surprise.
    """
    doc = document({'picard_1': {'CPG1': {'MEDIAN_COVERAGE': True}}})
    metric = extract_mod.extract(doc, SETTINGS, **PROVENANCE).metric('MEDIAN_COVERAGE')
    np.testing.assert_array_equal(metric.array, np.array([1.0]))
    assert metric.n_dropped == 0


def test_a_metric_absent_from_the_report_is_present_but_empty():
    """Absence must be recorded, not omitted - the report banner is built from this."""
    doc = document({'picard_1': {'CPG1': {'MEDIAN_COVERAGE': 30.0}}})
    result = extract_mod.extract(doc, SETTINGS, **PROVENANCE)
    assert 'FREEMIX' in result.metrics
    assert result.metric('FREEMIX').entries == ()


def test_a_document_that_is_not_an_object_is_an_error():
    with pytest.raises(extract_mod.ExtractError, match='not an object'):
        extract_mod.extract([], SETTINGS, **PROVENANCE)


def test_no_usable_general_stats_is_an_error():
    """The key being absent/malformed - we couldn't read the report - is one distinct failure.

    Refusing beats reporting a clean extraction on a report we could not read. This must
    not share wording with the "present but empty" case below: the two point an operator
    at different fixes (chase the file/credentials, versus exclude the dataset).
    """
    with pytest.raises(extract_mod.ExtractError, match=r'dataset-a: no usable report_general_stats_data.*gs://bucket'):
        extract_mod.extract(document(None), SETTINGS, **PROVENANCE)


def test_empty_general_stats_is_an_error_naming_the_dataset_and_uri():
    """The key being present but genuinely holding zero modules is the other distinct failure.

    Mirrors `check_multiqc.load_sections`'s two-branch split; this message must not match
    the "no usable" one above, or a future regression that re-merges the two would pass.
    """
    with pytest.raises(extract_mod.ExtractError, match=r'dataset-a.*present but empty.*gs://bucket'):
        extract_mod.extract(document({}), SETTINGS, **PROVENANCE)


def test_missing_config_version_is_recorded_as_unknown():
    result = extract_mod.extract(
        {'report_general_stats_data': {'picard_1': {'CPG1': {'MEDIAN_COVERAGE': 30.0}}}},
        SETTINGS,
        **PROVENANCE,
    )
    assert result.multiqc_version == 'unknown'


def test_explicit_null_config_version_is_recorded_as_unknown():
    """`str(None)` reads 'None', which is worse than not printing a version at all."""
    doc = document({'picard_1': {'CPG1': {'MEDIAN_COVERAGE': 30.0}}}, version=None)
    result = extract_mod.extract(doc, SETTINGS, **PROVENANCE)
    assert result.multiqc_version == 'unknown'


def test_the_result_round_trips_through_the_values_file(tmp_path):
    """extract's output is what the report stage reads, so it must survive save/load."""
    from align_genotype.qc_calibration import values as values_mod  # noqa: PLC0415

    doc = document({'picard_1': {'CPG1': {'MEDIAN_COVERAGE': 30.0}}, 'verifybamid': {'CPG1': {'FREEMIX': 0.001}}})
    result = extract_mod.extract(doc, SETTINGS, **PROVENANCE)
    path = tmp_path / 'values.json'
    values_mod.save(result, path)
    assert values_mod.load(path) == result
