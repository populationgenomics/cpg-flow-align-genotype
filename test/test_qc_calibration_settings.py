"""Unit tests for calibration settings read from the [qc_calibration] config block."""

import pytest

from cpg_utils import config as cpg_config

from align_genotype.qc_calibration import settings as settings_mod

METRICS = {
    'MEDIAN_COVERAGE': {'direction': 'min', 'unit': 'x'},
    'reads_duplicated_percent': {'direction': 'max', 'unit': '%', 'relative': True},
}

_MISSING = object()


@pytest.fixture
def patch_config(monkeypatch):
    """Wire config_retrieve to an in-memory nested dict.

    Mirrors the real `cpg_utils.config.config_retrieve`: a missing key with no default
    raises `ConfigError` rather than returning `None`, so tests can't pass for the wrong
    reason on a config read that omitted a default.
    """

    def _apply(tree: dict) -> None:
        def config_retrieve(keys, default=_MISSING):  # noqa: ANN202
            node = tree
            for key in keys:
                if not isinstance(node, dict) or key not in node:
                    if default is _MISSING:
                        raise cpg_config.ConfigError(f'missing config key: {list(keys)}')
                    return default
                node = node[key]
            return node

        monkeypatch.setattr(settings_mod.config, 'config_retrieve', config_retrieve)
        monkeypatch.setattr(settings_mod.check_multiqc.config, 'config_retrieve', config_retrieve)

    return _apply


def test_load_reads_metrics_for_the_runs_sequencing_type(patch_config):
    patch_config(
        {
            'workflow': {'sequencing_type': 'genome'},
            'qc_calibration': {'genome': {'metrics': METRICS}, 'exome': {'metrics': {}}},
        },
    )
    loaded = settings_mod.load()
    assert loaded.seq_type == 'genome'
    assert loaded.metric_keys == ('MEDIAN_COVERAGE', 'reads_duplicated_percent')
    assert loaded.metric('MEDIAN_COVERAGE').unit == 'x'
    assert loaded.relative_metrics == (loaded.metric('reads_duplicated_percent'),)


def test_load_applies_documented_defaults(patch_config):
    patch_config(
        {'workflow': {'sequencing_type': 'genome'}, 'qc_calibration': {'genome': {'metrics': METRICS}}},
    )
    loaded = settings_mod.load()
    assert loaded.k == pytest.approx(3.5)
    assert loaded.min_samples == 50
    assert loaded.bars.max_warn_rate == pytest.approx(0.10)
    assert loaded.bars.max_growth_churn == pytest.approx(0.02)
    assert loaded.bars.max_merge_churn == pytest.approx(0.05)


def test_load_honours_overrides(patch_config):
    patch_config(
        {
            'workflow': {'sequencing_type': 'exome'},
            'qc_calibration': {
                'k': 3.0,
                'min_samples': 25,
                'max_merge_churn': 0.5,
                'exome': {'metrics': METRICS},
            },
        },
    )
    loaded = settings_mod.load()
    assert loaded.k == pytest.approx(3.0)
    assert loaded.min_samples == 25
    assert loaded.bars.max_merge_churn == pytest.approx(0.5)


def test_load_without_metrics_for_this_seq_type_is_an_error(patch_config):
    patch_config(
        {'workflow': {'sequencing_type': 'exome'}, 'qc_calibration': {'genome': {'metrics': METRICS}}},
    )
    with pytest.raises(settings_mod.SettingsError, match="sequencing type 'exome'"):
        settings_mod.load()


def test_bad_direction_is_rejected():
    with pytest.raises(settings_mod.SettingsError, match='direction'):
        settings_mod.parse_metric('X', {'direction': 'up'})


def test_missing_direction_is_rejected():
    with pytest.raises(settings_mod.SettingsError, match='missing required key: direction'):
        settings_mod.parse_metric('X', {'unit': 'x'})


def test_bad_unit_is_rejected():
    with pytest.raises(settings_mod.SettingsError, match='unit'):
        settings_mod.parse_metric('X', {'direction': 'min', 'unit': 'furlongs'})


def test_unknown_metric_key_is_rejected():
    with pytest.raises(settings_mod.SettingsError, match='unknown key'):
        settings_mod.parse_metric('X', {'direction': 'min', 'gated': True})


def test_quoted_boolean_relative_is_rejected():
    """`bool("false")` is True, so a quoted boolean must not be coerced."""
    with pytest.raises(settings_mod.SettingsError, match='relative must be true or false'):
        settings_mod.parse_metric('X', {'direction': 'min', 'relative': 'false'})


def test_metric_defaults_to_frac_and_not_relative():
    metric = settings_mod.parse_metric('X', {'direction': 'min'})
    assert (metric.unit, metric.relative) == ('frac', False)


def test_metric_lookup_of_an_unknown_key_raises():
    loaded = settings_mod.CalibrationSettings(
        seq_type='genome',
        metrics=(settings_mod.MetricSpec(key='A', direction='min'),),
    )
    with pytest.raises(KeyError):
        loaded.metric('NOPE')


def test_enabled_defaults_to_false(patch_config):
    patch_config({'workflow': {'sequencing_type': 'genome'}, 'qc_calibration': {}})
    assert settings_mod.enabled() is False


def test_enabled_reads_the_flag(patch_config):
    patch_config({'workflow': {'sequencing_type': 'genome'}, 'qc_calibration': {'enabled': True}})
    assert settings_mod.enabled() is True


def test_enabled_rejects_a_quoted_boolean(patch_config):
    patch_config({'workflow': {'sequencing_type': 'genome'}, 'qc_calibration': {'enabled': 'false'}})
    with pytest.raises(settings_mod.SettingsError, match='must be true or false'):
        settings_mod.enabled()


def test_load_rejects_a_non_numeric_k(patch_config):
    patch_config(
        {
            'workflow': {'sequencing_type': 'genome'},
            'qc_calibration': {'k': 'abc', 'genome': {'metrics': METRICS}},
        },
    )
    with pytest.raises(settings_mod.SettingsError, match=r'qc_calibration\.k'):
        settings_mod.load()


def test_load_rejects_a_boolean_min_samples(patch_config):
    patch_config(
        {
            'workflow': {'sequencing_type': 'genome'},
            'qc_calibration': {'min_samples': True, 'genome': {'metrics': METRICS}},
        },
    )
    with pytest.raises(settings_mod.SettingsError, match=r'qc_calibration\.min_samples'):
        settings_mod.load()


def test_load_without_a_sequencing_type_raises(patch_config):
    patch_config({'qc_calibration': {'genome': {'metrics': METRICS}}})
    with pytest.raises(cpg_config.ConfigError):
        settings_mod.load()


def test_current_thresholds_reshapes_production_config(patch_config):
    patch_config(
        {
            'workflow': {'sequencing_type': 'genome'},
            'qc_thresholds': {
                'genome': {
                    'fail': {'min': {'MEDIAN_COVERAGE': 15}, 'max': {'FREEMIX': 0.04}},
                    'warn': {'min': {'MEDIAN_COVERAGE': 25}},
                },
            },
        },
    )
    assert settings_mod.current_thresholds('genome') == {
        'MEDIAN_COVERAGE': {'fail': 15, 'warn': 25},
        'FREEMIX': {'fail': 0.04},
    }
