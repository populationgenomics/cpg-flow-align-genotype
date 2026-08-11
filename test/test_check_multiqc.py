"""Unit tests for check_multiqc.run() with warn/fail severity tiers.

These run fully offline: config access and AR-GUID lookup are monkeypatched, so
no cpg-utils config file or Slack access is needed. The MultiQC JSON is a small
hand-built ``report_general_stats_data`` blob matching the real shape (see
testing_scripts/testing_data/*_multiqc*.json).
"""

import json

import pytest

from align_genotype.scripts import check_multiqc

# Threshold sets mirroring config_template.toml :: qc_thresholds, nested by
# severity tier then direction: {severity: {direction: {metric: threshold}}}.
GENOME_THRESHOLDS = {
    'fail': {
        'min': {'MEDIAN_COVERAGE': 10, 'reads_mapped_percent': 80},
        'max': {'FREEMIX': 0.04, 'reads_duplicated_percent': 25},
    },
}
EXOME_THRESHOLDS = {
    'fail': {
        'min': {'MEAN_TARGET_COVERAGE': 50, 'PCT_TARGET_BASES_20X': 0.90, 'reads_mapped_percent': 80},
        'max': {'FREEMIX': 0.04, 'FOLD_80_BASE_PENALTY': 3.0, 'ZERO_CVG_TARGETS_PCT': 0.05},
    },
    'warn': {
        'min': {'MEAN_TARGET_COVERAGE': 80, 'PCT_TARGET_BASES_20X': 0.95, 'reads_mapped_percent': 98},
        'max': {'FREEMIX': 0.01, 'FOLD_80_BASE_PENALTY': 2.0, 'ZERO_CVG_TARGETS_PCT': 0.01},
    },
}


@pytest.fixture
def patch_config(monkeypatch):
    """Return a helper that wires up config for a given seq type + threshold set."""

    def _apply(seq_type: str, thresholds: dict) -> None:
        def config_retrieve(keys, default=None):  # noqa: ANN202
            if keys == ['workflow', 'sequencing_type']:
                return seq_type
            # Nested tier lookup: ['qc_thresholds', <seq_type>, <severity>, <direction>]
            if len(keys) == 4 and keys[0] == 'qc_thresholds':
                _, st, severity, direction = keys
                empty = {} if default is None else default
                if st != seq_type:
                    return empty
                return thresholds.get(severity, {}).get(direction, empty)
            return default

        monkeypatch.setattr(check_multiqc.config, 'config_retrieve', config_retrieve)
        monkeypatch.setattr(check_multiqc.config, 'try_get_ar_guid', lambda: 'test-ar-guid')

    return _apply


def _write_json(tmp_path, sections: dict) -> str:
    path = tmp_path / 'multiqc_data.json'
    path.write_text(json.dumps({'report_general_stats_data': sections}))
    return str(path)


def _run(json_path, output_path) -> dict:
    return check_multiqc.run(
        multiqc_json_path=json_path,
        html_url='http://example/report.html',
        dataset='validation-test',
        title='CRAM check',
        send_to_slack=False,
        output_json_path=str(output_path),
    )


def _flags_by_metric(result: dict, sg: str) -> dict:
    return {f['flag']: f for f in result['qc_flags'].get(sg, [])}


# --- genome path (single fail tier) ----------------------------------------
def test_genome_flags_min_and_max(tmp_path, patch_config):
    patch_config('genome', GENOME_THRESHOLDS)
    sections = {
        'picard_4': {'CPG1': {'MEDIAN_COVERAGE': 8}},  # below min -> fail
        'samtools': {'CPG1': {'reads_mapped_percent': 99, 'reads_duplicated_percent': 30}},  # dup above max -> fail
        'verifybamid': {'CPG1': {'FREEMIX': 0.01}},  # ok
    }
    result = _run(_write_json(tmp_path, sections), tmp_path / 'out.json')

    assert result['sequencing_type'] == 'genome'
    assert result['n_samples_flagged'] == 1
    flags = _flags_by_metric(result, 'CPG1')
    assert set(flags) == {'MEDIAN_COVERAGE', 'reads_duplicated_percent'}
    assert all(f['severity'] == 'fail' for f in flags.values())  # no warn tier configured


def test_genome_all_pass(tmp_path, patch_config):
    patch_config('genome', GENOME_THRESHOLDS)
    sections = {
        'picard_4': {'CPG1': {'MEDIAN_COVERAGE': 35}},
        'samtools': {'CPG1': {'reads_mapped_percent': 99, 'reads_duplicated_percent': 5}},
        'verifybamid': {'CPG1': {'FREEMIX': 0.005}},
    }
    result = _run(_write_json(tmp_path, sections), tmp_path / 'out.json')
    assert result['n_samples_flagged'] == 0
    assert result['qc_flags'] == {}


# --- warn/fail tiers -------------------------------------------------------
def test_exome_warn_and_fail_tiers(tmp_path, patch_config):
    patch_config('exome', EXOME_THRESHOLDS)
    sections = {
        'picard': {
            'CPG2': {
                'MEAN_TARGET_COVERAGE': 40,  # < 50 fail
                'PCT_TARGET_BASES_20X': 0.99,  # ok
                'FOLD_80_BASE_PENALTY': 2.5,  # > 2.0 warn (<= 3.0 fail)
                'ZERO_CVG_TARGETS_PCT': 0.005,  # ok
            }
        },
        'samtools': {'CPG2': {'reads_mapped_percent': 99}},
        'verifybamid': {'CPG2': {'FREEMIX': 0.005}},
    }
    result = _run(_write_json(tmp_path, sections), tmp_path / 'out.json')
    flags = _flags_by_metric(result, 'CPG2')
    assert flags['MEAN_TARGET_COVERAGE']['severity'] == 'fail'
    assert flags['MEAN_TARGET_COVERAGE']['threshold'] == 50
    assert flags['FOLD_80_BASE_PENALTY']['severity'] == 'warn'
    assert flags['FOLD_80_BASE_PENALTY']['threshold'] == 2.0
    assert set(flags) == {'MEAN_TARGET_COVERAGE', 'FOLD_80_BASE_PENALTY'}


def test_value_breaching_both_tiers_recorded_once_as_fail(tmp_path, patch_config):
    patch_config('exome', EXOME_THRESHOLDS)
    # FOLD_80 = 5.0 breaches warn (>2.0) AND fail (>3.0): recorded once, as fail.
    sections = {'picard': {'CPG3': {'FOLD_80_BASE_PENALTY': 5.0}}}
    result = _run(_write_json(tmp_path, sections), tmp_path / 'out.json')
    flags = result['qc_flags']['CPG3']
    assert len(flags) == 1
    assert flags[0]['severity'] == 'fail'
    assert flags[0]['threshold'] == 3.0


def test_warn_only_metric_emits_warn(tmp_path, patch_config):
    # A metric configured with only a warn tier still flags (as warn).
    patch_config('genome', {'warn': {'max': {'FREEMIX': 0.01}}})
    sections = {'verifybamid': {'CPG4': {'FREEMIX': 0.02}}}
    result = _run(_write_json(tmp_path, sections), tmp_path / 'out.json')
    assert result['qc_flags']['CPG4'][0]['severity'] == 'warn'


# --- exome path (previously unchecked) -------------------------------------
def test_exome_with_no_thresholds_checks_nothing(tmp_path, patch_config):
    # Regression guard: before exome thresholds existed, this returned zero flags.
    patch_config('exome', {})
    sections = {'picard': {'CPG2': {'MEAN_TARGET_COVERAGE': 5}}}
    result = _run(_write_json(tmp_path, sections), tmp_path / 'out.json')
    assert result['n_samples_flagged'] == 0


# --- silent-inert-metric warning -------------------------------------------
def test_warns_on_configured_metric_absent_from_report(tmp_path, patch_config, caplog):
    # PCT_PF_READS_ALIGNED is the classic case: configured but never surfaced.
    patch_config('genome', {'fail': {'min': {'PCT_PF_READS_ALIGNED': 0.8}}})
    sections = {'samtools': {'CPG1': {'reads_mapped_percent': 99}}}
    with caplog.at_level('WARNING'):
        _run(_write_json(tmp_path, sections), tmp_path / 'out.json')
    assert any('PCT_PF_READS_ALIGNED' in r.message and 'not found' in r.message for r in caplog.records)


# --- non-numeric metric values ---------------------------------------------
def test_non_numeric_metric_value_is_skipped_not_crashed(tmp_path, patch_config, caplog):
    # Picard writes '?' for FOLD_80_BASE_PENALTY when coverage is ~0; comparing it
    # against a float used to crash the whole check (seen on a real 647-sample WES run).
    patch_config('exome', EXOME_THRESHOLDS)
    sections = {
        'picard': {
            'CPG_OK': {'FOLD_80_BASE_PENALTY': 4.0},
            'CPG_BAD': {'FOLD_80_BASE_PENALTY': '?'},
        },
    }
    with caplog.at_level('WARNING'):
        result = _run(_write_json(tmp_path, sections), tmp_path / 'out.json')

    assert result['qc_flags']['CPG_OK'][0]['flag'] == 'FOLD_80_BASE_PENALTY'
    assert 'CPG_BAD' not in result['qc_flags']
    assert any("non-numeric value '?'" in r.message for r in caplog.records)


def test_numeric_string_values_are_coerced(tmp_path, patch_config):
    # Values that arrive as numeric strings should still be compared, not skipped.
    patch_config('genome', GENOME_THRESHOLDS)
    sections = {'picard_4': {'CPG1': {'MEDIAN_COVERAGE': '8'}}}  # string "8" < 10 -> flag
    result = _run(_write_json(tmp_path, sections), tmp_path / 'out.json')
    assert result['qc_flags']['CPG1'][0]['flag'] == 'MEDIAN_COVERAGE'


def test_output_json_written_and_structured(tmp_path, patch_config):
    patch_config('genome', GENOME_THRESHOLDS)
    sections = {'picard_4': {'CPG1': {'MEDIAN_COVERAGE': 8}}}
    out = tmp_path / 'out.json'
    _run(_write_json(tmp_path, sections), out)
    written = json.loads(out.read_text())
    assert written['dataset'] == 'validation-test'
    flag = written['qc_flags']['CPG1'][0]
    assert flag['flag'] == 'MEDIAN_COVERAGE'
    assert flag['ar_guid'] == 'test-ar-guid'
    assert flag['severity'] == 'fail'


if __name__ == '__main__':
    pytest.main([__file__, '-v'])
