"""Unit tests for the two job entrypoint scripts, run in-process."""

import json
import logging

import pytest
from click.testing import CliRunner

from align_genotype.qc_calibration import settings as settings_mod
from align_genotype.qc_calibration import values as values_mod
from align_genotype.scripts import qc_calibration_extract

CALIBRATION_CONFIG = {
    'workflow': {'sequencing_type': 'genome'},
    'qc_calibration': {
        'genome': {
            'metrics': {
                'MEDIAN_COVERAGE': {'direction': 'min', 'unit': 'x'},
                'dup_pct': {'direction': 'max', 'unit': '%', 'relative': True},
            },
        },
    },
    'qc_thresholds': {'genome': {'fail': {'min': {'MEDIAN_COVERAGE': 15}}}},
}

REPORT = {
    'config_version': '1.33',
    'report_general_stats_data': {
        'picard_1': {f'CPG{i}': {'MEDIAN_COVERAGE': 30.0 + i} for i in range(6)},
        'samtools': {f'CPG{i}': {'dup_pct': 10.0 + i} for i in range(6)},
    },
}


@pytest.fixture
def patch_config(monkeypatch):
    """Point every config_retrieve used by the scripts at CALIBRATION_CONFIG."""

    def config_retrieve(keys, default=None):  # noqa: ANN202
        node = CALIBRATION_CONFIG
        for key in keys:
            if not isinstance(node, dict) or key not in node:
                return default
            node = node[key]
        return node

    from align_genotype.scripts import check_multiqc  # noqa: PLC0415

    monkeypatch.setattr(settings_mod.config, 'config_retrieve', config_retrieve)
    monkeypatch.setattr(check_multiqc.config, 'config_retrieve', config_retrieve)


def test_extract_writes_a_values_file(tmp_path, patch_config):  # noqa: ARG001
    report_path = tmp_path / 'multiqc_data.json'
    report_path.write_text(json.dumps(REPORT))
    output = tmp_path / 'values.json'

    result = CliRunner().invoke(
        qc_calibration_extract.main,
        [
            '--dataset',
            'ds-a',
            '--multiqc-json',
            str(report_path),
            '--analysis-id',
            '42',
            '--timestamp',
            '2026-06-01T00:00:00',
            '--uri',
            'gs://ds-a/multiqc_data.json',
            '--output',
            str(output),
        ],
    )
    assert result.exit_code == 0, result.output

    loaded = values_mod.load(output)
    assert loaded.dataset == 'ds-a'
    assert loaded.analysis_id == 42
    assert loaded.uri == 'gs://ds-a/multiqc_data.json'
    assert loaded.n_sequencing_groups == 6
    assert loaded.metric('MEDIAN_COVERAGE').n_values == 6


def test_extract_fails_loudly_on_an_unreadable_report(tmp_path, patch_config):  # noqa: ARG001
    report_path = tmp_path / 'multiqc_data.json'
    report_path.write_text(json.dumps({'config_version': '1.33'}))

    result = CliRunner().invoke(
        qc_calibration_extract.main,
        [
            '--dataset',
            'ds-a',
            '--multiqc-json',
            str(report_path),
            '--analysis-id',
            '42',
            '--timestamp',
            '2026-06-01T00:00:00',
            '--uri',
            'gs://ds-a/multiqc_data.json',
            '--output',
            str(tmp_path / 'values.json'),
        ],
    )
    assert result.exit_code != 0
    assert 'no usable report_general_stats_data' in str(result.exception)


def test_extract_warns_when_a_configured_metric_is_absent_from_the_report(tmp_path, patch_config, caplog):  # noqa: ARG001
    # No `samtools` section at all, so `dup_pct` (configured) has zero entries: the
    # `if not metric_values.entries` branch in main(), not exercised by the success test
    # above where every configured metric is present.
    report = {
        'config_version': '1.33',
        'report_general_stats_data': {
            'picard_1': {f'CPG{i}': {'MEDIAN_COVERAGE': 30.0 + i} for i in range(6)},
        },
    }
    report_path = tmp_path / 'multiqc_data.json'
    report_path.write_text(json.dumps(report))

    with caplog.at_level(logging.WARNING):
        result = CliRunner().invoke(
            qc_calibration_extract.main,
            [
                '--dataset',
                'ds-a',
                '--multiqc-json',
                str(report_path),
                '--analysis-id',
                '42',
                '--timestamp',
                '2026-06-01T00:00:00',
                '--uri',
                'gs://ds-a/multiqc_data.json',
                '--output',
                str(tmp_path / 'values.json'),
            ],
        )
    assert result.exit_code == 0, result.output
    assert 'ds-a: dup_pct is absent from this report' in caplog.text


def test_extract_warns_when_a_metric_lost_unusable_values(tmp_path, patch_config, caplog):  # noqa: ARG001
    # One non-numeric `dup_pct` cell among otherwise-usable ones: the
    # `elif metric_values.n_dropped` branch in main(), not exercised by either test above.
    report = {
        'config_version': '1.33',
        'report_general_stats_data': {
            'picard_1': {f'CPG{i}': {'MEDIAN_COVERAGE': 30.0 + i} for i in range(6)},
            'samtools': {
                'CPG0': {'dup_pct': '?'},
                **{f'CPG{i}': {'dup_pct': 10.0 + i} for i in range(1, 6)},
            },
        },
    }
    report_path = tmp_path / 'multiqc_data.json'
    report_path.write_text(json.dumps(report))

    with caplog.at_level(logging.WARNING):
        result = CliRunner().invoke(
            qc_calibration_extract.main,
            [
                '--dataset',
                'ds-a',
                '--multiqc-json',
                str(report_path),
                '--analysis-id',
                '42',
                '--timestamp',
                '2026-06-01T00:00:00',
                '--uri',
                'gs://ds-a/multiqc_data.json',
                '--output',
                str(tmp_path / 'values.json'),
            ],
        )
    assert result.exit_code == 0, result.output
    assert 'ds-a: dup_pct lost 1 unusable value(s)' in caplog.text
