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


def _write_values(tmp_path, name, coverage, dup):  # noqa: ANN202
    from align_genotype.qc_calibration import values as v  # noqa: PLC0415

    path = tmp_path / f'{name}.json'
    v.save(
        v.DatasetValues(
            dataset=name,
            seq_type='genome',
            analysis_id=1,
            timestamp='2026-06-01T00:00:00',
            uri=f'gs://{name}/multiqc_data.json',
            multiqc_version='1.33',
            generated='2026-08-12T00:00:00',
            n_sequencing_groups=len(coverage),
            section_sizes={'picard_1': len(coverage)},
            metrics={
                'MEDIAN_COVERAGE': v.MetricValues(
                    entries=tuple(('picard_1', f'{name}-CPG{i}', float(x)) for i, x in enumerate(coverage)),
                    n_dropped=0,
                ),
                'dup_pct': v.MetricValues(
                    entries=tuple(('samtools', f'{name}-CPG{i}', float(x)) for i, x in enumerate(dup)),
                    n_dropped=0,
                ),
            },
        ),
        path,
    )
    return path


def test_report_writes_json_and_html(tmp_path, patch_config, monkeypatch):  # noqa: ARG001
    from align_genotype.scripts import qc_calibration_report  # noqa: PLC0415

    monkeypatch.setattr(qc_calibration_report.config, 'try_get_ar_guid', lambda: 'test-ar-guid')
    a = _write_values(tmp_path, 'ds-a', [30, 32, 34, 36, 38, 10], [10, 10.5, 11, 11.5, 12, 40])
    b = _write_values(tmp_path, 'ds-b', [40, 42, 44, 46, 48, 50], [7, 7.5, 8, 8.5, 9, 9.5])
    out_json = tmp_path / 'calibration.json'
    out_html = tmp_path / 'calibration.html'

    result = CliRunner().invoke(
        qc_calibration_report.main,
        [
            '--values', str(a), '--values', str(b),
            '--skipped-dataset', 'ds-c',
            '--output-json', str(out_json),
            '--output-html', str(out_html),
        ],
    )
    assert result.exit_code == 0, result.output

    payload = json.loads(out_json.read_text())
    assert payload['sequencing_type'] == 'genome'
    assert [d['dataset'] for d in payload['datasets']] == ['ds-a', 'ds-b']
    assert payload['skipped_datasets'] == [
        {'dataset': 'ds-c', 'reason': 'no completed CramMultiQC qc analysis for genome'},
    ]
    assert payload['metrics']['MEDIAN_COVERAGE']['current'] == {'fail': 15}
    assert out_html.read_text().startswith('<!DOCTYPE html>')


def test_report_succeeds_even_when_a_metric_is_missing_everywhere(tmp_path, patch_config, monkeypatch, caplog):  # noqa: ARG001
    """Failing would destroy the HTML that explains the problem - Hail only copies
    write_output targets on job success."""
    from align_genotype.scripts import qc_calibration_report  # noqa: PLC0415

    monkeypatch.setattr(qc_calibration_report.config, 'try_get_ar_guid', lambda: 'x')
    path = tmp_path / 'ds-a.json'
    from align_genotype.qc_calibration import values as v  # noqa: PLC0415

    v.save(
        v.DatasetValues(
            dataset='ds-a', seq_type='genome', analysis_id=1, timestamp='t',
            uri='gs://a/multiqc_data.json', multiqc_version='1.33', generated='g',
            n_sequencing_groups=0, section_sizes={}, metrics={},
        ),
        path,
    )
    out_html = tmp_path / 'calibration.html'
    with caplog.at_level('ERROR'):
        result = CliRunner().invoke(
            qc_calibration_report.main,
            [
                '--values', str(path),
                '--output-json', str(tmp_path / 'calibration.json'),
                '--output-html', str(out_html),
            ],
        )
    assert result.exit_code == 0, result.output
    assert out_html.exists()
    assert 'checks nothing' in caplog.text


def test_report_warns_if_metric_missing_from_some_datasets(tmp_path, patch_config, monkeypatch, caplog):  # noqa: ARG001
    """The narrower problem - present in at least one dataset - is a WARNING, not an
    ERROR, and must not trip the same 'checks nothing' banner as the everywhere-missing
    case above."""
    from align_genotype.qc_calibration import values as v  # noqa: PLC0415
    from align_genotype.scripts import qc_calibration_report  # noqa: PLC0415

    monkeypatch.setattr(qc_calibration_report.config, 'try_get_ar_guid', lambda: 'x')
    a = _write_values(tmp_path, 'ds-a', [30, 32, 34, 36, 38, 40], [10, 10.5, 11, 11.5, 12, 12.5])
    b_path = tmp_path / 'ds-b.json'
    v.save(
        v.DatasetValues(
            dataset='ds-b',
            seq_type='genome',
            analysis_id=1,
            timestamp='2026-06-01T00:00:00',
            uri='gs://ds-b/multiqc_data.json',
            multiqc_version='1.33',
            generated='2026-08-12T00:00:00',
            n_sequencing_groups=6,
            section_sizes={'picard_1': 6},
            metrics={
                'MEDIAN_COVERAGE': v.MetricValues(
                    entries=tuple(('picard_1', f'ds-b-CPG{i}', float(30 + i)) for i in range(6)),
                    n_dropped=0,
                ),
                # dup_pct is entirely absent from this dataset's metrics.
            },
        ),
        b_path,
    )

    with caplog.at_level('WARNING'):
        result = CliRunner().invoke(
            qc_calibration_report.main,
            [
                '--values', str(a), '--values', str(b_path),
                '--output-json', str(tmp_path / 'calibration.json'),
                '--output-html', str(tmp_path / 'calibration.html'),
            ],
        )
    assert result.exit_code == 0, result.output
    assert 'dup_pct was absent from 1 dataset(s): ds-b' in caplog.text
    assert 'checks nothing' not in caplog.text
    assert not any(record.levelno >= logging.ERROR for record in caplog.records)


