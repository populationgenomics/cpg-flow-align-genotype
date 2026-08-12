"""Unit tests for the calibration stages' output paths and inert-by-default behaviour.

The stage classes are not instantiated here: `Stage.__init__` calls `get_workflow()`,
which raises unless a Workflow singleton exists. The `forced` flags and the DAG wiring
are verified by the dry run in Task 18, which prints `[forced]` beside a forced stage.
"""

from types import SimpleNamespace

import pytest

from cpg_utils import to_path

from align_genotype import qc_calibration_stages as stages_mod
from align_genotype.qc_calibration.discovery import MultiqcReport

REPORT = MultiqcReport(dataset='ds-a', uri='gs://x/multiqc_data.json', analysis_id=42, timestamp='t')


def fake_dataset(name: str) -> SimpleNamespace:
    return SimpleNamespace(
        name=name,
        prefix=lambda: to_path(f'gs://cpg-{name}-main'),
        web_prefix=lambda: to_path(f'gs://cpg-{name}-web'),
    )


@pytest.fixture
def patch_lookup(monkeypatch):
    """Control what discovery returns, the run's sequencing type, and the enabled flag."""

    def _apply(report, *, seq_type='genome', enabled=True):  # noqa: ANN202
        monkeypatch.setattr(stages_mod.settings, 'enabled', lambda: enabled)
        monkeypatch.setattr(stages_mod, 'sequencing_type', lambda: seq_type)
        monkeypatch.setattr(stages_mod, 'dataset_report', lambda _name: report)

    return _apply


def test_values_path_is_keyed_on_the_analysis_id():
    """A new MultiQC report means a new path, so extraction never reads a stale file."""
    path = stages_mod.dataset_values_path(fake_dataset('ds-a'), 'genome', REPORT)
    assert str(path) == 'gs://cpg-ds-a-main/qc_calibration/genome/values.42.json'


def test_values_path_is_prefixed_by_the_literal_sequencing_type():
    """Not sequencing_subdir(), which returns '' for genome and would look asymmetric."""
    path = stages_mod.dataset_values_path(fake_dataset('ds-a'), 'exome', REPORT)
    assert '/qc_calibration/exome/' in str(path)


def test_dataset_outputs_are_empty_when_calibration_is_disabled(patch_lookup):
    patch_lookup(REPORT, enabled=False)
    assert stages_mod.dataset_outputs(fake_dataset('ds-a')) == {}


def test_dataset_outputs_are_empty_without_a_report(patch_lookup):
    """cpg-flow treats a falsy expected output as reusable, so no job is queued."""
    patch_lookup(None)
    assert stages_mod.dataset_outputs(fake_dataset('ds-a')) == {}


def test_dataset_outputs_name_the_values_file_when_enabled(patch_lookup):
    patch_lookup(REPORT)
    outputs = stages_mod.dataset_outputs(fake_dataset('ds-a'))
    assert set(outputs) == {'values'}
    assert str(outputs['values']).endswith('values.42.json')


def test_report_outputs_are_empty_when_calibration_is_disabled(patch_lookup):
    patch_lookup(REPORT, enabled=False)
    assert stages_mod.report_outputs(fake_dataset('analysis')) == {}


def test_report_writes_json_to_main_and_html_to_web(patch_lookup):
    patch_lookup(REPORT)
    outputs = stages_mod.report_outputs(fake_dataset('analysis'))
    assert str(outputs['json']) == 'gs://cpg-analysis-main/qc_calibration/genome/calibration.json'
    assert str(outputs['html']) == 'gs://cpg-analysis-web/qc_calibration/genome/calibration.html'


def test_report_outputs_follow_the_sequencing_type(patch_lookup):
    patch_lookup(REPORT, seq_type='exome')
    assert '/qc_calibration/exome/' in str(stages_mod.report_outputs(fake_dataset('analysis'))['json'])


def test_collect_values_paths_separates_datasets_with_and_without_reports(monkeypatch):
    monkeypatch.setattr(stages_mod, 'dataset_report', lambda name: REPORT if name == 'ds-a' else None)
    found, skipped = stages_mod.collect_values_paths([fake_dataset('ds-a'), fake_dataset('ds-b')], 'genome')
    assert set(found) == {'ds-a'}
    assert skipped == ['ds-b']


def test_dataset_report_resolves_the_dataset_name_before_querying_discovery(monkeypatch):
    """`dataset_for_access_level` adds the -test suffix at test access level, and the
    resolved name - not the raw dataset name - is what goes into the GraphQL query.

    The other tests here monkeypatch `dataset_report` itself, so nothing exercises its
    own composition of `config.dataset_for_access_level` and `discovery.
    cached_latest_cram_multiqc` without this.
    """
    calls: list[tuple[str, str]] = []

    monkeypatch.setattr(stages_mod, 'sequencing_type', lambda: 'genome')
    monkeypatch.setattr(stages_mod.config, 'dataset_for_access_level', lambda name: f'{name}-test')
    monkeypatch.setattr(
        stages_mod.discovery,
        'cached_latest_cram_multiqc',
        lambda dataset, seq_type: calls.append((dataset, seq_type)),
    )

    stages_mod.dataset_report('ds-a')

    assert calls == [('ds-a-test', 'genome')]


def test_calibration_stages_are_wired_into_the_entrypoint():
    from align_genotype import run_workflow  # noqa: PLC0415

    assert stages_mod.QcCalibrationDatasetMetrics in run_workflow.STAGES
    assert stages_mod.QcCalibrationReport in run_workflow.STAGES
