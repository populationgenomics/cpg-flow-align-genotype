import json
from typing import NoReturn

import pytest

from align_genotype.scripts import create_cram_qc_analyses as script


def test_qc_prefix_genome_and_exome():
    assert script.qc_prefix('rdnow', 'genome') == 'gs://cpg-rdnow-main/qc'
    assert script.qc_prefix('rdnow', 'exome') == 'gs://cpg-rdnow-main/exome/qc'
    assert script.qc_prefix('rdnow-test', 'genome') == 'gs://cpg-rdnow-test/qc'


def test_expected_qc_files_matches_cpg_flow_layout():
    files = script.expected_qc_files('CPG1', 'rdnow', 'genome')
    assert files['samtools_stats'] == 'gs://cpg-rdnow-main/qc/samtools_stats/CPG1.samtools-stats'
    assert files['verifybamid'] == 'gs://cpg-rdnow-main/qc/verify_bamid/CPG1.verify-bamid.selfSM'
    assert files['picard_wgs_metrics'] == 'gs://cpg-rdnow-main/qc/picard_wgs_metrics/CPG1.picard-wgs-metrics'
    assert (
        files['markduplicates_metrics'] == 'gs://cpg-rdnow-main/qc/markduplicates_metrics/CPG1.markduplicates-metrics'
    )
    exome = script.expected_qc_files('CPG1', 'rdnow', 'exome')
    assert 'picard_wgs_metrics' not in exome
    assert exome['picard_hs_metrics'] == 'gs://cpg-rdnow-main/exome/qc/picard_hs_metrics/CPG1.picard-hs-metrics'


def test_relevant_projects_drops_test():
    assert script.relevant_projects(['rdnow', 'rdnow-test', 'acute-care'], True) == ['rdnow-test']
    assert script.relevant_projects(['rdnow', 'rdnow-test', 'acute-care'], False) == ['acute-care', 'rdnow']


@pytest.fixture
def metamist(monkeypatch):
    """
    SGs: CPG1 has a complete QC set, CPG5 a partial one, CPG2 already has cram_qc, CPG3 has no CRAM,
    CPG4 has no QC files, CPG6 everything but the optional duplicate metrics, CPG7 only the duplicate metrics.
    """
    sgs = [
        {
            'id': 'CPG1',
            'analyses': [
                {'id': 1, 'type': 'cram', 'outputs': {'path': 'gs://cpg-ds-main/cram/old.cram'}},
                {'id': 2, 'type': 'cram', 'outputs': {'path': 'gs://cpg-ds-main/cram/CPG1.cram'}},
            ],
        },
        {
            'id': 'CPG2',
            'analyses': [
                {'id': 3, 'type': 'cram', 'outputs': 'gs://cpg-ds-main/cram/CPG2.cram'},
                {'id': 4, 'type': 'cram_qc', 'outputs': 'gs://cpg-ds-main/cram/CPG2.cram'},
            ],
        },
        {'id': 'CPG3', 'analyses': []},
        {'id': 'CPG4', 'analyses': [{'id': 5, 'type': 'cram', 'outputs': 'gs://cpg-ds-main/cram/CPG4.cram'}]},
        {'id': 'CPG5', 'analyses': [{'id': 6, 'type': 'cram', 'outputs': 'gs://cpg-ds-main/cram/CPG5.cram'}]},
        {'id': 'CPG6', 'analyses': [{'id': 7, 'type': 'cram', 'outputs': 'gs://cpg-ds-main/cram/CPG6.cram'}]},
        {'id': 'CPG7', 'analyses': [{'id': 8, 'type': 'cram', 'outputs': 'gs://cpg-ds-main/cram/CPG7.cram'}]},
    ]
    cpg6_qc = script.expected_qc_files('CPG6', 'ds', 'genome')
    del cpg6_qc['markduplicates_metrics']
    monkeypatch.setattr(script, 'query', lambda *_a, **_k: {'project': {'sequencingGroups': sgs}})
    present = {
        'gs://cpg-ds-main/cram/CPG1.cram',
        'gs://cpg-ds-main/cram/CPG4.cram',
        'gs://cpg-ds-main/cram/CPG5.cram',
        *script.expected_qc_files('CPG1', 'ds', 'genome').values(),
        'gs://cpg-ds-main/qc/samtools_stats/CPG5.samtools-stats',
        'gs://cpg-ds-main/cram/CPG6.cram',
        *cpg6_qc.values(),
        'gs://cpg-ds-main/cram/CPG7.cram',
        'gs://cpg-ds-main/qc/markduplicates_metrics/CPG7.markduplicates-metrics',
    }
    monkeypatch.setattr(script, 'exists', lambda path: path in present)


@pytest.mark.usefixtures('metamist')
def test_find_records_keeps_only_unregistered_sgs_with_qc():
    records = {r.sg_id: r for r in script.find_records('ds', 'genome')}
    # CPG7's duplicate metrics alone don't count: they come from alignment, not QC
    assert set(records) == {'CPG1', 'CPG5', 'CPG6'}
    complete = records['CPG1']
    assert complete.cram == 'gs://cpg-ds-main/cram/CPG1.cram'  # the latest CRAM analysis
    assert complete.secondary == script.expected_qc_files('CPG1', 'ds', 'genome')
    assert not complete.missing
    assert not complete.missing_optional
    no_markdup = records['CPG6']
    assert not no_markdup.missing
    assert no_markdup.missing_optional == ['markduplicates_metrics']
    partial = records['CPG5']
    assert partial.secondary == {'samtools_stats': 'gs://cpg-ds-main/qc/samtools_stats/CPG5.samtools-stats'}
    assert 'verifybamid' in partial.missing


def test_no_mosdepth_keys_expected():
    for seq_type in script.SEQUENCING_TYPES:
        assert not any('mosdepth' in key for key in script.expected_qc_files('CPG1', 'ds', seq_type))


@pytest.mark.usefixtures('metamist')
def test_dry_run_writes_file_and_registers_nothing(tmp_path, monkeypatch):
    def fail(**_kwargs: object) -> NoReturn:
        raise AssertionError('dry run must not register')

    monkeypatch.setattr(script, 'create_new', fail)
    out = tmp_path / 'dry.json'
    incomplete = script.main(['ds'], ['genome'], str(out))
    assert [r.sg_id for r in incomplete] == ['CPG5']
    written = json.loads(out.read_text())
    assert [r['sg_id'] for r in written['to_register']] == ['CPG1', 'CPG6']
    assert written['to_register'][1]['missing_optional'] == ['markduplicates_metrics']
    assert written['to_register'][0]['meta']['sequencing_type'] == 'genome'
    assert [r['sg_id'] for r in written['errors']] == ['CPG5']
    assert 'verifybamid' in written['errors'][0]['missing']


@pytest.mark.usefixtures('metamist')
def test_registers_only_complete_qc_sets(monkeypatch, caplog):
    calls = []
    monkeypatch.setattr(script, 'create_new', lambda **kwargs: calls.append(kwargs))
    incomplete = script.main(['ds'], ['genome'], None)
    assert [c['sgs'] for c in calls] == [['CPG1'], ['CPG6']]
    assert calls[0]['project'] == 'ds'
    assert calls[0]['output'] == 'gs://cpg-ds-main/cram/CPG1.cram'
    assert calls[0]['analysis_type'] == 'cram_qc'
    assert calls[0]['secondary'] == script.expected_qc_files('CPG1', 'ds', 'genome')
    assert [r.sg_id for r in incomplete] == ['CPG5']
    assert any(r.levelname == 'ERROR' and 'CPG5' in r.message for r in caplog.records)
    # registered without its duplicate metrics, which is flagged but not an error
    assert 'markduplicates_metrics' not in calls[1]['secondary']
    assert any(
        r.levelname == 'WARNING' and 'CPG6' in r.message and 'markduplicates_metrics' in r.message
        for r in caplog.records
    )
    assert not any(r.levelname == 'ERROR' and 'CPG6' in r.message for r in caplog.records)


def test_register_refuses_incomplete_record(monkeypatch):
    monkeypatch.setattr(script, 'create_new', lambda **_k: pytest.fail('must not register'))
    record = script.CramQcRecord('CPG5', 'ds', 'genome', 'gs://x.cram', missing=['verifybamid'])
    with pytest.raises(ValueError, match='incomplete'):
        script.register(record)
