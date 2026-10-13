"""Tests for the alignment_recovery repair script.

The batch is built on a local backend and never run; batch creation, the reference resource group,
and config access are monkeypatched, so no cpg-utils config file or Hail credentials are needed.
All paths are local because hailtop looks up any gs:// bucket it reads from or writes to.
"""

import sys

import pytest

from hailtop.batch import Batch, LocalBackend

from align_genotype.scripts.repair_scripts import alignment_recovery

CONFIG = {('workflow', 'driver_image'): 'driver:latest'}


@pytest.fixture
def batch(monkeypatch, tmp_path):
    local_batch = Batch(backend=LocalBackend(tmp_dir=str(tmp_path)))
    monkeypatch.setattr(local_batch, 'run', lambda **_: None)
    monkeypatch.setattr(alignment_recovery.hail_batch, 'get_batch', lambda *_, **__: local_batch)
    monkeypatch.setattr(
        alignment_recovery.hail_batch,
        'fasta_res_group',
        lambda b: b.read_input_group(
            base=str(tmp_path / 'genome.fasta'),
            fai=str(tmp_path / 'genome.fasta.fai'),
            dict=str(tmp_path / 'genome.dict'),
        ),
    )
    monkeypatch.setattr(
        alignment_recovery.config,
        'config_retrieve',
        lambda key, default=None: CONFIG.get(tuple(key), default),
    )
    return local_batch


@pytest.fixture
def run_recovery(batch, monkeypatch, tmp_path):
    """Run the script's main() with the given extra CLI args and return the merge job."""

    def _run(*extra_args: str):  # noqa: ANN202
        monkeypatch.setattr(
            sys,
            'argv',
            [
                'alignment_recovery',
                *('-s', 'CPG1234', '-p', 'dataset', '--no-register'),
                *(
                    '-i',
                    str(tmp_path / 'Align_1_2_-abcde' / 'sorted_bam'),
                    str(tmp_path / 'Align_2_2_-fghij' / 'sorted_bam'),
                ),
                *extra_args,
            ],
        )
        alignment_recovery.main()
        [merge_job] = [job for job in batch._jobs if job.name == 'Merge BAMs for CPG1234']  # noqa: SLF001
        return merge_job

    return _run


def test_metrics_written_where_the_pipeline_writes_them(run_recovery, tmp_path):
    merge_job = run_recovery('-o', str(tmp_path / 'cram' / 'CPG1234.cram'))

    assert merge_job._command[-1].startswith('mv $BATCH_TMPDIR/metrics.duplicate-metrics.tsv ')  # noqa: SLF001
    assert merge_job.markdup_metrics._output_paths == {  # noqa: SLF001
        str(tmp_path / 'qc' / 'markduplicates_metrics' / 'CPG1234.markduplicates-metrics'),
    }


def test_metrics_out_overrides_the_derived_path(run_recovery, tmp_path):
    metrics_out = str(tmp_path / 'elsewhere' / 'CPG1234.markduplicates-metrics')
    merge_job = run_recovery('-o', str(tmp_path / 'cram' / 'CPG1234.cram'), '--metrics-out', metrics_out)

    assert merge_job.markdup_metrics._output_paths == {metrics_out}  # noqa: SLF001


def test_metrics_out_required_when_cram_is_outside_a_cram_dir(run_recovery, tmp_path):
    with pytest.raises(SystemExit):
        run_recovery('-o', str(tmp_path / 'somewhere' / 'CPG1234.cram'))
