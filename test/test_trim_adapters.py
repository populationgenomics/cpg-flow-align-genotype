"""Tests for the trim_adapters repair script's realignment job.

The batch is built on a local backend and never run. Config, the reference resource group and the
intermediate-file existence checks are monkeypatched, so both FASTQ intermediates are "reused" and
only the DRAGMAP realign job is created. Paths are local because hailtop looks up any gs:// bucket.
"""

import shlex

import pytest

from hailtop.batch import Batch, LocalBackend

from align_genotype.scripts.repair_scripts import trim_adapters


@pytest.fixture
def align_job(monkeypatch, tmp_path):
    config_values = {
        ('images', 'dragmap'): 'dragmap:latest',
        ('references', 'dragmap_prefix'): str(tmp_path / 'dragmap'),
    }
    monkeypatch.setattr(
        trim_adapters.config,
        'config_retrieve',
        lambda key, default=None: config_values.get(tuple(key), default),
    )
    monkeypatch.setattr(trim_adapters.config, 'dataset_path', lambda suffix: str(tmp_path / suffix))
    monkeypatch.setattr(trim_adapters, 'exists', lambda _: True)
    monkeypatch.setattr(
        trim_adapters.hail_batch,
        'fasta_res_group',
        lambda b: b.read_input_group(
            base=str(tmp_path / 'genome.fasta'),
            fai=str(tmp_path / 'genome.fasta.fai'),
            dict=str(tmp_path / 'genome.dict'),
        ),
    )

    batch = Batch(backend=LocalBackend(tmp_dir=str(tmp_path)))
    [job] = trim_adapters.trim_and_realign(batch, str(tmp_path / 'cram' / 'CPG1234.cram'), 'CPG1234', {})
    return job


def test_dupblaster_spills_to_batch_tmpdir(align_job):
    [dupblaster_line] = [line for line in align_job._command[-1].splitlines() if 'dupblaster' in line]  # noqa: SLF001
    tokens = shlex.split(dupblaster_line)
    assert tokens[tokens.index('--tmp-dir') + 1] == '$BATCH_TMPDIR'
