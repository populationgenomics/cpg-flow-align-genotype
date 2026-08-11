"""End-to-end dry run of the production check against a small report."""

import json
import os
from collections.abc import Iterator
from pathlib import Path
from typing import Any, NoReturn

import pytest

from cpg_utils import config

from align_genotype.qc_calibration import dryrun as dryrun_mod
from align_genotype.qc_calibration import spec as spec_mod
from align_genotype.qc_calibration.cache import CohortValues, ValueCache
from align_genotype.qc_calibration.manifest import Cohort

SPEC = spec_mod.loads(
    'seq_type = "genome"\ncache = "c.json"\n'
    '\n[metrics.MEDIAN_COVERAGE]\ndirection = "min"\nunit = "x"\nfail = 15\nwarn = 25\nreviewed = true\n'
    'rationale = "Depth gate."\n'
    '\n[metrics.DUP]\ndirection = "max"\nunit = "%"\nfail = 40\nreviewed = true\nrationale = "Dup gate."\n'
    '[metrics.DUP.relative]\nk = 3.5\nmin_cohort = 5\n',
)

CACHE = ValueCache(
    seq_type='genome',
    generated='x',
    complete=True,
    metrics=('MEDIAN_COVERAGE', 'DUP'),
    cohorts=(CohortValues('dataset-a', 6, '1.33', 'dict', 0, {'MEDIAN_COVERAGE': [30.0] * 6, 'DUP': [5.0] * 6}),),
)

# One sample below the fail line, one in the warn band, four healthy; DUP has one
# clear within-cohort outlier so the relative pass has something to find.
SECTIONS = {
    'picard': {
        'CPG1|S1': {'MEDIAN_COVERAGE': 10.0},
        'CPG2|S2': {'MEDIAN_COVERAGE': 20.0},
        'CPG3|S3': {'MEDIAN_COVERAGE': 30.0},
        'CPG4|S4': {'MEDIAN_COVERAGE': 31.0},
        'CPG5|S5': {'MEDIAN_COVERAGE': 32.0},
        'CPG6|S6': {'MEDIAN_COVERAGE': 33.0},
    },
    'samtools': {
        'CPG1|S1': {'DUP': 5.0},
        'CPG2|S2': {'DUP': 5.5},
        'CPG3|S3': {'DUP': 6.0},
        'CPG4|S4': {'DUP': 6.5},
        'CPG5|S5': {'DUP': 7.0},
        'CPG6|S6': {'DUP': 35.0},
    },
}


@pytest.fixture(autouse=True)
def _restore_config_paths() -> Iterator[None]:
    previous = os.environ.get('CPG_CONFIG_PATH', '')
    yield
    config.set_config_paths([p for p in previous.split(',') if p])


@pytest.fixture
def report(tmp_path):
    path = tmp_path / 'multiqc_data.json'
    path.write_text(json.dumps({'config_version': '1.33', 'report_general_stats_data': SECTIONS}))
    return Cohort('dataset-a', str(path))


def test_build_config_includes_sequencing_type_and_thresholds():
    text = dryrun_mod.build_config(SPEC, CACHE)
    assert 'sequencing_type = "genome"' in text
    assert '[qc_thresholds.genome.fail.min]' in text


def test_build_config_is_valid_loadable_toml(tmp_path):
    """A dry run that passes should prove the emitted config actually loads."""
    path = tmp_path / 'c.toml'
    path.write_text(dryrun_mod.build_config(SPEC, CACHE))
    config.set_config_paths([str(path)])
    assert config.config_retrieve(['workflow', 'sequencing_type']) == 'genome'


def test_dryrun_reports_absolute_fail_and_warn(report, tmp_path):
    result = dryrun_mod.execute(SPEC, CACHE, report, output_dir=tmp_path)
    assert result.counts[('MEDIAN_COVERAGE', 'fail', 'absolute')] == 1  # the 10.0
    assert result.counts[('MEDIAN_COVERAGE', 'warn', 'absolute')] == 1  # the 20.0


def test_dryrun_exercises_the_relative_pass(report, tmp_path):
    result = dryrun_mod.execute(SPEC, CACHE, report, output_dir=tmp_path)
    assert result.counts[('DUP', 'warn', 'relative')] == 1  # the 35.0 outlier


def test_dryrun_counts_flagged_samples_and_writes_structured_output(report, tmp_path):
    result = dryrun_mod.execute(SPEC, CACHE, report, output_dir=tmp_path)
    assert result.n_samples_flagged == 3  # CPG1, CPG2, CPG6
    written = json.loads((tmp_path / 'dryrun_dataset-a.json').read_text())
    assert written['sequencing_type'] == 'genome'


def test_dryrun_records_timing_and_memory(report, tmp_path):
    result = dryrun_mod.execute(SPEC, CACHE, report, output_dir=tmp_path)
    assert result.seconds >= 0
    assert result.peak_rss_gb > 0


def test_dryrun_restores_config_paths(report, tmp_path):
    existing = tmp_path / 'existing.toml'
    existing.write_text('[workflow]\nsequencing_type = "exome"\n')
    config.set_config_paths([str(existing)])
    dryrun_mod.execute(SPEC, CACHE, report, output_dir=tmp_path)
    assert config.get_config_paths() == [str(existing)]


@pytest.mark.usefixtures('report')
def test_dryrun_restores_config_paths_on_failure(tmp_path):
    """A crash inside the check must not leave global config state pointing at a temp file."""
    existing = tmp_path / 'existing.toml'
    existing.write_text('[workflow]\nsequencing_type = "exome"\n')
    config.set_config_paths([str(existing)])
    missing = Cohort('dataset-a', str(tmp_path / 'nope.json'))
    with pytest.raises(Exception):  # noqa: B017 - any failure is fine; the point is the restore
        dryrun_mod.execute(SPEC, CACHE, missing, output_dir=tmp_path)
    assert config.get_config_paths() == [str(existing)]


def test_restore_failure_does_not_mask_the_real_error(report, tmp_path, monkeypatch):
    """If the previous config vanished, the body's exception must still be what propagates."""
    doomed = tmp_path / 'doomed.toml'
    doomed.write_text('[workflow]\nsequencing_type = "exome"\n')
    config.set_config_paths([str(doomed)])

    def _boom(*_args: Any, **_kwargs: Any) -> NoReturn:
        doomed.unlink()  # previous config disappears mid-run
        raise RuntimeError('THE REAL ERROR')

    monkeypatch.setattr(dryrun_mod.check_multiqc, 'run', _boom)
    with pytest.raises(RuntimeError, match='THE REAL ERROR'):
        dryrun_mod.execute(SPEC, CACHE, report, output_dir=tmp_path)


def test_dryrun_does_not_send_to_slack(report, tmp_path, monkeypatch):
    """Flipping send_to_slack in execute() must fail this test - the spy records any post.

    There's a second, independent barrier behind this one: under `_config_from` the
    installed config carries only `[workflow]` and `[qc_thresholds...]`, so
    `send_message` -> `_get_token()` raises `ValueError('slack.token_secret_id and
    slack.token_project_id must be set in config')` before any HTTP call is made. A
    flipped flag fails loudly rather than silently posting.
    """
    calls = []
    monkeypatch.setattr(dryrun_mod.check_multiqc, 'send_message', lambda *a, **_k: calls.append(a))
    dryrun_mod.execute(SPEC, CACHE, report, output_dir=tmp_path)
    assert calls == []


def test_dryrun_works_from_a_cold_start_with_no_config_paths(report, tmp_path):
    config.set_config_paths([])
    os.environ.pop('CPG_CONFIG_PATH', None)
    result = dryrun_mod.execute(SPEC, CACHE, report, output_dir=tmp_path)
    assert result.n_samples_flagged == 3


def test_dryrun_creates_a_missing_output_dir(report, tmp_path):
    """A full run must not be thrown away just because output_dir doesn't exist yet."""
    missing_dir = tmp_path / 'nested' / 'output'
    result = dryrun_mod.execute(SPEC, CACHE, report, output_dir=missing_dir)
    assert Path(result.output_path).exists()


def test_dryrun_removes_stale_output_on_failure(tmp_path):
    """A pre-existing output file must not survive a failed run to be misread as current."""
    output_path = tmp_path / 'dryrun_dataset-a.json'
    output_path.write_text('{"stale": true}')
    missing = Cohort('dataset-a', str(tmp_path / 'nope.json'))
    with pytest.raises(Exception):  # noqa: B017 - any failure is fine; the point is the cleanup
        dryrun_mod.execute(SPEC, CACHE, missing, output_dir=tmp_path)
    assert not output_path.exists()
