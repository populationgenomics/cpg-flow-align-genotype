"""CLI behaviour, especially exit codes on the failure paths."""

import json
import os
from collections.abc import Iterator

import pytest
from click.testing import CliRunner, Result

from cpg_utils import config

from align_genotype.qc_calibration import cache as cache_mod
from align_genotype.qc_calibration import spec as spec_mod
from align_genotype.qc_calibration.cli import main
from align_genotype.qc_calibration.manifest import Cohort, Manifest

# Eight samples: MEDIAN_COVERAGE spans the fail line, FREEMIX is uniformly healthy, and
# DUP carries one clear within-cohort outlier so a relative tier has something to find.
SECTIONS = {
    'picard': {
        f'CPG{i}|S{i}': {'MEDIAN_COVERAGE': float(10 + i), 'DUP': 35.0 if i == 7 else 5.0 + i * 0.5} for i in range(8)
    },
    'verifybamid': {f'CPG{i}|S{i}': {'FREEMIX': 0.001} for i in range(8)},
}

SPEC_TEXT = (
    'seq_type = "genome"\ncache = "{cache}"\n'
    '\n[metrics.MEDIAN_COVERAGE]\ndirection = "min"\nunit = "x"\nfail = 15\nwarn = 25\nreviewed = true\n'
    'rationale = "Depth gate."\n'
    '\n[metrics.FREEMIX]\ndirection = "max"\nunit = "frac"\nfail = 0.04\nreviewed = true\n'
    'rationale = "Contamination gate."\n'
)

# The same spec plus a cohort-relative warn tier, for the `mad` command.
RELATIVE_SPEC_TEXT = (
    SPEC_TEXT + '\n[metrics.DUP]\ndirection = "max"\nunit = "%"\nfail = 40\nreviewed = true\n'
    'rationale = "Duplication gate."\n[metrics.DUP.relative]\nk = 3.5\nmin_cohort = 5\n'
)


@pytest.fixture(autouse=True)
def _restore_config_paths() -> Iterator[None]:
    """`dryrun` and `mad` install a throwaway config; put the caller's back afterwards."""
    previous = os.environ.get('CPG_CONFIG_PATH', '')
    yield
    config.set_config_paths([p for p in previous.split(',') if p])


@pytest.fixture
def workspace(tmp_path):
    """A spec, a manifest and one small report, all under tmp_path."""
    report = tmp_path / 'dataset-a.json'
    report.write_text(json.dumps({'config_version': '1.33', 'report_general_stats_data': SECTIONS}))

    manifest = tmp_path / 'manifest.toml'
    manifest.write_text(f'seq_type = "genome"\ngenerated = "x"\n\n[cohorts.dataset-a]\nuri = "{report}"\n')

    spec = tmp_path / 'spec.toml'
    spec.write_text(SPEC_TEXT.format(cache=tmp_path / 'values.json'))
    return {'spec': str(spec), 'manifest': str(manifest), 'cache': str(tmp_path / 'values.json'), 'dir': tmp_path}


def _run(*args: str) -> Result:
    return CliRunner().invoke(main, list(args))


def _handled_cleanly(result: Result) -> bool:
    """Whether click itself reported the failure, rather than an exception escaping.

    `SystemExit` is what a `ClickException` leaves behind once click has printed
    `Error: <msg>`; anything else means the operator got a traceback.
    """
    return result.exception is None or isinstance(result.exception, SystemExit)


def test_help_lists_every_subcommand():
    result = _run('--help')
    assert result.exit_code == 0
    for command in ('discover', 'collect', 'distributions', 'suggest', 'flagrates', 'mad', 'emit-config', 'dryrun'):
        assert command in result.output


def test_help_does_not_register_the_cmd_suffixed_names():
    """`suggest_cmd` would register as `suggest-cmd` without an explicit command name."""
    result = _run('--help')
    assert 'suggest-cmd' not in result.output
    assert 'dryrun-cmd' not in result.output
    assert 'emit-config-cmd' not in result.output


def test_discover_writes_a_manifest_and_says_to_review_it(tmp_path, monkeypatch):
    built = Manifest('genome', 'x', (Cohort('dataset-a', 'gs://x/a.json'), Cohort('dataset-b', 'gs://x/b.json')))
    monkeypatch.setattr('align_genotype.qc_calibration.cli.discovery.build_manifest', lambda *a, **k: built)  # noqa: ARG005
    out = tmp_path / 'manifest.toml'
    result = _run('discover', '--seq-type', 'genome', '--output', str(out))
    assert result.exit_code == 0, result.output
    assert '2 cohorts' in result.output
    assert 'eview' in result.output  # 'Review it before collecting'
    assert '[cohorts.dataset-b]' in out.read_text()


def test_collect_writes_a_complete_cache(workspace):
    result = _run('collect', '--spec', workspace['spec'], '--manifest', workspace['manifest'])
    assert result.exit_code == 0, result.output
    assert cache_mod.load(workspace['cache']).complete is True


def test_collect_exits_nonzero_when_a_gated_metric_is_missing(workspace):
    spec = workspace['dir'] / 'spec.toml'
    spec.write_text(spec.read_text() + '\n[metrics.NOT_A_REAL_KEY]\ndirection = "min"\nfail = 1\n')
    result = _run('collect', '--spec', str(spec), '--manifest', workspace['manifest'])
    assert result.exit_code == 1
    assert 'NOT_A_REAL_KEY' in result.output
    # The cache is still written - marked incomplete - so the survey above describes
    # something an operator can actually inspect.
    assert cache_mod.load(workspace['cache']).complete is False


def test_downstream_command_refuses_an_incomplete_cache(workspace):
    spec = workspace['dir'] / 'spec.toml'
    spec.write_text(spec.read_text() + '\n[metrics.NOT_A_REAL_KEY]\ndirection = "min"\nfail = 1\n')
    _run('collect', '--spec', str(spec), '--manifest', workspace['manifest'])
    result = _run('flagrates', '--spec', str(spec))
    assert result.exit_code == 1
    assert 'incomplete' in result.output


def test_flagrates_reports_rates(workspace):
    _run('collect', '--spec', workspace['spec'], '--manifest', workspace['manifest'])
    result = _run('flagrates', '--spec', workspace['spec'])
    assert result.exit_code == 0, result.output
    assert 'MEDIAN_COVERAGE' in result.output
    assert 'fail%' in result.output


def test_distributions_reports_percentiles(workspace):
    _run('collect', '--spec', workspace['spec'], '--manifest', workspace['manifest'])
    result = _run('distributions', '--spec', workspace['spec'])
    assert result.exit_code == 0, result.output
    assert 'p50' in result.output


def test_suggest_rewrites_the_spec_in_place(workspace):
    _run('collect', '--spec', workspace['spec'], '--manifest', workspace['manifest'])
    spec_path = workspace['dir'] / 'unreviewed.toml'
    spec_path.write_text(SPEC_TEXT.format(cache=workspace['cache']).replace('reviewed = true', 'reviewed = false'))
    result = _run('suggest', '--spec', str(spec_path))
    assert result.exit_code == 0, result.output
    reloaded = spec_mod.load(spec_path)
    assert reloaded.metric('MEDIAN_COVERAGE').reviewed is False
    # Seeded from the cohort's own low tail (values 10..17), replacing the spec's 15.
    assert reloaded.metric('MEDIAN_COVERAGE').fail == 10
    assert 'reviewed = false' in result.output


def test_suggest_can_write_to_a_different_spec(workspace):
    _run('collect', '--spec', workspace['spec'], '--manifest', workspace['manifest'])
    unreviewed = workspace['dir'] / 'unreviewed.toml'
    unreviewed.write_text(SPEC_TEXT.format(cache=workspace['cache']).replace('reviewed = true', 'reviewed = false'))
    out = workspace['dir'] / 'seeded.toml'
    result = _run('suggest', '--spec', str(unreviewed), '--output', str(out))
    assert result.exit_code == 0, result.output
    assert spec_mod.load(out).metric('MEDIAN_COVERAGE').fail == 10
    assert spec_mod.load(unreviewed).metric('MEDIAN_COVERAGE').fail == 15  # left alone


def test_emit_config_prints_the_block(workspace):
    _run('collect', '--spec', workspace['spec'], '--manifest', workspace['manifest'])
    result = _run('emit-config', '--spec', workspace['spec'], '--manifest', workspace['manifest'])
    assert result.exit_code == 0, result.output
    assert '[qc_thresholds.genome.fail.min]' in result.output
    # The generated header names the command that produced the block; if the registered
    # subcommand name ever drifts from what `emit` hardcodes, this is the guard.
    assert 'qc_calibrate emit-config' in result.output


def test_emit_config_writes_to_a_file_when_asked(workspace):
    _run('collect', '--spec', workspace['spec'], '--manifest', workspace['manifest'])
    out = workspace['dir'] / 'block.toml'
    result = _run('emit-config', '--spec', workspace['spec'], '--output', str(out))
    assert result.exit_code == 0, result.output
    assert '[qc_thresholds.genome.fail.min]' in out.read_text()


def test_emit_config_refuses_unreviewed_metrics(workspace):
    _run('collect', '--spec', workspace['spec'], '--manifest', workspace['manifest'])
    spec_path = workspace['dir'] / 'unreviewed.toml'
    spec_path.write_text(SPEC_TEXT.format(cache=workspace['cache']).replace('reviewed = true', 'reviewed = false'))
    result = _run('emit-config', '--spec', str(spec_path))
    assert result.exit_code == 1
    assert 'unreviewed' in result.output
    assert _handled_cleanly(result)


def test_mad_says_so_when_no_metric_has_a_relative_block(workspace):
    _run('collect', '--spec', workspace['spec'], '--manifest', workspace['manifest'])
    result = _run('mad', '--spec', workspace['spec'])
    # A spec with no relative candidate is a valid state, not an error.
    assert result.exit_code == 0, result.output
    assert 'no metrics' in result.output.lower()
    assert 'relative' in result.output


def test_mad_evaluates_a_relative_metric(workspace):
    spec_path = workspace['dir'] / 'relative.toml'
    spec_path.write_text(RELATIVE_SPEC_TEXT.format(cache=workspace['cache']))
    _run('collect', '--spec', str(spec_path), '--manifest', workspace['manifest'])
    result = _run('mad', '--spec', str(spec_path), '--metric', 'DUP')
    assert result.exit_code == 0, result.output
    assert 'Cohort-relative (MAD) evaluation: DUP' in result.output
    assert 'Verdict:' in result.output


def test_mad_rejects_an_unknown_metric(workspace):
    _run('collect', '--spec', workspace['spec'], '--manifest', workspace['manifest'])
    result = _run('mad', '--spec', workspace['spec'], '--metric', 'NOPE')
    assert result.exit_code == 1
    assert 'NOPE' in result.output
    assert _handled_cleanly(result)


def test_dryrun_runs_the_real_check(workspace):
    _run('collect', '--spec', workspace['spec'], '--manifest', workspace['manifest'])
    result = _run(
        'dryrun',
        '--spec',
        workspace['spec'],
        '--manifest',
        workspace['manifest'],
        '--cohort',
        'dataset-a',
        '--output-dir',
        str(workspace['dir']),
    )
    assert result.exit_code == 0, result.output
    assert 'Dry run: dataset-a' in result.output
    assert 'MEDIAN_COVERAGE' in result.output
    assert (workspace['dir'] / 'dryrun_dataset-a.json').exists()


def test_dryrun_rejects_a_cohort_absent_from_the_manifest(workspace):
    _run('collect', '--spec', workspace['spec'], '--manifest', workspace['manifest'])
    result = _run(
        'dryrun',
        '--spec',
        workspace['spec'],
        '--manifest',
        workspace['manifest'],
        '--cohort',
        'dataset-b',
        '--output-dir',
        str(workspace['dir']),
    )
    assert result.exit_code == 1
    assert 'dataset-b' in result.output
    assert _handled_cleanly(result)


def test_a_malformed_spec_is_a_clean_error_not_a_traceback(workspace):
    bad = workspace['dir'] / 'bad.toml'
    bad.write_text('seq_type = "genome"\ncache = "c.json"\n[metrics.M]\ndirection = "sideways"\nfail = 1\n')
    result = _run('flagrates', '--spec', str(bad))
    assert result.exit_code == 1
    assert 'direction' in result.output
    assert 'Traceback' not in result.output
    assert _handled_cleanly(result)


def test_a_toml_syntax_error_is_a_clean_error_not_a_traceback(workspace):
    """`TOMLDecodeError` subclasses `ValueError`, not any of the module error types."""
    bad = workspace['dir'] / 'broken.toml'
    bad.write_text('seq_type = "genome"\ncache = = "c.json"\n')
    result = _run('flagrates', '--spec', str(bad))
    assert result.exit_code == 1
    assert 'TOML syntax error' in result.output
    assert str(bad) in result.output  # which file, since TOMLDecodeError doesn't say
    assert 'Traceback' not in result.output
    assert _handled_cleanly(result)


def test_a_toml_syntax_error_in_the_manifest_is_a_clean_error(workspace):
    bad = workspace['dir'] / 'broken-manifest.toml'
    bad.write_text('seq_type = "genome"\n[cohorts.dataset-a\nuri = "x"\n')
    result = _run('collect', '--spec', workspace['spec'], '--manifest', str(bad))
    assert result.exit_code == 1
    assert 'TOML syntax error' in result.output
    assert str(bad) in result.output
    assert 'Traceback' not in result.output
    assert _handled_cleanly(result)
