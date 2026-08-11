"""End-to-end dry run of the production QC check against one real MultiQC report.

The old workflow had `dryrun_exome_check.py` and `dryrun_genome_check.py`: two
near-identical scripts kept in sync by hand, one per `seq_type`. Here the seq_type
comes from the calibration spec, so there is one command for both.

The old scripts also monkeypatched `config.config_retrieve`, so a passing dry run
proved nothing about whether the thresholds it used would actually load in
production - a monkeypatch answers whatever it's told to, regardless of whether the
real config shape would parse. Here the emitted `[qc_thresholds...]` block is
written to a temp TOML via `emit.render` and installed with
`cpg_utils.config.set_config_paths`, exactly as production config is installed. A
passing dry run therefore also proves the emitted block is what `load_thresholds`
(via `check_multiqc.run`) actually reads - the one thing monkeypatching could never
tell you.
"""

import logging
import sys
import time
from collections import Counter
from collections.abc import Iterator
from contextlib import contextmanager
from dataclasses import dataclass
from pathlib import Path
from tempfile import TemporaryDirectory
from typing import Any

from cpg_utils import config

from align_genotype.qc_calibration.cache import ValueCache
from align_genotype.qc_calibration.emit import render
from align_genotype.qc_calibration.manifest import Cohort
from align_genotype.qc_calibration.spec import CalibrationSpec
from align_genotype.scripts import check_multiqc


@dataclass(frozen=True)
class DryRunResult:
    """One cohort's dry-run outcome: what the production check flagged, and how."""

    cohort: str
    n_samples_flagged: int
    counts: Counter
    seconds: float
    peak_rss_gb: float
    output_path: str


def _peak_rss_gb() -> float:
    """Peak resident set size of this process, in GB.

    `resource.getrusage(...).ru_maxrss` is bytes on macOS/Darwin but KiB everywhere
    else (Linux), per the platform's `getrusage(2)`. Branching on `sys.platform` is
    the only way to interpret the same field correctly on both.
    """
    import resource  # noqa: PLC0415 - posix-only, and only needed here

    ru_maxrss = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    bytes_per_unit = 1 if sys.platform == 'darwin' else 1024
    return ru_maxrss * bytes_per_unit / 1e9


def build_config(spec: CalibrationSpec, cache: ValueCache) -> str:
    """The full config text a dry run installs: `[workflow]` plus the emitted thresholds.

    `check_multiqc.run` reads `['workflow', 'sequencing_type']` to pick which
    `qc_thresholds` table applies, so that key has to be present alongside the
    thresholds themselves for the check to do anything.
    """
    lines = ['[workflow]', f'sequencing_type = "{spec.seq_type}"', '']
    return '\n'.join(lines) + render(spec, cache, generated='dryrun')


@contextmanager
def _config_from(text: str) -> Iterator[None]:
    """Install `text` as the only config source for the duration of the block.

    Written as a real file and installed with `set_config_paths` rather than
    monkeypatching `config_retrieve`, so the dry run exercises the genuine config ->
    `load_thresholds` -> `check_multiqc.run` path (see module docstring).

    `get_config_paths` *raises* when nothing has ever been set, so a cold start is
    treated as "no previous paths" and restored to that. The restore itself is
    guarded so a failure there can't mask an exception raised by the body - see
    `_restore_config_paths`.
    """
    try:
        previous = list(config.get_config_paths())
    except config.ConfigError:
        previous = []
    with TemporaryDirectory(prefix='qc-calibration-dryrun-') as tmpdir:
        path = Path(tmpdir) / 'dryrun.toml'
        path.write_text(text)
        config.set_config_paths([str(path)])
        try:
            yield
        finally:
            _restore_config_paths(previous)


def _restore_config_paths(previous: list[str]) -> None:
    """Put the previous config paths back, without letting that failure win.

    `set_config_paths` re-validates every path it's given: it opens and parses each
    one. If a previous path became unreadable while the dry run was running, that
    failure would otherwise raise a `ValueError` naming the wrong problem *after* the
    body's own exception (or result) had already been produced. A failed restore is
    logged and the current (temp) paths are left in place instead, so the caller's
    real exception is what propagates.
    """
    try:
        config.set_config_paths(previous)
    except ValueError as exc:
        logging.warning(f'Could not restore previous config paths {previous}, leaving as-is - {exc}')


def execute(spec: CalibrationSpec, cache: ValueCache, cohort: Cohort, output_dir: Path) -> DryRunResult:
    """Run the production QC check against `cohort`'s report, under `spec`'s thresholds.

    `send_to_slack=False` always: a calibration dry run must never post to the lab's
    Slack channel. Structured output is written to
    `<output_dir>/dryrun_<cohort.label>.json`, exactly as a production run would
    write it with `--output-json`.
    """
    output_path = Path(output_dir) / f'dryrun_{cohort.label}.json'
    text = build_config(spec, cache)
    start = time.perf_counter()
    with _config_from(text):
        result: dict[str, Any] = check_multiqc.run(
            multiqc_json_path=cohort.uri,
            html_url=None,
            dataset=cohort.label,
            title=f'{spec.seq_type} calibration dry run',
            send_to_slack=False,
            output_json_path=str(output_path),
        )
    seconds = time.perf_counter() - start

    counts: Counter = Counter()
    for flags in result['qc_flags'].values():
        for flag in flags:
            counts[(flag['flag'], flag['severity'], flag.get('method', 'absolute'))] += 1

    return DryRunResult(
        cohort=cohort.label,
        n_samples_flagged=result['n_samples_flagged'],
        counts=counts,
        seconds=seconds,
        peak_rss_gb=_peak_rss_gb(),
        output_path=str(output_path),
    )
