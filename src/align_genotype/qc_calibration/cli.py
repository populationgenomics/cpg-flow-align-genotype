"""`qc_calibrate` - derive QC thresholds from a set of MultiQC reports.

Typical run:

    qc_calibrate discover      --seq-type genome --output calibration/manifest.genome.toml
    qc_calibrate collect       --spec calibration/spec.genome.toml --manifest calibration/manifest.genome.toml
    qc_calibrate distributions --spec calibration/spec.genome.toml
    qc_calibrate suggest       --spec calibration/spec.genome.toml
    qc_calibrate flagrates     --spec calibration/spec.genome.toml      # iterate here
    qc_calibrate mad           --spec calibration/spec.genome.toml
    qc_calibrate emit-config   --spec calibration/spec.genome.toml
    qc_calibrate dryrun        --spec ... --manifest ... --cohort <label>

A thin shell over the modules that do the work: every command parses its options, calls
one module, prints the string it hands back and picks an exit code. Nothing is rendered
here - `report` returns strings and this file echoes them - so the tables stay testable
against their exact content.

Two conventions worth knowing before editing:

- Command names are given explicitly wherever the function name differs from the
  command (`suggest_cmd`, `dryrun_cmd`, `emit_config`). Click derives a command's name
  from the function name, so `suggest_cmd` would register as `suggest-cmd`; the `_cmd`
  suffix exists only so the function doesn't shadow the module of the same name that
  this file imports. `emit.render`'s generated header names `qc_calibrate emit-config`
  in the block it writes into a committed file, so that one name in particular has to
  keep matching.
- Every expected failure becomes a `click.ClickException` via `@clean_errors`, so an
  operator gets one `Error: ...` line and exit 1 for a config problem rather than a
  stack trace.
"""

import functools
import logging
import sys
from collections.abc import Callable
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import click

from align_genotype.qc_calibration import cache as cache_mod
from align_genotype.qc_calibration import collect as collect_mod
from align_genotype.qc_calibration import discovery, dryrun, emit, relative, report, suggest
from align_genotype.qc_calibration import manifest as manifest_mod
from align_genotype.qc_calibration import spec as spec_mod

if sys.version_info >= (3, 11):
    import tomllib
else:  # pragma: no cover - 3.10 only
    import tomli as tomllib

logging.basicConfig(level=logging.INFO, format='%(levelname)s: %(message)s')

# Everything an operator can cause by handing this tool a bad file, a stale cache or a
# cohort set that doesn't match the spec. `TOMLDecodeError` is in here because a spec or
# manifest is hand-edited: it subclasses `ValueError` but none of the module error types,
# so a stray `=` in a spec would otherwise be the one config problem that still produced
# a traceback. `clean_errors` catches it in an earlier, more specific clause purely to
# name the file it came from; membership here is what guarantees it is caught at all.
CALIBRATION_ERRORS = (
    spec_mod.SpecError,
    manifest_mod.ManifestError,
    cache_mod.CacheError,
    collect_mod.CollectError,
    emit.EmitError,
    discovery.DiscoveryError,
    tomllib.TOMLDecodeError,
)

# The parameters naming a TOML file this tool *reads*. `TOMLDecodeError` reports a line
# and column but no filename, and `collect` and `dryrun` are each given two TOML paths,
# so the paths are named alongside it. Deliberately not every `*_path` parameter: an
# `--output` is written, never parsed, so citing it would send an operator to the wrong
# file.
_TOML_INPUT_PARAMS = ('spec_path', 'manifest_path')

_spec_option = click.option('--spec', 'spec_path', required=True, help='Path to the calibration spec TOML.')
_manifest_option = click.option('--manifest', 'manifest_path', required=True, help='Path to the cohort manifest TOML.')


def clean_errors(fn: Callable[..., Any]) -> Callable[..., Any]:
    """Turn expected calibration failures into a one-line error, not a traceback."""

    @functools.wraps(fn)
    def wrapper(*args: Any, **kwargs: Any) -> Any:
        try:
            return fn(*args, **kwargs)
        except tomllib.TOMLDecodeError as exc:
            paths = [str(kwargs[name]) for name in _TOML_INPUT_PARAMS if kwargs.get(name)]
            where = f' in {" or ".join(paths)}' if paths else ''
            raise click.ClickException(f'TOML syntax error{where} - {exc}') from exc
        except CALIBRATION_ERRORS as exc:
            raise click.ClickException(str(exc)) from exc

    return wrapper


def _load_spec_and_cache(spec_path: str) -> tuple[spec_mod.CalibrationSpec, cache_mod.ValueCache]:
    """The spec, plus the cache it points at - refused if it can't answer for that spec.

    Shared by every cache-consuming command so an incomplete or stale cache is rejected
    with the same message wherever it is met, rather than surfacing as an odd table or a
    `KeyError` several layers down.
    """
    spec = spec_mod.load(spec_path)
    values = cache_mod.load(spec.cache)
    cache_mod.require_usable(values, spec)
    return spec, values


def _now() -> str:
    return datetime.now(tz=timezone.utc).isoformat(timespec='seconds')


@click.group()
def main() -> None:
    """Derive warn/fail and cohort-relative QC thresholds from MultiQC reports."""


@main.command()
@click.option('--seq-type', required=True, help="Sequencing type, e.g. 'genome' or 'exome'.")
@click.option('--output', 'output_path', required=True, help='Where to write the manifest TOML.')
@clean_errors
def discover(seq_type: str, output_path: str) -> None:
    """Query Metamist for the latest MultiQC report per eligible dataset."""
    built = discovery.build_manifest(seq_type, generated=_now())
    manifest_mod.save(built, output_path)
    plural = '' if len(built.cohorts) == 1 else 's'
    click.echo(
        f'Wrote {len(built.cohorts)} cohort{plural} to {output_path}.\n'
        f'Review it before collecting: dropping a bad cohort, pinning an older analysis or substituting a local '
        f'file are all normal, and collect is the slow step.',
    )


@main.command()
@_spec_option
@_manifest_option
@clean_errors
def collect(spec_path: str, manifest_path: str) -> None:
    """Survey and extract every cohort - the one pass over the large reports."""
    spec = spec_mod.load(spec_path)
    result = collect_mod.collect_all(manifest_mod.load(manifest_path), spec, generated=_now())
    # Written before the exit check on purpose: a survey that describes a cache nobody
    # can open is useless, so an incomplete cache is still saved (marked incomplete, so
    # `require_usable` stops anything downstream running on it) and the survey printed.
    cache_mod.save(result.cache, spec.cache)
    click.echo(report.survey_report(result, spec))
    click.echo(f'\nCache written to {spec.cache} (complete={result.cache.complete}).')
    if not result.ok:
        raise click.ClickException('Collection incomplete - see the survey output above. Nothing downstream will run.')


@main.command()
@_spec_option
@clean_errors
def distributions(spec_path: str) -> None:
    """Percentile tables per metric, per cohort."""
    spec, values = _load_spec_and_cache(spec_path)
    click.echo(report.distributions_report(values, spec))


@main.command('suggest')
@_spec_option
@click.option('--output', 'output_path', default=None, help='Where to write the updated spec (default: in place).')
@clean_errors
def suggest_cmd(spec_path: str, output_path: str | None) -> None:
    """Seed candidate thresholds from the distributions, marked unreviewed."""
    spec, values = _load_spec_and_cache(spec_path)
    updated, seeded = suggest.seed(values, spec)
    spec_mod.save(updated, output_path or spec_path)
    click.echo(report.suggest_summary(seeded))


@main.command()
@_spec_option
@clean_errors
def flagrates(spec_path: str) -> None:
    """What the spec's candidate thresholds would flag, per cohort."""
    spec, values = _load_spec_and_cache(spec_path)
    click.echo(report.flagrates_report(values, spec))


@main.command()
@_spec_option
@click.option('--metric', 'metric_key', default=None, help='Evaluate one metric only.')
@clean_errors
def mad(spec_path: str, metric_key: str | None) -> None:
    """Cohort-relative (MAD) spread, warn rates and cohort-growth churn."""
    spec, values = _load_spec_and_cache(spec_path)
    if metric_key is not None and metric_key not in spec.metric_keys:
        raise click.ClickException(f'{metric_key!r} is not in {spec_path}; it defines {list(spec.metric_keys)}.')
    candidates = [m for m in spec.gated if m.relative is not None and metric_key in (None, m.key)]
    if not candidates:
        # A valid state, not an error: most metrics are gated absolutely, and a relative
        # tier is a candidate an operator proposes rather than something to be found.
        click.echo(
            'There are no metrics with a [metrics.<KEY>.relative] block to evaluate. To assess a candidate, add '
            'one (it needs an absolute `fail` gate behind it - cohort-relative flagging is warn-only) and re-run; '
            "the verdict then tells you whether the tier's warn rate and churn clear the adoption bar.",
        )
        return
    for metric in candidates:
        click.echo(report.mad_report(relative.evaluate(values, metric, spec.seq_type)))


@main.command('emit-config')
@_spec_option
@click.option('--manifest', 'manifest_path', default='', help='Manifest path, cited in the generated header.')
@click.option('--output', 'output_path', default=None, help='Write to a file instead of stdout.')
@clean_errors
def emit_config(spec_path: str, manifest_path: str, output_path: str | None) -> None:
    """Print the [qc_thresholds.<seq_type>...] block to paste into config_template.toml."""
    spec, values = _load_spec_and_cache(spec_path)
    block = emit.render(spec, values, spec_path=spec_path, manifest_path=manifest_path)
    if output_path:
        Path(output_path).write_text(block)
        click.echo(f'Wrote the config block to {output_path}.')
    else:
        click.echo(block)


@main.command('dryrun')
@_spec_option
@_manifest_option
@click.option('--cohort', 'cohort_label', required=True, help='Which cohort in the manifest to run against.')
@click.option('--output-dir', default='.', help='Where to write the structured QC flags JSON.')
@clean_errors
def dryrun_cmd(spec_path: str, manifest_path: str, cohort_label: str, output_dir: str) -> None:
    """Run the real check_multiqc against one cohort using the emitted config."""
    spec, values = _load_spec_and_cache(spec_path)
    loaded = manifest_mod.load(manifest_path)
    if cohort_label not in loaded.labels:
        raise click.ClickException(f'{cohort_label!r} is not in {manifest_path}; it lists {list(loaded.labels)}.')
    click.echo(report.dryrun_summary(dryrun.execute(spec, values, loaded.cohort(cohort_label), output_dir)))


if __name__ == '__main__':
    main()
