"""Render the `[qc_thresholds.<seq_type>...]` block for config_template.toml.

This module prints; it never writes the config. `config_template.toml` is hand-
maintained - it carries comments, ordering and judgement no generator can reproduce -
so the tool's job ends at handing an operator a block to diff and paste. That review
step is also the cheapest available guard against a calibration run quietly rewriting a
production gate.

Nothing here may name a cohort. The emitted block is pasted into a file that is
committed and pushed to a public repository, and cohort labels are real CPG dataset
names. The header therefore records the cohort *count* and cites the spec/manifest
*paths* (which live in a gitignored working directory) - never a label. Per-metric
evidence is aggregated across cohorts for the same reason: a range of medians says what
an operator needs without naming who contributed which end of it.

The review gate is the other hard rule. `suggest` seeds thresholds from percentiles and
marks them `reviewed = false`; `render` refuses to emit while any gated metric is still
unreviewed. That interlock is the only thing standing between a machine-generated
percentile and production config, so it fails closed.
"""

import textwrap
from datetime import datetime, timezone

import numpy as np

from align_genotype.qc_calibration import report, stats, tomlio
from align_genotype.qc_calibration.cache import ValueCache
from align_genotype.qc_calibration.spec import CalibrationSpec, MetricSpec, RelativeSpec

# The project line limit is 120 characters; comments carry a two-character '# ' prefix
# that `textwrap` does not count, so wrap the text itself two characters shorter.
COMMENT_PREFIX = '# '
COMMENT_WIDTH = 120 - len(COMMENT_PREFIX)

# (severity, direction) in the order config_template.toml writes them, so a diff between
# the emitted block and the committed one lines up rather than reading as a rewrite.
_ABSOLUTE_SECTIONS: tuple[tuple[str, str], ...] = (('fail', 'min'), ('fail', 'max'), ('warn', 'min'), ('warn', 'max'))

_RELATIVE_PREAMBLE = (
    'Cohort-relative (MAD / modified z-score) flagging - warn only. For metrics whose absolute level is strongly '
    'cohort-, kit- or library-prep dependent, a fixed warn line either floods every cohort or never fires; instead a '
    "sample is warned when it is a statistical outlier *within the current run's cohort*: "
    'mz = 0.6745 * (value - median) / MAD, so the threshold is median +/- k * MAD / 0.6745. Relative flags never fail '
    '(advisory only), and are skipped when the cohort has fewer than `min_cohort` samples or has zero MAD. The '
    'absolute fail gate above still applies as a hard stop; a sample caught there is not flagged twice.'
)


class EmitError(RuntimeError):
    """The config block cannot be emitted as things stand."""


def _cohort_count(n: int) -> str:
    """'1 cohort' / '3 cohorts' - the only cohort provenance this module may print."""
    return f'{n} cohort' if n == 1 else f'{n} cohorts'


def _comment(text: str) -> list[str]:
    """Wrap `text` into '# '-prefixed lines inside the project's line limit.

    Returns an empty list for empty text, so a metric with neither a rationale nor
    usable evidence simply gets no comment rather than a bare '#'.
    """
    return [f'{COMMENT_PREFIX}{line}' for line in textwrap.wrap(text, width=COMMENT_WIDTH)]


def _threshold(metric: MetricSpec, severity: str) -> float | None:
    """The absolute threshold `metric` defines at `severity`, if any."""
    return metric.fail if severity == 'fail' else metric.warn


def _require_reviewed(spec: CalibrationSpec) -> None:
    """Raise unless every gated metric has been signed off.

    Names all of them at once: an operator reviewing thresholds wants the full list to
    work through, not one metric per failed run. Un-gated metrics are surveyed rather
    than enforced, so their review state is irrelevant here.
    """
    unreviewed = sorted(m.key for m in spec.gated if not m.reviewed)
    if unreviewed:
        raise EmitError(
            f'refusing to emit {spec.seq_type} thresholds: unreviewed gated metric(s) {unreviewed}. '
            f'Check each one against `qc_calibrate flagrates`, then set `reviewed = true` in the spec. '
            f'A percentile-seeded threshold must not reach production config unlooked-at.',
        )


def _evidence(cache: ValueCache, metric: MetricSpec, severity: str) -> str:
    """The cohort median range and observed flag-rate range behind one threshold.

    Ranges, never per-cohort readings: this text is committed, so it must convey the
    spread without naming who sits at either end. Flag rates come from `stats.flag_rates`
    so they match what production would actually report - in particular the warn rate
    excludes samples already failing. Returns '' when no cohort has data for `metric`.
    """
    medians: list[float] = []
    rates: list[float] = []
    for label in cache.labels:
        values = cache.series(label, metric.key)
        if values.size == 0:
            continue
        medians.append(float(np.median(values)))
        fail_rate, warn_rate = stats.flag_rates(values, metric.direction, metric.fail, metric.warn)
        rates.append(fail_rate if severity == 'fail' else warn_rate)
    if not medians:
        return ''
    lo, hi = report.fmt_measure(min(medians), metric.unit), report.fmt_measure(max(medians), metric.unit)
    return (
        f'Cohort medians {lo}-{hi}; '
        f'{severity} rate {min(rates):.0%}-{max(rates):.0%} across {_cohort_count(len(medians))}.'
    )


def _header(spec: CalibrationSpec, cache: ValueCache, spec_path: str, manifest_path: str, generated: str) -> list[str]:
    """Provenance for the block: date, cohort count and artifact paths - no labels."""
    lines = _comment(
        f'{spec.seq_type.capitalize()} QC thresholds, calibrated from {_cohort_count(len(cache.cohorts))} '
        f'({generated}) by `qc_calibrate emit-config`. Review this against the block it replaces before pasting: '
        f'the tool prints and never edits this file.',
    )
    lines += _comment(
        'Cohort labels are omitted deliberately - this file is committed and pushed. The per-cohort detail behind '
        'these numbers lives in the calibration spec, manifest and value cache.',
    )
    if spec_path:
        lines += _comment(f'Spec: {spec_path}')
    if manifest_path:
        lines += _comment(f'Manifest: {manifest_path}')
    return lines


def _absolute_section(spec: CalibrationSpec, cache: ValueCache, severity: str, direction: str) -> list[str]:
    """One `[qc_thresholds.<seq_type>.<severity>.<direction>]` table, or nothing.

    A metric appears only in the tiers it actually defines, so a warn-only metric has no
    fail entry and a cohort-relative metric has no absolute warn entry.
    """
    metrics = [m for m in spec.gated if m.direction == direction and _threshold(m, severity) is not None]
    if not metrics:
        return []
    lines = ['', f'[qc_thresholds.{spec.seq_type}.{severity}.{direction}]']
    for metric in metrics:
        rationale_and_evidence = ' '.join(
            part for part in (metric.rationale, _evidence(cache, metric, severity)) if part
        )
        lines += _comment(rationale_and_evidence)
        lines.append(tomlio.fmt_kv(metric.key, _threshold(metric, severity), quote_key=True))
    return lines


def _relative_section(spec: CalibrationSpec) -> list[str]:
    """The `[qc_thresholds.<seq_type>.relative.<KEY>]` tables, or nothing."""
    relative: list[tuple[MetricSpec, RelativeSpec]] = [(m, m.relative) for m in spec.gated if m.relative is not None]
    if not relative:
        return []
    lines = ['', *_comment(_RELATIVE_PREAMBLE)]
    for index, (metric, settings) in enumerate(relative):
        if index:
            # Separate consecutive tables. The first sits directly under the shared
            # preamble, with no blank line, matching config_template.toml.
            lines.append('')
        # The metric key is written unquoted as part of a table header, so a key needing
        # quotes would silently reload as a further level of nesting. `spec` rejects
        # those on load; this is the same guard `spec.dumps` keeps for the same reason.
        tomlio.require_bare_key(metric.key, 'metric key')
        lines += [
            f'[qc_thresholds.{spec.seq_type}.relative.{metric.key}]',
            tomlio.fmt_kv('direction', metric.direction),
            tomlio.fmt_kv('k', settings.k),
            tomlio.fmt_kv('min_cohort', settings.min_cohort),
        ]
    return lines


def render(
    spec: CalibrationSpec,
    cache: ValueCache,
    spec_path: str = '',
    manifest_path: str = '',
    generated: str = '',
) -> str:
    """Render the config block for `spec`, with generated evidence drawn from `cache`.

    Raises `EmitError` if any gated metric is still unreviewed - nothing is rendered in
    that case, so a half-reviewed spec cannot produce a pasteable block. `spec_path` and
    `manifest_path` are cited in the header when supplied; `generated` defaults to
    today's date.
    """
    _require_reviewed(spec)
    stamp = generated or datetime.now(tz=timezone.utc).date().isoformat()
    lines = _header(spec, cache, spec_path, manifest_path, stamp)
    for severity, direction in _ABSOLUTE_SECTIONS:
        lines += _absolute_section(spec, cache, severity, direction)
    lines += _relative_section(spec)
    return '\n'.join(lines) + '\n'
