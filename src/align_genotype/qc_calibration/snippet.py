"""Render a `[qc_thresholds.<seq_type>...]` block for pasting into config_template.toml.

This produces a string for a human to diff and paste; nothing writes the config. That
file carries comments, ordering and judgement no generator reproduces, and the paste-
after-diff step is the cheapest guard against a calibration run quietly rewriting a
production gate. Section order matches the committed file so the diff reads as a change
rather than a rewrite.

Each threshold line is preceded by its `Candidate.basis` as a `# `-prefixed comment: a
pasted number with no stated evidence is exactly what calibration is meant to move away
from. The same metric's basis is repeated in both its `fail` and `warn` sections (and
again in its relative table, if any) rather than written once - the committed file
repeats context per section too, and a reader scanning `warn.min` should not have to
scroll to `fail.min` to find out where the number came from.
"""

import textwrap

from align_genotype.qc_calibration.settings import CalibrationSettings, MetricSpec
from align_genotype.qc_calibration.thresholds import Candidate

# Static glosses for the relative table's algorithm-level keys - what `k` and
# `min_samples` *mean*, not the per-metric statistical evidence `Candidate.basis`
# carries. Paving over these with only a basis comment would silently destroy
# documentation the committed file carries for every relative table.
_K_GLOSS = '# standard Iglewicz-Hoaglin outlier threshold'
_MIN_SAMPLES_GLOSS = '# below this, MAD is too noisy; skip relative flagging'

# (severity, direction) in the order config_template.toml writes them.
_SECTIONS: tuple[tuple[str, str], ...] = (('fail', 'min'), ('fail', 'max'), ('warn', 'min'), ('warn', 'max'))

# The project line limit is 120 characters; comments carry a two-character '# ' prefix
# that `textwrap` does not count, so wrap the text itself two characters shorter.
_LINE_LIMIT = 120
_COMMENT_PREFIX = '# '
_COMMENT_WIDTH = _LINE_LIMIT - len(_COMMENT_PREFIX)


def _comment(text: str) -> list[str]:
    """Wrap `text` into '# '-prefixed lines, or [] for empty text.

    Empty text renders no line at all rather than a bare '#': a candidate with no basis
    (or a relative metric whose basis is absent) must not leave a dangling comment marker
    with nothing after it.
    """
    if not text:
        return []
    return [f'{_COMMENT_PREFIX}{line}' for line in textwrap.wrap(text, width=_COMMENT_WIDTH)]


def _fmt(value: float, field: str) -> str:
    """Render a threshold as a TOML scalar. `field` names what's being rendered, for errors.

    `hasattr(value, 'item')` unwraps a numpy scalar first: `np.float64` reprs as
    `np.float64(0.75)`, which is not valid TOML, and `np.int64` would fall through to
    `str()` and be written as a quoted string - valid TOML of the wrong type, which is
    the worse failure. Candidates are cast in `thresholds._round_for_unit`, but this is
    the single sink for everything written, so it defends here too.

    `bool` subclasses `int`, so it is rejected explicitly rather than rendered as `true`/
    `false`: every value that reaches this sink is a threshold, a `k`, or a `min_samples`,
    never a genuine boolean, so a bool here is always an upstream bug, not a legitimate
    value to coerce and paste into production config. `field` (e.g.
    `qc_thresholds.genome.fail.max.FREEMIX`) says which one, since this only fires after
    the whole calibration run has already completed and `render` calls this in a loop -
    without it, whoever hits the error is bisecting `candidates` by hand.
    """
    if hasattr(value, 'item'):
        value = value.item()
    if isinstance(value, bool):
        raise TypeError(f'{field}: threshold value must be numeric, got a bool: {value!r}')
    return repr(value)


def _padded(base: str, comment: str, column: int) -> str:
    """`base` padded so `comment` starts at `column`, unless that would exceed the line limit."""
    line = f'{base.ljust(column)}{comment}'
    return line if len(line) <= _LINE_LIMIT else f'{base} {comment}'


def _relative_lines(metric: MetricSpec, settings: CalibrationSettings) -> list[str]:
    """The three annotated `key = value` lines a relative table needs, aligned like the committed file.

    `direction`'s annotation is derived from the value (`high`/`low`) rather than the
    committed file's metric-specific gloss ("more duplicate reads than the dataset norm")
    - a generic, accurate marker beats a generated guess at domain prose. `k` and
    `min_samples` get static boilerplate explaining what the constant means, which
    `Candidate.basis` - per-metric statistical evidence - does not and cannot substitute
    for; without it, a wholesale paste would silently delete that documentation.
    """
    field = f'qc_thresholds.{settings.seq_type}.relative.{metric.key}'
    bad = 'high' if metric.direction == 'max' else 'low'
    kvs = [
        ('direction', f'"{metric.direction}"', f'# bad = {bad}'),
        ('k', _fmt(settings.k, f'{field}.k'), _K_GLOSS),
        ('min_samples', _fmt(settings.min_samples, f'{field}.min_samples'), _MIN_SAMPLES_GLOSS),
    ]
    bases = [f'{key} = {value}' for key, value, _ in kvs]
    column = max(len(base) for base in bases) + 2
    return [_padded(base, comment, column) for base, (_, _, comment) in zip(bases, kvs, strict=True)]


def render(settings: CalibrationSettings, candidates: dict[str, Candidate]) -> str:
    """The config block for every metric that produced a candidate."""
    lines: list[str] = []
    for severity, direction in _SECTIONS:
        entries = [
            (metric.key, threshold, found.basis)
            for metric in settings.metrics
            if metric.direction == direction
            and (found := candidates.get(metric.key)) is not None
            and (threshold := getattr(found, severity)) is not None
        ]
        if not entries:
            continue
        lines += ['', f'[qc_thresholds.{settings.seq_type}.{severity}.{direction}]']
        for key, threshold, basis in entries:
            lines += _comment(basis)
            field = f'qc_thresholds.{settings.seq_type}.{severity}.{direction}.{key}'
            lines.append(f'"{key}" = {_fmt(threshold, field)}')

    for metric in settings.relative_metrics:
        found = candidates.get(metric.key)
        if found is None:
            continue
        lines += ['', f'[qc_thresholds.{settings.seq_type}.relative.{metric.key}]']
        lines += _comment(found.basis)
        lines += _relative_lines(metric, settings)

    return '\n'.join(lines).lstrip('\n') + ('\n' if lines else '')
