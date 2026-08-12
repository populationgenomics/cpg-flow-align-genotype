"""Render a `[qc_thresholds.<seq_type>...]` block for pasting into config_template.toml.

This produces a string for a human to diff and paste; nothing writes the config. That
file carries comments, ordering and judgement no generator reproduces, and the paste-
after-diff step is the cheapest guard against a calibration run quietly rewriting a
production gate. Section order matches the committed file so the diff reads as a change
rather than a rewrite.
"""

from align_genotype.qc_calibration.settings import CalibrationSettings
from align_genotype.qc_calibration.thresholds import Candidate

# (severity, direction) in the order config_template.toml writes them.
_SECTIONS: tuple[tuple[str, str], ...] = (('fail', 'min'), ('fail', 'max'), ('warn', 'min'), ('warn', 'max'))


def _fmt(value: float) -> str:
    """Render a threshold as a TOML scalar.

    `hasattr(value, 'item')` unwraps a numpy scalar first: `np.float64` reprs as
    `np.float64(0.75)`, which is not valid TOML, and `np.int64` would fall through to
    `str()` and be written as a quoted string - valid TOML of the wrong type, which is
    the worse failure. Candidates are cast in `thresholds._round_for_unit`, but this is
    the single sink for everything written, so it defends here too.

    `bool` subclasses `int`, so it is rejected explicitly rather than rendered as `true`/
    `false`: every value that reaches this sink is a threshold, a `k`, or a `min_samples`,
    never a genuine boolean, so a bool here is always an upstream bug, not a legitimate
    value to coerce and paste into production config.
    """
    if hasattr(value, 'item'):
        value = value.item()
    if isinstance(value, bool):
        raise TypeError(f'threshold value must be numeric, got a bool: {value!r}')
    return repr(value)


def render(settings: CalibrationSettings, candidates: dict[str, Candidate]) -> str:
    """The config block for every metric that produced a candidate."""
    lines: list[str] = []
    for severity, direction in _SECTIONS:
        entries = [
            (metric.key, threshold)
            for metric in settings.metrics
            if metric.direction == direction
            and (found := candidates.get(metric.key)) is not None
            and (threshold := getattr(found, severity)) is not None
        ]
        if not entries:
            continue
        lines += ['', f'[qc_thresholds.{settings.seq_type}.{severity}.{direction}]']
        lines += [f'"{key}" = {_fmt(threshold)}' for key, threshold in entries]

    for metric in settings.relative_metrics:
        if metric.key not in candidates:
            continue
        lines += [
            '',
            f'[qc_thresholds.{settings.seq_type}.relative.{metric.key}]',
            f'direction = "{metric.direction}"',
            f'k = {_fmt(settings.k)}',
            f'min_samples = {_fmt(settings.min_samples)}',
        ]

    return '\n'.join(lines).lstrip('\n') + ('\n' if lines else '')
