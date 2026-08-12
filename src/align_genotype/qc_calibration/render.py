"""Render the calibration summary to HTML.

The template directory and `autoescape=True` match `scripts/sg_qc_report.py`: dataset
names and skip reasons reach the page unmodified, so escaping is not optional.

The environment uses `StrictUndefined`: a typo in the template (`m.n_value` instead of
`m.n_values`) renders as an empty cell under the default `Undefined`, which is the worst
failure mode for a page people will use to decide production thresholds - it looks like
data, not an error. `StrictUndefined` turns that into a `jinja2.UndefinedError` at render
time. Fields that are legitimately optional (`m.current.fail` when nothing is shipped yet,
`m.candidate` when a metric has no data) are read through the `default` filter, which
still short-circuits on `Undefined` even in strict mode - only *unintended* gaps raise.
"""

from pathlib import Path
from typing import Any

import jinja2

JINJA_TEMPLATE_DIR = Path(__file__).absolute().parent.parent / 'templates'
TEMPLATE_NAME = 'qc_calibration_report.html.jinja'


def _pct(value: float | None) -> str:
    """A rate as a percentage, or an em dash when there was nothing to score."""
    return '—' if value is None else f'{value:.1%}'


def _num(value: float | None, unit: str = 'frac') -> str:
    """A metric value at a precision that reads sensibly for its unit."""
    if value is None:
        return '—'
    if unit == 'frac':
        return f'{value:.3f}'
    if unit in ('x', '%'):
        return f'{value:.1f}'
    return f'{value:.2f}'


def _stat(value: float | None) -> str:
    """A median, MAD or threshold at 4 dp, so the threshold explains the warn count."""
    return '—' if value is None else f'{value:.4f}'


def render(built: dict[str, Any]) -> str:
    """Render the dashboard for one assembled summary."""
    env = jinja2.Environment(
        loader=jinja2.FileSystemLoader(JINJA_TEMPLATE_DIR),
        autoescape=True,
        undefined=jinja2.StrictUndefined,
    )
    env.filters['pct'] = _pct
    env.filters['num'] = _num
    env.filters['stat'] = _stat
    return env.get_template(TEMPLATE_NAME).render(
        r=built,
        dataset_names=[d['dataset'] for d in built['datasets']],
        percentile_labels=[f'p{p}' for p in (1, 5, 10, 25, 50, 75, 90, 95, 99)],
    )
