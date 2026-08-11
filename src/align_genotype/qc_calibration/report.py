"""Human-readable formatting for the calibration commands.

Every renderer here *returns* a string rather than printing one. That keeps the tables
testable against their exact content - the survey report in particular is the thing an
operator reads to decide whether a gate is inert, so "does it actually say MISSING" is a
property worth asserting - and it keeps the CLI a thin shell that prints what it is
handed.

Nothing in this module imports the result types it renders at runtime. `emit` imports
`fmt_measure` from here and `dryrun` imports `emit`, so a runtime import of
`DryRunResult` would close the cycle `report -> dryrun -> emit -> report`. The result
types are therefore imported under `TYPE_CHECKING`, with `from __future__ import
annotations` in force so the annotations that mention them are never evaluated. Only
`stats` is imported for real, for the percentile list and the flag-rate guardrails the
tables are built from; it has no path back to here.
"""

from __future__ import annotations

import math
import textwrap
from typing import TYPE_CHECKING

from align_genotype.qc_calibration import stats

if TYPE_CHECKING:
    from collections.abc import Sequence

    from align_genotype.qc_calibration.cache import ValueCache
    from align_genotype.qc_calibration.collect import CollectResult
    from align_genotype.qc_calibration.dryrun import DryRunResult
    from align_genotype.qc_calibration.relative import MadEvaluation
    from align_genotype.qc_calibration.spec import CalibrationSpec, MetricSpec
    from align_genotype.qc_calibration.stats import ChurnResult
    from align_genotype.qc_calibration.suggest import Seeded

# Two spaces between columns: wide enough to read as a gap, narrow enough that a
# ten-cohort matrix still fits a wide terminal.
_GAP = '  '

# Prefix for a block an operator must not scroll past. Deliberately not colour: this
# output is as likely to be read in a CI log or a pasted paragraph as in a terminal.
_LOUD = '!! '

# Prose wrap width. Under the 120-column source limit so a wrapped note still has room
# for its indent, and narrow enough to stay readable next to the tables.
_WIDTH = 110

# Widest a table is allowed to get before its cohort columns are split into groups. The
# tables that grow with cohort *count* rather than cohort size are the ones that need
# this: the metric-presence matrix reaches ~250 columns for ten real dataset labels.
MAX_TABLE_WIDTH = 120

# Merged-cohort churn runs every unordered pair, which is quadratic in cohort count -
# 45 rows for ten cohorts. Only the worst few decide anything, so only those are shown,
# and the total is reported alongside so nothing looks hidden.
MAX_PAIRS_SHOWN = 5


def fmt_measure(value: float, unit: str) -> str:
    """Render a metric value at a precision that reads sensibly for its unit."""
    if unit == 'frac':
        return f'{value:.3f}'
    if unit in ('x', '%'):
        return f'{value:.1f}'
    return f'{value:.2f}'


def table(headers: Sequence[str], rows: Sequence[Sequence[str]]) -> str:
    """Fixed-width columns with a rule under the header.

    Every column is as wide as its widest cell, so values line up down the page and a
    missing entry is visible as a gap rather than shifting the row. Trailing padding is
    stripped: an empty last column (a marker column, say) must not leave invisible
    whitespace in output that gets pasted into a review.

    A row of the wrong width raises rather than rendering short - a silently dropped
    cohort column is exactly the kind of thing this report exists to make visible.
    """
    for row in rows:
        if len(row) != len(headers):
            raise ValueError(f'table: row {list(row)!r} has {len(row)} cells but there are {len(headers)} headers')
    widths = [max(len(h), max((len(row[i]) for row in rows), default=0)) for i, h in enumerate(headers)]
    lines = [_row(headers, widths), _GAP.join('-' * w for w in widths)]
    lines += [_row(row, widths) for row in rows]
    return '\n'.join(lines)


def _row(cells: Sequence[str], widths: Sequence[int]) -> str:
    return _GAP.join(c.ljust(w) for c, w in zip(cells, widths, strict=True)).rstrip()


def _titled(text: str, char: str = '=') -> str:
    return f'{text}\n{char * len(text)}'


def _join(blocks: Sequence[str]) -> str:
    return '\n\n'.join(b for b in blocks if b) + '\n'


def _sense(direction: str) -> str:
    """'min' means a low value is the problem, so higher is better."""
    return 'higher=better' if direction == 'min' else 'lower=better'


def _fmt_threshold(value: float | None, unit: str) -> str:
    return '-' if value is None else fmt_measure(value, unit)


def _fmt_rate(rate: float) -> str:
    """A flag rate as a percentage; '-' for the nan `flag_rates` returns on no samples."""
    return '-' if math.isnan(rate) else f'{rate:.1%}'


def _fmt_stat(value: float | None) -> str:
    """A median/MAD/threshold at 4 significant figures; '-' when there isn't one."""
    if value is None or math.isnan(value):
        return '-'
    return f'{value:.4g}'


def _wrap(text: str, initial: str = '', subsequent: str = '') -> list[str]:
    """Wrap prose to `_WIDTH`, never splitting a metric key across two lines."""
    return textwrap.wrap(
        text,
        width=_WIDTH,
        initial_indent=initial,
        subsequent_indent=subsequent,
        break_long_words=False,
        break_on_hyphens=False,
    )


def _fill(text: str) -> str:
    return '\n'.join(_wrap(text))


# ------------------------------------------------------------------------------ survey


def survey_report(result: CollectResult, spec: CalibrationSpec) -> str:
    """What each report looked like, and where (or whether) each metric was found.

    A cohort that failed to parse and a cohort missing a metric get separate blocks on
    purpose: the first is a URI, file or credentials problem and the cohort appears
    nowhere in the tables above, while the second is a spec-versus-MultiQC naming problem
    in a report that read perfectly well. Collapsing them into one "problems" list would
    send an operator to the wrong place half the time.
    """
    blocks = [_survey_cohorts(result), _survey_presence(result, spec)]
    if result.missing_gated:
        blocks.append(_survey_missing(result))
    if result.failures:
        blocks.append(_survey_failures(result))
    if result.ok:
        blocks.append('Survey OK: every gated metric was found in every surveyed cohort.')
    return _join(blocks)


def _survey_cohorts(result: CollectResult) -> str:
    headers = ['cohort', 'multiqc', 'shape', 'samples', 'sections', 'dropped']
    rows = [
        [
            row.label,
            row.multiqc_version,
            row.shape,
            str(row.n_samples),
            ', '.join(f'{name}={size}' for name, size in row.section_sizes.items()) or '-',
            str(row.n_dropped),
        ]
        for row in result.rows
    ]
    lines = [_titled('Cohorts surveyed'), table(headers, rows)]
    detail = [
        f'  {row.label}: ' + ', '.join(f'{key}={n}' for key, n in row.n_dropped_by_metric.items() if n)
        for row in result.rows
        if row.n_dropped
    ]
    if detail:
        # Per metric, not just the total: a metric where most samples are Picard's '?'
        # placeholder is a real signal, and one total can't be unpicked back into it.
        lines += ['', 'Dropped values by metric (non-numeric or non-finite):', *detail]
    return '\n'.join(lines)


def _column_width(header: str, cells: Sequence[str]) -> int:
    return max(len(header), max((len(cell) for cell in cells), default=0))


def _label_groups(labels: Sequence[str], columns: dict[str, list[str]], lead_width: int) -> list[list[str]]:
    """Greedily pack cohort columns into groups no wider than `MAX_TABLE_WIDTH`.

    A group can still overflow when one cohort's own column is wider than the budget -
    there is nowhere left to split - but that is one column too wide rather than ten.
    """
    groups: list[list[str]] = []
    current: list[str] = []
    width = lead_width
    for label in labels:
        needed = _column_width(label, columns[label]) + len(_GAP)
        if current and width + needed > MAX_TABLE_WIDTH:
            groups.append(current)
            current, width = [], lead_width
        current.append(label)
        width += needed
    groups.append(current)
    return groups


def _survey_presence(result: CollectResult, spec: CalibrationSpec) -> str:
    """The metric x cohort matrix, in column groups narrow enough to read.

    Ten cohorts with real dataset labels put a single matrix past 240 columns, which
    wraps into an unreadable block in any terminal. Splitting it into groups keeps every
    cell - truncating a section name or dropping a cohort column is exactly the silent
    loss this survey exists to prevent - and repeats the metric and gated columns in each
    group so every group reads on its own.
    """
    labels = [row.label for row in result.rows]
    columns = {
        row.label: [
            ','.join(where) if (where := row.where.get(metric.key, ())) else 'MISSING' for metric in spec.metrics
        ]
        for row in result.rows
    }
    lead_headers = ['metric', 'gated']
    lead_rows = [[metric.key, 'yes' if metric.gated else 'no'] for metric in spec.metrics]
    lead_width = sum(
        _column_width(header, [row[i] for row in lead_rows]) + len(_GAP) for i, header in enumerate(lead_headers)
    )

    groups = _label_groups(labels, columns, lead_width)
    blocks = [_titled('Metric presence (which general-stats section carries each metric)')]
    first = 1
    for group in groups:
        rows = [[*lead, *(columns[label][i] for label in group)] for i, lead in enumerate(lead_rows)]
        rendered = table([*lead_headers, *group], rows)
        if len(groups) > 1:
            last = first + len(group) - 1
            blocks.append(f'cohorts {first}-{last} of {len(labels)}\n{rendered}')
            first = last + 1
        else:
            blocks.append(rendered)
    return '\n\n'.join(blocks)


def _survey_missing(result: CollectResult) -> str:
    by_label = {row.label: row for row in result.rows}
    lines = [
        f'{_LOUD}MISSING GATED METRICS',
        f'{_LOUD}A gated metric absent from a cohort makes that gate inert there - it checks nothing',
        f'{_LOUD}at all - so the cache is marked incomplete and nothing downstream will run on it.',
    ]
    lines += [f'{_LOUD}  {label}: {", ".join(metrics)}' for label, metrics in result.missing_gated.items()]
    lines.append(f'{_LOUD}Keys present in each section of the affected cohorts - map a MultiQC rename from here:')
    for label in result.missing_gated:
        survey = by_label.get(label)
        if survey is None:  # a cohort can only be missing a metric if it was surveyed
            continue
        for section, keys in survey.section_keys.items():
            body = ', '.join(keys) or '(no keys)'
            lines += _wrap(body, initial=f'{_LOUD}  {label} / {section}: ', subsequent=f'{_LOUD}    ')
    return '\n'.join(lines)


def _survey_failures(result: CollectResult) -> str:
    lines = [
        f'{_LOUD}UNREADABLE COHORTS',
        f'{_LOUD}These reports could not be parsed at all, so they are absent from the tables above.',
        f'{_LOUD}That is a different problem from a cohort missing a metric, and it needs a different',
        f'{_LOUD}fix: check the URI, the file and your credentials, then re-run collect.',
    ]
    for label, error in result.failures:
        lines += _wrap(error, initial=f'{_LOUD}  {label}: ', subsequent=f'{_LOUD}    ')
    return '\n'.join(lines)


# ----------------------------------------------------------------------- distributions


def distributions_report(cache: ValueCache, spec: CalibrationSpec) -> str:
    """Per-metric percentile tables, one row per cohort."""
    blocks = [_titled('Value distributions by cohort')]
    blocks += [_distribution_block(cache, metric) for metric in spec.metrics]
    return _join(blocks)


def _distribution_block(cache: ValueCache, metric: MetricSpec) -> str:
    headers = ['cohort', 'n', *(f'p{pct}' for pct in stats.PERCENTILES)]
    rows = []
    for label in cache.labels:
        values = cache.series(label, metric.key)
        pcts = stats.percentiles(values)
        cells = [fmt_measure(v, metric.unit) for v in pcts] if pcts else ['-'] * len(stats.PERCENTILES)
        rows.append([label, str(values.size), *cells])
    title = f'{metric.key} ({_sense(metric.direction)}, unit {metric.unit})'
    return f'{_titled(title, "-")}\n{table(headers, rows)}'


# --------------------------------------------------------------------------- flagrates


def flagrates_report(cache: ValueCache, spec: CalibrationSpec) -> str:
    """What each candidate fail/warn pair would flag, per cohort."""
    title = _titled('Flag rates under the candidate thresholds')
    if not spec.gated:
        return _join([title, 'This spec has no gated metrics, so there are no thresholds to score.'])
    blocks = [title]
    blocks += [_flagrate_block(cache, metric) for metric in spec.gated]
    blocks.append(_flagrate_legend())
    return _join(blocks)


def _tiers(metric: MetricSpec) -> str:
    if metric.relative is not None:
        warn = f'warn = cohort-relative (MAD k={metric.relative.k:g}, min_cohort={metric.relative.min_cohort})'
    else:
        warn = f'warn = {_fmt_threshold(metric.warn, metric.unit)}'
    return f'fail = {_fmt_threshold(metric.fail, metric.unit)}  {warn}'


def _flagrate_block(cache: ValueCache, metric: MetricSpec) -> str:
    headers = ['cohort', 'n', 'fail%', 'warn%', 'review']
    rows = []
    review_anywhere = False
    for label in cache.labels:
        values = cache.series(label, metric.key)
        fail_rate, warn_rate = stats.flag_rates(values, metric.direction, metric.fail, metric.warn)
        review = stats.needs_review(fail_rate, warn_rate)
        review_anywhere = review_anywhere or review
        rows.append(
            [
                label,
                str(values.size),
                _fmt_rate(fail_rate) if metric.fail is not None else '-',
                _fmt_rate(warn_rate) if metric.warn is not None else '-',
                '*' if review else '',
            ],
        )
    title = f'{metric.key} ({_sense(metric.direction)}, unit {metric.unit})  {_tiers(metric)}'
    if review_anywhere:
        title += '  [REVIEW]'
    return f'{_titled(title, "-")}\n{table(headers, rows)}'


def _flagrate_legend() -> str:
    return '\n'.join(
        [
            _titled('Legend', '-'),
            *_wrap(
                'warn% excludes samples already failing, matching production, which evaluates fail before warn '
                'and records one flag per metric at the worst tier.',
            ),
            *_wrap(
                f'* marks a cohort whose fail rate exceeds {stats.FAIL_RATE_LIMIT:.0%} or whose warn rate exceeds '
                f'{stats.WARN_RATE_LIMIT:.0%}; [REVIEW] marks a metric where that happened in at least one cohort.',
            ),
            *_wrap(
                'A mark is a prompt to look, not a rejection - whether a cohort is healthy is your judgement, '
                'not a computable property, so read the threshold against the cohort and decide.',
            ),
            *_wrap(
                "A rate of '-' means there is no threshold to score it against: either the warn tier is "
                'cohort-relative (the header says so, and `qc_calibrate mad` scores it) or the metric has no '
                'tier at that severity at all.',
            ),
        ],
    )


# --------------------------------------------------------------------------------- mad


def mad_report(evaluation: MadEvaluation) -> str:
    """One candidate cohort-relative warn tier: the numbers, the churn, the verdict."""
    title = _titled(f'Cohort-relative (MAD) evaluation: {evaluation.metric} ({_sense(evaluation.direction)})')
    return _join(
        [
            title,
            _mad_cohorts(evaluation),
            _mad_homogeneous(evaluation),
            _mad_heterogeneous(evaluation),
            _mad_verdict(evaluation),
        ],
    )


def _mad_cohorts(evaluation: MadEvaluation) -> str:
    headers = ['cohort', 'values', 'samples', 'median', 'MAD', 'threshold', 'warn', 'warn rate', 'skipped']
    rows = [
        [
            cohort.label,
            str(cohort.n_values),
            str(cohort.n_samples),
            _fmt_stat(cohort.median),
            _fmt_stat(cohort.mad_raw),
            _fmt_stat(cohort.threshold),
            '-' if cohort.skipped else str(cohort.n_warn),
            '-' if cohort.skipped else f'{cohort.warn_rate:.1%}',
            cohort.skipped or '',
        ]
        for cohort in evaluation.cohorts
    ]
    lines = [_titled('Per-cohort median, MAD and derived threshold', '-'), table(headers, rows)]
    if duplicated := [c.label for c in evaluation.cohorts if c.duplicated]:
        lines.append(
            _fill(
                f'Note: some cohorts carry more values than samples ({", ".join(duplicated)}). The metric appears '
                'in more than one MultiQC section, so each affected sample contributes one value per section. The '
                'warn rate above is per value - the same basis the threshold is derived on - and overstates the '
                'per-sample rate, which the cache cannot answer because it stores values without sample identity.',
            ),
        )
    return '\n'.join(lines)


def _churn_cells(result: ChurnResult | None) -> list[str]:
    return ['-', '-'] if result is None else [str(result.flips), f'{result.flip_rate:.1%}']


def _mad_homogeneous(evaluation: MadEvaluation) -> str:
    title = _titled('Cohort-growth churn: a 60% before-slice scored against the full cohort threshold', '-')
    if not evaluation.homogeneous:
        return f'{title}\nNo cohort was large enough for a growth simulation.'
    headers = [
        'cohort',
        'before n',
        'ordered flips',
        'ordered churn',
        'shuffled flips',
        'shuffled churn',
        'worse',
        'note',
    ]
    rows = []
    for growth in evaluation.homogeneous:
        measured = growth.ordered or growth.shuffled
        rows.append(
            [
                growth.label,
                str(measured.n_initial) if measured else '-',
                *_churn_cells(growth.ordered),
                *_churn_cells(growth.shuffled),
                f'{growth.flip_rate:.1%}',
                'ORDERING-SENSITIVE' if growth.ordering_sensitive else '',
            ],
        )
    lines = [title, table(headers, rows)]
    if evaluation.ordering_sensitive:
        lines.append(
            _fill(
                f'ORDERING-SENSITIVE ({", ".join(evaluation.ordering_sensitive)}): the leading and the shuffled '
                'before-slice disagree on whether this cohort clears the churn bar. Both are shown because the '
                'answer depends on whether MultiQC key order tracks sequencing batches, which is plausible but '
                'nowhere guaranteed; the verdict below uses the worse of the two.',
            ),
        )
    return '\n'.join(lines)


def _mad_heterogeneous(evaluation: MadEvaluation) -> str:
    title = _titled('Merged-cohort churn: one cohort re-scored against the threshold a second one gives it', '-')
    if not evaluation.heterogeneous:
        return f'{title}\nNo cohort pair was large enough for a merge simulation.'
    worst = sorted(evaluation.heterogeneous, key=lambda pair: pair[2].flip_rate, reverse=True)[:MAX_PAIRS_SHOWN]
    headers = ['cohort', 'merged with', 'n', 'flagged before', 'flagged after', 'flips', 'churn']
    rows = [
        [
            label_a,
            label_b,
            str(r.n_initial),
            str(r.flagged_before),
            str(r.flagged_after),
            str(r.flips),
            f'{r.flip_rate:.1%}',
        ]
        for label_a, label_b, r in worst
    ]
    tally = f'Showing the worst {len(worst)} of {len(evaluation.heterogeneous)} pairs simulated.'
    return '\n'.join([title, table(headers, rows), tally])


def _mad_verdict(evaluation: MadEvaluation) -> str:
    reason = f' - {evaluation.verdict_reason}' if evaluation.verdict_reason else ''
    return '\n'.join(
        [
            f'Peak warn rate {evaluation.max_warn_rate:.1%}; peak cohort-growth churn {evaluation.max_churn:.1%}.',
            f'Verdict: {evaluation.verdict}{reason}',
            *_wrap(
                "Adoption is the operator's decision: this verdict is advice for you to sign off, not an "
                'automatic gate. Nothing here changes the spec.',
            ),
        ],
    )


# ----------------------------------------------------------------------------- suggest


def _fmt_seed(value: float | None) -> str:
    return '-' if value is None else f'{value:g}'


def suggest_summary(seeded: Sequence[Seeded]) -> str:
    """What `seed` proposed, on what evidence, and why none of it is a decision yet."""
    title = _titled('Seeded thresholds')
    if not seeded:
        return _join(
            [
                title,
                _fill(
                    'Nothing was seeded: every gated metric is either already reviewed, or has no values in any '
                    'cohort in the cache.',
                ),
            ],
        )
    rows = [[s.key, _fmt_seed(s.fail), _fmt_seed(s.warn)] for s in seeded]
    evidence = ['Evidence']
    for s in seeded:
        evidence += _wrap(s.evidence, initial=f'  {s.key}: ', subsequent='    ')
    caveat = _fill(
        'Every metric above was written with `reviewed = false`. These are starting points, not decisions: a '
        "percentile cannot know whether a threshold contradicts the lab's intent, or tracks ancestry rather than "
        'sample quality. `qc_calibrate emit-config` will refuse to emit a gated metric until you have checked its '
        'numbers and set `reviewed = true` in the spec.',
    )
    return _join([title, table(['metric', 'fail', 'warn'], rows), '\n'.join(evidence), caveat])


# ------------------------------------------------------------------------------ dryrun


def dryrun_summary(result: DryRunResult) -> str:
    """One cohort's dry run: what the production check flagged, how, and at what cost."""
    title = _titled(f'Dry run: {result.cohort}')
    facts = '\n'.join(
        [
            f'samples flagged: {result.n_samples_flagged}',
            f'elapsed: {result.seconds:.1f}s; peak RSS: {result.peak_rss_gb:.2f} GB',
            f'output: {result.output_path}',
        ],
    )
    if not result.counts:
        return _join([title, facts, 'No flags were raised.'])
    rows = [
        [metric, severity, method, str(count)] for (metric, severity, method), count in sorted(result.counts.items())
    ]
    tally = f'{sum(result.counts.values())} flags across {len(result.counts)} (metric, severity, method) groups.'
    return _join([title, facts, table(['metric', 'severity', 'method', 'count'], rows), tally])
