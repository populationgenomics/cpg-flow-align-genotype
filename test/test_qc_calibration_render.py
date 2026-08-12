"""Smoke tests for the calibration HTML report.

Fixtures are built locally rather than imported from `test.test_qc_calibration_summary`:
`test/` has no `__init__.py`, and the repo's Python distribution ships its own stdlib
`test` package in `site-packages` which wins the import over the local directory (a plain
namespace package only forms when no regular package of that name is found on any
earlier `sys.path` entry, and the stdlib one has an `__init__.py`). `from test.foo import
bar` therefore raises `ModuleNotFoundError: No module named 'test.test_qc_calibration_summary'`
even though the file exists. Adding `test/__init__.py` would fix it, but this task's scope
is limited to this file plus `render.py` and the template, so the fixture is rebuilt here
instead of introducing a package-wide change from an out-of-scope file.
"""

import pytest

from align_genotype.qc_calibration import render as render_mod
from align_genotype.qc_calibration import settings as settings_mod
from align_genotype.qc_calibration import summary as summary_mod
from align_genotype.qc_calibration import values as values_mod

SETTINGS = settings_mod.CalibrationSettings(
    seq_type='genome',
    metrics=(
        settings_mod.MetricSpec(key='MEDIAN_COVERAGE', direction='min', unit='x'),
        settings_mod.MetricSpec(key='dup_pct', direction='max', unit='%', relative=True),
        settings_mod.MetricSpec(key='ABSENT', direction='min', unit='x'),
    ),
    k=3.5,
    min_samples=4,
)

# No metric has `relative = True` - used to exercise the "nothing configured" fallback.
SETTINGS_NO_RELATIVE = settings_mod.CalibrationSettings(
    seq_type='genome',
    metrics=(settings_mod.MetricSpec(key='MEDIAN_COVERAGE', direction='min', unit='x'),),
    k=3.5,
    min_samples=4,
)


def dataset_values(name: str, coverage: list[float], dup: list[float], n_dropped: int = 0) -> values_mod.DatasetValues:
    """One dataset's values, mirroring the shape `qc_calibration/extract.py` produces."""
    return values_mod.DatasetValues(
        dataset=name,
        seq_type='genome',
        analysis_id=1,
        timestamp='2026-06-01T00:00:00',
        uri=f'gs://{name}/multiqc_data.json',
        multiqc_version='1.33',
        generated='2026-08-12T00:00:00',
        n_sequencing_groups=len(coverage),
        section_sizes={'picard_1': len(coverage)},
        metrics={
            'MEDIAN_COVERAGE': values_mod.MetricValues(
                entries=tuple(('picard_1', f'{name}-CPG{i}', float(v)) for i, v in enumerate(coverage)),
                n_dropped=n_dropped,
            ),
            'dup_pct': values_mod.MetricValues(
                entries=tuple(('samtools', f'{name}-CPG{i}', float(v)) for i, v in enumerate(dup)),
                n_dropped=0,
            ),
            'ABSENT': values_mod.MetricValues(entries=(), n_dropped=0),
        },
    )


@pytest.fixture
def html() -> str:
    built = summary_mod.build(
        [
            dataset_values('ds-a', [30, 32, 34, 36, 38, 10], [10, 10.5, 11, 11.5, 12, 40]),
            dataset_values('ds-b', [40, 42, 44, 46, 48, 50], [7, 7.5, 8, 8.5, 9, 9.5]),
        ],
        SETTINGS,
        current={'MEDIAN_COVERAGE': {'fail': 15, 'warn': 25}},
        skipped_datasets=[{'dataset': 'ds-c', 'reason': 'no completed CramMultiQC qc analysis for genome'}],
        generated='2026-08-12T00:00:00',
        ar_guid='test-ar-guid',
    )
    return render_mod.render(built)


def test_renders_a_complete_html_document(html: str) -> None:
    assert html.startswith('<!DOCTYPE html>')
    assert html.rstrip().endswith('</html>')


def test_carries_every_expected_section(html: str) -> None:
    for heading in (
        'Recommended fixed thresholds',
        'Dataset-relative (MAD) tiers',
        'Metric presence',
        'Percentile distributions',
        'Flag rates',
        'Churn detail',
        'Provenance',
        'Config block',
    ):
        assert heading in html, f'missing section: {heading}'


def test_headline_reports_the_run_shape(html: str) -> None:
    assert 'genome' in html
    assert 'test-ar-guid' in html


def test_shows_current_and_candidate_side_by_side(html: str) -> None:
    assert 'MEDIAN_COVERAGE' in html
    assert 'Current' in html
    assert 'Candidate' in html


def test_names_the_skipped_dataset_and_why(html: str) -> None:
    assert 'ds-c' in html
    assert 'no completed CramMultiQC qc analysis' in html


def test_banners_a_metric_absent_from_every_dataset(html: str) -> None:
    assert 'ABSENT' in html
    assert 'checks nothing' in html


def test_includes_the_pasteable_config_block(html: str) -> None:
    assert '[qc_thresholds.genome.fail.min]' in html


def test_states_that_candidates_are_not_decisions(html: str) -> None:
    assert 'not decisions' in html


def test_explains_the_merge_churn_tension(html: str) -> None:
    assert 'merge' in html.lower()
    assert 'shifts between datasets' in html


def test_frames_reject_as_a_live_question_not_a_tool_verdict(html: str) -> None:
    """Both shipped relative tiers are expected to REJECT on merge churn; the page must
    not present that as a tool-issued failure."""
    assert 'REJECT' in html
    assert 'live question for the team' in html
    assert 'not a tool verdict' in html


def test_escapes_html_in_dataset_names() -> None:
    """Dataset names reach the page unmodified, so autoescape must be on."""
    built = summary_mod.build(
        [dataset_values('<script>x</script>', [30, 32, 34, 36], [1, 2, 3, 4])],
        SETTINGS,
        current={},
        skipped_datasets=[],
        generated='2026-08-12T00:00:00',
        ar_guid='x',
    )
    assert '<script>x</script>' not in render_mod.render(built)


def test_escapes_html_in_skip_reasons() -> None:
    """Skip reasons are free text from discovery/extract and must also be escaped."""
    built = summary_mod.build(
        [dataset_values('ds-a', [30, 32, 34, 36], [1, 2, 3, 4])],
        SETTINGS,
        current={},
        skipped_datasets=[{'dataset': 'ds-c', 'reason': '<img src=x onerror=alert(1)>'}],
        generated='2026-08-12T00:00:00',
        ar_guid='x',
    )
    rendered = render_mod.render(built)
    assert '<img src=x onerror=alert(1)>' not in rendered
    assert '&lt;img' in rendered


# --- Tests beyond the plan's list -------------------------------------------------------
#
# The plan's single fixture (two 6-value datasets, min_samples=4) cannot exercise several
# branches the template renders differently:
#
# * Growth needs a 60%-sized before-slice >= min_samples (0.6 * 6 = 3.6 < 4), so growth is
#   *always* empty for the plan's fixture - the exact "simulation that could never run"
#   failure mode flagged as having recurred before. A dataset with >= 7 values is needed
#   to render a populated growth row at all.
# * Merge truncation to MAX_PAIRS=10 needs more than ~4 datasets (n*(n-1) pairs); the
#   plan's 2-dataset fixture only ever produces 2 pairs, so the "Worst 10 of N" wording
#   and the truncated table were never rendered.
# * `d.skipped` (a dataset excluded from its own per-dataset relative row) needs a
#   dataset with too few values, or zero values, for the relative metric - the plan's
#   fixture never puts a dataset below `min_samples` for `dup_pct`.
# * The relative section's for/else "no relative metric configured" fallback is only
#   reachable when settings.relative_metrics is empty, never true for the plan's SETTINGS.
# * A relative metric configured but with zero datasets (`build([], ...)`) takes a
#   different path through the same fallback in `summary.py` - the wording must not claim
#   "no metric is configured" when one is.
# * `r.skipped_datasets == []` (no "Datasets not included" heading) is the mirror of the
#   already-tested populated case.
# * A zero-dataset run end to end, to confirm the page still renders rather than raising
#   under `StrictUndefined` when every candidate is `None`.


def test_growth_churn_renders_when_a_dataset_is_large_enough() -> None:
    """A 10-value dataset clears the 0.6 * min_samples=4 threshold for a growth row."""
    coverage = list(range(20, 30))
    dup = [10, 10.5, 11, 11.5, 12, 12.5, 13, 13.5, 14, 14.5]
    built = summary_mod.build(
        [dataset_values('ds-big', coverage, dup)],
        SETTINGS,
        current={},
        skipped_datasets=[],
        generated='2026-08-12T00:00:00',
        ar_guid='x',
    )
    assert built['relative']['dup_pct']['n_growth_evaluated'] == 1  # sanity on the fixture itself
    rendered = render_mod.render(built)
    assert 'ds-big' in rendered
    assert '1 of 1 dataset(s) large enough' in rendered


def test_merge_pairs_truncate_to_ten_worst_in_the_rendered_table() -> None:
    """More than MAX_PAIRS=10 ordered pairs must render as a truncated 'worst N of M'."""
    dup_base = [10, 10.5, 11, 11.5, 12, 12.5, 13, 13.5, 14, 14.5]
    datasets = [dataset_values(f'ds-{i}', list(range(10)), [v + i * 0.3 for v in dup_base]) for i in range(5)]
    built = summary_mod.build(
        datasets,
        SETTINGS,
        current={},
        skipped_datasets=[],
        generated='2026-08-12T00:00:00',
        ar_guid='x',
    )
    evaluation = built['relative']['dup_pct']
    assert evaluation['merge_pairs_simulated'] == 20  # sanity: 5 * 4 ordered pairs
    rendered = render_mod.render(built)
    assert 'Worst 10 of 20 ordered pairs' in rendered
    assert '5 usable dataset(s)' in rendered


def test_a_dataset_too_small_for_the_relative_metric_shows_its_skip_reason() -> None:
    """A dataset below `min_samples` for the relative metric must show why, not a blank row."""
    built = summary_mod.build(
        [
            dataset_values('ds-ok', [30, 32, 34, 36], [10, 10.5, 11, 11.5, 12, 12.5, 13, 13.5]),
            dataset_values('ds-tiny', [30, 32], [1, 2]),  # only 2 dup_pct values, min_samples=4
        ],
        SETTINGS,
        current={},
        skipped_datasets=[],
        generated='2026-08-12T00:00:00',
        ar_guid='x',
    )
    rendered = render_mod.render(built)
    assert 'ds-tiny' in rendered
    assert 'min_samples 4' in rendered


def test_coverage_counts_are_shown_next_to_each_peak_when_datasets_are_skipped() -> None:
    """A skipped dataset must reduce the reported denominator, not just the peak figure."""
    big_x = [10, 10.5, 11, 11.5, 12, 12.5, 13, 13.5, 14, 14.5]
    big_y = [10.2, 10.7, 11.2, 11.7, 12.2, 12.7, 13.2, 13.7, 14.2, 14.7]
    tiny_z = [50.0, 55.0]  # below min_samples=4
    built = summary_mod.build(
        [
            dataset_values('ds-x', list(range(10)), big_x),
            dataset_values('ds-y', list(range(10)), big_y),
            dataset_values('ds-z', list(range(2)), tiny_z),
        ],
        SETTINGS,
        current={},
        skipped_datasets=[],
        generated='2026-08-12T00:00:00',
        ar_guid='x',
    )
    rendered = render_mod.render(built)
    assert '2 of 3 dataset(s)' in rendered  # warn-rate denominator: 3 total, 1 skipped
    assert '1 skipped' in rendered
    assert '2 of 3 dataset(s) large enough' in rendered  # growth denominator
    assert '2 usable dataset(s)' in rendered  # merge denominator


def test_relative_fallback_names_the_missing_metric_when_none_is_configured() -> None:
    """When no metric is `relative = true`, the page must say so, not just show an empty table."""
    built = summary_mod.build(
        [dataset_values('ds-a', [30, 32, 34, 36], [1, 2, 3, 4])],
        SETTINGS_NO_RELATIVE,
        current={},
        skipped_datasets=[],
        generated='2026-08-12T00:00:00',
        ar_guid='x',
    )
    assert built['relative'] == {}  # sanity on the fixture itself
    rendered = render_mod.render(built)
    assert 'No metric is configured with' in rendered


def test_relative_fallback_distinguishes_no_data_from_no_metric_configured() -> None:
    """A relative metric IS configured; the fallback text must not claim otherwise."""
    built = summary_mod.build(
        [],
        SETTINGS,
        current={},
        skipped_datasets=[],
        generated='2026-08-12T00:00:00',
        ar_guid='x',
    )
    assert built['relative'] == {}  # sanity: summary.py takes this path when there's no data at all
    rendered = render_mod.render(built)
    assert 'no dataset had data to evaluate it against' in rendered
    assert 'No metric is configured with' not in rendered


def test_omits_the_datasets_not_included_heading_when_nothing_was_skipped(html: str) -> None:  # noqa: ARG001
    built = summary_mod.build(
        [dataset_values('ds-a', [30, 32, 34, 36], [1, 2, 3, 4])],
        SETTINGS,
        current={},
        skipped_datasets=[],
        generated='2026-08-12T00:00:00',
        ar_guid='x',
    )
    rendered = render_mod.render(built)
    assert 'Datasets not included' not in rendered


def test_a_run_with_no_datasets_at_all_still_renders() -> None:
    """Every candidate is None and every table is empty; StrictUndefined must not fire."""
    built = summary_mod.build(
        [],
        SETTINGS,
        current={},
        skipped_datasets=[],
        generated='2026-08-12T00:00:00',
        ar_guid='x',
    )
    rendered = render_mod.render(built)
    assert rendered.startswith('<!DOCTYPE html>')
    assert rendered.rstrip().endswith('</html>')
