"""Unit tests for the stdout tables every calibration command prints.

Fixtures are built from the real dataclasses rather than stand-in dicts, so a field
rename in `collect`, `relative`, `suggest` or `dryrun` breaks these tests instead of
silently rendering an empty column.
"""

import subprocess
import sys
from collections import Counter

import pytest

from align_genotype.qc_calibration import report, stats
from align_genotype.qc_calibration import spec as spec_mod
from align_genotype.qc_calibration.cache import CohortValues, ValueCache
from align_genotype.qc_calibration.collect import CollectResult, SurveyRow
from align_genotype.qc_calibration.dryrun import DryRunResult
from align_genotype.qc_calibration.relative import CohortMad, HomogeneousChurn, MadEvaluation
from align_genotype.qc_calibration.stats import ChurnResult
from align_genotype.qc_calibration.suggest import Seeded

SPEC = spec_mod.loads(
    'seq_type = "genome"\ncache = "cache.json"\n'
    '\n[metrics.MEDIAN_COVERAGE]\ndirection = "min"\nunit = "x"\nfail = 15\nwarn = 25\n'
    '\n[metrics.FREEMIX]\ndirection = "max"\nunit = "frac"\nfail = 0.04\nwarn = 0.02\n'
    '\n[metrics.DUP]\ndirection = "max"\nunit = "%"\nfail = 40\n'
    '[metrics.DUP.relative]\nk = 3.5\nmin_cohort = 50\n'
    '\n[metrics.error_rate]\ndirection = "max"\nunit = "frac"\ngated = false\n',
)

CACHE = ValueCache(
    seq_type='genome',
    generated='2026-08-11T00:00:00+00:00',
    complete=True,
    metrics=('MEDIAN_COVERAGE', 'FREEMIX', 'DUP', 'error_rate'),
    cohorts=(
        CohortValues(
            label='dataset-a',
            n_samples=3,
            multiqc_version='1.33',
            shape='dict',
            n_dropped=0,
            values={
                'MEDIAN_COVERAGE': [30.0, 32.0, 34.0],
                'FREEMIX': [0.001, 0.002, 0.003],
                'DUP': [10.0, 12.0, 14.0],
                'error_rate': [],
            },
        ),
        CohortValues(
            label='dataset-b',
            n_samples=3,
            multiqc_version='1.21',
            shape='list',
            n_dropped=2,
            values={
                'MEDIAN_COVERAGE': [10.0, 12.0, 14.0],
                'FREEMIX': [0.001, 0.002, 0.003],
                'DUP': [11.0, 13.0, 15.0],
                'error_rate': [],
            },
        ),
    ),
)


def survey_rows():
    """Two surveyed cohorts; dataset-b carries no DUP key at all."""
    return (
        SurveyRow(
            label='dataset-a',
            uri='gs://bucket-a/multiqc.json',
            multiqc_version='1.33',
            shape='dict',
            n_samples=3,
            section_sizes={'general': 3, 'picard_1': 3},
            where={
                'MEDIAN_COVERAGE': ('general',),
                'FREEMIX': ('general',),
                'DUP': ('picard_1',),
                'error_rate': ('general',),
            },
            section_keys={
                'general': ('FREEMIX', 'MEDIAN_COVERAGE', 'error_rate'),
                'picard_1': ('DUP',),
            },
            n_dropped=0,
            n_dropped_by_metric={'MEDIAN_COVERAGE': 0, 'FREEMIX': 0, 'DUP': 0, 'error_rate': 0},
        ),
        SurveyRow(
            label='dataset-b',
            uri='gs://bucket-b/multiqc.json',
            multiqc_version='1.21',
            shape='list',
            n_samples=3,
            section_sizes={'general': 3},
            where={
                'MEDIAN_COVERAGE': ('general',),
                'FREEMIX': ('general',),
                'DUP': (),
                'error_rate': (),
            },
            section_keys={'general': ('FREEMIX', 'MEDIAN_COVERAGE', 'PERCENT_DUPLICATION')},
            n_dropped=2,
            n_dropped_by_metric={'MEDIAN_COVERAGE': 2, 'FREEMIX': 0, 'DUP': 0, 'error_rate': 0},
        ),
    )


def collect_result(missing_gated=None, failures=()):
    return CollectResult(
        cache=CACHE,
        rows=survey_rows(),
        missing_gated={} if missing_gated is None else missing_gated,
        failures=failures,
    )


# ---------------------------------------------------------------------------- table


def test_table_pads_every_column_to_its_widest_cell():
    out = report.table(['a', 'longheader'], [['xxx', 'y'], ['z', 'w']])
    assert out.splitlines() == ['a    longheader', '---  ----------', 'xxx  y', 'z    w']


def test_table_rule_spans_every_column():
    assert report.table(['metric', 'n'], [['DUP', '3']]).splitlines()[1] == '------  -'


def test_table_with_no_rows_is_header_and_rule():
    assert report.table(['a', 'b'], []).splitlines() == ['a  b', '-  -']


def test_table_rejects_a_row_of_the_wrong_width():
    with pytest.raises(ValueError, match='2 headers'):
        report.table(['a', 'b'], [['only-one']])


def test_table_cells_start_at_the_same_offset_in_every_row():
    lines = report.table(['cohort', 'n'], [['dataset-a', '3'], ['b', '11']]).splitlines()
    assert lines[2].index('3') == lines[3].index('11')
    assert lines[0].index('n') == lines[2].index('3')


# --------------------------------------------------------------------- survey_report


def test_survey_report_shows_version_shape_samples_and_drops():
    out = report.survey_report(collect_result(), SPEC)
    assert '1.33' in out
    assert '1.21' in out
    assert 'dict' in out
    assert 'list' in out
    assert 'general=3' in out
    assert 'picard_1=3' in out
    # the per-metric drop breakdown, not just the total
    assert 'MEDIAN_COVERAGE=2' in out


def test_survey_report_presence_matrix_names_the_carrying_section():
    out = report.survey_report(collect_result(), SPEC)
    dup_row = next(line for line in out.splitlines() if line.startswith('DUP'))
    assert 'picard_1' in dup_row
    assert 'MISSING' in dup_row  # dataset-b


def wide_survey_row(label):
    """A surveyed cohort whose label is as long as a real dataset label."""
    return SurveyRow(
        label=label,
        uri=f'gs://bucket/{label}.json',
        multiqc_version='1.33',
        shape='dict',
        n_samples=120,
        section_sizes={'general': 120, 'picard_1': 118},
        where={metric.key: ('picard_1',) for metric in SPEC.metrics},
        section_keys={'general': tuple(m.key for m in SPEC.metrics)},
        n_dropped=0,
        n_dropped_by_metric={},
    )


def matrix_headers(out):
    return [line for line in out.splitlines() if line.startswith('metric ') and 'gated' in line]


def test_survey_report_chunks_a_presence_matrix_too_wide_for_a_terminal():
    """Ten real dataset labels put one matrix past 240 columns; chunk, never truncate."""
    labels = [f'cohort-with-a-long-name-{i}' for i in range(10)]
    result = CollectResult(
        cache=CACHE,
        rows=tuple(wide_survey_row(label) for label in labels),
        missing_gated={},
        failures=(),
    )
    out = report.survey_report(result, SPEC)
    assert max(len(line) for line in out.splitlines()) <= report.MAX_TABLE_WIDTH
    for label in labels:
        assert label in out  # every cohort is still there
    assert len(matrix_headers(out)) > 1  # the matrix was split into groups
    assert f'of {len(labels)}' in out  # and each group says which cohorts it covers


def test_survey_report_does_not_chunk_a_narrow_presence_matrix():
    out = report.survey_report(collect_result(), SPEC)
    assert len(matrix_headers(out)) == 1


def test_survey_report_marks_gated_metrics():
    lines = report.survey_report(collect_result(), SPEC).splitlines()
    header = next(line for line in lines if line.startswith('metric '))
    assert 'gated' in header
    assert next(line for line in lines if line.startswith('error_rate')).split()[1] == 'no'
    assert next(line for line in lines if line.startswith('MEDIAN_COVERAGE')).split()[1] == 'yes'


def test_survey_report_dumps_candidate_keys_for_a_missing_gated_metric():
    """The whole point: an operator has to be able to spot the MultiQC rename."""
    out = report.survey_report(collect_result(missing_gated={'dataset-b': ('DUP',)}), SPEC)
    assert 'MISSING GATED METRICS' in out
    assert 'PERCENT_DUPLICATION' in out
    assert 'dataset-b' in out


def test_survey_report_no_missing_block_when_nothing_is_missing():
    out = report.survey_report(collect_result(), SPEC)
    assert 'MISSING GATED METRICS' not in out
    assert 'PERCENT_DUPLICATION' not in out


def test_survey_report_renders_unreadable_cohorts_distinctly_from_missing_metrics():
    out = report.survey_report(
        collect_result(missing_gated={'dataset-b': ('DUP',)}, failures=(('dataset-c', 'HTTP 403 on bucket'),)),
        SPEC,
    )
    assert 'UNREADABLE COHORTS' in out
    assert 'MISSING GATED METRICS' in out
    assert out.index('MISSING GATED METRICS') < out.index('UNREADABLE COHORTS')
    assert 'HTTP 403 on bucket' in out
    # the failed cohort was never surveyed, so it must not appear as a matrix column
    matrix_rows = [line for line in out.splitlines() if line.startswith('DUP')]
    assert matrix_rows
    assert all('dataset-c' not in line for line in matrix_rows)


def test_survey_report_says_so_when_everything_is_fine():
    out = report.survey_report(collect_result(), SPEC)
    assert 'UNREADABLE COHORTS' not in out
    assert 'OK' in out


def test_survey_report_ok_line_absent_when_there_are_failures():
    out = report.survey_report(collect_result(failures=(('dataset-c', 'boom'),)), SPEC)
    assert 'every gated metric' not in out


# -------------------------------------------------------------- distributions_report


def test_distributions_report_has_a_column_per_percentile():
    out = report.distributions_report(CACHE, SPEC)
    for pct in stats.PERCENTILES:
        assert f'p{pct}' in out


def test_distributions_report_shows_direction_sense_and_unit():
    out = report.distributions_report(CACHE, SPEC)
    assert 'higher=better' in out  # MEDIAN_COVERAGE is a 'min' metric
    assert 'lower=better' in out  # FREEMIX is a 'max' metric
    assert 'unit x' in out
    assert 'unit frac' in out


def test_distributions_report_formats_values_through_fmt_measure():
    out = report.distributions_report(CACHE, SPEC)
    assert '32.0' in out  # MEDIAN_COVERAGE p50, unit 'x' -> 1 dp
    assert '0.002' in out  # FREEMIX p50, unit 'frac' -> 3 dp


def test_distributions_report_covers_ungated_metrics_too():
    assert 'error_rate' in report.distributions_report(CACHE, SPEC)


def test_distributions_report_marks_a_cohort_with_no_values():
    block = report.distributions_report(CACHE, SPEC).split('error_rate')[1]
    row = next(line for line in block.splitlines() if line.startswith('dataset-a'))
    assert row.split()[1] == '0'
    assert '-' in row


# ------------------------------------------------------------------ flagrates_report


def test_flagrates_report_shows_candidate_thresholds():
    out = report.flagrates_report(CACHE, SPEC)
    assert 'fail = 15.0' in out
    assert 'warn = 25.0' in out
    assert 'fail = 0.040' in out


def test_flagrates_report_shows_per_cohort_fail_and_warn_rates():
    out = report.flagrates_report(CACHE, SPEC)
    # dataset-b MEDIAN_COVERAGE sits entirely below fail = 15
    block = out.split('MEDIAN_COVERAGE')[1].split('FREEMIX')[0]
    assert '100.0%' in next(line for line in block.splitlines() if line.startswith('dataset-b'))
    assert '0.0%' in next(line for line in block.splitlines() if line.startswith('dataset-a'))


def test_flagrates_report_marks_a_metric_needing_review():
    out = report.flagrates_report(CACHE, SPEC)
    assert 'REVIEW' in next(line for line in out.splitlines() if line.startswith('MEDIAN_COVERAGE'))


def test_flagrates_report_does_not_mark_a_healthy_metric():
    out = report.flagrates_report(CACHE, SPEC)
    assert 'REVIEW' not in next(line for line in out.splitlines() if line.startswith('FREEMIX'))


def test_flagrates_report_names_the_relative_warn_tier():
    out = report.flagrates_report(CACHE, SPEC)
    assert 'cohort-relative' in next(line for line in out.splitlines() if line.startswith('DUP '))


def test_flagrates_report_omits_ungated_metrics():
    assert 'error_rate' not in report.flagrates_report(CACHE, SPEC)


def test_flagrates_report_legend_explains_warn_exclusion_and_the_mark():
    out = report.flagrates_report(CACHE, SPEC)
    assert 'excludes samples already failing' in out
    assert f'{stats.FAIL_RATE_LIMIT:.0%}' in out
    assert f'{stats.WARN_RATE_LIMIT:.0%}' in out
    assert 'not a rejection' in out


def test_flagrates_report_handles_a_spec_with_no_gated_metrics():
    ungated = spec_mod.loads(
        'seq_type = "genome"\ncache = "cache.json"\n\n[metrics.error_rate]\ndirection = "max"\ngated = false\n',
    )
    assert 'no gated metrics' in report.flagrates_report(CACHE, ungated).lower()


# ------------------------------------------------------------------------ mad_report


def churn(n_initial=60, flips=0, before=1.0, after=1.1):
    return ChurnResult(
        threshold_before=before,
        threshold_after=after,
        n_initial=n_initial,
        flagged_before=flips,
        flagged_after=0,
        flips=flips,
    )


COHORT_MAD = CohortMad(
    label='dataset-a',
    n_values=100,
    n_samples=100,
    median=12.0,
    mad_raw=1.5,
    threshold=17.7825,
    n_warn=3,
    skipped=None,
)


def mad_evaluation(homogeneous=(), heterogeneous=(), cohorts=(COHORT_MAD,)):
    return MadEvaluation(
        metric='DUP',
        direction='max',
        cohorts=cohorts,
        homogeneous=homogeneous,
        heterogeneous=heterogeneous,
    )


def test_mad_report_shows_the_metric_and_direction_sense():
    out = report.mad_report(mad_evaluation())
    assert 'DUP' in out
    assert 'lower=better' in out


def test_mad_report_shows_per_cohort_numbers():
    row = next(line for line in report.mad_report(mad_evaluation()).splitlines() if line.startswith('dataset-a'))
    assert '12' in row  # median
    assert '1.5' in row  # MAD
    assert '17.78' in row  # derived threshold
    assert '3.0%' in row  # warn rate: 3 of 100 values


def test_mad_report_notes_a_skipped_cohort_and_its_reason():
    skipped = CohortMad(
        label='cohort-2',
        n_values=10,
        n_samples=10,
        median=12.0,
        mad_raw=1.0,
        threshold=None,
        n_warn=0,
        skipped='cohort 10 < min_cohort 50',
    )
    assert 'cohort 10 < min_cohort 50' in report.mad_report(mad_evaluation(cohorts=(skipped,)))


def test_mad_report_flags_a_cohort_carrying_more_values_than_samples():
    duplicated = CohortMad(
        label='dataset-a',
        n_values=120,
        n_samples=100,
        median=12.0,
        mad_raw=1.5,
        threshold=17.7825,
        n_warn=4,
        skipped=None,
    )
    out = report.mad_report(mad_evaluation(cohorts=(duplicated,)))
    assert '120' in out
    assert '100' in out
    assert 'more values than samples' in out


def test_mad_report_shows_both_ordered_and_shuffled_churn():
    homogeneous = (HomogeneousChurn(label='dataset-a', ordered=churn(flips=3), shuffled=churn(flips=0)),)
    out = report.mad_report(mad_evaluation(homogeneous=homogeneous))
    assert 'ordered' in out.lower()
    assert 'shuffled' in out.lower()
    rows = [line for line in out.splitlines() if line.startswith('dataset-a')]
    assert any('5.0%' in line and '0.0%' in line for line in rows)


def test_mad_report_flags_ordering_sensitive_cohorts():
    homogeneous = (HomogeneousChurn(label='dataset-a', ordered=churn(flips=3), shuffled=churn(flips=0)),)
    evaluation = mad_evaluation(homogeneous=homogeneous)
    assert evaluation.ordering_sensitive == ('dataset-a',)
    out = report.mad_report(evaluation)
    assert 'ORDERING-SENSITIVE' in out
    assert 'dataset-a' in out


def test_mad_report_handles_a_degenerate_slice():
    """One ordering can have a zero MAD while the other measures fine."""
    homogeneous = (HomogeneousChurn(label='dataset-a', ordered=None, shuffled=churn(flips=0)),)
    rows = [
        line
        for line in report.mad_report(mad_evaluation(homogeneous=homogeneous)).splitlines()
        if line.startswith('dataset-a')
    ]
    growth_row = rows[1].split()  # rows[0] is the per-cohort table
    assert growth_row[2:4] == ['-', '-']  # ordered flips and churn
    assert growth_row[4:6] == ['0', '0.0%']  # shuffled measured fine


def test_mad_report_says_when_no_cohort_could_be_simulated():
    assert 'no cohort' in report.mad_report(mad_evaluation()).lower()


def test_mad_report_shows_the_worst_heterogeneous_pairs_first():
    heterogeneous = tuple((f'cohort-{i}', f'cohort-{i + 1}', churn(n_initial=100, flips=i)) for i in range(1, 9))
    out = report.mad_report(mad_evaluation(heterogeneous=heterogeneous))
    pair_lines = [line for line in out.splitlines() if line.startswith('cohort-')]
    assert pair_lines[0].startswith('cohort-8')  # 8 flips of 100 is the worst pair
    assert len(pair_lines) == report.MAX_PAIRS_SHOWN
    assert '8 pairs' in out  # says how many were simulated


def test_mad_report_verdict_reject_carries_its_reason():
    homogeneous = (HomogeneousChurn(label='dataset-a', ordered=churn(flips=3), shuffled=churn(flips=3)),)
    evaluation = mad_evaluation(homogeneous=homogeneous)
    assert evaluation.verdict == 'REJECT'
    out = report.mad_report(evaluation)
    assert 'REJECT' in out
    assert evaluation.verdict_reason in out


def test_mad_report_verdict_recommend_needs_no_reason():
    evaluation = mad_evaluation(homogeneous=(HomogeneousChurn('dataset-a', churn(), churn()),))
    assert evaluation.verdict == 'RECOMMEND'
    assert 'RECOMMEND' in report.mad_report(evaluation)


def test_mad_report_says_adoption_is_the_operators_decision():
    out = report.mad_report(mad_evaluation())
    assert 'advice' in out.lower()
    assert 'operator' in out.lower()


# ------------------------------------------------------------------- suggest_summary


def test_suggest_summary_lists_what_was_seeded_and_its_evidence():
    seeded = (
        Seeded(key='FREEMIX', fail=0.04, warn=0.02, evidence='Seeded from p99 across 3 cohorts.'),
        Seeded(key='DUP', fail=40, warn=None, evidence='Seeded from p99 (39.8) across 3 cohorts.'),
    )
    out = report.suggest_summary(seeded)
    assert 'FREEMIX' in out
    assert '0.04' in out
    assert '0.02' in out
    assert '40' in out
    assert 'Seeded from p99 (39.8) across 3 cohorts.' in out


def test_suggest_summary_states_nothing_is_reviewed():
    out = report.suggest_summary((Seeded(key='FREEMIX', fail=0.04, warn=0.02, evidence='ev'),))
    assert 'reviewed = false' in out
    assert 'starting point' in out.lower()
    assert 'emit-config' in out
    assert 'refuse' in out.lower()


def test_suggest_summary_renders_an_absent_warn_tier():
    out = report.suggest_summary((Seeded(key='DUP', fail=40, warn=None, evidence='ev'),))
    assert next(line for line in out.splitlines() if line.startswith('DUP')).split()[2] == '-'


def test_suggest_summary_handles_the_empty_case():
    out = report.suggest_summary(())
    assert 'reviewed = false' not in out
    assert 'nothing' in out.lower()


# -------------------------------------------------------------------- dryrun_summary


def dryrun_result(counts=None):
    return DryRunResult(
        cohort='dataset-a',
        n_samples_flagged=12,
        counts=Counter() if counts is None else Counter(counts),
        seconds=41.234,
        peak_rss_gb=3.4567,
        output_path='out/dryrun_dataset-a.json',
    )


def test_dryrun_summary_groups_counts_by_metric_severity_and_method():
    result = dryrun_result(
        {
            ('DUP', 'warn', 'relative'): 4,
            ('MEDIAN_COVERAGE', 'fail', 'absolute'): 2,
            ('MEDIAN_COVERAGE', 'warn', 'absolute'): 7,
        },
    )
    out = report.dryrun_summary(result)
    assert 'dataset-a' in out
    assert '12' in out
    dup_row = next(line for line in out.splitlines() if line.startswith('DUP'))
    assert dup_row.split() == ['DUP', 'warn', 'relative', '4']
    coverage_rows = [line for line in out.splitlines() if line.startswith('MEDIAN_COVERAGE')]
    assert len(coverage_rows) == 2  # fail and warn are separate rows
    assert coverage_rows[0].split() == ['MEDIAN_COVERAGE', 'fail', 'absolute', '2']
    assert '13 flags' in out  # 4 + 2 + 7


def test_dryrun_summary_shows_timing_rss_and_output_path():
    out = report.dryrun_summary(dryrun_result())
    assert '41.2' in out
    assert '3.46' in out
    assert 'out/dryrun_dataset-a.json' in out


def test_dryrun_summary_says_so_when_nothing_was_flagged():
    assert 'no flags' in report.dryrun_summary(dryrun_result()).lower()


# ------------------------------------------------------------------------ provenance


def test_fmt_measure_is_unchanged():
    assert report.fmt_measure(0.1234, 'frac') == '0.123'
    assert report.fmt_measure(31.25, 'x') == '31.2'
    assert report.fmt_measure(31.25, '%') == '31.2'
    assert report.fmt_measure(31.257, 'other') == '31.26'


def probe(code):
    """Run `code` in a fresh interpreter and hand back the completed process."""
    return subprocess.run([sys.executable, '-c', code], capture_output=True, text=True, check=False)  # noqa: S603


def test_report_does_not_import_dryrun_or_emit_at_runtime():
    """report -> dryrun -> emit -> report would close an import cycle."""
    done = probe(
        'import sys\n'
        'import align_genotype.qc_calibration.report\n'
        "assert 'align_genotype.qc_calibration.dryrun' not in sys.modules, 'dryrun imported'\n"
        "assert 'align_genotype.qc_calibration.emit' not in sys.modules, 'emit imported'\n"
        "assert 'align_genotype.qc_calibration.relative' not in sys.modules, 'relative imported'\n",
    )
    assert done.returncode == 0, done.stderr


@pytest.mark.parametrize(('first', 'second'), [('dryrun', 'report'), ('report', 'dryrun')])
def test_report_and_dryrun_import_in_either_order(first, second):
    done = probe(f'import align_genotype.qc_calibration.{first}\nimport align_genotype.qc_calibration.{second}\n')
    assert done.returncode == 0, done.stderr
