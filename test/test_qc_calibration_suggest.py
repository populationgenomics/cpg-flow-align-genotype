"""Unit tests for threshold seeding."""

import pytest

from align_genotype.qc_calibration import spec as spec_mod
from align_genotype.qc_calibration import suggest as suggest_mod
from align_genotype.qc_calibration.cache import CohortValues, ValueCache

SPEC = spec_mod.loads(
    'seq_type = "genome"\ncache = "c.json"\n'
    # min metric, no thresholds yet beyond a placeholder fail
    '\n[metrics.MEDIAN_COVERAGE]\ndirection = "min"\nunit = "x"\nfail = 1\n'
    # max metric
    '\n[metrics.FREEMIX]\ndirection = "max"\nunit = "frac"\nfail = 1\n'
    # relative metric - warn tier is cohort-relative, so only fail should be seeded
    '\n[metrics.DUP]\ndirection = "max"\nunit = "%"\nfail = 1\n[metrics.DUP.relative]\nk = 3.5\nmin_cohort = 5\n'
    # already signed off - must not be touched
    '\n[metrics.PCT_20X]\ndirection = "min"\nunit = "frac"\nfail = 0.75\nwarn = 0.85\nreviewed = true\n'
    # not gated - must not be touched
    '\n[metrics.error_rate]\ndirection = "max"\nunit = "frac"\ngated = false\n',
)

# Cohort A: 10..109, Cohort B: 20..119. Lowest p1 = 10.99, lowest p5 = 14.95.
# FREEMIX 0.000..0.099 in both: highest p99 = 0.09801, highest p95 = 0.09405.
_COHORT_A = {
    'MEDIAN_COVERAGE': [float(v) for v in range(10, 110)],
    'FREEMIX': [v / 1000 for v in range(100)],
    'DUP': [float(v) for v in range(100)],
    'PCT_20X': [0.9] * 100,
    'error_rate': [0.01] * 100,
}
_COHORT_B = {**_COHORT_A, 'MEDIAN_COVERAGE': [float(v) for v in range(20, 120)]}

CACHE = ValueCache(
    seq_type='genome',
    generated='x',
    complete=True,
    metrics=('MEDIAN_COVERAGE', 'FREEMIX', 'DUP', 'PCT_20X', 'error_rate'),
    cohorts=(
        CohortValues('dataset-a', 100, '1.33', 'dict', 0, _COHORT_A),
        CohortValues('dataset-b', 100, '1.33', 'dict', 0, _COHORT_B),
    ),
)


def test_min_metric_seeded_from_the_low_tail():
    updated, _ = suggest_mod.seed(CACHE, SPEC)
    metric = updated.metric('MEDIAN_COVERAGE')
    assert metric.fail == 11  # lowest p1 across cohorts (10.99), rounded for unit 'x'
    assert metric.warn == 15  # lowest p5 across cohorts (14.95)


def test_max_metric_seeded_from_the_high_tail():
    updated, _ = suggest_mod.seed(CACHE, SPEC)
    metric = updated.metric('FREEMIX')
    assert metric.fail == pytest.approx(0.1)  # highest p99 (0.09801), 2 dp for unit 'frac'
    assert metric.warn == pytest.approx(0.09)  # highest p95 (0.09405)


def test_relative_metric_seeds_fail_only():
    """Its warn tier is cohort-relative; an absolute warn alongside would double-flag."""
    updated, _ = suggest_mod.seed(CACHE, SPEC)
    assert updated.metric('DUP').fail is not None
    assert updated.metric('DUP').warn is None


def test_seeded_metrics_are_marked_unreviewed():
    updated, _ = suggest_mod.seed(CACHE, SPEC)
    assert updated.metric('MEDIAN_COVERAGE').reviewed is False


def test_seeded_metrics_record_their_evidence():
    updated, seeded = suggest_mod.seed(CACHE, SPEC)
    assert 'p1' in updated.metric('MEDIAN_COVERAGE').rationale
    assert {s.key for s in seeded} == {'MEDIAN_COVERAGE', 'FREEMIX', 'DUP'}


def test_seeded_evidence_names_the_cohort_count():
    _, seeded = suggest_mod.seed(CACHE, SPEC)
    assert '2 cohorts' in next(s for s in seeded if s.key == 'MEDIAN_COVERAGE').evidence


def test_reviewed_metric_is_never_overwritten():
    updated, seeded = suggest_mod.seed(CACHE, SPEC)
    assert updated.metric('PCT_20X') == SPEC.metric('PCT_20X')
    assert 'PCT_20X' not in {s.key for s in seeded}


def test_ungated_metric_is_never_seeded():
    updated, seeded = suggest_mod.seed(CACHE, SPEC)
    assert updated.metric('error_rate') == SPEC.metric('error_rate')
    assert 'error_rate' not in {s.key for s in seeded}


def test_seed_does_not_mutate_the_input_spec():
    suggest_mod.seed(CACHE, SPEC)
    assert SPEC.metric('MEDIAN_COVERAGE').warn is None
    assert SPEC.metric('MEDIAN_COVERAGE').fail == 1


def test_metric_with_no_data_is_reported_not_seeded():
    empty = ValueCache(
        seq_type='genome',
        generated='x',
        complete=True,
        metrics=('MEDIAN_COVERAGE',),
        cohorts=(CohortValues('dataset-a', 0, '1.33', 'dict', 0, {'MEDIAN_COVERAGE': []}),),
    )
    spec = spec_mod.loads(
        'seq_type = "genome"\ncache = "c.json"\n[metrics.MEDIAN_COVERAGE]\ndirection = "min"\nfail = 1\n',
    )
    updated, seeded = suggest_mod.seed(empty, spec)
    assert seeded == ()
    assert updated == spec


def test_metric_absent_from_one_cohort_still_seeds_from_the_others():
    """A cohort missing an un-gated metric shouldn't block seeding the rest."""
    partial = ValueCache(
        seq_type='genome',
        generated='x',
        complete=True,
        metrics=('MEDIAN_COVERAGE',),
        cohorts=(
            CohortValues('dataset-a', 100, '1.33', 'dict', 0, {'MEDIAN_COVERAGE': [float(v) for v in range(10, 110)]}),
            CohortValues('dataset-b', 0, '1.33', 'dict', 0, {'MEDIAN_COVERAGE': []}),
        ),
    )
    spec = spec_mod.loads(
        'seq_type = "genome"\ncache = "c.json"\n[metrics.MEDIAN_COVERAGE]\ndirection = "min"\nunit = "x"\nfail = 1\n',
    )
    updated, seeded = suggest_mod.seed(partial, spec)
    assert [s.key for s in seeded] == ['MEDIAN_COVERAGE']
    assert updated.metric('MEDIAN_COVERAGE').fail == 11


@pytest.mark.parametrize(
    ('unit', 'expected_fail'),
    [('x', 11), ('%', 11), ('frac', 10.99)],
)
def test_rounding_follows_the_unit(unit, expected_fail):
    """Integers read sensibly for coverage and percentages; fractions need 2 dp."""
    spec = spec_mod.loads(
        f'seq_type = "genome"\ncache = "c.json"\n[metrics.M]\ndirection = "min"\nunit = "{unit}"\nfail = 1\n',
    )
    cache = ValueCache(
        seq_type='genome',
        generated='x',
        complete=True,
        metrics=('M',),
        cohorts=(CohortValues('dataset-a', 100, '1.33', 'dict', 0, {'M': [float(v) for v in range(10, 110)]}),),
    )
    updated, _ = suggest_mod.seed(cache, spec)
    assert updated.metric('M').fail == pytest.approx(expected_fail)


def test_seeded_values_are_python_scalars_not_numpy():
    """These get written to TOML; numpy scalars would render wrong or raise."""
    updated, _ = suggest_mod.seed(CACHE, SPEC)
    fail = updated.metric('FREEMIX').fail
    assert type(fail) in (int, float)


def test_seed_result_survives_a_spec_round_trip():
    """with_metric re-validates, but the seeded spec must also dump and reload."""
    updated, _ = suggest_mod.seed(CACHE, SPEC)
    assert spec_mod.loads(spec_mod.dumps(updated)) == updated


def test_seed_preserves_a_hand_written_rationale():
    """An operator who annotates a still-unreviewed metric shouldn't lose it on re-run."""
    annotated = SPEC.with_metric('MEDIAN_COVERAGE', rationale='Operator note: legacy assay, keep loose for now.')
    updated, seeded = suggest_mod.seed(CACHE, annotated)
    assert updated.metric('MEDIAN_COVERAGE').rationale == 'Operator note: legacy assay, keep loose for now.'
    assert updated.metric('MEDIAN_COVERAGE').fail == 11  # thresholds still refresh
    assert 'MEDIAN_COVERAGE' in {s.key for s in seeded}


def test_seed_overwrites_its_own_previous_rationale():
    """A prior seed run's machine-written text refreshes rather than accumulating."""
    once, _ = suggest_mod.seed(CACHE, SPEC)
    twice, _ = suggest_mod.seed(CACHE, once)
    assert twice.metric('MEDIAN_COVERAGE').rationale == once.metric('MEDIAN_COVERAGE').rationale
    assert twice.metric('MEDIAN_COVERAGE').rationale.startswith(suggest_mod.RATIONALE_MARKER)
