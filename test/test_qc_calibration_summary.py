"""Unit tests for the assembled calibration summary."""

import pytest

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

CURRENT = {'MEDIAN_COVERAGE': {'fail': 15, 'warn': 25}}


def dataset_values(name, coverage, dup, n_dropped=0) -> values_mod.DatasetValues:
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
def built():
    return summary_mod.build(
        [
            dataset_values('ds-a', [30, 32, 34, 36, 38, 10], [10, 10.5, 11, 11.5, 12, 40], n_dropped=2),
            dataset_values('ds-b', [40, 42, 44, 46, 48, 50], [7, 7.5, 8, 8.5, 9, 9.5]),
        ],
        SETTINGS,
        current=CURRENT,
        skipped_datasets=[{'dataset': 'ds-c', 'reason': 'no completed CramMultiQC qc analysis for genome'}],
        generated='2026-08-12T00:00:00',
        ar_guid='test-ar-guid',
    )


def test_records_run_level_provenance(built):
    assert built['sequencing_type'] == 'genome'
    assert built['ar_guid'] == 'test-ar-guid'
    assert built['settings']['k'] == pytest.approx(3.5)
    assert built['settings']['min_samples'] == 4
    assert [d['dataset'] for d in built['datasets']] == ['ds-a', 'ds-b']
    assert built['datasets'][0]['uri'] == 'gs://ds-a/multiqc_data.json'
    assert built['skipped_datasets'][0]['dataset'] == 'ds-c'


def test_counts_sequencing_groups_and_values_separately(built):
    metric = built['metrics']['MEDIAN_COVERAGE']
    assert metric['n_values'] == 12
    assert metric['n_groups_with_values'] == 12
    assert metric['n_datasets'] == 2
    assert metric['n_dropped'] == 2


def test_carries_the_shipped_thresholds_for_comparison(built):
    assert built['metrics']['MEDIAN_COVERAGE']['current'] == {'fail': 15, 'warn': 25}
    assert built['metrics']['dup_pct']['current'] == {}


def test_proposes_a_candidate_per_metric_with_data(built):
    assert built['metrics']['MEDIAN_COVERAGE']['candidate']['fail'] is not None
    assert built['metrics']['dup_pct']['candidate']['warn'] is None
    assert built['metrics']['ABSENT']['candidate'] is None


def test_scores_flag_rates_for_both_current_and_candidate(built):
    rates = built['metrics']['MEDIAN_COVERAGE']['flag_rates']
    # ds-a has one value of 10, below the shipped fail of 15.
    assert rates['current']['ds-a']['fail'] == pytest.approx(1 / 6)
    assert rates['current']['ds-b']['fail'] == pytest.approx(0.0)
    assert set(rates['candidate']) == {'ds-a', 'ds-b'}


def test_reports_percentiles_per_dataset(built):
    percentiles = built['metrics']['MEDIAN_COVERAGE']['percentiles']
    assert set(percentiles) == {'ds-a', 'ds-b'}
    assert percentiles['ds-b']['p50'] == pytest.approx(45.0)


def test_records_where_each_metric_was_found(built):
    assert built['metrics']['MEDIAN_COVERAGE']['present_in'] == ['picard_1']
    assert built['metrics']['MEDIAN_COVERAGE']['missing_from'] == []
    assert built['metrics']['ABSENT']['missing_from'] == ['ds-a', 'ds-b']


def test_a_metric_missing_everywhere_produces_a_loud_warning(built):
    assert any('every dataset' in w for w in built['warnings'])
    assert any('ABSENT' in w for w in built['warnings'])


def test_a_metric_present_everywhere_produces_no_warning(built):
    assert not any('MEDIAN_COVERAGE' in w for w in built['warnings'])


def test_evaluates_only_the_relative_metrics(built):
    assert set(built['relative']) == {'dup_pct'}
    evaluation = built['relative']['dup_pct']
    assert evaluation['verdict'] in ('RECOMMEND', 'REJECT')
    assert {d['dataset'] for d in evaluation['datasets']} == {'ds-a', 'ds-b'}


def test_relative_evaluation_reports_all_three_bars(built):
    evaluation = built['relative']['dup_pct']
    assert 'max_warn_rate' in evaluation
    assert 'max_growth_churn' in evaluation
    assert 'max_merge_churn' in evaluation
    assert evaluation['bars']['max_merge_churn'] == pytest.approx(0.05)


def test_includes_a_pasteable_config_snippet(built):
    assert '[qc_thresholds.genome.fail.min]' in built['config_snippet']


def test_is_json_serialisable(built):
    import json  # noqa: PLC0415

    json.dumps(built, allow_nan=False)


def test_no_datasets_at_all_still_builds(built):  # noqa: ARG001
    empty = summary_mod.build(
        [],
        SETTINGS,
        current={},
        skipped_datasets=[],
        generated='2026-08-12T00:00:00',
        ar_guid='x',
    )
    assert empty['datasets'] == []
    assert empty['metrics']['MEDIAN_COVERAGE']['candidate'] is None
    assert any('every dataset' in w for w in empty['warnings'])


# --- Tests added beyond the plan's list ------------------------------------------------
#
# The plan's `built` fixture never exercises several branches that matter for a report
# meant to inform a real threshold decision:
#
# * The shipped `warn` tier for MEDIAN_COVERAGE never fires without `fail` also firing
#   (the only sub-`fail` value in `ds-a` is also sub-`warn`), so `_rates`'s warn branch
#   with a genuinely non-zero rate was untested.
# * The relative metric's candidate always had `tiers.get('warn') is None`, but that was
#   only ever exercised on `built`, which never separately confirmed `flag_rates` maps it
#   to `None` rather than a computed-then-discarded number.
# * `missing_from` had only the all-or-nothing cases (present everywhere / absent
#   everywhere); the narrower "absent from *some* datasets" warning branch - explicitly
#   required to read differently from the "absent from every dataset" one - was never
#   triggered.
# * `DatasetMad.median`/`mad_raw` are `nan` only when a dataset has zero values for the
#   relative metric; `built` never has such a dataset, so the `nan -> None` guard in
#   `_relative_block` was never actually exercised, only trivially satisfied.
# * `built`'s datasets have 6 values each with `min_samples=4`; growth needs
#   `0.6 * n >= min_samples`, i.e. `n >= 7`, so growth churn is *always* empty for
#   `built` and the `_churn_result` mapping over `evaluation.growth` was never run.
# * `built` only has 2 datasets, so merge simulation only ever produces 2 ordered pairs -
#   far short of `MAX_PAIRS = 10` - so truncation and "worst first" sorting were never
#   checked.
# * `n_values` and `n_groups_with_values` were equal in every fixture, so the divergence
#   the module's docstring says the two fields exist to show was never actually produced.
# * The mid-task addition (coverage counts alongside each relative peak) needs its own
#   tests: a mix of usable/skipped datasets, and the "field present but zero" case for a
#   single-dataset run.


def test_flag_rate_reports_a_real_warn_only_breach():
    """A value below the shipped `warn` line but not the `fail` line must show up in `warn`.

    `ds-a` in the shared fixture only ever breaches `warn` via a value that also breaches
    `fail`, so the warn-only branch of `_rates` was never exercised by the plan's tests.
    """
    built = summary_mod.build(
        [dataset_values('ds-only', [30, 32, 34, 36, 38, 20], [10, 10.5, 11, 11.5, 12, 12.5])],
        SETTINGS,
        current=CURRENT,
        skipped_datasets=[],
        generated='2026-08-12T00:00:00',
        ar_guid='x',
    )
    rate = built['metrics']['MEDIAN_COVERAGE']['flag_rates']['current']['ds-only']
    assert rate['fail'] == pytest.approx(0.0)
    assert rate['warn'] == pytest.approx(1 / 6)


def test_relative_candidate_reports_no_warn_flag_rate(built):
    """A relative metric's candidate has no absolute `warn`; the flag rate must be `None`.

    Not `0.0` - `0.0` would claim "nothing warned", but nothing was scored at all.
    """
    rates = built['metrics']['dup_pct']['flag_rates']['candidate']['ds-a']
    assert rates['warn'] is None
    assert rates['fail'] is not None


def test_metric_missing_from_some_datasets_produces_a_narrower_warning():
    """Present in one dataset, absent in another - distinct from absent everywhere."""
    built = summary_mod.build(
        [
            # 6 points, not 3: a 3-point candidate here collapses fail and warn onto the
            # same rounded integer and trips the inert-tier warning too - see the
            # dedicated inert-tier tests below for that case in isolation.
            dataset_values('ds-p', [30, 32, 34, 36, 38, 10], [10, 10.5, 11, 11.5, 12, 12.5]),
            dataset_values('ds-q', [], [12, 12.5, 13]),
        ],
        SETTINGS,
        current=CURRENT,
        skipped_datasets=[],
        generated='2026-08-12T00:00:00',
        ar_guid='x',
    )
    assert built['metrics']['MEDIAN_COVERAGE']['missing_from'] == ['ds-q']
    narrow = [w for w in built['warnings'] if 'MEDIAN_COVERAGE' in w]
    assert len(narrow) == 1
    assert 'ds-q' in narrow[0]
    assert 'every dataset' not in narrow[0]


def test_relative_dataset_with_no_values_reports_none_for_median_and_mad():
    """A dataset with zero values for the relative metric has nan median/MAD - guarded to None."""
    built = summary_mod.build(
        [
            dataset_values('ds-r', [30, 32, 34, 36], [10, 10.5, 11, 11.5]),
            dataset_values('ds-s', [40, 42, 44, 46], []),
        ],
        SETTINGS,
        current=CURRENT,
        skipped_datasets=[],
        generated='2026-08-12T00:00:00',
        ar_guid='x',
    )
    entries = {d['dataset']: d for d in built['relative']['dup_pct']['datasets']}
    assert entries['ds-s']['median'] is None
    assert entries['ds-s']['mad'] is None
    assert entries['ds-s']['skipped'] is not None
    assert entries['ds-r']['median'] is not None
    assert entries['ds-r']['mad'] is not None


def test_relative_growth_churn_populates_for_large_enough_datasets():
    """Growth needs a 60% before-slice at least `min_samples` big; 6-value fixtures never clear that."""
    dup = [10, 10.5, 11, 11.5, 12, 12.5, 13, 13.5, 14, 14.5]  # median 12.25, MAD 1.25: not degenerate
    built = summary_mod.build(
        [dataset_values('ds-big', list(range(10)), dup)],
        SETTINGS,
        current=CURRENT,
        skipped_datasets=[],
        generated='2026-08-12T00:00:00',
        ar_guid='x',
    )
    growth = built['relative']['dup_pct']['growth']
    assert len(growth) == 1
    entry = growth[0]
    assert entry['dataset'] == 'ds-big'
    assert entry['ordered'] is not None
    assert isinstance(entry['ordered']['threshold_before'], float)
    assert isinstance(entry['ordered']['flip_rate'], float)
    assert isinstance(entry['worse'], float)
    assert isinstance(entry['ordering_sensitive'], bool)


def test_relative_merge_churn_truncates_and_reports_total_and_worst_first():
    """More than MAX_PAIRS ordered pairs must be truncated, worst first, with the true total kept."""
    datasets = [
        dataset_values(f'ds-{i}', list(range(6)), [10 + i * 0.3 + step * 0.5 for step in range(6)]) for i in range(4)
    ]
    built = summary_mod.build(
        datasets,
        SETTINGS,
        current=CURRENT,
        skipped_datasets=[],
        generated='2026-08-12T00:00:00',
        ar_guid='x',
    )
    evaluation = built['relative']['dup_pct']
    assert evaluation['merge_pairs_simulated'] == 12  # 4 datasets, n*(n-1) ordered pairs
    assert len(evaluation['merge_worst']) == 10  # truncated to MAX_PAIRS
    flip_rates = [entry['flip_rate'] for entry in evaluation['merge_worst']]
    assert flip_rates == sorted(flip_rates, reverse=True)


def test_metric_with_duplicated_entries_reports_the_divergence():
    """A metric split across two MultiQC sections carries more values than sequencing groups.

    The per-value rate must still use `n_values`, never `n_groups_with_values` - that is
    the whole point of keeping both fields.
    """
    duplicated = values_mod.DatasetValues(
        dataset='ds-dup',
        seq_type='genome',
        analysis_id=1,
        timestamp='2026-06-01T00:00:00',
        uri='gs://ds-dup/multiqc_data.json',
        multiqc_version='1.33',
        generated='2026-08-12T00:00:00',
        n_sequencing_groups=2,
        section_sizes={'picard_1': 2, 'picard_4': 1},
        metrics={
            'MEDIAN_COVERAGE': values_mod.MetricValues(
                entries=(
                    ('picard_1', 'ds-dup-CPG0', 30.0),
                    ('picard_1', 'ds-dup-CPG1', 32.0),
                    ('picard_4', 'ds-dup-CPG0', 31.0),  # same sequencing group, second section
                ),
                n_dropped=0,
            ),
            'dup_pct': values_mod.MetricValues(entries=(), n_dropped=0),
            'ABSENT': values_mod.MetricValues(entries=(), n_dropped=0),
        },
    )
    built = summary_mod.build(
        [duplicated],
        SETTINGS,
        current=CURRENT,
        skipped_datasets=[],
        generated='2026-08-12T00:00:00',
        ar_guid='x',
    )
    metric = built['metrics']['MEDIAN_COVERAGE']
    assert metric['n_values'] == 3
    assert metric['n_groups_with_values'] == 2
    assert metric['duplicated_in'] == ['ds-dup']
    rate = built['metrics']['MEDIAN_COVERAGE']['flag_rates']['current']['ds-dup']['fail']
    assert rate == pytest.approx(0.0)  # none of 30/32/31 breach fail=15; denominator would not change this either way


def test_no_datasets_relative_block_is_empty():
    built = summary_mod.build(
        [],
        SETTINGS,
        current={},
        skipped_datasets=[],
        generated='2026-08-12T00:00:00',
        ar_guid='x',
    )
    assert built['relative'] == {}


def test_relative_block_reports_dataset_and_simulation_coverage_counts():
    """The counts behind each peak figure must reflect the mix of usable and skipped datasets."""
    big_x = [10, 10.5, 11, 11.5, 12, 12.5, 13, 13.5, 14, 14.5]
    big_y = [10.2, 10.7, 11.2, 11.7, 12.2, 12.7, 13.2, 13.7, 14.2, 14.7]
    tiny_z = [50.0, 55.0]  # below min_samples=4: skipped before MAD is even considered
    built = summary_mod.build(
        [
            dataset_values('ds-x', list(range(10)), big_x),
            dataset_values('ds-y', list(range(10)), big_y),
            dataset_values('ds-z', list(range(2)), tiny_z),
        ],
        SETTINGS,
        current=CURRENT,
        skipped_datasets=[],
        generated='2026-08-12T00:00:00',
        ar_guid='x',
    )
    evaluation = built['relative']['dup_pct']
    assert evaluation['n_datasets'] == 3
    assert evaluation['n_datasets_skipped'] == 1
    assert evaluation['n_datasets_evaluated'] == 2
    assert evaluation['n_growth_evaluated'] == 2
    assert evaluation['n_merge_datasets_evaluated'] == 2
    assert evaluation['merge_pairs_simulated'] == 2


def test_relative_block_reports_zero_merge_pairs_for_a_single_dataset():
    """A single-dataset run must report `merge_pairs_simulated: 0`, not omit the field."""
    dup = [10, 10.5, 11, 11.5, 12, 12.5, 13, 13.5, 14, 14.5]
    built = summary_mod.build(
        [dataset_values('ds-solo', list(range(10)), dup)],
        SETTINGS,
        current=CURRENT,
        skipped_datasets=[],
        generated='2026-08-12T00:00:00',
        ar_guid='x',
    )
    evaluation = built['relative']['dup_pct']
    assert evaluation['merge_pairs_simulated'] == 0
    assert evaluation['n_merge_datasets_evaluated'] == 0
    assert evaluation['merge_worst'] == []


# --- Inert warn tier: a tier that looks configured and checks nothing ------------------
#
# Found by running the whole pipeline end to end, not by any test: `check_multiqc.worst_breach`
# evaluates fail before warn and returns the first breach, so a `min` metric's warn threshold
# must sit strictly above its fail threshold (the mirror for `max`) or every value that would
# warn is already recorded as a fail. Rounding for 'x'/'%' units can collapse two distinct raw
# percentiles onto the same integer and produce exactly this - `fail=27, warn=27` was observed
# in a real smoke run.


def test_candidate_warn_tier_equal_to_fail_produces_an_inert_warning():
    """p1=26.04 and p5=26.2 both round to 26 for a two-point 'x'-unit dataset - fail == warn."""
    built = summary_mod.build(
        [dataset_values('ds-only', [26.0, 30.0], [10, 11])],
        SETTINGS,
        current={},
        skipped_datasets=[],
        generated='2026-08-12T00:00:00',
        ar_guid='x',
    )
    candidate = built['metrics']['MEDIAN_COVERAGE']['candidate']
    assert candidate['fail'] == candidate['warn'] == 26
    inert = [w for w in built['warnings'] if 'unreachable' in w and 'MEDIAN_COVERAGE' in w]
    assert len(inert) == 1
    assert 'candidate' in inert[0]
    assert 'checks nothing' not in inert[0]  # exclusive to the absent-everywhere warning


def test_candidate_warn_at_or_past_fail_for_a_max_metric_produces_an_inert_warning():
    """The mirrored case: a 'max' metric is inert when warn >= fail (p99=29.96, p95=29.8, both round to 30)."""
    metric = settings_mod.MetricSpec(key='FREEMIX', direction='max', unit='%')
    max_settings = settings_mod.CalibrationSettings(seq_type='genome', metrics=(metric,), k=3.5, min_samples=4)
    values = values_mod.DatasetValues(
        dataset='ds-only',
        seq_type='genome',
        analysis_id=1,
        timestamp='2026-06-01T00:00:00',
        uri='gs://ds-only/multiqc_data.json',
        multiqc_version='1.33',
        generated='2026-08-12T00:00:00',
        n_sequencing_groups=2,
        section_sizes={'picard_1': 2},
        metrics={
            'FREEMIX': values_mod.MetricValues(
                entries=(('picard_1', 'ds-only-CPG0', 26.0), ('picard_1', 'ds-only-CPG1', 30.0)),
                n_dropped=0,
            ),
        },
    )
    built = summary_mod.build(
        [values],
        max_settings,
        current={},
        skipped_datasets=[],
        generated='2026-08-12T00:00:00',
        ar_guid='x',
    )
    candidate = built['metrics']['FREEMIX']['candidate']
    assert candidate['fail'] == candidate['warn'] == 30
    inert = [w for w in built['warnings'] if 'unreachable' in w and 'FREEMIX' in w]
    assert len(inert) == 1
    assert 'candidate' in inert[0]


def test_correctly_ordered_tiers_produce_no_inert_warning(built):
    """MEDIAN_COVERAGE's candidate in the shared fixture has a real gap between fail and warn."""
    candidate = built['metrics']['MEDIAN_COVERAGE']['candidate']
    assert candidate['fail'] != candidate['warn']
    assert not any('unreachable' in w for w in built['warnings'])


def test_inert_current_tier_produces_the_shipped_config_variant_of_the_warning():
    """An inert pair already in `current` is a live production defect, not just a proposal to review."""
    built = summary_mod.build(
        [dataset_values('ds-only', [30, 32, 34, 36, 38, 10], [10, 10.5, 11, 11.5, 12, 12.5])],
        SETTINGS,
        current={'MEDIAN_COVERAGE': {'fail': 27, 'warn': 27}},
        skipped_datasets=[],
        generated='2026-08-12T00:00:00',
        ar_guid='x',
    )
    # The candidate for this fixture is fail=11, warn=15 (a real gap) - only the shipped pair
    # is inert, so exactly one warning should fire, and it must be the shipped-config variant.
    inert = [w for w in built['warnings'] if 'unreachable' in w and 'MEDIAN_COVERAGE' in w]
    assert len(inert) == 1
    assert 'shipped' in inert[0]
    assert 'config_template.toml' in inert[0]
    assert 'checks nothing' not in inert[0]
