"""Unit tests for cohort-relative (MAD) evaluation."""

import os
from collections.abc import Iterator

import pytest

from cpg_utils import config

from align_genotype.qc_calibration import relative as relative_mod
from align_genotype.qc_calibration import spec as spec_mod
from align_genotype.qc_calibration.cache import CohortValues, ValueCache
from align_genotype.qc_calibration.relative import CohortMad, MadEvaluation
from align_genotype.qc_calibration.stats import ChurnResult

METRIC = spec_mod.loads(
    'seq_type = "genome"\ncache = "c.json"\n'
    '[metrics.DUP]\ndirection = "max"\nunit = "%"\nfail = 1000.0\n'
    '[metrics.DUP.relative]\nk = 3.5\nmin_cohort = 5\n',
).metric('DUP')

# 49 evenly spread values plus one extreme: median 25.5, raw MAD 12.5,
# so the modified-z threshold lands at 25.5 + 3.5*12.5/0.6745 = 90.3629.
OUTLIER_COHORT = [*range(1, 50), 500.0]


@pytest.fixture(autouse=True)
def _restore_config_paths() -> Iterator[None]:
    """Keep the global cpg-utils config state from leaking between tests."""
    previous = os.environ.get('CPG_CONFIG_PATH', '')
    yield
    config.set_config_paths([p for p in previous.split(',') if p])


def _cache(**cohorts: list[float]) -> ValueCache:
    return ValueCache(
        seq_type='genome',
        generated='x',
        complete=True,
        metrics=('DUP',),
        cohorts=tuple(
            CohortValues(label, len(values), '1.33', 'dict', 0, {'DUP': [float(v) for v in values]})
            for label, values in cohorts.items()
        ),
    )


# --- per-cohort numbers -----------------------------------------------------------


def test_evaluate_reports_median_mad_and_threshold():
    evaluation = relative_mod.evaluate(_cache(**{'dataset-a': OUTLIER_COHORT}), METRIC, 'genome')
    cohort = evaluation.cohorts[0]
    assert cohort.label == 'dataset-a'
    assert cohort.n == 50
    assert cohort.median == pytest.approx(25.5)
    assert cohort.mad_raw == pytest.approx(12.5)
    assert cohort.threshold == pytest.approx(90.363, abs=0.001)


def test_warn_count_comes_from_the_production_relative_flags_path():
    evaluation = relative_mod.evaluate(_cache(**{'dataset-a': OUTLIER_COHORT}), METRIC, 'genome')
    assert evaluation.cohorts[0].n_warn == 1
    assert evaluation.cohorts[0].warn_rate == pytest.approx(0.02)


def test_cohort_below_min_cohort_is_skipped():
    evaluation = relative_mod.evaluate(_cache(**{'dataset-a': [1.0, 2.0, 3.0]}), METRIC, 'genome')
    cohort = evaluation.cohorts[0]
    assert cohort.skipped is not None
    assert 'min_cohort' in cohort.skipped
    assert cohort.n_warn == 0
    assert cohort.threshold is None


def test_degenerate_mad_cohort_is_skipped():
    evaluation = relative_mod.evaluate(_cache(**{'dataset-a': [5.0] * 10}), METRIC, 'genome')
    cohort = evaluation.cohorts[0]
    assert cohort.threshold is None
    assert 'zero MAD' in cohort.skipped
    assert cohort.n_warn == 0


def test_evaluate_restores_previous_config_paths(tmp_path):
    existing = tmp_path / 'existing.toml'
    existing.write_text('[workflow]\nsequencing_type = "genome"\n')
    config.set_config_paths([str(existing)])
    relative_mod.evaluate(_cache(**{'dataset-a': OUTLIER_COHORT}), METRIC, 'genome')
    assert config.get_config_paths() == [str(existing)]


def test_evaluate_restores_config_paths_even_on_failure(tmp_path, monkeypatch):
    """The context manager must not leak global config state when the body raises."""
    existing = tmp_path / 'existing.toml'
    existing.write_text('[workflow]\nsequencing_type = "genome"\n')
    config.set_config_paths([str(existing)])
    monkeypatch.setattr(relative_mod.check_multiqc, 'relative_flags', lambda *_a, **_k: 1 / 0)
    with pytest.raises(ZeroDivisionError):
        relative_mod.evaluate(_cache(**{'dataset-a': OUTLIER_COHORT}), METRIC, 'genome')
    assert config.get_config_paths() == [str(existing)]


def test_metric_without_a_relative_block_raises():
    plain = spec_mod.loads(
        'seq_type = "genome"\ncache = "c.json"\n[metrics.M]\ndirection = "max"\nfail = 1\n',
    ).metric('M')
    with pytest.raises(ValueError, match=r'no .* relative'):
        relative_mod.evaluate(_cache(**{'dataset-a': OUTLIER_COHORT}), plain, 'genome')


# --- churn --------------------------------------------------------------------


def test_evaluate_simulates_homogeneous_and_heterogeneous_growth():
    evaluation = relative_mod.evaluate(
        _cache(**{'dataset-a': OUTLIER_COHORT, 'dataset-b': [float(v) for v in range(200, 250)]}),
        METRIC,
        'genome',
    )
    assert {label for label, _ in evaluation.homogeneous} == {'dataset-a', 'dataset-b'}
    assert {(a, b) for a, b, _ in evaluation.heterogeneous} == {('dataset-a', 'dataset-b')}


def test_cohorts_too_small_are_excluded_from_churn():
    evaluation = relative_mod.evaluate(_cache(**{'dataset-a': [1.0, 2.0, 3.0]}), METRIC, 'genome')
    assert evaluation.homogeneous == ()
    assert evaluation.heterogeneous == ()


def test_heterogeneous_pairs_are_unordered_and_not_self_paired():
    evaluation = relative_mod.evaluate(
        _cache(
            **{
                'dataset-a': OUTLIER_COHORT,
                'dataset-b': [float(v) for v in range(200, 250)],
                'dataset-c': [float(v) for v in range(300, 350)],
            }
        ),
        METRIC,
        'genome',
    )
    pairs = {(a, b) for a, b, _ in evaluation.heterogeneous}
    assert pairs == {('dataset-a', 'dataset-b'), ('dataset-a', 'dataset-c'), ('dataset-b', 'dataset-c')}


# --- verdict ------------------------------------------------------------------


def _evaluation(warn_rates: list[float], flip_rates: list[float]) -> MadEvaluation:
    cohorts = tuple(
        CohortMad(
            f'cohort-{i}',
            n=1000,
            median=1.0,
            mad_raw=1.0,
            threshold=5.0,
            n_warn=round(rate * 1000),
            skipped=None,
        )
        for i, rate in enumerate(warn_rates)
    )
    churns = tuple(
        (f'cohort-{i}', ChurnResult(1.0, 1.0, 1000, 0, 0, round(rate * 1000))) for i, rate in enumerate(flip_rates)
    )
    return MadEvaluation('DUP', 'max', cohorts, churns, ())


def test_verdict_recommends_when_warn_and_churn_are_low():
    evaluation = _evaluation(warn_rates=[0.0, 0.042], flip_rates=[0.01, 0.02])
    assert evaluation.max_warn_rate == pytest.approx(0.042)
    assert evaluation.max_churn == pytest.approx(0.02)
    assert evaluation.verdict == 'RECOMMEND'
    assert evaluation.verdict_reason == ''


def test_verdict_rejects_on_high_churn():
    evaluation = _evaluation(warn_rates=[0.02], flip_rates=[0.245])
    assert evaluation.verdict == 'REJECT'
    assert 'churn' in evaluation.verdict_reason


def test_verdict_rejects_on_high_warn_rate():
    evaluation = _evaluation(warn_rates=[0.4], flip_rates=[0.0])
    assert evaluation.verdict == 'REJECT'
    assert 'warn' in evaluation.verdict_reason


def test_verdict_reports_both_reasons_when_both_fail():
    evaluation = _evaluation(warn_rates=[0.4], flip_rates=[0.245])
    assert evaluation.verdict == 'REJECT'
    assert 'warn' in evaluation.verdict_reason
    assert 'churn' in evaluation.verdict_reason


def test_skipped_cohorts_are_excluded_from_max_warn_rate():
    """A skipped cohort has n_warn 0, which would otherwise drag the peak down."""
    cohorts = (
        CohortMad('dataset-a', n=100, median=1.0, mad_raw=1.0, threshold=5.0, n_warn=40, skipped=None),
        CohortMad(
            'dataset-b',
            n=3,
            median=1.0,
            mad_raw=0.0,
            threshold=None,
            n_warn=0,
            skipped='cohort 3 < min_cohort 5',
        ),
    )
    evaluation = MadEvaluation('DUP', 'max', cohorts, (), ())
    assert evaluation.max_warn_rate == pytest.approx(0.4)


def test_max_churn_of_no_churn_data_is_zero():
    assert MadEvaluation('DUP', 'max', (), (), ()).max_churn == 0.0


def test_warn_count_tracks_the_k_written_into_the_production_config():
    """Proof the warn count really comes from the temp config, not a silent no-op.

    With k=3.5 only the 500 outlier clears the line; a k that small pulls the
    threshold down below most of the cohort. If `relative_flags` were reading an
    empty spec, both counts would be 0.
    """
    loose = relative_mod.evaluate(_cache(**{'dataset-a': OUTLIER_COHORT}), METRIC, 'genome')
    tight_metric = spec_mod.loads(
        'seq_type = "genome"\ncache = "c.json"\n'
        '[metrics.DUP]\ndirection = "max"\nunit = "%"\nfail = 1000.0\n'
        '[metrics.DUP.relative]\nk = 0.1\nmin_cohort = 5\n',
    ).metric('DUP')
    tight = relative_mod.evaluate(_cache(**{'dataset-a': OUTLIER_COHORT}), tight_metric, 'genome')
    assert loose.cohorts[0].n_warn == 1
    assert tight.cohorts[0].n_warn > 1


def test_evaluate_works_from_a_cold_start_with_no_config_paths_set(monkeypatch):
    """`get_config_paths()` raises when nothing was ever set - the restore must cope."""
    monkeypatch.delenv('CPG_CONFIG_PATH', raising=False)
    config.set_config_paths([])
    evaluation = relative_mod.evaluate(_cache(**{'dataset-a': OUTLIER_COHORT}), METRIC, 'genome')
    assert evaluation.cohorts[0].n_warn == 1
    with pytest.raises(config.ConfigError):
        config.get_config_paths()
