"""Golden test for the emitted config block."""

import os
from collections.abc import Iterator

import pytest

from cpg_utils import config

from align_genotype.qc_calibration import emit as emit_mod
from align_genotype.qc_calibration import spec as spec_mod
from align_genotype.qc_calibration import suggest as suggest_mod
from align_genotype.qc_calibration import tomlio
from align_genotype.qc_calibration.cache import CohortValues, ValueCache
from align_genotype.qc_calibration.emit import EmitError
from align_genotype.scripts import check_multiqc

SPEC = spec_mod.loads(
    'seq_type = "genome"\ncache = "c.json"\n'
    '\n[metrics.MEDIAN_COVERAGE]\ndirection = "min"\nunit = "x"\nfail = 15\nwarn = 25\nreviewed = true\n'
    'rationale = "Primary depth gate."\n'
    '\n[metrics.FREEMIX]\ndirection = "max"\nunit = "frac"\nfail = 0.04\nwarn = 0.02\nreviewed = true\n'
    'rationale = "Contamination safety gate."\n'
    '\n[metrics.DUP]\ndirection = "max"\nunit = "%"\nfail = 40\nreviewed = true\n'
    'rationale = "Library-prep dependent."\n[metrics.DUP.relative]\nk = 3.5\nmin_cohort = 50\n'
    '\n[metrics.error_rate]\ndirection = "max"\nunit = "frac"\ngated = false\n',
)

CACHE = ValueCache(
    seq_type='genome',
    generated='2026-08-11T00:00:00',
    complete=True,
    metrics=('MEDIAN_COVERAGE', 'FREEMIX', 'DUP', 'error_rate'),
    cohorts=(
        CohortValues(
            'dataset-a',
            5,
            '1.33',
            'dict',
            0,
            {
                'MEDIAN_COVERAGE': [10.0, 20.0, 30.0, 40.0, 50.0],
                'FREEMIX': [0.001] * 5,
                'DUP': [5.0, 6.0, 7.0, 8.0, 9.0],
                'error_rate': [0.01] * 5,
            },
        ),
    ),
)


@pytest.fixture(autouse=True)
def _restore_config_paths() -> Iterator[None]:
    """Keep the global cpg-utils config state from leaking between tests."""
    previous = os.environ.get('CPG_CONFIG_PATH', '')
    yield
    config.set_config_paths([p for p in previous.split(',') if p])


def _rendered() -> str:
    return emit_mod.render(SPEC, CACHE, spec_path='calibration/spec.genome.toml', generated='2026-08-11')


# --- the golden assertion: the parsed block, exactly --------------------------


def test_emitted_block_parses_to_exactly_the_expected_config():
    assert tomlio.loads(_rendered()) == {
        'qc_thresholds': {
            'genome': {
                'fail': {'min': {'MEDIAN_COVERAGE': 15}, 'max': {'FREEMIX': 0.04, 'DUP': 40}},
                'warn': {'min': {'MEDIAN_COVERAGE': 25}, 'max': {'FREEMIX': 0.02}},
                'relative': {'DUP': {'direction': 'max', 'k': 3.5, 'min_cohort': 50}},
            },
        },
    }


def test_ungated_metric_is_not_emitted():
    assert 'error_rate' not in _rendered()


def test_relative_metric_has_no_absolute_warn_entry():
    warn_max = tomlio.loads(_rendered())['qc_thresholds']['genome']['warn']['max']
    assert 'DUP' not in warn_max


def test_round_trips_through_production_load_thresholds(tmp_path):
    """The block must be loadable by the code that will actually enforce it."""
    path = tmp_path / 'emitted.toml'
    path.write_text(_rendered())
    config.set_config_paths([str(path)])
    assert check_multiqc.load_thresholds('genome') == {
        'min': {'MEDIAN_COVERAGE': {'fail': 15, 'warn': 25}},
        'max': {'FREEMIX': {'fail': 0.04, 'warn': 0.02}, 'DUP': {'fail': 40}},
    }


def test_relative_block_is_readable_by_production_relative_flags(tmp_path):
    """The relative table must match the shape relative_flags expects."""
    path = tmp_path / 'emitted.toml'
    path.write_text(_rendered())
    config.set_config_paths([str(path)])
    spec = config.config_retrieve(['qc_thresholds', 'genome', 'relative'])
    assert spec == {'DUP': {'direction': 'max', 'k': 3.5, 'min_cohort': 50}}


# --- structure and provenance --------------------------------------------------


def test_sections_appear_in_config_template_order():
    text = _rendered()
    order = [
        '[qc_thresholds.genome.fail.min]',
        '[qc_thresholds.genome.fail.max]',
        '[qc_thresholds.genome.warn.min]',
        '[qc_thresholds.genome.warn.max]',
        '[qc_thresholds.genome.relative.DUP]',
    ]
    positions = [text.index(section) for section in order]
    assert positions == sorted(positions)


def test_header_records_provenance_without_naming_datasets():
    """Cohort labels are real dataset names; this block gets committed and pushed."""
    text = _rendered()
    assert '1 cohort' in text
    assert 'calibration/spec.genome.toml' in text
    assert 'dataset-a' not in text


def test_no_cohort_label_appears_anywhere_in_the_output():
    """Belt and braces on the hard rule - check every label in the cache."""
    text = _rendered()
    for label in CACHE.labels:
        assert label not in text


# --- composed with `suggest`: the actual documented workflow ------------------
#
# `test_no_cohort_label_appears_anywhere_in_the_output` above uses a hand-written spec
# with prose rationales, so it can't catch a label `suggest` itself writes into
# `rationale`. The workflow an operator actually runs is seed -> review -> emit; only a
# test that runs all three proves the composition is safe.

_TWO_COHORT_CACHE = ValueCache(
    seq_type='genome',
    generated='x',
    complete=True,
    metrics=('MEDIAN_COVERAGE',),
    cohorts=(
        CohortValues('dataset-a', 100, '1.33', 'dict', 0, {'MEDIAN_COVERAGE': [float(v) for v in range(10, 110)]}),
        CohortValues('dataset-b', 100, '1.33', 'dict', 0, {'MEDIAN_COVERAGE': [float(v) for v in range(20, 120)]}),
    ),
)


def test_seed_then_review_then_emit_never_leaks_a_cohort_label():
    """The documented workflow: seed, operator reviews, emit. No label survives it."""
    raw_spec = spec_mod.loads(
        'seq_type = "genome"\ncache = "c.json"\n[metrics.MEDIAN_COVERAGE]\ndirection = "min"\nunit = "x"\nfail = 1\n',
    )
    updated, seeded = suggest_mod.seed(_TWO_COHORT_CACHE, raw_spec)
    assert seeded  # sanity: the seed actually produced something to review
    for result in seeded:
        updated = updated.with_metric(result.key, reviewed=True)
    text = emit_mod.render(updated, _TWO_COHORT_CACHE, generated='2026-08-11')
    for label in _TWO_COHORT_CACHE.labels:
        assert label not in text


def test_render_rejects_a_hand_written_rationale_naming_a_cohort():
    """The backstop: a label reaching `rationale` by hand-edit must still be caught."""
    tainted = SPEC.with_metric('MEDIAN_COVERAGE', rationale='Elevated in dataset-a, see lab notes.')
    with pytest.raises(EmitError, match=r"metric 'MEDIAN_COVERAGE'.*dataset-a"):
        emit_mod.render(tainted, CACHE, generated='2026-08-11')


def test_metric_comment_carries_rationale_and_generated_evidence():
    text = _rendered()
    assert '# Primary depth gate.' in text
    assert 'Cohort medians 30.0-30.0' in text
    assert '20%' in text  # 1 of 5 samples below the fail line of 15


def test_comment_lines_stay_within_the_line_length():
    for line in _rendered().splitlines():
        assert len(line) <= 120


def test_empty_section_is_omitted():
    spec = spec_mod.loads(
        'seq_type = "genome"\ncache = "c.json"\n'
        '[metrics.MEDIAN_COVERAGE]\ndirection = "min"\nunit = "x"\nfail = 15\nreviewed = true\n',
    )
    text = emit_mod.render(spec, CACHE, generated='2026-08-11')
    assert '[qc_thresholds.genome.fail.min]' in text
    assert 'warn' not in text.split('[qc_thresholds.genome.fail.min]')[1]


def test_manifest_path_cited_when_given():
    text = emit_mod.render(SPEC, CACHE, manifest_path='calibration/manifest.genome.toml', generated='2026-08-11')
    assert 'calibration/manifest.genome.toml' in text


# --- the review gate -----------------------------------------------------------


def test_refuses_to_emit_an_unreviewed_gated_metric():
    unreviewed = SPEC.with_metric(replace_key='MEDIAN_COVERAGE', reviewed=False)
    with pytest.raises(EmitError, match=r"unreviewed.*\['MEDIAN_COVERAGE'\]"):
        emit_mod.render(unreviewed, CACHE, generated='2026-08-11')


def test_refusal_names_every_unreviewed_metric():
    unreviewed = SPEC.with_metric(replace_key='MEDIAN_COVERAGE', reviewed=False).with_metric(
        replace_key='FREEMIX',
        reviewed=False,
    )
    with pytest.raises(EmitError, match='FREEMIX'):
        emit_mod.render(unreviewed, CACHE, generated='2026-08-11')


def test_unreviewed_ungated_metric_does_not_block_emission():
    """error_rate is gated = false and unreviewed; only gated metrics matter."""
    assert 'qc_thresholds' in _rendered()
