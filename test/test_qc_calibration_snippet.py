"""Unit tests for the copy-pasteable qc_thresholds block."""

import sys
from typing import Any

import numpy as np
import pytest

from align_genotype.qc_calibration import settings as settings_mod
from align_genotype.qc_calibration import snippet as snippet_mod
from align_genotype.qc_calibration.thresholds import Candidate

if sys.version_info >= (3, 11):
    import tomllib
else:
    import tomli as tomllib

SETTINGS = settings_mod.CalibrationSettings(
    seq_type='genome',
    metrics=(
        settings_mod.MetricSpec(key='MEDIAN_COVERAGE', direction='min', unit='x'),
        settings_mod.MetricSpec(key='FREEMIX', direction='max', unit='frac'),
        settings_mod.MetricSpec(key='reads_duplicated_percent', direction='max', unit='%', relative=True),
    ),
)

CANDIDATES = {
    'MEDIAN_COVERAGE': Candidate('MEDIAN_COVERAGE', fail=17, warn=26, basis='...'),
    'FREEMIX': Candidate('FREEMIX', fail=0.03, warn=0.015, basis='...'),
    'reads_duplicated_percent': Candidate('reads_duplicated_percent', fail=38, warn=None, basis='...'),
}


def test_block_parses_as_toml():
    parsed = tomllib.loads(snippet_mod.render(SETTINGS, CANDIDATES))
    assert parsed['qc_thresholds']['genome']['fail']['min'] == {'MEDIAN_COVERAGE': 17}
    assert parsed['qc_thresholds']['genome']['fail']['max'] == {
        'FREEMIX': 0.03,
        'reads_duplicated_percent': 38,
    }
    assert parsed['qc_thresholds']['genome']['warn']['min'] == {'MEDIAN_COVERAGE': 26}
    assert parsed['qc_thresholds']['genome']['warn']['max'] == {'FREEMIX': 0.015}


def test_relative_metric_has_no_absolute_warn_entry():
    parsed = tomllib.loads(snippet_mod.render(SETTINGS, CANDIDATES))
    assert 'reads_duplicated_percent' not in parsed['qc_thresholds']['genome']['warn']['max']


def test_relative_table_carries_direction_k_and_min_samples():
    parsed = tomllib.loads(snippet_mod.render(SETTINGS, CANDIDATES))
    relative = parsed['qc_thresholds']['genome']['relative']['reads_duplicated_percent']
    assert relative == {'direction': 'max', 'k': 3.5, 'min_samples': 50}


def test_block_round_trips_through_the_production_loader(monkeypatch):
    """Proof the emitted shape is what load_thresholds actually reads."""
    from align_genotype.scripts import check_multiqc  # noqa: PLC0415

    parsed = tomllib.loads(snippet_mod.render(SETTINGS, CANDIDATES))
    monkeypatch.setattr(
        check_multiqc.config,
        'config_retrieve',
        lambda keys, default=None: _dig(parsed, keys, default),
    )
    loaded = check_multiqc.load_thresholds('genome')
    assert loaded['min']['MEDIAN_COVERAGE'] == {'fail': 17, 'warn': 26}
    assert loaded['max']['reads_duplicated_percent'] == {'fail': 38}


def _dig(tree: dict, keys: list, default: Any) -> Any:
    node = tree
    for key in keys:
        if not isinstance(node, dict) or key not in node:
            return default
        node = node[key]
    return node


def test_section_order_matches_the_committed_config():
    """So the block diffs against config_template.toml as a change, not a rewrite."""
    body = snippet_mod.render(SETTINGS, CANDIDATES)
    order = [line for line in body.splitlines() if line.startswith('[qc_thresholds')]
    assert order == [
        '[qc_thresholds.genome.fail.min]',
        '[qc_thresholds.genome.fail.max]',
        '[qc_thresholds.genome.warn.min]',
        '[qc_thresholds.genome.warn.max]',
        '[qc_thresholds.genome.relative.reads_duplicated_percent]',
    ]


def test_a_metric_with_no_candidate_is_omitted():
    body = snippet_mod.render(SETTINGS, {'MEDIAN_COVERAGE': CANDIDATES['MEDIAN_COVERAGE']})
    assert 'FREEMIX' not in body


def test_empty_candidates_render_an_empty_block():
    assert snippet_mod.render(SETTINGS, {}).strip() == ''


def test_metric_keys_are_quoted_like_the_committed_config():
    assert '"MEDIAN_COVERAGE" = 17' in snippet_mod.render(SETTINGS, CANDIDATES)


def test_integers_render_without_a_decimal_point():
    """Pasting must not turn a shipped `= 15` into `= 15.0`."""
    body = snippet_mod.render(SETTINGS, CANDIDATES)
    assert '= 17' in body
    assert '= 17.0' not in body


def test_exome_seq_type_is_honoured():
    settings = settings_mod.CalibrationSettings(
        seq_type='exome',
        metrics=(settings_mod.MetricSpec(key='FREEMIX', direction='max', unit='frac'),),
    )
    body = snippet_mod.render(settings, {'FREEMIX': Candidate('FREEMIX', fail=0.03, warn=None, basis='')})
    assert '[qc_thresholds.exome.fail.max]' in body


# --- Additional tests beyond the plan's list ---
#
# The plan's fixtures only exercise a fail-only metric (FREEMIX has both, MEDIAN_COVERAGE
# has both, reads_duplicated_percent is relative-so-warn-is-None). None of them cover a
# *non-relative* metric that has a candidate but is missing one tier - the walrus-filtered
# branch in `render`'s list comprehension when `found` exists but `getattr(found, severity)`
# is None for a reason other than "this metric is relative". Nor do they cover the
# numpy-scalar or bool branches of `_fmt`. `_fmt` is exercised only through `render` (the
# public entry point), rather than by importing the private helper directly.


def test_a_non_relative_warn_only_candidate_is_omitted_from_fail():
    """A metric with no fail candidate (e.g. a correlated warn-only signal) skips fail."""
    settings = settings_mod.CalibrationSettings(
        seq_type='genome',
        metrics=(settings_mod.MetricSpec(key='reads_properly_paired_percent', direction='min', unit='%'),),
    )
    candidates = {
        'reads_properly_paired_percent': Candidate('reads_properly_paired_percent', fail=None, warn=92, basis=''),
    }
    body = snippet_mod.render(settings, candidates)
    assert '[qc_thresholds.genome.fail.min]' not in body
    assert '[qc_thresholds.genome.warn.min]' in body
    assert '"reads_properly_paired_percent" = 92' in body


def test_a_non_relative_fail_only_candidate_is_omitted_from_warn():
    """A metric with no warn candidate skips the warn section, independent of `relative`."""
    settings = settings_mod.CalibrationSettings(
        seq_type='genome',
        metrics=(settings_mod.MetricSpec(key='reads_mapped_percent', direction='min', unit='%'),),
    )
    candidates = {'reads_mapped_percent': Candidate('reads_mapped_percent', fail=80, warn=None, basis='')}
    body = snippet_mod.render(settings, candidates)
    assert '[qc_thresholds.genome.fail.min]' in body
    assert '[qc_thresholds.genome.warn.min]' not in body


def test_a_numpy_float_threshold_renders_without_the_np_repr_wrapper():
    """`repr(np.float64(0.03))` is `np.float64(0.03)`, not valid TOML - must be unwrapped."""
    settings = settings_mod.CalibrationSettings(
        seq_type='genome',
        metrics=(settings_mod.MetricSpec(key='FREEMIX', direction='max', unit='frac'),),
    )
    candidates = {'FREEMIX': Candidate('FREEMIX', fail=np.float64(0.03), warn=None, basis='')}
    body = snippet_mod.render(settings, candidates)
    assert '"FREEMIX" = 0.03' in body
    assert 'np.float64' not in body


def test_a_numpy_int_threshold_renders_without_the_np_repr_wrapper():
    settings = settings_mod.CalibrationSettings(
        seq_type='genome',
        metrics=(settings_mod.MetricSpec(key='MEDIAN_COVERAGE', direction='min', unit='x'),),
    )
    candidates = {'MEDIAN_COVERAGE': Candidate('MEDIAN_COVERAGE', fail=np.int64(17), warn=None, basis='')}
    body = snippet_mod.render(settings, candidates)
    assert '"MEDIAN_COVERAGE" = 17' in body
    assert 'np.int64' not in body


def test_a_bool_threshold_is_rejected_rather_than_rendered_as_true_false():
    """`bool` subclasses `int`; a threshold must never silently become a TOML boolean."""
    settings = settings_mod.CalibrationSettings(
        seq_type='genome',
        metrics=(settings_mod.MetricSpec(key='FREEMIX', direction='max', unit='frac'),),
    )
    candidates = {'FREEMIX': Candidate('FREEMIX', fail=True, warn=None, basis='')}
    with pytest.raises(TypeError):
        snippet_mod.render(settings, candidates)
