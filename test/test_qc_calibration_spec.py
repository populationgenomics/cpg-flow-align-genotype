"""Unit tests for the calibration spec: load, validate, dump, round-trip."""

import pytest

from align_genotype.qc_calibration import spec as spec_mod
from align_genotype.qc_calibration.spec import CalibrationSpec, MetricSpec, RelativeSpec, SpecError

VALID = """
seq_type = "genome"
cache = "calibration/genome_values.json"

[metrics.MEDIAN_COVERAGE]
direction = "min"
unit = "x"
fail = 15
warn = 25
reviewed = true
rationale = "Primary depth gate."

[metrics.reads_duplicated_percent]
direction = "max"
unit = "%"
fail = 40
reviewed = true
[metrics.reads_duplicated_percent.relative]
k = 3.5
min_cohort = 50

[metrics.error_rate]
direction = "max"
unit = "frac"
gated = false
rationale = "Varies with ancestry; not a defensible hard gate."
"""


def test_loads_valid_spec():
    spec = spec_mod.loads(VALID)
    assert spec.seq_type == 'genome'
    assert spec.cache == 'calibration/genome_values.json'
    assert [m.key for m in spec.metrics] == ['MEDIAN_COVERAGE', 'reads_duplicated_percent', 'error_rate']


def test_metric_defaults():
    depth = spec_mod.loads(VALID).metric('MEDIAN_COVERAGE')
    assert depth == MetricSpec(
        key='MEDIAN_COVERAGE',
        direction='min',
        unit='x',
        gated=True,
        fail=15,
        warn=25,
        relative=None,
        reviewed=True,
        rationale='Primary depth gate.',
    )


def test_relative_block_parsed():
    dup = spec_mod.loads(VALID).metric('reads_duplicated_percent')
    assert dup.relative == RelativeSpec(k=3.5, min_cohort=50)
    assert dup.warn is None


def test_gated_property_excludes_ungated():
    assert [m.key for m in spec_mod.loads(VALID).gated] == ['MEDIAN_COVERAGE', 'reads_duplicated_percent']


def test_metric_keys_property():
    assert spec_mod.loads(VALID).metric_keys == ('MEDIAN_COVERAGE', 'reads_duplicated_percent', 'error_rate')


def test_unknown_metric_lookup_raises():
    with pytest.raises(KeyError, match='NOPE'):
        spec_mod.loads(VALID).metric('NOPE')


# --- validation ---------------------------------------------------------------


@pytest.mark.parametrize(
    ('body', 'match'),
    [
        ('[metrics.M]\ndirection = "sideways"\nfail = 1\n', "direction must be 'min' or 'max'"),
        ('[metrics.M]\nfail = 1\n', 'missing required key: direction'),
        ('[metrics.M]\ndirection = "min"\nunit = "furlongs"\nfail = 1\n', 'unit must be one of'),
        ('[metrics.M]\ndirection = "min"\n', 'gated metric defines no fail, warn or relative tier'),
        (
            '[metrics.M]\ndirection = "max"\n[metrics.M.relative]\nk = 3.5\n',
            'relative tier requires an absolute fail',
        ),
        ('[metrics.M]\ndirection = "max"\ngated = false\nfail = 1\n', 'non-gated metric must not define'),
        ('[metrics.M]\ndirection = "max"\ngated = false\nwarn = 1\n', 'non-gated metric must not define'),
        (
            '[metrics.M]\ndirection = "max"\ngated = false\n[metrics.M.relative]\nk = 3.5\n',
            'non-gated metric must not define',
        ),
    ],
)
def test_validation_errors(body, match):
    text = f'seq_type = "genome"\ncache = "c.json"\n\n{body}'
    with pytest.raises(SpecError, match=match):
        spec_mod.loads(text)


def test_relative_with_absolute_warn_rejected():
    text = (
        'seq_type = "genome"\ncache = "c.json"\n\n'
        '[metrics.M]\ndirection = "max"\nfail = 40\nwarn = 30\n'
        '[metrics.M.relative]\nk = 3.5\n'
    )
    with pytest.raises(SpecError, match='relative tier replaces the absolute warn'):
        spec_mod.loads(text)


def test_missing_seq_type_rejected():
    with pytest.raises(SpecError, match='missing required key: seq_type'):
        spec_mod.loads('cache = "c.json"\n[metrics.M]\ndirection = "min"\nfail = 1\n')


def test_missing_cache_rejected():
    with pytest.raises(SpecError, match='missing required key: cache'):
        spec_mod.loads('seq_type = "genome"\n[metrics.M]\ndirection = "min"\nfail = 1\n')


def test_no_metrics_rejected():
    with pytest.raises(SpecError, match='defines no metrics'):
        spec_mod.loads('seq_type = "genome"\ncache = "c.json"\n')


@pytest.mark.parametrize(
    'text',
    [
        'seq_type = "genome"\ncache = "c.json"\n\n[[metrics]]\ndirection = "min"\nfail = 1\n',
        'seq_type = "genome"\ncache = "c.json"\nmetrics = "oops"\n',
    ],
)
def test_non_table_metrics_rejected(text):
    """A hand-edit to `[[metrics]]` (an array of tables) is a plausible typo for a table of tables."""
    with pytest.raises(SpecError, match=r'\[metrics\] must be a table'):
        spec_mod.loads(text)


def test_metric_key_needing_quoting_rejected():
    """dumps() writes metric keys unquoted, so a key that needs quoting can't round-trip."""
    text = 'seq_type = "genome"\ncache = "c.json"\n\n[metrics."odd.key"]\ndirection = "min"\nfail = 1\n'
    with pytest.raises(SpecError, match='not a bare TOML key'):
        spec_mod.loads(text)


def test_unknown_metric_key_rejected():
    text = 'seq_type = "genome"\ncache = "c.json"\n\n[metrics.M]\ndirection = "min"\nfail = 1\nnope = 2\n'
    with pytest.raises(SpecError, match=r"unknown key\(s\) \['nope'\]"):
        spec_mod.loads(text)


def test_unknown_relative_key_rejected():
    text = (
        'seq_type = "genome"\ncache = "c.json"\n\n'
        '[metrics.M]\ndirection = "max"\nfail = 1\n[metrics.M.relative]\nk = 3.5\nnope = 2\n'
    )
    with pytest.raises(SpecError, match=r"unknown relative key\(s\) \['nope'\]"):
        spec_mod.loads(text)


def test_relative_defaults():
    text = 'seq_type = "genome"\ncache = "c.json"\n\n[metrics.M]\ndirection = "max"\nfail = 1\n[metrics.M.relative]\n'
    assert spec_mod.loads(text).metric('M').relative == RelativeSpec(k=3.5, min_cohort=50)


@pytest.mark.parametrize('field', ['reviewed', 'gated'])
def test_quoted_boolean_rejected(field):
    """A quoted 'false' is a plausible typo; bool() coercion would silently flip it true."""
    text = f'seq_type = "genome"\ncache = "c.json"\n\n[metrics.M]\ndirection = "max"\nfail = 1\n{field} = "false"\n'
    with pytest.raises(SpecError, match=f'{field} must be true or false'):
        spec_mod.loads(text)


@pytest.mark.parametrize('field', ['fail', 'warn'])
def test_string_threshold_rejected(field):
    text = f'seq_type = "genome"\ncache = "c.json"\n\n[metrics.M]\ndirection = "max"\n{field} = "high"\n'
    with pytest.raises(SpecError, match=f'{field} must be a number'):
        spec_mod.loads(text)


@pytest.mark.parametrize('field', ['fail', 'warn'])
def test_boolean_threshold_rejected(field):
    """bool subclasses int, so this must be excluded explicitly or `fail = true` becomes 1."""
    text = f'seq_type = "genome"\ncache = "c.json"\n\n[metrics.M]\ndirection = "max"\n{field} = true\n'
    with pytest.raises(SpecError, match=f'{field} must be a number'):
        spec_mod.loads(text)


@pytest.mark.parametrize(
    'relative_body',
    [
        'k = "big"',
        'min_cohort = "many"',
        'k = [1, 2]',
    ],
)
def test_relative_non_numeric_rejected(relative_body):
    text = (
        'seq_type = "genome"\ncache = "c.json"\n\n'
        f'[metrics.M]\ndirection = "max"\nfail = 1\n[metrics.M.relative]\n{relative_body}\n'
    )
    with pytest.raises(SpecError, match='must be numeric'):
        spec_mod.loads(text)


# --- dump / round-trip ----------------------------------------------------------


def test_dump_round_trips():
    original = spec_mod.loads(VALID)
    assert spec_mod.loads(spec_mod.dumps(original)) == original


def test_dump_places_relative_subtable_after_its_parent_keys():
    text = spec_mod.dumps(spec_mod.loads(VALID))
    parent = text.index('[metrics.reads_duplicated_percent]')
    sub = text.index('[metrics.reads_duplicated_percent.relative]')
    nxt = text.index('[metrics.error_rate]')
    assert parent < sub < nxt


def test_with_metric_returns_a_new_spec():
    """Specs are immutable; suggest() builds an updated copy."""
    original = spec_mod.loads(VALID)
    updated = original.with_metric(
        replace_key='MEDIAN_COVERAGE',
        fail=20,
        warn=30,
        reviewed=False,
        rationale='seeded',
    )
    assert original.metric('MEDIAN_COVERAGE').fail == 15
    assert updated.metric('MEDIAN_COVERAGE').fail == 20
    assert updated.metric('MEDIAN_COVERAGE').reviewed is False
    assert [m.key for m in updated.metrics] == [m.key for m in original.metrics]


def test_with_metric_unknown_key_raises():
    """metric() raises KeyError for a typo'd key; with_metric must agree, not silently no-op."""
    original = spec_mod.loads(VALID)
    with pytest.raises(KeyError, match='NOPE'):
        original.with_metric(replace_key='NOPE', fail=1)


def test_with_metric_rename_blocked():
    original = spec_mod.loads(VALID)
    with pytest.raises(SpecError, match='cannot rename'):
        original.with_metric(replace_key='MEDIAN_COVERAGE', key='RENAMED')


def test_with_metric_validates_the_result():
    """A change that would violate rule 3 (relative forbids warn) must raise, not build silently."""
    original = spec_mod.loads(VALID)
    with pytest.raises(SpecError, match='relative tier replaces the absolute warn'):
        original.with_metric(replace_key='reads_duplicated_percent', warn=30)


def test_with_metric_result_is_valid():
    """The updated spec must itself be valid - i.e. usable and round-trippable, not just built."""
    original = spec_mod.loads(VALID)
    updated = original.with_metric(replace_key='error_rate', gated=True, fail=0.05, reviewed=True)
    assert spec_mod.loads(spec_mod.dumps(updated)) == updated


def test_dumps_on_directly_constructed_spec():
    """suggest() builds specs directly rather than via loads(); dumps must handle that path too."""
    built = CalibrationSpec(
        seq_type='exome',
        cache='calibration/exome_values.json',
        metrics=(MetricSpec(key='DEPTH', direction='min', unit='x', fail=10, warn=20, reviewed=True, rationale='r'),),
    )
    assert spec_mod.loads(spec_mod.dumps(built)) == built


def test_save_and_load_path_round_trip(tmp_path):
    path = tmp_path / 'spec.toml'
    original = spec_mod.loads(VALID)
    spec_mod.save(original, path)
    assert spec_mod.load(path) == original
