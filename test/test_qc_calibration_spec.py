"""Unit tests for the calibration spec: load, validate, dump, round-trip."""

import pytest

from align_genotype.qc_calibration import spec as spec_mod
from align_genotype.qc_calibration.spec import MetricSpec, RelativeSpec, SpecError

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
            '[metrics.M]\ndirection = "max"\nwarn = 1\n[metrics.M.relative]\nk = 3.5\n',
            'relative tier requires an absolute fail',
        ),
        ('[metrics.M]\ndirection = "max"\ngated = false\nfail = 1\n', 'non-gated metric must not define'),
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


def test_no_metrics_rejected():
    with pytest.raises(SpecError, match='defines no metrics'):
        spec_mod.loads('seq_type = "genome"\ncache = "c.json"\n')


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


def test_save_and_load_path_round_trip(tmp_path):
    path = tmp_path / 'spec.toml'
    original = spec_mod.loads(VALID)
    spec_mod.save(original, path)
    assert spec_mod.load(path) == original
