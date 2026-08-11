"""Unit tests for the calibration TOML helpers."""

import pytest

from align_genotype.qc_calibration import tomlio


@pytest.mark.parametrize(
    ('value', 'expected'),
    [
        (True, 'true'),
        (False, 'false'),
        (15, '15'),
        (3.5, '3.5'),
        (0.75, '0.75'),
        ('genome', '"genome"'),
    ],
)
def test_fmt_value(value, expected):
    assert tomlio.fmt_value(value) == expected


def test_fmt_value_bool_wins_over_int():
    """bool is a subclass of int; it must not render as 1/0."""
    assert tomlio.fmt_value(True) == 'true'


def test_fmt_value_escapes_quotes_and_backslashes():
    assert tomlio.fmt_value('a "quoted" c:\\path') == '"a \\"quoted\\" c:\\\\path"'


def test_fmt_kv():
    assert tomlio.fmt_kv('k', 3.5) == 'k = 3.5'


def test_fmt_kv_quoted_key():
    assert tomlio.fmt_kv('MEDIAN_COVERAGE', 15, quote_key=True) == '"MEDIAN_COVERAGE" = 15'


@pytest.mark.parametrize('key', ['MEDIAN_COVERAGE', 'reads_mapped_percent', 'dataset-a', 'PCT_20X'])
def test_require_bare_key_accepts_identifier_like_keys(key):
    tomlio.require_bare_key(key, 'metric')  # does not raise


@pytest.mark.parametrize('key', ['odd.key', 'has space', 'quote"d', '', 'sl/ash'])
def test_require_bare_key_rejects_anything_needing_quoting(key):
    with pytest.raises(ValueError, match='not a bare TOML key'):
        tomlio.require_bare_key(key, 'metric')


def test_loads_round_trips_written_scalars():
    text = '\n'.join(
        [
            '[table]',
            tomlio.fmt_kv('name', 'dataset-a'),
            tomlio.fmt_kv('k', 3.5),
            tomlio.fmt_kv('gated', True),
        ]
    )
    assert tomlio.loads(text) == {'table': {'name': 'dataset-a', 'k': 3.5, 'gated': True}}
