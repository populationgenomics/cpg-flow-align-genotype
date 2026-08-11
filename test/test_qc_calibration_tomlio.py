"""Unit tests for the calibration TOML helpers."""

import numpy as np
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


@pytest.mark.parametrize(
    ('value', 'expected'),
    [
        (np.float64(0.75), '0.75'),  # subclasses float, but reprs as np.float64(0.75)
        (np.int64(15), '15'),  # subclasses nothing we check - would become the string "15"
        (np.bool_(True), 'true'),  # likewise - would become the string "True"
        (np.float32(0.5), '0.5'),
    ],
)
def test_fmt_value_unwraps_numpy_scalars(value, expected):
    """Every threshold this tool writes comes out of numpy; the sink must not mangle them."""
    assert tomlio.fmt_value(value) == expected


def test_numpy_scalars_round_trip_with_the_right_type():
    text = '\n'.join(
        [
            '[table]',
            tomlio.fmt_kv('depth', np.float64(0.75)),
            tomlio.fmt_kv('count', np.int64(15)),
            tomlio.fmt_kv('gated', np.bool_(True)),
        ]
    )
    parsed = tomlio.loads(text)['table']
    assert parsed == {'depth': 0.75, 'count': 15, 'gated': True}
    assert isinstance(parsed['count'], int) and not isinstance(parsed['count'], bool)
    assert isinstance(parsed['gated'], bool)


@pytest.mark.parametrize('value', [None, [1, 2], {'a': 1}, object()])
def test_fmt_value_rejects_non_scalars(value):
    """Better a loud TypeError than a plausible-looking quoted string in a config file."""
    with pytest.raises(TypeError, match='as a TOML scalar'):
        tomlio.fmt_value(value)


def test_require_bare_key_rejects_a_trailing_newline():
    r"""`$` matches before a final \n, so this needs fullmatch, not match."""
    with pytest.raises(ValueError, match='not a bare TOML key'):
        tomlio.require_bare_key('abc\n', 'metric')


def test_require_bare_key_names_the_kind():
    """Task 3 passes 'cohort label' to distinguish error sources; that half must work."""
    with pytest.raises(ValueError, match='cohort label'):
        tomlio.require_bare_key('odd.name', 'cohort label')


def test_load_path_reads_a_local_file(tmp_path):
    path = tmp_path / 'sample.toml'
    path.write_text('[table]\nname = "dataset-a"\nk = 3.5\n')
    assert tomlio.load_path(path) == {'table': {'name': 'dataset-a', 'k': 3.5}}


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
