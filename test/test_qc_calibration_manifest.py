"""Unit tests for the cohort manifest."""

import pytest

from align_genotype.qc_calibration import manifest as manifest_mod
from align_genotype.qc_calibration.manifest import Cohort, Manifest, ManifestError

VALID = """
seq_type = "genome"
generated = "2026-08-11T13:40:00"

[cohorts.dataset-a]
uri = "gs://cpg-dataset-a-main/qc/multiqc_data.json"
analysis_id = 84213
timestamp = "2026-06-02T04:11:09"

[cohorts.dataset-b]
uri = "file:///tmp/dataset-b_multiqc_data.json"
"""


def test_loads_valid_manifest():
    m = manifest_mod.loads(VALID)
    assert m.seq_type == 'genome'
    assert m.generated == '2026-08-11T13:40:00'
    assert m.labels == ('dataset-a', 'dataset-b')


def test_cohort_fields_and_defaults():
    m = manifest_mod.loads(VALID)
    assert m.cohort('dataset-a') == Cohort(
        label='dataset-a',
        uri='gs://cpg-dataset-a-main/qc/multiqc_data.json',
        analysis_id=84213,
        timestamp='2026-06-02T04:11:09',
    )
    assert m.cohort('dataset-b') == Cohort(label='dataset-b', uri='file:///tmp/dataset-b_multiqc_data.json')


def test_unknown_cohort_raises():
    with pytest.raises(KeyError, match='dataset-z'):
        manifest_mod.loads(VALID).cohort('dataset-z')


def test_unknown_cohort_error_lists_known_labels():
    """The operator hand-edits these files; a typo should say what was available."""
    with pytest.raises(KeyError, match='dataset-a'):
        manifest_mod.loads(VALID).cohort('dataset-z')


@pytest.mark.parametrize(
    ('text', 'match'),
    [
        ('generated = "x"\n[cohorts.a]\nuri = "u"\n', 'missing required key: seq_type'),
        ('seq_type = "genome"\ngenerated = "x"\n', 'lists no cohorts'),
        (
            'seq_type = "genome"\ngenerated = "x"\n[cohorts.a]\nanalysis_id = 1\n',
            "cohort 'a': missing required key: uri",
        ),
        ('seq_type = "genome"\ngenerated = "x"\n[cohorts."a.b"]\nuri = "u"\n', 'not a bare TOML key'),
        ('seq_type = "genome"\ngenerated = "x"\n[cohorts.a]\nuri = "u"\nnope = 1\n', r"unknown key\(s\) \['nope'\]"),
    ],
)
def test_validation_errors(text, match):
    with pytest.raises(ManifestError, match=match):
        manifest_mod.loads(text)


def test_non_numeric_analysis_id_rejected():
    """int() on a bad value would raise ValueError, which the CLI's except ManifestError misses."""
    text = 'seq_type = "genome"\ngenerated = "x"\n[cohorts.a]\nuri = "u"\nanalysis_id = "many"\n'
    with pytest.raises(ManifestError, match='analysis_id'):
        manifest_mod.loads(text)


def test_empty_uri_rejected():
    text = 'seq_type = "genome"\ngenerated = "x"\n[cohorts.a]\nuri = ""\n'
    with pytest.raises(ManifestError, match='uri'):
        manifest_mod.loads(text)


def test_dump_round_trips():
    original = manifest_mod.loads(VALID)
    assert manifest_mod.loads(manifest_mod.dumps(original)) == original


def test_dump_omits_absent_optional_fields():
    text = manifest_mod.dumps(manifest_mod.loads(VALID))
    dataset_b_block = text[text.index('[cohorts.dataset-b]') :]
    assert 'analysis_id' not in dataset_b_block
    assert 'timestamp' not in dataset_b_block


def test_dumps_guards_bare_keys_on_a_directly_constructed_manifest():
    """dumps writes labels unquoted, so a label needing quotes can't round-trip."""
    bad = Manifest(seq_type='genome', generated='x', cohorts=(Cohort(label='odd.name', uri='u'),))
    with pytest.raises(ValueError, match='not a bare TOML key'):
        manifest_mod.dumps(bad)


def test_save_and_load_path_round_trip(tmp_path):
    path = tmp_path / 'manifest.toml'
    original = manifest_mod.loads(VALID)
    manifest_mod.save(original, path)
    assert manifest_mod.load(path) == original
