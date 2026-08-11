"""Unit tests for Metamist cohort discovery, using a fake query function."""

import pytest

from align_genotype.qc_calibration import discovery as discovery_mod
from align_genotype.qc_calibration import manifest as manifest_mod

PROJECTS = {
    'myProjects': [
        {'name': 'dataset-a', 'dataset': 'dataset-a', 'meta': {'is_seqr': True}},
        {'name': 'dataset-b', 'dataset': 'dataset-b', 'meta': {'is_seqr': True}},
        {'name': 'dataset-a-test', 'dataset': 'dataset-a-test', 'meta': {'is_seqr': True}},
        {'name': 'dataset-c-training', 'dataset': 'dataset-c-training', 'meta': {'is_seqr': True}},
        {'name': 'dataset-d', 'dataset': 'dataset-d', 'meta': {'is_seqr': False}},
    ],
}

ANALYSES = {
    'dataset-a': [
        {
            'id': 1,
            'meta': {'sequencing_type': 'genome'},
            'output': 'gs://a/old.json',
            'timestampCompleted': '2026-01-01T00:00:00',
        },
        {
            'id': 2,
            'meta': {'sequencing_type': 'genome'},
            'output': 'gs://a/new.json',
            'timestampCompleted': '2026-06-01T00:00:00',
        },
        {
            'id': 3,
            'meta': {'sequencing_type': 'exome'},
            'output': 'gs://a/exome.json',
            'timestampCompleted': '2026-07-01T00:00:00',
        },
    ],
    'dataset-b': [
        {
            'id': 4,
            'meta': {'sequencing_type': 'genome'},
            'output': 'gs://b/g.json',
            'timestampCompleted': '2026-03-01T00:00:00',
        },
    ],
}


def _fake_query(query_text, variables=None) -> dict:
    if 'myProjects' in query_text:
        return PROJECTS
    return {'project': {'analyses': ANALYSES.get(variables['datasetName'], [])}}


def test_selects_eligible_datasets_only():
    manifest = discovery_mod.build_manifest('genome', query_fn=_fake_query, generated='2026-08-11')
    assert manifest.labels == ('dataset-a', 'dataset-b')


@pytest.mark.parametrize(
    ('project', 'eligible'),
    [
        ({'name': 'dataset-a', 'meta': {'is_seqr': True}}, True),
        ({'name': 'dataset-a-test', 'meta': {'is_seqr': True}}, False),
        ({'name': 'dataset-c-training', 'meta': {'is_seqr': True}}, False),
        ({'name': 'dataset-seqr-x', 'meta': {'is_seqr': True}}, False),
        ({'name': 'dataset-d', 'meta': {'is_seqr': False}}, False),
        ({'name': 'dataset-e', 'meta': {}}, False),
        ({'name': 'dataset-f'}, False),
        # Malformed rows from a live external service must not crash is_eligible.
        ({'name': 'dataset-a', 'meta': 'oops'}, False),
        ({'name': 'dataset-a', 'meta': ['is_seqr']}, False),
        ({'name': None, 'meta': {'is_seqr': True}}, True),
    ],
)
def test_dataset_eligibility(project, eligible):
    assert discovery_mod.is_eligible(project) is eligible


def test_excluded_dataset_name_is_logged_at_info(caplog):
    caplog.set_level('INFO')
    assert discovery_mod.is_eligible({'name': 'dataset-a-test', 'meta': {'is_seqr': True}}) is False
    assert 'dataset-a-test' in caplog.text


def test_picks_the_latest_analysis_for_the_requested_seq_type():
    manifest = discovery_mod.build_manifest('genome', query_fn=_fake_query, generated='2026-08-11')
    cohort = manifest.cohort('dataset-a')
    assert cohort.uri == 'gs://a/new.json'
    assert cohort.analysis_id == 2
    assert cohort.timestamp == '2026-06-01T00:00:00'


def test_latest_is_by_timestamp_not_by_id_or_order():
    """A higher id with an older timestamp must lose."""
    analyses = [
        {
            'id': 99,
            'meta': {'sequencing_type': 'genome'},
            'output': 'gs://a/old.json',
            'timestampCompleted': '2020-01-01T00:00:00',
        },
        {
            'id': 1,
            'meta': {'sequencing_type': 'genome'},
            'output': 'gs://a/new.json',
            'timestampCompleted': '2026-01-01T00:00:00',
        },
    ]
    assert discovery_mod.latest_analysis(analyses, 'genome', 'dataset-a')['id'] == 1


def test_other_seq_types_are_ignored():
    manifest = discovery_mod.build_manifest('exome', query_fn=_fake_query, generated='2026-08-11')
    assert manifest.labels == ('dataset-a',)
    assert manifest.cohort('dataset-a').uri == 'gs://a/exome.json'


def test_dataset_with_no_matching_analysis_is_omitted():
    manifest = discovery_mod.build_manifest('exome', query_fn=_fake_query, generated='2026-08-11')
    assert 'dataset-b' not in manifest.labels


def test_dataset_with_no_matching_analysis_logs_at_info(caplog):
    caplog.set_level('INFO')
    discovery_mod.build_manifest('exome', query_fn=_fake_query, generated='2026-08-11')
    assert 'dataset-b' in caplog.text


def test_project_with_no_name_or_dataset_is_skipped_with_a_warning(caplog):
    """is_eligible tolerates a None name, but build_manifest still needs a label."""

    def query_fn(query_text, variables=None) -> dict:  # noqa: ARG001 - fake must match QueryFn signature
        if 'myProjects' in query_text:
            return {'myProjects': [{'name': None, 'meta': {'is_seqr': True}}]}
        return {'project': {'analyses': []}}

    with pytest.raises(discovery_mod.DiscoveryError, match='no cohorts'):
        discovery_mod.build_manifest('genome', query_fn=query_fn, generated='x')
    assert 'no dataset/name' in caplog.text


def test_null_project_response_is_tolerated():
    """A present-but-null `project` (e.g. no read access) must not crash discovery."""

    def query_fn(query_text, variables=None) -> dict:  # noqa: ARG001 - fake must match QueryFn signature
        if 'myProjects' in query_text:
            return {'myProjects': [{'name': 'dataset-a', 'dataset': 'dataset-a', 'meta': {'is_seqr': True}}]}
        return {'project': None}

    with pytest.raises(discovery_mod.DiscoveryError, match='no cohorts'):
        discovery_mod.build_manifest('genome', query_fn=query_fn, generated='x')


def test_analysis_with_null_timestamp_is_skipped_but_valid_one_wins(caplog):
    def query_fn(query_text, variables=None) -> dict:  # noqa: ARG001 - fake must match QueryFn signature
        if 'myProjects' in query_text:
            return {'myProjects': [{'name': 'dataset-a', 'dataset': 'dataset-a', 'meta': {'is_seqr': True}}]}
        return {
            'project': {
                'analyses': [
                    {
                        'id': 1,
                        'meta': {'sequencing_type': 'genome'},
                        'output': 'gs://a/no-ts.json',
                        'timestampCompleted': None,
                    },
                    {
                        'id': 2,
                        'meta': {'sequencing_type': 'genome'},
                        'output': 'gs://a/g.json',
                        'timestampCompleted': '2026-01-01T00:00:00',
                    },
                ],
            },
        }

    manifest = discovery_mod.build_manifest('genome', query_fn=query_fn, generated='x')
    cohort = manifest.cohort('dataset-a')
    assert cohort.analysis_id == 2
    assert cohort.uri == 'gs://a/g.json'
    assert 'dataset-a' in caplog.text
    assert '1' in caplog.text
    assert '1 analyses across 1 datasets had unusable timestamps' in caplog.text


def test_all_null_timestamps_omits_dataset_and_raises_when_none_left():
    def query_fn(query_text, variables=None) -> dict:  # noqa: ARG001 - fake must match QueryFn signature
        if 'myProjects' in query_text:
            return {'myProjects': [{'name': 'dataset-a', 'dataset': 'dataset-a', 'meta': {'is_seqr': True}}]}
        return {
            'project': {
                'analyses': [
                    {
                        'id': 1,
                        'meta': {'sequencing_type': 'genome'},
                        'output': 'gs://a/no-ts.json',
                        'timestampCompleted': None,
                    },
                    {
                        'id': 2,
                        'meta': {'sequencing_type': 'genome'},
                        'output': 'gs://a/no-ts-2.json',
                        'timestampCompleted': 12345,
                    },
                ],
            },
        }

    with pytest.raises(discovery_mod.DiscoveryError, match='no cohorts'):
        discovery_mod.build_manifest('genome', query_fn=query_fn, generated='x')


def test_analysis_without_an_output_path_is_ignored():
    def query_fn(query_text, variables=None) -> dict:  # noqa: ARG001 - fake must match QueryFn signature
        if 'myProjects' in query_text:
            return {'myProjects': [{'name': 'dataset-a', 'dataset': 'dataset-a', 'meta': {'is_seqr': True}}]}
        return {
            'project': {
                'analyses': [
                    {
                        'id': 1,
                        'meta': {'sequencing_type': 'genome'},
                        'output': None,
                        'timestampCompleted': '2026-09-01T00:00:00',
                    },
                    {
                        'id': 2,
                        'meta': {'sequencing_type': 'genome'},
                        'output': 'gs://a/g.json',
                        'timestampCompleted': '2026-01-01T00:00:00',
                    },
                ],
            },
        }

    manifest = discovery_mod.build_manifest('genome', query_fn=query_fn, generated='x')
    assert manifest.cohort('dataset-a').analysis_id == 2


def test_analysis_with_no_meta_is_ignored():
    def query_fn(query_text, variables=None) -> dict:  # noqa: ARG001 - fake must match QueryFn signature
        if 'myProjects' in query_text:
            return {'myProjects': [{'name': 'dataset-a', 'dataset': 'dataset-a', 'meta': {'is_seqr': True}}]}
        return {'project': {'analyses': [{'id': 1, 'output': 'gs://a/g.json', 'timestampCompleted': 'x'}]}}

    with pytest.raises(discovery_mod.DiscoveryError, match='no cohorts'):
        discovery_mod.build_manifest('genome', query_fn=query_fn, generated='x')


def test_dataset_label_that_is_not_a_bare_toml_key_is_skipped(caplog):
    def query_fn(query_text, variables=None) -> dict:  # noqa: ARG001 - fake must match QueryFn signature
        if 'myProjects' in query_text:
            return {'myProjects': [{'name': 'odd.name', 'dataset': 'odd.name', 'meta': {'is_seqr': True}}]}
        return {
            'project': {
                'analyses': [
                    {
                        'id': 1,
                        'meta': {'sequencing_type': 'genome'},
                        'output': 'gs://x/g.json',
                        'timestampCompleted': '2026-01-01T00:00:00',
                    },
                ]
            }
        }

    with pytest.raises(discovery_mod.DiscoveryError, match='no cohorts'):
        discovery_mod.build_manifest('genome', query_fn=query_fn, generated='x')
    assert 'odd.name' in caplog.text


def test_no_eligible_datasets_raises():
    def query_fn(query_text, variables=None) -> dict:  # noqa: ARG001 - fake must match QueryFn signature
        return {'myProjects': []}

    with pytest.raises(discovery_mod.DiscoveryError, match='no cohorts'):
        discovery_mod.build_manifest('genome', query_fn=query_fn, generated='x')


def test_error_names_the_eligibility_filter_so_it_can_be_diagnosed():
    def query_fn(query_text, variables=None) -> dict:  # noqa: ARG001 - fake must match QueryFn signature
        return {'myProjects': []}

    with pytest.raises(discovery_mod.DiscoveryError, match='is_seqr'):
        discovery_mod.build_manifest('genome', query_fn=query_fn, generated='x')


def test_error_also_names_the_name_exclusion_filter():
    def query_fn(query_text, variables=None) -> dict:  # noqa: ARG001 - fake must match QueryFn signature
        return {'myProjects': []}

    with pytest.raises(discovery_mod.DiscoveryError, match='training'):
        discovery_mod.build_manifest('genome', query_fn=query_fn, generated='x')


def test_manifest_round_trips_through_toml():
    original = discovery_mod.build_manifest('genome', query_fn=_fake_query, generated='2026-08-11')
    assert manifest_mod.loads(manifest_mod.dumps(original)) == original


def test_analysis_ids_are_ints_not_strings():
    """Manifest validation requires a genuine int; GraphQL can return either."""

    def query_fn(query_text, variables=None) -> dict:  # noqa: ARG001 - fake must match QueryFn signature
        if 'myProjects' in query_text:
            return {'myProjects': [{'name': 'dataset-a', 'dataset': 'dataset-a', 'meta': {'is_seqr': True}}]}
        return {
            'project': {
                'analyses': [
                    {
                        'id': '7',
                        'meta': {'sequencing_type': 'genome'},
                        'output': 'gs://a/g.json',
                        'timestampCompleted': 'x',
                    },
                ]
            }
        }

    manifest = discovery_mod.build_manifest('genome', query_fn=query_fn, generated='x')
    assert manifest.cohort('dataset-a').analysis_id == 7
    assert isinstance(manifest.cohort('dataset-a').analysis_id, int)
