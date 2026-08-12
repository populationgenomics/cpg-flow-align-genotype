"""Unit tests for per-dataset CramMultiQC report discovery, using a fake query function."""

from align_genotype.qc_calibration import discovery as discovery_mod


def analysis(analysis_id, output, timestamp) -> dict:
    return {'id': analysis_id, 'output': output, 'timestampCompleted': timestamp}


def fake_query(analyses):
    """A query function returning `analyses`, recording the variables it was given."""
    calls = []

    def _query(query_text, variables=None):  # noqa: ANN202, ARG001
        calls.append(variables)
        return {'project': {'analyses': analyses}}

    _query.calls = calls
    return _query


def test_returns_the_newest_completed_report():
    query = fake_query(
        [
            analysis(1, 'gs://a/old/multiqc_data.json', '2026-01-01T00:00:00'),
            analysis(2, 'gs://a/new/multiqc_data.json', '2026-06-01T00:00:00'),
        ],
    )
    found = discovery_mod.latest_cram_multiqc('ds-a', 'genome', query_fn=query)
    assert found.dataset == 'ds-a'
    assert found.analysis_id == 2
    assert found.uri == 'gs://a/new/multiqc_data.json'
    assert found.timestamp == '2026-06-01T00:00:00'


def test_selects_on_timestamp_not_analysis_id():
    """A higher id with an older timestamp must lose."""
    query = fake_query(
        [
            analysis(99, 'gs://a/old/multiqc_data.json', '2026-01-01T00:00:00'),
            analysis(2, 'gs://a/new/multiqc_data.json', '2026-06-01T00:00:00'),
        ],
    )
    assert discovery_mod.latest_cram_multiqc('ds-a', 'genome', query_fn=query).analysis_id == 2


def test_filters_server_side_on_stage_and_sequencing_type():
    """CramMultiQC and GvcfMultiQC are both analysis_type='qc'."""
    query = fake_query([analysis(1, 'gs://a/multiqc_data.json', '2026-01-01T00:00:00')])
    discovery_mod.latest_cram_multiqc('ds-a', 'genome', query_fn=query)
    assert query.calls[0] == {
        'dataset': 'ds-a',
        'analysisType': 'qc',
        'metaFilter': {'stage': {'eq': 'CramMultiQC'}, 'sequencing_type': {'eq': 'genome'}},
    }


def test_ignores_the_html_analysis_of_the_same_stage():
    """analysis_keys=['json','html'] registers two analyses per run."""
    query = fake_query(
        [
            analysis(1, 'gs://a/multiqc.html', '2026-06-02T00:00:00'),
            analysis(2, 'gs://a/multiqc_data.json', '2026-06-01T00:00:00'),
        ],
    )
    found = discovery_mod.latest_cram_multiqc('ds-a', 'genome', query_fn=query)
    assert found.uri == 'gs://a/multiqc_data.json'


def test_no_matching_analysis_returns_none():
    assert discovery_mod.latest_cram_multiqc('ds-a', 'genome', query_fn=fake_query([])) is None


def test_only_html_analyses_returns_none():
    query = fake_query([analysis(1, 'gs://a/multiqc.html', '2026-06-01T00:00:00')])
    assert discovery_mod.latest_cram_multiqc('ds-a', 'genome', query_fn=query) is None


def test_analysis_without_an_output_is_skipped():
    query = fake_query(
        [
            analysis(1, None, '2026-06-02T00:00:00'),
            analysis(2, 'gs://a/multiqc_data.json', '2026-06-01T00:00:00'),
        ],
    )
    assert discovery_mod.latest_cram_multiqc('ds-a', 'genome', query_fn=query).analysis_id == 2


def test_unrankable_timestamp_is_skipped_with_a_warning(caplog):
    """A row we cannot order must be dropped, not compared - and never silently."""
    query = fake_query(
        [
            analysis(1, 'gs://a/broken/multiqc_data.json', None),
            analysis(2, 'gs://a/good/multiqc_data.json', '2026-01-01T00:00:00'),
        ],
    )
    with caplog.at_level('WARNING'):
        found = discovery_mod.latest_cram_multiqc('ds-a', 'genome', query_fn=query)
    assert found.analysis_id == 2
    assert 'timestampCompleted' in caplog.text


def test_all_timestamps_unrankable_returns_none():
    query = fake_query([analysis(1, 'gs://a/multiqc_data.json', None)])
    assert discovery_mod.latest_cram_multiqc('ds-a', 'genome', query_fn=query) is None


def test_a_null_project_returns_none():
    """`project` can be present-but-null, e.g. with no read access."""

    def _query(query_text, variables=None):  # noqa: ANN202, ARG001
        return {'project': None}

    assert discovery_mod.latest_cram_multiqc('ds-a', 'genome', query_fn=_query) is None


def test_a_null_analyses_list_returns_none():
    def _query(query_text, variables=None):  # noqa: ANN202, ARG001
        return {'project': {'analyses': None}}

    assert discovery_mod.latest_cram_multiqc('ds-a', 'genome', query_fn=_query) is None


def test_exome_sequencing_type_is_passed_through():
    query = fake_query([analysis(1, 'gs://a/multiqc_data.json', '2026-01-01T00:00:00')])
    discovery_mod.latest_cram_multiqc('ds-a', 'exome', query_fn=query)
    assert query.calls[0]['metaFilter']['sequencing_type'] == {'eq': 'exome'}


def test_cached_latest_cram_multiqc_memoizes_per_dataset_and_seq_type(monkeypatch):
    """The DAG-assembly driver calls this repeatedly; it must not re-query per call."""
    discovery_mod.cached_latest_cram_multiqc.cache_clear()
    calls = []

    def fake(dataset, seq_type):  # noqa: ANN202
        calls.append((dataset, seq_type))
        return discovery_mod.MultiqcReport(dataset=dataset, uri='gs://x', analysis_id=1, timestamp='t')

    monkeypatch.setattr(discovery_mod, 'latest_cram_multiqc', fake)
    discovery_mod.cached_latest_cram_multiqc('ds-a', 'genome')
    discovery_mod.cached_latest_cram_multiqc('ds-a', 'genome')
    discovery_mod.cached_latest_cram_multiqc('ds-b', 'genome')
    discovery_mod.cached_latest_cram_multiqc('ds-a', 'exome')

    assert calls == [('ds-a', 'genome'), ('ds-b', 'genome'), ('ds-a', 'exome')]
    discovery_mod.cached_latest_cram_multiqc.cache_clear()


def test_cached_latest_cram_multiqc_caches_a_none_result(monkeypatch):
    """A dataset with no report yet must not be re-queried within the same process."""
    discovery_mod.cached_latest_cram_multiqc.cache_clear()
    calls = []

    def fake(dataset, seq_type):  # noqa: ANN202
        calls.append((dataset, seq_type))

    monkeypatch.setattr(discovery_mod, 'latest_cram_multiqc', fake)
    assert discovery_mod.cached_latest_cram_multiqc('ds-a', 'genome') is None
    assert discovery_mod.cached_latest_cram_multiqc('ds-a', 'genome') is None

    assert calls == [('ds-a', 'genome')]
    discovery_mod.cached_latest_cram_multiqc.cache_clear()
