"""Unit tests for the SG QC report builders and template rendering.

These run fully offline — no Metamist/DB access — by feeding the builders the
same shape of data the DB query returns.
"""

import pytest

from align_genotype.scripts.sg_qc_report import (
    SGInfo,
    SGReport,
    _extract_reads,
    _fmt_num,
    _section_label,
    _value_display,
    build_sections,
    collect_qc_flags,
    metric_histogram,
    render_report,
    severity_histogram,
    source_histogram,
    summarise_flags,
)

# Mirrors the ``meta`` blob returned by Metamist for a flagged SG.
SAMPLE_SG = {
    'id': 'CPG276402',
    'meta': {
        'cram_qc_flags': [
            {
                'flag': 'reads_mapped_percent',
                'value': 79,
                'comparison': '<',
                'threshold': 80,
                'section': 'samtools',
                'date': '2026-07-01T15:27:41.484818',
                'ar_guid': '12345-abcde',
                'resolved': True,
                'resolution_date': '2026-07-03T02:17:23.609252',
            },
            {
                'flag': 'FREEMIX',
                'value': 0.12,
                'comparison': '>',
                'threshold': 0.04,
                'section': 'verifybamid',
                'date': '2026-07-01T15:27:41.484818',
                'ar_guid': '12345-abcde',
                'resolved': False,
                'resolution_date': None,
            },
            {
                'flag': 'MEDIAN_COVERAGE',
                'value': 22,
                'comparison': '<',
                'threshold': 30,
                'section': 'picard_4',
                'date': '2026-07-03T02:17:14.395432',
                'ar_guid': 'fccc1e42-791c-4e61-a221-d5fca0c7b072',
                'resolved': False,
                'resolution_date': None,
            },
        ],
        'gvcf_qc_flags': [
            {
                'flag': 'gvcf_issue',
                'value': 99,
                'comparison': '<',
                'threshold': 100,
                'section': 'vcfcheck',
                'date': '2026-07-02T15:27:41.484818',
                'ar_guid': '12345-abcde',
                'resolved': False,
                'resolution_date': None,
            },
        ],
    },
}


def _sg_info(sg_id: str = 'CPG276402') -> SGInfo:
    return SGInfo(
        sg_id=sg_id,
        sg_type='genome',
        sg_technology='short-read',
        sg_platform='illumina',
        crams=[],
        fastq_pairs=[('R1_001.fastq.gz', 'R2_001.fastq.gz'), ('R1_002.fastq.gz', 'R2_002.fastq.gz')],
        other_reads=[],
        sample_external_id='HG003_NA24149',
        sample_type='blood',
        participant_external_id='HG003_NA24149',
        family_external_id='GIAB_ASHKENAZI',
    )


def _report(sg_id: str = 'CPG276402') -> SGReport:
    collected = collect_qc_flags([{**SAMPLE_SG, 'id': sg_id}])[0]
    return SGReport(
        sg_info=_sg_info(sg_id),
        cram_flags=collected['cram_qc_flags'],
        gvcf_flags=collected['gvcf_qc_flags'],
    )


# --- flag collection -------------------------------------------------------
def test_collect_qc_flags_parses_meta():
    collected = collect_qc_flags([SAMPLE_SG])
    assert len(collected) == 1
    assert len(collected[0]['cram_qc_flags']) == 3
    assert len(collected[0]['gvcf_qc_flags']) == 1
    assert collected[0]['cram_qc_flags'][1].flag == 'FREEMIX'


def test_collect_qc_flags_handles_missing_meta():
    collected = collect_qc_flags([{'id': 'CPG000000', 'meta': None}])
    assert collected[0]['cram_qc_flags'] == []
    assert collected[0]['gvcf_qc_flags'] == []


# --- section split (active vs resolved) ------------------------------------
def test_build_sections_splits_active_and_resolved():
    unresolved, resolved = build_sections([_report()])
    assert len(unresolved) == 1
    assert len(resolved) == 1
    # The one SG appears in both sections, with flags partitioned.
    assert unresolved[0]['info'].sg_id == 'CPG276402'
    active_metrics = {f['flag'] for f in unresolved[0]['flags']}
    resolved_metrics = {f['flag'] for f in resolved[0]['flags']}
    assert active_metrics == {'FREEMIX', 'MEDIAN_COVERAGE', 'gvcf_issue'}
    assert resolved_metrics == {'reads_mapped_percent'}


def test_unresolved_row_summary():
    unresolved, _ = build_sections([_report()])
    row = unresolved[0]
    assert row['count_summary'] == '3 active flags'
    assert row['source_summary'] == '2 CRAM · 1 GVCF'


def test_resolved_row_summary_singular():
    _, resolved = build_sections([_report()])
    row = resolved[0]
    assert row['count_summary'] == '1 resolved flag'  # singular, no trailing 's'


def test_clean_sg_absent_from_both_sections():
    clean = SGReport(sg_info=_sg_info('CPG_CLEAN'), cram_flags=[], gvcf_flags=[])
    unresolved, resolved = build_sections([clean])
    assert unresolved == []
    assert resolved == []


def test_sections_sorted_by_family_then_participant():
    a = _report('CPG_A')
    b = _report('CPG_B')
    object.__setattr__(b.sg_info, 'family_external_id', 'AAA_FIRST')
    unresolved, _ = build_sections([a, b])
    assert [r['info'].sg_id for r in unresolved] == ['CPG_B', 'CPG_A']


# --- flag formatting -------------------------------------------------------
def test_resolved_flag_uses_resolution_date_truncated():
    _, resolved = build_sections([_report()])
    f = resolved[0]['flags'][0]
    assert f['date_full'] == '2026-07-03T02:17:23.609252'
    assert f['date_short'] == '2026-07-03'


def test_active_flag_uses_detection_date():
    unresolved, _ = build_sections([_report()])
    freemix = next(f for f in unresolved[0]['flags'] if f['flag'] == 'FREEMIX')
    assert freemix['date_short'] == '2026-07-01'


def test_metric_and_section_labels():
    unresolved, _ = build_sections([_report()])
    median = next(f for f in unresolved[0]['flags'] if f['flag'] == 'MEDIAN_COVERAGE')
    assert median['metric_label'] == 'Median coverage'
    assert median['section_label'] == 'Picard'  # 'picard_4' -> 'Picard'


def test_value_display_direction():
    # '<' means it fell below a minimum, '>' means it exceeded a maximum.
    assert _value_display(22, '<', 30, '×') == '22× (below minimum 30×)'
    assert _value_display(0.12, '>', 0.04, '') == '0.12 (above maximum 0.04)'


def test_fmt_num_rounding():
    assert _fmt_num(30) == '30'  # int stays int
    assert _fmt_num(30.0) == '30'  # integer-valued float stays int
    assert _fmt_num(64.579124) == '64.58'  # >= 1 -> 2 decimal places
    assert _fmt_num(64.5) == '64.5'  # trailing zeros stripped
    assert _fmt_num(0.0616722) == '0.062'  # < 1 -> 2 significant figures
    assert _fmt_num(0.12) == '0.12'
    assert _fmt_num(0.04) == '0.04'


def test_value_display_rounds_long_floats():
    assert _value_display(64.579124, '<', 80, '%') == '64.58% (below minimum 80%)'
    assert _value_display(0.0616722, '>', 0.04, '') == '0.062 (above maximum 0.04)'


def test_section_label_normalises_multiqc_suffix():
    assert _section_label('picard_4') == 'Picard'
    assert _section_label('verifybamid') == 'VerifyBamID'
    assert _section_label('unknown_tool_2') == 'Unknown Tool'


# --- read grouping ---------------------------------------------------------
def test_extract_reads_pairs_fastqs():
    assays = [
        {
            'meta': {
                'reads_type': 'fastq',
                'reads': [
                    {'basename': 'S_L001_R1_001.fastq.gz'},
                    {'basename': 'S_L001_R2_001.fastq.gz'},
                ],
            }
        },
        {
            'meta': {
                'reads_type': 'fastq',
                'reads': [
                    {'location': 'gs://b/S_L002_R1.fastq.gz'},
                    {'location': 'gs://b/S_L002_R2.fastq.gz'},
                ],
            }
        },
    ]
    crams, pairs, other = _extract_reads(assays)
    assert crams == []
    assert other == []
    assert pairs == [
        ('S_L001_R1_001.fastq.gz', 'S_L001_R2_001.fastq.gz'),
        ('S_L002_R1.fastq.gz', 'S_L002_R2.fastq.gz'),
    ]


def test_extract_reads_handles_cram():
    assays = [{'meta': {'reads_type': 'cram', 'reads': [{'basename': 'sample.cram'}]}}]
    crams, pairs, other = _extract_reads(assays)
    assert crams == ['sample.cram']
    assert pairs == []


def test_extract_reads_classifies_unknown_type_by_extension():
    assays = [{'meta': {'reads': [{'basename': 'x.bam'}, {'basename': 'y.fq.gz'}, {'basename': 'notes.txt'}]}}]
    crams, pairs, other = _extract_reads(assays)
    assert crams == ['x.bam']
    assert pairs == [('y.fq.gz', '')]
    assert other == ['notes.txt']


# --- summary counts --------------------------------------------------------
def test_summarise_flags_counts():
    sg_data = collect_qc_flags([SAMPLE_SG])
    summary = summarise_flags(sg_data)
    assert summary['total_sgs'] == 1
    assert summary['active_flags'] == 3
    assert summary['active_cram'] == 2
    assert summary['active_gvcf'] == 1
    assert summary['active_fail'] == 3  # SAMPLE_SG flags have no severity -> default 'fail'
    assert summary['active_warn'] == 0
    assert summary['sgs_affected'] == 1
    assert summary['resolved_flags'] == 1


def test_summarise_flags_all_resolved_is_not_affected():
    sg = {
        'id': 'CPG_R',
        'meta': {'cram_qc_flags': [{**SAMPLE_SG['meta']['cram_qc_flags'][0]}], 'gvcf_qc_flags': []},
    }
    summary = summarise_flags(collect_qc_flags([sg]))
    assert summary['active_flags'] == 0
    assert summary['sgs_affected'] == 0
    assert summary['resolved_flags'] == 1


# --- end-to-end render -----------------------------------------------------
def test_render_report_smoke():
    summary = summarise_flags(collect_qc_flags([SAMPLE_SG]))
    html = render_report('dataset', [_report()], summary=summary)
    assert 'Unresolved flags' in html
    assert 'Resolved — past incidents' in html
    assert 'QC Flag Summary' in html  # new title
    assert 'id-family">GIAB_ASHKENAZI' in html  # family id, no "Family" prefix, distinct styling
    assert 'id-participant">HG003_NA24149' in html  # participant styled separately
    assert 'CPG276402' in html  # SG id still present (muted)
    assert 'Median coverage' in html  # human metric label
    assert 'below minimum 30×' in html  # value display
    assert 'FASTQ pair' in html  # grouped reads


# --- filter histograms -----------------------------------------------------
def test_metric_histogram_counts_sgs_per_metric():
    # Two SGs each flagged for reads_duplicated_percent; one also for FREEMIX.
    dup = {**SAMPLE_SG['meta']['cram_qc_flags'][1], 'flag': 'reads_duplicated_percent'}
    free = {**SAMPLE_SG['meta']['cram_qc_flags'][1], 'flag': 'FREEMIX'}
    sg_a = collect_qc_flags([{'id': 'A', 'meta': {'cram_qc_flags': [dup], 'gvcf_qc_flags': []}}])[0]
    sg_b = collect_qc_flags([{'id': 'B', 'meta': {'cram_qc_flags': [dup, free], 'gvcf_qc_flags': []}}])[0]
    reports = [
        SGReport(_sg_info('A'), sg_a['cram_qc_flags'], []),
        SGReport(_sg_info('B'), sg_b['cram_qc_flags'], []),
    ]
    unresolved, _ = build_sections(reports)
    hist = metric_histogram(unresolved)
    # Sorted by count desc: duplicated reads (2 SGs) before contamination (1 SG).
    assert [(h['key'], h['count']) for h in hist] == [
        ('reads_duplicated_percent', 2),
        ('FREEMIX', 1),
    ]
    assert hist[0]['label'] == 'Duplicated reads'


def test_source_histogram():
    unresolved, _ = build_sections([_report()])
    hist = source_histogram(unresolved)
    assert {h['key']: h['count'] for h in hist} == {'CRAM': 1, 'GVCF': 1}


def test_render_report_includes_filter_bar_and_row_data():
    summary = summarise_flags(collect_qc_flags([SAMPLE_SG]))
    html = render_report('dataset', [_report()], summary=summary)
    assert 'filter-bar' in html
    assert 'metric-chip' in html
    assert 'data-metrics="' in html
    assert 'data-search="' in html
    # metric key present as a chip value, human label as its text
    assert 'data-metric="MEDIAN_COVERAGE"' in html
    assert 'Median coverage (1)' in html


def test_render_report_shows_flags_at_top_level():
    # The 3 active flags render inline (2 shown) with a "+1 more" affordance,
    # without needing to expand the detail row.
    summary = summarise_flags(collect_qc_flags([SAMPLE_SG]))
    html = render_report('dataset', [_report()], summary=summary)
    assert 'flag-glance' in html
    assert '+1 more' in html  # 3 active flags -> 2 shown inline + 1 more


def test_render_report_all_clear_banner():
    clean = SGReport(sg_info=_sg_info('CPG_CLEAN'), cram_flags=[], gvcf_flags=[])
    summary = {
        'total_sgs': 5,
        'active_flags': 0,
        'active_cram': 0,
        'active_gvcf': 0,
        'active_fail': 0,
        'active_warn': 0,
        'sgs_affected': 0,
        'resolved_flags': 0,
    }
    html = render_report('dataset', [clean], summary=summary)
    assert 'All clear' in html


# --- severity (warn vs fail) -----------------------------------------------
def _sg_with_severities() -> dict:
    """A flagged SG carrying one warn and one fail CRAM flag."""
    return {
        'id': 'CPG_SEV',
        'meta': {
            'cram_qc_flags': [
                {
                    'flag': 'FOLD_80_BASE_PENALTY',
                    'value': 2.5,
                    'comparison': '>',
                    'threshold': 2.0,
                    'section': 'picard',
                    'date': '2026-07-01T00:00:00',
                    'ar_guid': 'x',
                    'severity': 'warn',
                    'resolved': False,
                    'resolution_date': None,
                },
                {
                    'flag': 'MEAN_TARGET_COVERAGE',
                    'value': 40,
                    'comparison': '<',
                    'threshold': 50,
                    'section': 'picard',
                    'date': '2026-07-01T00:00:00',
                    'ar_guid': 'x',
                    'severity': 'fail',
                    'resolved': False,
                    'resolution_date': None,
                },
            ],
            'gvcf_qc_flags': [],
        },
    }


def _sev_report() -> SGReport:
    collected = collect_qc_flags([_sg_with_severities()])[0]
    return SGReport(_sg_info('CPG_SEV'), collected['cram_qc_flags'], [])


def test_missing_severity_defaults_to_fail():
    # A stored flag without a 'severity' key (pre-tier data) loads as 'fail'.
    collected = collect_qc_flags([SAMPLE_SG])[0]
    assert all(f.severity == 'fail' for f in collected['cram_qc_flags'])


def test_summarise_flags_splits_warn_and_fail():
    summary = summarise_flags(collect_qc_flags([_sg_with_severities()]))
    assert summary['active_fail'] == 1
    assert summary['active_warn'] == 1


def test_row_orders_fail_before_warn():
    unresolved, _ = build_sections([_sev_report()])
    row = unresolved[0]
    assert row['n_fail'] == 1
    assert row['n_warn'] == 1
    assert row['row_severity'] == 'fail'
    # The fail flag leads the at-a-glance list.
    assert row['flags'][0]['severity'] == 'fail'
    assert row['flags'][0]['flag'] == 'MEAN_TARGET_COVERAGE'


def test_render_report_shows_severity_badges_and_chips():
    summary = summarise_flags(collect_qc_flags([_sg_with_severities()]))
    html = render_report('dataset', [_sev_report()], summary=summary)
    assert 'badge badge-warn' in html
    assert 'badge badge-fail' in html
    assert 'Failing flags' in html  # summary card
    assert 'Warning flags' in html
    assert 'data-severities="' in html
    assert 'severity-chip' in html  # both tiers present -> severity filter shown


def test_severity_histogram_counts():
    unresolved, _ = build_sections([_sev_report()])
    hist = severity_histogram(unresolved)
    assert {h['key']: h['count'] for h in hist} == {'fail': 1, 'warn': 1}


if __name__ == '__main__':
    pytest.main([__file__, '-v'])
