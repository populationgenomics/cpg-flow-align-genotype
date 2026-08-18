"""
Queries Metamist for all QC flags across a dataset's sequencing groups
of a given type (exome or genome) and renders them into a report using
the sg_qc_overview.html.jinja template.
"""

import re
from argparse import ArgumentParser
from dataclasses import dataclass
from datetime import datetime
from pathlib import Path
from time import perf_counter

import jinja2
from loguru import logger

from cpg_utils import to_path
from cpg_utils.config import config_retrieve, dataset_for_access_level
from cpg_utils.metamist_registration import create_new
from cpg_utils.slack import send_message
from metamist.graphql import gql, query

from align_genotype.utils import QcFlag

STAGE_NAME = 'GenerateSgQcReport'
JINJA_TEMPLATE_DIR = Path(__file__).absolute().parent.parent / 'templates'

DATASET_SGS_QUERY = gql(
    """
    query datasetSgs($dataset: String!, $seqType: String!, $seqTech: String!) {
        project(name: $dataset) {
            sequencingGroups(type: {eq: $seqType}, technology: {eq: $seqTech}) {
                id
                meta
                type
            }
        }
    }
    """
)

SGS_INFO_QUERY = gql(
    """
    query sgInfo($sgIds: [String!]!) {
        sequencingGroups(id: {in_: $sgIds}) {
            id
            meta
            type
            technology
            platform
            assays {
                id
                meta
            }
            sample {
                id
                externalIds
                type
                participant {
                    id
                    externalIds
                    families {
                        id
                        externalIds
                    }
                }
            }
        }
    }
    """
)

EXISTING_ANALYSES_QUERY = gql(
    """
    query existingAnalyses($dataset: String!, $metaFilter: JSON!) {
        project(name: $dataset) {
            analyses(type: {eq: "web"}, meta: $metaFilter) {
                id
                outputs
                timestampCompleted
                meta
            }
        }
    }
    """
)

# ---------------------------------------------------------------------------
# Human-readable labels for MultiQC metric keys and module sections.
#
# The `flag`/`section` on a QcFlag are the raw MultiQC metric key and module
# key (see check_multiqc.py). MultiQC's friendly headers aren't captured in the
# JSON we parse, so we maintain a small map here. The metric set is small and
# stable (see config_template.toml :: qc_thresholds & README.md).
#
# Note that some metrics use integers for percentages, others use floats in [0,1].
# ---------------------------------------------------------------------------
METRIC_LABELS: dict[str, tuple[str, str, int]] = {
    # Metric source
    # Metric key: (human label, unit suffix, multiplier to convert to unit)
    # samtools stats metrics
    'reads_mapped_percent': ('Reads mapped', '%', 1),
    'reads_duplicated_percent': ('Duplicated reads', '%', 1),
    'reads_properly_paired_percent': ('Reads properly paired', '%', 1),
    # Picard CollectWgsMetrics (Genome)
    'PCT_PF_READS_ALIGNED': ('Passing Filter Reads aligned', '%', 100),
    'PCT_20X': ('Bases ≥20× coverage', '%', 100),
    'MEDIAN_COVERAGE': ('Median coverage', '×', 1),
    'MEAN_COVERAGE': ('Mean coverage', '×', 1),
    # Picard CollectHsMetrics (Exome) target-coverage metrics
    'ZERO_CVG_TARGETS_PCT': ('Zero-coverage targets', '%', 100),
    'PCT_TARGET_BASES_20X': ('Target bases ≥20×', '%', 100),
    'FOLD_80_BASE_PENALTY': ('Fold-80 base penalty', '×', 1),
    'MEAN_TARGET_COVERAGE': ('Mean target coverage', '×', 1),
    # VerifyBamID2 contamination metric
    'FREEMIX': ('Contamination (FreeMix)', '%', 100),
}


SECTION_LABELS: dict[str, str] = {
    'samtools': 'Samtools',
    'verifybamid': 'VerifyBamID',
    'picard': 'Picard',
}


@dataclass
class SGInfo:
    sg_id: str
    sg_type: str
    sg_technology: str
    sg_platform: str
    crams: list[str]
    fastq_pairs: list[tuple[str, str]]
    other_reads: list[str]
    sample_external_id: str
    sample_type: str
    participant_external_id: str
    family_external_id: str


@dataclass(frozen=True)
class SGReport:
    """All QC flags for a single sequencing group, plus its metadata."""

    sg_info: SGInfo
    cram_flags: list[QcFlag]
    gvcf_flags: list[QcFlag]


def _has_active(flags: list[QcFlag]) -> bool:
    return any(not f.resolved for f in flags)


# ---------------------------------------------------------------------------
# Flag label / value / date formatting
# ---------------------------------------------------------------------------
def _metric_label(metric: str) -> tuple[str, str, int]:
    """(human label, unit, multiplier) for a MultiQC metric key; falls back to the key."""
    return METRIC_LABELS.get(metric, (metric, '', 1))


def _section_label(section: str) -> str:
    """Human tool name for a MultiQC module key (e.g. 'picard_4' -> 'Picard')."""
    base = re.sub(r'_\d+$', '', section or '')  # strip MultiQC module-instance suffix
    return SECTION_LABELS.get(base.lower(), base.replace('_', ' ').title() or '—')


def _fmt_num(n: float | str) -> str:
    """Format a metric number for display.

    Integers stay integers (30 -> '30'); values >= 1 get 2 decimal places
    (64.579124 -> '64.58'); small values < 1 get 2 significant figures
    (0.0616722 -> '0.062'). Trailing zeros are stripped (64.50 -> '64.5').
    """
    if isinstance(n, bool) or not isinstance(n, (int, float)):
        return str(n)
    f = float(n)
    if f.is_integer():
        return str(int(f))
    s = f'{f:.2f}' if abs(f) >= 1 else f'{f:.2g}'
    if '.' in s and 'e' not in s.lower():
        s = s.rstrip('0').rstrip('.')
    return s


def _value_display(value: float, comparison: str, threshold: float, unit: str) -> str:
    """Plain-language value vs threshold, using the comparison direction."""
    v = f'{_fmt_num(value)}{unit}'
    t = f'{_fmt_num(threshold)}{unit}'
    if comparison == '<':
        return f'{v} (below minimum {t})'
    if comparison == '>':
        return f'{v} (above maximum {t})'
    return f'{v} {comparison} {t}'


# Severity ordering for display: failures rank ahead of warnings.
SEVERITY_RANK: dict[str, int] = {'fail': 0, 'warn': 1}


def _flag_to_dict(flag: QcFlag, source: str) -> dict:
    """Flatten a QcFlag into a template-ready, human-readable dict."""
    label, unit, multiplier = _metric_label(flag.flag)
    # Active flags carry their detection date; resolved carry the resolution date.
    date_full = (flag.resolution_date if flag.resolved else flag.date) or ''
    severity = flag.severity or 'fail'
    return {
        'source': source,
        'flag': flag.flag,
        'metric_label': label,
        'section': flag.section,
        'section_label': _section_label(flag.section),
        'resolved': flag.resolved,
        'severity': severity,
        'severity_label': 'Fail' if severity == 'fail' else 'Warn',
        'value_display': _value_display(flag.value * multiplier, flag.comparison, flag.threshold * multiplier, unit),
        'ar_guid': flag.ar_guid,
        'date_full': date_full,
        'date_short': date_full[:10],  # YYYY-MM-DD
    }


# ---------------------------------------------------------------------------
# Row / section building
# ---------------------------------------------------------------------------
def _make_row(report: SGReport, flags: list[dict], *, resolved: bool) -> dict:
    """Build a template-ready row for one SG within one section."""
    info = report.sg_info
    n_cram = sum(1 for f in flags if f['source'] == 'CRAM')
    n_gvcf = sum(1 for f in flags if f['source'] == 'GVCF')
    parts = []
    if n_cram:
        parts.append(f'{n_cram} CRAM')
    if n_gvcf:
        parts.append(f'{n_gvcf} GVCF')

    noun = 'resolved' if resolved else 'active'
    count_summary = f'{len(flags)} {noun} flag' + ('' if len(flags) == 1 else 's')

    n_fail = sum(1 for f in flags if f['severity'] == 'fail')
    n_warn = sum(1 for f in flags if f['severity'] == 'warn')
    # A row's overall severity is its worst flag - drives row-level styling/sorting.
    row_severity = 'fail' if n_fail else 'warn'
    sev_parts = []
    if n_fail:
        sev_parts.append(f'{n_fail} fail')
    if n_warn:
        sev_parts.append(f'{n_warn} warn')

    return {
        'info': info,
        'flags': flags,
        'source_summary': ' · '.join(parts),
        'count_summary': count_summary,
        'severity_summary': ' · '.join(sev_parts),
        'n_fail': n_fail,
        'n_warn': n_warn,
        'row_severity': row_severity,
        # Group by family, then participant, then SG (collaborator-facing order),
        # but surface failing SGs above warn-only ones within the active section.
        'sort_key': (
            0 if n_fail else 1,
            info.family_external_id or '~',
            info.participant_external_id or '~',
            info.sg_id,
        ),
    }


def build_sections(reports: list[SGReport]) -> tuple[list[dict], list[dict]]:
    """Split reports into (unresolved rows, resolved rows).

    An SG with both active and resolved flags appears in both lists, showing
    only the flags relevant to each section.
    """
    unresolved, resolved = [], []
    for report in reports:
        all_flags = [_flag_to_dict(f, 'CRAM') for f in report.cram_flags]
        all_flags += [_flag_to_dict(f, 'GVCF') for f in report.gvcf_flags]

        active = [f for f in all_flags if not f['resolved']]
        past = [f for f in all_flags if f['resolved']]

        if active:
            # Failures first, then by source and detection date, so the two flags
            # shown in the "at a glance" row lead with the most serious.
            active.sort(key=lambda f: (SEVERITY_RANK.get(f['severity'], 0), f['source'], f['date_full']))
            unresolved.append(_make_row(report, active, resolved=False))
        if past:
            past.sort(key=lambda f: f['date_full'], reverse=True)  # most recent first
            resolved.append(_make_row(report, past, resolved=True))

    unresolved.sort(key=lambda row: row['sort_key'])
    resolved.sort(key=lambda row: row['sort_key'])
    return unresolved, resolved


# ---------------------------------------------------------------------------
# Metamist querying
# ---------------------------------------------------------------------------
def _primary_external_id(obj: dict | None) -> str:
    """Metamist keys the primary external ID under ''. Fall back gracefully."""
    ext = (obj or {}).get('externalIds') or {}
    return ext.get('') or next(iter(ext.values()), '')


def _basename(entry: dict | str | None) -> str | None:
    """Pull a file basename from a reads entry (dict or path string)."""
    if isinstance(entry, dict):
        return entry.get('basename') or ((entry.get('location') or '').rsplit('/', 1)[-1] or None)
    if isinstance(entry, str):
        return entry.rsplit('/', 1)[-1] or None
    return None


def _extract_reads(assays: list[dict]) -> tuple[list[str], list[tuple[str, str]], list[str]]:
    """Group an SG's assay read files into (crams, fastq_pairs, other).

    Each assay's ``meta.reads`` is a list of file entries; ``meta.reads_type``
    tells us whether they're fastq (R1/R2 pairs) or aligned (bam/cram).
    """
    crams: list[str] = []
    fastq_pairs: list[tuple[str, str]] = []
    other: list[str] = []

    for assay in assays:
        meta = assay.get('meta') or {}
        reads = meta.get('reads')
        reads_type = (meta.get('reads_type') or '').lower()
        entries = reads if isinstance(reads, list) else ([reads] if reads else [])
        names = [n for n in (_basename(e) for e in entries) if n]
        if not names:
            continue

        if reads_type == 'fastq':
            # A fastq assay is one (or more) R1/R2 pair(s); pair sequentially.
            for i in range(0, len(names), 2):
                r2 = names[i + 1] if i + 1 < len(names) else ''
                fastq_pairs.append((names[i], r2))
        elif reads_type in ('bam', 'cram'):
            crams.extend(names)
        else:
            # Unknown reads_type — classify by extension.
            for name in names:
                low = name.lower()
                if low.endswith(('.cram', '.bam')):
                    crams.append(name)
                elif low.endswith(('.fastq.gz', '.fq.gz', '.fastq', '.fq')):
                    fastq_pairs.append((name, ''))
                else:
                    other.append(name)

    return crams, fastq_pairs, other


def get_sg_infos(sg_ids: list[str]) -> dict[str, SGInfo]:
    """Query Metamist for detailed SG info, keyed by SG id."""
    if not sg_ids:
        logger.warning('No SGs flagged')
        return {}
    logger.info(f'Querying Metamist for detailed info on {len(sg_ids)} SG(s): {", ".join(sg_ids)}')
    started = perf_counter()
    response = query(SGS_INFO_QUERY, variables={'sgIds': sg_ids})
    logger.info(f'Received SG info for {len(sg_ids)} SG(s) in {perf_counter() - started:.1f}s')
    infos: dict[str, SGInfo] = {}
    for sg in response['sequencingGroups']:
        sample = sg.get('sample') or {}
        participant = sample.get('participant') or {}
        families = participant.get('families') or []
        family = families[0] if families else None

        crams, fastq_pairs, other = _extract_reads(sg.get('assays', []))

        infos[sg['id']] = SGInfo(
            sg_id=sg['id'],
            sg_type=sg.get('type') or '',
            sg_technology=sg.get('technology') or '',
            sg_platform=sg.get('platform') or '',
            crams=crams,
            fastq_pairs=fastq_pairs,
            other_reads=other,
            sample_external_id=_primary_external_id(sample),
            sample_type=sample.get('type') or '',
            participant_external_id=_primary_external_id(participant),
            family_external_id=_primary_external_id(family),
        )
    return infos


def get_previous_analysis(dataset: str, meta_filter: dict) -> dict | None:
    """Query Metamist for existing web analyses matching a meta filter and take the most recent one, if any."""
    response = query(EXISTING_ANALYSES_QUERY, variables={'dataset': dataset, 'metaFilter': meta_filter})
    existing_analyses = response['project']['analyses']
    if not existing_analyses:
        return None
    existing_analyses.sort(key=lambda a: a.get('timestampCompleted') or '', reverse=True)
    if len(existing_analyses) == 1:
        # Only one analysis is the current one, so there are no previous analyses to compare against.
        return None
    previous_analysis = existing_analyses[1]  # The second most recent analysis is the previous one
    logger.info(f'Found previous analysis {previous_analysis["id"]} from {previous_analysis["timestampCompleted"]}')
    if not previous_analysis['meta'].get('summary'):
        logger.warning(f'Previous analysis {previous_analysis["id"]} has no summary in meta; skipping')
        return None
    return previous_analysis


def collect_qc_flags(sequencing_groups: list[dict]) -> list[dict]:
    """Extract QC flags from each sequencing group's metadata."""
    results = []
    for sg in sequencing_groups:
        meta = sg.get('meta') or {}
        results.append(
            {
                'id': sg['id'],
                'cram_qc_flags': [QcFlag(**flag) for flag in meta.get('cram_qc_flags', [])],
                'gvcf_qc_flags': [QcFlag(**flag) for flag in meta.get('gvcf_qc_flags', [])],
            }
        )
    return results


def summarise_flags(sg_data: list[dict]) -> dict:
    """Dataset-wide, flag-centric summary counts for the header cards."""
    all_flags = [(f, 'CRAM') for sg in sg_data for f in sg['cram_qc_flags']]
    all_flags += [(f, 'GVCF') for sg in sg_data for f in sg['gvcf_qc_flags']]

    active = [(f, src) for f, src in all_flags if not f.resolved]
    active_cram = sum(1 for f, src in active if src == 'CRAM')
    active_gvcf = sum(1 for f, src in active if src == 'GVCF')
    active_fail = sum(1 for f, _ in active if (f.severity or 'fail') == 'fail')
    active_warn = sum(1 for f, _ in active if (f.severity or 'fail') == 'warn')
    resolved_flags = sum(1 for f, _ in all_flags if f.resolved)
    sgs_affected = sum(1 for sg in sg_data if _has_active(sg['cram_qc_flags']) or _has_active(sg['gvcf_qc_flags']))

    return {
        'total_sgs': len(sg_data),
        'active_flags': active_cram + active_gvcf,
        'active_cram': active_cram,
        'active_gvcf': active_gvcf,
        'active_fail': active_fail,
        'active_warn': active_warn,
        'sgs_affected': sgs_affected,
        'resolved_flags': resolved_flags,
    }


def metric_histogram(rows: list[dict]) -> list[dict]:
    """Per-metric SG counts for the filter chips (one count per SG per metric)."""
    counts: dict[str, dict] = {}
    for row in rows:
        for key in {f['flag'] for f in row['flags']}:
            entry = counts.setdefault(key, {'key': key, 'label': _metric_label(key)[0], 'count': 0})
            entry['count'] += 1
    return sorted(counts.values(), key=lambda d: (-d['count'], d['label']))


def source_histogram(rows: list[dict]) -> list[dict]:
    """Per-source (CRAM/GVCF) SG counts for the filter toggles."""
    counts: dict[str, int] = {}
    for row in rows:
        for source in {f['source'] for f in row['flags']}:
            counts[source] = counts.get(source, 0) + 1
    return [{'key': k, 'count': counts[k]} for k in ('CRAM', 'GVCF') if k in counts]


def severity_histogram(rows: list[dict]) -> list[dict]:
    """Per-severity SG counts for the filter toggles (one count per SG per tier)."""
    counts: dict[str, int] = {}
    for row in rows:
        for severity in {f['severity'] for f in row['flags']}:
            counts[severity] = counts.get(severity, 0) + 1
    labels = {'fail': 'Failing', 'warn': 'Warnings'}
    return [{'key': k, 'label': labels[k], 'count': counts[k]} for k in ('fail', 'warn') if k in counts]


def render_report(dataset: str, reports: list[SGReport], *, summary: dict) -> str:
    """Build template context and render the Jinja template.

    Only SGs with at least one flag appear in ``reports``; ``summary`` holds the
    dataset-wide, flag-centric counts for the header cards.
    """
    unresolved, resolved = build_sections(reports)

    env = jinja2.Environment(loader=jinja2.FileSystemLoader(JINJA_TEMPLATE_DIR), autoescape=True)
    template = env.get_template('sg_qc_overview.html.jinja')
    return template.render(
        dataset=dataset,
        generated_at=datetime.now().strftime('%Y-%m-%d %H:%M:%S'),  # noqa: DTZ005
        summary=summary,
        unresolved=unresolved,
        resolved=resolved,
        active_metrics=metric_histogram(unresolved),
        active_sources=source_histogram(unresolved),
        active_severities=severity_histogram(unresolved),
    )


def construct_summary_message(
    dataset: str, out_html_url: str, seq_type: str, seq_tech: str, summary: dict, previous_analysis: dict | None
):
    """Construct a Slack message with a concise summary and a link to the report."""
    report_title = f'SG QC report ({seq_type} | {seq_tech})'
    messages = [f'*[{dataset}]* <{out_html_url}|{report_title}>']
    if summary['sgs_affected'] == 0:
        messages.append('✅ No sequencing groups flagged')
    else:
        messages.append(f'{summary["sgs_affected"]} / {summary["total_sgs"]} sequencing groups flagged')
        messages.append(f'🚩 *{summary["active_flags"]} Total active flags*')
        if summary['active_warn']:
            messages.append(f'⚠️ {summary["active_warn"]} warning flags')
        if summary['active_fail']:
            messages.append(f'❗ {summary["active_fail"]} failure flags')

    if previous_analysis:
        previous_summary = previous_analysis['meta']['summary']  # This exists because we already checked it did
        additional_flags = summary['active_flags'] - previous_summary.get('active_flags', 0)
        additional_fail = summary['active_fail'] - previous_summary.get('active_fail', 0)
        additional_warn = summary['active_warn'] - previous_summary.get('active_warn', 0)
        additional_flagged_sgs = summary['sgs_affected'] - previous_summary.get('sgs_affected', 0)
        additional_sgs = summary['total_sgs'] - previous_summary.get('total_sgs', 0)
        timestamp_str = previous_analysis['timestampCompleted'].split('T')[0]  # Extract date portion
        if additional_flags > 0 or additional_flagged_sgs > 0 or additional_sgs > 0:
            messages.append(f'📢 *New since last report on {timestamp_str}*')
            if additional_sgs > 0:
                messages.append(f'+{additional_sgs} new sequencing groups added to the dataset')
            if additional_flagged_sgs > 0:
                messages.append(f'+{additional_flagged_sgs} additional sequencing groups flagged')
            if additional_flags > 0:
                messages.append(f'+{additional_flags} new flags ({additional_fail} fail, {additional_warn} warn)')
        else:
            messages.append(f'No new flags since last report on {timestamp_str}')

    text = '\n'.join(messages)
    logger.info(text)
    if config_retrieve(['workflow', 'sg_qc_report', 'send_to_slack']):
        send_message(text)


def main(dataset: str, output: str, timestamped_output: str, out_html_url: str):
    """Query Metamist for QC flags and generate a SG QC HTML report."""

    dataset = dataset_for_access_level(dataset)
    seq_type = config_retrieve(['workflow', 'sequencing_type'])
    seq_tech = config_retrieve(['workflow', 'sequencing_technology'])

    logging_prefix = f'{dataset} ({seq_type} | {seq_tech})'

    logger.info(f'{logging_prefix} :: Querying Metamist for QC flags')
    started = perf_counter()
    response = query(DATASET_SGS_QUERY, variables={'dataset': dataset, 'seqType': seq_type, 'seqTech': seq_tech})
    sequencing_groups = response['project']['sequencingGroups']
    logger.info(
        f'{logging_prefix} :: Found {len(sequencing_groups)} sequencing groups in {perf_counter() - started:.1f}s'
    )

    sg_data = collect_qc_flags(sequencing_groups)
    summary = summarise_flags(sg_data)

    # Only fetch rich SG metadata for SGs that actually have flags to report.
    flagged = [sg for sg in sg_data if sg['cram_qc_flags'] or sg['gvcf_qc_flags']]
    logger.info(f'{logging_prefix} :: {len(flagged)} sequencing groups have QC flags.')

    infos = get_sg_infos([sg['id'] for sg in flagged])
    reports = [
        SGReport(
            sg_info=infos[sg['id']],
            cram_flags=sg['cram_qc_flags'],
            gvcf_flags=sg['gvcf_qc_flags'],
        )
        for sg in flagged
        if sg['id'] in infos
    ]

    logger.info(f'{logging_prefix} :: Rendering report for {len(reports)} flagged SG(s)')
    started = perf_counter()
    html = render_report(dataset, reports, summary=summary)
    logger.info(f'{logging_prefix} :: Rendered report in {perf_counter() - started:.1f}s')

    with to_path(output).open('w') as f:
        f.write(html)
    logger.info(f'{logging_prefix} :: Wrote SG QC report to {output}')

    with to_path(timestamped_output).open('w') as f:
        f.write(html)
    logger.info(f'{logging_prefix} :: Wrote timestamped SG QC report to {timestamped_output}')

    # Register results in Metamist manually to capture all dataset SGs in scope, not just the input_cohorts SGs
    meta = {
        'stage': 'GenerateSgQcReport',
        'dataset': dataset,
        'sequencing_type': seq_type,
        'sequencing_technology': seq_tech,
        'summary': summary,
    }

    create_new(
        project=dataset,
        output=timestamped_output,
        analysis_type='web',
        sgs=[sg['id'] for sg in sequencing_groups],
        meta=meta,
    )
    logger.info(f'{logging_prefix} :: Registered web analysis for {len(sequencing_groups)} SG(s)')

    meta.pop('summary')
    construct_summary_message(dataset, out_html_url, seq_type, seq_tech, summary, get_previous_analysis(dataset, meta))


if __name__ == '__main__':
    parser = ArgumentParser()
    parser.add_argument('--dataset', required=True, help='Metamist dataset/project name')
    parser.add_argument('--fixed-output', required=True, help='Path to write the HTML report')
    parser.add_argument('--timestamped-output', required=True, help='Path to write the timestamped HTML report')
    parser.add_argument('--html-url', required=True, help='Clickable URL for the HTML report')
    args = parser.parse_args()
    main(
        dataset=args.dataset,
        output=args.fixed_output,
        timestamped_output=args.timestamped_output,
        out_html_url=args.html_url,
    )
