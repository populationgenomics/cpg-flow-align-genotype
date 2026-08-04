"""
Queries Metamist for all QC flags across a dataset's sequencing groups
and renders the sg_qc_overview.html.jinja template.
"""

import asyncio
from argparse import ArgumentParser
from dataclasses import dataclass
from datetime import datetime
from pathlib import Path

import jinja2
from loguru import logger

from cpg_utils.config import dataset_for_access_level
from metamist.graphql import gql, query, query_async

from align_genotype.utils import QcFlag

JINJA_TEMPLATE_DIR = Path(__file__).absolute().parent.parent / 'templates'

DATASET_SGS_QUERY = gql(
    """
    query datasetSgs($dataset: String!) {
        project(name: $dataset) {
            sequencingGroups {
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
    query sgInfo($sgId: String!) {
        sequencingGroup(id: $sgId) {
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


@dataclass
class SGInfo:
    sg_id: str
    sg_type: str
    sg_technology: str
    sg_platform: str
    read_files: list[str]
    sample_external_id: str
    sample_type: str
    participant_external_id: str
    family_external_id: str


@dataclass
class SGQCFlag:
    sg_info: SGInfo
    flag: QcFlag
    source: str  # 'CRAM' or 'GVCF'


def _has_active(flags: list[QcFlag]) -> bool:
    return any(not f.resolved for f in flags)


def _prepare_sg_rows(sg_data: list[dict[str, SGInfo | list[QcFlag]]]) -> list[SGQCFlag]:
    """
    Extracts each SG's info and QC flags into a SGQCFlag dataclass
    """
    rows = []
    for sg in sg_data:
        sg_info = sg['sg_info']
        cram_flags = sg['cram_qc_flags']
        gvcf_flags = sg['gvcf_qc_flags']
        for flag in cram_flags:
            rows.append(
                SGQCFlag(
                    sg_info=sg_info,
                    flag=flag,
                    source='CRAM',
                )
            )
        for flag in gvcf_flags:
            rows.append(
                SGQCFlag(
                    sg_info=sg_info,
                    flag=flag,
                    source='GVCF',
                )
            )
    return sorted(rows, key=lambda r: (r.sg_info.sg_id, r.source, r.flag.date or ''))


def _prepare_sg_row(sg: SGQCFlag) -> dict:
    """Transform raw SG data into a template-ready row dict."""
    # We need to rework this function to handle the SGQCFlag dataclass instead of the previous dict structure.
    # From this, we want to render the SG info and the QC flags in a way that is suitable for the HTML template.
    # Ideally we should display the SG ID, type, technology, platform, sample external ID, participant external ID,
    # family external ID, and then the QC flags with their details.
    cram_flags = sg.flag if sg.source == 'CRAM' else []
    gvcf_flags = sg.flag if sg.source == 'GVCF' else []
    has_active_cram = _has_active(cram_flags)
    has_active_gvcf = _has_active(gvcf_flags)

    if has_active_cram or has_active_gvcf:
        status_class = 'status-fail'
        status_text = 'Flagged'
    elif cram_flags or gvcf_flags:
        status_class = 'status-resolved'
        status_text = 'Resolved'
    else:
        status_class = 'status-pass'
        status_text = 'Clean'

    display_id = sg.sg_info.sg_id + (
        f' ({sg.sg_info.sample_external_id})' if sg.sg_info.participant_external_id else ''
    )

    flags = []
    for f in cram_flags:
        flags.append({**f, 'source': 'CRAM', 'date': f.get('resolution_date') or f.get('date', '')})
    for f in gvcf_flags:
        flags.append({**f, 'source': 'GVCF', 'date': f.get('resolution_date') or f.get('date', '')})

    n_active_cram = sum(1 for f in cram_flags if not f.get('resolved', False))
    n_active_gvcf = sum(1 for f in gvcf_flags if not f.get('resolved', False))
    n_resolved = (len(cram_flags) + len(gvcf_flags)) - (n_active_cram + n_active_gvcf)

    summary_parts = []
    if n_active_cram:
        summary_parts.append(f'{n_active_cram} CRAM')
    if n_active_gvcf:
        summary_parts.append(f'{n_active_gvcf} GVCF')

    return {
        'display_id': display_id,
        'status_class': status_class,
        'status_text': status_text,
        'active_summary': ', '.join(summary_parts) if summary_parts else '—',
        'resolved_summary': str(n_resolved) if n_resolved else '—',
        'flags': flags,
    }


async def get_sg_info(sg_id: str) -> SGInfo:
    """Query Metamist for detailed SG info."""
    response = query_async(SGS_INFO_QUERY, variables={'sgId': sg_id})
    sg = response['sequencingGroup']
    sample = sg['sample']
    participant = sample['participant']
    family = participant['families'][0] if participant['families'] else None

    read_files = []
    for assay in sg.get('assays', []):
        assay_meta = assay.get('meta') or {}
        reads = assay_meta.get('reads', [])
        if isinstance(reads, dict):
            read_files.append(reads.get('basename'))
        elif isinstance(reads, list):
            for r in reads:
                if isinstance(r, dict):
                    read_files.append(r.get('basename'))
                elif isinstance(r, str):
                    read_files.append(r)

    return SGInfo(
        sg_id=sg['id'],
        sg_type=sg['type'],
        sg_technology=sg['technology'],
        sg_platform=sg['platform'],
        read_files=read_files,
        sample_external_id=sample['externalIds'][''],
        sample_type=sample['type'],
        participant_external_id=participant['externalIds'][''],
        family_external_id=family['externalIds'][''],
    )


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


def render_report(dataset: str, sg_data: list[dict[str, SGInfo | list[QcFlag]]]) -> str:
    """Build template context and render the Jinja template."""
    total = len(sg_data)
    flagged_cram = sum(1 for sg in sg_data if _has_active(sg['cram_qc_flags']))
    flagged_gvcf = sum(1 for sg in sg_data if _has_active(sg['gvcf_qc_flags']))
    flagged_any = sum(1 for sg in sg_data if _has_active(sg['cram_qc_flags']) or _has_active(sg['gvcf_qc_flags']))

    rows = [_prepare_sg_row(sg) for sg in _prepare_sg_rows(sg_data)]

    env = jinja2.Environment(loader=jinja2.FileSystemLoader(JINJA_TEMPLATE_DIR), autoescape=True)
    template = env.get_template('sg_qc_overview.html.jinja')
    return template.render(
        dataset=dataset,
        generated_at=datetime.now().strftime('%Y-%m-%d %H:%M:%S'),  # noqa: DTZ005
        total=total,
        flagged_any=flagged_any,
        flagged_cram=flagged_cram,
        flagged_gvcf=flagged_gvcf,
        sequencing_groups=rows,
    )


def main(dataset: str, output: str):
    """Query Metamist for QC flags and generate a SG QC HTML report."""

    dataset = dataset_for_access_level(dataset)

    logger.info(f'{dataset} :: Querying Metamist for QC flags')
    response = query(DATASET_SGS_QUERY, variables={'dataset': dataset})
    sequencing_groups = response['project']['sequencingGroups']
    logger.info(f'{dataset} :: Found {len(sequencing_groups)} sequencing groups.')

    sg_data = collect_qc_flags(sequencing_groups)
    flagged = [sg for sg in sg_data if sg['cram_qc_flags'] or sg['gvcf_qc_flags']]
    logger.info(f'{dataset} :: {len(flagged)} sequencing groups have QC flags.')
    flagged_sgs_info = [
        {
            'sg_info': asyncio.run(get_sg_info(sg['id'])),
            'cram_qc_flags': sg['cram_qc_flags'],
            'gvcf_qc_flags': sg['gvcf_qc_flags'],
        }
        for sg in flagged
    ]

    html = render_report(dataset, flagged_sgs_info)

    with open(output, 'w') as f:
        f.write(html)
    logger.info(f'{dataset} :: Wrote SG QC report to {output}')


if __name__ == '__main__':
    parser = ArgumentParser()
    parser.add_argument('--dataset', required=True, help='Metamist dataset/project name')
    parser.add_argument('--output', required=True, help='Path to write the HTML report')
    args = parser.parse_args()
    main(dataset=args.dataset, output=args.output)
