"""
Queries Metamist for all QC flags across a dataset's sequencing groups
and renders the master_qc.html.jinja template.
"""

from argparse import ArgumentParser
from datetime import datetime
from pathlib import Path

import jinja2
from loguru import logger
from metamist.graphql import gql, query

from cpg_utils.config import config_retrieve


JINJA_TEMPLATE_DIR = Path(__file__).absolute().parent.parent / 'templates'

DATASET_QC_FLAGS_QUERY = gql(
    """
    query datasetQcFlags($dataset: String!) {
        project(name: $dataset) {
            sequencingGroups {
                id
                meta
                sample {
                    externalId
                }
            }
        }
    }
    """
)


def _has_active(flags: list[dict]) -> bool:
    return any(not f.get('resolved', False) for f in flags)


def _prepare_sg_row(sg: dict) -> dict:
    """Transform raw SG data into a template-ready row dict."""
    cram_flags = sg['cram_flags']
    gvcf_flags = sg['gvcf_flags']
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

    display_id = sg['id'] + (f' ({sg["external_id"]})' if sg['external_id'] else '')

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


def collect_qc_flags(sequencing_groups: list[dict]) -> list[dict]:
    """Extract QC flags from each sequencing group's metadata."""
    results = []
    for sg in sequencing_groups:
        meta = sg.get('meta') or {}
        external_id = (sg.get('sample') or {}).get('externalId', '')
        results.append(
            {
                'id': sg['id'],
                'external_id': external_id,
                'cram_flags': meta.get('cram_qc_flags', []),
                'gvcf_flags': meta.get('gvcf_qc_flags', []),
            }
        )
    return results


def render_report(dataset: str, sg_data: list[dict]) -> str:
    """Build template context and render the Jinja template."""
    total = len(sg_data)
    flagged_cram = sum(1 for sg in sg_data if _has_active(sg['cram_flags']))
    flagged_gvcf = sum(1 for sg in sg_data if _has_active(sg['gvcf_flags']))
    flagged_any = sum(1 for sg in sg_data if _has_active(sg['cram_flags']) or _has_active(sg['gvcf_flags']))

    rows = [_prepare_sg_row(sg) for sg in sorted(sg_data, key=lambda s: s['id'])]

    env = jinja2.Environment(loader=jinja2.FileSystemLoader(JINJA_TEMPLATE_DIR), autoescape=True)
    template = env.get_template('master_qc.html.jinja')
    return template.render(
        dataset=dataset,
        generated_at=datetime.now().strftime('%Y-%m-%d %H:%M:%S'),  # noqa: DTZ005
        total=total,
        flagged_any=flagged_any,
        flagged_cram=flagged_cram,
        flagged_gvcf=flagged_gvcf,
        sequencing_groups=rows,
    )


def main(dataset: str, output_html: str):
    """Query Metamist for QC flags and generate a master QC HTML report."""

    if config_retrieve(['workflow', 'access_level'], 'standard') == 'test' and '-test' not in dataset:
        dataset = f'{dataset}-test'

    logger.info(f'Querying Metamist for QC flags in dataset: {dataset}')
    response = query(DATASET_QC_FLAGS_QUERY, variables={'dataset': dataset})
    sequencing_groups = response['project']['sequencingGroups']
    logger.info(f'Found {len(sequencing_groups)} sequencing groups.')

    sg_data = collect_qc_flags(sequencing_groups)
    flagged = [sg for sg in sg_data if sg['cram_flags'] or sg['gvcf_flags']]
    logger.info(f'{len(flagged)} sequencing groups have QC flags.')

    html = render_report(dataset, sg_data)

    with open(output_html, 'w') as f:
        f.write(html)
    logger.info(f'Wrote master QC report to {output_html}')


if __name__ == '__main__':
    parser = ArgumentParser()
    parser.add_argument('--dataset', required=True, help='Metamist dataset/project name')
    parser.add_argument('--output', required=True, help='Path to write the HTML report')
    args = parser.parse_args()
    main(args.dataset, args.output_html)
