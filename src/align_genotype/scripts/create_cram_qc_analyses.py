"""
Retrospectively register `cram_qc` analyses in Metamist for CRAM QC files that cpg-flow wrote but never registered.

cpg-flow (stages.py) writes each QC files to a path like <dataset prefix>/qc/<tool>/<sg>.<suffix>, e.g.
    gs://cpg-rdnow-main/qc/samtools_stats/CPGNNN.samtools-stats

where the dataset prefix is the dataset's bucket (gs://cpg-<dataset>-main), plus exome/ for exomes.

nf-align (nf-align/main.nf) registers one `cram_qc` analysis per SG: the CRAM as the primary output, and every QC
file as a secondary file, keyed by its type. This script builds the same analysis for the QC files cpg-flow left
behind, using the same keys.

For each project, every active SG of the requested sequencing types that has a registered CRAM, and no `cram_qc`
analysis yet, is checked for its QC files:
- every required QC file exists: the analysis is created. The optional duplicate metrics (written at alignment, and
  not by every aligner) are registered if they exist, and a warning is logged if they don't
- some required files exist, but not all: an error. Nothing is registered for the SG, and the script exits non-zero
  at the end, after registering the complete ones
- no required files exist: the SG was never QC'd, and is skipped

    python3 -m align_genotype.scripts.create_cram_qc_analyses --dry-run cram_qc_dry_run.json
    python3 -m align_genotype.scripts.create_cram_qc_analyses --projects rdnow acute-care
"""

import argparse
import json
import logging
import sys
from dataclasses import asdict, dataclass, field

# a functools.cache decorated existence checker, running one query per directory
from cpg_utils.existence_checks import exists
from cpg_utils.metamist_registration import create_new
from metamist.apis import ProjectApi
from metamist.graphql import gql, query

ANALYSIS_TYPE = 'cram_qc'

logger = logging.getLogger(__name__)

SEQUENCING_TYPES = ('genome', 'exome')

# QC files an analysis is registered without: flagged when missing, but never an error
OPTIONAL_KEYS = frozenset({'markduplicates_metrics'})

# every active SG of one type in a project, with its active CRAM and CRAM QC analyses
PROJECT_SGS_QUERY = gql(
    """
    query ProjectSgAnalyses($project: String!, $seq_type: String!) {
        project(name: $project) {
            sequencingGroups(type: {eq: $seq_type}, activeOnly: {eq: true}) {
                id
                analyses(type: {in_: ["cram", "cram_qc"]}, status: {eq: COMPLETED}, active: {eq: true}) {
                    id
                    type
                    outputs
                }
            }
        }
    }
    """
)


@dataclass
class CramQcRecord:
    """One SG's cram_qc analysis: the CRAM it describes, and the QC files found for it."""

    sg_id: str
    project: str
    sequencing_type: str
    cram: str
    secondary: dict[str, str] = field(default_factory=dict)
    missing: list[str] = field(default_factory=list)  # required QC files not found: the record can't be registered
    missing_optional: list[str] = field(default_factory=list)  # optional QC files not found: flagged only

    @property
    def meta(self) -> dict[str, str]:
        return {'sequencing_type': self.sequencing_type, 'source': 'create_cram_qc_analyses'}


def bucket_name(project: str) -> str:
    """A Metamist project's main bucket: cpg-<dataset>-main, or cpg-<dataset>-test for a <dataset>-test project."""
    if project.endswith('-test'):
        return f'gs://cpg-{project.removesuffix("-test")}-test'
    return f'gs://cpg-{project}-main'


def qc_prefix(project: str, sequencing_type: str) -> str:
    """Where cpg-flow writes QC for a project: the dataset prefix (with exome/ for exomes), then qc/."""
    subdir = 'exome/' if sequencing_type == 'exome' else ''
    return f'{bucket_name(project)}/{subdir}qc'


def expected_qc_files(sg_id: str, project: str, sequencing_type: str) -> dict[str, str]:
    """
    Every QC file cpg-flow can write for an SG, keyed by its secondary-file key in a cram_qc analysis.
    Paths match stages.py's expected_outputs; keys match nf-align's.
    """
    prefix = qc_prefix(project, sequencing_type)
    coverage = 'wgs' if sequencing_type == 'genome' else 'hs'
    # key -> (directory under qc/, file suffix)
    layout = {
        'verifybamid': ('verify_bamid', 'verify-bamid.selfSM'),
        'samtools_stats': ('samtools_stats', 'samtools-stats'),
        'alignment_summary_metrics': ('alignment_summary_metrics', 'alignment_summary_metrics'),
        'base_distribution_by_cycle_metrics': (
            'base_distribution_by_cycle_metrics',
            'base_distribution_by_cycle_metrics',
        ),
        'insert_size_metrics': ('insert_size_metrics', 'insert_size_metrics'),
        'quality_by_cycle_metrics': ('quality_by_cycle_metrics', 'quality_by_cycle_metrics'),
        'quality_yield_metrics': ('quality_yield_metrics', 'quality_yield_metrics'),
        f'picard_{coverage}_metrics': (f'picard_{coverage}_metrics', f'picard-{coverage}-metrics'),
        'markduplicates_metrics': ('markduplicates_metrics', 'markduplicates-metrics'),
    }
    return {key: f'{prefix}/{directory}/{sg_id}.{suffix}' for key, (directory, suffix) in layout.items()}


def output_path(outputs: dict | str | None) -> str | None:
    """The primary path from Metamist's polymorphic analysis outputs field."""
    if isinstance(outputs, dict):
        return outputs.get('path')
    return outputs if isinstance(outputs, str) else None


def relevant_projects(projects: list[str], test_only: bool = False) -> list[str]:
    """Drop any project whose name contains an excluded name."""
    def conditional(project: str, test_only: bool) -> bool:
        if (test_only and '-test' in project) or (not test_only and '-test' not in project):
            return True
        return False
    return sorted(p for p in projects if conditional(p, test_only))


def find_records(project: str, sequencing_type: str) -> list[CramQcRecord]:
    """A record for each SG in the project that has a CRAM, no cram_qc analysis yet, and at least one QC file."""
    result = query(PROJECT_SGS_QUERY, variables={'project': project, 'seq_type': sequencing_type})
    records = []
    for sg in result['project']['sequencingGroups']:
        sg_id = sg['id']
        analyses = sg['analyses']

        if any(a['type'] == ANALYSIS_TYPE for a in analyses):
            logger.debug(f'{project} {sg_id}: already has a {ANALYSIS_TYPE} analysis')
            continue

        crams = sorted((a for a in analyses if a['type'] == 'cram'), key=lambda a: a['id'])
        cram = output_path(crams[-1]['outputs']) if crams else None
        if not cram:
            logger.debug(f'{project} {sg_id}: no CRAM analysis')
            continue
        if not exists(cram):
            logger.warning(f'{project} {sg_id}: registered CRAM {cram} does not exist, skipping')
            continue

        record = CramQcRecord(sg_id=sg_id, project=project, sequencing_type=sequencing_type, cram=cram)
        for key, path in expected_qc_files(sg_id, project, sequencing_type).items():
            if exists(path):
                record.secondary[key] = path
            elif key in OPTIONAL_KEYS:
                record.missing_optional.append(key)
            else:
                record.missing.append(key)

        # duplicate metrics are optional in CPG-Flow
        if not record.secondary.keys() - OPTIONAL_KEYS:
            logger.debug(f"{project} {sg_id}: no required QC files found, never QC'd")
            continue
        records.append(record)
    return records


def register(record: CramQcRecord) -> None:
    """Create the cram_qc analysis in Metamist. create_new fails if any file doesn't exist."""
    if record.missing:
        raise ValueError(f'{record.project} {record.sg_id}: refusing to register an incomplete QC set')
    create_new(
        project=record.project,
        output=record.cram,
        analysis_type=ANALYSIS_TYPE,
        meta=record.meta,
        secondary=record.secondary,
        sgs=[record.sg_id],
    )


def main(projects: list[str], sequencing_types: list[str], dry_run: str | None) -> list[CramQcRecord]:
    """
    Registers (or, for a dry run, writes out) every SG with a complete QC set.
    Returns the SGs with an incomplete set, which are never registered.
    """
    records: list[CramQcRecord] = []
    for project in projects:
        for sequencing_type in sequencing_types:
            found = find_records(project, sequencing_type)
            logger.info(f'{project} {sequencing_type}: {len(found)} SGs with unregistered QC files')
            records.extend(found)

    complete = [r for r in records if not r.missing]
    incomplete = [r for r in records if r.missing]
    for record in incomplete:
        logger.error(
            f'{record.project} {record.sg_id} ({record.sequencing_type}): incomplete QC set, not registering. '
            f'Missing: {", ".join(record.missing)}'
        )
    without_optional = [r for r in complete if r.missing_optional]
    for record in without_optional:
        logger.warning(
            f'{record.project} {record.sg_id} ({record.sequencing_type}): no optional '
            f'{", ".join(record.missing_optional)}, registering without'
        )
    logger.info(
        f'{len(complete)} cram_qc analyses to create ({len(without_optional)} without optional files), '
        f'{len(incomplete)} SGs with incomplete QC sets'
    )

    if dry_run:
        with open(dry_run, 'w') as handle:
            json.dump(
                {
                    'to_register': [asdict(r) | {'meta': r.meta} for r in complete],
                    'errors': [asdict(r) for r in incomplete],
                },
                handle,
                indent=2,
            )
        logger.info(f'Dry run: wrote {len(complete)} records and {len(incomplete)} errors to {dry_run}')
        return incomplete

    for record in complete:
        register(record)
        logger.info(f'Registered {ANALYSIS_TYPE} for {record.project} {record.sg_id}')
    return incomplete


def cli_main() -> None:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument(
        '--projects',
        nargs='+',
        help='Metamist projects to scan. Default: every project you can access, minus the excluded ones.',
    )
    parser.add_argument(
        '--sequencing-types',
        nargs='+',
        choices=SEQUENCING_TYPES,
        default=list(SEQUENCING_TYPES),
    )
    parser.add_argument('--test-only', action='store_true', help='Only run test datasets')
    parser.add_argument(
        '--dry-run',
        metavar='JSON',
        help='Write the SGs and QC files found to this file, instead of creating analyses in Metamist.',
    )
    parser.add_argument('--verbose', action='store_true', help='Log every SG skipped, and why.')
    args = parser.parse_args()

    logging.basicConfig(level=logging.INFO, format='%(levelname)s: %(message)s')
    logger.setLevel(logging.DEBUG if args.verbose else logging.INFO)
    projects = args.projects or relevant_projects(ProjectApi().get_my_projects(), args.test_only)
    logger.info(f'Scanning {len(projects)} projects: {", ".join(projects)}')
    if incomplete := main(projects, args.sequencing_types, args.dry_run):
        sys.exit(f'{len(incomplete)} SGs have incomplete QC sets and were not registered, see the errors above')


if __name__ == '__main__':
    cli_main()
