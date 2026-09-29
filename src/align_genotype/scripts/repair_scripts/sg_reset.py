"""Retire the stale derived analyses for a sequencing group.

Finds every active analysis for a sequencing group, copies its outputs to an archive
prefix, then inactivates the analysis. Intended to run after a repair script has
produced a corrected CRAM, so the outputs derived from the bad CRAM stop being picked
up by the pipeline.

The sequencing group entry itself is never touched - only analyses.

`cram` analyses are never touched either. Swapping the CRAM a sequencing group points
at is the consequential step - inactivating one without a registered replacement makes
the pipeline realign from scratch - so it is left to be done deliberately by hand. That
also means this is safe to run after a repair has registered its new CRAM.
"""

import argparse
import subprocess

from loguru import logger

from cpg_flow.utils import exists
from cpg_utils import to_path
from metamist.graphql import gql, query

# exact analysis types this tool refuses to touch; swapping a CRAM is done by hand.
# note this is an exact match, so derived types like 'mito-cram' are still reset.
PROTECTED_TYPES = ('cram',)

# index/checksum companions are not recorded in metamist, so probe for them
SIBLING_SUFFIXES = ('.tbi', '.crai', '.bai', '.idx', '.md5')

ACTIVE_ANALYSES_QUERY = gql(
    """
    query SgActiveAnalyses($sg_id: String!) {
        sequencingGroups(id: {eq: $sg_id}) {
            id
            analyses(active: {eq: true}) {
                id
                type
                status
                outputs
                sequencingGroups { id }
            }
        }
    }
    """
)


def _output_path(outputs: dict | str | None) -> str | None:
    """Pull a single usable path out of the polymorphic metamist outputs field."""
    if isinstance(outputs, dict):
        return outputs.get('path')
    return outputs if isinstance(outputs, str) else None


def collect(sg_id: str, include_multi_sample: bool) -> tuple[list[dict], list[tuple[dict, str]]]:
    """Return (targets, skipped) where skipped carries a reason per analysis."""

    result = query(ACTIVE_ANALYSES_QUERY, variables={'sg_id': sg_id})
    targets: list[dict] = []
    skipped: list[tuple[dict, str]] = []

    for sg in result['sequencingGroups']:
        for analysis in sg['analyses']:
            path = _output_path(analysis['outputs'])
            n_sgs = len(analysis['sequencingGroups'])

            if analysis['type'] in PROTECTED_TYPES:
                skipped.append((analysis, 'cram type, handle manually'))
            elif not path:
                skipped.append((analysis, 'no output path'))
            elif n_sgs != 1 and not include_multi_sample:
                skipped.append((analysis, f'multi-sample ({n_sgs} sgs)'))
            else:
                targets.append(analysis | {'path': path})

    return targets, skipped


def files_for(path: str) -> list[str]:
    """The analysis output plus any index/checksum companions that exist beside it."""
    found = [path] if exists(path) else []
    found += [f'{path}{suffix}' for suffix in SIBLING_SUFFIXES if exists(f'{path}{suffix}')]
    return found


def archive_dest(path: str, dest: str) -> str:
    """Archive under <dest>/<original parent dir>/, preserving structure to avoid collisions."""
    p = to_path(path)
    return f'{dest.rstrip("/")}/{p.parent.name}/'


Plan = list[tuple[dict, list[str], str]]


def build_plan(targets: list[dict], dest: str) -> Plan:
    """Pair each target analysis with the files to archive and where they are going."""
    return [(a, files_for(a['path']), archive_dest(a['path'], dest)) for a in targets]


def report(sg_id: str, skipped: list[tuple[dict, str]], plan: Plan) -> None:
    """Print the full set of intended actions, for dry runs and for the record."""
    if skipped:
        print(f'\nSkipping {len(skipped)} analyses:')
        for analysis, reason in skipped:
            print(f'  {analysis["type"]:<12} id={analysis["id"]:<8} [{reason}]')

    print(f'\n{len(plan)} analyses to reset for {sg_id}:')
    for analysis, files, dest_dir in plan:
        print(f'  {analysis["type"]:<12} id={analysis["id"]}')
        if not files:
            print(f'      ! output missing on disk: {analysis["path"]}')
        for f in files:
            print(f'      copy {f}')
            print(f'        -> {dest_dir}{to_path(f).name}')
        print('      inactivate')


def apply_plan(plan: Plan) -> None:
    """Archive each analysis output, verify it landed, then inactivate the record."""
    from metamist.apis import AnalysisApi  # noqa: PLC0415
    from metamist.models import AnalysisUpdateModel  # noqa: PLC0415

    aapi = AnalysisApi()

    for analysis, files, dest_dir in plan:
        if files:
            subprocess.run(['gcloud', 'storage', 'cp', *files, dest_dir], check=True)  # noqa: S603, S607
            for f in files:
                copied = f'{dest_dir}{to_path(f).name}'
                if not exists(copied):
                    msg = f'archive missing after copy: {copied}'
                    raise RuntimeError(msg)
            logger.info(f'archived {len(files)} file(s) for analysis {analysis["id"]} to {dest_dir}')
        else:
            logger.warning(f'analysis {analysis["id"]} has no files on disk, inactivating record only')

        aapi.update_analysis(
            analysis_id=analysis['id'],
            analysis_update_model=AnalysisUpdateModel(
                active=False,
                meta={'archived_to': dest_dir} if files else {},
            ),
        )
        logger.info(f'inactivated analysis {analysis["id"]} ({analysis["type"]})')


def main() -> None:
    parser = argparse.ArgumentParser(
        description='Archive and inactivate the active analyses for a sequencing group.',
    )
    parser.add_argument('--sg-id', required=True, help='Sequencing group ID to reset.')
    parser.add_argument('--dest', required=True, help='GCS prefix to copy analysis outputs to.')
    parser.add_argument(
        '--include-multi-sample',
        action='store_true',
        help='Also reset cohort-level analyses. Off by default: their outputs are shared '
        'with other sequencing groups, so archiving and inactivating them affects other samples.',
    )
    parser.add_argument('--dry-run', action='store_true', help='Report what would be done, change nothing.')
    args = parser.parse_args()

    targets, skipped = collect(args.sg_id, args.include_multi_sample)
    plan = build_plan(targets, args.dest)
    report(args.sg_id, skipped, plan)

    if not plan:
        print(f'\nNothing to reset for {args.sg_id}.')
        return

    if args.dry_run:
        print('\n[dry run] nothing changed.')
        return

    apply_plan(plan)
    print(f'\nReset {len(plan)} analyses for {args.sg_id}.')


if __name__ == '__main__':
    main()
