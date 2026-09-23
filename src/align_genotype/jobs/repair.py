"""
Entrypoint for CRAM repair workflows.

Dispatches to the appropriate repair script based on --repair-type, handles
metamist registration of the repaired CRAM and inactivation of the old entry.

Usage via analysis-runner:
    analysis-runner \\
        --image <driver-image> \\
        --config src/align_genotype/scripts/repair_scripts/repair_config.toml \\
        --dataset seqr --access-level full \\
        --output-dir "repair/trim-adapters" \\
        repair_crams -- --repair-type trim-adapters --dataset seqr --sg-ids CPG123
"""

import argparse
from collections.abc import Callable

from loguru import logger

from hailtop.batch.job import Job

from cpg_utils import config, hail_batch, to_path

from align_genotype.scripts.repair_scripts import repair_utils, strip_qnames, trim_adapters

REPAIR_FUNCTIONS: dict[str, Callable[..., list[Job]]] = {
    'strip-qnames': strip_qnames.run,
    'trim-adapters': trim_adapters.run,
}


def _output_cram_path(cram_path: str) -> str:
    """Derive the repaired CRAM output path from the original CRAM path.

    Replaces the /cram/ directory with /cram_repaired/ in the same bucket.
    """
    p = to_path(cram_path)
    parts = list(p.parts)
    try:
        idx = parts.index('cram')
        parts[idx] = 'cram_repaired'
    except ValueError:
        parts.insert(-1, 'cram_repaired')
    return str(to_path('/'.join(parts)))


def main() -> None:
    parser = argparse.ArgumentParser(description='Repair CRAM files with known defects.')
    parser.add_argument(
        '--repair-type',
        choices=list(REPAIR_FUNCTIONS),
        required=True,
        help='Type of repair to apply.',
    )
    parser.add_argument(
        '--dataset',
        required=True,
        help='Metamist dataset/project name for analysis registration.',
    )
    parser.add_argument(
        '--sg-ids',
        nargs='+',
        required=True,
        help='Sequencing group IDs to repair (queries metamist for CRAM paths).',
    )
    parser.add_argument(
        '--skip-jobs',
        nargs='+',
        choices=['extract', 'trim', 'realign', 'strip'],
        default=[],
        help='Force-skip specific jobs, assuming their outputs already exist.',
    )
    parser.add_argument('--dry-run', action='store_true', help='Print the CRAM files that would be repaired.')
    args = parser.parse_args()

    sg_crams = repair_utils.get_cram_paths_for_sgs(args.sg_ids)

    if not sg_crams:
        print('No CRAM files found.')
        raise SystemExit(1)

    repair_fn = REPAIR_FUNCTIONS[args.repair_type]

    if args.dry_run:
        for sg_id, cram_path, _analysis_id in sg_crams:
            output = _output_cram_path(cram_path)
            print(f'[{args.repair_type}] {sg_id} {cram_path} -> {output}')
        return

    batch = hail_batch.get_batch()
    for sg_id, cram_path, old_analysis_id in sg_crams:
        output_cram = _output_cram_path(cram_path)

        kwargs: dict = {
            'job_attrs': {'repair_type': args.repair_type},
            'skip_jobs': set(args.skip_jobs),
        }

        repair_jobs = repair_fn(batch, cram_path, sg_id, output_cram, **kwargs)

        if repair_jobs:
            reg_job = batch.new_python_job(
                f'register_repair_{sg_id}',
                attributes={'tool': 'metamist'},
            )
            reg_job.image(config.config_retrieve(['workflow', 'driver_image']))
            for j in repair_jobs:
                reg_job.depends_on(j)
            reg_job.call(
                repair_utils.register_and_inactivate,
                new_cram_path=output_cram,
                sg_id=sg_id,
                dataset=args.dataset,
                old_analysis_id=old_analysis_id,
                repair_type=args.repair_type,
            )
        else:
            logger.info(f'All repair jobs for {sg_id} already completed, skipping registration')

    batch.run(wait=False)


if __name__ == '__main__':
    main()
