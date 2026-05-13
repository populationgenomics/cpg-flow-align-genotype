"""
Hail Batch job to diagnose CRAM issues causing Manta's "0 high-confidence read pairs" failure.

Submits a batch job that localises the CRAM+CRAI and runs diagnose_cram_manta.sh.
Results appear in the Hail Batch job log.

Usage:
    analysis-runner ... diagnose_cram_manta --cram_path gs://bucket/sample.cram
"""

import argparse

from cpg_utils import config, hail_batch

SCRIPT_PATH = '/align_genotype/src/align_genotype/scripts/diagnose_cram_manta.sh'


def main(cram_path: str) -> None:
    batch = hail_batch.get_batch()

    job = batch.new_job('Diagnose CRAM for Manta')
    job.image(config.config_retrieve(['workflow', 'driver_image']))
    job.cpu(2)
    job.memory('8Gi')
    job.storage('100Gi')

    cram = batch.read_input_group(
        cram=cram_path,
        crai=f'{cram_path}.crai',
    ).cram

    job.command(f'bash {SCRIPT_PATH} {cram}')

    batch.run(wait=False)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Diagnose CRAM issues causing Manta failure')
    parser.add_argument('--cram_path', type=str, required=True, help='GCS path to CRAM file (.crai must exist alongside)')
    args = parser.parse_args()
    main(cram_path=args.cram_path)