"""
Verify what causes Manta's "0 high-confidence read pairs" failure.

Checks reference compatibility, extracts a region to BAM, and attempts
Manta on both CRAM and BAM to isolate the issue. Results in Hail Batch job log.

Usage:
    analysis-runner ... verify_fixmate_manta --cram_path gs://bucket/sample.cram
"""

import argparse

from cpg_utils import config, hail_batch

SCRIPT_PATH = '/align_genotype/src/align_genotype/scripts/verify_fixmate_manta.sh'


def main(cram_path: str) -> None:
    batch = hail_batch.get_batch()

    job = batch.new_job('Verify CRAM Manta compatibility')
    job.image(config.config_retrieve(['workflow', 'driver_image']))
    job.cpu(4)
    job.memory('16Gi')
    job.storage('50Gi')

    cram = batch.read_input_group(
        cram=cram_path,
        crai=f'{cram_path}.crai',
    ).cram

    reference = hail_batch.fasta_res_group(batch)

    job.command(f'bash {SCRIPT_PATH} {cram} {reference.base}')

    batch.run(wait=False)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Verify CRAM compatibility with Manta')
    parser.add_argument('--cram_path', type=str, required=True, help='GCS path to CRAM file (.crai must exist alongside)')
    args = parser.parse_args()
    main(cram_path=args.cram_path)