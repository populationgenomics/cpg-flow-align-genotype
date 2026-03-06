"""
Runs Picard FixMateInformation on input BAMs to clean them prior to running Picard MarkDuplicates,
which is used in the workflow for both genome and exome samples.

This is required for some BAMs which may have missing or incorrect mate information that
causes MarkDuplicates to fail.
"""

import argparse

from cpg_utils import config, hail_batch
from hailtop.batch.job import Job


def run_fixmate(batch: hail_batch.Batch, bam_path: str, out_bam_path: str) -> Job:
    """Run Picard FixMateInformation on BAM file to clean mate information."""

    job = batch.new_job(
        'Picard FixMateInformation',
        attributes={'tool': 'picard'},
    )

    job.image(config.image_path('picard', '3.4.0-3'))
    job.memory('16Gi')
    job.storage('100Gi')
    job.cpu(4)

    # read in the BAM and index
    bam_localised = batch.read_input_group(
        bam=bam_path,
        bai=f'{bam_path}.bai',
    ).bam

    job.command(f"""\
    set -e

    picard FixMateInformation \\
        -I {bam_localised} \\
        -O {job.output_bam} \\
        --ADD_MATE_CIGAR true \\
        --IGNORE_MISSING_MATES false
    """)
    batch.write_output(job.output_bam, out_bam_path.removesuffix('.bam'))
    return job


def main(args):
    batch = hail_batch.get_batch()

    run_fixmate(batch, bam_path=args.bam_path, out_bam_path=args.out_bam_path)

    batch.run(wait=False)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--bam_path', type=str, required=True, help='Path to input BAM file')
    parser.add_argument('--out_bam_path', type=str, required=True, help='Path to output BAM file')
    args = parser.parse_args()
    main(args)
