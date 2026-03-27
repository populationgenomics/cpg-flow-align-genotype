"""
Pipeline script to fix duplicate primary alignments in BAM files using Hail Batch.

This script orchestrates a 3-step pipeline:
1. Name sort the BAM file
2. Fix duplicate primary alignments (using fix_primary_alignments.py)
3. Coordinate sort and index the output

Usage:
    python picard_fixprim.py --bam_path <input.bam> --out_bam_path <output.bam>
"""

import argparse

from cpg_utils import hail_batch
from hailtop.batch.job import Job

DOCKER_IMAGE = 'australia-southeast1-docker.pkg.dev/cpg-common/images-dev/cpg-flow-align-genotype:0.4.6-1'


def job_name_sort(batch: hail_batch.Batch, bam_path: str) -> Job:
    """Job 1: Sort BAM by name to group duplicate alignments together."""
    job = batch.new_job('Name Sort BAM')

    job.image(DOCKER_IMAGE)
    job.memory('16Gi')
    job.storage('500Gi')

    cpu_count = 8
    job.cpu(cpu_count)

    bam_localised = batch.read_input_group(bam=bam_path, bai=f'{bam_path}.bai')

    job.declare_resource_group(name_sorted={'bam': '{root}.bam'})

    job.command(f"""
        set -ex
        samtools sort -n -@ {cpu_count} -o {job.name_sorted.bam} {bam_localised.bam} -T $BATCH_TMPDIR
    """)

    return job


def job_fix_primary_alignments(batch: hail_batch.Batch, name_sorted_bam) -> Job:
    """Job 2: Process name-sorted BAM to fix duplicate primary alignments."""
    job = batch.new_job('Fix Primary Alignments')

    job.image(DOCKER_IMAGE)
    job.memory('16Gi')
    job.storage('500Gi')
    job.cpu(4)

    job.declare_resource_group(fixed={'bam': '{root}.bam'})

    job.command(f"""
        set -ex
        python -m align_genotype.scripts.fix_primary_alignments --input {name_sorted_bam} --output {job.fixed.bam}
    """)

    return job


def job_coordinate_sort_and_index(batch: hail_batch.Batch, fixed_bam, out_bam_path: str) -> Job:
    """Job 3: Sort BAM by coordinate and create index."""
    job = batch.new_job('Coordinate Sort and Index')

    job.image(DOCKER_IMAGE)
    job.memory('16Gi')
    job.storage('500Gi')

    cpu_count = 8
    job.cpu(cpu_count)

    job.declare_resource_group(output_bam={'bam': '{root}.bam', 'bam.bai': '{root}.bam.bai'})

    job.command(f"""
        set -ex
        samtools sort -@ {cpu_count} -o {job.output_bam.bam} {fixed_bam} -T $BATCH_TMPDIR
        samtools index {job.output_bam.bam}
    """)

    batch.write_output(job.output_bam, out_bam_path.removesuffix('.bam'))
    return job


def main(args):
    batch = hail_batch.get_batch()

    # Job 1: Name sort
    job1 = job_name_sort(batch, bam_path=args.bam_path)

    # Job 2: Fix primary alignments (depends on job1)
    job2 = job_fix_primary_alignments(batch, name_sorted_bam=job1.name_sorted.bam)

    # Job 3: Coordinate sort and index (depends on job2)
    job_coordinate_sort_and_index(batch, fixed_bam=job2.fixed.bam, out_bam_path=args.out_bam_path)

    batch.run(wait=False)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Fix duplicate primary alignments in BAM files using Hail Batch')
    parser.add_argument('--bam_path', type=str, required=True, help='Path to input BAM file')
    parser.add_argument('--out_bam_path', type=str, required=True, help='Path to output BAM file')
    args = parser.parse_args()
    main(args)
