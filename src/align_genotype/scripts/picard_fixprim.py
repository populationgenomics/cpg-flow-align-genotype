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

from cpg_utils import config, hail_batch


def main(bam_in: str, bam_out: str) -> None:
    batch = hail_batch.get_batch()

    job = batch.new_job('Name Sort BAM')

    # use the _current_ image, hard coding would let this get out of sync
    job.image(config.config_retrieve(['workflow', 'driver_image']))
    job.memory('16Gi')
    job.storage('500Gi')

    cpu_count = 8
    job.cpu(cpu_count)

    bam_localised = batch.read_input_group(bam=bam_in, bai=f'{bam_in}.bai')

    job.declare_resource_group(output_bam={'bam': '{root}.bam', 'bam.bai': '{root}.bam.bai'})

    job.command(f"""
        set -ex
        samtools sort -n -@ {cpu_count} -o $BATCH_TMPDIR/name_sorted.bam {bam_localised.bam} -T $BATCH_TMPDIR
        python -m align_genotype.scripts.fix_primary_alignments --input $BATCH_TMPDIR/name_sorted.bam --output $BATCH_TMPDIR/fixed.bam
        rm $BATCH_TMPDIR/name_sorted.bam
        samtools sort -@ {cpu_count} -o {job.output_bam.bam} $BATCH_TMPDIR/name_sorted.bam -T $BATCH_TMPDIR
        rm $BATCH_TMPDIR/fixed.bam
        samtools index -@ {cpu_count} {job.output_bam.bam}
    """)  # noqa: E501

    batch.write_output(job.output_bam, bam_out.removesuffix('.bam'))

    batch.run(wait=False)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Fix duplicate primary alignments in BAM files using Hail Batch')
    parser.add_argument('--bam_path', type=str, required=True, help='Path to input BAM file')
    parser.add_argument('--out_bam_path', type=str, required=True, help='Path to output BAM file')
    args = parser.parse_args()
    main(bam_in=args.bam_path, bam_out=args.out_bam_path)
