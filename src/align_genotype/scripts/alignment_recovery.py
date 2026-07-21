"""
The current alignment process can fail during merging, which leaves the workflow in limbo, having to
repeat the (expensive) alignment of every shard from scratch.

This script rescues those alignments: the per-shard name-sorted BAMs are usually still sitting in the
batch tmp bucket, so we can pick up where the pipeline failed - merge them, deduplicate, coordinate-sort,
and convert to a final CRAM - without re-running DRAGMAP.

The BAM shards MUST be provided in the same order the pipeline emitted them (Align_1, Align_2, ...), as
`samtools merge -n` relies on a consistent RG-ordered stream for dupblaster to catch cross-shard duplicates.

Example:
    python alignment_recovery.py \\
        -s CPG1234 \\
        -o gs://cpg-dataset-main/cram/CPG1234.cram \\
        -p seqr \\
        -i gs://tmp-bucket/c2dfbe/Align_1_6_-P167e/sorted_bam \\
           gs://tmp-bucket/c2dfbe/Align_2_6_-qPVi8/sorted_bam \\
           gs://tmp-bucket/c2dfbe/Align_3_6_-tBItI/sorted_bam \\
           gs://tmp-bucket/c2dfbe/Align_4_6_-Slx2L/sorted_bam \\
           gs://tmp-bucket/c2dfbe/Align_5_6_-M1pXr/sorted_bam \\
           gs://tmp-bucket/c2dfbe/Align_6_6_-51gSX/sorted_bam
"""

import argparse

from cpg_utils import config, hail_batch, to_path

from align_genotype.jobs.align import dedup_sort_cmd


def parser_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument('-s', '--sgid', type=str, required=True, help='Sequencing Group ID')
    parser.add_argument('-o', '--out', type=str, required=True, help='Path to final CRAM')
    parser.add_argument('-p', '--project', type=str, required=True, help='Project/dataset for registration')
    parser.add_argument(
        '-i',
        '--input',
        nargs='+',
        type=str,
        required=True,
        help='Paths to the per-shard name-sorted BAMs, in shard order (Align_1, Align_2, ...)',
    )
    parser.add_argument('--stage', type=str, default='AlignWithDragmap', help='Stage name, recorded in Metamist meta')
    parser.add_argument('--seq_type', type=str, default='genome', help='Sequencing Type')
    parser.add_argument('--storage', type=int, default=400, help='Storage in GiB to use for merge/sort')
    parser.add_argument(
        '--no-register',
        action='store_true',
        help='Skip creating a Metamist analysis entry (only produce the CRAM)',
    )
    return parser.parse_args()


def register_cram(output: str, sgid: str, dataset: str, stage: str, seq_type: str) -> None:
    """
    Create a completed `cram` analysis entry in Metamist. Runs inside the batch as a PythonJob that
    depends on the merge job, so registration only happens if the CRAM was produced successfully.
    """
    from cpg_flow.metamist import Metamist  # noqa: PLC0415

    Metamist().create_analysis(
        output=output,
        type_='cram',
        status='completed',
        sequencing_group_ids=[sgid],
        dataset=dataset,
        meta={'stage': stage, 'sequencing_type': seq_type, 'source': 'alignment_recovery'},
    )


if __name__ == '__main__':
    args = parser_args()

    batch_instance = hail_batch.get_batch(name=f'Merge alignments for {args.sgid}')

    # Read the BAM shards into the batch, preserving CLI order so `samtools merge -n` sees a
    # consistent RG-ordered stream. The list comprehension keeps ordering 1:1 with args.input.
    local_bams = [batch_instance.read_input(bam) for bam in args.input]

    fasta_reference = hail_batch.fasta_res_group(batch_instance)

    nthreads = config.config_retrieve(['workflow', 'merge_threads'], 16)

    merge_job = batch_instance.new_bash_job(name=f'Merge BAMs for {args.sgid}', attributes={'tool': 'samtools_merge'})
    merge_job.image(config.config_retrieve(['workflow', 'driver_image']))
    merge_job.storage(f'{args.storage}GiB')
    merge_job.cpu(nthreads)
    merge_job.memory('highmem')
    merge_job.declare_resource_group(
        output_cram={
            'cram': '{root}.cram',
            'cram.crai': '{root}.cram.crai',
        },
    )

    # Reproduce the pipeline's merge path exactly (see jobs/align.py):
    #   1. name-order merge of the shards -> RG-ordered stream
    #   2. dupblaster + coordinate-sort   -> catches cross-shard duplicates
    #   4. samtools view                  -> indexed CRAM (v3.0) written to GCS
    bam_string = ' '.join(map(str, local_bams))
    CMD = (
        f'samtools merge -n -@5 - {bam_string} '
        f'{dedup_sort_cmd(nthreads, merge_job.markdup_metrics)} '
        f'| samtools view --write-index -@5 '
        f'-T {fasta_reference.base} -O cram,version=3.0 -o {merge_job.output_cram.cram} -'
    )
    merge_job.command(CMD)

    # write_output takes the resource-group root (no suffix); .cram / .cram.crai are appended
    out_path = to_path(args.out)
    batch_instance.write_output(merge_job.output_cram, str(out_path.with_suffix('')))
    batch_instance.write_output(
        merge_job.markdup_metrics,
        f'{out_path.with_suffix("")}.markduplicates-metrics',
    )

    if not args.no_register:
        register_job = batch_instance.new_python_job(name=f'Register CRAM for {args.sgid}')
        register_job.image(config.config_retrieve(['workflow', 'driver_image']))
        register_job.depends_on(merge_job)
        register_job.call(register_cram, str(args.out), args.sgid, args.project, args.stage, args.seq_type)

    batch_instance.run(wait=False)
