"""
Create Hail Batch jobs to run Picard tools (marking duplicates, QC).
"""

import hailtop.batch as hb
from cpg_flow import resources, utils
from cpg_utils import Path, config, hail_batch
from hailtop.batch.job import BashJob


def markdup(
    batch_instance: hb.Batch,
    sorted_bam: hb.Resource,
    job_attrs: dict,
    output_path: Path,
    out_markdup_metrics_path: Path,
) -> BashJob:
    """
    Make job that runs Picard MarkDuplicates and converts the result to CRAM.
    """
    job = batch_instance.new_bash_job(
        f'MarkDuplicates {"mito" if "mito" in str(output_path) else ""}',
        attributes=job_attrs | {'tool': 'picard_MarkDuplicates'},
    )

    job.image(config.config_retrieve(['images', 'picard']))

    # check for a memory override for impossible sequencing groups
    # if RAM is overridden, update the memory resource setting
    # check for a storage override for unreasonably large sequencing groups
    resource = resources.HIGHMEM.request_resources(
        ncpu=4,
        mem_gb=config.config_retrieve(['workflow', 'picard_mem_gb'], None),
        storage_gb=config.config_retrieve(['workflow', 'picard_storage_gb']),
    )

    resource.set_to_job(job)

    job.declare_resource_group(
        output_cram={
            'cram': '{root}.cram',
            'cram.crai': '{root}.cram.crai',
        },
    )

    fasta_reference = hail_batch.fasta_res_group(batch_instance)

    cmd = f"""
    picard {resource.java_mem_options()} MarkDuplicates \\
    I={sorted_bam} O={job.temp_bam} M={job.markdup_metrics} \\
    TMP_DIR=$(dirname {job.output_cram.cram})/picard-tmp \\
    ASSUME_SORT_ORDER=coordinate
    echo "MarkDuplicates finished successfully"

    rm {sorted_bam}

    samtools view --write-index -@{resource.get_nthreads() - 1} \\
    -T {fasta_reference.base} \\
    -O cram \\
    -o {job.output_cram.cram} \\
    {job.temp_bam}
    echo "samtools view finished successfully"
    """
    job.command(hail_batch.command(cmd, monitor_space=True))

    batch_instance.write_output(job.output_cram, str(output_path.with_suffix('')))
    batch_instance.write_output(
        job.markdup_metrics,
        out_markdup_metrics_path,
    )

    return job


def get_intervals(
    b: hb.Batch,
    scatter_count: int,
    job_attrs: dict[str, str] | None = None,
    output_prefix: Path | None = None,
) -> tuple[BashJob | None, list[hb.Resource]]:
    """
    Add a job that splits genome/exome intervals into sub-intervals to be used to
    parallelize variant calling.

    @param b: Hail Batch object,
    @param scatter_count: number of target sub-intervals,
    @param job_attrs: attributes for Hail Batch job,
    @param output_prefix: path optionally to save split subintervals.

    The job calls picard IntervalListTools to scatter the input interval list
    into scatter_count sub-interval lists, inspired by this WARP task :
    https://github.com/broadinstitute/warp/blob/bc90b0db0138747685b459c83ce52c8576ce03cd/tasks/broad/Utilities.wdl

    Note that we use the mode INTERVAL_SUBDIVISION instead of
    BALANCING_WITHOUT_INTERVAL_SUBDIVISION_WITH_OVERFLOW. Modes other than
    INTERVAL_SUBDIVISION produce an unpredictable number of intervals. WDL can
    handle that, but Hail Batch is not dynamic and expects a certain number
    of output files.
    """
    sequencing_type = config.config_retrieve(['workflow', 'sequencing_type'])
    source_intervals_path = config.config_retrieve(['references', f'{sequencing_type}_calling_interval_lists'])
    exclude_intervals_path = config.config_retrieve(['references', 'hg38_telomeres_and_centromeres_intervals_list'])

    if scatter_count == 1:
        # Special case when we don't need to split
        return None, [b.read_input(source_intervals_path)]

    if output_prefix:
        interval_lists_exist = all(
            utils.exists(output_prefix / f'{idx}.interval_list') for idx in range(1, scatter_count + 1)
        )
        if interval_lists_exist:
            return None, [b.read_input(str(output_prefix / f'{idx + 1}.interval_list')) for idx in range(scatter_count)]

    job = b.new_bash_job(
        f'Make {scatter_count} intervals for {sequencing_type}',
        attributes=(job_attrs or {}) | {'tool': 'picard IntervalListTools'},
    )
    job.image(config.config_retrieve(['images', 'picard']))
    resources.STANDARD.set_resources(job, storage_gb=16, mem_gb=2)

    break_bands_at_multiples_of = {
        'genome': 100000,
        'exome': 0,
    }.get(sequencing_type, 0)

    # If there are intervals to exclude, subtract them from the source intervals
    extra_cmd = f'-ACTION SUBTRACT -SI {b.read_input(str(exclude_intervals_path))}' if exclude_intervals_path else ''

    cmd = f"""
    set -o pipefail
    set -ex

    mkdir $BATCH_TMPDIR/out

    picard -Xms1000m -Xmx1500m \
    IntervalListTools \
    -SCATTER_COUNT {scatter_count} \
    -SUBDIVISION_MODE INTERVAL_SUBDIVISION \
    -UNIQUE true \
    -SORT true \
    -BREAK_BANDS_AT_MULTIPLES_OF {break_bands_at_multiples_of} \
    -I {b.read_input(source_intervals_path)} {extra_cmd} \
    -OUTPUT $BATCH_TMPDIR/out
    ls $BATCH_TMPDIR/out
    ls $BATCH_TMPDIR/out/*
    """

    for idx in range(scatter_count):
        name = f'temp_{str(idx + 1).zfill(4)}_of_{scatter_count}'
        cmd += f"""
        ln $BATCH_TMPDIR/out/{name}/scattered.interval_list {job[f'{idx + 1}.interval_list']}
        """

    job.command(cmd)

    if output_prefix:
        for idx in range(scatter_count):
            b.write_output(
                job[f'{idx + 1}.interval_list'],
                str(output_prefix / f'{idx + 1}.interval_list'),
            )

    intervals = [job[f'{idx + 1}.interval_list'] for idx in range(scatter_count)]

    return job, intervals
