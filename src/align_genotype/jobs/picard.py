"""
Create Hail Batch jobs to run Picard tools (marking duplicates, QC).
"""

import hailtop.batch as hb
from cpg_flow import resources
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
        storage_gb=config.config_retrieve(['workflow', 'picard_storage_gb'], 250),
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
