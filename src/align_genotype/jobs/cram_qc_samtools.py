"""
Create Hail Batch jobs for samtools.
"""

from cpg_flow import resources
from cpg_utils import Path, config, hail_batch
from hailtop.batch.job import Job


def samtools_stats(
    cram_path: str,
    output: Path,
    job_attrs: dict,
) -> Job:
    """Run `samtools stats` for alignment QC."""

    batch_instance = hail_batch.get_batch()

    job = batch_instance.new_job(
        'samtools stats',
        attributes=job_attrs | {'tool': 'samtools'},
    )

    job.image(config.config_retrieve(['images', 'samtools']))
    res = resources.STANDARD.set_resources(job, fraction=1)
    reference = hail_batch.fasta_res_group(batch_instance)

    # read in the CRAM and index
    cram_localised = batch_instance.read_input_group(
        cram=cram_path,
        crai=f'{cram_path}.crai',
    ).cram

    job.command(f"""\
    samtools stats \\
        -@{res.get_nthreads() - 1} \\
        --reference {reference.base} \\
        {cram_localised} > {job.output_stats}
    """)
    batch_instance.write_output(job.output_stats, output)

    return job
