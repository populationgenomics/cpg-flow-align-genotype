"""Strip /1 and /2 QNAME suffixes from a CRAM.

These suffixes, left by older sequencers, break mate pairing in samtools fastq
without collation.
"""

from hailtop.batch.job import Job

from cpg_utils import config, hail_batch, to_path


def run(
    batch: hail_batch.Batch,
    cram_path: str,
    sg_id: str,
    output_cram: str,
    job_attrs: dict,
) -> list[Job]:
    """Strip QNAME suffixes and write repaired CRAM to output_cram."""

    job = batch.new_job(
        'repair CRAM: strip QNAME suffixes',
        attributes=job_attrs | {'tool': 'samtools'},
    )

    job.image(config.config_retrieve(['images', 'samtools']))
    job.memory('standard')
    job.storage(f'{config.config_retrieve(["workflow", "genome_cram_gb"], "100")}Gi')

    cram_localised = batch.read_input_group(
        cram=cram_path,
        crai=f'{cram_path}.crai',
    ).cram

    job.declare_resource_group(
        output_cram={
            'cram': '{root}.cram',
            'cram.crai': '{root}.cram.crai',
        },
    )

    reference = hail_batch.fasta_res_group(batch)

    awk_strip = r"""awk 'BEGIN{FS=OFS="\t"} !/^@/{sub(/\/[12]$/,"",$1)} {print}'"""

    job.command(f"""\
    set -eo pipefail

    samtools view -h -T {reference.base} -@ 3 {cram_localised} | \
    {awk_strip} | \
    samtools view --write-index \
        -C -T {reference.base} -@ 3 \
        -o {job.output_cram.cram} -
    """)

    batch.write_output(job.output_cram, to_path(output_cram).with_suffix('').as_posix())
    return [job]