"""
A wrapper script that queues Hail batch jobs to repair CRAM files with known
defects that break downstream pipelines (e.g. GATK-SV, Mito, etc.).

Supported repair types:
  strip-qnames   — Remove /1 /2 QNAME suffixes left by older sequencers.
                    These break mate pairing in samtools fastq without collation.
  trim-adapters  — Trim Illumina adapters and poly-G artefacts with fastp,
                    then realign with BWA. For samples where adapters were not
                    stripped before alignment, causing spurious soft-clipping.

Queries metamist for CRAM paths by --sg-ids.
"""

import argparse

from hailtop.batch.job import Job
from loguru import logger

from cpg_flow.status import complete_analysis_job
from cpg_flow.utils import exists
from cpg_utils import config, hail_batch, to_path
from metamist.graphql import gql, query


def strip_qname_suffixes(
    batch: hail_batch.Batch,
    cram_path: str,
    sg_id: str,
    job_attrs: dict,
    staging_prefix: str,
) -> list[Job]:
    """Strip /1 and /2 QNAME suffixes from a CRAM, overwriting in place."""
    staging_cram = to_path(staging_prefix) / sg_id / 'stripped.cram'

    if exists(staging_cram):
        logger.info(f'Skipping strip QNAME suffixes for {sg_id}: output exists at {staging_cram}')
        return []

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

    batch.write_output(job.output_cram, str(staging_cram).removesuffix('.cram'))
    batch.write_output(job.output_cram, cram_path.removesuffix('.cram'))
    return [job]


def trim_adapters(
    batch: hail_batch.Batch,
    cram_path: str,
    sg_id: str,
    job_attrs: dict,
    staging_prefix: str,
    fastq_path: str | None = None,
) -> list[Job]:
    """Trim adapters and poly-G, then realign with BWA. Overwrites in place."""

    staging = to_path(staging_prefix) / sg_id

    bwa_image = config.config_retrieve(['images', 'bwa'])
    fastp_image = config.config_retrieve(['images', 'fastp'])
    storage = f'{config.config_retrieve(["workflow", "genome_cram_gb"], "400")}Gi'

    reference = hail_batch.fasta_res_group(batch, indices=['amb', 'ann', 'bwt', 'pac', 'sa'])

    jobs: list[Job] = []

    # Job 1: CRAM → interleaved FASTQ
    if fastq_path:
        fastq_input = batch.read_input(fastq_path)
    else:
        fastq_out = staging / 'interleaved.fastq.gz'
        if exists(fastq_out):
            logger.info(f'Skipping FASTQ extraction for {sg_id}: output exists at {fastq_out}')
            fastq_input = batch.read_input(str(fastq_out))
        else:
            cram_localised = batch.read_input_group(
                cram=cram_path,
                crai=f'{cram_path}.crai',
            ).cram

            extract_fastq = batch.new_job('repair CRAM: CRAM to FASTQ', attributes=job_attrs | {'tool': 'samtools'})
            extract_fastq.image(bwa_image)
            extract_fastq.memory('32Gi')
            extract_fastq.storage('1000Gi')

            extract_fastq.command(f"""\
            set -eo pipefail

            samtools collate -u -O \
                --reference {reference.base} {cram_localised} $BATCH_TMPDIR/collate_tmp | \
            samtools fastq -n -@ 3 - | \
            gzip > {extract_fastq.fastq_gz}
            """)
            batch.write_output(extract_fastq.fastq_gz, str(fastq_out))
            fastq_input = extract_fastq.fastq_gz
            jobs.append(extract_fastq)

    # Job 2: fastp adapter + poly-G trimming
    trimmed_out = staging / 'trimmed.fastq'
    if exists(trimmed_out):
        logger.info(f'Skipping fastp trim for {sg_id}: output exists at {trimmed_out}')
        trimmed_input = batch.read_input(str(trimmed_out))
    else:
        trim_reads = batch.new_job(
            'repair CRAM: fastp trim',
            attributes=job_attrs | {'tool': 'fastp'},
        )
        trim_reads.image(fastp_image)
        trim_reads.memory('standard')
        trim_reads.storage(storage)

        trim_reads.command(f"""\
        set -eo pipefail

        gunzip -c {fastq_input} | \
        fastp --stdin --interleaved_in \
            --stdout \
            --detect_adapter_for_pe \
            --trim_poly_g \
            --thread 4 \
            --json /dev/null --html /dev/null \
            > {trim_reads.trimmed_fastq}
        """)
        batch.write_output(trim_reads.trimmed_fastq, str(trimmed_out))
        trimmed_input = trim_reads.trimmed_fastq
        jobs.append(trim_reads)

    # Job 3: BWA realign → sorted CRAM
    staging_cram = staging / 'realigned.cram'
    if exists(staging_cram):
        logger.info(f'Skipping BWA realign for {sg_id}: output exists at {staging_cram}')
    else:
        bwa_realign = batch.new_job(
            'repair CRAM: BWA realign',
            attributes=job_attrs | {'tool': 'bwa'},
        )
        bwa_realign.image(bwa_image)
        bwa_realign.memory('highmem')
        bwa_realign.storage(storage)
        bwa_realign.spot(False)

        bwa_realign.declare_resource_group(
            output_cram={
                'cram': '{root}.cram',
                'cram.crai': '{root}.cram.crai',
            },
        )

        bwa_realign.command(f"""\
        set -eo pipefail

        bwa mem -K 100000000 -p -v 3 -t 8 -Y \
            -R '@RG\\tID:{sg_id}\\tLB:LB0\\tPL:PL0\\tPU:PU0\\tSM:{sg_id}' \
            {reference.base} {trimmed_input} | \
        samtools view -C -T {reference.base} - | \
        samtools sort --write-index \
            -@ 4 \
            -o {bwa_realign.output_cram.cram}
        """)

        batch.write_output(bwa_realign.output_cram, str(staging_cram).removesuffix('.cram'))
        batch.write_output(bwa_realign.output_cram, cram_path.removesuffix('.cram'))
        jobs.append(bwa_realign)

    return jobs


def _register_repair(
    batch: hail_batch.Batch,
    cram_path: str,
    sg_id: str,
    repair_type: str,
    dataset: str,
    depends_on: list[Job],
) -> Job:
    """Register a cram-repair analysis entry in metamist after the repair job completes."""

    reg_job = batch.new_python_job(
        f'register_repair_{sg_id}',
        attributes={'tool': 'metamist'},
    )
    reg_job.image(config.config_retrieve(['workflow', 'driver_image']))
    reg_job.call(
        complete_analysis_job,
        output=cram_path,
        analysis_type='custom',
        cohort_ids=[],
        sg_ids=[sg_id],
        project_name=dataset,
        meta={'stage': 'cram-repair', 'repair_type': repair_type},
    )
    for j in depends_on:
        reg_job.depends_on(j)
    return reg_job


REPAIR_FUNCTIONS = {
    'strip-qnames': strip_qname_suffixes,
    'trim-adapters': trim_adapters,
}


def get_cram_paths_for_sgs(sg_ids: list[str]) -> list[tuple[str, str]]:
    """Query metamist for CRAM paths by sequencing group IDs. Returns (sg_id, cram_path) tuples."""
    query_str = gql(
        """
        query GetCramPaths($sg_id: String!) {
            sequencingGroups(id: {eq: $sg_id}) {
                id
                analyses(type: {eq: "cram"}, status: {eq: COMPLETED}) {
                    timestampCompleted
                    outputs
                }
            }
        }
        """
    )
    results = []
    for sg_id in sg_ids:
        result = query(query_str, variables={'sg_id': sg_id})
        for sg in result['sequencingGroups']:
            analyses = sorted(
                sg['analyses'],
                key=lambda a: a['timestampCompleted'] or '',
                reverse=True,
            )
            if analyses:
                outputs = analyses[0]['outputs']
                path = outputs.get('path') if isinstance(outputs, dict) else outputs
                if path:
                    results.append((sg['id'], path))
    return results


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Repair CRAM files with known defects.')
    parser.add_argument(
        '--repair-type',
        choices=list(REPAIR_FUNCTIONS),
        required=True,
        help='Type of repair to apply.',
    )
    parser.add_argument(
        '--dataset',
        required=True,
        help='Metamist dataset/project name for analysis registration.',
    )
    parser.add_argument(
        '--sg-ids',
        nargs='+',
        required=True,
        help='Sequencing group IDs to repair (queries metamist for CRAM paths).',
    )
    parser.add_argument(
        '--staging-prefix',
        help='GCS prefix for intermediate outputs (enables job skipping on re-runs). Required for non-dry runs.',
    )
    parser.add_argument(
        '--fastq-path',
        help='Skip CRAM-to-FASTQ extraction and use this existing GCS FASTQ path instead.',
    )
    parser.add_argument('--dry-run', action='store_true', help='Print the CRAM files that would be repaired.')
    args = parser.parse_args()

    sg_crams = get_cram_paths_for_sgs(args.sg_ids)

    if not sg_crams:
        print('No CRAM files found.')
        raise SystemExit(1)

    repair_fn = REPAIR_FUNCTIONS[args.repair_type]

    if args.dry_run:
        for sg_id, cram_path in sg_crams:
            print(f'[{args.repair_type}] {sg_id} {cram_path}')
    else:
        if not args.staging_prefix:
            parser.error('--staging-prefix is required for non-dry runs')

        batch = hail_batch.get_batch()
        for sg_id, cram_path in sg_crams:
            kwargs: dict = {
                'job_attrs': {'repair_type': args.repair_type},
                'staging_prefix': args.staging_prefix,
            }
            if args.fastq_path and repair_fn == trim_adapters:
                kwargs['fastq_path'] = args.fastq_path
            repair_jobs = repair_fn(batch, cram_path, sg_id, **kwargs)
            if repair_jobs:
                _register_repair(batch, cram_path, sg_id, args.repair_type, args.dataset, depends_on=repair_jobs)
            else:
                logger.info(f'All repair jobs for {sg_id} already completed, skipping registration')
        batch.run(wait=False)
