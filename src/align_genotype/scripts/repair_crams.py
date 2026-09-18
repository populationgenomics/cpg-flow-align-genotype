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

from cpg_flow.status import complete_analysis_job
from cpg_utils import config, hail_batch
from metamist.graphql import gql, query


def strip_qname_suffixes(
    batch: hail_batch.Batch,
    cram_path: str,
    sg_id: str,
    job_attrs: dict,
) -> list[Job]:
    """Strip /1 and /2 QNAME suffixes from a CRAM, overwriting in place."""

    job = batch.new_job(
        'repair CRAM: strip QNAME suffixes',
        attributes=job_attrs | {'tool': 'samtools'},
    )

    job.image(config.config_retrieve(['images', 'samtools']))
    job.memory('16Gi')
    job.storage(f'{config.config_retrieve(["workflow", "genome_cram_gb"], "100")}Gi')
    job.cpu(4)

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

    batch.write_output(job.output_cram, cram_path.removesuffix('.cram'))
    return [job]


def trim_adapters(
    batch: hail_batch.Batch,
    cram_path: str,
    sg_id: str,
    job_attrs: dict,
) -> list[Job]:
    """Trim adapters and poly-G, then realign with BWA. Overwrites in place."""

    bwa_image = config.config_retrieve(['images', 'bwa'])
    fastp_image = config.config_retrieve(['images', 'fastp'])
    storage = f'{config.config_retrieve(["workflow", "genome_cram_gb"], "100")}Gi'

    reference = hail_batch.fasta_res_group(batch, indices=['amb', 'ann', 'bwt', 'pac', 'sa'])

    cram_localised = batch.read_input_group(
        cram=cram_path,
        crai=f'{cram_path}.crai',
    ).cram

    # Job 1: CRAM → interleaved FASTQ
    extract_fastq = batch.new_job(
        'repair CRAM: CRAM to FASTQ',
        attributes=job_attrs | {'tool': 'samtools'},
    )
    extract_fastq.image(bwa_image)
    extract_fastq.cpu(4)
    extract_fastq.memory('16Gi')
    extract_fastq.storage(storage)

    extract_fastq.command(f"""\
    set -eo pipefail

    samtools collate -u -O -T /tmp/collate_tmp \
        --reference {reference.base} {cram_localised} | \
    samtools fastq -n -@ 3 - > {extract_fastq.fastq}
    """)

    # Job 2: fastp adapter + poly-G trimming
    trim_reads = batch.new_job(
        'repair CRAM: fastp trim',
        attributes=job_attrs | {'tool': 'fastp'},
    )
    trim_reads.image(fastp_image)
    trim_reads.cpu(4)
    trim_reads.memory('16Gi')
    trim_reads.storage(storage)

    trim_reads.command(f"""\
    set -eo pipefail

    fastp --in1 {extract_fastq.fastq} --interleaved_in \
        --stdout \
        --detect_adapter_for_pe \
        --trim_poly_g \
        --thread 4 \
        --json /dev/null --html /dev/null \
        > {trim_reads.trimmed_fastq}
    """)

    # Job 3: BWA realign → sorted CRAM
    bwa_realign = batch.new_job(
        'repair CRAM: BWA realign',
        attributes=job_attrs | {'tool': 'bwa'},
    )
    bwa_realign.image(bwa_image)
    bwa_realign.cpu(8)
    bwa_realign.memory('highmem')
    bwa_realign.storage(storage)

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
        {reference.base} {trim_reads.trimmed_fastq} | \
    samtools view -C -T {reference.base} - | \
    samtools sort --write-index \
        -@ 4 \
        -o {bwa_realign.output_cram.cram}
    """)

    batch.write_output(bwa_realign.output_cram, cram_path.removesuffix('.cram'))
    return [extract_fastq, trim_reads, bwa_realign]


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
        analysis_type='cram-repair',
        cohort_ids=[],
        sg_ids=[sg_id],
        project_name=dataset,
        meta={'repair_type': repair_type},
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
            if analyses and (path := analyses[0]['outputs'].get('path')):
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
        batch = hail_batch.get_batch()
        for sg_id, cram_path in sg_crams:
            repair_jobs = repair_fn(batch, cram_path, sg_id, job_attrs={'repair_type': args.repair_type})
            _register_repair(batch, cram_path, sg_id, args.repair_type, args.dataset, depends_on=repair_jobs)
        batch.run(wait=False)