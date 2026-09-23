"""Trim Illumina adapters and poly-G artefacts with fastp, then realign with DRAGMAP.

For samples where adapters were not stripped before alignment, causing spurious
soft-clipping. The alignment pipeline matches the production DRAGMAP workflow in
align_genotype/jobs/align.py: DRAGMAP → dupblaster → coordinate sort → CRAM v3.0.
"""

import os.path

from loguru import logger

from hailtop.batch.job import Job

from cpg_flow.utils import exists
from cpg_utils import config, hail_batch, to_path

DRAGMAP_INDEX_FILES = ['hash_table.cfg.bin', 'hash_table.cmp', 'reference.bin']


def run(  # noqa: PLR0915
    batch: hail_batch.Batch,
    cram_path: str,
    sg_id: str,
    output_cram: str,
    job_attrs: dict,
    skip_jobs: set[str] | None = None,
) -> list[Job]:
    """Trim adapters/poly-G and realign. Writes repaired CRAM to output_cram."""

    _skip = skip_jobs or set()
    output_dir = to_path(output_cram).parent

    fastp_image = config.config_retrieve(['images', 'fastp'])
    dragmap_image = config.config_retrieve(['images', 'dragmap'])
    storage = f'{config.config_retrieve(["workflow", "genome_cram_gb"], "400")}Gi'

    nthreads = 16
    sort_threads = min(nthreads, 6) - 1

    reference = hail_batch.fasta_res_group(batch)
    dragmap_index = batch.read_input_group(
        **{
            k.replace('.', '_'): os.path.join(config.config_retrieve(['references', 'dragmap_prefix']), k)
            for k in DRAGMAP_INDEX_FILES
        },
    )

    jobs: list[Job] = []

    # --- Job 1: CRAM → interleaved FASTQ ---
    fastq_out = output_dir / f'{sg_id}_interleaved.fastq.gz'
    if 'extract' in _skip or exists(fastq_out):
        logger.info(f'Skipping FASTQ extraction for {sg_id}: output exists at {fastq_out}')
        fastq_input = batch.read_input(str(fastq_out))
    else:
        cram_localised = batch.read_input_group(
            cram=cram_path,
            crai=f'{cram_path}.crai',
        ).cram

        extract_fastq = batch.new_job('repair CRAM: CRAM to FASTQ', attributes=job_attrs | {'tool': 'samtools'})
        extract_fastq.image(dragmap_image)
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

    # --- Job 2: fastp adapter + poly-G trimming ---
    trimmed_out = output_dir / f'{sg_id}_trimmed.fastq.gz'
    if 'trim' in _skip or exists(trimmed_out):
        logger.info(f'Skipping fastp trim for {sg_id}: output exists at {trimmed_out}')
        trimmed_input = batch.read_input(str(trimmed_out))
    else:
        trim_reads = batch.new_job(
            'repair CRAM: fastp trim',
            attributes=job_attrs | {'tool': 'fastp'},
        )
        trim_reads.image(fastp_image)
        trim_reads.cpu(8)
        trim_reads.memory('16Gi')
        trim_reads.storage(storage)

        trim_reads.command(f"""\
        set -eo pipefail

        pigz -dc {fastq_input} | \
        fastp --stdin --interleaved_in \
            --stdout \
            --detect_adapter_for_pe \
            --trim_poly_g \
            --thread 4 \
            --json /dev/null --html /dev/null | \
        pigz -p 4 > {trim_reads.trimmed_fastq}
        """)
        batch.write_output(trim_reads.trimmed_fastq, str(trimmed_out))
        trimmed_input = trim_reads.trimmed_fastq
        jobs.append(trim_reads)

    # --- Job 3: DRAGMAP align → dupblaster dedup → coordinate sort → CRAM v3.0 ---
    if 'realign' in _skip or exists(output_cram):
        logger.info(f'Skipping DRAGMAP realign for {sg_id}: output exists at {output_cram}')
    else:
        align_job = batch.new_job(
            'repair CRAM: DRAGMAP realign',
            attributes=job_attrs | {'tool': 'dragmap'},
        )
        align_job.image(dragmap_image)
        align_job.cpu(nthreads)
        align_job.memory('highmem')
        align_job.storage(storage)
        align_job.spot(False)

        align_job.declare_resource_group(
            output_cram={
                'cram': '{root}.cram',
                'cram.crai': '{root}.cram.crai',
            },
        )

        storage_buffer = config.config_retrieve(['workflow', 'align_buffer_kb'], 2097152)

        align_job.command(f"""\
        set -eo pipefail

        watch_disk() {{
        local min_free_kb={storage_buffer}
        while true; do
          local avail
          avail=$(df --output=avail "$BATCH_TMPDIR" | tail -1)
          if (( avail < min_free_kb )); then
            echo "FATAL: $BATCH_TMPDIR has ${{avail}}KB free (< ${{min_free_kb}}KB) — aborting before disk fills" >&2
            kill -TERM -$$ 2>/dev/null
            exit 1
          fi
          sleep 15
        done
        }}
        watch_disk &
        WATCHDOG_PID=$!
        trap 'kill "$WATCHDOG_PID" 2>/dev/null' EXIT

        mkfifo r1
        pigz -dc {trimmed_input} > r1 &
        pid_r1=$!

        dragen-os -r {dragmap_index} --interleaved=1 -b r1 \
            --RGID {sg_id} --RGSM {sg_id} \
            --num-threads {nthreads - 1} \
        | dupblaster --stats {align_job.markdup_metrics} \
        | samtools sort -@{sort_threads} -T $BATCH_TMPDIR/samtools-dd-tmp -Obam \
        | samtools view --write-index -@{sort_threads} \
            -T {reference.base} -O cram,version=3.0 \
            -o {align_job.output_cram.cram} -

        if wait $pid_r1; then
            echo "Background decompression finished successfully"
        else
            echo "Background decompression failed" >&2
            exit 1
        fi
        """)

        batch.write_output(align_job.output_cram, to_path(output_cram).with_suffix('').as_posix())
        jobs.append(align_job)

    return jobs
