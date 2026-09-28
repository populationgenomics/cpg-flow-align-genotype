"""Trim Illumina adapters and poly-G artefacts with fastp, then realign with DRAGMAP.

For samples where adapters were not stripped before alignment, leaving spurious
soft-clipping. The alignment matches the production DRAGMAP workflow in
align_genotype/jobs/align.py: DRAGMAP -> dupblaster -> coordinate sort -> CRAM v3.0.

Archives the original CRAM before realigning, and optionally registers the result as a
new `cram` analysis. Retiring the analyses derived from the old CRAM is a separate step
- see sg_reset.py.

Intermediates (interleaved FASTQ, trimmed FASTQ) are written to declared GCS paths
beside the output CRAM, so a re-run picks up whatever a previous run completed.

Standalone by design: no shared imports, so it can be lifted into another workflow
engine as a single process with an input path, an output path and a sample ID.
"""

import argparse
import os.path

from loguru import logger

from hailtop.batch.job import Job

from cpg_flow.utils import exists
from cpg_utils import config, hail_batch, to_path

DRAGMAP_INDEX_FILES = ['hash_table.cfg.bin', 'hash_table.cmp', 'reference.bin']
BACKUP_DIR = 'bad_cram'


def archive_original(batch: hail_batch.Batch, cram_path: str, job_attrs: dict) -> Job | None:
    """Copy the original CRAM and index to bad_cram/ before it is superseded.

    Server-side GCS copy, so the data never transits the worker and the job needs no
    attached storage. Returns None if the archive already exists.
    """

    dest = config.dataset_path(f'{BACKUP_DIR}/{to_path(cram_path).name}')
    dest_dir = str(to_path(dest).parent)

    if exists(dest) and exists(f'{dest}.crai'):
        logger.info(f'Skipping CRAM archive: already present at {dest}')
        return None

    job = batch.new_bash_job('repair CRAM: archive original', attributes=job_attrs | {'tool': 'gcloud'})
    job.image(config.config_retrieve(['workflow', 'driver_image']))
    hail_batch.authenticate_cloud_credentials_in_job(job)

    job.command(f"""\
    set -eo pipefail

    gcloud storage cp {cram_path} {cram_path}.crai {dest_dir}/

    # fail loudly rather than let a repair proceed on an unverified archive
    for f in {dest} {dest}.crai; do
        gcloud storage objects describe "$f" --format='value(size)' > /dev/null \
            || {{ echo "FATAL: archive missing after copy: $f" >&2; exit 1; }}
    done
    echo "archived original CRAM to {dest_dir}/"
    """)

    return job


def trim_and_realign(
    batch: hail_batch.Batch,
    cram_path: str,
    sg_id: str,
    output_cram: str,
    job_attrs: dict,
) -> list[Job]:
    """CRAM -> interleaved FASTQ -> fastp -> DRAGMAP -> CRAM. Skips whatever already exists."""

    output_dir = to_path(output_cram).parent
    storage = f'{config.config_retrieve(["workflow", "genome_cram_gb"], "400")}Gi'
    dragmap_image = config.config_retrieve(['images', 'dragmap'])

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

    # --- CRAM to interleaved FASTQ ---
    fastq_out = output_dir / f'{sg_id}_interleaved.fastq.gz'
    if exists(fastq_out):
        logger.info(f'Reusing extracted FASTQ: {fastq_out}')
        fastq_input = batch.read_input(str(fastq_out))
    else:
        cram_localised = batch.read_input_group(cram=cram_path, crai=f'{cram_path}.crai').cram

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

    # --- fastp adapter and poly-G trimming ---
    trimmed_out = output_dir / f'{sg_id}_trimmed.fq.gz'
    if exists(trimmed_out):
        logger.info(f'Reusing trimmed FASTQ: {trimmed_out}')
        trimmed_input = batch.read_input(str(trimmed_out))
    else:
        trim_reads = batch.new_job('repair CRAM: fastp trim', attributes=job_attrs | {'tool': 'fastp'})
        trim_reads.image(config.config_retrieve(['images', 'fastp']))
        trim_reads.cpu(8)
        trim_reads.memory('16Gi')
        trim_reads.storage(storage)

        trim_reads.command(f"""\
        set -eo pipefail

        fastp --in1 {fastq_input} --interleaved_in \
            --stdout \
            --detect_adapter_for_pe \
            --trim_poly_g \
            --thread 4 \
            --json /dev/null --html /dev/null | \
        gzip > {trim_reads.trimmed_fastq}
        """)
        batch.write_output(trim_reads.trimmed_fastq, str(trimmed_out))
        trimmed_input = trim_reads.trimmed_fastq
        jobs.append(trim_reads)

    # --- DRAGMAP align -> dupblaster dedup -> coordinate sort -> CRAM v3.0 ---
    align_job = batch.new_job('repair CRAM: DRAGMAP realign', attributes=job_attrs | {'tool': 'dragmap'})
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

    # dragen-os reads gzipped FASTQ via -1, but a batch-tmp input path has no extension
    # for it to detect, hence the symlink.
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

    ln -s {trimmed_input} $BATCH_TMPDIR/trimmed_interleaved.fq.gz

    dragen-os -r {dragmap_index} --interleaved=1 \
        -1 $BATCH_TMPDIR/trimmed_interleaved.fq.gz \
        --RGID {sg_id} --RGSM {sg_id} \
        --num-threads {nthreads - 1} \
    | dupblaster --stats {align_job.markdup_metrics} \
    | samtools sort -@{sort_threads} -T $BATCH_TMPDIR/samtools-dd-tmp -Obam \
    | samtools view --write-index -@{sort_threads} \
        -T {reference.base} -O cram,version=3.0 \
        -o {align_job.output_cram.cram} -
    """)

    batch.write_output(align_job.output_cram, str(to_path(output_cram).with_suffix('')))
    jobs.append(align_job)

    return jobs


def register_cram(output: str, sg_id: str, dataset: str) -> None:
    """Create a completed `cram` analysis for the repaired CRAM.

    Runs inside the batch as a PythonJob depending on the alignment job, so registration
    only happens if the CRAM was actually produced. Does not inactivate anything.
    """
    from cpg_flow.metamist import Metamist  # noqa: PLC0415

    Metamist().create_analysis(
        output=output,
        type_='cram',
        status='completed',
        sequencing_group_ids=[sg_id],
        dataset=dataset,
        meta={'source': 'cram-repair', 'repair_type': 'trim-adapters'},
    )


def main() -> None:
    parser = argparse.ArgumentParser(description='Trim adapters and poly-G, then realign with DRAGMAP.')
    parser.add_argument('--cram-path', required=True, help='Input CRAM (GCS path). Needs a .crai beside it.')
    parser.add_argument('--output-path', required=True, help='Output CRAM (GCS path).')
    parser.add_argument('--sg-id', required=True, help='Sequencing group ID, used for the read group and metamist.')
    parser.add_argument('--dataset', help='Metamist dataset for registration. Required with --register.')
    parser.add_argument(
        '--register',
        action='store_true',
        help='Register the repaired CRAM as a new `cram` analysis. Off by default.',
    )
    parser.add_argument('--dry-run', action='store_true', help='Print what would be done without submitting.')
    args = parser.parse_args()

    if args.register and not args.dataset:
        parser.error('--register requires --dataset')

    output_dir = to_path(args.output_path).parent

    if args.dry_run:
        print(f'[trim-adapters] {args.sg_id}  {args.cram_path} -> {args.output_path}')
        print(f'  archive -> {config.dataset_path(f"{BACKUP_DIR}/{to_path(args.cram_path).name}")}')
        for label, path in (
            ('interleaved FASTQ', output_dir / f'{args.sg_id}_interleaved.fastq.gz'),
            ('trimmed FASTQ', output_dir / f'{args.sg_id}_trimmed.fq.gz'),
            ('output CRAM', to_path(args.output_path)),
        ):
            print(f'  {label}: {path} [{"exists, will reuse" if exists(path) else "will create"}]')
        print(f'  register: {args.register}')
        return

    if exists(args.output_path):
        logger.info(f'Output already exists, nothing to do: {args.output_path}')
        return

    batch = hail_batch.get_batch()
    job_attrs = {'repair_type': 'trim-adapters', 'sequencing_group': args.sg_id}

    archive_job = archive_original(batch, args.cram_path, job_attrs)
    repair_jobs = trim_and_realign(batch, args.cram_path, args.sg_id, args.output_path, job_attrs)

    # nothing references the archive job's outputs, so the ordering must be explicit.
    # applied to every job, since which one runs first depends on what already exists.
    if archive_job:
        for job in repair_jobs:
            job.depends_on(archive_job)

    if args.register:
        reg_job = batch.new_python_job(f'register repaired CRAM {args.sg_id}', attributes={'tool': 'metamist'})
        reg_job.image(config.config_retrieve(['workflow', 'driver_image']))
        reg_job.depends_on(repair_jobs[-1])
        reg_job.call(register_cram, args.output_path, args.sg_id, args.dataset)

    batch.run(wait=False)


if __name__ == '__main__':
    main()
