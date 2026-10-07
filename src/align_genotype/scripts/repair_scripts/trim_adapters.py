"""Trim Illumina adapters and poly-G artefacts with fastp, then realign with DRAGMAP.

For samples where adapters were not stripped before alignment, leaving spurious
soft-clipping. The alignment matches the production DRAGMAP workflow in
align_genotype/jobs/align.py: DRAGMAP -> dupblaster -> coordinate sort -> CRAM v3.0.

The realigned CRAM replaces the original in place, so the archived copy is the only
remaining copy of the original - the archive job verifies it by CRC32C before the
realignment is allowed to run.

The archive location is noted on that CRAM's existing `cram` analysis. No new analysis
is created. Retiring the analyses derived from the old CRAM is a separate step - see
sg_reset.py.

Intermediates (interleaved FASTQ, trimmed FASTQ) are written under cram_repair/ so a
re-run picks up whatever a previous run completed - extraction and trimming together
are most of the runtime.

Standalone by design: no shared imports, so it can be lifted into another workflow
engine as a single process with a CRAM path and a sample ID.
"""

import argparse
import os.path

from loguru import logger

from hailtop.batch.job import Job

from cpg_flow.utils import exists
from cpg_utils import config, hail_batch, to_path

DRAGMAP_INDEX_FILES = ['hash_table.cfg.bin', 'hash_table.cmp', 'reference.bin']
BACKUP_DIR = 'bad_cram'
INTERMEDIATE_DIR = 'cram_repair'


def archive_original(batch: hail_batch.Batch, cram_path: str, job_attrs: dict) -> Job | None:
    """Copy the original CRAM and index to bad_cram/ before it is overwritten in place.

    Server-side GCS copy, so the data never transits the worker and the job needs no
    attached storage. Compares CRC32C before and after: the realignment overwrites the
    original, so an unverified archive would mean the only copy is unaccounted for.
    Returns None if the archive already exists.
    """

    dest = config.dataset_path(f'{BACKUP_DIR}/{to_path(cram_path).name}')
    dest_dir = str(to_path(dest).parent)

    # An existing archive means a repair has run before, so cram_path may already hold a
    # repaired CRAM. Copying again would overwrite the preserved original with it.
    if exists(dest):
        if exists(f'{dest}.crai'):
            logger.warning(f'Archive already present at {dest} - this CRAM has been repaired before')
            return None
        msg = (
            f'Archive at {dest} exists but its index does not. Re-archiving would overwrite the '
            f'preserved original with whatever is at {cram_path} now, which may already be repaired. '
            f'Reindex or remove the partial archive by hand before re-running.'
        )
        raise RuntimeError(msg)

    job = batch.new_bash_job('repair CRAM: archive original', attributes=job_attrs | {'tool': 'gcloud'})
    job.image(config.config_retrieve(['workflow', 'driver_image']))
    hail_batch.authenticate_cloud_credentials_in_job(job)

    job.command(f"""\
    set -eo pipefail

    crc() {{ gcloud storage objects describe "$1" --format="value(crc32c_hash)"; }}

    before_cram=$(crc {cram_path})
    before_crai=$(crc {cram_path}.crai)

    # --no-clobber so the preserved original can never be overwritten, even if the
    # existence check above raced or was bypassed
    gcloud storage cp --no-clobber {cram_path} {cram_path}.crai {dest_dir}/

    after_cram=$(crc {dest})
    after_crai=$(crc {dest}.crai)

    # the realignment overwrites the original, so refuse to continue on an unverified archive
    if [[ -z "$before_cram" || "$before_cram" != "$after_cram" ]]; then
        echo "FATAL: CRAM archive checksum mismatch ($before_cram vs $after_cram)" >&2
        exit 1
    fi
    if [[ -z "$before_crai" || "$before_crai" != "$after_crai" ]]; then
        echo "FATAL: index archive checksum mismatch ($before_crai vs $after_crai)" >&2
        exit 1
    fi

    echo "archived original CRAM to {dest_dir}/ (crc32c $after_cram)"
    """)

    return job


def trim_and_realign(
    batch: hail_batch.Batch,
    cram_path: str,
    sg_id: str,
    job_attrs: dict,
) -> list[Job]:
    """CRAM -> interleaved FASTQ -> fastp -> DRAGMAP -> CRAM, in place. Reuses what exists."""

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
    fastq_out = config.dataset_path(f'{INTERMEDIATE_DIR}/{sg_id}_interleaved.fastq.gz')
    if exists(fastq_out):
        logger.info(f'Reusing extracted FASTQ: {fastq_out}')
        fastq_input = batch.read_input(fastq_out)
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
        batch.write_output(extract_fastq.fastq_gz, fastq_out)
        fastq_input = extract_fastq.fastq_gz
        jobs.append(extract_fastq)

    # --- fastp adapter and poly-G trimming ---
    trimmed_out = config.dataset_path(f'{INTERMEDIATE_DIR}/{sg_id}_trimmed.fq.gz')
    if exists(trimmed_out):
        logger.info(f'Reusing trimmed FASTQ: {trimmed_out}')
        trimmed_input = batch.read_input(trimmed_out)
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
        batch.write_output(trim_reads.trimmed_fastq, trimmed_out)
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
    # --metrics-prefix is mandatory, but we don't need the result
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
    | dupblaster --metrics-prefix markdup_metrics \
    | samtools sort -@{sort_threads} -T $BATCH_TMPDIR/samtools-dd-tmp -Obam \
    | samtools view --write-index -@{sort_threads} \
        -T {reference.base} -O cram,version=3.0 \
        -o {align_job.output_cram.cram} -
    """)

    # input localisation happens before the job runs and the output is copied back after,
    # so writing to the input path replaces it rather than racing it
    batch.write_output(align_job.output_cram, str(to_path(cram_path).with_suffix('')))
    jobs.append(align_job)

    return jobs


def record_archive(cram_path: str, sg_id: str, archived_cram: str, repair_image: str) -> None:
    """Note the archive location on the existing `cram` analysis for this CRAM.

    Finds the active `cram` analysis whose output is cram_path and patches its meta.
    Metamist meta updates merge, so existing keys are preserved. No new analysis is
    created - the repaired CRAM supersedes the file this record already describes.

    Runs inside the batch as a PythonJob depending on the alignment job, so the meta is
    only written if the repaired CRAM was actually produced.
    """
    from datetime import datetime, timezone  # noqa: PLC0415

    from metamist.apis import AnalysisApi  # noqa: PLC0415
    from metamist.graphql import gql, query  # noqa: PLC0415
    from metamist.models import AnalysisUpdateModel  # noqa: PLC0415

    find = gql(
        """
        query CramAnalysis($sg_id: String!) {
            sequencingGroups(id: {eq: $sg_id}) {
                analyses(type: {eq: "cram"}, active: {eq: true}) { id outputs }
            }
        }
        """
    )

    for sg in query(find, variables={'sg_id': sg_id})['sequencingGroups']:
        for analysis in sg['analyses']:
            outputs = analysis['outputs']
            path = outputs.get('path') if isinstance(outputs, dict) else outputs
            if path != cram_path:
                continue
            AnalysisApi().update_analysis(
                analysis_id=analysis['id'],
                analysis_update_model=AnalysisUpdateModel(
                    meta={
                        'repair_type': 'trim-adapters',
                        'old_cram_path': archived_cram,
                        'old_cram_index_path': f'{archived_cram}.crai',
                        'repair_script_used': 'src/align_genotype/scripts/repair_scripts/trim_adapters.py',
                        'repair_image': repair_image,
                        'repair_date': datetime.now(timezone.utc).date().isoformat(),
                    },
                ),
            )
            return

    msg = f'No active cram analysis for {sg_id} with output {cram_path}, cannot record the archive'
    raise RuntimeError(msg)


def main() -> None:
    parser = argparse.ArgumentParser(
        description='Trim adapters and poly-G, then realign with DRAGMAP, replacing the CRAM in place.',
    )
    parser.add_argument(
        '--cram-path',
        required=True,
        help='CRAM to repair (GCS path). Needs a .crai beside it. Replaced in place.',
    )
    parser.add_argument('--sg-id', required=True, help='Sequencing group ID, used for the read group and metamist.')
    parser.add_argument('--dry-run', action='store_true', help='Print what would be done without submitting.')
    args = parser.parse_args()

    archived_cram = config.dataset_path(f'{BACKUP_DIR}/{to_path(args.cram_path).name}')

    if args.dry_run:
        print(f'[trim-adapters] {args.sg_id}')
        state = 'already archived' if exists(archived_cram) else 'will archive'
        print(f'  archive {args.cram_path} -> {archived_cram} [{state}]')
        for label, path in (
            ('interleaved FASTQ', config.dataset_path(f'{INTERMEDIATE_DIR}/{args.sg_id}_interleaved.fastq.gz')),
            ('trimmed FASTQ', config.dataset_path(f'{INTERMEDIATE_DIR}/{args.sg_id}_trimmed.fq.gz')),
        ):
            print(f'  {label}: {path} [{"exists, will reuse" if exists(path) else "will create"}]')
        print(f'  repair in place: {args.cram_path}')
        print(f'  note the archive on the existing cram analysis for {args.cram_path}')
        return

    batch = hail_batch.get_batch()
    job_attrs = {'repair_type': 'trim-adapters', 'sequencing_group': args.sg_id}

    archive_job = archive_original(batch, args.cram_path, job_attrs)
    repair_jobs = trim_and_realign(batch, args.cram_path, args.sg_id, job_attrs)

    # nothing references the archive job's outputs, so the ordering must be explicit.
    # applied to every job, since which one runs first depends on what already exists.
    if archive_job:
        for job in repair_jobs:
            job.depends_on(archive_job)

    driver_image = config.config_retrieve(['workflow', 'driver_image'])
    meta_job = batch.new_python_job(f'note archive on cram analysis {args.sg_id}', attributes={'tool': 'metamist'})
    meta_job.image(driver_image)
    meta_job.depends_on(repair_jobs[-1])
    meta_job.call(record_archive, args.cram_path, args.sg_id, archived_cram, driver_image)

    batch.run(wait=False)


if __name__ == '__main__':
    main()
