"""Strip /1 and /2 QNAME suffixes from a CRAM.

These suffixes, left by older sequencers, break mate pairing in samtools fastq
without collation.

The repaired CRAM replaces the original in place, so the archived copy is the only
remaining copy of the original - the archive job verifies it by CRC32C before the
repair is allowed to run.

The archive location is noted on that CRAM's existing `cram` analysis. No new analysis
is created. Retiring the analyses derived from the old CRAM is a separate step - see
sg_reset.py.

Standalone by design: no shared imports, so it can be lifted into another workflow
engine as a single process with a CRAM path and a sample ID.
"""

import argparse

from loguru import logger

from hailtop.batch.job import Job

from cpg_flow.utils import exists
from cpg_utils import config, hail_batch, to_path

BACKUP_DIR = 'bad_cram'


def archive_original(batch: hail_batch.Batch, cram_path: str, job_attrs: dict) -> Job | None:
    """Copy the original CRAM and index to bad_cram/ before it is overwritten in place.

    Server-side GCS copy, so the data never transits the worker and the job needs no
    attached storage. Compares CRC32C before and after: the repair overwrites the
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

    # the repair overwrites the original, so refuse to continue on an unverified archive
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


def strip_qnames(batch: hail_batch.Batch, cram_path: str, job_attrs: dict) -> Job:
    """Rewrite the CRAM in place with /1 and /2 stripped from every QNAME."""

    job = batch.new_job('repair CRAM: strip QNAME suffixes', attributes=job_attrs | {'tool': 'samtools'})
    job.image(config.config_retrieve(['images', 'samtools']))
    job.memory('standard')
    job.storage(f'{config.config_retrieve(["workflow", "genome_cram_gb"], "100")}Gi')

    cram_localised = batch.read_input_group(cram=cram_path, crai=f'{cram_path}.crai').cram

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
    samtools view --write-index -@ 3 \
        -T {reference.base} -O cram,version=3.0 \
        -o {job.output_cram.cram} -
    """)

    # input localisation happens before the job runs and the output is copied back after,
    # so writing to the input path replaces it rather than racing it
    batch.write_output(job.output_cram, str(to_path(cram_path).with_suffix('')))
    return job


def record_archive(cram_path: str, sg_id: str, archived_cram: str, repair_image: str) -> None:
    """Note the archive location on the existing `cram` analysis for this CRAM.

    Finds the active `cram` analysis whose output is cram_path and patches its meta.
    Metamist meta updates merge, so existing keys are preserved. No new analysis is
    created - the repaired CRAM supersedes the file this record already describes.

    Runs inside the batch as a PythonJob depending on the repair job, so the meta is
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
                        'repair_type': 'strip-qnames',
                        'old_cram_path': archived_cram,
                        'old_cram_index_path': f'{archived_cram}.crai',
                        'repair_script_used': 'src/align_genotype/scripts/repair_scripts/strip_qnames.py',
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
        description='Strip /1 and /2 QNAME suffixes from a CRAM, replacing it in place.',
    )
    parser.add_argument(
        '--cram-path',
        required=True,
        help='CRAM to repair (GCS path). Needs a .crai beside it. Replaced in place.',
    )
    parser.add_argument('--sg-id', required=True, help='Sequencing group ID, used to find the cram analysis.')
    parser.add_argument('--dry-run', action='store_true', help='Print what would be done without submitting.')
    args = parser.parse_args()

    archived_cram = config.dataset_path(f'{BACKUP_DIR}/{to_path(args.cram_path).name}')

    if args.dry_run:
        print(f'[strip-qnames] {args.sg_id}')
        state = 'already archived' if exists(archived_cram) else 'will archive'
        print(f'  archive {args.cram_path} -> {archived_cram} [{state}]')
        print(f'  repair in place: {args.cram_path}')
        print(f'  note the archive on the existing cram analysis for {args.cram_path}')
        return

    batch = hail_batch.get_batch()
    job_attrs = {'repair_type': 'strip-qnames', 'sequencing_group': args.sg_id}

    archive_job = archive_original(batch, args.cram_path, job_attrs)
    repair_job = strip_qnames(batch, args.cram_path, job_attrs)

    # nothing references the archive job's outputs, so the ordering must be explicit
    if archive_job:
        repair_job.depends_on(archive_job)

    driver_image = config.config_retrieve(['workflow', 'driver_image'])
    meta_job = batch.new_python_job(f'note archive on cram analysis {args.sg_id}', attributes={'tool': 'metamist'})
    meta_job.image(driver_image)
    meta_job.depends_on(repair_job)
    meta_job.call(record_archive, args.cram_path, args.sg_id, archived_cram, driver_image)

    batch.run(wait=False)


if __name__ == '__main__':
    main()
