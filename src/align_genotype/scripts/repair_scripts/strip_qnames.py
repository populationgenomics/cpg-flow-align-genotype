"""Strip /1 and /2 QNAME suffixes from a CRAM.

These suffixes, left by older sequencers, break mate pairing in samtools fastq
without collation.

Archives the original CRAM before rewriting it, and optionally registers the result
as a new `cram` analysis. Retiring the analyses derived from the old CRAM is a
separate step - see sg_reset.py.

Standalone by design: no shared imports, so it can be lifted into another workflow
engine as a single process with an input path, an output path and a sample ID.
"""

import argparse

from loguru import logger

from hailtop.batch.job import Job

from cpg_flow.utils import exists
from cpg_utils import config, hail_batch, to_path

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


def strip_qnames(
    batch: hail_batch.Batch,
    cram_path: str,
    output_cram: str,
    job_attrs: dict,
) -> Job:
    """Rewrite the CRAM with /1 and /2 stripped from every QNAME."""

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
    samtools view --write-index \
        -C -T {reference.base} -@ 3 \
        -o {job.output_cram.cram} -
    """)

    batch.write_output(job.output_cram, str(to_path(output_cram).with_suffix('')))
    return job


def register_cram(output: str, sg_id: str, dataset: str) -> None:
    """Create a completed `cram` analysis for the repaired CRAM.

    Runs inside the batch as a PythonJob depending on the repair job, so registration
    only happens if the CRAM was actually produced. Does not inactivate anything.
    """
    from cpg_flow.metamist import Metamist  # noqa: PLC0415

    Metamist().create_analysis(
        output=output,
        type_='cram',
        status='completed',
        sequencing_group_ids=[sg_id],
        dataset=dataset,
        meta={'source': 'cram-repair', 'repair_type': 'strip-qnames'},
    )


def main() -> None:
    parser = argparse.ArgumentParser(description='Strip /1 and /2 QNAME suffixes from a CRAM.')
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

    if args.dry_run:
        print(f'[strip-qnames] {args.sg_id}  {args.cram_path} -> {args.output_path}')
        print(f'  archive -> {config.dataset_path(f"{BACKUP_DIR}/{to_path(args.cram_path).name}")}')
        print(f'  register: {args.register}')
        return

    if exists(args.output_path):
        logger.info(f'Output already exists, nothing to do: {args.output_path}')
        return

    batch = hail_batch.get_batch()
    job_attrs = {'repair_type': 'strip-qnames', 'sequencing_group': args.sg_id}

    archive_job = archive_original(batch, args.cram_path, job_attrs)
    repair_job = strip_qnames(batch, args.cram_path, args.output_path, job_attrs)

    # nothing references the archive job's outputs, so the ordering must be explicit
    if archive_job:
        repair_job.depends_on(archive_job)

    if args.register:
        reg_job = batch.new_python_job(f'register repaired CRAM {args.sg_id}', attributes={'tool': 'metamist'})
        reg_job.image(config.config_retrieve(['workflow', 'driver_image']))
        reg_job.depends_on(repair_job)
        reg_job.call(register_cram, args.output_path, args.sg_id, args.dataset)

    batch.run(wait=False)


if __name__ == '__main__':
    main()
