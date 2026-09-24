"""Copy a CRAM and rewrite its @RG read group to a different sequencing group ID.

Used to promote a CRAM repaired under a test sequencing group to the main
sequencing group for the same sample. CRAM stores RG as an index into the
header's @RG list, so rewriting the header re-labels every record without
re-encoding the alignment data.
"""

import argparse

from hailtop.batch.job import Job

from cpg_utils import config, hail_batch, to_path


def reheader(
    batch: hail_batch.Batch,
    cram_path: str,
    old_sg: str,
    new_sg: str,
    output_cram: str,
) -> Job:
    """Rewrite @RG ID/SM from old_sg to new_sg, writing an indexed CRAM to output_cram."""

    job = batch.new_job(
        f'reheader CRAM {old_sg} -> {new_sg}',
        attributes={'tool': 'samtools', 'old_sg': old_sg, 'new_sg': new_sg},
    )

    nthreads = 4
    job.image(config.config_retrieve(['images', 'samtools']))
    job.cpu(nthreads)
    job.memory('standard')
    job.storage(f'{config.config_retrieve(["workflow", "reheader_storage_gb"], 100)}Gi')

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

    # Only touch @RG ID/SM. @PG lines are left alone so the original command
    # lines remain an accurate record of how the data was actually produced.
    awk_rg = (
        'awk \'BEGIN{OFS="\\t"} $1=="@RG"{for(i=2;i<=NF;i++){'
        f'if($i ~ /^ID:/) $i="ID:{new_sg}"; if($i ~ /^SM:/) $i="SM:{new_sg}"'
        '}} {print}\''
    )

    job.command(f"""\
    set -eo pipefail

    samtools view -H {cram_localised} > old_header.sam

    grep -c "{old_sg}" old_header.sam || true
    {awk_rg} old_header.sam > new_header.sam
    printf '@CO\\treheadered {old_sg} -> {new_sg} from {cram_path}\\n' >> new_header.sam

    echo "--- @RG before ---"; grep '^@RG' old_header.sam
    echo "--- @RG after  ---"; grep '^@RG' new_header.sam

    samtools reheader new_header.sam {cram_localised} > {job.output_cram.cram}
    samtools index -@ {nthreads} {job.output_cram.cram}

    echo "--- verifying written header ---"
    samtools view -H {job.output_cram.cram} | grep '^@RG'
    if samtools view -H {job.output_cram.cram} | grep -q '^@RG.*{old_sg}'; then
        echo "FATAL: {old_sg} still present in @RG after reheader" >&2
        exit 1
    fi
    """)

    batch.write_output(job.output_cram, str(to_path(output_cram).with_suffix('')))
    return job


def main() -> None:
    parser = argparse.ArgumentParser(description='Copy a CRAM and rewrite its @RG to a new sequencing group ID.')
    parser.add_argument('--cram-path', required=True, help='Source CRAM (GCS path).')
    parser.add_argument('--old-sg', required=True, help='Sequencing group ID currently in the @RG lines.')
    parser.add_argument('--new-sg', required=True, help='Sequencing group ID to write into @RG ID/SM.')
    parser.add_argument('--output-path', required=True, help='Destination CRAM (GCS path).')
    parser.add_argument('--dry-run', action='store_true', help='Print what would be done without submitting.')
    args = parser.parse_args()

    if args.dry_run:
        print(f'{args.cram_path} -> {args.output_path}  (@RG {args.old_sg} -> {args.new_sg})')
        return

    batch = hail_batch.get_batch()
    reheader(batch, args.cram_path, args.old_sg, args.new_sg, args.output_path)
    batch.run(wait=False)


if __name__ == '__main__':
    main()