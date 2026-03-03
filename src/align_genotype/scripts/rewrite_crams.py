"""
A wrapper script that queues Hail batch jobs to rewrite CRAM files
with samtools from v3.1 to v3.0 format for compatibility with downstream tools.

Queries datasets from CRAMs written after 25th February 2026.
"""

import argparse
from datetime import datetime, timezone

from cpg_utils import config, hail_batch
from hailtop.batch.job import Job
from metamist.graphql import gql, query

# The date the new Picard image with updated samtools was deployed, which produced CRAM v3.1 by default.
FROM_DATE = datetime(2026, 2, 25, tzinfo=timezone.utc)


def rewrite_cram(
    batch: hail_batch.Batch,
    cram_path: str,
    job_attrs: dict,
) -> Job:
    """Rewrite CRAM file to version 3.0 format using samtools."""

    job = batch.new_job(
        'rewrite CRAM to v3.0',
        attributes=job_attrs | {'tool': 'samtools'},
    )

    job.image(config.image_path('samtools', '1.18-1'))

    # read in the CRAM and index
    cram_localised = batch.read_input_group(
        cram=cram_path,
        crai=f'{cram_path}.crai',
    ).cram

    reference = hail_batch.fasta_res_group(batch)

    job.command(f"""\
    samtools view \\
        -@{job.attrs.get('nthreads', 4) - 1} \\
        --reference {reference.base} \\
        -O cram,version=3.0 \\
        -o {job.output_cram.cram} \\
        {cram_localised}
    """)
    # Overwrite the original CRAM file with the rewritten version.
    batch.write_output(job.output_cram.cram, cram_path)
    return job


def get_crams_to_rewrite(datasets: list[str]) -> list[str]:
    """Query the database for CRAM files that need to be rewritten."""
    query_str = gql(
        """
        query GetCramsToRewrite(dataset: String!) {
            project(name: $dataset) {
                sequencingGroups(type: {eq: "genome"}) {
                    id
                    analyses(type: {eq: "cram"}, status: {eq: "completed"}) {
                        id
                        timestampCompleted
                        outputs
                    }
                }
            }
        }
        """
    )
    crams_to_rewrite = []
    for dataset in datasets:
        variables = {'dataset': dataset}
        result = query(query_str, variables)
        for sg in result['project']['sequencingGroups']:
            for analysis in sg['analyses']:
                timestamp = datetime.fromisoformat(analysis['timestampCompleted']).replace(tzinfo=timezone.utc)
                if timestamp > FROM_DATE:
                    crams_to_rewrite.append(analysis['outputs']['path'])

    return crams_to_rewrite


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Rewrite CRAM files to version 3.0 format.')
    parser.add_argument(
        '--datasets',
        nargs='+',
        required=True,
        help='List of datasets to query for CRAM files to rewrite.',
    )
    parser.add_argument('--dry-run', action='store_true', help='Print the CRAM files that would be rewritten.')
    args = parser.parse_args()

    crams_to_rewrite = get_crams_to_rewrite(args.datasets)

    if args.dry_run:
        for cram_path in crams_to_rewrite:
            print(cram_path)
    else:
        batch = hail_batch.get_batch()
        for cram_path in crams_to_rewrite:
            rewrite_cram(batch, cram_path, job_attrs={'nthreads': 4})
        batch.run()
