"""Shared utilities for CRAM repair scripts."""

from loguru import logger

from hailtop.batch.job import Job

from cpg_flow.utils import exists
from cpg_utils import config, hail_batch, to_path
from metamist.graphql import gql, query

BACKUP_DIR = 'bad_cram'


def backup_cram_path(cram_path: str) -> str:
    """Derive the archive path for the pre-repair CRAM, under bad_cram/.

    Built with cpg_utils dataset_path, so the bucket namespace follows the
    analysis-runner --dataset and --access-level rather than being derived from
    the input path: test runs archive to -test, full runs to -main.
    """
    return config.dataset_path(f'{BACKUP_DIR}/{to_path(cram_path).name}')


def backup_original_cram(
    batch: hail_batch.Batch,
    cram_path: str,
    job_attrs: dict,
) -> Job | None:
    """Archive the original CRAM and its index to bad_cram/ before any repair runs.

    Uses a server-side GCS copy, so the data never transits the worker and the job
    needs no attached storage. Returns None if the archive already exists.
    """

    dest = backup_cram_path(cram_path)
    dest_dir = str(to_path(dest).parent)

    if exists(dest) and exists(f'{dest}.crai'):
        logger.info(f'Skipping CRAM backup: already archived at {dest}')
        return None

    job = batch.new_bash_job(
        'repair CRAM: archive original',
        attributes=job_attrs | {'tool': 'gcloud'},
    )
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


def get_cram_paths_for_sgs(sg_ids: list[str]) -> list[tuple[str, str, int]]:
    """Query metamist for CRAM paths by sequencing group IDs.

    Returns (sg_id, cram_path, analysis_id) tuples.
    """
    query_str = gql(
        """
        query GetCramPaths($sg_id: String!) {
            sequencingGroups(id: {eq: $sg_id}) {
                id
                analyses(type: {eq: "cram"}, status: {eq: COMPLETED}) {
                    id
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
                latest = analyses[0]
                outputs = latest['outputs']
                path = outputs.get('path') if isinstance(outputs, dict) else outputs
                if path:
                    results.append((sg['id'], path, latest['id']))
    return results


def register_and_inactivate(
    new_cram_path: str,
    sg_id: str,
    dataset: str,
    old_analysis_id: int,
    repair_type: str,
) -> None:
    """Create a new CRAM analysis entry and inactivate the old one.

    Intended to run inside a Hail Batch PythonJob.
    """
    from cpg_flow.metamist import Metamist  # noqa: PLC0415
    from metamist.apis import AnalysisApi  # noqa: PLC0415
    from metamist.models import AnalysisUpdateModel  # noqa: PLC0415

    m = Metamist()
    aid = m.create_analysis(
        output=new_cram_path,
        type_='cram',
        status='completed',
        sequencing_group_ids=[sg_id],
        dataset=dataset,
        meta={'source': 'cram-repair', 'repair_type': repair_type},
    )

    if aid is None:
        logger.error(f'Failed to create new CRAM entry for {sg_id}, skipping inactivation')
        return

    aapi = AnalysisApi()
    aapi.update_analysis(
        analysis_id=old_analysis_id,
        analysis_update_model=AnalysisUpdateModel(active=False),
    )
    logger.info(f'Inactivated old analysis {old_analysis_id} for {sg_id}')
