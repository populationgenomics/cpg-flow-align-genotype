"""Shared utilities for CRAM repair scripts."""

from loguru import logger

from metamist.graphql import gql, query


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
    from metamist.api import AnalysisApi  # noqa: PLC0415
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
