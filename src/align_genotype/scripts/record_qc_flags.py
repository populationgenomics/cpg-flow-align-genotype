import json
from dataclasses import asdict
from datetime import datetime

from loguru import logger
from metamist.graphql import gql, query

from align_genotype.utils import QcFlag

DATASET_SG_META_QUERY = gql(
    """
    query datasetSgMeta($dataset: String!) {
        project(name: $dataset) {
            sequencingGroups {
                id
                meta
            }
        }
    }
    """
)

SG_META_MUTATION = gql(
    """
    mutation updateSgMeta($dataset: String!, $sgId: String!, $sgMeta: JSON!) {
        sequencingGroup {
            updateSequencingGroup(
                project: $dataset
                sequencingGroup: {id: $sgId, meta: $sgMeta}
            ) {
                id
                meta
            }
        }
    }
    """
)

def compare_qc_flag(current_flag: dict, new_flag: dict) -> bool:
    """
    Compares two QC flags and returns True if they are the same, False otherwise.
    """
    return (
        current_flag.get('flag') == new_flag.get('flag') and
        current_flag.get('value') == new_flag.get('value') and
        current_flag.get('comparison') == new_flag.get('comparison') and
        current_flag.get('threshold') == new_flag.get('threshold') and
        current_flag.get('section') == new_flag.get('section') and
        not current_flag.get('resolved', False)  # We only want to compare unresolved flags
    )


def main(
    dataset: str,
    qc_flags_json_path: str,
):
    """
    Reads the qc flags JSON file and updates the sequencing group meta in Metamist.

    If the SG meta already has a 'qc_flags' key, it will be updated with the new flags, including
    recording resolution information. If the 'qc_flags' key does not exist, it will be created.
    """
    today = datetime.now()  # noqa: DTZ005

    # Load the QC flags from the JSON file
    with open(qc_flags_json_path) as f:
        qc_flags_data = json.load(f)

    if not qc_flags_data:
        logger.info(f"No QC flags found in {qc_flags_json_path}. No updates will be made.")
        return

    # Query the sequencing groups for the given dataset
    response = query(DATASET_SG_META_QUERY, variables={'dataset': dataset})
    sequencing_groups = response['project']['sequencingGroups']

    # Update each sequencing group's meta with the new QC flags
    for sg in sequencing_groups:
        sg_id = sg['id']
        current_qc_flags: list[dict] = sg['meta'].get('qc_flags', [])
        new_qc_flags: list[dict] = qc_flags_data['qc_flags'].get(sg_id, [])
        new_qc_flags_by_flag = {flag['flag']: flag for flag in new_qc_flags}

        if not current_qc_flags and not new_qc_flags:
            continue  # No existing or new QC flags for this SG, skip

        stats = {"resolved": 0, "retained": 0}
        final_flags: list[QcFlag] = []
        if current_qc_flags:
            logger.info(f"{sg_id} :: Found {len(current_qc_flags)} existing QC flags. Updating existing flags.")
            # Compare the current and new QC flags to determine if an update is needed

            for flag in current_qc_flags:
                if flag['flag'] not in new_qc_flags_by_flag and not flag['resolved']:
                    # If an unresolved current flag is not in the new flags, mark it resolved on the current date
                    flag['resolved'] = True
                    flag['resolution_date'] = today.isoformat()
                    stats['resolved'] += 1
                    logger.info(f"{sg_id} :: Marking QC flag '{flag['flag']}' for SG as resolved.")
                elif compare_qc_flag(flag, new_qc_flags_by_flag[flag['flag']]):
                    # If the current flag matches the new flag, we retain it
                    stats['retained'] += 1
                    logger.info(f"{sg_id} :: QC flag '{flag['flag']}' for SG remains unresolved.")
                else:
                    # If the current flag does not match the new flag, we update it
                    flag.update(new_qc_flags_by_flag[flag['flag']])
                    stats['retained'] += 1
                    logger.info(f"{sg_id} :: QC flag '{flag['flag']}' for SG updated with new information.")
                final_flags.append(QcFlag(**flag))
        else:
            logger.info(f"{sg_id} :: No existing QC flags found, adding {len(new_qc_flags)} new flags.")

        for flag in new_qc_flags:
            flag['resolved'] = False
            flag['resolution_date'] = None
            final_flags.append(QcFlag(**flag))

        # Perform the mutation to update the SG meta
        mutation_response = query(
            SG_META_MUTATION,
            variables={ 'sgId': sg_id, 'sgMeta': {'qc_flags': [asdict(flag) for flag in final_flags]}}
        )
        logger.info(f"{sg_id} :: {len(final_flags)} total flags. Mutation response: {mutation_response}")
        logger.info(f"{sg_id} :: Resolved: {stats['resolved']}, Retained: {stats['retained']}")
