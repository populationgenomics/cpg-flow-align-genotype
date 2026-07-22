import json
from dataclasses import asdict
from datetime import datetime

import click
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
    Returns True if the two flags refer to the same QC issue and the current flag is still
    unresolved. Flag identity is (flag + comparison + threshold + section); the measured
    `value` is intentionally excluded because it can drift slightly between MultiQC runs.
    """
    return (
        current_flag.get('flag') == new_flag.get('flag')
        and current_flag.get('comparison') == new_flag.get('comparison')
        and current_flag.get('threshold') == new_flag.get('threshold')
        and current_flag.get('section') == new_flag.get('section')
        and not current_flag.get('resolved', False)
    )


def reconcile_sg_qc_flags(
    sg: dict,
    new_flags_by_sg: dict[str, list[dict]],
    dataset: str,
    label: str,
    today: datetime,
) -> None:
    """
    Reconcile current and new QC flags for a single SG and update its meta in Metamist.

    Only considers flags of the type specified by the label (CRAM or GVCF).

    Marks absent flags as resolved, retains unchanged flags, updates differing flags,
    and adds new flags that didn't previously exist.
    """
    sg_id = sg['id']
    report = 'CramMultiQC' if label == 'CRAM' else 'GvcfMultiQC' if label == 'GVCF' else 'MultiQC'

    # Get all the existing QC flags of the specified type for this SG
    qc_flags_key = f'{label.lower()}_qc_flags'
    current_qc_flags: list[dict] = (sg['meta'] or {}).get(qc_flags_key, [])
    unresolved_current_flags = [flag for flag in current_qc_flags if not flag.get('resolved', False)]

    new_qc_flags: list[dict] = new_flags_by_sg.get(sg_id, [])

    if not current_qc_flags and not new_qc_flags:
        logger.info(f'{sg_id} :: No existing or new {report} flags for this SG, skipping.')
        return  # No existing or new QC flags for this SG, skip

    if not unresolved_current_flags and not new_qc_flags:
        logger.info(f'{sg_id} :: No unresolved existing or new {report} flags for this SG, skipping.')
        return  # No unresolved existing or new QC flags for this SG, skip

    # Key flags by (section, flag) so the same metric in different MultiQC sections
    # is treated as a distinct QC issue.
    new_qc_flags_by_key = {(flag['section'], flag['flag']): flag for flag in new_qc_flags}
    existing_flag_keys = {(flag['section'], flag['flag']) for flag in current_qc_flags}

    # Track the final set of flags to be recorded in Metamist, including resolved, retained, updated, and added flags
    final_flags: list[QcFlag] = []
    stats = {'resolved': 0, 'retained': 0, 'updated': 0, 'added': 0}

    if current_qc_flags:
        logger.info(f'{sg_id} :: Found {len(current_qc_flags)} existing {report} flags. Reconciling.')
        for flag in current_qc_flags:
            key = (flag['section'], flag['flag'])
            present_in_new = key in new_qc_flags_by_key
            if not present_in_new:
                if not flag['resolved']:
                    # Previously-flagged issue no longer present: mark resolved
                    flag['resolved'] = True
                    flag['resolution_date'] = today.isoformat(timespec='seconds')
                    logger.info(f"{sg_id} :: Marking {report} flag '{flag['flag']}' as resolved.")
                    stats['resolved'] += 1
                else:
                    # Already resolved and still absent: keep as-is
                    logger.info(f"{sg_id} :: {report} flag '{flag['flag']}' remains resolved.")
            elif compare_qc_flag(flag, new_qc_flags_by_key[key]):
                # Same unresolved issue is still present: refresh the measured value and
                # but keep resolution status. Identity (metric/threshold/section/
                # comparison) is unchanged so this counts as 'retained', not 'updated'.
                new_flag = new_qc_flags_by_key[key]
                flag['value'] = new_flag['value']
                logger.info(f"{sg_id} :: {report} flag '{flag['flag']}' remains unresolved (value refreshed).")
                stats['retained'] += 1
            else:
                # Current flag exists in new run but differs (or was resolved and has reappeared):
                # overwrite with new flag data (which sets resolved=False)
                flag.update(new_qc_flags_by_key[key])
                logger.info(f"{sg_id} :: {report} flag '{flag['flag']}' updated with new information.")
                stats['updated'] += 1
            final_flags.append(QcFlag(**flag))
    else:
        logger.info(f'{sg_id} :: No existing {report} flags found, adding {len(new_qc_flags)} new flags.')

    # Add new flags that aren't already represented in current_qc_flags
    for flag in new_qc_flags:
        if (flag['section'], flag['flag']) in existing_flag_keys:
            continue
        logger.info(f"{sg_id} :: Adding new {report} flag '{flag['flag']}'.")
        final_flags.append(QcFlag(**flag))
        stats['added'] += 1

    # Perform the mutation to update the SG meta
    query(
        SG_META_MUTATION,
        variables={
            'dataset': dataset,
            'sgId': sg_id,
            'sgMeta': {qc_flags_key: [asdict(flag) for flag in final_flags]},
        },
    )
    logger.info(
        f'{sg_id} :: Recorded {len(final_flags)} {report} flags in Metamist. '
        f'Resolved: {stats["resolved"]}, Retained: {stats["retained"]}, '
        f'Updated: {stats["updated"]}, Added: {stats["added"]}'
    )


@click.command()
@click.option('--dataset', required=True, help='Dataset name')
@click.option('--label', required=True, help='Report type (CRAM or GVCF)')
@click.option(
    '--qc-flags-json',
    'qc_flags_json_path',
    required=True,
    help='Path to the QC flags JSON file',
)
@click.option(
    '--sequencing-group-ids-map',
    'sg_id_mapping_file',
    required=True,
    help='Path to the sequencing group IDs mapping file',
)
def main(
    dataset: str,
    label: str,
    qc_flags_json_path: str,
    sg_id_mapping_file: str,
):
    """
    Reads the qc flags JSON file and the SG mapping file, and updates any flagged QC issues in the
    sequencing group meta in Metamist. The mapping file is required to ensure that the SG was in
    scope for the dataset being processed and to avoid updating unrelated SGs.

    If the SG meta already has a 'qc_flags' key, it will be updated with the new flags, including
    recording resolution information. If the 'qc_flags' key does not exist, it will be created.
    """
    today = datetime.now()  # noqa: DTZ005

    # Load the sequencing group IDs mapping file tsv
    with open(sg_id_mapping_file) as f:
        sg_id_map = dict(line.strip().split('\t') for line in f if line.strip())

    # Load the QC flags from the JSON file
    with open(qc_flags_json_path) as f:
        qc_flags_data = json.load(f)

    # Query the sequencing groups for the given dataset
    response = query(DATASET_SG_META_QUERY, variables={'dataset': dataset})
    sequencing_groups = response['project']['sequencingGroups']

    # Reconcile each sequencing group's QC flags
    new_flags_by_sg = qc_flags_data['qc_flags']
    for sg in sequencing_groups:
        if sg['id'] in sg_id_map:
            reconcile_sg_qc_flags(sg, new_flags_by_sg, dataset, label, today)


if __name__ == '__main__':
    main()  # pylint: disable=E1120
