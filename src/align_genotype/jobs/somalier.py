"""
Adding jobs for fingerprinting and pedigree checks. Mostly using Somalier.
"""

from hailtop.batch import Resource
from hailtop.batch.job import BashJob

from cpg_flow import targets, utils
from cpg_utils import Path, hail_batch, config

# We want to exclude contaminated sequencing groups from relatedness checks. Somalier is not
# designed to work with contaminated sequencing groups, and in a presence of contamination it
# can generate a lot of false positive families.
MAX_FREEMIX = 0.04


def pedigree(
    dataset: targets.Dataset,
    outputs: dict[str, Path],
    somalier_path_by_sgid: dict[str, Path],
    out_html_url: str | None = None,
    verifybamid_by_sgid: dict[str, Path | str] | None = None,
    label: str | None = None,
    job_attrs: dict | None = None,
) -> list[BashJob]:
    """
    Add somalier and peddy jobs that infer relatedness and sex, compare that
    to the provided PED file, and attempt to recover it. If unable to recover, cancel
    the further workflow jobs.

    Returns a job, a path to a fixed PED file if able to recover, and a path to a file
    with relatedness information for each sample pair

    somalier_path_by_sgid should be a dict of paths to .somalier fingerprints.
    """

    batch_instance = hail_batch.get_batch()

    relate_j = _relate(
        somalier_path_by_sgid=somalier_path_by_sgid,
        verifybamid_by_sgid=verifybamid_by_sgid,
        rich_id_map=dataset.rich_id_map(),
        inputs=outputs,
        label=label,
        job_attrs=job_attrs,
    )

    expected_ped = dataset.write_ped_file(outputs['expected_ped'])

    check_j = _check_pedigree(
        samples_file=relate_j.output_samples,
        pairs_file=relate_j.output_pairs,
        expected_ped=batch_instance.read_input(expected_ped),
        somalier_html_url=out_html_url,
        rich_id_map=dataset.rich_id_map(),
        dataset_name=dataset.name,
        label=label,
        job_attrs=job_attrs,
        check_file=outputs['checks'],
    )
    check_j.depends_on(relate_j)

    return [relate_j, check_j]


def _check_pedigree(
    samples_file: Resource,
    pairs_file: Resource,
    expected_ped: Resource,
    dataset_name: str,
    somalier_html_url: str,
    rich_id_map: dict[str, str],
    label: str,
    job_attrs: dict,
    check_file: Path,
) -> BashJob:
    """
    Run job that checks pedigree and batch correctness. The job will send a Slack message about any mismatches.
    """
    batch_instance = hail_batch.get_batch()

    title = f'Pedigree check [{label}]'

    check_j = batch_instance.new_bash_job(title, job_attrs)
    check_j.image(config.config_retrieve(['workflow', 'driver_image']))

    send_to_slack = config.config_retrieve(['workflow', 'send_to_slack'], default=True)
    cmd = f"""\
    {
        utils.rich_sequencing_group_id_seds(rich_id_map, [str(samples_file), str(pairs_file), str(expected_ped)])
        if rich_id_map
        else ''
    }
    python3 -m align_genotype.scripts.check_pedigree \\
    --somalier-samples {samples_file} \\
    --somalier-pairs {pairs_file} \\
    --expected-ped {expected_ped} \\
    --html-url {somalier_html_url} \\
    --dataset {dataset_name} \\
    --title "{title}" \\
    {'--slack' if send_to_slack else ''}
    touch {check_j.output}
    """
    if somalier_html_url:
        cmd += f'echo "HTML URL: {somalier_html_url}"'

    hail_batch.copy_common_env(check_j)
    hail_batch.authenticate_cloud_credentials_in_job(check_j)
    check_j.command(cmd)
    batch_instance.write_output(check_j.output, check_file)
    return check_j


def _relate(
    somalier_path_by_sgid: dict[str, Path],
    verifybamid_by_sgid: dict[str, Path],
    inputs: dict[str, Path],
    rich_id_map: dict[str, str],
    job_attrs: dict,
    label: str,
) -> BashJob:
    batch_instance = hail_batch.get_batch()

    j = batch_instance.new_bash_job(f'Somalier relate {label}', job_attrs | {'tool': 'somalier'})
    j.image(config.config_retrieve(['images', 'somalier']))
    # Size of one somalier file is 212K, so we add another G only if the number of sequencing groups is >4k
    storage_gb = 1 + len(somalier_path_by_sgid) // 4000
    j.storage(f'{storage_gb}Gi')

    cmd = ''
    input_files_file = '$BATCH_TMPDIR/input_files.list'
    sequencing_groups_ids_file = '$BATCH_TMPDIR/sample_ids.list'
    for sgid, somalier_path in somalier_path_by_sgid.items():
        somalier_file = batch_instance.read_input(somalier_path)
        verify_file = batch_instance.read_input(verifybamid_by_sgid[sgid])
        cmd += f"""
        FREEMIX=$(cat {verify_file} | tail -n1 | cut -f7)
        if [[ $(echo "$FREEMIX > {MAX_FREEMIX}" | bc) -eq 0 ]]; then \
        echo "{somalier_file}" >> {input_files_file}; \
        echo "{sgid}" >> {sequencing_groups_ids_file}; \
        fi
        """

    cmd += f"""
    cat {batch_instance.read_input(inputs['expected_ped'])} | \
    grep -v Family.ID | \
    grep -f {sequencing_groups_ids_file} > expected.ped
    """

    cmd += f"""
    somalier relate \\
    $(cat {input_files_file}) \\
    --ped expected.ped \\
    -o related \\
    --infer
    
    mv related.pairs.tsv {j.output_pairs}
    mv related.samples.tsv {j.output_samples}
    {utils.rich_sequencing_group_id_seds(rich_id_map, ['related.html'])}
    mv related.html {j.output_html}
    """

    j.command(cmd)
    # Copy somalier outputs to final destination.
    batch_instance.write_output(j.output_samples, inputs['samples'])
    batch_instance.write_output(j.output_pairs, inputs['pairs'])
    batch_instance.write_output(j.output_html, inputs['html'])
    return j
