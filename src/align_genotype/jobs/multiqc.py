"""
Batch jobs to run MultiQC.
"""

from cpg_flow import resources, targets
from cpg_utils import Path, config, hail_batch
from hailtop.batch import Batch, ResourceFile
from hailtop.batch.job import Job


def multiqc(
    dataset: targets.Dataset,
    tmp_prefix: Path,
    paths: list[Path],
    outputs: dict[str, Path],
    out_html_url: str | None,
    out_checks_path: Path | None,
    label: str,
    ending_to_trim: set[str],
    modules_to_trim_endings: set[str],
    job_attrs: dict,
    sequencing_group_id_map: dict[str, str],
    extra_config: dict,
    send_to_slack: bool = True,
) -> list[Job]:
    """
    Run MultiQC for the files in `qc_paths`
    @param tmp_prefix: bucket for tmp files
    @param paths: file bucket paths to pass into MultiQC
    @param dataset: Dataset object
    @param out_json_path: where to write MultiQC-generated JSON file
    @param out_html_path: where to write the HTML report
    @param out_html_url: URL corresponding to the HTML report
    @param out_checks_path: flag indicating that QC checks were done
    @param label: To add to the report's, Batch job's, and Slack message's titles
    @param ending_to_trim: trim these endings from input files to get sequencing group names
    @param modules_to_trim_endings: list of modules for which trim the endings
    @param job_attrs: attributes to add to Hail Batch job
    @param sequencing_group_id_map: sequencing group ID map for bulk sequencing group renaming:
        (https://multiqc.info/docs/#bulk-sample-renaming-in-reports)
    @param extra_config: additional configuration options for MultiQC
    @param send_to_slack: whether or not to send a Slack message to the qc channel
    """

    batch_instance = hail_batch.get_batch()

    title = f'MultiQC [{label}]'

    mqc_j = batch_instance.new_bash_job(title, job_attrs | {'tool': 'MultiQC'})
    mqc_j.image(config.config_retrieve(['images', 'multiqc']))
    mqc_j.cpu(8)
    mqc_j.storage('20Gi')

    file_list_path = tmp_prefix / f'{dataset.get_alignment_inputs_hash()}_multiqc-file-list.txt'
    dry_run = config.config_retrieve(['workflow', 'dry_run'], False)
    if not dry_run:
        with file_list_path.open('w') as f:
            f.writelines([f'{p}\n' for p in paths])

    file_list = batch_instance.read_input(file_list_path)

    sg_id_mapping_file_path = tmp_prefix / f'{dataset.get_alignment_inputs_hash()}_rename-sg-map.tsv'
    if not dry_run:
        with sg_id_mapping_file_path.open('w') as fh:
            for sgid, new_sgid in sequencing_group_id_map.items():
                fh.write('\t'.join([sgid, new_sgid]) + '\n')

    sg_id_mapping_file = batch_instance.read_input(sg_id_mapping_file_path)

    joined_endings = ', '.join(ending_to_trim)
    joined_modules = ', '.join(modules_to_trim_endings)

    if extra_config:
        extra_config_param = ''
        for k, v in extra_config.items():
            serialised = f'{k}: {v}'
            extra_config_param += f'''--cl-config "{serialised}" \\
            '''  # noqa: Q001
    else:
        extra_config_param = ''

    mqc_j.command(
        f"""
        mkdir inputs
        cat {file_list} | gcloud storage cp -I inputs/

        multiqc -f inputs -o output \\
        --replace-names {sg_id_mapping_file} \\
        --title "{title} for dataset <b>{dataset.name}</b>" \\
        --filename report.html \\
        --cl-config "extra_fn_clean_exts: [{joined_endings}]" \\
        --cl-config "max_table_rows: 10000" \\
        --cl-config "use_filename_as_sample_name: [{joined_modules}]" \\
        {extra_config_param}

        cp output/report.html {mqc_j.html}
        cp output/report_data/multiqc_data.json {mqc_j.json}
        """
    )
    if out_html_url:
        mqc_j.command(f'echo "HTML URL: {out_html_url}"')

    batch_instance.write_output(mqc_j.html, outputs['html'])
    batch_instance.write_output(mqc_j.json, outputs['json'])
    jobs: list[Job] = [mqc_j]
    if config.config_retrieve(['qc_thresholds']):
        check_j = check_report_job(
            b=batch_instance,
            multiqc_json_file=mqc_j.json,
            multiqc_html_url=out_html_url,
            rich_id_map=dataset.rich_id_map(),
            dataset_name=dataset.name,
            label=label,
            out_checks_path=out_checks_path,
            job_attrs=job_attrs,
            send_to_slack=send_to_slack,
        )
        check_j.depends_on(mqc_j)
        jobs.append(check_j)

        # Important to add the -test suffix as dataset_name is used in GraphQL queries
        test = config.config_retrieve(['workflow', 'access_level'], None) == 'test'
        dataset_name = dataset.name + '-test' if test else dataset.name
        record_j = record_qc_flags_job(
            b=batch_instance,
            dataset_name=dataset_name,
            label=label,
            sg_id_mapping_file=sg_id_mapping_file,
            check_multiqc_json_file=check_j.output,
            job_attrs=job_attrs,
        )
        record_j.depends_on(check_j)
        jobs.append(record_j)

    return jobs


def check_report_job(
    b: Batch,
    multiqc_json_file: ResourceFile,
    dataset_name: str,
    multiqc_html_url: str | None = None,
    label: str | None = None,
    rich_id_map: dict[str, str] | None = None,
    out_checks_path: Path | None = None,
    job_attrs: dict | None = None,
    send_to_slack: bool = True,
) -> Job:
    """
    Run job that checks MultiQC JSON result and sends a Slack notification
    about failed samples.
    """
    title = 'MultiQC'
    if label:
        title += f' [{label}]'
    check_j = b.new_job(f'{title} check', (job_attrs or {}) | {'tool': 'python'})
    resources.STANDARD.set_resources(j=check_j, ncpu=2)
    check_j.image(config.config_retrieve(['workflow', 'driver_image']))

    cmd = f"""\
    {rich_sequencing_group_id_seds(rich_id_map, [multiqc_json_file]) if rich_id_map else ''}

    python3 -m align_genotype.scripts.check_multiqc \\
    --multiqc-json {multiqc_json_file} \\
    --html-url {multiqc_html_url} \\
    --dataset {dataset_name} \\
    --title "{title}" \\
    --output-json {check_j.output} \\
    --{'no-' if not send_to_slack else ''}send-to-slack

    echo "HTML URL: {multiqc_html_url}"
    """

    check_j.command(cmd)

    if out_checks_path:
        b.write_output(check_j.output, str(out_checks_path))
    return check_j


def record_qc_flags_job(
    b: Batch,
    dataset_name: str,
    label: str,
    sg_id_mapping_file: ResourceFile,
    check_multiqc_json_file: ResourceFile,
    job_attrs: dict | None = None,
) -> Job:
    """
    Run job that records QC flags in Metamist.
    """
    report_name = 'CramMultiQC' if label == 'CRAM' else 'GvcfMultiQC' if label == 'GVCF' else 'MultiQC'
    record_j = b.new_job('Record QC flags', (job_attrs or {}) | {'tool': 'python'})
    resources.STANDARD.set_resources(j=record_j, ncpu=2)
    record_j.image(config.config_retrieve(['workflow', 'driver_image']))

    cmd = f"""\
    python3 -m align_genotype.scripts.record_qc_flags \\
    --dataset {dataset_name} \\
    --report {report_name} \\
    --qc-flags-json {check_multiqc_json_file} \\
    --sequencing-group-ids-map {sg_id_mapping_file}
    """

    record_j.command(cmd)
    return record_j


def rich_sequencing_group_id_seds(
    rich_id_map: dict[str, str],
    file_names: list[str | ResourceFile],
) -> str:
    """
    Helper function to add seds into a command that would extend sequencing group IDs
    in each file in `file_names` with an external ID, only if external ID is
    different from the original.

    @param rich_id_map: map used to replace sequencing groups, e.g. {'CPGAA': 'CPGAA|EXTID'}
    @param file_names: file names and Hail Batch Resource files where to replace IDs
    @return: bash command that does replacement
    """
    cmd = ''
    for sgid, rich_sgid in rich_id_map.items():
        for fname in file_names:
            cmd += f"sed -iBAK 's/{sgid}/{rich_sgid}/g' {fname}"
            cmd += '\n'
    return cmd
