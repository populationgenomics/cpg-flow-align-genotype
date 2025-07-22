"""
Batch jobs to run MultiQC.
"""

from cpg_flow import targets
from cpg_utils import Path, config, hail_batch
from hailtop.batch.job import Job


def multiqc(  # noqa: PLR0913
    dataset: targets.Dataset,
    tmp_prefix: Path,
    paths: list[Path],
    outputs: dict[str, Path],
    label: str,
    ending_to_trim: set[str],
    modules_to_trim_endings: set[str],
    job_attrs: dict,
    sequencing_group_id_map: dict[str, str],
    extra_config: str,
) -> Job:
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
    """

    batch_instance = hail_batch.get_batch()

    title = f'MultiQC [{label}]'

    mqc_j = batch_instance.new_bash_job(title, job_attrs | {'tool': 'MultiQC'})
    mqc_j.image(config.config_retrieve(['images', 'multiqc']))
    mqc_j.cpu(16)

    file_list_path = tmp_prefix / f'{dataset.get_alignment_inputs_hash()}_multiqc-file-list.txt'
    dry_run = config.config_retrieve(['workflow', 'dry_run'], False)
    if not dry_run:
        with file_list_path.open('w') as f:
            f.writelines([f'{p}\n' for p in paths])

    file_list = batch_instance.read_input(file_list_path)

    sample_map_path = tmp_prefix / f'{dataset.get_alignment_inputs_hash()}_rename-sample-map.tsv'
    if not dry_run:
        with sample_map_path.open('w') as fh:
            for sgid, new_sgid in sequencing_group_id_map.items():
                fh.write('\t'.join([sgid, new_sgid]) + '\n')

    sample_map_file = batch_instance.read_input(sample_map_path)

    joined_endings = ', '.join(ending_to_trim)
    joined_modules = ', '.join(modules_to_trim_endings)

    mqc_j.command(
        f"""
        mkdir inputs
        cat {file_list} | gcloud storage cp -I inputs/

        multiqc -f inputs -o output \\
        --replace-names {sample_map_file} \\
        --title "{title} for dataset <b>{dataset.name}</b>" \\
        --filename report.html \\
        --cl-config "extra_fn_clean_exts: [{joined_endings}]" \\
        --cl-config "max_table_rows: 10000" \\
        --cl-config "use_filename_as_sample_name: [{joined_modules}]" \\
        --cl-config "table_columns_visible: {extra_config}"

        cp output/report.html {mqc_j.html}
        cp output/report_data/multiqc_data.json {mqc_j.json}
        """
    )

    batch_instance.write_output(mqc_j.html, outputs['html'])
    batch_instance.write_output(mqc_j.json, outputs['json'])
    return mqc_j
