#!/usr/bin/env python3

"""
Batch jobs to run MultiQC.
"""

from typing import cast

from hailtop.batch import Batch, ResourceFile
from hailtop.batch.job import Job

from cpg_flow import targets, utils
from cpg_utils import Path, hail_batch, config
from cpg_utils.config import get_config, image_path


def multiqc(
    dataset: targets.Dataset,
    tmp_prefix: Path,
    paths: list[Path],
    out_json_path: Path,
    out_html_path: Path,
    out_html_url: str,
    label: str,
    ending_to_trim: set[str] | None = None,
    modules_to_trim_endings: set[str] | None = None,
    job_attrs: dict | None = None,
    sequencing_group_id_map: dict[str, str] | None = None,
    extra_config: dict | None = None,
) -> list[Job]:
    """
    Run MultiQC for the files in `qc_paths`
    @param b: batch object
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
    @param send_to_slack: whether or not to send a Slack message to the qc channel
    @param extra_config: extra config to pass to MultiQC
    @return: job objects
    """

    batch_instance = hail_batch.get_batch()

    title = f'MultiQC [{label}]'

    mqc_j = batch_instance.new_job(title, job_attrs | {'tool': 'MultiQC'})
    mqc_j.image(image_path('multiqc'))
    mqc_j.cpu(16)

    file_list_path = tmp_prefix / f'{dataset.get_alignment_inputs_hash()}_multiqc-file-list.txt'
    dry_run = config.config_retrieve(['workflow', 'dry_run'], False)
    if not dry_run:
        with file_list_path.open('w') as f:
            f.writelines([f'{p}\n' for p in paths])
    file_list = batch_instance.read_input(str(file_list_path))

    endings_conf = ', '.join(list(ending_to_trim)) if ending_to_trim else ''
    modules_conf = ', '.join(list(modules_to_trim_endings)) if modules_to_trim_endings else ''

    if sequencing_group_id_map:
        sample_map_path = tmp_prefix / f'{dataset.get_alignment_inputs_hash()}_rename-sample-map.tsv'
        if not dry_run:
            with sample_map_path.open('w') as fh:
                for sgid, new_sgid in sequencing_group_id_map.items():
                    fh.write('\t'.join([sgid, new_sgid]) + '\n')

        sample_map_file = batch_instance.read_input(str(sample_map_path))
    else:
        sample_map_file = None

    if extra_config:
        extra_config_param = ''
        for k, v in extra_config.items():
            serialised = f'{k}: {v}'
            extra_config_param += f'''--cl-config "{serialised}" \\
            '''
    else:
        extra_config_param = ''

    report_filename = 'report'
    mqc_j.command(f"""\
        mkdir inputs
        cat {file_list} | gsutil -m cp -I inputs/
    
        multiqc -f inputs -o output \\
        {f'--replace-names {sample_map_file} ' if sample_map_file else ''} \\
        --title "{title} for dataset <b>{dataset.name}</b>" \\
        --filename {report_filename}.html \\
        --cl-config "extra_fn_clean_exts: [{endings_conf}]" \\
        --cl-config "max_table_rows: 10000" \\
        --cl-config "use_filename_as_sample_name: [{modules_conf}]" \\
        {extra_config_param}
    
        ls output/{report_filename}_data
        cp output/{report_filename}.html {mqc_j.html}
        cp output/{report_filename}_data/multiqc_data.json {mqc_j.json}
    """)
    batch_instance.write_output(mqc_j.html, str(out_html_path))
    batch_instance.write_output(mqc_j.json, str(out_json_path))

    assert isinstance(mqc_j.json, ResourceFile)
    jobs: list[Job] = [mqc_j]
    if get_config().get('qc_thresholds'):
        check_j = check_report_job(
            b=batch_instance,
            multiqc_json_file=mqc_j.json,
            multiqc_html_url=out_html_url,
            rich_id_map=dataset.rich_id_map(),
            dataset_name=dataset.name,
            label=label,
            job_attrs=job_attrs,
        )
        check_j.depends_on(mqc_j)
        jobs.append(check_j)
    return jobs


def check_report_job(
    b: Batch,
    multiqc_json_file: ResourceFile,
    dataset_name: str,
    multiqc_html_url: str,
    label: str,
    rich_id_map: dict[str, str],
    job_attrs: dict,
) -> Job:
    """
    Run job that checks MultiQC JSON result and sends a Slack notification about failed samples.
    """
    title = f'MultiQC [{label}]'
    check_j = b.new_bash_job(f'{title} check', job_attrs)
    check_j.image(config.config_retrieve(['workflow', 'driver_image']))
    send_to_slack = config.config_retrieve(['workflow', 'send_to_slack'], default=True)

    check_j.command(
        f"""\
    {utils.rich_sequencing_group_id_seds(rich_id_map, [multiqc_json_file])}

    python3 -m align_genotype.scripts.check_multiqc \\
    --multiqc-json {multiqc_json_file} \\
    --html-url {multiqc_html_url} \\
    --dataset {dataset_name} \\
    --title "{title}" \\
    --{'no-' if not send_to_slack else ''}send-to-slack

    touch {check_j.output}
    echo "HTML URL: {multiqc_html_url}"
    """
    )

    return check_j
