"""Batch jobs for QC threshold calibration."""

from hailtop.batch.job import Job

from cpg_utils import Path, config, hail_batch

from align_genotype.qc_calibration.discovery import MultiqcReport


def extract_dataset_metrics(report: MultiqcReport, output: Path, job_attrs: dict) -> Job:
    """Localise one dataset's MultiQC report and distil it to a values file.

    `read_input` rather than a `gcloud storage cp` in the command, matching every other
    job here. Reports run to hundreds of megabytes and parse into several GB of Python
    objects, so this asks for highmem and enough disk for one localised report. The peak
    is dominated by the entry script's `json.load`, not by extraction itself:
    `normalise_sections` shares the leaf per-sample dicts rather than copying them, so
    `extract` adds only cheap scaffolding and the values file it returns is far smaller
    than the report.
    """
    batch = hail_batch.get_batch()

    job = batch.new_bash_job(f'QC calibration extract: {report.dataset}', job_attrs | {'tool': 'python'})
    job.image(config.config_retrieve(['workflow', 'driver_image']))
    job.cpu(2).memory('highmem').storage(config.config_retrieve(['qc_calibration', 'extract_storage'], '20Gi'))

    localised = batch.read_input(report.uri)

    job.command(
        f"""\
    python3 -m align_genotype.scripts.qc_calibration_extract \\
        --dataset {report.dataset} \\
        --multiqc-json {localised} \\
        --analysis-id {report.analysis_id} \\
        --timestamp {report.timestamp} \\
        --uri {report.uri} \\
        --output {job.values}
    """
    )

    batch.write_output(job.values, output)
    return job


def calibration_report(
    values_paths: dict[str, Path],
    skipped_datasets: list[str],
    outputs: dict[str, Path],
    job_attrs: dict,
) -> Job:
    """Read every dataset's values file and render the calibration report.

    Small by construction: the values files are a few hundred KB each, because the
    expensive parse already happened one dataset at a time in `extract_dataset_metrics`.
    """
    batch = hail_batch.get_batch()

    job = batch.new_bash_job('QC calibration report', job_attrs | {'tool': 'python'})
    job.image(config.config_retrieve(['workflow', 'driver_image']))
    job.cpu(2).memory('standard')

    localised = [batch.read_input(str(path)) for path in values_paths.values()]
    values_args = ' '.join(f'--values {resource}' for resource in localised)
    skipped_args = ' '.join(f'--skipped-dataset {name}' for name in skipped_datasets)

    job.command(
        f"""\
    python3 -m align_genotype.scripts.qc_calibration_report \\
        {values_args} {skipped_args} \\
        --output-json {job.report_json} \\
        --output-html {job.report_html}
    """
    )

    batch.write_output(job.report_json, outputs['json'])
    batch.write_output(job.report_html, outputs['html'])
    return job
