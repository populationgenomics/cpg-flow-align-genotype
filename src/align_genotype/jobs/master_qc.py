"""
Batch job to generate the master QC HTML report from Metamist QC flags.
"""

from cpg_flow import resources
from cpg_utils import Path, config, hail_batch
from hailtop.batch.job import Job


def master_qc(
    dataset_name: str,
    outputs: dict[str, Path],
    out_html_url: str | None,
    job_attrs: dict,
) -> list[Job]:
    """
    Create a Hail Batch job that queries Metamist for all QC flags
    in the dataset and generates a summary HTML report.
    """
    batch = hail_batch.get_batch()

    j = batch.new_bash_job('Master QC Report', job_attrs | {'tool': 'python'})
    resources.STANDARD.set_resources(j=j, ncpu=2)
    j.image(config.config_retrieve(['workflow', 'driver_image']))

    j.command(
        f"""\
    python3 -m align_genotype.scripts.master_qc \\
    --dataset {dataset_name} \\
    --output-html {j.html}
    """
    )

    if out_html_url:
        j.command(f'echo "Master QC HTML URL: {out_html_url}"')

    batch.write_output(j.html, str(outputs['html']))
    return [j]
