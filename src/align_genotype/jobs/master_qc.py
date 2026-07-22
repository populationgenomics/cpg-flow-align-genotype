"""
Batch job to generate the master QC HTML report from Metamist QC flags.
"""

from hailtop.batch.job import Job

from cpg_utils import Path, config, hail_batch


def master_qc(
    dataset: str,
    outputs: dict[str, Path],
    job_attrs: dict,
) -> Job:
    """
    Create a Hail Batch job that queries Metamist for all QC flags in the dataset.
    Generates a summary HTML report showing all open and resolved flags.
    """
    batch = hail_batch.get_batch()

    j = batch.new_bash_job('Master QC Report', job_attrs | {'tool': 'python'})
    j.image(config.config_retrieve(['workflow', 'driver_image'])).memory('standard').cpu(2)

    j.command(
        f"""\
    python3 -m align_genotype.scripts.master_qc \\
    --dataset {dataset} \\
    --output-html {j.html}
    """
    )

    batch.write_output(j.html, outputs['html'])
    return j
