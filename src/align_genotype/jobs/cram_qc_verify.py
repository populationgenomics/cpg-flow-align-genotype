"""
Create Hail Batch jobs for VerifyBAMID.
"""

from cpg_flow import resources
from cpg_utils import Path, config, hail_batch
from hailtop.batch.job import Job


def verifybamid(
    cram_path: str,
    output: Path,
    job_attrs: dict,
) -> Job:
    """
    Run `VerifyBamID` contamination checks.
    Based on https://github.com/broadinstitute/warp/blob/57edec5591182d120b7d288b4b409e92a6539871/tasks/broad/BamProcessing.wdl#L395

    Creates a *.selfSM file, a TSV file with 2 rows, 19 columns.
    First row are the keys (e.g., SEQ_SM, RG, FREEMIX), second row are the associated values
    """
    batch_instance = hail_batch.get_batch()

    job = batch_instance.new_bash_job('VerifyBamID', attributes=job_attrs | {'tool': 'VerifyBamID'})
    job.image(config.config_retrieve(['images', 'verifybamid']))
    res = resources.STANDARD.request_resources(ncpu=4, storage_gb=resources.storage_for_cram_qc_job())
    res.set_to_job(job)

    reference = hail_batch.fasta_res_group(batch_instance)

    sequencing_type = config.config_retrieve(['workflow', 'sequencing_type'])
    contam_ud = batch_instance.read_input(config.config_retrieve(['cramqc', f'{sequencing_type}_contam_ud']))
    contam_bed = batch_instance.read_input(config.config_retrieve(['cramqc', f'{sequencing_type}_contam_bed']))
    contam_mu = batch_instance.read_input(config.config_retrieve(['cramqc', f'{sequencing_type}_contam_mu']))

    # define number of PCs used in estimation of contamination
    num_pcs = config.config_retrieve(['cramqc', 'num_pcs'])

    # read in the CRAM and index
    cram_localised = batch_instance.read_input_group(
        cram=cram_path,
        crai=f'{cram_path}.crai',
    ).cram

    job.command(
        f"""\
        /root/micromamba/share/verifybamid2-2.0.1-8/VerifyBamID \
            --NumThread {res.get_nthreads()} \
            --Verbose \
            --NumPC {num_pcs} \
            --Output OUTPUT \
            --BamFile {cram_localised} \
            --Reference {reference.base} \
            --UDPath {contam_ud} \
            --MeanPath {contam_mu} \
            --BedPath {contam_bed} \
            {'--max-depth 1000 ' if sequencing_type == 'exome' else ''} 1>/dev/null

        cp OUTPUT.selfSM {job.out_selfsm}
        """
    )
    batch_instance.write_output(job.out_selfsm, output)
    return job
