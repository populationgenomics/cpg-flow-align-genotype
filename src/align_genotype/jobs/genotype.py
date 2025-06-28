"""
CRAM to GVCF: create Hail Batch jobs to genotype individual sequencing groups.
"""

import logging

import hailtop.batch as hb
from cpg_flow import filetypes, resources, utils
from cpg_utils import Path, config, hail_batch
from hailtop.batch.job import BashJob

from align_genotype.jobs.picard import get_intervals


def genotype(
    sequencing_group_name: str,
    tmp_prefix: Path,
    cram_path: str,
    output_path: Path,
    job_attrs: dict[str, str],
) -> list[BashJob]:
    """
    Takes a CRAM file and runs GATK tools to make a GVCF.
    """
    hc_gvcf_path = tmp_prefix / 'haplotypecaller' / f'{sequencing_group_name}.g.vcf.gz'

    # collect all the real jobs - [reuse] jobs substituted for None
    jobs = [
        job
        for job in haplotype_caller(
            sequencing_group_name=sequencing_group_name,
            job_attrs=job_attrs,
            output_path=hc_gvcf_path,
            cram_path=cram_path,
            tmp_prefix=tmp_prefix,
            scatter_count=config.config_retrieve(['workflow', 'scatter_count_genotype']),
        )
        if job is not None
    ]
    postproc_j = postprocess_gvcf(
        gvcf_path=filetypes.GvcfPath(hc_gvcf_path),
        sequencing_group_name=sequencing_group_name,
        job_attrs=job_attrs,
        output_path=output_path,
    )

    # only keep elements which are not None
    # if both exist, set the dependency
    if postproc_j and jobs:
        postproc_j.depends_on(*jobs)

    if postproc_j:
        jobs.append(postproc_j)

    return jobs


intervals: list[hb.ResourceFile] | None = None


def haplotype_caller(  # noqa: PLR0913
    sequencing_group_name: str,
    cram_path: str,
    tmp_prefix: Path,
    scatter_count: int,
    job_attrs: dict[str, str],
    output_path: Path,
) -> list[BashJob]:
    """
    Run GATK Haplotype Caller in parallel, split by intervals.
    """

    batch_instance = hail_batch.get_batch()
    jobs = []

    if scatter_count > 1:
        global intervals
        if intervals is None:
            intervals_j, intervals = get_intervals(
                b=batch_instance,
                scatter_count=scatter_count,
                job_attrs=job_attrs,
                output_prefix=tmp_prefix / f'intervals_{scatter_count}',
            )
            if intervals_j:
                jobs.append(intervals_j)

        hc_fragments = []
        # Splitting variant calling by intervals
        for idx in range(scatter_count):
            # give each fragment a tmp location
            fragment = tmp_prefix / 'haplotypecaller' / f'{idx}_of_{scatter_count}_{sequencing_group_name}.g.vcf.gz'
            j, result = _haplotype_caller_one(
                sequencing_group_name=sequencing_group_name,
                cram_path=cram_path,
                job_attrs=job_attrs | {'part': f'{idx + 1}/{scatter_count}'},
                interval=intervals[idx],
                out_gvcf_path=fragment,
            )
            hc_fragments.append(result)

            # only consider jobs which weren't scheduled for [reuse]
            if j:
                jobs.append(j)

        jobs.append(
            merge_gvcfs_job(
                batch_instance=batch_instance,
                sequencing_group_name=sequencing_group_name,
                gvcf_groups=hc_fragments,
                job_attrs=job_attrs,
                out_gvcf_path=output_path,
            )
        )
    else:
        hc_j, _result = _haplotype_caller_one(
            sequencing_group_name=sequencing_group_name,
            job_attrs=job_attrs,
            cram_path=cram_path,
            out_gvcf_path=output_path,
        )
        if hc_j:
            jobs.append(hc_j)

    return jobs


def _haplotype_caller_one(
    sequencing_group_name: str,
    cram_path: str,
    job_attrs: dict,
    out_gvcf_path: Path,
    interval: hb.Resource | None = None,
) -> tuple[BashJob | None, hb.Resource]:
    """
    Add one GATK HaplotypeCaller job on an interval.
    """

    batch_instance = hail_batch.get_batch()

    if utils.can_reuse(out_gvcf_path):
        logging.info(f'Reusing HaplotypeCaller {out_gvcf_path}, job attrs: {job_attrs}')
        # localise the output GVCF and return it
        gvcf_group = batch_instance.read_input_group(
            **{'g.vcf.gz': str(out_gvcf_path), 'g.vcf.gz.tbi': f'{out_gvcf_path!s}.tbi'}
        )
        return None, gvcf_group

    job = batch_instance.new_bash_job('HaplotypeCaller', job_attrs | {'tool': 'gatk HaplotypeCaller'})
    job.image(config.config_retrieve(['images', 'gatk']))

    # Enough storage to localize CRAMs (can't pass GCS URL to CRAM to gatk directly
    # because we will hit GCP egress bandwidth limit:
    # https://batch.hail.populationgenomics.org.au/batches/7493/jobs/2)
    # CRAMs can be as big as 80G:
    # https://batch.hail.populationgenomics.org.au/batches/74042/jobs/3346
    # plus we need enough space to fit output GVCF (1G) and reference data (5G).
    # HaplotypeCaller is not parallelized, so we request the minimal possible chunk
    # of a Hail worker (2 cores), plus with the `highmem` instance that would
    # give enough memory: 13G. That's not going to give enough disk storage, so we
    # are explicitly requesting more storage.
    #
    # Based on an audit of RD crams on 19/05/23, 99% of crams are <34Gb. Will set the
    # default to 40Gb for genomes then use a run specific confg to run the rare
    # sequencing group that will fail from this limit.
    storage_default = 40 if config.config_retrieve(['workflow', 'sequencing_type']) == 'genome' else None
    # enough for input CRAM and output GVCF
    job_res = resources.HIGHMEM.request_resources(
        ncpu=2,
        storage_gb=config.config_retrieve(['workflow', 'haplotypecaller_storage'], storage_default),
    )
    job_res.set_to_job(job)

    job.declare_resource_group(
        output_gvcf={
            'g.vcf.gz': '{root}-' + sequencing_group_name + '.g.vcf.gz',
            'g.vcf.gz.tbi': '{root}-' + sequencing_group_name + '.g.vcf.gz.tbi',
        },
    )

    reference = hail_batch.fasta_res_group(batch_instance)

    cmd = f"""\
    CRAM=$BATCH_TMPDIR/{sequencing_group_name}.cram
    CRAI=$BATCH_TMPDIR/{sequencing_group_name}.cram.crai

    # Retrying copying to avoid google bandwidth limits
    retry_gs_cp {cram_path} $CRAM
    retry_gs_cp {f'{cram_path}.crai'} $CRAI

    gatk --java-options \
    "{job_res.java_mem_options()} \
    -XX:GCTimeLimit=50 \
    -XX:GCHeapFreeLimit=10" \\
    HaplotypeCaller \\
    -R {reference.base} \\
    -I $CRAM \\
    --read-index $CRAI \\
    {f'-L {interval} ' if interval is not None else ''} \\
    --disable-spanning-event-genotyping \\
    --dragen-mode \\
    -O {job.output_gvcf['g.vcf.gz']} \\
    -G AS_StandardAnnotation \\
    -GQB 10 -GQB 20 -GQB 30 -GQB 40 -GQB 50 -GQB 60 -GQB 70 -GQB 80 -GQB 90 \\
    -ERC GVCF \\
    --create-output-variant-index
    """
    job.command(hail_batch.command(cmd, monitor_space=True, setup_gcp=True, define_retry_function=True))
    if out_gvcf_path:
        batch_instance.write_output(job.output_gvcf, str(out_gvcf_path).replace('.g.vcf.gz', ''))
    return job, job.output_gvcf


def merge_gvcfs_job(
    batch_instance: hb.Batch,
    sequencing_group_name: str,
    gvcf_groups: list[hb.ResourceGroup],
    job_attrs: dict,
    out_gvcf_path: Path,
) -> BashJob:
    """
    Combine by-interval GVCFs into a single sequencing group wide GVCF file.
    """

    job = batch_instance.new_job(f'Merge {len(gvcf_groups)} GVCFs', job_attrs | {'tool': 'picard MergeVcfs'})

    job.image(config.config_retrieve(['images', 'picard']))
    job.cpu(2)
    job.memory('highmem')  # ~ 6G/core ~ 12G
    job.storage(f'{len(gvcf_groups) * 1.5 + 2}G')
    job.declare_resource_group(
        output_gvcf={
            'g.vcf.gz': '{root}-' + sequencing_group_name + '.g.vcf.gz',
            'g.vcf.gz.tbi': '{root}-' + sequencing_group_name + '.g.vcf.gz.tbi',
        },
    )

    input_cmd = ' '.join([f'INPUT={gvcf_group["g.vcf.gz"]} ' for gvcf_group in gvcf_groups])

    job.command(
        f"""
        set -o pipefail
        set -ex

        picard -Xms11g \
        MergeVcfs {input_cmd} OUTPUT={job.output_gvcf['g.vcf.gz']}
        """
    )
    batch_instance.write_output(job.output_gvcf, str(out_gvcf_path).replace('.g.vcf.gz', ''))
    return job


def postprocess_gvcf(
    gvcf_path: filetypes.GvcfPath,
    sequencing_group_name: str,
    job_attrs: dict,
    output_path: Path,
) -> BashJob:
    """
    1. Runs ReblockGVCF to annotate with allele-specific VCF INFO fields
    required for recalibration, and reduce the number of GVCF blocking bins to 2.
    2. Subsets GVCF to main, not-alt chromosomes to avoid downstream errors.
    3. Removes the DS INFO field that is added to some HGDP GVCFs to avoid errors
       from Hail about mismatched INFO annotations
    4. Renames the GVCF sequencing group name to use CPG ID.
    """
    batch_instance = hail_batch.get_batch()

    logging.info(f'Adding GVCF postproc job for sequencing group {sequencing_group_name}, gvcf {gvcf_path}')

    job = batch_instance.new_bash_job('Postproc GVCF', job_attrs | {'tool': 'gatk ReblockGVCF'})
    job.image(config.config_retrieve(['images', 'gatk']))

    # We need at least 2 CPU, so on 16-core instance it would be 8 jobs,
    # meaning we have more than enough disk (265/8=33.125G).
    # Enough to fit a pre-reblocked GVCF, which can be as big as 10G,
    # the reblocked result (1G), and ref data (5G).
    storage_gb = config.config_retrieve(['workflow', 'postproc_gvcf_storage'])
    job_res = resources.STANDARD.set_resources(job, ncpu=2, storage_gb=storage_gb)

    job.declare_resource_group(
        output_gvcf={
            'g.vcf.gz': '{root}.g.vcf.gz',
            'g.vcf.gz.tbi': '{root}.g.vcf.gz.tbi',
        },
    )

    reference = hail_batch.fasta_res_group(batch_instance)
    noalt_regions = batch_instance.read_input(config.config_retrieve(['references', 'noalt_bed']))
    gvcf = batch_instance.read_input(gvcf_path.path)
    gq_bands = config.config_retrieve(['workflow', 'reblock_gq_bands'])

    cmd = f"""\
    GVCF={gvcf}
    GVCF_NODP=$BATCH_TMPDIR/{sequencing_group_name}-nodp.g.vcf.gz
    REBLOCKED=$BATCH_TMPDIR/{sequencing_group_name}-reblocked.g.vcf.gz

    # Reindexing just to make sure the index is not corrupted
    bcftools index --tbi $GVCF

    # Remove INFO/DP field, which contradicts the FORMAT/DP, in the way that
    # it has _all_ reads, not just variant-calling-usable reads. If we keep INFO/DP,
    # ReblockGVCF would prioritize it over FORMAT/DP and change FORMAT/DP to INFO/DP
    # in the resulting merged blocks. It would pick the highest INFO/DP when merging
    # multiple blocks, so a variant in a small homopolymer region (surrounded by
    # long DP=0 areas), that attracted piles of low-MQ reads with INFO/DP=1000
    # will translate into a long GQ<20 block with the same FORMAT/DP=1000,
    # which is wrong, because most of this block has no reads.
    bcftools view $GVCF \\
    | bcftools annotate -x INFO/DP \\
    | bcftools view -Oz -o $GVCF_NODP
    tabix -p vcf $GVCF_NODP

    gatk --java-options "{job_res.java_mem_options()}" \\
    ReblockGVCF \\
    --reference {reference.base} \\
    -V $GVCF_NODP \\
    -do-qual-approx \\
    --floor-blocks {' '.join(f'-GQB {b}' for b in gq_bands)} \\
    -O $REBLOCKED \\
    --create-output-variant-index true

    EXISTING_SN=$(bcftools query -l $GVCF)

    bcftools view $REBLOCKED -T {noalt_regions} \\
    | bcftools annotate -x INFO/DS \\
    | bcftools reheader -s <(echo "$EXISTING_SN {sequencing_group_name}") \\
    | bcftools view -Oz -o {job.output_gvcf['g.vcf.gz']}

    tabix -p vcf {job.output_gvcf['g.vcf.gz']}
    """
    job.command(
        hail_batch.command(
            cmd,
            setup_gcp=True,
            monitor_space=True,
            define_retry_function=True,
        )
    )

    batch_instance.write_output(job.output_gvcf, str(output_path).replace('.g.vcf.gz', ''))

    return job
