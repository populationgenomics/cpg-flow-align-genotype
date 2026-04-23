"""
CRAM to VNtyper results: create Hail Batch jobs to run VNtyper on individual sequencing groups.
"""

import hailtop.batch as hb
from cpg_flow import resources
from cpg_utils import Path, config, hail_batch
from hailtop.batch.job import BashJob

"""
Possible outputs:

results/
├── pipeline_summary.json        # Machine-readable pipeline summary
├── pipeline_summary.csv         # Optional (--summary-formats csv)
├── pipeline_summary.tsv         # Optional (--summary-formats tsv)
├── pipeline.log                 # Pipeline execution log
├── summary_report.html          # HTML report with IGV visualization
├── predefined_regions_<assembly>.bed  # Region BED file (e.g., hg19, hg38)
├── kestrel/
│   ├── kestrel_result.tsv       # Final genotyping result
│   ├── kestrel_pre_result.tsv   # Pre-filter variants (all candidates)
│   ├── output.vcf               # Raw Kestrel VCF
│   ├── output_indel.vcf         # Filtered INDEL VCF
│   ├── output_indel.vcf.gz      # Compressed INDEL VCF (if bcftools available)
│   ├── output.bam               # Kestrel alignment BAM
│   ├── output.bam.bai           # BAM index
│   └── output.bed               # BED file for coverage visualization
├── fastq_bam_processing/
│   ├── output_R1.fastq.gz       # Extracted/processed R1 reads
│   ├── output_R2.fastq.gz       # Extracted/processed R2 reads
│   └── pipeline_info.json       # BAM header metadata (BAM/CRAM input)
├── alignment_processing/
│   └── output_sorted.bam        # BWA-aligned BAM (FASTQ input only)
├── coverage/
│   └── coverage_summary.tsv     # VNTR region coverage statistics
└── advntr/                      # Only when --extra-modules advntr
    ├── output_adVNTR.tsv         # Raw adVNTR output
    ├── output_adVNTR_result.tsv  # Processed adVNTR result
    └── cross_match_results.tsv   # Kestrel vs adVNTR comparison
"""


def vntyper(
    sequencing_group_name: str,
    cram_path: str,
    out_paths: dict[str, Path],
    job_attrs: dict[str, str],
) -> tuple[BashJob | None, hb.Resource]:
    """
    Add one VNtyper job.
    """

    batch_instance = hail_batch.get_batch()

    sequencing_type = config.config_retrieve(['workflow', 'sequencing_type'])
    job = batch_instance.new_bash_job(f'VNtyper {sequencing_group_name}', job_attrs | {'tool': 'vntyper'})
    job.image(config.config_retrieve(['images', 'vntyper']))

    storage_default = 40 if sequencing_type == 'genome' else None

    job_res = resources.HIGHMEM.request_resources(
        ncpu=config.config_retrieve(['workflow', 'vntyper_cpu'], 4),
        storage_gb=config.config_retrieve(['workflow', 'vntyper_storage'], storage_default),
    )
    job_res.set_to_job(job)

    reference = hail_batch.fasta_res_group(batch_instance)

    cram_resource_group = batch_instance.read_input_group(**{'cram': cram_path, 'cram.crai': f'{cram_path!s}.crai'})

    job.command(f"""\
    export REF_PATH={reference['fasta']} && \

    vntyper pipeline \
        --cram {cram_resource_group.cram} \
        --reference-assembly hg38 \
        -o ./results \
        --threads 4
    """)

    job.command(f"""\
        mv ./results/summary_report.html {job.html} && \
        mv ./results/kestrel/kestrel_result.tsv {job.kestrel}
    """)

    batch_instance.write_output(job.html, str(out_paths['html']))
    batch_instance.write_output(job.kestrel, str(out_paths['kestrel']))

    return job, job.output_files
