# Alignment and Genotyping

A [CPG-Flow](https://github.com/populationgenomics/cpg-flow) migration of the Alignment and Genotyping workflow from Production Pipelines.

These workflows start with the single-sample assay data we've received from our collaborators, FastQ, BAM, or FQ.ora, and carry out the alignment and genotyping steps of the analysis, including:

- Aligning the reads to the reference genome using DRAGMAP (Dragen-OS)
- Generating alignment quality metrics using Samtools, Picard, and Somalier fingerprinting
- Genotyping the aligned reads using GATK HaplotypeCaller
- Generating gVCF quality metrics using Picard
- Running VerifyBamID to check for contamination and sample swaps
- Running VNtyper to genotype the VNTR regions of the genome

## QC stages

QC metrics are collected at the single sample level, and then aggregated at the Dataset level to produce a MultiQC report. The QC metrics are compared against thresholds defined in the configuration file, and any failures or warnings are reported to Slack.

### QC thresholds - genome

| Metric Key | Source | Description | Warn | Fail | Justification |
| --- | --- | --- | --- | --- | --- |
| `MEDIAN_COVERAGE` | Picard `CollectWgsMetrics` | Median Coverage | < 25x | < 15x | Fails below 15x catches the under-sequenced tail (~0-1% of every cohort) without tripping on normal depth variation. |
| `PCT_20X` | Picard `CollectWgsMetrics` | % Bases >= 20x | < 85% | < 75% | Typical cohort medians ~95%; fail below 75% flags seriously inadequate breadth. |
| `reads_mapped_percent` | `samtools stats` | % Reads mapped | < 97% | < 80% | Typical cohort medians ~98%; fail below 80% flags disastrous mapping issues. |
| `reads_properly_paired_percent` | `samtools stats` | % Reads properly paired | < 92% | < 90% | Typical cohort medians ~96-98%; fail below 90% flags potential issues with library preparation. |
| `reads_duplicated_percent` | `samtools stats` | % Reads duplicated | > 30% | > 40% | Typical cohort medians ~5-10%; fail above 40% flags serious issues with the sequencing. |
| `FREEMIX` | VerifyBamID2 | Contamination / Freemix | > 2% | > 4% | Typical cohort medians ~1%; fail above 4% flags potential contamination. |

### QC thresholds - exome

Exome thresholds. Whole-genome coverage metrics from `CollectWgsMetrics` are meaningless for exomes, so we gate on Picard `CollectHsMetrics` target-coverage metrics instead.

Note that the cohort medians for these metrics are highly dependent on the capture kit used and vary greatly between datasets.

| Metric Key | Source | Description | Warn | Fail | Justification |
| --- | --- | --- | --- | --- | --- |
| `MEAN_TARGET_COVERAGE` | Picard `CollectHsMetrics` | Mean coverage over target region | < 50x | < 30x | Typical cohort medians ~100x; fail below 30x flags inadequate coverage. |
| `PCT_TARGET_BASES_20X` | Picard `CollectHsMetrics` | % Target bases >= 20x | < 90% | < 80% | Typical cohort medians ~90-95%; fail below 75% flags inadequate breadth. |
| `FOLD_80_BASE_PENALTY` | Picard `CollectHsMetrics` | Fold 80 base penalty | > 2.0 | > 3.0 | Typical cohort medians ~1.8; fail above 3.0 flags inadequate uniformity. |
| `ZERO_CVG_TARGETS_PCT` | Picard `CollectHsMetrics` | % Target bases with zero coverage | > 7% | > 10% | Typical cohort medians ~0-1%; fail above 10% flags inadequate coverage. |
| `reads_mapped_percent` | `samtools stats` | % Reads mapped | < 95% | < 80% | Typical cohort medians ~98%; fail below 80% flags disastrous mapping issues. |
| `reads_duplicated_percent` | `samtools stats` | % Reads duplicated | > 30% | > 50% | Typical cohort medians ~5-10%; fail above 50% flags serious issues with the sequencing. |
| `FREEMIX` | VerifyBamID2 | Contamination / Freemix | > 2% | > 4% | Typical cohort medians ~1%; fail above 4% flags potential contamination. |

## Running the workflows

This single-sample workflow has a dedicated entrypoint, and can be operated through analysis runner as follows:

```bash
analysis-runner --skip-repo-checkout \
    --image australia-southeast1-docker.pkg.dev/cpg-common/images/cpg-flow-align-genotype:0.5.4-1 \
    --config src/align_genotype/config_template.toml \
    --dataset seqr \
    --description 'Single-Sample data generation' \
    --access-level full \
    --output-dir OUTPUT_DIR \
    run_workflow
```

A secondary workflow continues on from the single-sample steps, and produces Dataset-level outputs, including:

- Runs Somalier Relate on the Dataset's Somalier fingerprints to generate a Dataset-level pedigree
- Uses the somalier outputs to check relationships against the expected pedigree
- Runs Dataset-level MultiQC for the CRAM and gVCF metrics, publishing an HTML report, and writing the results to Slack

This Dataset-level workflow can be run in a similar way, but with a different entrypoint:

```bash
analysis-runner --skip-repo-checkout \
    --image australia-southeast1-docker.pkg.dev/cpg-common/images/cpg-flow-align-genotype:0.5.4-1 \
    --config src/align_genotype/config_template.toml \
    --dataset seqr \
    --description 'Dataset-Level QC workflow' \
    --access-level full \
    --output-dir OUTPUT_DIR \
    dataset_workflow
```

A third workflow is available to run VNtyper on a set of samples, which can be used to genotype the VNTR regions of the genome. This workflow requires the Align stage to have completed and can be run with the following command:

```bash
analysis-runner --skip-repo-checkout \
    --image australia-southeast1-docker.pkg.dev/cpg-common/images/cpg-flow-align-genotype:0.5.4-1 \
    --config src/align_genotype/config_template.toml \
    --config src/align_genotype/vntyper.toml \
    --dataset seqr \
    --description 'VNtyper workflow' \
    --access-level full \
    --output-dir OUTPUT_DIR \
    run_workflow
```

## Configuration

A [configuration template file](src/align_genotype/config_template.toml) is provided, which contains all the settings and references required to run the workflow. This file can be copied and modified to create a specific configuration for a given run. This config file should contain all values required by the workflow, meaning there is no reliance on the default global configuration generated by `Analysis-Runner`, except for relating to storage locations.

Two entries in the config template should be modified for each run:

- `workflow.input_cohorts`:  a list of Cohort IDs to be used as input
- `workflow.sequencing_type`: 'exome' or 'genome', depending on the type of sequencing data being processed

A secondary config file is provided for the VNtyper workflow, which contains settings specific to only this stage. Without this config, the stage is always skipped. It should be last in the config list when running the workflow, so that it can override any relevant settings in the main config template.

## Structure

This repository has the following structure:

```text
src
├── align_genotype
│   ├── __init__.py
│   ├── config_template.toml
│   ├── dataset_stages.py
│   ├── dataset_workflow.py
│   ├── jobs
│   │   ├── __init__.py
│   │   ├── align.py
│   │   ├── cram_qc_samtools.py
│   │   ├── cram_qc_somalier.py
│   │   ├── cram_qc_verify.py
│   │   ├── genotype.py
│   │   ├── multiqc.py
│   │   ├── picard.py
│   │   └── somalier.py
│   ├── run_workflow.py
│   ├── scripts
│   │   ├── __init__.py
│   │   └── check_pedigree.py
│   ├── stages.py
│   └── utils.py
```

## Original code

The original code for the alignment and genotyping workflow can be found in the Production Pipelines repository, in the following files:

- [Alignment](https://github.com/populationgenomics/production-pipelines/blob/ca8741c9d34c85f3f3e0811f081e67d56086d831/cpg_workflows/stages/align.py)
- [Genotyping](https://github.com/populationgenomics/production-pipelines/blob/ca8741c9d34c85f3f3e0811f081e67d56086d831/cpg_workflows/stages/genotype.py)
- [Cram QC](https://github.com/populationgenomics/production-pipelines/blob/main/cpg_workflows/stages/cram_qc.py)
- [gVCF QC](https://github.com/populationgenomics/production-pipelines/blob/main/cpg_workflows/stages/gvcf_qc.py)
