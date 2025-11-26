from cpg_flow import stage, targets
from cpg_utils import Path, config

from align_genotype.jobs import multiqc, somalier
from align_genotype.stages import (
    CramQcPicardCollectMetrics,
    CramQcPicardMultiMetrics,
    CramQcSamtoolsStats,
    CramQcSomalier,
    CramQcVerifyBamId,
    RunGvcfQc,
)


def filter_to_dataset_sgids(
    inputs_by_sgid: dict,
    dataset: targets.Dataset,
):
    """Filter a dict of inputs by sequencing group ID to only those in the dataset."""
    return {k: v for k, v in inputs_by_sgid.items() if k in dataset.get_sequencing_group_ids()}


@stage.stage(required_stages=[CramQcVerifyBamId, CramQcSomalier])
class SomalierPedigree(stage.DatasetStage):
    """
    Checks pedigree from CRAM fingerprints.
    """

    def expected_outputs(self, dataset: targets.Dataset) -> dict[str, Path]:
        """
        Return the report for MultiQC, plus putting an HTML into the web bucket.
        MultiQC expects the following patterns:
        * *.samples.tsv
        * *.pairs.tsv
        https://github.com/ewels/MultiQC/blob/master/multiqc/utils/search_patterns.yaml#L472-L481
        """

        prefix = dataset.prefix() / 'somalier' / 'cram' / dataset.get_alignment_inputs_hash()
        web_prefix = dataset.web_prefix() / 'somalier' / 'cram' / dataset.get_alignment_inputs_hash()
        return {
            'samples': prefix / f'{dataset.name}.samples.tsv',
            'expected_ped': prefix / f'{dataset.name}.expected.ped',
            'pairs': prefix / f'{dataset.name}.pairs.tsv',
            'html': web_prefix / 'cram-somalier-pedigree.html',
            'checks': prefix / f'{dataset.name}-checks.done',
        }

    def queue_jobs(self, dataset: targets.Dataset, inputs: stage.StageInput) -> stage.StageOutput:
        """
        Checks calls job from the pedigree module
        """

        outputs = self.expected_outputs(dataset)

        verifybamid_by_sgid = filter_to_dataset_sgids(inputs.as_path_by_target(CramQcVerifyBamId), dataset)
        somalier_path_by_sgid = filter_to_dataset_sgids(inputs.as_path_by_target(CramQcSomalier), dataset)

        html_url = str(outputs['html']).replace(str(dataset.web_prefix()), dataset.web_url())

        jobs = somalier.pedigree(
            dataset=dataset,
            somalier_path_by_sgid=somalier_path_by_sgid,
            verifybamid_by_sgid=verifybamid_by_sgid,
            outputs=outputs,
            out_html_url=html_url,
            label='Somalier',
            job_attrs=self.get_job_attrs(dataset),
        )
        return self.make_outputs(dataset, data=outputs, jobs=jobs)


@stage.stage(
    required_stages=[
        CramQcVerifyBamId,
        CramQcSomalier,
        CramQcSamtoolsStats,
        CramQcPicardCollectMetrics,
        CramQcPicardMultiMetrics,
        SomalierPedigree,
    ],
    analysis_type='qc',
    analysis_keys=['json'],
)
class CramMultiQC(stage.DatasetStage):
    """
    Run MultiQC to aggregate CRAM QC stats.
    """

    def expected_outputs(self, dataset: targets.Dataset) -> dict[str, Path]:
        """
        Expected to produce an HTML and a corresponding JSON file.
        """

        # get the unique hash for these Sequencing Groups
        sg_hash = dataset.get_alignment_inputs_hash()
        return {
            'html': dataset.web_prefix() / 'qc' / 'cram' / sg_hash / 'multiqc.html',
            'json': dataset.prefix() / 'qc' / 'cram' / sg_hash / 'multiqc_data.json',
            'checks': dataset.prefix() / 'qc' / 'cram' / sg_hash / '.checks',
        }

    def queue_jobs(self, dataset: targets.Dataset, inputs: stage.StageInput) -> stage.StageOutput | None:
        """
        Call a function from the `jobs` module using inputs from `cramqc` and `somalier` stages.
        """
        outputs = self.expected_outputs(dataset)

        paths = [str(each_path) for each_path in inputs.as_dict(dataset, SomalierPedigree).values()]

        ending_to_trim = {
            '.alignment_summary_metrics',
            '.base_distribution_by_cycle_metrics',
            '.insert_size_metrics',
            '.quality_by_cycle_metrics',
            '.quality_yield_metrics',
            '.verify-bamid.selfSM',
            '.samtools-stats',
        }
        modules_to_trim_endings = {
            'picard/alignment_metrics',
            'picard/insertsize',
            'picard/basedistributionbycycle',
            'picard/quality_by_cycle',
            'picard/quality_yield_metrics',
            'verifybamid/selfsm',
            'samtools/stats',
        }

        # add in alignemnt metrics
        for content in filter_to_dataset_sgids(inputs.as_dict_by_target(CramQcPicardMultiMetrics), dataset).values():
            paths.extend(
                [
                    content['summary'],
                    content['base_dist'],
                    content['insert_size'],
                    content['qual_by_cycle'],
                    content['yield'],
                ]
            )
        paths.extend(filter_to_dataset_sgids(inputs.as_path_by_target(CramQcSomalier), dataset).values())
        paths.extend(filter_to_dataset_sgids(inputs.as_path_by_target(CramQcVerifyBamId), dataset).values())
        paths.extend(filter_to_dataset_sgids(inputs.as_path_by_target(CramQcSamtoolsStats), dataset).values())
        paths.extend(filter_to_dataset_sgids(inputs.as_path_by_target(CramQcPicardCollectMetrics), dataset).values())

        if base_url := dataset.web_url():
            html_url = str(outputs['html']).replace(str(dataset.web_prefix()), base_url)
        else:
            html_url = None

        sequencing_type = config.config_retrieve(['workflow', 'sequencing_type'])
        if sequencing_type == 'genome':
            ending_to_trim.add('.picard-wgs-metrics')
            modules_to_trim_endings.add('picard/wgs_metrics')
        if sequencing_type == 'exome':
            ending_to_trim.add('.picard-hs-metrics')
            modules_to_trim_endings.add('picard/hsmetrics')

        send_to_slack = config.config_retrieve(['workflow', 'gvcf_multiqc', 'send_to_slack'], default=True)
        extra_config = config.config_retrieve(['workflow', 'gvcf_multiqc', 'extra_config'], default={})
        extra_config['table_columns_visible'] = {'FastQC': False}

        jobs = multiqc.multiqc(
            tmp_prefix=dataset.tmp_prefix() / 'multiqc' / 'cram',
            paths=paths,
            ending_to_trim=ending_to_trim,
            modules_to_trim_endings=modules_to_trim_endings,
            dataset=dataset,
            outputs=outputs,
            out_checks_path=outputs['checks'],
            out_html_url=html_url,
            job_attrs=self.get_job_attrs(dataset),
            sequencing_group_id_map=dataset.rich_id_map(),
            label='CRAM',
            extra_config=extra_config,
            send_to_slack=send_to_slack,
        )
        return self.make_outputs(dataset, data=outputs, jobs=jobs)


def _update_meta(output_path: str) -> dict:
    import json  # noqa: PLC0415

    from cloudpathlib import CloudPath  # noqa: PLC0415

    with CloudPath(output_path).open() as f:
        d = json.load(f)
    return {'multiqc': d['report_general_stats_data']}


@stage.stage(
    required_stages=[RunGvcfQc],
    analysis_type='qc',
    analysis_keys=['json'],
    update_analysis_meta=_update_meta,
)
class GvcfMultiQC(stage.DatasetStage):
    """Run MultiQC to summarise all GVCF QC."""

    def expected_outputs(self, dataset: stage.Dataset) -> dict[str, Path]:
        """Expected to produce an HTML and a corresponding JSON file."""
        sg_hash = dataset.get_alignment_inputs_hash()
        return {
            'html': dataset.web_prefix() / 'qc' / 'gvcf' / sg_hash / 'multiqc.html',
            'json': dataset.prefix() / 'qc' / 'gvcf' / sg_hash / 'multiqc_data.json',
            'checks': dataset.prefix() / 'qc' / 'gvcf' / sg_hash / '.checks',
        }

    def queue_jobs(self, dataset: targets.Dataset, inputs: stage.StageInput) -> stage.StageOutput:
        """Collect gVCF QC."""
        outputs = self.expected_outputs(dataset)

        paths = [
            content['detail']
            for content in filter_to_dataset_sgids(inputs.as_dict_by_target(RunGvcfQc), dataset).values()
        ]

        if base_url := dataset.web_url():
            html_url = str(outputs['html']).replace(str(dataset.web_prefix()), base_url)
        else:
            html_url = None

        send_to_slack = config.config_retrieve(['workflow', 'cram_multiqc', 'send_to_slack'], default=True)
        extra_config = config.config_retrieve(['workflow', 'cram_multiqc', 'extra_config'], default={})
        extra_config['table_columns_visible'] = {'Picard': True}

        jobs = multiqc.multiqc(
            dataset=dataset,
            outputs=outputs,
            out_html_url=html_url,
            out_checks_path=outputs['checks'],
            tmp_prefix=dataset.tmp_prefix() / 'multiqc' / 'gvcf',
            paths=paths,
            ending_to_trim={'.variant_calling_detail_metrics'},
            modules_to_trim_endings={'picard/variant_calling_metrics'},
            job_attrs=self.get_job_attrs(dataset),
            sequencing_group_id_map=dataset.rich_id_map(),
            label='GVCF',
            extra_config=extra_config,
            send_to_slack=send_to_slack,
        )
        return self.make_outputs(dataset, data=outputs, jobs=jobs)
