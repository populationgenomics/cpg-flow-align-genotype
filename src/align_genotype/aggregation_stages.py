from cpg_flow import stage, targets
from cpg_utils import Path, config

from align_genotype.jobs import somalier, multiqc
from align_genotype.stages import CramQcVerifyBamId, CramQcSomalier


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
        https://github.com/ewels/MultiQC/blob/master/multiqc/utils/search_patterns
        .yaml#L472-L481
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

        verifybamid_by_sgid = inputs.as_path_by_target(CramQcVerifyBamId)
        somalier_path_by_sgid = inputs.as_path_by_target(CramQcSomalier)

        html_url = str(outputs['html']).replace(str(dataset.web_prefix()), dataset.web_url())

        if any(sg.pedigree.dad or sg.pedigree.mom for sg in dataset.get_sequencing_groups()):
            jobs = somalier.pedigree(
                dataset=dataset,
                somalier_path_by_sgid=somalier_path_by_sgid,
                verifybamid_by_sgid=verifybamid_by_sgid,
                outputs=outputs,
                out_html_url=html_url,
                job_attrs=self.get_job_attrs(dataset),
            )
            return self.make_outputs(dataset, data=outputs, jobs=jobs)
        else:
            return self.make_outputs(dataset, skipped=True)


# @stage.stage(required_stages=[CramQC, SomalierPedigree], analysis_type='qc', analysis_keys=['json'])
# class CramMultiQC(stage.DatasetStage):
#     """
#     Run MultiQC to aggregate CRAM QC stats.
#     """
#
#     def expected_outputs(self, dataset: targets.Dataset) -> dict[str, Path]:
#         """
#         Expected to produce an HTML and a corresponding JSON file.
#         """
#
#         # get the unique hash for these Sequencing Groups
#         sg_hash = dataset.get_alignment_inputs_hash()
#         return {
#             'html': dataset.web_prefix() / 'qc' / 'cram' / sg_hash / 'multiqc.html',
#             'json': dataset.prefix() / 'qc' / 'cram' / sg_hash / 'multiqc_data.json',
#             'checks': dataset.prefix() / 'qc' / 'cram' / sg_hash / '.checks',
#         }
#
#     def queue_jobs(self, dataset: targets.Dataset, inputs: stage.StageInput) -> stage.StageOutput | None:
#         """
#         Call a function from the `jobs` module using inputs from `cramqc`
#         and `somalier` stages.
#         """
#         outputs = self.expected_outputs(dataset)
#         html_path = outputs['html']
#         checks_path = outputs['checks']
#         html_url = str(html_path).replace(str(dataset.web_prefix()), dataset.web_url())
#
#         paths = []
#         try:
#             somalier_samples = inputs.as_path(dataset, SomalierPedigree, key='samples')
#             somalier_pairs = inputs.as_path(dataset, SomalierPedigree, key='pairs')
#         except StageInputNotFoundError:
#             pass
#         else:
#             paths = [
#                 somalier_samples,
#                 somalier_pairs,
#             ]
#
#         ending_to_trim = set()  # endings to trim to get sample names
#         modules_to_trim_endings = set()
#
#         for sequencing_group in dataset.get_sequencing_groups():
#             for qc in qc_functions():
#                 for key, out in qc.outs.items():
#                     if not out:
#                         continue
#                     try:
#                         path = inputs.as_path(sequencing_group, CramQC, key)
#                     except StageInputNotFoundError:  # allow missing inputs
#                         logging.warning(
#                             f'Output CramQc/"{key}" not found for {sequencing_group}, '
#                             f'it will be silently excluded from MultiQC',
#                         )
#                         continue
#                     modules_to_trim_endings.add(out.multiqc_key)
#                     paths.append(path)
#                     ending_to_trim.add(path.name.replace(sequencing_group.id, ''))
#
#         if not paths:
#             logging.warning('No CRAM QC found to aggregate with MultiQC')
#             return self.make_outputs(dataset)
#
#         send_to_slack = config_retrieve(['workflow', 'cram_multiqc', 'send_to_slack'], default=True)
#         extra_config = config_retrieve(['workflow', 'cram_multiqc', 'extra_config'], default={})
#         extra_config['table_columns_visible'] = {'FastQC': False}
#
#         jobs = multiqc.multiqc(
#             tmp_prefix=dataset.tmp_prefix() / 'multiqc' / 'cram',
#             paths=paths,
#             ending_to_trim=ending_to_trim,
#             modules_to_trim_endings=modules_to_trim_endings,
#             dataset=dataset,
#             out_json_path=outputs['json'],
#             out_html_path=html_path,
#             out_html_url=html_url,
#             out_checks_path=checks_path,
#             job_attrs=self.get_job_attrs(dataset),
#             sequencing_group_id_map=dataset.rich_id_map(),
#             label='CRAM',
#             send_to_slack=send_to_slack,
#             extra_config=extra_config,
#         )
#         return self.make_outputs(dataset, data=outputs, jobs=jobs)
#
#
# @stage(
#     required_stages=[
#         GvcfQC,
#         GvcfHappy,
#     ],
#     analysis_type='qc',
#     analysis_keys=['json'],
#     update_analysis_meta=_update_meta,
# )
# class GvcfMultiQC(DatasetStage):
#     """
#     Run MultiQC to summarise all GVCF QC.
#     """
#
#     def expected_outputs(self, dataset: Dataset) -> dict[str, Path]:
#         """
#         Expected to produce an HTML and a corresponding JSON file.
#         """
#         if config_retrieve(['workflow', 'skip_qc'], False):
#             return {}
#
#         # get the unique hash for these Sequencing Groups
#         sg_hash = dataset.get_alignment_inputs_hash()
#         return {
#             'html': dataset.web_prefix() / 'qc' / 'gvcf' / sg_hash / 'multiqc.html',
#             'json': dataset.prefix() / 'qc' / 'gvcf' / sg_hash / 'multiqc_data.json',
#             'checks': dataset.prefix() / 'qc' / 'gvcf' / sg_hash / '.checks',
#         }
#
#     def queue_jobs(self, dataset: Dataset, inputs: StageInput) -> StageOutput | None:
#         """
#         Collect QC.
#         """
#         outputs = self.expected_outputs(dataset)
#         json_path = outputs['json']
#         html_path = outputs['html']
#         checks_path = outputs['checks']
#         if base_url := dataset.web_url():
#             html_url = str(html_path).replace(str(dataset.web_prefix()), base_url)
#         else:
#             html_url = None
#
#         paths = []
#         ending_to_trim = set()  # endings to trim to get sample names
#
#         for sequencing_group in dataset.get_sequencing_groups():
#             for _stage, key in [
#                 (GvcfQC, 'qc_detail'),
#                 (GvcfHappy, None),
#             ]:
#                 try:
#                     path = inputs.as_path(sequencing_group, _stage, key)
#                 except StageInputNotFoundError:  # allow missing inputs
#                     if _stage != GvcfHappy:
#                         logging.warning(
#                             f'Output {_stage.__name__}/"{key}" not found for {sequencing_group}, '
#                             f'it will be silently excluded from MultiQC',
#                         )
#                 else:
#                     paths.append(path)
#                     ending_to_trim.add(path.name.replace(sequencing_group.id, ''))
#
#         if not paths:
#             logging.warning('No GVCF QC found to aggregate with MultiQC')
#             return self.make_outputs(dataset)
#
#         modules_to_trim_endings = {
#             'picard/variant_calling_metrics',
#             'happy',
#         }
#
#         send_to_slack = config_retrieve(['workflow', 'gvcf_multiqc', 'send_to_slack'], default=True)
#         extra_config = config_retrieve(['workflow', 'gvcf_multiqc', 'extra_config'], default={})
#         extra_config['table_columns_visible'] = {'Picard': True}
#
#         jobs = multiqc(
#             get_batch(),
#             tmp_prefix=dataset.tmp_prefix() / 'multiqc' / 'gvcf',
#             paths=paths,
#             ending_to_trim=ending_to_trim,
#             modules_to_trim_endings=modules_to_trim_endings,
#             dataset=dataset,
#             out_json_path=json_path,
#             out_html_path=html_path,
#             out_html_url=html_url,
#             out_checks_path=checks_path,
#             job_attrs=self.get_job_attrs(dataset),
#             sequencing_group_id_map=dataset.rich_id_map(),
#             extra_config=extra_config,
#             send_to_slack=send_to_slack,
#             label='GVCF',
#         )
#         return self.make_outputs(dataset, data=outputs, jobs=jobs)
