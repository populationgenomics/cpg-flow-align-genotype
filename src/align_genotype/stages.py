"""
Alignment and genotyping... oof. Where to start.
"""

from cpg_flow import stage, targets, workflow
from cpg_utils import Path, config, to_path

from align_genotype.jobs.align import align
from align_genotype.jobs.cram_qc_samtools import samtools_stats
from align_genotype.jobs.cram_qc_somalier import extract_somalier
from align_genotype.jobs.cram_qc_verify import verifybamid
from align_genotype.jobs.genotype import genotype
from align_genotype.jobs.picard import collect_metrics, generate_intervals, hs_metrics, vcf_qc, wgs_metrics


@stage.stage
class GenerateIntervalsOnce(stage.MultiCohortStage):
    def expected_outputs(self, multicohort: targets.MultiCohort) -> dict[str, list[Path]]:
        scatter_count = config.config_retrieve(['workflow', 'scatter_count_genotype'])
        prefix = multicohort.analysis_dataset.prefix(category='tmp') / 'genotype' / 'intervals'
        return {'intervals': [prefix / to_path(f'{idx}.interval_list') for idx in range(1, scatter_count + 1)]}

    def queue_jobs(
        self,
        multicohort: targets.MultiCohort,
        inputs: stage.StageInput,  # noqa: ARG002
    ) -> stage.StageOutput:
        outputs = self.expected_outputs(multicohort)
        job = generate_intervals(outputs['intervals'], self.get_job_attrs(multicohort))
        return self.make_outputs(multicohort, jobs=job, data=outputs)


@stage.stage(
    analysis_type='cram',
    analysis_keys=['cram'],
)
class AlignWithDragmap(stage.SequencingGroupStage):
    """
    This is a generic stage that runs a bash command.
    """

    def expected_outputs(self, sequencing_group: targets.SequencingGroup) -> dict[str, Path | str]:
        """
        n.b. markduplicates-metrics are a String here so their existence isn't detected during DAG assembly.
        This means we will not accidentally restart alignment where this optional accessory file doesn't exist.
        """
        cram_path = sequencing_group.cram if sequencing_group.cram else sequencing_group.make_cram_path()
        return {
            'cram': cram_path.path,
            'sorted_bam': str(
                sequencing_group.dataset.prefix(category='tmp') / 'align' / f'{sequencing_group.id}.sorted.bam'
            ),
            'markdup': str(
                sequencing_group.dataset.prefix()
                / 'qc'
                / 'markduplicates_metrics'
                / f'{sequencing_group.id}.markduplicates-metrics'
            ),
        }

    def queue_jobs(self, sequencing_group: targets.SequencingGroup, inputs: stage.StageInput) -> stage.StageOutput:  # noqa: ARG002
        """
        This is where we generate jobs for this stage.
        """

        outputs = self.expected_outputs(sequencing_group)

        jobs = align(
            sequencing_group=sequencing_group,
            job_attrs=self.get_job_attrs(sequencing_group),
            output_path=outputs['cram'],
            sorted_bam_path=outputs['sorted_bam'],
            markdup_metrics_path=outputs['markdup'],
        )

        # return the jobs and outputs
        return self.make_outputs(target=sequencing_group, data=outputs, jobs=jobs)


@stage.stage(
    required_stages=[AlignWithDragmap, GenerateIntervalsOnce],
    analysis_type='gvcf',
)
class GenotypeWithGatk(stage.SequencingGroupStage):
    """
    Use HaplotypeCaller to genotype individual sequencing groups (i.e. CRAM -> GVCF).
    """

    def expected_outputs(self, sequencing_group: targets.SequencingGroup) -> Path:
        """
        Generate a GVCF and corresponding TBI index.
        """
        gvcf = sequencing_group.gvcf or sequencing_group.make_gvcf_path()
        return gvcf.path

    def queue_jobs(self, sequencing_group: targets.SequencingGroup, inputs: stage.StageInput) -> stage.StageOutput:
        """
        Use function from the jobs module
        """

        output = self.expected_outputs(sequencing_group)

        cram = inputs.as_str(sequencing_group, AlignWithDragmap, 'cram')
        intervals = inputs.as_dict(workflow.get_multicohort(), GenerateIntervalsOnce)['intervals']

        jobs = genotype(
            output_path=output,
            sequencing_group_name=sequencing_group.id,
            cram_path=cram,
            tmp_prefix=self.tmp_prefix / sequencing_group.id,
            intervals=intervals,
            job_attrs=self.get_job_attrs(sequencing_group),
        )
        return self.make_outputs(sequencing_group, data=output, jobs=jobs)


@stage.stage(required_stages=[AlignWithDragmap])
class CramQcSomalier(stage.SequencingGroupStage):
    """Run somalier extract on a CRAM file."""

    def expected_outputs(self, sequencing_group: targets.SequencingGroup) -> Path:
        return sequencing_group.make_cram_path().somalier_path

    def queue_jobs(self, sequencing_group: targets.SequencingGroup, inputs: stage.StageInput) -> stage.StageOutput:
        output = self.expected_outputs(sequencing_group)

        cram = inputs.as_str(sequencing_group, AlignWithDragmap, 'cram')

        jobs = extract_somalier(
            cram_path=cram,
            output=output,
            job_attrs=self.get_job_attrs(sequencing_group),
        )

        return self.make_outputs(sequencing_group, data=output, jobs=jobs)


@stage.stage(required_stages=[AlignWithDragmap])
class CramQcVerifyBamId(stage.SequencingGroupStage):
    """Run verifyBamId on a CRAM file."""

    def expected_outputs(self, sequencing_group: targets.SequencingGroup) -> Path:
        return sequencing_group.dataset.prefix() / 'qc' / 'verify_bamid' / f'{sequencing_group.id}.verify-bamid.selfSM'

    def queue_jobs(self, sequencing_group: targets.SequencingGroup, inputs: stage.StageInput) -> stage.StageOutput:
        output = self.expected_outputs(sequencing_group)

        cram = inputs.as_str(sequencing_group, AlignWithDragmap, 'cram')

        jobs = verifybamid(
            cram_path=cram,
            output=output,
            job_attrs=self.get_job_attrs(sequencing_group),
        )

        return self.make_outputs(sequencing_group, data=output, jobs=jobs)


@stage.stage(required_stages=[AlignWithDragmap])
class CramQcSamtoolsStats(stage.SequencingGroupStage):
    """Run Samtools  on a CRAM file."""

    def expected_outputs(self, sequencing_group: targets.SequencingGroup) -> Path:
        return sequencing_group.dataset.prefix() / 'qc' / 'samtools_stats' / f'{sequencing_group.id}.samtools-stats'

    def queue_jobs(self, sequencing_group: targets.SequencingGroup, inputs: stage.StageInput) -> stage.StageOutput:
        output = self.expected_outputs(sequencing_group)

        cram = inputs.as_str(sequencing_group, AlignWithDragmap, 'cram')

        jobs = samtools_stats(
            cram_path=cram,
            output=output,
            job_attrs=self.get_job_attrs(sequencing_group),
        )

        return self.make_outputs(sequencing_group, data=output, jobs=jobs)


@stage.stage(required_stages=[AlignWithDragmap])
class CramQcPicardMultiMetrics(stage.SequencingGroupStage):
    """
    Run Picard CollectMetrics  on a CRAM file.

    Severe technical debt warning - in production-pipelines this one Stage ran to completion, and wrote each
    individual output to a different folder, instead of segregating the results by Stage.

    To prevent re-running this Stage on old data, incurring substantial costs to remove data from cold storage, that
    same output convention is used here. i.e., this stage will write 4 outputs into 4 different folders, all under the
    qc top-level.
    """

    def expected_outputs(self, sequencing_group: targets.SequencingGroup) -> dict[str, Path]:
        qc_prefix = sequencing_group.dataset.prefix() / 'qc'
        return {
            'summary': qc_prefix / 'alignment_summary_metrics' / f'{sequencing_group.id}.alignment_summary_metrics',
            'base_dist': qc_prefix
            / 'base_distribution_by_cycle_metrics'
            / f'{sequencing_group.id}.base_distribution_by_cycle_metrics',
            'insert_size': qc_prefix / 'insert_size_metrics' / f'{sequencing_group.id}.insert_size_metrics',
            'qual_by_cycle': qc_prefix / 'quality_by_cycle_metrics' / f'{sequencing_group.id}.quality_by_cycle_metrics',
            'yield': qc_prefix / 'quality_yield_metrics' / f'{sequencing_group.id}.quality_yield_metrics',
        }

    def queue_jobs(self, sequencing_group: targets.SequencingGroup, inputs: stage.StageInput) -> stage.StageOutput:
        outputs = self.expected_outputs(sequencing_group)

        cram = inputs.as_str(sequencing_group, AlignWithDragmap, 'cram')

        jobs = collect_metrics(
            cram_path=cram,
            outputs=outputs,
            job_attrs=self.get_job_attrs(sequencing_group),
        )

        return self.make_outputs(sequencing_group, data=outputs, jobs=jobs)


@stage.stage(required_stages=[AlignWithDragmap])
class CramQcPicardCollectMetrics(stage.SequencingGroupStage):
    """
    Run Picard CollectMetrics on a CRAM file.

    Severe technical debt warning - in production-pipelines this one Stage ran to completion, and wrote each
    individual output to a different folder, instead of segregating the results by Stage.
    To prevent re-running this Stage on old data, incurring substantial costs to remove data from cold storage, that
    same output convention is used here.
    """

    def expected_outputs(self, sequencing_group: targets.SequencingGroup) -> Path:
        qc_prefix = sequencing_group.dataset.prefix() / 'qc'
        # genomes use wgs, exome uses hs - this is in the path and file extension
        qc_type = 'wgs' if sequencing_group.sequencing_type == 'genome' else 'hs'
        return qc_prefix / f'picard_{qc_type}_metrics/{sequencing_group.id}.picard-{qc_type}-metrics'

    def queue_jobs(self, sequencing_group: targets.SequencingGroup, inputs: stage.StageInput) -> stage.StageOutput:
        output = self.expected_outputs(sequencing_group)

        cram = inputs.as_str(sequencing_group, AlignWithDragmap, 'cram')

        if sequencing_group.sequencing_type == 'genome':
            job = wgs_metrics(
                cram_path=cram,
                output=output,
                job_attrs=self.get_job_attrs(sequencing_group),
            )
        elif sequencing_group.sequencing_type == 'exome':
            job = hs_metrics(
                cram_path=cram,
                output=output,
                job_attrs=self.get_job_attrs(sequencing_group),
            )

        else:
            raise ValueError(f'Unsupported sequencing type: {sequencing_group.sequencing_type}')

        return self.make_outputs(sequencing_group, data=output, jobs=job)


@stage.stage(required_stages=[GenotypeWithGatk])
class RunGvcfQc(stage.SequencingGroupStage):
    """
    Run GVCF QC tools on a GVCF file.
    """

    def expected_outputs(self, sequencing_group: targets.SequencingGroup) -> dict[str, Path]:
        qc_prefix = sequencing_group.dataset.prefix() / 'qc'
        return {
            'summary': qc_prefix / f'{sequencing_group.id}.variant_calling_summary_metrics',
            'detail': qc_prefix / f'{sequencing_group.id}.variant_calling_detail_metrics',
        }

    def queue_jobs(self, sequencing_group: targets.SequencingGroup, inputs: stage.StageInput) -> stage.StageOutput:
        """
        Use function from the jobs module
        """

        gvcf_path = inputs.as_str(sequencing_group, GenotypeWithGatk)

        outputs = self.expected_outputs(sequencing_group)

        job = vcf_qc(
            gvcf=gvcf_path,
            output_prefix=str(outputs['summary']).removesuffix('.variant_calling_summary_metrics'),
            job_attrs=self.get_job_attrs(sequencing_group),
        )
        return self.make_outputs(sequencing_group, data=outputs, jobs=job)
