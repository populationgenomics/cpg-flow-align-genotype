"""
Alignment and genotyping... oof. Where to start.
"""

from cpg_flow import stage, targets
from cpg_utils import Path

from align_genotype.jobs.align import align
from align_genotype.jobs.genotype import genotype
from align_genotype.jobs.CramQcSomalier import extract_somalier
from align_genotype.jobs.CramQcVerifyBamId import verifybamid
from align_genotype.jobs.CramQcSamtoolsStats import samtools_stats


@stage.stage(
    analysis_type='cram',
    analysis_keys=['cram'],
)
class AlignWithDragmap(stage.SequencingGroupStage):
    """
    This is a generic stage that runs a bash command.
    """

    def expected_outputs(self, sequencing_group: targets.SequencingGroup) -> Path:
        cram_path = sequencing_group.cram if sequencing_group.cram else sequencing_group.make_cram_path()
        return {
            'cram': cram_path.path,
            'sorted_bam': str(
                sequencing_group.dataset.prefix(category='tmp') / 'align' / f'{sequencing_group.id}.sorted.bam'
            ),
            'markdup': (
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
    required_stages=AlignWithDragmap,
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

        jobs = genotype(
            output_path=output,
            sequencing_group_name=sequencing_group.id,
            cram_path=cram,
            tmp_prefix=self.tmp_prefix / sequencing_group.id,
            job_attrs=self.get_job_attrs(sequencing_group),
        )
        return self.make_outputs(sequencing_group, data=output, jobs=jobs)


@stage.stage(required_stages=AlignWithDragmap)
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


@stage.stage(required_stages=AlignWithDragmap)
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


@stage.stage(required_stages=AlignWithDragmap)
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
