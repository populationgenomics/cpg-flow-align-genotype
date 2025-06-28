"""
Alignment and genotyping... oof. Where to start.
"""

from cpg_flow import stage, targets
from cpg_utils import Path

from align_genotype.jobs.align import align
from align_genotype.jobs.genotype import genotype


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
        return sequencing_group.gvcf or sequencing_group.make_gvcf_path()

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
