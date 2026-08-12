"""QC threshold calibration stages.

An isolated branch of the DAG: neither stage declares a `required_stages` dependency on
a production stage, so adding them to `run_workflow` never pulls alignment or genotyping
work into a calibration run, and never adds calibration work to a production one.

Both stages return no jobs unless `qc_calibration.enabled` is set. That flag ships false,
because a dependency-free stage whose outputs do not exist would otherwise queue a job
per dataset and register Metamist analyses on every production invocation. A calibration
run sets it true and passes
`only_stages = ['QcCalibrationDatasetMetrics', 'QcCalibrationReport']` - both names, since
`only_stages` marks every unlisted stage skipped, so naming just the report would make
the extract stage check for outputs rather than produce them.

Every decision lives in the free functions below rather than in the stage methods:
`Stage.__init__` calls `get_workflow()`, so a stage cannot be constructed outside a real
run and its methods cannot be unit tested.
"""

from cpg_flow import stage, targets
from cpg_utils import Path, config

from align_genotype.jobs import qc_calibration
from align_genotype.qc_calibration import discovery, settings
from align_genotype.qc_calibration.discovery import MultiqcReport


def sequencing_type() -> str:
    return config.config_retrieve(['workflow', 'sequencing_type'])


def dataset_report(dataset_name: str) -> MultiqcReport | None:
    """That dataset's latest CramMultiQC report, memoised across DAG assembly.

    The Metamist project name carries a `-test` suffix at test access level, and this
    name goes into a GraphQL query, so it must be resolved rather than used raw.
    """
    return discovery.cached_latest_cram_multiqc(
        config.dataset_for_access_level(dataset_name),
        sequencing_type(),
    )


def dataset_values_path(dataset: targets.Dataset, seq_type: str, report: MultiqcReport) -> Path:
    """Where one dataset's values file lives.

    Keyed on the analysis ID, so a new MultiQC report produces a new path and extraction
    re-runs exactly when the underlying data changed - never on a stale file.

    The literal sequencing type, not `targets.sequencing_subdir()`, which returns '' for
    genome and would leave the two runs' outputs asymmetric.
    """
    return dataset.prefix() / 'qc_calibration' / seq_type / f'values.{report.analysis_id}.json'


def dataset_outputs(dataset: targets.Dataset) -> dict[str, Path]:
    """The extract stage's expected outputs, empty when there is nothing to do.

    cpg-flow treats a falsy expected output as reusable, so the action becomes REUSE,
    `queue_jobs` is never called, and no Metamist analysis is created.
    """
    if not settings.enabled():
        return {}
    report = dataset_report(dataset.name)
    if report is None:
        return {}
    return {'values': dataset_values_path(dataset, sequencing_type(), report)}


def report_outputs(analysis_dataset: targets.Dataset) -> dict[str, Path]:
    """The report stage's expected outputs: JSON to main, HTML to web."""
    if not settings.enabled():
        return {}
    seq_type = sequencing_type()
    return {
        'json': analysis_dataset.prefix() / 'qc_calibration' / seq_type / 'calibration.json',
        'html': analysis_dataset.web_prefix() / 'qc_calibration' / seq_type / 'calibration.html',
    }


def collect_values_paths(
    datasets: list[targets.Dataset],
    seq_type: str,
) -> tuple[dict[str, Path], list[str]]:
    """Split the multicohort's datasets into those with a values file and those without.

    Paths are recomputed from the same memoised discovery the extract stage used rather
    than read out of `StageInput`: `StageInput._each` raises when *no* dataset produced an
    output, which is a legitimate state here. `required_stages` still supplies the job
    ordering, which is the part `inputs` is actually needed for.
    """
    found: dict[str, Path] = {}
    skipped: list[str] = []
    for dataset in datasets:
        report = dataset_report(dataset.name)
        if report is None:
            skipped.append(dataset.name)
        else:
            found[dataset.name] = dataset_values_path(dataset, seq_type, report)
    return found, skipped


@stage.stage
class QcCalibrationDatasetMetrics(stage.DatasetStage):
    """Distil one dataset's latest CramMultiQC report to a small values file."""

    def expected_outputs(self, dataset: targets.Dataset) -> dict[str, Path]:
        return dataset_outputs(dataset)

    def queue_jobs(self, dataset: targets.Dataset, inputs: stage.StageInput) -> stage.StageOutput:  # noqa: ARG002
        outputs = self.expected_outputs(dataset)
        if not outputs:
            return self.make_outputs(dataset, data=None, jobs=None)

        job = qc_calibration.extract_dataset_metrics(
            report=dataset_report(dataset.name),
            output=outputs['values'],
            job_attrs=self.get_job_attrs(dataset),
        )
        return self.make_outputs(dataset, data=outputs, jobs=job)


@stage.stage(
    required_stages=[QcCalibrationDatasetMetrics],
    analysis_type='web',
    analysis_keys=['html'],
    forced=True,
)
class QcCalibrationReport(stage.MultiCohortStage):
    """Cross-dataset analysis and the HTML dashboard the teams review."""

    def expected_outputs(self, multicohort: targets.MultiCohort) -> dict[str, Path]:
        return report_outputs(multicohort.analysis_dataset)

    def queue_jobs(self, multicohort: targets.MultiCohort, inputs: stage.StageInput) -> stage.StageOutput:  # noqa: ARG002
        outputs = self.expected_outputs(multicohort)
        if not outputs:
            return self.make_outputs(multicohort, data=None, jobs=None)

        values_paths, skipped = collect_values_paths(multicohort.get_datasets(), sequencing_type())
        if not values_paths:
            raise ValueError(
                f'No dataset in this multicohort has a completed CramMultiQC {sequencing_type()} analysis, '
                f'so there is nothing to calibrate from. Check that CramMultiQC has run for these '
                f'datasets at this sequencing type.',
            )

        job = qc_calibration.calibration_report(
            values_paths=values_paths,
            skipped_datasets=skipped,
            outputs=outputs,
            job_attrs=self.get_job_attrs(multicohort),
        )
        return self.make_outputs(multicohort, data=outputs, jobs=job)
