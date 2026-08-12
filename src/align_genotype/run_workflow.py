#!/usr/bin/env python3

from argparse import ArgumentParser

from cpg_flow.workflow import run_workflow

from align_genotype.qc_calibration_stages import QcCalibrationDatasetMetrics, QcCalibrationReport
from align_genotype.stages import (
    CramQcPicardCollectMetrics,
    CramQcPicardMultiMetrics,
    CramQcSamtoolsStats,
    CramQcSomalier,
    CramQcVerifyBamId,
    GenotypeWithGatk,
    RunGvcfQc,
    VntyperIndexPage,
)

STAGES = [
    GenotypeWithGatk,
    CramQcPicardMultiMetrics,
    CramQcPicardCollectMetrics,
    CramQcSomalier,
    CramQcSamtoolsStats,
    CramQcVerifyBamId,
    RunGvcfQc,
    VntyperIndexPage,
    # An isolated branch: no dependency on the stages above, and inert unless
    # `qc_calibration.enabled` is set. See qc_calibration/README.md.
    QcCalibrationDatasetMetrics,
    QcCalibrationReport,
]


def cli_main():
    """
    CLI entrypoint - starts up the workflow
    """
    parser = ArgumentParser()
    parser.add_argument('--dry_run', action='store_true', help='Dry run')
    args = parser.parse_args()

    run_workflow(name='align_genotype', stages=STAGES, dry_run=args.dry_run)


if __name__ == '__main__':
    cli_main()
