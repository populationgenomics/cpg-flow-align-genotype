"""Distil one dataset's MultiQC report into a values file.

Run as a Hail Batch job by `QcCalibrationDatasetMetrics`, with the report already
localised by `batch.read_input`. Also runnable off a checkout against a downloaded
report, which is the quickest way to check what a new MultiQC version surfaces.
"""

import json
import logging
from datetime import datetime, timezone

import click

from cpg_utils import to_path

from align_genotype.qc_calibration import extract, settings, values

logging.basicConfig(level=logging.INFO, format='%(levelname)s: %(message)s')


@click.command()
@click.option('--dataset', required=True, help='Dataset name, used to label the values file.')
@click.option('--multiqc-json', 'multiqc_json_path', required=True, help='Path to the MultiQC JSON.')
@click.option('--analysis-id', type=int, required=True, help='Metamist analysis ID the report came from.')
@click.option('--timestamp', required=True, help='timestampCompleted of that analysis.')
@click.option('--uri', required=True, help='The report URI, recorded for provenance.')
@click.option('--output', 'output_path', required=True, help='Where to write the values file.')
def main(
    dataset: str,
    multiqc_json_path: str,
    analysis_id: int,
    timestamp: str,
    uri: str,
    output_path: str,
) -> None:
    """Extract every configured calibration metric from one MultiQC report."""
    loaded = settings.load()
    logging.info(f'{dataset}: extracting {len(loaded.metrics)} {loaded.seq_type} metrics from {multiqc_json_path}')

    with to_path(multiqc_json_path).open() as f:
        document = json.load(f)

    result = extract.extract(
        document,
        loaded,
        dataset=dataset,
        analysis_id=analysis_id,
        timestamp=timestamp,
        uri=uri,
        generated=datetime.now(tz=timezone.utc).isoformat(timespec='seconds'),
    )
    # `document` is the multi-hundred-MB parsed report; `result` is the small distillate
    # `extract` built from it. Drop the reference before the write below so the parsed
    # document's memory isn't held alive across it - peak memory in this job is the
    # `json.load` above, and there is no reason for the write to extend that peak.
    del document

    values.save(result, output_path)

    for key, metric_values in result.metrics.items():
        if not metric_values.entries:
            logging.warning(f'{dataset}: {key} is absent from this report')
        elif metric_values.n_dropped:
            logging.warning(f'{dataset}: {key} lost {metric_values.n_dropped} unusable value(s)')
    logging.info(f'{dataset}: wrote {output_path} ({result.n_sequencing_groups} sequencing groups)')


if __name__ == '__main__':
    main()
