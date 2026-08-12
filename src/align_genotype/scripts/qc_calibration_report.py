"""Turn every dataset's values file into the calibration report.

Run as a Hail Batch job by `QcCalibrationReport`, with the values files already localised
by `batch.read_input`. Also runnable off a checkout against downloaded values files,
which is the quickest way to try a different `k` or metric list without re-parsing any
MultiQC report:

    python -m align_genotype.scripts.qc_calibration_report \\
        --values a.json --values b.json --output-json out.json --output-html out.html

A metric absent from every dataset, or a shipped `qc_thresholds` warn tier that is inert,
is logged at ERROR and banners the report, but neither fails the job: Hail only copies
`write_output` targets on success, so failing here would destroy the page that explains
the problem.
"""

import json
import logging
from datetime import datetime, timezone

import click

from cpg_utils import config, to_path

from align_genotype.qc_calibration import render, settings, summary, values

logging.basicConfig(level=logging.INFO, format='%(levelname)s: %(message)s')

SKIP_REASON = 'no completed CramMultiQC qc analysis for {seq_type}'


@click.command()
@click.option('--values', 'value_paths', multiple=True, required=True, help='A dataset values file. Repeatable.')
@click.option('--skipped-dataset', 'skipped', multiple=True, help='A dataset with no report. Repeatable.')
@click.option('--output-json', 'output_json', required=True, help='Where to write calibration.json.')
@click.option('--output-html', 'output_html', required=True, help='Where to write calibration.html.')
def main(value_paths: tuple[str, ...], skipped: tuple[str, ...], output_json: str, output_html: str) -> None:
    """Assemble the cross-dataset calibration report."""
    loaded = settings.load()
    datasets = [values.load(path) for path in value_paths]
    logging.info(f'Read {len(datasets)} dataset values file(s) for {loaded.seq_type}')

    built = summary.build(
        datasets,
        loaded,
        current=settings.current_thresholds(loaded.seq_type),
        skipped_datasets=[
            {'dataset': name, 'reason': SKIP_REASON.format(seq_type=loaded.seq_type)} for name in skipped
        ],
        generated=datetime.now(tz=timezone.utc).isoformat(timespec='seconds'),
        ar_guid=config.try_get_ar_guid() or 'unknown',
    )

    with to_path(output_json).open('w') as f:
        json.dump(built, f, indent=2, allow_nan=False)
    with to_path(output_html).open('w') as f:
        f.write(render.render(built))

    for warning in built['warnings']:
        if warning['severity'] == 'error':
            logging.error(warning['message'])
        else:
            logging.warning(warning['message'])
    for key, evaluation in built['relative'].items():
        logging.info(f'{key}: {evaluation["verdict"]} {evaluation["reason"]}')
    logging.info(f'Wrote {output_json} and {output_html}')


if __name__ == '__main__':
    main()
