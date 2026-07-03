"""
Checks metrics in MultiQC output, based on thresholds in the qc_thresholds
config section.

Script can send a report to a Slack channel. To enable that, set SLACK_TOKEN
and SLACK_CHANNEL environment variables, and add "Seqr Loader" app into
a channel with:

/invite @Seqr Loader
"""

import json
import logging
import pprint
from collections import defaultdict
from dataclasses import asdict
from datetime import datetime
from typing import Any

import click
from cpg_utils import config, to_path
from cpg_utils.slack import send_message

from align_genotype.utils import QcFlag

logging.basicConfig()
logging.getLogger().setLevel(logging.DEBUG)


@click.command()
@click.option(
    '--multiqc-json',
    'multiqc_json_path',
    required=True,
    help='Path to MultiQC JSON output',
)
@click.option(
    '--html-url',
    'html_url',
    help='MultiQC HTML URL',
)
@click.option('--dataset', 'dataset', help='Dataset name')
@click.option('--title', 'title', help='Report title')
@click.option(
    '--send-to-slack/--no-send-to-slack',
    'send_to_slack',
    help='Send log to Slack message, according to environment variables SLACK_CHANNEL and SLACK_TOKEN',
)
@click.option(
    '--output-json',
    'output_json_path',
    help='Path to write structured QC flags JSON output',
)
def main(
    multiqc_json_path: str,
    html_url: str | None = None,
    dataset: str | None = None,
    title: str | None = None,
    send_to_slack: bool = True,
    output_json_path: str | None = None,
):
    """
    Check metrics in MultiQC json and compare them against thresholds, then send info about failed samples
    as a Slack message and save structured QC flags JSON output to a file if specified.
    """
    run(
        multiqc_json_path=multiqc_json_path,
        html_url=html_url,
        dataset=dataset,
        title=title,
        send_to_slack=send_to_slack,
        output_json_path=output_json_path,
    )


def run(  # noqa: C901
    multiqc_json_path: str,
    html_url: str | None = None,
    dataset: str | None = None,
    title: str | None = None,
    send_to_slack: bool = True,
    output_json_path: str | None = None,
) -> dict[str, Any]:
    seq_type = config.config_retrieve(['workflow', 'sequencing_type'])

    today = datetime.now()  # noqa: DTZ005

    with to_path(multiqc_json_path).open() as f:
        d = json.load(f)
        sections = d['report_general_stats_data']
        logging.info(f'report_general_stats_data: {pprint.pformat(sections)}')

    bad_lines_by_sample: dict[str, list[str]] = defaultdict(list)
    qc_flags_by_sample: dict[str, list[QcFlag]] = defaultdict(list)
    for config_key, fail_sign, good_sign, is_fail in [
        (
            'min',
            '<',
            '≥',
            lambda val_, thresh_: val_ < thresh_,
        ),
        (
            'max',
            '>',
            '≤',
            lambda val_, thresh_: val_ > thresh_,
        ),
    ]:
        threshold_d = config.config_retrieve(['qc_thresholds', seq_type, config_key], {})
        for section_name, section in sections.items():
            for sample, val_by_metric in section.items():
                for metric, threshold in threshold_d.items():
                    if metric in val_by_metric:
                        val = val_by_metric[metric]
                        if is_fail(val, threshold):
                            line = f'{metric}={val:0.2f}{fail_sign}{threshold:0.2f}'
                            bad_lines_by_sample[sample].append(line)
                            sg_id = sample.split('|', 1)[0]
                            qc_flags_by_sample[sg_id].append(
                                QcFlag(
                                    flag=metric,
                                    value=val,
                                    comparison=fail_sign,
                                    threshold=threshold,
                                    section=section_name,
                                    date=today.isoformat(timespec='seconds'),
                                    ar_guid=config.try_get_ar_guid(),
                                ),
                            )
                            logging.info(f'❗ {sample}: {line}')
                        else:
                            line = f'{metric}={val:0.2f}{good_sign}{threshold:0.2f}'
                            logging.info(f'✅ {sample}: {line}')
    logging.info('')

    # Constructing Slack message
    report_title = title or 'MultiQC report'
    title = f'*[{dataset}]* <{html_url}|{report_title}>' if dataset and html_url else report_title
    messages = []
    if bad_lines_by_sample:
        messages.append(f'{title}. {len(bad_lines_by_sample)} samples are flagged:')
        for sample, bad_lines in bad_lines_by_sample.items():
            messages.append(f'❗ {sample}: ' + ', '.join(bad_lines))
    else:
        messages.append(f'✅ {title}')
    text = '\n'.join(messages)
    logging.info(text)

    if send_to_slack:
        send_message(text)

    result: dict[str, Any] = {
        'title': report_title,
        'dataset': dataset,
        'html_url': html_url,
        'sequencing_type': seq_type,
        'n_samples_flagged': len(qc_flags_by_sample),
        'qc_flags': {sample: [asdict(flag) for flag in flags] for sample, flags in qc_flags_by_sample.items()},
    }

    if output_json_path:
        with to_path(output_json_path).open('w') as f:
            json.dump(result, f, indent=2)

    return result


if __name__ == '__main__':
    main()  # pylint: disable=E1120
