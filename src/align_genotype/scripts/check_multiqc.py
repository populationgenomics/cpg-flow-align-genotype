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

import click
from cpg_utils import config, to_path
from cpg_utils.slack import send_message

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
def main(
    multiqc_json_path: str,
    html_url: str | None = None,
    dataset: str | None = None,
    title: str | None = None,
    send_to_slack: bool = True,
):
    """
    Check metrics in MultiQC json and send info about failed samples
    as a Slack message.
    """
    run(
        multiqc_json_path=multiqc_json_path,
        html_url=html_url,
        dataset=dataset,
        title=title,
        send_to_slack=send_to_slack,
    )


def run(  # noqa: C901, PLR0912
    multiqc_json_path: str,
    html_url: str | None = None,
    dataset: str | None = None,
    title: str | None = None,
    send_to_slack: bool = True,
):
    seq_type = config.config_retrieve(['workflow', 'sequencing_type'])

    with to_path(multiqc_json_path).open() as f:
        d = json.load(f)
        sections = d['report_general_stats_data']
        logging.info(f'report_general_stats_data: {pprint.pformat(sections)}')

    bad_lines_by_sample = defaultdict(list)
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
        for section in sections:
            for sample, val_by_metric in section.items():
                for metric, threshold in threshold_d.items():
                    if metric in val_by_metric:
                        val = val_by_metric[metric]
                        if is_fail(val, threshold):
                            line = f'{metric}={val:0.2f}{fail_sign}{threshold:0.2f}'
                            bad_lines_by_sample[sample].append(line)
                            logging.info(f'❗ {sample}: {line}')
                        else:
                            line = f'{metric}={val:0.2f}{good_sign}{threshold:0.2f}'
                            logging.info(f'✅ {sample}: {line}')
    logging.info('')

    # Constructing Slack message
    if dataset and html_url:
        title = f'*[{dataset}]* <{html_url}|{title or "MultiQC report"}>'
    elif not title:
        title = 'MultiQC report'
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


if __name__ == '__main__':
    main()  # pylint: disable=E1120
