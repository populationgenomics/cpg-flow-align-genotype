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


# Direction -> (comparison sign written into the flag, predicate for "breaches this threshold").
DIRECTIONS: dict[str, tuple[str, Any]] = {
    'under': ('<', lambda val, thresh: val < thresh),
    'over': ('>', lambda val, thresh: val > thresh),
}
# Severity tiers, worst first. `fail` is evaluated before `warn`, so a value that
# breaches both is recorded once, as a fail.
SEVERITIES: tuple[str, ...] = ('fail', 'warn')


def load_thresholds(seq_type: str) -> dict[str, dict[str, dict[str, float]]]:
    """Read the nested qc_thresholds config into {direction: {metric: {severity: threshold}}}.

    Config layout is ``[qc_thresholds.<seq_type>.<severity>.<direction>]`` (see
    config_template.toml). A metric may define only some tiers.
    """
    thresholds: dict[str, dict[str, dict[str, float]]] = {direction: {} for direction in DIRECTIONS}
    for severity in SEVERITIES:
        for direction in DIRECTIONS:
            configured = config.config_retrieve(['qc_thresholds', seq_type, severity, direction], {})
            for metric, threshold in configured.items():
                thresholds[direction].setdefault(metric, {})[severity] = threshold
    return thresholds


def worst_breach(val: float, tiers: dict[str, float], direction: str) -> tuple[str, float] | None:
    """Return (severity, threshold) of the most severe tier `val` breaches, else None.

    `tiers` maps severity -> threshold for one metric/direction. `fail` is checked
    before `warn` so a value breaching both is reported once as a fail.
    """
    _, breaches = DIRECTIONS[direction]
    for severity in SEVERITIES:
        if severity in tiers and breaches(val, tiers[severity]):
            return severity, tiers[severity]
    return None


def warn_unmatched_metrics(sections: dict[str, Any], seq_type: str) -> None:
    """Log a warning for any configured threshold metric MultiQC never surfaced.

    Without this, a typo'd or unsurfaced metric key silently checks nothing -
    which is exactly how the old ``PCT_PF_READS_ALIGNED`` genome gate (present only
    in ``report_saved_raw_data``, not ``report_general_stats_data``) went unnoticed.
    Scans every severity tier and direction.
    """
    present_metrics = {
        metric for section in sections.values() for val_by_metric in section.values() for metric in val_by_metric
    }
    thresholds = load_thresholds(seq_type)
    configured_metrics = {metric for by_metric in thresholds.values() for metric in by_metric}

    if not configured_metrics:
        logging.warning(f'No qc_thresholds configured for sequencing_type={seq_type!r}; nothing will be checked.')
    for metric in sorted(configured_metrics - present_metrics):
        logging.warning(
            f'Configured threshold metric {metric!r} not found in any MultiQC section for '
            f'sequencing_type={seq_type!r}; this threshold will not be checked.',
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

    # Log a compact structural summary rather than pprint-ing the whole blob: on a
    # large cohort (e.g. 647 WES samples) the full dump is a multi-MB string built
    # on every run. The per-sample detail is still logged as checks are evaluated.
    sections_summary = ', '.join(f'{name}={len(section)} samples' for name, section in sections.items())
    logging.info(f'report_general_stats_data: {sections_summary}')

    warn_unmatched_metrics(sections, seq_type)

    thresholds = load_thresholds(seq_type)
    bad_lines_by_sample: dict[str, list[str]] = defaultdict(list)
    qc_flags_by_sample: dict[str, list[QcFlag]] = defaultdict(list)
    for direction, metric_tiers in thresholds.items():
        sign = DIRECTIONS[direction][0]
        for section_name, section in sections.items():
            for sample, val_by_metric in section.items():
                for metric, tiers in metric_tiers.items():
                    if metric not in val_by_metric:
                        continue
                    # MultiQC/Picard can emit non-numeric placeholders (e.g. Picard
                    # writes '?' for FOLD_80_BASE_PENALTY when coverage is ~0). Coerce
                    # to float and skip anything we can't compare, rather than crashing
                    # the whole check.
                    try:
                        val = float(val_by_metric[metric])
                    except (TypeError, ValueError):
                        logging.warning(
                            f'{sample}: metric {metric!r} has non-numeric value '
                            f'{val_by_metric[metric]!r}; skipping threshold check.',
                        )
                        continue
                    verdict = worst_breach(val, tiers, direction)
                    if verdict is None:
                        logging.info(f'✅ {sample}: {metric}={val:0.2f} within thresholds')
                        continue
                    severity, threshold = verdict
                    icon = '❗' if severity == 'fail' else '⚠️'
                    line = f'{metric}={val:0.2f}{sign}{threshold:0.2f} [{severity}]'
                    bad_lines_by_sample[sample].append(f'{icon} {line}')
                    sg_id = sample.split('|', 1)[0]
                    qc_flags_by_sample[sg_id].append(
                        QcFlag(
                            flag=metric,
                            value=val,
                            comparison=sign,
                            threshold=threshold,
                            section=section_name,
                            date=today.isoformat(timespec='seconds'),
                            ar_guid=config.try_get_ar_guid(),
                            severity=severity,
                        ),
                    )
                    logging.info(f'{icon} {sample}: {line}')
    logging.info('')

    # Constructing Slack message
    report_title = title or 'MultiQC report'
    title = f'*[{dataset}]* <{html_url}|{report_title}>' if dataset and html_url else report_title
    messages = []
    if bad_lines_by_sample:
        n_fail = sum(1 for flags in qc_flags_by_sample.values() for f in flags if f.severity == 'fail')
        n_warn = sum(1 for flags in qc_flags_by_sample.values() for f in flags if f.severity == 'warn')
        messages.append(
            f'{title}. {len(bad_lines_by_sample)} samples flagged ({n_fail} failing, {n_warn} warnings):',
        )
        for sample, bad_lines in bad_lines_by_sample.items():
            messages.append(f'{sample}: ' + ', '.join(bad_lines))
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
