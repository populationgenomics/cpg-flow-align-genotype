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
import statistics
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
    'min': ('<', lambda val, thresh: val < thresh),
    'max': ('>', lambda val, thresh: val > thresh),
}
# Severity tiers, worst first. `fail` is evaluated before `warn`, so a value that
# breaches both is recorded once, as a fail.
SEVERITIES: tuple[str, ...] = ('fail', 'warn')

# 0.75 quantile of the standard normal; the scaling constant in the Iglewicz-Hoaglin
# modified z-score used for cohort-relative (MAD) flagging.
MODIFIED_Z_CONST = 0.6745


def robust_threshold(values: list[float], direction: str, k: float) -> float | None:
    """Cohort-relative outlier threshold from the modified z-score.

    Returns ``median +/- k*MAD/0.6745`` on the bad side for ``direction``
    ('max' = high is bad, 'min' = low is bad). Returns None when MAD == 0
    (degenerate cohort, e.g. all-identical values) so callers can skip.
    """
    if not values:
        return None
    med = statistics.median(values)
    mad = statistics.median([abs(v - med) for v in values])
    if mad == 0:
        return None
    delta = k * mad / MODIFIED_Z_CONST
    return med + delta if direction == 'max' else med - delta


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


def _gather_metric_values(sections: dict[str, Any], metric: str) -> list[tuple[str, str, float]]:
    """(section, sample, float value) for every sample carrying `metric`; non-numeric skipped."""
    entries: list[tuple[str, str, float]] = []
    for section_name, section in sections.items():
        for sample, val_by_metric in section.items():
            if metric not in val_by_metric:
                continue
            try:
                entries.append((section_name, sample, float(val_by_metric[metric])))
            except (TypeError, ValueError):
                logging.warning(
                    f'{sample}: relative metric {metric!r} non-numeric {val_by_metric[metric]!r}; skipping.',
                )
    return entries


def _relative_flags_for_metric(
    metric: str,
    cfg: dict,
    entries: list[tuple[str, str, float]],
    today: datetime,
    already_flagged: dict[str, set[tuple[str, str]]],
) -> list[tuple[str, str, QcFlag]]:
    """Warn flags for one metric's cohort; empty if MAD is degenerate."""
    direction = cfg['direction']
    sign, breaches = DIRECTIONS[direction]
    threshold = robust_threshold([v for _, _, v in entries], direction, cfg['k'])
    if threshold is None:
        logging.warning(f'Relative flagging skipped for {metric!r}: zero MAD (degenerate cohort).')
        return []
    threshold = round(threshold, 4)  # damp sub-0.0001 jitter -> fewer spurious "updated" churns
    logging.info(f'Relative threshold for {metric!r} ({direction}, k={cfg["k"]}): {threshold}')

    results: list[tuple[str, str, QcFlag]] = []
    for section_name, sample, val in entries:
        if not breaches(val, threshold):
            continue
        sg_id = sample.split('|', 1)[0]
        if (section_name, metric) in already_flagged.get(sg_id, set()):
            continue  # already flagged absolutely (e.g. fail) - don't double up
        results.append((
            sample,
            sg_id,
            QcFlag(
                flag=metric,
                value=val,
                comparison=sign,
                threshold=threshold,
                section=section_name,
                date=today.isoformat(timespec='seconds'),
                ar_guid=config.try_get_ar_guid(),
                severity='warn',
                method='relative',
            ),
        ))
    return results


def relative_flags(
    sections: dict[str, Any],
    seq_type: str,
    today: datetime,
    already_flagged: dict[str, set[tuple[str, str]]],
) -> list[tuple[str, str, QcFlag]]:
    """Cohort-relative (MAD) warn-only flags for metrics in the `relative` config.

    For each configured metric, gathers every sample's value across the current run
    (the "cohort"), derives a robust median/MAD outlier threshold, and warns samples
    beyond it. Relative flags are always ``severity='warn'`` / ``method='relative'``,
    and are skipped when the cohort is smaller than ``min_cohort`` or MAD is zero.
    ``already_flagged`` maps sg_id -> {(section, metric)} flagged by the absolute pass;
    those are not double-flagged, so the absolute fail gate takes precedence.

    Returns (sample, sg_id, QcFlag) tuples.
    """
    spec = config.config_retrieve(['qc_thresholds', seq_type, 'relative'], {})
    results: list[tuple[str, str, QcFlag]] = []
    for metric, cfg in spec.items():
        entries = _gather_metric_values(sections, metric)
        min_cohort = cfg.get('min_cohort', 0)
        if len(entries) < min_cohort:
            logging.info(f'Relative flagging skipped for {metric!r}: cohort {len(entries)} < min_cohort {min_cohort}.')
            continue
        results.extend(_relative_flags_for_metric(metric, cfg, entries, today, already_flagged))
    return results


def apply_relative_flags(
    sections: dict[str, Any],
    seq_type: str,
    today: datetime,
    qc_flags_by_sample: dict[str, list[QcFlag]],
    bad_lines_by_sample: dict[str, list[str]],
) -> None:
    """Run the relative pass and merge its warn flags into the accumulators in place.

    Absolute flags are computed first; `already_flagged` lets the relative pass defer
    to them so a sample isn't flagged twice for the same metric.
    """
    already_flagged = {sg_id: {(f.section, f.flag) for f in flags} for sg_id, flags in qc_flags_by_sample.items()}
    for sample, sg_id, flag in relative_flags(sections, seq_type, today, already_flagged):
        qc_flags_by_sample[sg_id].append(flag)
        line = f'{flag.flag}={flag.value:0.2f}{flag.comparison}{flag.threshold:0.2f} [warn·relative]'
        bad_lines_by_sample[sample].append(f'⚠️ {line}')
        logging.info(f'⚠️ {sample}: {line}')


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

    # Cohort-relative (MAD) warn-only pass. Runs after the absolute pass so it can
    # skip metrics already flagged (absolute fail takes precedence over relative warn).
    apply_relative_flags(sections, seq_type, today, qc_flags_by_sample, bad_lines_by_sample)
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
