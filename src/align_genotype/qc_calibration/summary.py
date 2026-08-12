"""Assemble every calibration result into the one structure the outputs are built from.

This dict is both `calibration.json` and the template context. Keeping one
representation means the HTML can never show a number the JSON does not carry.
"""

import math
from typing import Any

import numpy as np

from align_genotype.qc_calibration import relative, snippet, stats, thresholds
from align_genotype.qc_calibration.settings import CalibrationSettings, MetricSpec
from align_genotype.qc_calibration.values import DatasetValues, MetricValues


def _percentiles(values: np.ndarray) -> dict[str, float]:
    measured = stats.percentiles(values)
    return {f'p{pct}': value for pct, value in zip(stats.PERCENTILES, measured, strict=True)}


def _rates(values: np.ndarray, metric: MetricSpec, tiers: dict[str, float | None]) -> dict[str, float | None]:
    """Fail and warn rates under one pair of thresholds, `None` where nothing is scored.

    `stats.flag_rates` returns nan for an empty array; nan is not JSON-serialisable under
    `allow_nan=False`, and "no values" is a different statement from "0%", so it maps to
    None rather than being silently zeroed.
    """
    fail_rate, warn_rate = stats.flag_rates(values, metric.direction, tiers.get('fail'), tiers.get('warn'))
    return {
        'fail': None if math.isnan(fail_rate) or tiers.get('fail') is None else fail_rate,
        'warn': None if math.isnan(warn_rate) or tiers.get('warn') is None else warn_rate,
    }


def _churn_result(result: stats.ChurnResult | None) -> dict[str, Any] | None:
    if result is None:
        return None
    return {
        'threshold_before': result.threshold_before,
        'threshold_after': result.threshold_after,
        'n_initial': result.n_initial,
        'flagged_before': result.flagged_before,
        'flagged_after': result.flagged_after,
        'flips': result.flips,
        'flip_rate': result.flip_rate,
    }


def _relative_block(evaluation: relative.MadEvaluation) -> dict[str, Any]:
    """One relative metric's evaluation, with the denominator behind every peak figure.

    `max_warn_rate`, `max_growth_churn` and `max_merge_churn` are each computed over a
    different subset of datasets - skipped datasets drop out of the warn-rate peak,
    growth additionally drops any dataset whose 60% before-slice would itself fall below
    `min_samples`, and merge needs at least two usable datasets to produce any pairs at
    all. A reader who sees only "peak warn rate 4.2%" has no way to tell whether that
    describes every dataset in the run or a fraction of them, so the counts behind each
    figure are reported alongside it rather than left to a docstring nobody reading the
    HTML will see.
    """
    n_skipped = sum(1 for d in evaluation.datasets if d.skipped is not None)
    return {
        'metric': evaluation.metric,
        'direction': evaluation.direction,
        'verdict': evaluation.verdict,
        'reason': evaluation.verdict_reason,
        'bars': {
            'max_warn_rate': evaluation.bars.max_warn_rate,
            'max_growth_churn': evaluation.bars.max_growth_churn,
            'max_merge_churn': evaluation.bars.max_merge_churn,
        },
        'max_warn_rate': evaluation.max_warn_rate,
        'max_growth_churn': evaluation.max_growth_churn,
        'max_merge_churn': evaluation.max_merge_churn,
        # The denominators behind the three peaks above - see the docstring.
        'n_datasets': len(evaluation.datasets),
        'n_datasets_skipped': n_skipped,
        'n_datasets_evaluated': len(evaluation.datasets) - n_skipped,
        'n_growth_evaluated': len(evaluation.growth),
        'n_merge_datasets_evaluated': len({name for a, b, _ in evaluation.merge for name in (a, b)}),
        'ordering_sensitive': list(evaluation.ordering_sensitive),
        'datasets': [
            {
                'dataset': d.dataset,
                'n_values': d.n_values,
                'n_groups_with_values': d.n_groups_with_values,
                'median': None if math.isnan(d.median) else d.median,
                'mad': None if math.isnan(d.mad_raw) else d.mad_raw,
                'threshold': d.threshold,
                'n_warn': d.n_warn,
                'warn_rate': d.warn_rate,
                'duplicated': d.duplicated,
                'skipped': d.skipped,
            }
            for d in evaluation.datasets
        ],
        'growth': [
            {
                'dataset': g.dataset,
                'ordered': _churn_result(g.ordered),
                'shuffled': _churn_result(g.shuffled),
                'worse': g.flip_rate,
                'ordering_sensitive': g.ordering_sensitive(evaluation.bars.max_growth_churn),
            }
            for g in evaluation.growth
        ],
        # Merge is quadratic in dataset count - 240 ordered pairs at 16 datasets - and
        # only the worst few decide anything. The total is reported so nothing is hidden.
        'merge_worst': [
            {'dataset': a, 'merged_with': b, **_churn_result(r)}
            for a, b, r in sorted(evaluation.merge, key=lambda pair: pair[2].flip_rate, reverse=True)[:MAX_PAIRS]
        ],
        'merge_pairs_simulated': len(evaluation.merge),
    }


# How many merge pairs the report lists, worst first.
MAX_PAIRS = 10


def build(
    values: list[DatasetValues],
    settings: CalibrationSettings,
    *,
    current: dict[str, dict[str, float]],
    skipped_datasets: list[dict[str, str]],
    generated: str,
    ar_guid: str,
) -> dict[str, Any]:
    """Assemble the calibration result from every dataset's extracted values."""
    by_name: dict[str, DatasetValues] = {v.dataset: v for v in values}
    warnings: list[str] = []

    metrics: dict[str, Any] = {}
    for metric in settings.metrics:
        per_dataset: dict[str, MetricValues] = {name: v.metric(metric.key) for name, v in by_name.items()}
        with_data = {name: mv for name, mv in per_dataset.items() if mv.array.size}
        missing_from = sorted(name for name, mv in per_dataset.items() if not mv.array.size)

        candidate = thresholds.candidate(per_dataset, metric) if per_dataset else None
        candidate_tiers: dict[str, float | None] = (
            {'fail': candidate.fail, 'warn': candidate.warn} if candidate else {'fail': None, 'warn': None}
        )
        shipped = current.get(metric.key, {})

        metrics[metric.key] = {
            'direction': metric.direction,
            'unit': metric.unit,
            'relative': metric.relative,
            'n_values': sum(mv.n_values for mv in per_dataset.values()),
            'n_groups_with_values': len({sg for mv in per_dataset.values() for _, sg, _ in mv.entries}),
            'n_datasets': len(with_data),
            'n_dropped': sum(mv.n_dropped for mv in per_dataset.values()),
            'present_in': sorted({section for mv in per_dataset.values() for section in mv.sections}),
            'missing_from': missing_from,
            'duplicated_in': sorted(name for name, mv in per_dataset.items() if mv.duplicated),
            'current': dict(shipped),
            'candidate': (
                {'fail': candidate.fail, 'warn': candidate.warn, 'basis': candidate.basis} if candidate else None
            ),
            'flag_rates': {
                'current': {name: _rates(mv.array, metric, shipped) for name, mv in with_data.items()},
                'candidate': {name: _rates(mv.array, metric, candidate_tiers) for name, mv in with_data.items()},
            },
            'percentiles': {name: _percentiles(mv.array) for name, mv in with_data.items()},
        }

        # A metric absent everywhere is the PCT_PF_READS_ALIGNED class of bug: a key that
        # checks nothing at all. Absent from some datasets is a narrower question.
        if not with_data:
            warnings.append(
                f'{metric.key} was absent from every dataset - the key checks nothing. '
                f'Either it is misspelled for this sequencing type, or MultiQC writes it '
                f'only to report_saved_raw_data, which the production check never reads.',
            )
        elif missing_from:
            warnings.append(f'{metric.key} was absent from {len(missing_from)} dataset(s): {", ".join(missing_from)}')

    relative_blocks: dict[str, Any] = {}
    for metric in settings.relative_metrics:
        per_dataset = {name: v.metric(metric.key) for name, v in by_name.items()}
        if not per_dataset:
            continue
        relative_blocks[metric.key] = _relative_block(relative.evaluate(per_dataset, metric, settings))

    candidates = {
        key: thresholds.Candidate(
            metric=key,
            fail=body['candidate']['fail'],
            warn=body['candidate']['warn'],
            basis=body['candidate']['basis'],
        )
        for key, body in metrics.items()
        if body['candidate'] is not None
    }

    return {
        'sequencing_type': settings.seq_type,
        'generated': generated,
        'ar_guid': ar_guid,
        'settings': {
            'k': settings.k,
            'min_samples': settings.min_samples,
            'max_warn_rate': settings.bars.max_warn_rate,
            'max_growth_churn': settings.bars.max_growth_churn,
            'max_merge_churn': settings.bars.max_merge_churn,
        },
        'datasets': [
            {
                'dataset': v.dataset,
                'analysis_id': v.analysis_id,
                'timestamp': v.timestamp,
                'uri': v.uri,
                'n_sequencing_groups': v.n_sequencing_groups,
                'multiqc_version': v.multiqc_version,
                'section_sizes': v.section_sizes,
            }
            for v in values
        ],
        'skipped_datasets': list(skipped_datasets),
        'metrics': metrics,
        'relative': relative_blocks,
        'warnings': warnings,
        'config_snippet': snippet.render(settings, candidates),
    }
