"""Dataset-relative (MAD) tier evaluation.

Some metrics have no defensible fixed warn line because their normal level shifts by
dataset or protocol - duplication rate on whole genomes is the canonical case. A fixed
warn line either floods the high-duplication datasets or never fires on the low ones, so
the warn tier is derived per run from that run's own median and MAD (an Iglewicz-Hoaglin
modified z-score line) with an absolute `fail` gate behind it.

Nothing here re-implements the modified z-score. Thresholds come from
`check_multiqc.robust_threshold` - production's own function - and a warn is counted by
breaching the same 4-dp-rounded threshold `_relative_flags_for_metric` compares against.

This module answers two questions about a candidate tier: how many values it would warn
on per dataset, and whether the flag set stays stable as the dataset changes. The second
decides adoption. Every flip is a value that stops or starts being flagged purely because
the dataset's median moved, which is a spurious "updated" flag in the database, so a tier
that churns is worse than no tier at all.
"""

import itertools
from dataclasses import dataclass

import numpy as np

from align_genotype.qc_calibration import stats
from align_genotype.qc_calibration.settings import Bars, CalibrationSettings, MetricSpec
from align_genotype.qc_calibration.values import MetricValues
from align_genotype.scripts import check_multiqc

# The "before" slice for growth: 60% of a dataset, scored against the threshold the whole
# dataset produces.
_GROWTH_FRACTION = 0.6

# Which 60% changes the answer. Entry order is inherited from MultiQC's JSON key order,
# and the same 50 values differing only in order have measured 16.7% churn on the leading
# slice against 0.0% on a shuffled one - a REJECT and a RECOMMEND for one metric. The
# leading slice models real batch growth *if* report order tracks sequencing batches,
# which is plausible for sequentially-assigned CPG IDs but nowhere guaranteed; the
# shuffled slice models an arbitrary smaller dataset. Neither is safe to assume, so both
# are simulated and the verdict takes the worse. Seeded so a re-run is reproducible.
_SHUFFLE_SEED = 0


@dataclass(frozen=True)
class DatasetMad:
    """One dataset's median, MAD, derived threshold and warn count.

    `threshold` is None exactly when `skipped` is set - a dataset below `min_samples`, one
    with a degenerate (zero) MAD, or one with no values for this metric has no relative
    line and therefore no warn count.
    """

    dataset: str
    n_values: int
    n_groups_with_values: int
    median: float
    mad_raw: float
    threshold: float | None
    n_warn: int
    skipped: str | None

    @property
    def warn_rate(self) -> float:
        """Warned *values* over total values - deliberately not a per-group rate.

        The threshold is derived per value, exactly as production derives it, so this is
        the rate consistent with the threshold shown beside it. Where a metric appears in
        two MultiQC sections it overstates the per-sequencing-group rate; `duplicated` is
        the signal that this is in play, and the values file carries the identity needed
        to say by how much.
        """
        return self.n_warn / self.n_values if self.n_values else 0.0

    @property
    def duplicated(self) -> bool:
        return self.n_values > self.n_groups_with_values


@dataclass(frozen=True)
class GrowthChurn:
    """One dataset's growth simulation under both before-slice orderings.

    Both slices are the same size and are scored against the same whole-dataset
    threshold; only which values they contain differs.
    """

    dataset: str
    ordered: stats.ChurnResult | None
    shuffled: stats.ChurnResult | None

    @property
    def flip_rate(self) -> float:
        """The worse of the two orderings - the number the verdict uses."""
        rates = [r.flip_rate for r in (self.ordered, self.shuffled) if r is not None]
        return max(rates) if rates else 0.0

    def ordering_sensitive(self, bar: float) -> bool:
        """Whether the two orderings disagree about clearing `bar`.

        Diagnostic only. When true, the churn number depends on an assumption about
        MultiQC key order that a reader should be told about rather than have silently
        resolved for them.
        """
        rates = [r.flip_rate for r in (self.ordered, self.shuffled) if r is not None]
        return any(r > bar for r in rates) and any(r <= bar for r in rates)


@dataclass(frozen=True)
class MadEvaluation:
    """A candidate relative tier's per-dataset numbers, churn simulations and verdict."""

    metric: str
    direction: str
    bars: Bars
    datasets: tuple[DatasetMad, ...]
    growth: tuple[GrowthChurn, ...]
    merge: tuple[tuple[str, str, stats.ChurnResult], ...]

    @property
    def max_warn_rate(self) -> float:
        """Peak warn rate across datasets, ignoring skipped ones."""
        rates = [d.warn_rate for d in self.datasets if d.skipped is None]
        return max(rates) if rates else 0.0

    @property
    def max_growth_churn(self) -> float:
        rates = [g.flip_rate for g in self.growth]
        return max(rates) if rates else 0.0

    @property
    def max_merge_churn(self) -> float:
        rates = [r.flip_rate for _, _, r in self.merge]
        return max(rates) if rates else 0.0

    @property
    def ordering_sensitive(self) -> tuple[str, ...]:
        return tuple(g.dataset for g in self.growth if g.ordering_sensitive(self.bars.max_growth_churn))

    @property
    def verdict_reason(self) -> str:
        """Why this tier misses the bar, or '' when it clears every one.

        Each simulation is named against its own bar. A reader given only a conflated
        "peak churn" cannot tell whether to re-scope the dataset set or reconsider the
        metric, and those have different answers.
        """
        reasons = []
        if self.max_warn_rate > self.bars.max_warn_rate:
            reasons.append(f'peak warn rate {self.max_warn_rate:.1%} exceeds {self.bars.max_warn_rate:.0%}')
        # Churn to 2 dp, warn rate to 1: churn values sit close to their bars, and
        # "2.0% exceeds 2%" reads as a contradiction where "2.01% exceeds 2%" does not.
        if self.max_growth_churn > self.bars.max_growth_churn:
            reasons.append(
                f'peak dataset-growth churn {self.max_growth_churn:.2%} exceeds {self.bars.max_growth_churn:.0%}',
            )
        if self.max_merge_churn > self.bars.max_merge_churn:
            reasons.append(
                f'peak cross-dataset merge churn {self.max_merge_churn:.2%} exceeds {self.bars.max_merge_churn:.0%}',
            )
        return '; '.join(reasons)

    @property
    def verdict(self) -> str:
        return 'REJECT' if self.verdict_reason else 'RECOMMEND'


def _evaluate_dataset(
    dataset: str,
    metric_values: MetricValues,
    metric: MetricSpec,
    min_samples: int,
    k: float,
) -> DatasetMad:
    """One dataset's relative numbers, mirroring the skips production makes."""
    values = metric_values.array
    n_values = int(values.size)
    counts = {'n_values': n_values, 'n_groups_with_values': metric_values.n_groups_with_values}
    if n_values == 0:
        # Distinguished from "too small": nothing was extracted for this metric here,
        # which is a collection problem, not a size one. Reporting it as
        # "0 values < min_samples 50" would send a reader to the wrong place.
        return DatasetMad(
            dataset,
            **counts,
            median=float('nan'),
            mad_raw=float('nan'),
            threshold=None,
            n_warn=0,
            skipped=f'metric {metric.key!r} has no values in this dataset',
        )
    median = float(np.median(values))
    mad_raw = float(np.median(np.abs(values - median)))
    if n_values < min_samples:
        return DatasetMad(
            dataset,
            **counts,
            median=median,
            mad_raw=mad_raw,
            threshold=None,
            n_warn=0,
            skipped=f'{n_values} values < min_samples {min_samples}',
        )
    threshold = check_multiqc.robust_threshold(list(values), metric.direction, k)
    if threshold is None:
        return DatasetMad(
            dataset,
            **counts,
            median=median,
            mad_raw=mad_raw,
            threshold=None,
            n_warn=0,
            skipped='zero MAD (degenerate dataset); use the absolute gate here',
        )
    # Production rounds before comparing, so the displayed threshold always explains the
    # displayed count.
    threshold = round(threshold, 4)
    return DatasetMad(
        dataset,
        **counts,
        median=median,
        mad_raw=mad_raw,
        threshold=threshold,
        n_warn=int(stats.breach(values, threshold, metric.direction).sum()),
        skipped=None,
    )


def _growth_churn(
    usable: dict[str, np.ndarray],
    direction: str,
    k: float,
    min_samples: int,
) -> tuple[GrowthChurn, ...]:
    """Each dataset's 60% before-slice re-scored against the whole dataset's threshold.

    A before-slice below `min_samples` is not simulated: production would have skipped a
    dataset that size outright and emitted no flags, so deriving a threshold from it and
    counting flips measures churn against a state that cannot occur.
    """
    results = []
    for dataset, values in usable.items():
        size = int(values.size * _GROWTH_FRACTION)
        if size < min_samples:
            continue
        shuffled = np.random.default_rng(_SHUFFLE_SEED).permutation(values)
        results.append(
            GrowthChurn(
                dataset=dataset,
                ordered=stats.churn(values[:size], values, direction, k),
                shuffled=stats.churn(shuffled[:size], values, direction, k),
            ),
        )
    # An entry where both orderings were degenerate carries no measurement.
    return tuple(g for g in results if g.ordered is not None or g.shuffled is not None)


def _merge_churn(
    usable: dict[str, np.ndarray],
    direction: str,
    k: float,
) -> tuple[tuple[str, str, stats.ChurnResult], ...]:
    """Each dataset re-scored against the threshold it gets once another joins it.

    An entry `(a, b, result)` is *a's* flag set after b merges in. A merge disturbs both
    projects' flag sets by different amounts, so both directions of every pair run -
    `n*(n-1)` entries, not `n*(n-1)/2`. One direction per pair would leave the headline
    figure depending on insertion order, the same hazard `_SHUFFLE_SEED` guards against.

    The harsher simulation by construction, and there is a tension worth naming: a metric
    earns a relative tier precisely *because* its normal level shifts between datasets,
    and a merge punishes exactly that. Read a high figure as "how much would pooling two
    projects disturb this", not as a defect count.
    """
    results = []
    for (label_a, values_a), (label_b, values_b) in itertools.permutations(usable.items(), 2):
        result = stats.churn(values_a, np.concatenate([values_a, values_b]), direction, k)
        if result is not None:
            results.append((label_a, label_b, result))
    return tuple(results)


def evaluate(
    by_dataset: dict[str, MetricValues],
    metric: MetricSpec,
    settings: CalibrationSettings,
) -> MadEvaluation:
    """Evaluate `metric`'s dataset-relative warn tier across every dataset supplied."""
    if not metric.relative:
        raise ValueError(
            f'metric {metric.key!r} has no configured relative tier to evaluate; '
            f'set relative = true on [qc_calibration.{settings.seq_type}.metrics.{metric.key}]',
        )
    datasets = tuple(
        _evaluate_dataset(name, metric_values, metric, settings.min_samples, settings.k)
        for name, metric_values in by_dataset.items()
    )
    # Datasets production would skip outright cannot churn, so they are excluded rather
    # than contributing a misleading zero flip rate.
    usable = {
        name: metric_values.array
        for name, metric_values in by_dataset.items()
        if metric_values.array.size >= settings.min_samples
    }
    return MadEvaluation(
        metric=metric.key,
        direction=metric.direction,
        bars=settings.bars,
        datasets=datasets,
        growth=_growth_churn(usable, metric.direction, settings.k, settings.min_samples),
        merge=_merge_churn(usable, metric.direction, settings.k),
    )
