"""Cohort-relative (MAD) evaluation, driven through the production code path.

Some metrics have no defensible fixed warn line because their normal level shifts by
cohort or protocol - duplication rate on whole genomes is the canonical case, where
cohort medians span roughly 7% to 18% depending on library prep. A fixed warn line
either floods the high-duplication cohorts or never fires on the low-duplication ones,
so the warn tier is derived per-run from the cohort's own median and MAD (an
Iglewicz-Hoaglin modified z-score line) with an absolute `fail` gate behind it.

This module answers two questions about a candidate relative tier:

1. How many samples would it warn on, per cohort?
2. Does the flag set stay stable as the cohort grows?

The second question is the one that decides adoption. Every flip - a sample that stops
being flagged purely because the cohort's median moved - is a spurious "updated" flag in
the database, so a relative tier that churns is worse than no relative tier at all.

There is deliberately no modified z-score implementation here. Warn counts come from
calling ``check_multiqc.relative_flags`` for real, against a throwaway config file, and
churn comes from ``stats.churn``, which calls ``check_multiqc.robust_threshold``
directly. The manual workflow this replaces had a MAD prototype plus a second script
that verified the prototype against production; keeping both invites exactly the
divergence the verification existed to catch. Routing through the genuine config ->
``load_thresholds`` -> ``relative_flags`` path also proves, incidentally, that the config
shape being emitted is loadable.
"""

import itertools
import tempfile
from collections.abc import Iterator
from contextlib import contextmanager
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path

import numpy as np

from cpg_utils import config

from align_genotype.qc_calibration import stats, tomlio
from align_genotype.qc_calibration.cache import ValueCache
from align_genotype.qc_calibration.spec import MetricSpec, RelativeSpec
from align_genotype.scripts import check_multiqc

# The adoption bar. These are the values that historically admitted exome
# ZERO_CVG_TARGETS_PCT and genome reads_duplicated_percent, and rejected
# PCT_SELECTED_BASES and PCT_OFF_BAIT - both of which churned up to 24.5% of their flag
# set on simulated cohort growth. Both are advisory: `verdict` is a recommendation for
# an operator to sign off, not an automatic gate.
MAX_WARN_RATE = 0.10
MAX_CHURN = 0.02

# The "before" slice for homogeneous growth: score the first 60% of a cohort against the
# threshold the full cohort produces, i.e. what an operator would have seen if this
# cohort had been sequenced in two batches.
_GROWTH_FRACTION = 0.6

# `relative_flags` stamps every flag with a date. A warn *rate* doesn't depend on it, so
# it's fixed rather than `now()` - an evaluation re-run must be reproducible.
_FIXED_DATE = datetime(2000, 1, 1, tzinfo=timezone.utc)

# Section name for the synthetic MultiQC-shaped cohort handed to `relative_flags`.
_SECTION = 'calibration'


@dataclass(frozen=True)
class CohortMad:
    """One cohort's median, MAD and the warn count the production path produced.

    `threshold` is None exactly when `skipped` is set - a cohort below `min_cohort`, or
    one with a degenerate (zero) MAD, has no relative line and therefore no warn count.
    """

    label: str
    n: int
    median: float
    mad_raw: float
    threshold: float | None
    n_warn: int
    skipped: str | None

    @property
    def warn_rate(self) -> float:
        return self.n_warn / self.n if self.n else 0.0


@dataclass(frozen=True)
class MadEvaluation:
    """A candidate relative tier's per-cohort numbers, churn simulations and verdict."""

    metric: str
    direction: str
    cohorts: tuple[CohortMad, ...]
    homogeneous: tuple[tuple[str, stats.ChurnResult], ...]
    heterogeneous: tuple[tuple[str, str, stats.ChurnResult], ...]

    @property
    def max_warn_rate(self) -> float:
        """Peak warn rate across cohorts, ignoring skipped ones.

        A skipped cohort's `n_warn` is 0 by construction, so including it would mask a
        real peak elsewhere rather than reporting "no data" for that cohort.
        """
        rates = [c.warn_rate for c in self.cohorts if c.skipped is None]
        return max(rates) if rates else 0.0

    @property
    def max_churn(self) -> float:
        """Peak flip rate across both growth simulations; 0.0 when there is no data."""
        rates = [r.flip_rate for _, r in self.homogeneous]
        rates += [r.flip_rate for _, _, r in self.heterogeneous]
        return max(rates) if rates else 0.0

    @property
    def verdict_reason(self) -> str:
        """Why this tier misses the adoption bar, or '' when it clears it."""
        reasons = []
        if self.max_warn_rate > MAX_WARN_RATE:
            reasons.append(f'peak warn rate {self.max_warn_rate:.1%} exceeds {MAX_WARN_RATE:.0%}')
        if self.max_churn > MAX_CHURN:
            reasons.append(f'peak cohort-growth churn {self.max_churn:.1%} exceeds {MAX_CHURN:.0%}')
        return '; '.join(reasons)

    @property
    def verdict(self) -> str:
        return 'REJECT' if self.verdict_reason else 'RECOMMEND'


@contextmanager
def _production_config(seq_type: str, metric: MetricSpec, relative: RelativeSpec) -> Iterator[None]:
    """Point cpg-utils at a throwaway config carrying just this metric's relative spec.

    Written as a real file and installed with `set_config_paths` rather than
    monkeypatching `config_retrieve`, so the numbers come from the genuine config ->
    `load_thresholds` -> `relative_flags` path. `set_config_paths` only validates that
    the files exist, end in `.toml` and parse, so no other keys are needed.

    `get_config_paths` *raises* when nothing has ever been set, so a cold start is
    treated as "no previous paths" and restored to that.
    """
    tomlio.require_bare_key(seq_type, 'sequencing type')
    tomlio.require_bare_key(metric.key, 'metric key')
    try:
        previous = list(config.get_config_paths())
    except config.ConfigError:
        previous = []
    text = '\n'.join(
        [
            '[workflow]',
            tomlio.fmt_kv('sequencing_type', seq_type),
            '',
            f'[qc_thresholds.{seq_type}.relative.{metric.key}]',
            tomlio.fmt_kv('direction', metric.direction),
            tomlio.fmt_kv('k', relative.k),
            tomlio.fmt_kv('min_cohort', relative.min_cohort),
            '',
        ],
    )
    with tempfile.TemporaryDirectory(prefix='qc-calibration-') as tmpdir:
        path = Path(tmpdir) / 'relative.toml'
        path.write_text(text)
        config.set_config_paths([str(path)])
        try:
            yield
        finally:
            # Unconditional: leaking a deleted temp path into global state would break
            # every later `config_retrieve` in the process, including an operator's
            # subsequent commands in the same CLI invocation.
            config.set_config_paths(previous)


def _warn_count(values: np.ndarray, metric: MetricSpec, seq_type: str) -> int:
    """Warn flags production would raise on `values`, via `check_multiqc.relative_flags`.

    The cohort is reshaped into the ``{section: {sample: {metric: value}}}`` form
    `relative_flags` consumes. Sample names are positional (`S0`, `S1`, ...): the count
    is all that's wanted, and a synthetic name can't be mistaken for a real sample ID.
    `already_flagged` is empty because an absolute fail gate suppressing a relative warn
    would understate the tier's warn rate, which is the number being calibrated.
    """
    sections = {_SECTION: {f'S{i}': {metric.key: float(v)} for i, v in enumerate(values)}}
    return len(check_multiqc.relative_flags(sections, seq_type, _FIXED_DATE, already_flagged={}))


def _evaluate_cohort(
    label: str,
    values: np.ndarray,
    metric: MetricSpec,
    relative: RelativeSpec,
    seq_type: str,
) -> CohortMad:
    """One cohort's relative numbers, mirroring the two skips `relative_flags` makes.

    The reported `threshold` is `robust_threshold`'s value rounded to 4 dp - the same
    rounding, on the same number, that `relative_flags` applies before comparing. The
    displayed threshold therefore always explains the displayed count.
    """
    n = int(values.size)
    median = float(np.median(values)) if n else float('nan')
    mad_raw = float(np.median(np.abs(values - median))) if n else float('nan')
    if n < relative.min_cohort:
        return CohortMad(label, n, median, mad_raw, None, 0, f'cohort {n} < min_cohort {relative.min_cohort}')
    threshold = check_multiqc.robust_threshold(list(values), metric.direction, relative.k)
    if threshold is None:
        return CohortMad(label, n, median, mad_raw, None, 0, 'zero MAD (degenerate cohort)')
    return CohortMad(label, n, median, mad_raw, round(threshold, 4), _warn_count(values, metric, seq_type), None)


def _homogeneous_churn(
    usable: dict[str, np.ndarray],
    direction: str,
    k: float,
) -> tuple[tuple[str, stats.ChurnResult], ...]:
    """Each cohort's first 60% re-scored against the whole cohort's threshold."""
    results = []
    for label, values in usable.items():
        result = stats.churn(values[: int(values.size * _GROWTH_FRACTION)], values, direction, k)
        if result is not None:
            results.append((label, result))
    return tuple(results)


def _heterogeneous_churn(
    usable: dict[str, np.ndarray],
    direction: str,
    k: float,
) -> tuple[tuple[str, str, stats.ChurnResult], ...]:
    """Each cohort re-scored against the threshold it gets once a second one joins it.

    The harsher of the two simulations, and the realistic one: cohorts differ by
    protocol, so a merged run shifts the median further than growth within one cohort.
    Pairs are unordered and never self-paired.
    """
    results = []
    for (label_a, values_a), (label_b, values_b) in itertools.combinations(usable.items(), 2):
        result = stats.churn(values_a, np.concatenate([values_a, values_b]), direction, k)
        if result is not None:
            results.append((label_a, label_b, result))
    return tuple(results)


def evaluate(cache: ValueCache, metric: MetricSpec, seq_type: str) -> MadEvaluation:
    """Evaluate `metric`'s cohort-relative warn tier across every cohort in `cache`.

    Per-cohort warn counts are computed inside a single throwaway-config block - one
    temp config for the whole metric, since the spec written into it doesn't vary by
    cohort. Churn is computed outside it: `stats.churn` calls `robust_threshold`
    directly and needs no config at all.
    """
    relative = metric.relative
    if relative is None:
        raise ValueError(
            f'metric {metric.key!r} has no configured relative tier to evaluate; '
            f'add a [metrics.{metric.key}.relative] block to the calibration spec.',
        )
    series = {label: cache.series(label, metric.key) for label in cache.labels}
    with _production_config(seq_type, metric, relative):
        cohorts = tuple(_evaluate_cohort(label, values, metric, relative, seq_type) for label, values in series.items())
    # Cohorts production would skip outright can't churn, so they're excluded rather
    # than contributing a misleading zero flip rate.
    usable = {label: values for label, values in series.items() if values.size >= relative.min_cohort}
    return MadEvaluation(
        metric=metric.key,
        direction=metric.direction,
        cohorts=cohorts,
        homogeneous=_homogeneous_churn(usable, metric.direction, relative.k),
        heterogeneous=_heterogeneous_churn(usable, metric.direction, relative.k),
    )
