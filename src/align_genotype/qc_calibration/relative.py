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
the database, so a relative tier that churns is worse than no relative tier at all. It
is simulated conservatively: two cohorts merging, and one cohort growing under both
plausible readings of the cache's value order, with the worst result deciding.

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
import logging
import os
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

# The "before" slice for homogeneous growth: 60% of a cohort, scored against the
# threshold the full cohort produces.
_GROWTH_FRACTION = 0.6

# Which 60%, though, changes the answer. `cache.series` inherits its order from MultiQC's
# JSON key order, and the same 50 values differing only in order have been measured at
# 16.7% churn (leading slice) versus 0.0% (shuffled slice) - a REJECT and a RECOMMEND for
# one metric. The leading slice models real batch growth *if* report order tracks
# sequencing batches, which is plausible for sequentially-assigned CPG IDs but is
# nowhere guaranteed; the shuffled slice models an arbitrary smaller cohort. Neither
# interpretation is safe to assume, so both are simulated and the verdict uses the worse.
# The shuffle is seeded so an evaluation re-run is reproducible.
_SHUFFLE_SEED = 0

# `relative_flags` stamps every flag with a date. A warn *rate* doesn't depend on it, so
# it's fixed rather than `now()` - an evaluation re-run must be reproducible.
_FIXED_DATE = datetime(2000, 1, 1, tzinfo=timezone.utc)

# Section name for the synthetic MultiQC-shaped cohort handed to `relative_flags`.
_SECTION = 'calibration'


@dataclass(frozen=True)
class CohortMad:
    """One cohort's median, MAD and the warn count the production path produced.

    `threshold` is None exactly when `skipped` is set - a cohort below `min_cohort`, one
    with a degenerate (zero) MAD, or one the cache holds no values for, has no relative
    line and therefore no warn count.

    `n_values` and `n_samples` are both reported because they can differ. `n_values` is
    what the threshold is actually computed over, matching production; `n_samples` is
    what the cohort survey counted. A metric appearing in two MultiQC sections (v1.33 can
    carry both `picard_1` and `picard_4` in the Picard namespace) contributes one value
    per section per sample, so `n_values > n_samples`. Surfacing both makes that
    duplication visible instead of silent - see `warn_rate`.
    """

    label: str
    n_values: int
    n_samples: int
    median: float
    mad_raw: float
    threshold: float | None
    n_warn: int
    skipped: str | None

    @property
    def warn_rate(self) -> float:
        """Warned *values* over total values - deliberately not a per-sample rate.

        The threshold is derived per value, exactly as production derives it, so this is
        the rate consistent with the threshold shown next to it. Under *uniform*
        cross-section duplication it also equals the per-sample rate (both numerator and
        denominator scale together). Under *partial* duplication it overstates it: 26
        samples with 2 outliers present in both sections report 4/28 = 14.3% where the
        true per-sample rate is 2/26 = 7.7% - enough to cross `MAX_WARN_RATE` and flip
        the verdict. `n_values != n_samples` is the signal that this is in play; a true
        per-sample rate is not computable here, since the cache stores values without
        sample identity.
        """
        return self.n_warn / self.n_values if self.n_values else 0.0

    @property
    def duplicated(self) -> bool:
        """Whether this cohort carries more values than samples (see `warn_rate`)."""
        return self.n_values > self.n_samples


@dataclass(frozen=True)
class HomogeneousChurn:
    """One cohort's growth simulation under both before-slice orderings.

    Both slices are the same size and are scored against the same full-cohort threshold;
    only which samples they contain differs. Either can be None when that slice has a
    degenerate MAD. See `_SHUFFLE_SEED` for why both are run.
    """

    label: str
    ordered: stats.ChurnResult | None
    shuffled: stats.ChurnResult | None

    @property
    def flip_rate(self) -> float:
        """The worse of the two orderings - the number the verdict uses."""
        rates = [r.flip_rate for r in (self.ordered, self.shuffled) if r is not None]
        return max(rates) if rates else 0.0

    @property
    def ordering_sensitive(self) -> bool:
        """Whether the two orderings disagree on whether this cohort clears `MAX_CHURN`.

        Diagnostic only. When true, the churn number depends on an assumption about
        MultiQC's key order that the operator should be told about rather than have
        silently resolved for them.
        """
        rates = [r.flip_rate for r in (self.ordered, self.shuffled) if r is not None]
        return any(r > MAX_CHURN for r in rates) and any(r <= MAX_CHURN for r in rates)


@dataclass(frozen=True)
class MadEvaluation:
    """A candidate relative tier's per-cohort numbers, churn simulations and verdict."""

    metric: str
    direction: str
    cohorts: tuple[CohortMad, ...]
    homogeneous: tuple[HomogeneousChurn, ...]
    heterogeneous: tuple[tuple[str, str, stats.ChurnResult], ...]

    @property
    def max_warn_rate(self) -> float:
        """Peak warn rate across cohorts, ignoring skipped ones.

        The skip filter is defensive rather than load-bearing: a skipped cohort has
        `n_warn == 0` and therefore `warn_rate == 0.0`, and adding zeros to `max()`
        cannot change the result. It is here so that relaxing that invariant later -
        recording a would-be count for a cohort production skips, say - can't silently
        start feeding a rate production would never produce into the adoption bar.
        """
        rates = [c.warn_rate for c in self.cohorts if c.skipped is None]
        return max(rates) if rates else 0.0

    @property
    def max_churn(self) -> float:
        """Peak flip rate across every growth simulation; 0.0 when there is no data."""
        rates = [h.flip_rate for h in self.homogeneous]
        rates += [r.flip_rate for _, _, r in self.heterogeneous]
        return max(rates) if rates else 0.0

    @property
    def ordering_sensitive(self) -> tuple[str, ...]:
        """Cohorts whose churn verdict depends on the before-slice ordering."""
        return tuple(h.label for h in self.homogeneous if h.ordering_sensitive)

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
        # Copied because `get_config_paths` hands back the module-level list itself.
        # `set_config_paths` rebinds rather than mutating, so nothing currently writes
        # through the alias - the copy just keeps that a local detail.
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
            # Always attempted: leaking a deleted temp path into global state would break
            # every later `config_retrieve` in the process, including an operator's
            # subsequent commands in the same CLI invocation.
            _restore_config_paths(previous)


def _restore_config_paths(previous: list[str]) -> None:
    """Put the previous config paths back, without letting that failure win.

    `set_config_paths` re-validates: it opens and parses every path. Under
    analysis-runner those are `gs://` paths, so a failure here is not hypothetical - a
    transient GCS or credential error on the restore would replace whatever the caller
    was actually doing with a `ValueError` naming the wrong problem, and would do so
    *after* the evaluation had already succeeded. A failed restore is logged and the
    current paths left in place instead; the caller's own exception, or result, survives.
    """
    try:
        config.set_config_paths(previous)
    except ValueError as exc:
        logging.warning(f'Could not restore previous config paths {previous}, leaving as-is - {exc}')
        return
    if not previous:
        # `set_config_paths([])` writes an empty CPG_CONFIG_PATH, turning "absent" into
        # "present but empty". Identical to cpg-utils, but visible to any subprocess that
        # tests for the key, so restore the absence exactly.
        os.environ.pop('CPG_CONFIG_PATH', None)


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
    n_samples: int,
    metric: MetricSpec,
    relative: RelativeSpec,
    seq_type: str,
) -> CohortMad:
    """One cohort's relative numbers, mirroring the skips `relative_flags` makes.

    The reported `threshold` is `robust_threshold`'s value rounded to 4 dp - the same
    rounding, on the same number, that `relative_flags` applies before comparing. The
    displayed threshold therefore always explains the displayed count.
    """
    n_values = int(values.size)
    if n_values == 0:
        # Distinguished from "too small": the cache holds nothing for this metric here,
        # which is a collection problem, not a cohort-size one. Reporting it as
        # `cohort 0 < min_cohort 50` would send an operator to the wrong place.
        return CohortMad(
            label,
            0,
            n_samples,
            float('nan'),
            float('nan'),
            None,
            0,
            f'metric {metric.key!r} has no values in this cohort',
        )
    median = float(np.median(values))
    mad_raw = float(np.median(np.abs(values - median)))
    if n_values < relative.min_cohort:
        return CohortMad(
            label,
            n_values,
            n_samples,
            median,
            mad_raw,
            None,
            0,
            f'cohort {n_values} < min_cohort {relative.min_cohort}',
        )
    threshold = check_multiqc.robust_threshold(list(values), metric.direction, relative.k)
    if threshold is None:
        return CohortMad(
            label,
            n_values,
            n_samples,
            median,
            mad_raw,
            None,
            0,
            'zero MAD (degenerate cohort); MAD gives no usable line here, use the absolute gate',
        )
    return CohortMad(
        label,
        n_values,
        n_samples,
        median,
        mad_raw,
        round(threshold, 4),
        _warn_count(values, metric, seq_type),
        None,
    )


def _homogeneous_churn(
    usable: dict[str, np.ndarray],
    direction: str,
    k: float,
    min_cohort: int,
) -> tuple[HomogeneousChurn, ...]:
    """Each cohort's 60% before-slice re-scored against the whole cohort's threshold.

    Run twice per cohort - leading slice and seeded-shuffled slice - because the answer
    depends on which 60% is chosen (see `_SHUFFLE_SEED`).

    A before-slice below `min_cohort` is not simulated at all. Production would have
    skipped a cohort that size outright and emitted no flags, so deriving a threshold
    from it and counting flips against it measures churn against a state that cannot
    occur. With the shipped `min_cohort = 50` this excludes every cohort of 50-83
    samples, whose 60% slice lands at 30-49.
    """
    results = []
    for label, values in usable.items():
        size = int(values.size * _GROWTH_FRACTION)
        if size < min_cohort:
            continue
        shuffled = np.random.default_rng(_SHUFFLE_SEED).permutation(values)
        results.append(
            HomogeneousChurn(
                label=label,
                ordered=stats.churn(values[:size], values, direction, k),
                shuffled=stats.churn(shuffled[:size], values, direction, k),
            ),
        )
    # An entry where both orderings were degenerate carries no measurement, so it is
    # dropped exactly as a single degenerate result was before.
    return tuple(h for h in results if h.ordered is not None or h.shuffled is not None)


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
    n_samples = {c.label: c.n_samples for c in cache.cohorts}
    with _production_config(seq_type, metric, relative):
        cohorts = tuple(
            _evaluate_cohort(label, values, n_samples[label], metric, relative, seq_type)
            for label, values in series.items()
        )
    # Cohorts production would skip outright can't churn, so they're excluded rather
    # than contributing a misleading zero flip rate.
    usable = {label: values for label, values in series.items() if values.size >= relative.min_cohort}
    return MadEvaluation(
        metric=metric.key,
        direction=metric.direction,
        cohorts=cohorts,
        homogeneous=_homogeneous_churn(usable, metric.direction, relative.k, relative.min_cohort),
        heterogeneous=_heterogeneous_churn(usable, metric.direction, relative.k),
    )
