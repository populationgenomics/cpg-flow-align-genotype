"""The per-metric value cache - the small artifact everything downstream reads.

MultiQC reports run to hundreds of megabytes each, so `collect` parses every report
exactly once and distils it to this. Threshold tuning then iterates against numbers
already in memory instead of re-parsing gigabytes, which is what makes the
flagrates/mad loop usable interactively.
"""

import json
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import numpy as np

from align_genotype.qc_calibration.spec import CalibrationSpec


class CacheError(RuntimeError):
    """The cache can't answer the question being asked of it."""


@dataclass(frozen=True)
class CohortValues:
    """One cohort's extracted numbers, plus the provenance the survey reported."""

    label: str
    n_samples: int
    multiqc_version: str
    shape: str
    n_dropped: int
    values: dict[str, list[float]]


@dataclass(frozen=True)
class ValueCache:
    seq_type: str
    generated: str
    complete: bool
    metrics: tuple[str, ...]
    cohorts: tuple[CohortValues, ...]

    @property
    def labels(self) -> tuple[str, ...]:
        return tuple(c.label for c in self.cohorts)

    def cohort(self, label: str) -> CohortValues:
        for c in self.cohorts:
            if c.label == label:
                return c
        raise KeyError(f'{label!r} is not in this cache; known cohorts: {list(self.labels)}')

    def series(self, label: str, metric: str) -> np.ndarray:
        """This cohort's finite values for `metric`, as a float array (empty if absent).

        Non-finite values are filtered here as well as during collect: the cache is a
        plain JSON file an operator may hand-edit, and every downstream percentile and
        MAD calculation assumes finite input.
        """
        raw = self.cohort(label).values.get(metric, [])
        if not raw:
            return np.array([], dtype=float)
        values = np.array([v for v in raw if v is not None], dtype=float)
        return values[np.isfinite(values)]


def save(cache: ValueCache, path: str | Path) -> None:
    """Write the cache as JSON. Raises rather than emitting NaN, which isn't valid JSON."""
    payload: dict[str, Any] = {
        'seq_type': cache.seq_type,
        'generated': cache.generated,
        'complete': cache.complete,
        'metrics': list(cache.metrics),
        'cohorts': {
            c.label: {
                'n_samples': c.n_samples,
                'multiqc_version': c.multiqc_version,
                'shape': c.shape,
                'n_dropped': c.n_dropped,
                'values': {metric: [float(v) for v in values] for metric, values in c.values.items()},
            }
            for c in cache.cohorts
        },
    }
    # Lazy: cpg_utils pulls in cloudpathlib and the S3 client stack, ~1.7s of import
    # time. Keeping it in here lets the rest of qc_calibration - and its tests - run
    # without paying that, and without coupling to cpg-utils config state.
    from cpg_utils import to_path  # noqa: PLC0415

    with to_path(path).open('w') as f:
        json.dump(payload, f, indent=2, allow_nan=False)


def load(path: str | Path) -> ValueCache:
    """Read a cache written by `save`."""
    from cpg_utils import to_path  # noqa: PLC0415

    with to_path(path).open() as f:
        raw = json.load(f)
    return ValueCache(
        seq_type=raw['seq_type'],
        generated=raw['generated'],
        complete=bool(raw['complete']),
        metrics=tuple(raw['metrics']),
        cohorts=tuple(
            CohortValues(
                label=label,
                n_samples=body['n_samples'],
                multiqc_version=body['multiqc_version'],
                shape=body['shape'],
                n_dropped=body['n_dropped'],
                values=body['values'],
            )
            for label, body in raw['cohorts'].items()
        ),
    )


def require_usable(cache: ValueCache, spec: CalibrationSpec) -> None:
    """Raise unless `cache` can answer every question `spec` will ask of it.

    Narrowing the spec's metric list is free (a subset is fine); adding a metric means
    the reports have to be parsed again. Incompleteness is checked first, since it is
    the more fundamental problem - the cache was built against a metric MultiQC never
    surfaced, so its numbers can't be trusted regardless of the metric list.
    """
    if not cache.complete:
        raise CacheError(
            'Cache is incomplete: a gated metric was missing from at least one cohort when it was collected. '
            'Fix the spec or the cohort set, then re-run `qc_calibrate collect`.',
        )
    if missing := sorted(set(spec.metric_keys) - set(cache.metrics)):
        raise CacheError(
            f'Cache is missing {missing} - it was collected for a different metric list. '
            f'Re-run `qc_calibrate collect` to pick the new metric(s) up.',
        )
