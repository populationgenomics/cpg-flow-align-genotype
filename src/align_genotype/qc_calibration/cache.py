"""The per-metric value cache - the small artifact everything downstream reads.

MultiQC reports run to hundreds of megabytes each, so `collect` parses every report
exactly once and distils it to this. Threshold tuning then iterates against numbers
already in memory instead of re-parsing gigabytes, which is what makes the
flagrates/mad loop usable interactively.
"""

from __future__ import annotations

import json
import os
import tempfile
from dataclasses import dataclass
from pathlib import Path
from typing import TYPE_CHECKING, Any

import numpy as np

if TYPE_CHECKING:
    # `from __future__ import annotations` above means annotations are never evaluated
    # at runtime, so both of these - like the lazy `to_path` imports below - cost
    # nothing outside of type checking. `cpg_utils.Path` is the CloudPath | pathlib.Path
    # union `to_path` can hand back.
    from cpg_utils import Path as CpgPath

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
        MAD calculation assumes finite input. A value that isn't numeric at all (a
        typo'd string, a stray table) is the same hand-edit risk, so it raises
        `CacheError` naming the cohort and metric rather than a bare numpy exception.
        """
        raw = self.cohort(label).values.get(metric, [])
        if not raw:
            return np.array([], dtype=float)
        try:
            values = np.array([v for v in raw if v is not None], dtype=float)
        except (TypeError, ValueError) as exc:
            raise CacheError(f'cache: cohort {label!r} metric {metric!r} has a non-numeric value - {exc}') from exc
        return values[np.isfinite(values)]


def _coerce_metric_values(label: str, metric: str, values: list) -> list[float]:
    """Cast every value to float for writing, or fail loudly naming cohort and metric.

    `series` tolerates a hand-edited `None` or junk value on read (see its docstring)
    because a hand-edited cache is expected to need defending against. Writing one back
    out is a different situation: `collect` never emits a non-numeric value, so if one
    reaches `save` it's either a bug or a hand-edit that already broke the contract this
    function exists to keep. Dropping it silently would hide that; a plain `TypeError`
    from `float()` would report it with no cohort or metric context.
    """
    coerced = []
    for v in values:
        try:
            coerced.append(float(v))
        except (TypeError, ValueError) as exc:
            raise CacheError(
                f'cache: cohort {label!r} metric {metric!r} has a non-numeric value {v!r} - {exc}',
            ) from exc
    return coerced


def save(cache: ValueCache, path: str | Path | CpgPath) -> None:
    """Write the cache as JSON, all-or-nothing.

    This is the one artifact in the workflow that costs ten minutes of report-parsing
    to regenerate, so a failed or interrupted `save` must never destroy a previously
    good cache at the same path. Guards, in order:

    - Every value is cast with `float()` first (`_coerce_metric_values`); a `None` or
      otherwise non-numeric value raises `CacheError` naming the cohort and metric.
    - The payload is then serialised with `json.dumps` before anything touches disk, so
      a stray NaN (`allow_nan=False`) raises with nothing written, rather than partway
      through streaming into an already-open target file.
    - On local disk, the write lands in a sibling temp file - via `tempfile.mkstemp` in
      the same directory, so two concurrent saves can't collide on one temp name and
      the final rename stays same-filesystem - and is only moved onto `path` via
      `Path.replace`, an atomic rename. The temp file is removed on any failure so a
      crash doesn't leave a full-size unexplained file next to the cache.

    `CloudPath.replace` is deliberately not used as the non-local equivalent of that
    temp-then-replace dance: cloudpathlib's implementation unlinks the destination
    *before* copying the source (`cloudpathlib/cloudpath.py`, `Path.replace`: `if
    target.exists(): target.unlink()`), so a tmp-then-replace write on a CloudPath
    would delete a good cache and then perform a non-transactional copy - strictly
    worse than writing directly, whether or not that copy happens to be server-side
    (e.g. GCS's `copy_blob`). A direct `open('w')` write is already a single atomic
    per-object upload on the object stores we use, so that's what the non-local branch
    does. Caches are written locally in practice, so there's nothing to gain chasing
    this further.
    """
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
                'values': {
                    metric: _coerce_metric_values(c.label, metric, values) for metric, values in c.values.items()
                },
            }
            for c in cache.cohorts
        },
    }
    text = json.dumps(payload, indent=2, allow_nan=False)

    # Lazy: cpg_utils pulls in cloudpathlib and the S3 client stack, ~1.7s of import
    # time. Keeping it in here lets the rest of qc_calibration - and its tests - run
    # without paying that, and without coupling to cpg-utils config state.
    from cpg_utils import to_path  # noqa: PLC0415

    target = to_path(path)
    if isinstance(target, Path):
        fd, tmp_name = tempfile.mkstemp(dir=target.parent, prefix=f'{target.name}.', suffix='.tmp')
        tmp = Path(tmp_name)
        try:
            with os.fdopen(fd, 'w') as f:
                f.write(text)
            tmp.replace(target)
        finally:
            # A no-op on the success path (`replace` already consumed `tmp`); on
            # failure this is what stops a crash from leaving a full-size unexplained
            # file next to the cache.
            tmp.unlink(missing_ok=True)
    else:
        with target.open('w') as f:
            f.write(text)


def load(path: str | Path | CpgPath) -> ValueCache:
    """Read a cache written by `save`.

    Wraps the parse in a single `try` so a missing key, a wrong-shaped value, or
    malformed JSON all surface as one `CacheError` naming the path, instead of a bare
    `KeyError`/`TypeError`/`json.JSONDecodeError` with no indication of which file or
    that it's the cache format that's the problem.
    """
    from cpg_utils import to_path  # noqa: PLC0415

    try:
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
    except (KeyError, TypeError, AttributeError, json.JSONDecodeError) as exc:
        raise CacheError(f'{path}: not a usable value cache - {exc}') from exc


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
