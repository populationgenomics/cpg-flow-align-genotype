"""The per-dataset values file - the small artifact the report stage reads.

One MultiQC report runs to hundreds of megabytes; this is what it distils to, and it is
the only thing the cross-dataset analysis ever loads. Values carry the sequencing group
they came from, because MultiQC 1.33 can split one tool across general-stats sections
(`picard_1` and `picard_4` both sit in the Picard namespace) and a sequencing group then
contributes one value per section. Threshold derivation keeps that duplication because
production does - `check_multiqc._relative_flags_for_metric` feeds the same doubled list
to `robust_threshold` - but reporting a per-sequencing-group rate needs the identity, and
storing bare floats made that uncomputable.
"""

import json
import math
from dataclasses import dataclass, field
from typing import Any

import numpy as np

from cpg_utils import Path, to_path

# (section, sequencing group, value) - the shape `check_multiqc.gather_metric_values`
# returns, kept verbatim so nothing has to be reshaped on the way in.
Entry = tuple[str, str, float]


class ValuesError(RuntimeError):
    """A values file could not be read, or holds something that must never be written."""


def _require_int(value: Any, what: str) -> int:
    """Reject a bool masquerading as a count.

    `bool` subclasses `int`, so plain `int(value)` accepts `True`/`False` as `1`/`0`
    with no error - a hand-edit or an upstream bug that put a boolean where a count
    belongs would be silently accepted rather than raising.
    """
    if isinstance(value, bool) or not isinstance(value, int):
        raise TypeError(f'{what} must be an integer, got {value!r}')
    return value


def _require_number(value: Any, what: str) -> float:
    """Reject a bool masquerading as a metric value; see `_require_int`."""
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        raise TypeError(f'{what} must be a number, got {value!r}')
    return float(value)


@dataclass(frozen=True)
class MetricValues:
    """One metric's usable values in one dataset, plus how many were unusable."""

    entries: tuple[Entry, ...]
    n_dropped: int

    @property
    def array(self) -> np.ndarray:
        """The values as a float array, in entry order.

        Non-finite values are filtered here as well as in `extract`: this is a plain JSON
        file, and every percentile and MAD downstream assumes finite input.
        """
        if not self.entries:
            return np.array([], dtype=float)
        values = np.array([value for _, _, value in self.entries], dtype=float)
        return values[np.isfinite(values)]

    @property
    def n_values(self) -> int:
        return len(self.entries)

    @property
    def n_sequencing_groups(self) -> int:
        return len({sg for _, sg, _ in self.entries})

    @property
    def sections(self) -> tuple[str, ...]:
        return tuple(sorted({section for section, _, _ in self.entries}))

    @property
    def duplicated(self) -> bool:
        """Whether this metric carries more values than sequencing groups."""
        return self.n_values > self.n_sequencing_groups


@dataclass(frozen=True)
class DatasetValues:
    """One dataset's extracted values, with the provenance of the report they came from."""

    dataset: str
    seq_type: str
    analysis_id: int
    timestamp: str
    uri: str
    multiqc_version: str
    generated: str
    n_sequencing_groups: int
    section_sizes: dict[str, int] = field(default_factory=dict)
    metrics: dict[str, MetricValues] = field(default_factory=dict)

    def metric(self, key: str) -> MetricValues:
        """This dataset's values for `key`, empty rather than missing if absent."""
        return self.metrics.get(key, MetricValues(entries=(), n_dropped=0))


def save(values: DatasetValues, path: str | Path) -> None:
    """Write a values file as JSON.

    Every value is checked finite before anything is opened, so a rejected save leaves no
    partial file behind. A non-finite value reaching this file would poison every
    percentile and MAD in the report stage.
    """
    for key, metric in values.metrics.items():
        for _, sg, value in metric.entries:
            if not math.isfinite(value):
                raise ValuesError(f'{values.dataset}: metric {key!r} sequencing group {sg!r} has a non-finite value')

    payload: dict[str, Any] = {
        'dataset': values.dataset,
        'seq_type': values.seq_type,
        'analysis_id': values.analysis_id,
        'timestamp': values.timestamp,
        'uri': values.uri,
        'multiqc_version': values.multiqc_version,
        'generated': values.generated,
        'n_sequencing_groups': values.n_sequencing_groups,
        'section_sizes': dict(values.section_sizes),
        'metrics': {
            key: {
                'n_dropped': metric.n_dropped,
                'entries': [[section, sg, value] for section, sg, value in metric.entries],
            }
            for key, metric in values.metrics.items()
        },
    }
    # allow_nan=False is belt-and-braces behind the finite check above.
    text = json.dumps(payload, indent=2, allow_nan=False)
    with to_path(path).open('w') as f:
        f.write(text)


def load(path: str | Path) -> DatasetValues:
    """Read a values file written by `save`.

    Every parse failure surfaces as one `ValuesError` naming the path: a missing key, a
    wrong-shaped value and malformed JSON otherwise raise three unrelated exception types
    with no indication of which of N files was at fault.
    """
    try:
        with to_path(path).open() as f:
            raw = json.load(f)
        return DatasetValues(
            dataset=raw['dataset'],
            seq_type=raw['seq_type'],
            analysis_id=_require_int(raw['analysis_id'], 'analysis_id'),
            timestamp=raw['timestamp'],
            uri=raw['uri'],
            multiqc_version=raw['multiqc_version'],
            generated=raw['generated'],
            n_sequencing_groups=_require_int(raw['n_sequencing_groups'], 'n_sequencing_groups'),
            section_sizes=dict(raw['section_sizes']),
            metrics={
                key: MetricValues(
                    entries=tuple(
                        (section, sg, _require_number(value, f'metric {key!r} entry value'))
                        for section, sg, value in body['entries']
                    ),
                    n_dropped=_require_int(body['n_dropped'], f'metric {key!r} n_dropped'),
                )
                for key, body in raw['metrics'].items()
            },
        )
    except (KeyError, TypeError, ValueError, AttributeError) as exc:
        raise ValuesError(f'{path}: not a usable values file - {exc}') from exc
