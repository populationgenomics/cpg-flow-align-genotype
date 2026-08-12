# QC Calibration Workflow Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Convert the `qc_calibrate` operator CLI into two CPG Flow stages that derive candidate QC thresholds from every dataset's latest CramMultiQC report and render an HTML dashboard for the genomics and bioinformatics teams to review.

**Architecture:** A `DatasetStage` fans out one Hail Batch job per dataset, each localising that dataset's `multiqc_data.json` and distilling it to a small values file in the dataset's own bucket. A `MultiCohortStage` then reads all of those, does the cross-dataset analysis (percentile-tail threshold candidates, per-dataset flag rates, MAD tiers, growth and merge churn), and writes `calibration.json` to the analysis dataset's main bucket plus `calibration.html` to its web bucket. Neither stage depends on a production stage, so they form an isolated branch of the DAG.

**Tech Stack:** Python 3.10–3.11, cpg-flow ~1.3, cpg-utils, Hail Batch, numpy, jinja2, click, pytest. Design spec: `docs/superpowers/specs/2026-08-12-qc-calibration-workflow-design.md`.

---

## File Structure

**Deleted** (the CLI surface):

| Path | Why |
|---|---|
| `src/align_genotype/qc_calibration/cli.py` | the stage is the entrypoint |
| `src/align_genotype/qc_calibration/tomlio.py` | nothing reads or writes TOML any more |
| `src/align_genotype/qc_calibration/manifest.py` | no manifest artifact; provenance goes in the report JSON |
| `src/align_genotype/qc_calibration/dryrun.py` | proved the emitted block loads; production proves that every run |
| `src/align_genotype/qc_calibration/spec.py` | replaced by `settings.py` |
| `src/align_genotype/qc_calibration/cache.py` | replaced by `values.py` |
| `src/align_genotype/qc_calibration/collect.py` | replaced by `extract.py` |
| `src/align_genotype/qc_calibration/suggest.py` | replaced by `thresholds.py` |
| `src/align_genotype/qc_calibration/emit.py` | replaced by `snippet.py` |
| `src/align_genotype/qc_calibration/report.py` | replaced by `summary.py` + a Jinja template |
| `test/test_qc_calibration_{tomlio,manifest,cli,dryrun,spec,cache,collect,suggest,emit,report}.py` | with their modules |

**Created:**

| Path | Responsibility |
|---|---|
| `src/align_genotype/qc_calibration/settings.py` | read `[qc_calibration]` config into `MetricSpec` / `CalibrationSettings`; read shipped `qc_thresholds` for comparison |
| `src/align_genotype/qc_calibration/values.py` | the per-dataset values file: dataclasses plus JSON read/write |
| `src/align_genotype/qc_calibration/extract.py` | one MultiQC document to one `DatasetValues` |
| `src/align_genotype/qc_calibration/thresholds.py` | candidate fixed thresholds from percentile tails |
| `src/align_genotype/qc_calibration/snippet.py` | render the `[qc_thresholds.<seq_type>...]` block |
| `src/align_genotype/qc_calibration/summary.py` | assemble every result into the one dict that serves both the JSON output and the template |
| `src/align_genotype/qc_calibration/render.py` | Jinja environment and template render |
| `src/align_genotype/templates/qc_calibration_report.html.jinja` | the dashboard |
| `src/align_genotype/scripts/qc_calibration_extract.py` | per-dataset job entrypoint |
| `src/align_genotype/scripts/qc_calibration_report.py` | report job entrypoint |
| `src/align_genotype/jobs/qc_calibration.py` | the two job builders |
| `src/align_genotype/qc_calibration_stages.py` | the two stages |

**Modified:** `src/align_genotype/scripts/check_multiqc.py` (`min_cohort` → `min_samples`), `src/align_genotype/config_template.toml`, `src/align_genotype/qc_calibration/{stats,relative,discovery}.py`, `src/align_genotype/run_workflow.py`, `pyproject.toml`.

---

## Task 1: Rename `min_cohort` to `min_samples` in shipped config

This touches production behaviour, so it lands as its own commit reviewable in isolation. In CPG Flow terms the value means "the minimum number of sequencing groups in this MultiQC run", and since `CramMultiQC` is a `DatasetStage` that run is one dataset — `cohort` is exactly the word this whole piece of work is removing.

**Files:**
- Modify: `src/align_genotype/scripts/check_multiqc.py:307`
- Modify: `src/align_genotype/config_template.toml:151`, `:212`
- Test: `test/test_check_multiqc.py`

- [ ] **Step 1: Write the failing test**

Add to `test/test_check_multiqc.py`, after `test_relative_skipped_below_min_cohort` (find it with `grep -n min_cohort test/test_check_multiqc.py`):

```python
def test_relative_uses_min_samples_key(patch_config):
    """The config key is `min_samples`; a run of 3 must be skipped by a bar of 4."""
    patch_config(
        'genome',
        {
            'relative': {
                'reads_duplicated_percent': {'direction': 'max', 'k': 3.5, 'min_samples': 4},
            },
        },
    )
    sections = {
        'samtools': {
            f'CPG{i}': {'reads_duplicated_percent': value}
            for i, value in enumerate([10.0, 11.0, 50.0])
        },
    }
    flags = check_multiqc.relative_flags(sections, 'genome', datetime(2026, 1, 1), already_flagged={})
    assert flags == []
```

Note `test_check_multiqc.py` already imports `datetime`; confirm with `grep -n '^from datetime' test/test_check_multiqc.py` and add `from datetime import datetime` if absent.

**Three data points is the minimum that makes this test mean anything.** With exactly two
points the modified z-score is always exactly `±0.6745`, because MAD equals half the
range - so `|mz| > k` is false for every `k` this codebase uses, and no two-point set can
ever be flagged whether or not the size guard runs. A two-point version of this test
passes identically before and after the rename: a vacuous regression guard on the only
production behaviour change in this plan. The same trap applies to any small MAD fixture
you are tempted to build later; use at least three points with a genuine outlier.

- [ ] **Step 2: Run the test to verify it fails**

Run: `uv run pytest test/test_check_multiqc.py::test_relative_uses_min_samples_key -v`

Expected: FAIL. Before the rename `min_samples` is ignored, `cfg.get('min_cohort', 0)`
defaults the bar to `0`, so the three-value run is not skipped, `CPG2` (50.0) is flagged
as an outlier, and `flags` is non-empty. After the rename the bar is 4, `3 < 4`, and the
metric is skipped.

- [ ] **Step 3: Make the rename**

In `src/align_genotype/scripts/check_multiqc.py`, inside `relative_flags`, change:

```python
        min_cohort = cfg.get('min_cohort', 0)
        if len(entries) < min_cohort:
            logging.info(f'Relative flagging skipped for {metric!r}: cohort {len(entries)} < min_cohort {min_cohort}.')
            continue
```

to:

```python
        min_samples = cfg.get('min_samples', 0)
        if len(entries) < min_samples:
            logging.info(
                f'Relative flagging skipped for {metric!r}: '
                f'{len(entries)} values < min_samples {min_samples}.',
            )
            continue
```

Also update the `relative_flags` docstring: replace `` ``min_cohort`` `` with `` ``min_samples`` ``.

- [ ] **Step 4: Update the config template**

In `src/align_genotype/config_template.toml`, under `[qc_thresholds.genome.relative.reads_duplicated_percent]` and `[qc_thresholds.exome.relative.ZERO_CVG_TARGETS_PCT]`, change both occurrences of:

```toml
min_cohort = 50    # below this, MAD is too noisy; skip relative flagging
```

to:

```toml
min_samples = 50   # below this, MAD is too noisy; skip relative flagging
```

- [ ] **Step 5: Update any pre-existing test that still spells the old key**

Run: `grep -rn min_cohort src/ test/`

Every remaining hit in `test/test_check_multiqc.py` must be changed to `min_samples`. Hits in `src/align_genotype/qc_calibration/` are the old CLI modules deleted in Tasks 6b and 17 — leave them, they still import cleanly until then.

- [ ] **Step 6: Run the full check_multiqc suite**

Run: `uv run pytest test/test_check_multiqc.py -v`

Expected: PASS, all tests.

- [ ] **Step 7: Commit**

```bash
git add src/align_genotype/scripts/check_multiqc.py src/align_genotype/config_template.toml test/test_check_multiqc.py
git commit -m "refactor(qc): rename the relative min_cohort config key to min_samples

Under CPG Flow a 'cohort' is a collection of datasets and sequencing groups.
This key bounds the number of sequencing groups in one MultiQC run, which for
a DatasetStage is one dataset - so min_samples is what it actually means."
```

---

## Task 2: `settings.py` — read the `[qc_calibration]` config block

Replaces `spec.py`. The old module loaded an operator-authored TOML carrying candidate `fail`/`warn` values, a `gated` flag, a `reviewed` interlock and free-text `rationale`. None of that survives: the tool now derives candidates and writes them to a report rather than to config, so there is nothing to hand-author and nothing to sign off inside the tool.

**Files:**
- Create: `src/align_genotype/qc_calibration/settings.py`
- Test: `test/test_qc_calibration_settings.py`

- [ ] **Step 1: Write the failing test**

Create `test/test_qc_calibration_settings.py`:

```python
"""Unit tests for calibration settings read from the [qc_calibration] config block."""

import pytest

from cpg_utils import config as cpg_config

from align_genotype.qc_calibration import settings as settings_mod

_MISSING = object()

METRICS = {
    'MEDIAN_COVERAGE': {'direction': 'min', 'unit': 'x'},
    'reads_duplicated_percent': {'direction': 'max', 'unit': '%', 'relative': True},
}


@pytest.fixture
def patch_config(monkeypatch):
    """Wire config_retrieve to an in-memory nested dict."""

    def _apply(tree: dict) -> None:
        def config_retrieve(keys, default=_MISSING):  # noqa: ANN202
            node = tree
            for key in keys:
                if not isinstance(node, dict) or key not in node:
                    # The real config_retrieve raises when a key is missing and no
                    # default was passed. `load()` relies on that for
                    # workflow.sequencing_type, so the fake has to do it too.
                    if default is _MISSING:
                        raise cpg_config.ConfigError(f'missing config key: {list(keys)}')
                    return default
                node = node[key]
            return node

        monkeypatch.setattr(settings_mod.config, 'config_retrieve', config_retrieve)

    return _apply


def test_load_reads_metrics_for_the_runs_sequencing_type(patch_config):
    patch_config(
        {
            'workflow': {'sequencing_type': 'genome'},
            'qc_calibration': {'genome': {'metrics': METRICS}, 'exome': {'metrics': {}}},
        },
    )
    loaded = settings_mod.load()
    assert loaded.seq_type == 'genome'
    assert loaded.metric_keys == ('MEDIAN_COVERAGE', 'reads_duplicated_percent')
    assert loaded.metric('MEDIAN_COVERAGE').unit == 'x'
    assert loaded.relative_metrics == (loaded.metric('reads_duplicated_percent'),)


def test_load_applies_documented_defaults(patch_config):
    patch_config(
        {'workflow': {'sequencing_type': 'genome'}, 'qc_calibration': {'genome': {'metrics': METRICS}}},
    )
    loaded = settings_mod.load()
    assert loaded.k == pytest.approx(3.5)
    assert loaded.min_samples == 50
    assert loaded.bars.max_warn_rate == pytest.approx(0.10)
    assert loaded.bars.max_growth_churn == pytest.approx(0.02)
    assert loaded.bars.max_merge_churn == pytest.approx(0.05)


def test_load_honours_overrides(patch_config):
    patch_config(
        {
            'workflow': {'sequencing_type': 'exome'},
            'qc_calibration': {
                'k': 3.0,
                'min_samples': 25,
                'max_merge_churn': 0.5,
                'exome': {'metrics': METRICS},
            },
        },
    )
    loaded = settings_mod.load()
    assert loaded.k == pytest.approx(3.0)
    assert loaded.min_samples == 25
    assert loaded.bars.max_merge_churn == pytest.approx(0.5)


def test_load_without_metrics_for_this_seq_type_is_an_error(patch_config):
    patch_config(
        {'workflow': {'sequencing_type': 'exome'}, 'qc_calibration': {'genome': {'metrics': METRICS}}},
    )
    with pytest.raises(settings_mod.SettingsError, match="sequencing type 'exome'"):
        settings_mod.load()


def test_bad_direction_is_rejected():
    with pytest.raises(settings_mod.SettingsError, match='direction'):
        settings_mod.parse_metric('X', {'direction': 'up'})


def test_missing_direction_is_rejected():
    with pytest.raises(settings_mod.SettingsError, match='missing required key: direction'):
        settings_mod.parse_metric('X', {'unit': 'x'})


def test_bad_unit_is_rejected():
    with pytest.raises(settings_mod.SettingsError, match='unit'):
        settings_mod.parse_metric('X', {'direction': 'min', 'unit': 'furlongs'})


def test_quoted_boolean_relative_is_rejected():
    """`bool("false")` is True, so a quoted boolean must not be coerced."""
    with pytest.raises(settings_mod.SettingsError, match='relative must be true or false'):
        settings_mod.parse_metric('X', {'direction': 'min', 'relative': 'false'})


def test_metric_defaults_to_frac_and_not_relative():
    metric = settings_mod.parse_metric('X', {'direction': 'min'})
    assert (metric.unit, metric.relative) == ('frac', False)


def test_enabled_defaults_to_false(patch_config):
    patch_config({'workflow': {'sequencing_type': 'genome'}, 'qc_calibration': {}})
    assert settings_mod.enabled() is False


def test_enabled_reads_the_flag(patch_config):
    patch_config({'workflow': {'sequencing_type': 'genome'}, 'qc_calibration': {'enabled': True}})
    assert settings_mod.enabled() is True


def test_current_thresholds_reshapes_production_config(patch_config):
    patch_config(
        {
            'workflow': {'sequencing_type': 'genome'},
            'qc_thresholds': {
                'genome': {
                    'fail': {'min': {'MEDIAN_COVERAGE': 15}, 'max': {'FREEMIX': 0.04}},
                    'warn': {'min': {'MEDIAN_COVERAGE': 25}},
                },
            },
        },
    )
    assert settings_mod.current_thresholds('genome') == {
        'MEDIAN_COVERAGE': {'fail': 15, 'warn': 25},
        'FREEMIX': {'fail': 0.04},
    }


def test_enabled_rejects_a_quoted_boolean(patch_config):
    """The highest-consequence setting in the module: a truthy string would turn every
    production run into a job-per-dataset calibration run."""
    patch_config({'workflow': {'sequencing_type': 'genome'}, 'qc_calibration': {'enabled': 'false'}})
    with pytest.raises(settings_mod.SettingsError, match='must be true or false'):
        settings_mod.enabled()


def test_load_rejects_a_non_numeric_k(patch_config):
    patch_config(
        {
            'workflow': {'sequencing_type': 'genome'},
            'qc_calibration': {'k': 'abc', 'genome': {'metrics': METRICS}},
        },
    )
    with pytest.raises(settings_mod.SettingsError, match='qc_calibration.k'):
        settings_mod.load()


def test_load_rejects_a_boolean_min_samples(patch_config):
    """`int(True)` is 1, which would silently enable relative tiers on tiny datasets."""
    patch_config(
        {
            'workflow': {'sequencing_type': 'genome'},
            'qc_calibration': {'min_samples': True, 'genome': {'metrics': METRICS}},
        },
    )
    with pytest.raises(settings_mod.SettingsError, match='must be an integer'):
        settings_mod.load()


def test_load_without_a_sequencing_type_raises(patch_config):
    """No usable config at all is a different problem from a bad [qc_calibration] block,
    so cpg-utils' own error is left to propagate rather than rewrapped."""
    patch_config({'qc_calibration': {'genome': {'metrics': METRICS}}})
    with pytest.raises(cpg_config.ConfigError):
        settings_mod.load()
```

`current_thresholds` delegates to `check_multiqc.load_thresholds`, which calls `config.config_retrieve` on the `check_multiqc` module object — so that test needs the patch applied there too. Add to the fixture, after the `settings_mod.config` line:

```python
        monkeypatch.setattr(settings_mod.check_multiqc.config, 'config_retrieve', config_retrieve)
```

- [ ] **Step 2: Run the tests to verify they fail**

Run: `uv run pytest test/test_qc_calibration_settings.py -v`

Expected: FAIL — `ModuleNotFoundError: No module named 'align_genotype.qc_calibration.settings'`.

- [ ] **Step 3: Write the implementation**

Create `src/align_genotype/qc_calibration/settings.py`:

```python
"""Calibration settings, read from the `[qc_calibration]` config block.

Metric lists are per sequencing type because the Picard module differs: genome uses
`CollectWgsMetrics`, exome `CollectHsMetrics`. That is why the candidate metric set is a
config parameter and not a Python constant - an exome spec copied onto a genome run would
gate keys that do not exist.
"""

from dataclasses import dataclass
from typing import Any

from cpg_utils import config

from align_genotype.scripts import check_multiqc

DIRECTIONS = ('min', 'max')
UNITS = ('x', 'frac', '%')

_METRIC_KEYS = frozenset({'direction', 'unit', 'relative'})


class SettingsError(ValueError):
    """The `[qc_calibration]` config block is missing something or is inconsistent."""


@dataclass(frozen=True)
class MetricSpec:
    """One metric's direction, display unit and whether to evaluate a relative tier.

    `direction`: 'min' = higher is better, so a value below the threshold is flagged;
    'max' is the mirror. `unit` affects rounding and display only.
    """

    key: str
    direction: str
    unit: str = 'frac'
    relative: bool = False


@dataclass(frozen=True)
class Bars:
    """The advisory bars a relative tier's verdict is judged against.

    Growth and merge churn have separate bars because they model different things.
    Growth - a before-slice re-scored against the threshold the whole dataset produces -
    is a forecast: datasets accrete sequencing groups over time, which is what a shipped
    tier actually faces. Merge - two whole projects pooled into one run - is a stress
    test, not a prediction of anything scheduled. Hence the looser merge bar.
    """

    max_warn_rate: float = 0.10
    max_growth_churn: float = 0.02
    max_merge_churn: float = 0.05


@dataclass(frozen=True)
class CalibrationSettings:
    seq_type: str
    metrics: tuple[MetricSpec, ...]
    k: float = 3.5
    min_samples: int = 50
    bars: Bars = Bars()

    @property
    def metric_keys(self) -> tuple[str, ...]:
        return tuple(m.key for m in self.metrics)

    @property
    def relative_metrics(self) -> tuple[MetricSpec, ...]:
        return tuple(m for m in self.metrics if m.relative)

    def metric(self, key: str) -> MetricSpec:
        for m in self.metrics:
            if m.key == key:
                return m
        raise KeyError(f'{key!r} is not a configured calibration metric; known: {list(self.metric_keys)}')


def parse_metric(key: str, raw: Any) -> MetricSpec:
    """Validate one `[qc_calibration.<seq_type>.metrics.<KEY>]` table."""
    if not isinstance(raw, dict):
        raise SettingsError(f'metric {key!r}: must be a table of direction/unit/relative keys')
    if unknown := sorted(set(raw) - _METRIC_KEYS):
        raise SettingsError(f'metric {key!r}: unknown key(s) {unknown}; expected {sorted(_METRIC_KEYS)}')
    if 'direction' not in raw:
        raise SettingsError(f'metric {key!r}: missing required key: direction')
    direction = raw['direction']
    if direction not in DIRECTIONS:
        raise SettingsError(f"metric {key!r}: direction must be 'min' or 'max', got {direction!r}")
    unit = raw.get('unit', 'frac')
    if unit not in UNITS:
        raise SettingsError(f'metric {key!r}: unit must be one of {list(UNITS)}, got {unit!r}')
    relative = raw.get('relative', False)
    # `bool('false')` is True, so a quoted boolean must be rejected rather than coerced.
    if not isinstance(relative, bool):
        raise SettingsError(f'metric {key!r}: relative must be true or false, got {relative!r}')
    return MetricSpec(key=key, direction=direction, unit=unit, relative=relative)


def _require_bool(key: str, value: Any) -> bool:
    """Reject anything that is not already a boolean rather than coercing it.

    `bool('false')` is True. For `enabled` in particular that typo would turn every
    ordinary production run into a job-per-dataset calibration run.
    """
    if not isinstance(value, bool):
        raise SettingsError(f'qc_calibration.{key} must be true or false, got {value!r}')
    return value


def _require_number(key: str, value: Any) -> float:
    """Reject a non-numeric setting with a message naming the key.

    `bool` subclasses `int`, so it is excluded explicitly. A bare `float()` here would
    raise `ValueError: could not convert string to float: 'abc'`, which never mentions
    which setting was wrong.
    """
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        raise SettingsError(f'qc_calibration.{key} must be a number, got {value!r}')
    return float(value)


def _require_int(key: str, value: Any) -> int:
    if isinstance(value, bool) or not isinstance(value, int):
        raise SettingsError(f'qc_calibration.{key} must be an integer, got {value!r}')
    return value


def enabled() -> bool:
    """Whether the calibration stages should queue any jobs at all.

    Defaults to False so that the stages sit inert in the DAG of an ordinary production
    run. They carry no `required_stages` dependency, so without this they would queue a
    job per dataset and register Metamist analyses on every invocation. Validated rather
    than coerced - see `_require_bool`.
    """
    return _require_bool('enabled', config.config_retrieve(['qc_calibration', 'enabled'], False))


def load() -> CalibrationSettings:
    """Read the settings for the run's sequencing type."""
    seq_type = config.config_retrieve(['workflow', 'sequencing_type'])
    raw_metrics = config.config_retrieve(['qc_calibration', seq_type, 'metrics'], {})
    if not raw_metrics:
        raise SettingsError(
            f'no calibration metrics configured for sequencing type {seq_type!r}; '
            f'add [qc_calibration.{seq_type}.metrics.<KEY>] tables',
        )
    return CalibrationSettings(
        seq_type=seq_type,
        metrics=tuple(parse_metric(key, value) for key, value in raw_metrics.items()),
        k=_require_number('k', config.config_retrieve(['qc_calibration', 'k'], 3.5)),
        min_samples=_require_int('min_samples', config.config_retrieve(['qc_calibration', 'min_samples'], 50)),
        bars=Bars(
            max_warn_rate=_require_number(
                'max_warn_rate', config.config_retrieve(['qc_calibration', 'max_warn_rate'], 0.10)
            ),
            max_growth_churn=_require_number(
                'max_growth_churn', config.config_retrieve(['qc_calibration', 'max_growth_churn'], 0.02)
            ),
            max_merge_churn=_require_number(
                'max_merge_churn', config.config_retrieve(['qc_calibration', 'max_merge_churn'], 0.05)
            ),
        ),
    )


def current_thresholds(seq_type: str) -> dict[str, dict[str, float]]:
    """The shipped absolute thresholds, as `{metric: {severity: value}}`.

    Read through `check_multiqc.load_thresholds` rather than a second config reader, so
    the "current" column in the report is exactly what production enforces today.
    """
    by_direction = check_multiqc.load_thresholds(seq_type)
    return {metric: dict(tiers) for by_metric in by_direction.values() for metric, tiers in by_metric.items()}
```

- [ ] **Step 4: Run the tests to verify they pass**

Run: `uv run pytest test/test_qc_calibration_settings.py -v`

Expected: PASS, 18 tests.

- [ ] **Step 5: Commit**

```bash
git add src/align_genotype/qc_calibration/settings.py test/test_qc_calibration_settings.py
git commit -m "feat(qc_calibration): read calibration settings from config

Replaces the operator-authored spec TOML. The tool now derives candidate
thresholds rather than scoring hand-written ones, so gated/fail/warn/reviewed
/rationale all go; what remains is direction, unit and whether to evaluate a
relative tier."
```

---

## Task 3: Add the `[qc_calibration]` block to `config_template.toml`

**Files:**
- Modify: `src/align_genotype/config_template.toml`
- Test: `test/test_qc_calibration_settings.py`

- [ ] **Step 1: Write the failing test**

Append to `test/test_qc_calibration_settings.py`:

```python
def test_config_template_ships_a_parseable_calibration_block():
    """Every metric in the shipped template must pass validation, for both seq types."""
    import sys  # noqa: PLC0415
    from pathlib import Path  # noqa: PLC0415

    if sys.version_info >= (3, 11):
        import tomllib  # noqa: PLC0415
    else:
        import tomli as tomllib  # noqa: PLC0415

    template = Path('src/align_genotype/config_template.toml')
    parsed = tomllib.loads(template.read_text())
    block = parsed['qc_calibration']

    assert block['enabled'] is False, 'must ship disabled so production runs stay inert'
    for seq_type in ('genome', 'exome'):
        metrics = block[seq_type]['metrics']
        assert metrics, f'no calibration metrics shipped for {seq_type}'
        for key, raw in metrics.items():
            settings_mod.parse_metric(key, raw)

    # Every metric marked relative must have a shipped absolute fail gate behind it:
    # relative flagging is warn-only, so without one a uniformly poor dataset is ungated.
    for seq_type in ('genome', 'exome'):
        fails = set(parsed['qc_thresholds'][seq_type].get('fail', {}).get('min', {})) | set(
            parsed['qc_thresholds'][seq_type].get('fail', {}).get('max', {}),
        )
        for key, raw in block[seq_type]['metrics'].items():
            if raw.get('relative'):
                assert key in fails, f'{seq_type} {key} is relative but has no absolute fail gate'
```

- [ ] **Step 2: Run the test to verify it fails**

Run: `uv run pytest test/test_qc_calibration_settings.py::test_config_template_ships_a_parseable_calibration_block -v`

Expected: FAIL with `KeyError: 'qc_calibration'`.

- [ ] **Step 3: Add the config block**

In `src/align_genotype/config_template.toml`, immediately **before** the `[references.broad]` line (find it with `grep -n '^\[references.broad\]' src/align_genotype/config_template.toml`), insert:

```toml
# QC threshold calibration. Derives candidate warn/fail thresholds and dataset-relative
# (MAD) tiers from every dataset's latest CramMultiQC report, and writes an HTML report
# to the analysis dataset's web bucket. See qc_calibration/README.md.
#
# `enabled` is false so that these stages sit inert in an ordinary production run: they
# carry no stage dependencies, so without it they would queue a job per dataset every
# time. A calibration run sets it true and passes
#   only_stages = ['QcCalibrationDatasetMetrics', 'QcCalibrationReport']
[qc_calibration]
enabled = false
k = 3.5                  # Iglewicz-Hoaglin outlier constant for relative tiers
min_samples = 50         # below this, MAD is too noisy to derive a relative tier
max_warn_rate = 0.10     # advisory bars behind the RECOMMEND / REJECT verdict
max_growth_churn = 0.02  # a forecast: datasets accrete sequencing groups over time
max_merge_churn = 0.05   # a stress test: nothing schedules two projects into one run
extract_storage = "20Gi" # disk for one localised multiqc_data.json

# Genome candidates - Picard CollectWgsMetrics, plus samtools and verifybamid.
# direction: min = higher is better. unit: x | % | frac, display and rounding only.
[qc_calibration.genome.metrics.MEDIAN_COVERAGE]
direction = "min"
unit = "x"

[qc_calibration.genome.metrics.PCT_20X]
direction = "min"
unit = "frac"

[qc_calibration.genome.metrics.FREEMIX]
direction = "max"
unit = "frac"

[qc_calibration.genome.metrics.reads_mapped_percent]
direction = "min"
unit = "%"

[qc_calibration.genome.metrics.reads_properly_paired_percent]
direction = "min"
unit = "%"

# Duplication is strongly library-prep dependent on WGS, so its warn tier is evaluated
# dataset-relatively rather than as a fixed line. The absolute fail gate stays.
[qc_calibration.genome.metrics.reads_duplicated_percent]
direction = "max"
unit = "%"
relative = true

# Exome candidates - Picard CollectHsMetrics. There are no *_TARGET_* or FOLD_80 keys on
# a genome report, and no MEDIAN_COVERAGE on an exome one, so the two lists are disjoint
# where the capture matters.
[qc_calibration.exome.metrics.MEAN_TARGET_COVERAGE]
direction = "min"
unit = "x"

[qc_calibration.exome.metrics.PCT_TARGET_BASES_20X]
direction = "min"
unit = "frac"

[qc_calibration.exome.metrics.FOLD_80_BASE_PENALTY]
direction = "max"
unit = "frac"

[qc_calibration.exome.metrics.FREEMIX]
direction = "max"
unit = "frac"

[qc_calibration.exome.metrics.reads_mapped_percent]
direction = "min"
unit = "%"

[qc_calibration.exome.metrics.reads_duplicated_percent]
direction = "max"
unit = "%"

# Zero-coverage target rate is capture-kit dependent, so its warn tier is dataset-relative.
[qc_calibration.exome.metrics.ZERO_CVG_TARGETS_PCT]
direction = "max"
unit = "frac"
relative = true
```

- [ ] **Step 4: Strip the local-scratch citations from the `qc_thresholds` comments**

Still in `config_template.toml`, three edits. Find them with
`grep -n 'testing_scripts\|mad_relative_prototype\|verify_relative_6cohorts\|qc_metric_distributions' src/align_genotype/config_template.toml`.

Change:

```
# stop; a sample caught there is not double-flagged here. Validated across all 10 WGS
# cohorts via the production relative_flags() path (testing_scripts/mad_relative_prototype.py,
# verify_relative_6cohorts.py).
```

to:

```
# stop; a sample caught there is not double-flagged here.
```

Change:

```
# Calibrated against two real WES cohorts (647- and 522-sample MultiQC runs; see
# testing_scripts/qc_metric_distributions.py + qc_threshold_flagrates.py). Values
# are set so a healthy cohort flags at ~0% fail / single-digit % warn, while a
# genuinely poor batch still surfaces. These metrics are capture-kit dependent -
# re-run those scripts on a known-good batch when onboarding a new kit.
```

to:

```
# Values are set so a healthy dataset flags at ~0% fail / single-digit % warn, while
# a genuinely poor batch still surfaces. These metrics are capture-kit dependent -
# re-run the QC calibration workflow on a known-good batch when onboarding a new kit.
```

Then replace the remaining bare uses of "cohort" in the `qc_thresholds` comments with "dataset", since under CPG Flow each MultiQC run covers one dataset. Find them with
`grep -n 'cohort' src/align_genotype/config_template.toml` and change every hit inside a `qc_thresholds` comment (e.g. "Cohort medians ~32-37x" → "Dataset medians ~32-37x", "~25-37% of legitimately higher-dup preps in the two highest-duplication cohorts" → "... datasets").

- [ ] **Step 5: Run the tests to verify they pass**

Run: `uv run pytest test/test_qc_calibration_settings.py -v`

Expected: PASS, 19 tests (the 18 from Task 2, plus this one).

- [ ] **Step 6: Verify no scratch references remain**

Run: `grep -rn 'testing_scripts' src/`

Expected: no output.

- [ ] **Step 7: Commit**

```bash
git add src/align_genotype/config_template.toml test/test_qc_calibration_settings.py
git commit -m "feat(qc_calibration): ship the calibration metric lists in the config template

Genome and exome candidate sets are disjoint where the Picard module differs, so
they are config rather than constants. Ships disabled. Also drops the local
scratch-script citations from the qc_thresholds comments and switches their
cohort/dataset wording to CPG Flow's."
```

---

## Task 4: `values.py` — the per-dataset values file

Replaces `cache.py`. Two changes of substance. First, the sequencing group ID is kept alongside every value, so value counts and sequencing-group counts are both exact — the old cache stored bare floats and therefore could not answer "what is the per-sequencing-group warn rate", which matters because MultiQC 1.33 can split one tool across sections and a sequencing group then contributes one value per section. Second, the file is written by a Hail Batch job to a local path that Hail uploads, so the old `mkstemp`-then-`replace` dance and the `CloudPath` branch are gone.

**Files:**
- Create: `src/align_genotype/qc_calibration/values.py`
- Test: `test/test_qc_calibration_values.py`

- [ ] **Step 1: Write the failing test**

Create `test/test_qc_calibration_values.py`:

```python
"""Unit tests for the per-dataset values file."""

import json

import numpy as np
import pytest

from align_genotype.qc_calibration import values as values_mod


def make_values(**overrides) -> values_mod.DatasetValues:
    defaults = {
        'dataset': 'dataset-a',
        'seq_type': 'genome',
        'analysis_id': 42,
        'timestamp': '2026-06-01T00:00:00',
        'uri': 'gs://bucket/multiqc_data.json',
        'multiqc_version': '1.33',
        'generated': '2026-08-12T00:00:00',
        'n_sequencing_groups': 3,
        'section_sizes': {'picard_1': 3},
        'metrics': {
            'MEDIAN_COVERAGE': values_mod.MetricValues(
                entries=(('picard_1', 'CPG1', 30.0), ('picard_1', 'CPG2', 34.0), ('picard_1', 'CPG3', 38.0)),
                n_dropped=1,
            ),
        },
    }
    return values_mod.DatasetValues(**{**defaults, **overrides})


def test_metric_values_exposes_an_array_in_entry_order():
    metric = make_values().metrics['MEDIAN_COVERAGE']
    np.testing.assert_array_equal(metric.array, np.array([30.0, 34.0, 38.0]))


def test_metric_values_counts_values_and_sequencing_groups_separately():
    """One sequencing group in two sections yields two values but one group."""
    metric = values_mod.MetricValues(
        entries=(('picard_1', 'CPG1', 1.0), ('picard_4', 'CPG1', 1.0), ('picard_1', 'CPG2', 2.0)),
        n_dropped=0,
    )
    assert metric.n_values == 3
    assert metric.n_sequencing_groups == 2
    assert metric.duplicated is True
    assert metric.sections == ('picard_1', 'picard_4')


def test_metric_values_not_duplicated_when_one_section():
    assert make_values().metrics['MEDIAN_COVERAGE'].duplicated is False


def test_empty_metric_values_gives_an_empty_array():
    metric = values_mod.MetricValues(entries=(), n_dropped=0)
    assert metric.array.size == 0
    assert metric.n_values == 0
    assert metric.n_sequencing_groups == 0
    assert metric.sections == ()


def test_array_filters_non_finite_defensively():
    """extract already drops these; the file is JSON and could be hand-edited."""
    metric = values_mod.MetricValues(
        entries=(('s', 'CPG1', 1.0), ('s', 'CPG2', float('nan')), ('s', 'CPG3', float('inf'))),
        n_dropped=0,
    )
    np.testing.assert_array_equal(metric.array, np.array([1.0]))


def test_round_trips_through_json(tmp_path):
    path = tmp_path / 'values.json'
    original = make_values()
    values_mod.save(original, path)
    assert values_mod.load(path) == original


def test_saved_json_is_readable_and_shaped_as_documented(tmp_path):
    path = tmp_path / 'values.json'
    values_mod.save(make_values(), path)
    payload = json.loads(path.read_text())
    assert payload['dataset'] == 'dataset-a'
    assert payload['analysis_id'] == 42
    assert payload['metrics']['MEDIAN_COVERAGE']['n_dropped'] == 1
    assert payload['metrics']['MEDIAN_COVERAGE']['entries'][0] == ['picard_1', 'CPG1', 30.0]


def test_save_refuses_a_non_finite_value(tmp_path):
    """A NaN reaching the file would poison every percentile downstream."""
    broken = make_values(
        metrics={'X': values_mod.MetricValues(entries=(('s', 'CPG1', float('nan')),), n_dropped=0)},
    )
    with pytest.raises(values_mod.ValuesError, match='non-finite'):
        values_mod.save(broken, tmp_path / 'values.json')


def test_load_names_the_file_when_the_shape_is_wrong(tmp_path):
    path = tmp_path / 'values.json'
    path.write_text('{"dataset": "a"}')
    with pytest.raises(values_mod.ValuesError, match=str(path)):
        values_mod.load(path)


def test_load_names_the_file_when_the_json_is_malformed(tmp_path):
    path = tmp_path / 'values.json'
    path.write_text('{not json')
    with pytest.raises(values_mod.ValuesError, match=str(path)):
        values_mod.load(path)
```

- [ ] **Step 2: Run the tests to verify they fail**

Run: `uv run pytest test/test_qc_calibration_values.py -v`

Expected: FAIL — `ModuleNotFoundError: No module named 'align_genotype.qc_calibration.values'`.

- [ ] **Step 3: Write the implementation**

Create `src/align_genotype/qc_calibration/values.py`:

```python
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

    `allow_nan=False` is the guard that matters: a non-finite value written here would
    reach every percentile and MAD in the report stage. Serialising before opening the
    target means such a value raises with nothing written.
    """
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
    for key, metric in values.metrics.items():
        for _, sg, value in metric.entries:
            if not math.isfinite(value):
                raise ValuesError(f'{values.dataset}: metric {key!r} sequencing group {sg!r} has a non-finite value')
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
            analysis_id=int(raw['analysis_id']),
            timestamp=raw['timestamp'],
            uri=raw['uri'],
            multiqc_version=raw['multiqc_version'],
            generated=raw['generated'],
            n_sequencing_groups=int(raw['n_sequencing_groups']),
            section_sizes=dict(raw['section_sizes']),
            metrics={
                key: MetricValues(
                    entries=tuple((section, sg, float(value)) for section, sg, value in body['entries']),
                    n_dropped=int(body['n_dropped']),
                )
                for key, body in raw['metrics'].items()
            },
        )
    except (KeyError, TypeError, ValueError, AttributeError) as exc:
        raise ValuesError(f'{path}: not a usable values file - {exc}') from exc
```

`json.JSONDecodeError` subclasses `ValueError`, so it is caught by the clause above.

- [ ] **Step 4: Run the tests to verify they pass**

Run: `uv run pytest test/test_qc_calibration_values.py -v`

Expected: PASS, 10 tests.

- [ ] **Step 5: Commit**

```bash
git add src/align_genotype/qc_calibration/values.py test/test_qc_calibration_values.py
git commit -m "feat(qc_calibration): per-dataset values file keyed by sequencing group

Replaces the multi-dataset value cache. Keeping the sequencing group alongside
each value makes value and group counts both exact, which the old bare-float
cache could not do when MultiQC splits one tool across general-stats sections.
Hail writes the file, so the atomic-rename and CloudPath branches go."
```

---

## Task 5: `extract.py` — one MultiQC document to one `DatasetValues`

Replaces `collect.py`. The old module looped over a manifest, called `gc.collect()` between reports and wrapped each in a broad `except Exception` so one bad dataset did not discard the successes behind it. None of that is needed: one Hail job handles one report, and Hail isolates per-job failure.

**Files:**
- Create: `src/align_genotype/qc_calibration/extract.py`
- Test: `test/test_qc_calibration_extract.py`

- [ ] **Step 1: Write the failing test**

Create `test/test_qc_calibration_extract.py`:

```python
"""Unit tests for extracting one MultiQC report into a values file."""

import numpy as np
import pytest

from align_genotype.qc_calibration import extract as extract_mod
from align_genotype.qc_calibration import settings as settings_mod

SETTINGS = settings_mod.CalibrationSettings(
    seq_type='genome',
    metrics=(
        settings_mod.MetricSpec(key='MEDIAN_COVERAGE', direction='min', unit='x'),
        settings_mod.MetricSpec(key='FREEMIX', direction='max', unit='frac'),
    ),
)

PROVENANCE = {
    'dataset': 'dataset-a',
    'analysis_id': 42,
    'timestamp': '2026-06-01T00:00:00',
    'uri': 'gs://bucket/multiqc_data.json',
    'generated': '2026-08-12T00:00:00',
}


def document(general_stats, version='1.33') -> dict:
    return {'config_version': version, 'report_general_stats_data': general_stats}


def test_extracts_values_from_the_v133_dict_shape():
    doc = document(
        {
            'picard_1': {'CPG1': {'MEDIAN_COVERAGE': 30.0}, 'CPG2': {'MEDIAN_COVERAGE': 34.0}},
            'verifybamid': {'CPG1': {'FREEMIX': 0.001}},
        },
    )
    result = extract_mod.extract(doc, SETTINGS, **PROVENANCE)
    assert result.dataset == 'dataset-a'
    assert result.seq_type == 'genome'
    assert result.multiqc_version == '1.33'
    assert result.n_sequencing_groups == 2
    assert result.section_sizes == {'picard_1': 2, 'verifybamid': 1}
    np.testing.assert_array_equal(result.metric('MEDIAN_COVERAGE').array, np.array([30.0, 34.0]))
    np.testing.assert_array_equal(result.metric('FREEMIX').array, np.array([0.001]))


def test_extracts_values_from_the_v114_list_shape():
    """v1.14 stores general stats as a positional list; both shapes exist in real reports."""
    doc = document([{'CPG1': {'MEDIAN_COVERAGE': 30.0}}, {'CPG1': {'FREEMIX': 0.001}}], version='1.14')
    result = extract_mod.extract(doc, SETTINGS, **PROVENANCE)
    assert sorted(result.section_sizes) == ['section_0', 'section_1']
    np.testing.assert_array_equal(result.metric('MEDIAN_COVERAGE').array, np.array([30.0]))


def test_strips_the_rich_id_suffix_to_get_the_sequencing_group():
    """MultiQC runs with --replace-names, so sample keys can read CPG1|EXTID."""
    doc = document({'picard_1': {'CPG1|EXT1': {'MEDIAN_COVERAGE': 30.0}}})
    result = extract_mod.extract(doc, SETTINGS, **PROVENANCE)
    assert result.metric('MEDIAN_COVERAGE').entries == (('picard_1', 'CPG1', 30.0),)


def test_one_group_in_two_sections_yields_two_values_and_one_group():
    doc = document(
        {
            'picard_1': {'CPG1': {'MEDIAN_COVERAGE': 30.0}},
            'picard_4': {'CPG1': {'MEDIAN_COVERAGE': 30.0}},
        },
    )
    metric = extract_mod.extract(doc, SETTINGS, **PROVENANCE).metric('MEDIAN_COVERAGE')
    assert (metric.n_values, metric.n_sequencing_groups) == (2, 1)


def test_picard_question_mark_placeholder_is_dropped_and_counted():
    doc = document({'picard_1': {'CPG1': {'MEDIAN_COVERAGE': '?'}, 'CPG2': {'MEDIAN_COVERAGE': 34.0}}})
    metric = extract_mod.extract(doc, SETTINGS, **PROVENANCE).metric('MEDIAN_COVERAGE')
    np.testing.assert_array_equal(metric.array, np.array([34.0]))
    assert metric.n_dropped == 1


def test_non_finite_value_is_dropped_and_counted():
    doc = document({'picard_1': {'CPG1': {'MEDIAN_COVERAGE': float('nan')}, 'CPG2': {'MEDIAN_COVERAGE': 34.0}}})
    metric = extract_mod.extract(doc, SETTINGS, **PROVENANCE).metric('MEDIAN_COVERAGE')
    np.testing.assert_array_equal(metric.array, np.array([34.0]))
    assert metric.n_dropped == 1


def test_numeric_strings_are_coerced():
    doc = document({'picard_1': {'CPG1': {'MEDIAN_COVERAGE': '30.5'}}})
    metric = extract_mod.extract(doc, SETTINGS, **PROVENANCE).metric('MEDIAN_COVERAGE')
    np.testing.assert_array_equal(metric.array, np.array([30.5]))
    assert metric.n_dropped == 0


def test_a_metric_absent_from_the_report_is_present_but_empty():
    """Absence must be recorded, not omitted - the report banner is built from this."""
    doc = document({'picard_1': {'CPG1': {'MEDIAN_COVERAGE': 30.0}}})
    result = extract_mod.extract(doc, SETTINGS, **PROVENANCE)
    assert 'FREEMIX' in result.metrics
    assert result.metric('FREEMIX').entries == ()


def test_a_document_that_is_not_an_object_is_an_error():
    with pytest.raises(extract_mod.ExtractError, match='not an object'):
        extract_mod.extract([], SETTINGS, **PROVENANCE)


def test_no_usable_general_stats_is_an_error():
    """Refusing beats reporting a clean extraction on a report we could not read."""
    with pytest.raises(extract_mod.ExtractError, match='no usable report_general_stats_data'):
        extract_mod.extract(document(None), SETTINGS, **PROVENANCE)


def test_empty_general_stats_is_an_error_naming_the_dataset_and_uri():
    with pytest.raises(extract_mod.ExtractError, match='dataset-a.*gs://bucket'):
        extract_mod.extract(document({}), SETTINGS, **PROVENANCE)
```

- [ ] **Step 2: Run the tests to verify they fail**

Run: `uv run pytest test/test_qc_calibration_extract.py -v`

Expected: FAIL — `ModuleNotFoundError: No module named 'align_genotype.qc_calibration.extract'`.

- [ ] **Step 3: Write the implementation**

Create `src/align_genotype/qc_calibration/extract.py`:

```python
"""Distil one MultiQC report into one dataset's values file.

Extraction goes through `check_multiqc.normalise_sections` and
`check_multiqc.gather_metric_values`, so calibration sees precisely what enforcement
sees. A metric configured but absent from the report is recorded present-and-empty
rather than omitted: that distinction is what the report's presence matrix is built
from, and a silently absent metric is how a gate goes inert.
"""

import math
from typing import Any

from align_genotype.qc_calibration.settings import CalibrationSettings
from align_genotype.qc_calibration.values import DatasetValues, MetricValues
from align_genotype.scripts import check_multiqc


class ExtractError(RuntimeError):
    """A MultiQC report could not be read, or holds nothing to extract."""


def extract(
    document: Any,
    settings: CalibrationSettings,
    *,
    dataset: str,
    analysis_id: int,
    timestamp: str,
    uri: str,
    generated: str,
) -> DatasetValues:
    """Extract every configured metric from one parsed MultiQC document."""
    # `[]`, `null`, `42` and `"str"` are all valid JSON, so a truncated or wrong-file URI
    # can parse cleanly and then fail on `.get` with a bare AttributeError naming neither
    # the dataset nor the path.
    if not isinstance(document, dict):
        raise ExtractError(f'{dataset}: report is a {type(document).__name__}, not an object, in {uri}')

    version = str(document.get('config_version', 'unknown'))
    sections = check_multiqc.normalise_sections(document.get('report_general_stats_data'))
    if not sections:
        raise ExtractError(f'{dataset}: no usable report_general_stats_data (multiqc {version}) in {uri}')

    metrics: dict[str, MetricValues] = {}
    for metric in settings.metrics:
        entries, n_non_numeric = check_multiqc.gather_metric_values(sections, metric.key)
        finite = tuple(
            # Production derives a sequencing group ID the same way: MultiQC runs with
            # --replace-names against the dataset's rich ID map, so a sample key can read
            # `CPG1|EXTID`.
            (section, sample.split('|', 1)[0], value)
            for section, sample, value in entries
            if math.isfinite(value)
        )
        # Two disjoint kinds of loss - values float() refused, and values it accepted that
        # came back nan/inf - so summing them cannot double-count.
        metrics[metric.key] = MetricValues(
            entries=finite,
            n_dropped=n_non_numeric + (len(entries) - len(finite)),
        )

    return DatasetValues(
        dataset=dataset,
        seq_type=settings.seq_type,
        analysis_id=analysis_id,
        timestamp=timestamp,
        uri=uri,
        multiqc_version=version,
        generated=generated,
        n_sequencing_groups=len(
            {sample.split('|', 1)[0] for section in sections.values() for sample in section},
        ),
        section_sizes={name: len(section) for name, section in sections.items()},
        metrics=metrics,
    )
```

- [ ] **Step 4: Run the tests to verify they pass**

Run: `uv run pytest test/test_qc_calibration_extract.py -v`

Expected: PASS, 11 tests.

- [ ] **Step 5: Commit**

```bash
git add src/align_genotype/qc_calibration/extract.py test/test_qc_calibration_extract.py
git commit -m "feat(qc_calibration): extract one MultiQC report per dataset

One report per call, so the gc.collect() loop and the broad per-dataset except
that kept one bad URI from discarding the successes behind it both go - Hail
isolates per-job failure. Sequencing group IDs are stripped of their rich-ID
suffix exactly as production does."
```

---

## Task 6: Retire the "cohort" wording in `stats.py`

`stats.py` needs no logic change — percentiles, flag rates and the churn simulation are
already pure functions over arrays. Only the wording is wrong.

**Files:**
- Modify: `src/align_genotype/qc_calibration/stats.py`
- Modify: `test/test_qc_calibration_stats.py`

- [ ] **Step 1: Confirm the suite is green before touching anything**

Run: `uv run pytest test/test_qc_calibration_stats.py -v`

Expected: PASS, 20 tests.

- [ ] **Step 2: Rewrite the module docstring and the two comments**

In `src/align_genotype/qc_calibration/stats.py`, replace the module docstring with:

```python
"""Numeric analysis over extracted values - percentiles, flag rates, dataset-growth churn.

No I/O and no config: everything here is a pure function of arrays already in memory.
"""
```

Replace the `FAIL_RATE_LIMIT` / `WARN_RATE_LIMIT` comment with:

```python
# A healthy dataset should sit near 0% fail and single-digit % warn. Beyond these, a
# candidate threshold is marked for a second look - advice, not a rejection: "healthy
# dataset" is a human judgement, not a computable property.
```

In `ChurnResult`, replace the class docstring with:

```python
    """How a dataset-relative threshold moved, and who changed status because of it."""
```

In `churn`, replace the docstring with:

```python
    """Re-score the *initial* values against the *grown* set's threshold.

    `flips` counts values whose flag status changes purely because the dataset grew -
    each one a spurious "updated" flag in the database, which is the cost that decides
    whether a dataset-relative tier is safe to adopt. Returns None when either set has a
    degenerate (zero) MAD, since no threshold exists to compare.

    Thresholds are rounded to 4 dp exactly as production does, so the simulation measures
    the churn operators would actually see rather than sub-0.0001 jitter.
    """
```

- [ ] **Step 3: Update the test file's wording**

Run: `grep -n cohort test/test_qc_calibration_stats.py`

Change every hit to `dataset`. These are docstrings and comments only; no assertion changes.

- [ ] **Step 4: Verify nothing broke**

Run: `uv run pytest test/test_qc_calibration_stats.py -v`

Expected: PASS, 20 tests.

- [ ] **Step 5: Verify no stray wording remains**

Run: `grep -n cohort src/align_genotype/qc_calibration/stats.py`

Expected: no output.

- [ ] **Step 6: Commit**

```bash
git add src/align_genotype/qc_calibration/stats.py test/test_qc_calibration_stats.py
git commit -m "refactor(qc_calibration): use dataset wording in stats

Comments and docstrings only - the arithmetic is unchanged."
```

---

## Task 6b: Delete the first half of the CLI surface

Deletion is split in two because the old modules import each other. `cli.py` imports
almost everything; `emit.py` imports `report`, `spec`, `cache` and `tomlio`; `dryrun.py`
imports `emit` and `manifest`. Rewriting `relative.py` in Task 7 removes its `spec` /
`cache` / `tomlio` imports, but the modules deleted here would still be importing the old
versions — so they have to go first, or the full suite is red from Task 7 to Task 17.

What stays for now: `spec.py`, `cache.py`, `tomlio.py` and `manifest.py`, because the
un-rewritten `discovery.py` still imports the last two. They go in Task 17, once Task 12
has rewritten it.

**Files:**
- Delete: `src/align_genotype/qc_calibration/{cli,dryrun,emit,report,suggest,collect}.py`
- Delete: `test/test_qc_calibration_{cli,dryrun,emit,report,suggest,collect}.py`
- Modify: `pyproject.toml`

- [ ] **Step 1: Delete the modules and their tests**

```bash
git rm src/align_genotype/qc_calibration/{cli,dryrun,emit,report,suggest,collect}.py
git rm test/test_qc_calibration_{cli,dryrun,emit,report,suggest,collect}.py
```

- [ ] **Step 2: Drop the console script**

In `pyproject.toml`, remove these three lines — the entry point now names a module that
no longer exists:

```toml
# the QC threshold calibration tool - operator-facing, run by hand off a checkout; the
# pipeline never invokes it
qc_calibrate = 'align_genotype.qc_calibration.cli:main'
```

- [ ] **Step 3: Verify nothing still imports them**

Run: `grep -rn 'qc_calibration import \(cli\|dryrun\|emit\|report\|suggest\|collect\)\|qc_calibration\.\(cli\|dryrun\|emit\|report\|suggest\|collect\)\|qc_calibrate' src/ test/ pyproject.toml`

Expected: no output. `report.py`'s `fmt_measure` was imported by `emit.py`; both are gone,
and the new `render.py` carries its own formatting filters.

- [ ] **Step 4: Run the whole suite**

Run: `uv run pytest -q`

Expected: PASS. The remaining old tests (`spec`, `cache`, `tomlio`, `manifest`,
`discovery`, `relative`, `stats`) still cover unchanged modules, alongside the new
`settings`, `values` and `extract` suites.

- [ ] **Step 5: Commit**

```bash
git add -A src/align_genotype/qc_calibration test/ pyproject.toml
git commit -m "refactor(qc_calibration): remove the CLI and its report renderers

cli, dryrun, emit, report, suggest and collect have no consumer now that the
stage is the entrypoint and settings/values/extract cover the data layer. The
qc_calibrate console script goes with them. spec, cache, tomlio and manifest
follow once discovery no longer imports them."
```

---

## Task 7: Rewrite `relative.py` without the global-config detour

The current module writes a throwaway TOML, installs it with `config.set_config_paths()`,
calls `check_multiqc.relative_flags` to count warns, then tries to restore the previous
paths. That was defensible in a local CLI. Inside a Hail Batch job, where CPG Flow has
already installed the run's config and every later `config_retrieve` depends on it,
mutating global config paths to obtain a warn count is a hazard.

The "no second MAD implementation anywhere" invariant is preserved: `robust_threshold`
is production's own function, `stats.churn` already calls it directly, and counting a
breach against a 4-dp-rounded threshold is exactly what `_relative_flags_for_metric`
does after it derives one.

**Files:**
- Modify: `src/align_genotype/qc_calibration/relative.py` (full rewrite)
- Replace: `test/test_qc_calibration_relative.py` (full rewrite)

- [ ] **Step 1: Write the failing test**

Replace the whole of `test/test_qc_calibration_relative.py` with:

```python
"""Unit tests for dataset-relative (MAD) tier evaluation."""

import numpy as np
import pytest

from align_genotype.qc_calibration import relative as relative_mod
from align_genotype.qc_calibration import settings as settings_mod
from align_genotype.qc_calibration import values as values_mod

METRIC = settings_mod.MetricSpec(key='dup_pct', direction='max', unit='%', relative=True)


def settings(min_samples=4, **bars) -> settings_mod.CalibrationSettings:
    return settings_mod.CalibrationSettings(
        seq_type='genome',
        metrics=(METRIC,),
        k=3.5,
        min_samples=min_samples,
        bars=settings_mod.Bars(**bars),
    )


def metric_values(values, section='samtools') -> values_mod.MetricValues:
    return values_mod.MetricValues(
        entries=tuple((section, f'CPG{i}', float(v)) for i, v in enumerate(values)),
        n_dropped=0,
    )


def test_derives_median_mad_and_threshold_per_dataset():
    by_dataset = {'ds-a': metric_values([10.0, 10.0, 12.0, 12.0, 40.0])}
    result = relative_mod.evaluate(by_dataset, METRIC, settings())
    (dataset,) = result.datasets
    assert dataset.dataset == 'ds-a'
    assert dataset.median == pytest.approx(12.0)
    assert dataset.mad_raw == pytest.approx(2.0)
    # median + k*MAD/0.6745 = 12 + 3.5*2/0.6745 = 22.3795...
    assert dataset.threshold == pytest.approx(22.3795, abs=1e-4)
    assert dataset.n_warn == 1  # only the 40.0
    assert dataset.skipped is None


def test_threshold_matches_production_rounded_to_four_dp():
    """The displayed threshold must be the number production compares against."""
    from align_genotype.scripts import check_multiqc  # noqa: PLC0415

    values = [10.0, 10.0, 12.0, 12.0, 40.0]
    expected = round(check_multiqc.robust_threshold(values, 'max', 3.5), 4)
    result = relative_mod.evaluate({'ds-a': metric_values(values)}, METRIC, settings())
    assert result.datasets[0].threshold == expected


def test_dataset_below_min_samples_is_skipped_with_a_size_reason():
    result = relative_mod.evaluate({'ds-a': metric_values([1.0, 2.0])}, METRIC, settings(min_samples=50))
    (dataset,) = result.datasets
    assert dataset.threshold is None
    assert dataset.n_warn == 0
    assert '2 values < min_samples 50' in dataset.skipped


def test_dataset_with_no_values_says_so_rather_than_too_small():
    """A collection problem and a size problem send an operator to different places."""
    result = relative_mod.evaluate({'ds-a': metric_values([])}, METRIC, settings())
    assert 'no values' in result.datasets[0].skipped


def test_zero_mad_dataset_is_skipped_as_degenerate():
    result = relative_mod.evaluate({'ds-a': metric_values([5.0] * 10)}, METRIC, settings())
    (dataset,) = result.datasets
    assert dataset.threshold is None
    assert 'zero MAD' in dataset.skipped


def test_warn_rate_is_per_value_and_duplication_is_flagged():
    by_dataset = {
        'ds-a': values_mod.MetricValues(
            entries=(
                ('picard_1', 'CPG0', 10.0),
                ('picard_4', 'CPG0', 10.0),
                ('picard_1', 'CPG1', 12.0),
                ('picard_4', 'CPG1', 12.0),
                ('picard_1', 'CPG2', 12.0),
                ('picard_4', 'CPG2', 12.0),
            ),
            n_dropped=0,
        ),
    }
    (dataset,) = relative_mod.evaluate(by_dataset, METRIC, settings()).datasets
    assert dataset.n_values == 6
    assert dataset.n_sequencing_groups == 3
    assert dataset.duplicated is True


def test_growth_churn_reports_both_orderings_and_uses_the_worse():
    rng = np.random.default_rng(1)
    # Leading 60% is tight, trailing 40% is high - the batch-ordering effect.
    ordered = [*list(rng.normal(10, 0.5, 30)), *list(rng.normal(20, 0.5, 20))]
    result = relative_mod.evaluate({'ds-a': metric_values(ordered)}, METRIC, settings(min_samples=10))
    (growth,) = result.growth
    assert growth.ordered is not None
    assert growth.shuffled is not None
    assert growth.flip_rate == max(growth.ordered.flip_rate, growth.shuffled.flip_rate)


def test_growth_is_not_simulated_when_the_before_slice_is_below_min_samples():
    """Production would have skipped a dataset that size, so churn against it is fiction."""
    result = relative_mod.evaluate({'ds-a': metric_values(range(1, 21))}, METRIC, settings(min_samples=20))
    assert result.growth == ()


def test_merge_churn_simulates_both_directions_of_every_pair():
    by_dataset = {
        'ds-a': metric_values([10.0, 10.5, 11.0, 11.5, 12.0, 40.0]),
        'ds-b': metric_values([30.0, 30.5, 31.0, 31.5, 32.0, 60.0]),
    }
    result = relative_mod.evaluate(by_dataset, METRIC, settings())
    assert {(a, b) for a, b, _ in result.merge} == {('ds-a', 'ds-b'), ('ds-b', 'ds-a')}


def test_merge_churn_never_self_pairs():
    by_dataset = {'ds-a': metric_values([10.0, 11.0, 12.0, 13.0, 40.0])}
    assert relative_mod.evaluate(by_dataset, METRIC, settings()).merge == ()


def test_verdict_recommends_when_every_bar_is_cleared():
    by_dataset = {'ds-a': metric_values([10.0, 10.5, 11.0, 11.5, 12.0, 11.2, 10.8, 11.4])}
    result = relative_mod.evaluate(by_dataset, METRIC, settings())
    assert result.verdict == 'RECOMMEND'
    assert result.verdict_reason == ''


def test_verdict_names_the_bar_that_was_missed():
    by_dataset = {'ds-a': metric_values([10.0, 10.5, 11.0, 11.5, 12.0, 40.0])}
    result = relative_mod.evaluate(by_dataset, METRIC, settings(max_warn_rate=0.0))
    assert result.verdict == 'REJECT'
    assert 'warn rate' in result.verdict_reason


def test_verdict_distinguishes_growth_from_merge_churn():
    """Conflating them would send an operator to the wrong table."""
    by_dataset = {
        'ds-a': metric_values([10.0, 10.5, 11.0, 11.5, 12.0, 40.0]),
        'ds-b': metric_values([30.0, 30.5, 31.0, 31.5, 32.0, 60.0]),
    }
    result = relative_mod.evaluate(by_dataset, METRIC, settings(max_merge_churn=0.0))
    assert 'merge churn' in result.verdict_reason
    assert 'growth churn' not in result.verdict_reason


def test_evaluate_refuses_a_metric_with_no_relative_tier():
    plain = settings_mod.MetricSpec(key='dup_pct', direction='max', unit='%', relative=False)
    with pytest.raises(ValueError, match='no configured relative tier'):
        relative_mod.evaluate({'ds-a': metric_values([1.0, 2.0, 3.0])}, plain, settings())


def test_max_warn_rate_ignores_skipped_datasets():
    by_dataset = {
        'ds-a': metric_values([10.0, 10.5, 11.0, 11.5, 12.0, 40.0]),
        'ds-tiny': metric_values([1.0]),
    }
    result = relative_mod.evaluate(by_dataset, METRIC, settings())
    skipped = next(d for d in result.datasets if d.dataset == 'ds-tiny')
    assert skipped.skipped is not None
    assert result.max_warn_rate == max(d.warn_rate for d in result.datasets if d.skipped is None)
```

- [ ] **Step 2: Run the tests to verify they fail**

Run: `uv run pytest test/test_qc_calibration_relative.py -v`

Expected: FAIL — `AttributeError` or `TypeError` on `relative_mod.evaluate`, whose current
signature is `(cache, metric, seq_type)`.

- [ ] **Step 3: Write the implementation**

Replace the whole of `src/align_genotype/qc_calibration/relative.py` with:

```python
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
    n_sequencing_groups: int
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
        return self.n_values > self.n_sequencing_groups


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
                f'peak dataset-growth churn {self.max_growth_churn:.2%} '
                f'exceeds {self.bars.max_growth_churn:.0%}',
            )
        if self.max_merge_churn > self.bars.max_merge_churn:
            reasons.append(
                f'peak cross-dataset merge churn {self.max_merge_churn:.2%} '
                f'exceeds {self.bars.max_merge_churn:.0%}',
            )
        return '; '.join(reasons)

    @property
    def verdict(self) -> str:
        return 'REJECT' if self.verdict_reason else 'RECOMMEND'


def _evaluate_dataset(dataset: str, metric_values: MetricValues, metric: MetricSpec, min_samples: int) -> DatasetMad:
    """One dataset's relative numbers, mirroring the skips production makes."""
    values = metric_values.array
    n_values = int(values.size)
    counts = {'n_values': n_values, 'n_sequencing_groups': metric_values.n_sequencing_groups}
    if n_values == 0:
        # Distinguished from "too small": nothing was extracted for this metric here,
        # which is a collection problem, not a size one. Reporting it as
        # "0 values < min_samples 50" would send a reader to the wrong place.
        return DatasetMad(
            dataset, **counts, median=float('nan'), mad_raw=float('nan'),
            threshold=None, n_warn=0,
            skipped=f'metric {metric.key!r} has no values in this dataset',
        )
    median = float(np.median(values))
    mad_raw = float(np.median(np.abs(values - median)))
    if n_values < min_samples:
        return DatasetMad(
            dataset, **counts, median=median, mad_raw=mad_raw, threshold=None, n_warn=0,
            skipped=f'{n_values} values < min_samples {min_samples}',
        )
    threshold = check_multiqc.robust_threshold(list(values), metric.direction, _k_for(metric))
    if threshold is None:
        return DatasetMad(
            dataset, **counts, median=median, mad_raw=mad_raw, threshold=None, n_warn=0,
            skipped='zero MAD (degenerate dataset); use the absolute gate here',
        )
    # Production rounds before comparing, so the displayed threshold always explains the
    # displayed count.
    threshold = round(threshold, 4)
    return DatasetMad(
        dataset, **counts, median=median, mad_raw=mad_raw, threshold=threshold,
        n_warn=int(stats.breach(values, threshold, metric.direction).sum()),
        skipped=None,
    )


# `k` is a single run-wide setting rather than per metric, but reading it through one
# function keeps the call sites honest if that ever changes.
_K: dict[str, float] = {}


def _k_for(metric: MetricSpec) -> float:
    return _K[metric.key]


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
    _K[metric.key] = settings.k
    datasets = tuple(
        _evaluate_dataset(name, metric_values, metric, settings.min_samples)
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
```

- [ ] **Step 4: Remove the module-level `_K` hack**

The `_K` dict above is process-global mutable state, which is exactly the class of thing
this task exists to remove. Replace it: delete the `_K` dict and `_k_for`, add a `k`
parameter to `_evaluate_dataset`, and pass `settings.k` at the call site.

```python
def _evaluate_dataset(
    dataset: str,
    metric_values: MetricValues,
    metric: MetricSpec,
    min_samples: int,
    k: float,
) -> DatasetMad:
```

with `threshold = check_multiqc.robust_threshold(list(values), metric.direction, k)`, and
in `evaluate`:

```python
    datasets = tuple(
        _evaluate_dataset(name, metric_values, metric, settings.min_samples, settings.k)
        for name, metric_values in by_dataset.items()
    )
```

Then delete the `_K[metric.key] = settings.k` line.

- [ ] **Step 5: Run the tests to verify they pass**

Run: `uv run pytest test/test_qc_calibration_relative.py -v`

Expected: PASS, 15 tests.

- [ ] **Step 6: Verify the global-config detour is gone**

Run: `grep -n 'set_config_paths\|tomlio\|relative_flags' src/align_genotype/qc_calibration/relative.py`

Expected: no output.

- [ ] **Step 7: Commit**

```bash
git add src/align_genotype/qc_calibration/relative.py test/test_qc_calibration_relative.py
git commit -m "refactor(qc_calibration): evaluate MAD tiers without mutating global config

Warn counts came from relative_flags driven through a throwaway TOML installed
with set_config_paths. Inside a Hail Batch job, where cpg-flow has already
installed the run config, that is a hazard for a warn count. robust_threshold
is production's own function, so counting a breach against the same 4-dp
threshold keeps calibration and enforcement identical without the detour."
```

---

## Task 8: `thresholds.py` — candidate fixed thresholds from percentile tails

Replaces `suggest.py`. The old module wrote seeds back into a spec file behind a
`reviewed = false` interlock; there is no spec file now, so it simply returns candidates.

**Files:**
- Create: `src/align_genotype/qc_calibration/thresholds.py`
- Test: `test/test_qc_calibration_thresholds.py`

- [ ] **Step 1: Write the failing test**

Create `test/test_qc_calibration_thresholds.py`:

```python
"""Unit tests for candidate fixed thresholds."""

import pytest

from align_genotype.qc_calibration import settings as settings_mod
from align_genotype.qc_calibration import thresholds as thresholds_mod
from align_genotype.qc_calibration import values as values_mod


def metric_values(values) -> values_mod.MetricValues:
    return values_mod.MetricValues(
        entries=tuple(('s', f'CPG{i}', float(v)) for i, v in enumerate(values)),
        n_dropped=0,
    )


def test_min_metric_seeds_from_the_lowest_per_dataset_tails():
    """Worst-case per dataset, so a candidate does not already flag the marginal one."""
    metric = settings_mod.MetricSpec(key='COV', direction='min', unit='x')
    by_dataset = {'ds-a': metric_values(range(30, 131)), 'ds-b': metric_values(range(10, 111))}
    candidate = thresholds_mod.candidate(by_dataset, metric)
    # ds-b is the marginal dataset: its p1 is ~11 and its p5 is ~15.
    assert candidate.fail == 11
    assert candidate.warn == 15


def test_max_metric_seeds_from_the_highest_per_dataset_tails():
    metric = settings_mod.MetricSpec(key='DUP', direction='max', unit='%')
    by_dataset = {'ds-a': metric_values(range(1, 102)), 'ds-b': metric_values(range(20, 121))}
    candidate = thresholds_mod.candidate(by_dataset, metric)
    assert candidate.fail == 119
    assert candidate.warn == 115


def test_x_and_percent_units_round_to_whole_numbers():
    metric = settings_mod.MetricSpec(key='COV', direction='min', unit='x')
    candidate = thresholds_mod.candidate({'ds': metric_values([10.4, 20.0, 30.0])}, metric)
    assert isinstance(candidate.fail, int)


def test_frac_unit_keeps_two_decimal_places():
    """A fraction rounded to an integer would collapse to 0 or 1."""
    metric = settings_mod.MetricSpec(key='PCT_20X', direction='min', unit='frac')
    candidate = thresholds_mod.candidate({'ds': metric_values([0.9012, 0.95, 0.97])}, metric)
    assert candidate.fail == pytest.approx(0.9, abs=0.01)
    assert isinstance(candidate.fail, float)


def test_a_relative_metric_gets_no_absolute_warn():
    """The relative tier is the warn tier; both would double-flag the same values."""
    metric = settings_mod.MetricSpec(key='DUP', direction='max', unit='%', relative=True)
    candidate = thresholds_mod.candidate({'ds': metric_values(range(1, 102))}, metric)
    assert candidate.fail is not None
    assert candidate.warn is None
    assert 'relative' in candidate.basis


def test_a_metric_with_no_values_anywhere_yields_no_candidate():
    metric = settings_mod.MetricSpec(key='ABSENT', direction='min', unit='x')
    assert thresholds_mod.candidate({'ds-a': metric_values([]), 'ds-b': metric_values([])}, metric) is None


def test_datasets_without_values_are_skipped_not_counted():
    metric = settings_mod.MetricSpec(key='COV', direction='min', unit='x')
    by_dataset = {'ds-a': metric_values(range(30, 131)), 'ds-empty': metric_values([])}
    candidate = thresholds_mod.candidate(by_dataset, metric)
    assert '1 dataset' in candidate.basis


def test_basis_names_the_percentiles_used():
    metric = settings_mod.MetricSpec(key='COV', direction='min', unit='x')
    candidate = thresholds_mod.candidate({'ds': metric_values(range(30, 131))}, metric)
    assert 'p1' in candidate.basis
    assert 'p5' in candidate.basis
```

- [ ] **Step 2: Run the tests to verify they fail**

Run: `uv run pytest test/test_qc_calibration_thresholds.py -v`

Expected: FAIL — `ModuleNotFoundError: No module named 'align_genotype.qc_calibration.thresholds'`.

- [ ] **Step 3: Write the implementation**

Create `src/align_genotype/qc_calibration/thresholds.py`:

```python
"""Candidate fixed thresholds, seeded from the per-dataset distribution tails.

These are starting points, not answers, and the report says so. What a percentile cannot
supply is the judgement a threshold actually turns on: preserving the lab's intent for a
hard gate unless the data clearly contradicts it, refusing to hard-fail on metrics that
track ancestry, biology or chemistry rather than sample quality, and overriding the lab
only with a measured count of good sequencing groups their line would have discarded.

Which tail seeds which tier: for `min` metrics (higher is better) the bad values sit in
the low tail, so `fail` comes from the lowest per-dataset p1 and `warn` from the lowest
p5. For `max` metrics it is the high tail - highest p99 for `fail`, highest p95 for
`warn`. Taking the worst per-dataset value rather than pooling means a candidate does not
already flag a large slice of the most marginal dataset in the set.
"""

from dataclasses import dataclass

import numpy as np

from align_genotype.qc_calibration.settings import MetricSpec
from align_genotype.qc_calibration.values import MetricValues

_SEED_PERCENTILES: dict[str, dict[str, int]] = {
    'min': {'fail': 1, 'warn': 5},
    'max': {'fail': 99, 'warn': 95},
}

# 'x' (coverage) and '%' round to whole numbers for readability; 'frac' keeps 2 dp, since
# a fraction rounded to an integer would collapse to 0 or 1.
_INTEGER_UNITS = {'x', '%'}
_FRACTION_DP = 2


@dataclass(frozen=True)
class Candidate:
    """One metric's candidate thresholds and the evidence behind them."""

    metric: str
    fail: float | None
    warn: float | None
    basis: str


def _round_for_unit(value: float, unit: str) -> float | int:
    """Round a raw percentile to a sensible reading for `unit`, as a Python builtin.

    `round(np.float64(...), 2)` returns another `np.float64`, which serialises to JSON
    only by accident, so the fraction branch casts first. The integer branch calls
    `round` with no `ndigits` on a plain float, which returns a genuine `int`.
    """
    if unit in _INTEGER_UNITS:
        return round(float(value))
    return round(float(value), _FRACTION_DP)


def _tail(by_dataset: dict[str, MetricValues], metric: MetricSpec, percentile: int) -> float | None:
    """The worst per-dataset value of `percentile`, or None if no dataset has data.

    'Worst' is direction-dependent: the lowest per-dataset percentile for a `min` metric,
    the highest for a `max` one.
    """
    readings = [
        float(np.percentile(values, percentile))
        for metric_values in by_dataset.values()
        if (values := metric_values.array).size
    ]
    if not readings:
        return None
    return min(readings) if metric.direction == 'min' else max(readings)


def candidate(by_dataset: dict[str, MetricValues], metric: MetricSpec) -> Candidate | None:
    """Candidate thresholds for one metric, or None when no dataset has any values.

    A metric with a relative tier gets a `fail` candidate only: its warn tier is derived
    per run from the dataset's own spread, and an absolute warn alongside it would
    double-flag the same values and disagree about which threshold was breached.
    """
    pcts = _SEED_PERCENTILES[metric.direction]
    fail_raw = _tail(by_dataset, metric, pcts['fail'])
    if fail_raw is None:
        return None

    n_with_data = sum(1 for v in by_dataset.values() if v.array.size)
    plural = '' if n_with_data == 1 else 's'
    basis = f'worst per-dataset p{pcts["fail"]}={fail_raw:.4g} across {n_with_data} dataset{plural}'

    warn: float | None = None
    if metric.relative:
        basis += '; warn tier is dataset-relative, so no absolute warn is proposed'
    else:
        warn_raw = _tail(by_dataset, metric, pcts['warn'])
        # Unreachable: `_tail` walks the same per-dataset arrays at a different
        # percentile, so a dataset with fail-percentile data always has warn-percentile
        # data. Raised rather than asserted, since ruff bans `assert` outside test/**.
        if warn_raw is None:
            raise AssertionError('a dataset with fail-percentile data must also have warn-percentile data')
        warn = _round_for_unit(warn_raw, metric.unit)
        basis += f', warn p{pcts["warn"]}={warn_raw:.4g}'

    return Candidate(
        metric=metric.key,
        fail=_round_for_unit(fail_raw, metric.unit),
        warn=warn,
        basis=basis,
    )
```

- [ ] **Step 4: Run the tests to verify they pass**

Run: `uv run pytest test/test_qc_calibration_thresholds.py -v`

Expected: PASS, 8 tests.

- [ ] **Step 5: Commit**

```bash
git add src/align_genotype/qc_calibration/thresholds.py test/test_qc_calibration_thresholds.py
git commit -m "feat(qc_calibration): derive candidate thresholds from distribution tails

Replaces suggest.py. There is no spec file to seed and no reviewed interlock to
enforce, so this returns candidates for the report rather than mutating a spec."
```

---

## Task 9: `snippet.py` — render the `qc_thresholds` block

Replaces `emit.py`. The `reviewed` interlock and the dataset-label guard both go: nothing
here is pasted into a public file by the tool, and the report lives in a private web
bucket where dataset names are the point.

**Files:**
- Create: `src/align_genotype/qc_calibration/snippet.py`
- Test: `test/test_qc_calibration_snippet.py`

- [ ] **Step 1: Write the failing test**

Create `test/test_qc_calibration_snippet.py`:

```python
"""Unit tests for the copy-pasteable qc_thresholds block."""

import sys

from align_genotype.qc_calibration import settings as settings_mod
from align_genotype.qc_calibration import snippet as snippet_mod
from align_genotype.qc_calibration.thresholds import Candidate

if sys.version_info >= (3, 11):
    import tomllib
else:
    import tomli as tomllib

SETTINGS = settings_mod.CalibrationSettings(
    seq_type='genome',
    metrics=(
        settings_mod.MetricSpec(key='MEDIAN_COVERAGE', direction='min', unit='x'),
        settings_mod.MetricSpec(key='FREEMIX', direction='max', unit='frac'),
        settings_mod.MetricSpec(key='reads_duplicated_percent', direction='max', unit='%', relative=True),
    ),
)

CANDIDATES = {
    'MEDIAN_COVERAGE': Candidate('MEDIAN_COVERAGE', fail=17, warn=26, basis='...'),
    'FREEMIX': Candidate('FREEMIX', fail=0.03, warn=0.015, basis='...'),
    'reads_duplicated_percent': Candidate('reads_duplicated_percent', fail=38, warn=None, basis='...'),
}


def test_block_parses_as_toml():
    parsed = tomllib.loads(snippet_mod.render(SETTINGS, CANDIDATES))
    assert parsed['qc_thresholds']['genome']['fail']['min'] == {'MEDIAN_COVERAGE': 17}
    assert parsed['qc_thresholds']['genome']['fail']['max'] == {
        'FREEMIX': 0.03,
        'reads_duplicated_percent': 38,
    }
    assert parsed['qc_thresholds']['genome']['warn']['min'] == {'MEDIAN_COVERAGE': 26}
    assert parsed['qc_thresholds']['genome']['warn']['max'] == {'FREEMIX': 0.015}


def test_relative_metric_has_no_absolute_warn_entry():
    parsed = tomllib.loads(snippet_mod.render(SETTINGS, CANDIDATES))
    assert 'reads_duplicated_percent' not in parsed['qc_thresholds']['genome']['warn']['max']


def test_relative_table_carries_direction_k_and_min_samples():
    parsed = tomllib.loads(snippet_mod.render(SETTINGS, CANDIDATES))
    relative = parsed['qc_thresholds']['genome']['relative']['reads_duplicated_percent']
    assert relative == {'direction': 'max', 'k': 3.5, 'min_samples': 50}


def test_block_round_trips_through_the_production_loader(monkeypatch):
    """Proof the emitted shape is what load_thresholds actually reads."""
    from align_genotype.scripts import check_multiqc  # noqa: PLC0415

    parsed = tomllib.loads(snippet_mod.render(SETTINGS, CANDIDATES))
    monkeypatch.setattr(
        check_multiqc.config,
        'config_retrieve',
        lambda keys, default=None: _dig(parsed, keys, default),
    )
    loaded = check_multiqc.load_thresholds('genome')
    assert loaded['min']['MEDIAN_COVERAGE'] == {'fail': 17, 'warn': 26}
    assert loaded['max']['reads_duplicated_percent'] == {'fail': 38}


def _dig(tree, keys, default):
    node = tree
    for key in keys:
        if not isinstance(node, dict) or key not in node:
            return default
        node = node[key]
    return node


def test_section_order_matches_the_committed_config():
    """So the block diffs against config_template.toml as a change, not a rewrite."""
    body = snippet_mod.render(SETTINGS, CANDIDATES)
    order = [line for line in body.splitlines() if line.startswith('[qc_thresholds')]
    assert order == [
        '[qc_thresholds.genome.fail.min]',
        '[qc_thresholds.genome.fail.max]',
        '[qc_thresholds.genome.warn.min]',
        '[qc_thresholds.genome.warn.max]',
        '[qc_thresholds.genome.relative.reads_duplicated_percent]',
    ]


def test_a_metric_with_no_candidate_is_omitted():
    body = snippet_mod.render(SETTINGS, {'MEDIAN_COVERAGE': CANDIDATES['MEDIAN_COVERAGE']})
    assert 'FREEMIX' not in body


def test_empty_candidates_render_an_empty_block():
    assert snippet_mod.render(SETTINGS, {}).strip() == ''


def test_metric_keys_are_quoted_like_the_committed_config():
    assert '"MEDIAN_COVERAGE" = 17' in snippet_mod.render(SETTINGS, CANDIDATES)
```

- [ ] **Step 2: Run the tests to verify they fail**

Run: `uv run pytest test/test_qc_calibration_snippet.py -v`

Expected: FAIL — `ModuleNotFoundError: No module named 'align_genotype.qc_calibration.snippet'`.

- [ ] **Step 3: Write the implementation**

Create `src/align_genotype/qc_calibration/snippet.py`:

```python
"""Render a `[qc_thresholds.<seq_type>...]` block for pasting into config_template.toml.

This produces a string for a human to diff and paste; nothing writes the config. That
file carries comments, ordering and judgement no generator reproduces, and the paste-
after-diff step is the cheapest guard against a calibration run quietly rewriting a
production gate. Section order matches the committed file so the diff reads as a change
rather than a rewrite.
"""

from align_genotype.qc_calibration.settings import CalibrationSettings
from align_genotype.qc_calibration.thresholds import Candidate

# (severity, direction) in the order config_template.toml writes them.
_SECTIONS: tuple[tuple[str, str], ...] = (('fail', 'min'), ('fail', 'max'), ('warn', 'min'), ('warn', 'max'))


def _fmt(value: float) -> str:
    """Render a threshold as a TOML scalar.

    `hasattr(value, 'item')` unwraps a numpy scalar first: `np.float64` reprs as
    `np.float64(0.75)`, which is not valid TOML, and `np.int64` would fall through to
    `str()` and be written as a quoted string - valid TOML of the wrong type, which is
    the worse failure. Candidates are cast in `thresholds._round_for_unit`, but this is
    the single sink for everything written, so it defends here too.
    """
    if hasattr(value, 'item'):
        value = value.item()
    return repr(value)


def render(settings: CalibrationSettings, candidates: dict[str, Candidate]) -> str:
    """The config block for every metric that produced a candidate."""
    lines: list[str] = []
    for severity, direction in _SECTIONS:
        entries = [
            (metric.key, threshold)
            for metric in settings.metrics
            if metric.direction == direction
            and (found := candidates.get(metric.key)) is not None
            and (threshold := getattr(found, severity)) is not None
        ]
        if not entries:
            continue
        lines += ['', f'[qc_thresholds.{settings.seq_type}.{severity}.{direction}]']
        lines += [f'"{key}" = {_fmt(threshold)}' for key, threshold in entries]

    for metric in settings.relative_metrics:
        if metric.key not in candidates:
            continue
        lines += [
            '',
            f'[qc_thresholds.{settings.seq_type}.relative.{metric.key}]',
            f'direction = "{metric.direction}"',
            f'k = {_fmt(settings.k)}',
            f'min_samples = {_fmt(settings.min_samples)}',
        ]

    return '\n'.join(lines).lstrip('\n') + ('\n' if lines else '')
```

- [ ] **Step 4: Run the tests to verify they pass**

Run: `uv run pytest test/test_qc_calibration_snippet.py -v`

Expected: PASS, 8 tests.

- [ ] **Step 5: Commit**

```bash
git add src/align_genotype/qc_calibration/snippet.py test/test_qc_calibration_snippet.py
git commit -m "feat(qc_calibration): render the qc_thresholds block for pasting

Replaces emit.py. The reviewed interlock has nothing left to gate and the
dataset-label guard protected a public file this no longer writes to, so both
go; what remains is the section ordering that keeps the diff readable."
```

---

## Task 10: `summary.py` — assemble the one result structure

Replaces `report.py`. The old module rendered fixed-width ASCII tables, and roughly 250
of its 587 lines were column-packing machinery (`_label_groups`, `MAX_TABLE_WIDTH`
chunking) that exists only because ten dataset labels overflow a terminal. HTML has no
such problem.

`build` returns a plain dict rather than a dataclass tree, because that dict *is* both
the `calibration.json` output and the template context. A parallel dataclass hierarchy
would be ~80 lines mirroring a schema the spec already fixes.

**Files:**
- Create: `src/align_genotype/qc_calibration/summary.py`
- Test: `test/test_qc_calibration_summary.py`

- [ ] **Step 1: Write the failing test**

Create `test/test_qc_calibration_summary.py`:

```python
"""Unit tests for the assembled calibration summary."""

import pytest

from align_genotype.qc_calibration import settings as settings_mod
from align_genotype.qc_calibration import summary as summary_mod
from align_genotype.qc_calibration import values as values_mod

SETTINGS = settings_mod.CalibrationSettings(
    seq_type='genome',
    metrics=(
        settings_mod.MetricSpec(key='MEDIAN_COVERAGE', direction='min', unit='x'),
        settings_mod.MetricSpec(key='dup_pct', direction='max', unit='%', relative=True),
        settings_mod.MetricSpec(key='ABSENT', direction='min', unit='x'),
    ),
    k=3.5,
    min_samples=4,
)

CURRENT = {'MEDIAN_COVERAGE': {'fail': 15, 'warn': 25}}


def dataset_values(name, coverage, dup, n_dropped=0) -> values_mod.DatasetValues:
    return values_mod.DatasetValues(
        dataset=name,
        seq_type='genome',
        analysis_id=1,
        timestamp='2026-06-01T00:00:00',
        uri=f'gs://{name}/multiqc_data.json',
        multiqc_version='1.33',
        generated='2026-08-12T00:00:00',
        n_sequencing_groups=len(coverage),
        section_sizes={'picard_1': len(coverage)},
        metrics={
            'MEDIAN_COVERAGE': values_mod.MetricValues(
                entries=tuple(('picard_1', f'{name}-CPG{i}', float(v)) for i, v in enumerate(coverage)),
                n_dropped=n_dropped,
            ),
            'dup_pct': values_mod.MetricValues(
                entries=tuple(('samtools', f'{name}-CPG{i}', float(v)) for i, v in enumerate(dup)),
                n_dropped=0,
            ),
            'ABSENT': values_mod.MetricValues(entries=(), n_dropped=0),
        },
    )


@pytest.fixture
def built():
    return summary_mod.build(
        [
            dataset_values('ds-a', [30, 32, 34, 36, 38, 10], [10, 10.5, 11, 11.5, 12, 40], n_dropped=2),
            dataset_values('ds-b', [40, 42, 44, 46, 48, 50], [7, 7.5, 8, 8.5, 9, 9.5]),
        ],
        SETTINGS,
        current=CURRENT,
        skipped_datasets=[{'dataset': 'ds-c', 'reason': 'no completed CramMultiQC qc analysis for genome'}],
        generated='2026-08-12T00:00:00',
        ar_guid='test-ar-guid',
    )


def test_records_run_level_provenance(built):
    assert built['sequencing_type'] == 'genome'
    assert built['ar_guid'] == 'test-ar-guid'
    assert built['settings']['k'] == pytest.approx(3.5)
    assert built['settings']['min_samples'] == 4
    assert [d['dataset'] for d in built['datasets']] == ['ds-a', 'ds-b']
    assert built['datasets'][0]['uri'] == 'gs://ds-a/multiqc_data.json'
    assert built['skipped_datasets'][0]['dataset'] == 'ds-c'


def test_counts_sequencing_groups_and_values_separately(built):
    metric = built['metrics']['MEDIAN_COVERAGE']
    assert metric['n_values'] == 12
    assert metric['n_sequencing_groups'] == 12
    assert metric['n_datasets'] == 2
    assert metric['n_dropped'] == 2


def test_carries_the_shipped_thresholds_for_comparison(built):
    assert built['metrics']['MEDIAN_COVERAGE']['current'] == {'fail': 15, 'warn': 25}
    assert built['metrics']['dup_pct']['current'] == {}


def test_proposes_a_candidate_per_metric_with_data(built):
    assert built['metrics']['MEDIAN_COVERAGE']['candidate']['fail'] is not None
    assert built['metrics']['dup_pct']['candidate']['warn'] is None
    assert built['metrics']['ABSENT']['candidate'] is None


def test_scores_flag_rates_for_both_current_and_candidate(built):
    rates = built['metrics']['MEDIAN_COVERAGE']['flag_rates']
    # ds-a has one value of 10, below the shipped fail of 15.
    assert rates['current']['ds-a']['fail'] == pytest.approx(1 / 6)
    assert rates['current']['ds-b']['fail'] == pytest.approx(0.0)
    assert set(rates['candidate']) == {'ds-a', 'ds-b'}


def test_reports_percentiles_per_dataset(built):
    percentiles = built['metrics']['MEDIAN_COVERAGE']['percentiles']
    assert set(percentiles) == {'ds-a', 'ds-b'}
    assert percentiles['ds-b']['p50'] == pytest.approx(45.0)


def test_records_where_each_metric_was_found(built):
    assert built['metrics']['MEDIAN_COVERAGE']['present_in'] == ['picard_1']
    assert built['metrics']['MEDIAN_COVERAGE']['missing_from'] == []
    assert built['metrics']['ABSENT']['missing_from'] == ['ds-a', 'ds-b']


def test_a_metric_missing_everywhere_produces_a_loud_warning(built):
    assert any('every dataset' in w for w in built['warnings'])
    assert any('ABSENT' in w for w in built['warnings'])


def test_a_metric_present_everywhere_produces_no_warning(built):
    assert not any('MEDIAN_COVERAGE' in w for w in built['warnings'])


def test_evaluates_only_the_relative_metrics(built):
    assert set(built['relative']) == {'dup_pct'}
    evaluation = built['relative']['dup_pct']
    assert evaluation['verdict'] in ('RECOMMEND', 'REJECT')
    assert {d['dataset'] for d in evaluation['datasets']} == {'ds-a', 'ds-b'}


def test_relative_evaluation_reports_all_three_bars(built):
    evaluation = built['relative']['dup_pct']
    assert 'max_warn_rate' in evaluation
    assert 'max_growth_churn' in evaluation
    assert 'max_merge_churn' in evaluation
    assert evaluation['bars']['max_merge_churn'] == pytest.approx(0.05)


def test_includes_a_pasteable_config_snippet(built):
    assert '[qc_thresholds.genome.fail.min]' in built['config_snippet']


def test_is_json_serialisable(built):
    import json  # noqa: PLC0415

    json.dumps(built, allow_nan=False)


def test_no_datasets_at_all_still_builds(built):  # noqa: ARG001
    empty = summary_mod.build(
        [], SETTINGS, current={}, skipped_datasets=[], generated='2026-08-12T00:00:00', ar_guid='x',
    )
    assert empty['datasets'] == []
    assert empty['metrics']['MEDIAN_COVERAGE']['candidate'] is None
    assert any('every dataset' in w for w in empty['warnings'])
```

- [ ] **Step 2: Run the tests to verify they fail**

Run: `uv run pytest test/test_qc_calibration_summary.py -v`

Expected: FAIL — `ModuleNotFoundError: No module named 'align_genotype.qc_calibration.summary'`.

- [ ] **Step 3: Write the implementation**

Create `src/align_genotype/qc_calibration/summary.py`:

```python
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
        'ordering_sensitive': list(evaluation.ordering_sensitive),
        'datasets': [
            {
                'dataset': d.dataset,
                'n_values': d.n_values,
                'n_sequencing_groups': d.n_sequencing_groups,
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
            'n_sequencing_groups': len({sg for mv in per_dataset.values() for _, sg, _ in mv.entries}),
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
```

- [ ] **Step 4: Run the tests to verify they pass**

Run: `uv run pytest test/test_qc_calibration_summary.py -v`

Expected: PASS, 14 tests.

- [ ] **Step 5: Commit**

```bash
git add src/align_genotype/qc_calibration/summary.py test/test_qc_calibration_summary.py
git commit -m "feat(qc_calibration): assemble the calibration summary

Replaces report.py. One dict serves both calibration.json and the template, so
the HTML cannot show a number the JSON does not carry. The ~250 lines of ASCII
column-packing that existed because ten dataset labels overflow a terminal go
with it."
```

---

## Task 11: The HTML dashboard

**Files:**
- Create: `src/align_genotype/qc_calibration/render.py`
- Create: `src/align_genotype/templates/qc_calibration_report.html.jinja`
- Test: `test/test_qc_calibration_render.py`

- [ ] **Step 1: Write the failing test**

Create `test/test_qc_calibration_render.py`:

```python
"""Smoke tests for the calibration HTML report."""

import pytest

from align_genotype.qc_calibration import render as render_mod
from test.test_qc_calibration_summary import SETTINGS, dataset_values  # noqa: PT013
from align_genotype.qc_calibration import summary as summary_mod


@pytest.fixture
def html() -> str:
    built = summary_mod.build(
        [
            dataset_values('ds-a', [30, 32, 34, 36, 38, 10], [10, 10.5, 11, 11.5, 12, 40]),
            dataset_values('ds-b', [40, 42, 44, 46, 48, 50], [7, 7.5, 8, 8.5, 9, 9.5]),
        ],
        SETTINGS,
        current={'MEDIAN_COVERAGE': {'fail': 15, 'warn': 25}},
        skipped_datasets=[{'dataset': 'ds-c', 'reason': 'no completed CramMultiQC qc analysis for genome'}],
        generated='2026-08-12T00:00:00',
        ar_guid='test-ar-guid',
    )
    return render_mod.render(built)


def test_renders_a_complete_html_document(html):
    assert html.startswith('<!DOCTYPE html>')
    assert html.rstrip().endswith('</html>')


def test_carries_every_expected_section(html):
    for heading in (
        'Recommended fixed thresholds',
        'Dataset-relative (MAD) tiers',
        'Metric presence',
        'Percentile distributions',
        'Flag rates',
        'Churn detail',
        'Provenance',
        'Config block',
    ):
        assert heading in html, f'missing section: {heading}'


def test_headline_reports_the_run_shape(html):
    assert 'genome' in html
    assert 'test-ar-guid' in html


def test_shows_current_and_candidate_side_by_side(html):
    assert 'MEDIAN_COVERAGE' in html
    assert 'Current' in html
    assert 'Candidate' in html


def test_names_the_skipped_dataset_and_why(html):
    assert 'ds-c' in html
    assert 'no completed CramMultiQC qc analysis' in html


def test_banners_a_metric_absent_from_every_dataset(html):
    assert 'ABSENT' in html
    assert 'checks nothing' in html


def test_includes_the_pasteable_config_block(html):
    assert '[qc_thresholds.genome.fail.min]' in html


def test_states_that_candidates_are_not_decisions(html):
    assert 'not decisions' in html


def test_explains_the_merge_churn_tension(html):
    assert 'merge' in html.lower()
    assert 'shifts between datasets' in html


def test_escapes_html_in_dataset_names():
    """Dataset names reach the page unmodified, so autoescape must be on."""
    built = summary_mod.build(
        [dataset_values('<script>x</script>', [30, 32, 34, 36], [1, 2, 3, 4])],
        SETTINGS,
        current={},
        skipped_datasets=[],
        generated='2026-08-12T00:00:00',
        ar_guid='x',
    )
    assert '<script>x</script>' not in render_mod.render(built)
```

- [ ] **Step 2: Run the tests to verify they fail**

Run: `uv run pytest test/test_qc_calibration_render.py -v`

Expected: FAIL — `ModuleNotFoundError: No module named 'align_genotype.qc_calibration.render'`.

- [ ] **Step 3: Write the render module**

Create `src/align_genotype/qc_calibration/render.py`:

```python
"""Render the calibration summary to HTML.

The template directory and `autoescape=True` match `scripts/sg_qc_report.py`: dataset
names and skip reasons reach the page unmodified, so escaping is not optional.
"""

from pathlib import Path
from typing import Any

import jinja2

JINJA_TEMPLATE_DIR = Path(__file__).absolute().parent.parent / 'templates'
TEMPLATE_NAME = 'qc_calibration_report.html.jinja'


def _pct(value: float | None) -> str:
    """A rate as a percentage, or an em dash when there was nothing to score."""
    return '—' if value is None else f'{value:.1%}'


def _num(value: float | None, unit: str = 'frac') -> str:
    """A metric value at a precision that reads sensibly for its unit."""
    if value is None:
        return '—'
    if unit == 'frac':
        return f'{value:.3f}'
    if unit in ('x', '%'):
        return f'{value:.1f}'
    return f'{value:.2f}'


def _stat(value: float | None) -> str:
    """A median, MAD or threshold at 4 dp, so the threshold explains the warn count."""
    return '—' if value is None else f'{value:.4f}'


def render(built: dict[str, Any]) -> str:
    """Render the dashboard for one assembled summary."""
    env = jinja2.Environment(loader=jinja2.FileSystemLoader(JINJA_TEMPLATE_DIR), autoescape=True)
    env.filters['pct'] = _pct
    env.filters['num'] = _num
    env.filters['stat'] = _stat
    return env.get_template(TEMPLATE_NAME).render(
        r=built,
        dataset_names=[d['dataset'] for d in built['datasets']],
        percentile_labels=[f'p{p}' for p in (1, 5, 10, 25, 50, 75, 90, 95, 99)],
    )
```

- [ ] **Step 4: Write the template**

Create `src/align_genotype/templates/qc_calibration_report.html.jinja`:

```jinja
<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8">
<meta name="viewport" content="width=device-width, initial-scale=1.0">
<title>QC threshold calibration — {{ r.sequencing_type }}</title>
<style>
  * { margin: 0; padding: 0; box-sizing: border-box; }
  body { font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', Roboto, sans-serif;
         background: #f5f6fa; color: #2d3436; padding: 24px; }
  .container { max-width: 1400px; margin: 0 auto; }
  h1 { font-size: 1.6rem; margin-bottom: 4px; }
  h2 { font-size: 1.15rem; margin: 32px 0 12px; }
  h3 { font-size: 0.95rem; margin: 20px 0 8px; color: #2d3436; }
  .subtitle { color: #636e72; font-size: 0.9rem; margin-bottom: 16px; }
  .blurb { color: #636e72; font-size: 0.85rem; max-width: 900px; line-height: 1.5; margin-bottom: 16px; }
  .cards { display: flex; gap: 16px; flex-wrap: wrap; margin-bottom: 24px; }
  .card { background: #fff; border-radius: 8px; padding: 16px 24px;
          box-shadow: 0 1px 3px rgba(0,0,0,0.08); min-width: 140px; }
  .card .label { font-size: 0.75rem; color: #636e72; text-transform: uppercase; letter-spacing: 0.5px; }
  .card .value { font-size: 1.7rem; font-weight: 700; margin-top: 4px; }
  .banner { border-radius: 8px; padding: 14px 18px; margin-bottom: 16px; line-height: 1.5;
            background: #fdecea; color: #a4271b; border: 1px solid #f5c2bd; font-size: 0.88rem; }
  .banner ul { margin: 8px 0 0 20px; }
  table { width: 100%; border-collapse: collapse; background: #fff; border-radius: 8px;
          overflow: hidden; box-shadow: 0 1px 3px rgba(0,0,0,0.08); font-size: 0.85rem;
          margin-bottom: 16px; }
  th, td { padding: 8px 12px; text-align: left; border-bottom: 1px solid #eef1f4; white-space: nowrap; }
  th { background: #f1f3f5; font-weight: 600; font-size: 0.78rem; text-transform: uppercase;
       letter-spacing: 0.4px; color: #636e72; }
  td.metric { font-family: ui-monospace, SFMono-Regular, Menlo, monospace; font-weight: 600; }
  td.num { text-align: right; font-variant-numeric: tabular-nums; }
  .verdict { font-weight: 700; padding: 2px 8px; border-radius: 10px; font-size: 0.75rem; }
  .verdict.RECOMMEND { background: #e8f8f2; color: #0b8a63; }
  .verdict.REJECT { background: #fdecea; color: #a4271b; }
  .over { color: #a4271b; font-weight: 600; }
  .muted { color: #95a5a6; }
  .skip { color: #b9770e; font-size: 0.8rem; white-space: normal; }
  details { background: #fff; border-radius: 8px; padding: 14px 18px; margin-bottom: 10px;
            box-shadow: 0 1px 3px rgba(0,0,0,0.08); }
  details summary { cursor: pointer; font-weight: 600; font-size: 0.95rem; }
  details[open] summary { margin-bottom: 12px; }
  details .scroll { overflow-x: auto; }
  pre { background: #2d3436; color: #dfe6e9; padding: 16px; border-radius: 8px;
        overflow-x: auto; font-size: 0.8rem; line-height: 1.5; }
</style>
</head>
<body>
<div class="container">

  <h1>QC threshold calibration — {{ r.sequencing_type }}</h1>
  <p class="subtitle">
    Generated {{ r.generated }} · ar-guid <code>{{ r.ar_guid }}</code> ·
    k = {{ r.settings.k }} · min_samples = {{ r.settings.min_samples }}
  </p>

  <div class="cards">
    <div class="card"><div class="label">Datasets</div><div class="value">{{ r.datasets | length }}</div></div>
    <div class="card"><div class="label">Sequencing groups</div>
      <div class="value">{{ r.datasets | sum(attribute='n_sequencing_groups') }}</div></div>
    <div class="card"><div class="label">Metrics</div><div class="value">{{ r.metrics | length }}</div></div>
    <div class="card"><div class="label">Relative tiers</div><div class="value">{{ r.relative | length }}</div></div>
  </div>

  {% if r.warnings %}
  <div class="banner">
    <strong>Metric coverage problems</strong>
    <ul>{% for warning in r.warnings %}<li>{{ warning }}</li>{% endfor %}</ul>
  </div>
  {% endif %}

  <h2>Recommended fixed thresholds</h2>
  <p class="blurb">
    Candidates are the worst per-dataset percentile tails — descriptions of where the data sits, <strong>not
    decisions</strong>. A percentile cannot know whether a line contradicts the lab's intent, or whether a metric
    tracks ancestry and chemistry rather than sample quality. Read the flag rates below against what is shipped
    today, and decide as a team.
  </p>
  <div class="scroll">
  <table>
    <tr>
      <th>Metric</th><th>Dir</th><th>Unit</th>
      <th>Current fail</th><th>Current warn</th>
      <th>Candidate fail</th><th>Candidate warn</th>
      <th>Values</th><th>Groups</th><th>Datasets</th><th>Dropped</th>
    </tr>
    {% for key, m in r.metrics.items() %}
    <tr>
      <td class="metric">{{ key }}</td>
      <td>{{ m.direction }}</td>
      <td>{{ m.unit }}</td>
      <td class="num">{{ m.current.fail | default(none) | num(m.unit) }}</td>
      <td class="num">{{ m.current.warn | default(none) | num(m.unit) }}</td>
      <td class="num">{% if m.candidate %}{{ m.candidate.fail | num(m.unit) }}{% else %}—{% endif %}</td>
      <td class="num">{% if m.candidate %}{{ m.candidate.warn | num(m.unit) }}{% else %}—{% endif %}</td>
      <td class="num">{{ m.n_values }}</td>
      <td class="num">{{ m.n_sequencing_groups }}</td>
      <td class="num">{{ m.n_datasets }}</td>
      <td class="num">{{ m.n_dropped }}</td>
    </tr>
    {% endfor %}
  </table>
  </div>

  <h2>Dataset-relative (MAD) tiers</h2>
  <p class="blurb">
    A relative tier ships only <code>direction</code>, <code>k</code> and <code>min_samples</code> — production
    recomputes median and MAD from whatever sequencing groups are in each run, so the per-dataset thresholds below
    are illustrative of today's data and will move as datasets grow. That movement is what churn measures.
    Growth churn is a forecast; merge churn is a stress test, which is why its bar is looser. Note the tension: a
    metric earns a relative tier precisely because its normal level <em>shifts between datasets</em>, and the merge
    simulation punishes exactly that property — read a high figure as "how much would pooling two projects disturb
    this", not as a defect count.
  </p>

  {% for key, rel in r.relative.items() %}
  <h3>{{ key }} <span class="verdict {{ rel.verdict }}">{{ rel.verdict }}</span></h3>
  {% if rel.reason %}<p class="blurb">{{ rel.reason }}</p>{% endif %}
  <table>
    <tr><th>Measure</th><th>Observed</th><th>Bar</th></tr>
    <tr><td>Peak warn rate</td>
      <td class="num {% if rel.max_warn_rate > rel.bars.max_warn_rate %}over{% endif %}">
        {{ rel.max_warn_rate | pct }}</td>
      <td class="num">{{ rel.bars.max_warn_rate | pct }}</td></tr>
    <tr><td>Peak growth churn</td>
      <td class="num {% if rel.max_growth_churn > rel.bars.max_growth_churn %}over{% endif %}">
        {{ rel.max_growth_churn | pct }}</td>
      <td class="num">{{ rel.bars.max_growth_churn | pct }}</td></tr>
    <tr><td>Peak merge churn</td>
      <td class="num {% if rel.max_merge_churn > rel.bars.max_merge_churn %}over{% endif %}">
        {{ rel.max_merge_churn | pct }}</td>
      <td class="num">{{ rel.bars.max_merge_churn | pct }}</td></tr>
  </table>
  <div class="scroll">
  <table>
    <tr><th>Dataset</th><th>Values</th><th>Groups</th><th>Median</th><th>MAD</th>
        <th>Threshold</th><th>Warned</th><th>Warn rate</th><th>Skipped</th></tr>
    {% for d in rel.datasets %}
    <tr>
      <td>{{ d.dataset }}</td>
      <td class="num">{{ d.n_values }}</td>
      <td class="num">{{ d.n_sequencing_groups }}</td>
      <td class="num">{{ d.median | stat }}</td>
      <td class="num">{{ d.mad | stat }}</td>
      <td class="num">{{ d.threshold | stat }}</td>
      <td class="num">{% if d.skipped %}—{% else %}{{ d.n_warn }}{% endif %}</td>
      <td class="num">{% if d.skipped %}—{% else %}{{ d.warn_rate | pct }}{% endif %}</td>
      <td class="skip">{{ d.skipped or '' }}</td>
    </tr>
    {% endfor %}
  </table>
  </div>
  {% else %}
  <p class="blurb muted">No metric is configured with <code>relative = true</code> for this sequencing type.</p>
  {% endfor %}

  <h2>Evidence</h2>

  <details>
    <summary>Metric presence</summary>
    <p class="blurb">Which general-stats section carried each metric. An empty cell means the key is absent or
    renamed in that dataset — a metric absent everywhere is a key that checks nothing.</p>
    <div class="scroll">
    <table>
      <tr><th>Metric</th>{% for name in dataset_names %}<th>{{ name }}</th>{% endfor %}</tr>
      {% for key, m in r.metrics.items() %}
      <tr>
        <td class="metric">{{ key }}</td>
        {% for name in dataset_names %}
        <td>{% if name in m.missing_from %}<span class="over">MISSING</span>
            {% else %}{{ m.present_in | join(', ') }}{% endif %}</td>
        {% endfor %}
      </tr>
      {% endfor %}
    </table>
    </div>
  </details>

  <details>
    <summary>Percentile distributions</summary>
    {% for key, m in r.metrics.items() %}
    <h3>{{ key }} <span class="muted">({{ m.direction }}, {{ m.unit }})</span></h3>
    <div class="scroll">
    <table>
      <tr><th>Dataset</th>{% for label in percentile_labels %}<th>{{ label }}</th>{% endfor %}</tr>
      {% for name, pcts in m.percentiles.items() %}
      <tr><td>{{ name }}</td>
        {% for label in percentile_labels %}<td class="num">{{ pcts[label] | num(m.unit) }}</td>{% endfor %}
      </tr>
      {% endfor %}
    </table>
    </div>
    {% endfor %}
  </details>

  <details>
    <summary>Flag rates</summary>
    <p class="blurb">Warn excludes sequencing groups already failing, matching production, which evaluates fail
    before warn and records one flag per metric at the worst tier.</p>
    {% for key, m in r.metrics.items() %}
    <h3>{{ key }}</h3>
    <div class="scroll">
    <table>
      <tr><th>Dataset</th><th>Current fail</th><th>Current warn</th><th>Candidate fail</th><th>Candidate warn</th></tr>
      {% for name in m.flag_rates.candidate %}
      <tr>
        <td>{{ name }}</td>
        <td class="num">{{ m.flag_rates.current.get(name, {}).get('fail') | pct }}</td>
        <td class="num">{{ m.flag_rates.current.get(name, {}).get('warn') | pct }}</td>
        <td class="num">{{ m.flag_rates.candidate[name].fail | pct }}</td>
        <td class="num">{{ m.flag_rates.candidate[name].warn | pct }}</td>
      </tr>
      {% endfor %}
    </table>
    </div>
    {% endfor %}
  </details>

  <details>
    <summary>Churn detail</summary>
    {% for key, rel in r.relative.items() %}
    <h3>{{ key }} — growth</h3>
    <p class="blurb">A 60% before-slice re-scored against the whole dataset's threshold, run on the leading slice and
    on a seeded shuffle. The verdict uses the worse of the two, because whether MultiQC key order tracks sequencing
    batches is plausible but nowhere guaranteed.</p>
    <div class="scroll">
    <table>
      <tr><th>Dataset</th><th>Before n</th><th>Ordered flips</th><th>Ordered churn</th>
          <th>Shuffled flips</th><th>Shuffled churn</th><th>Worse</th><th>Note</th></tr>
      {% for g in rel.growth %}
      <tr>
        <td>{{ g.dataset }}</td>
        <td class="num">{{ (g.ordered or g.shuffled).n_initial }}</td>
        <td class="num">{% if g.ordered %}{{ g.ordered.flips }}{% else %}—{% endif %}</td>
        <td class="num">{% if g.ordered %}{{ g.ordered.flip_rate | pct }}{% else %}—{% endif %}</td>
        <td class="num">{% if g.shuffled %}{{ g.shuffled.flips }}{% else %}—{% endif %}</td>
        <td class="num">{% if g.shuffled %}{{ g.shuffled.flip_rate | pct }}{% else %}—{% endif %}</td>
        <td class="num">{{ g.worse | pct }}</td>
        <td>{% if g.ordering_sensitive %}<span class="over">ORDERING-SENSITIVE</span>{% endif %}</td>
      </tr>
      {% endfor %}
    </table>
    </div>
    <h3>{{ key }} — merge</h3>
    <p class="blurb">Worst {{ rel.merge_worst | length }} of {{ rel.merge_pairs_simulated }} ordered pairs. Both
    directions of every pair run, because a merge disturbs both projects' flag sets by different amounts.</p>
    <div class="scroll">
    <table>
      <tr><th>Dataset</th><th>Merged with</th><th>n</th><th>Flagged before</th>
          <th>Flagged after</th><th>Flips</th><th>Churn</th></tr>
      {% for m in rel.merge_worst %}
      <tr>
        <td>{{ m.dataset }}</td><td>{{ m.merged_with }}</td>
        <td class="num">{{ m.n_initial }}</td><td class="num">{{ m.flagged_before }}</td>
        <td class="num">{{ m.flagged_after }}</td><td class="num">{{ m.flips }}</td>
        <td class="num">{{ m.flip_rate | pct }}</td>
      </tr>
      {% endfor %}
    </table>
    </div>
    {% endfor %}
  </details>

  <details>
    <summary>Provenance</summary>
    <div class="scroll">
    <table>
      <tr><th>Dataset</th><th>Analysis ID</th><th>Completed</th><th>Groups</th><th>MultiQC</th><th>Report</th></tr>
      {% for d in r.datasets %}
      <tr>
        <td>{{ d.dataset }}</td><td class="num">{{ d.analysis_id }}</td><td>{{ d.timestamp }}</td>
        <td class="num">{{ d.n_sequencing_groups }}</td><td>{{ d.multiqc_version }}</td>
        <td><code>{{ d.uri }}</code></td>
      </tr>
      {% endfor %}
    </table>
    </div>
    {% if r.skipped_datasets %}
    <h3>Datasets not included</h3>
    <table>
      <tr><th>Dataset</th><th>Reason</th></tr>
      {% for d in r.skipped_datasets %}
      <tr><td>{{ d.dataset }}</td><td class="skip">{{ d.reason }}</td></tr>
      {% endfor %}
    </table>
    {% endif %}
  </details>

  <h2>Config block</h2>
  <p class="blurb">The candidate thresholds as a <code>qc_thresholds</code> block. Diff this against
  <code>config_template.toml</code> before pasting — that file carries comments and judgement no generator
  reproduces, and the diff is the cheapest guard against silently rewriting a production gate.</p>
  <pre>{{ r.config_snippet }}</pre>

</div>
</body>
</html>
```

- [ ] **Step 5: Run the tests to verify they pass**

Run: `uv run pytest test/test_qc_calibration_render.py -v`

Expected: PASS, 10 tests.

If `test_carries_every_expected_section` fails on `Metric presence`, check that the
`<summary>` text matches the assertion exactly.

- [ ] **Step 6: Eyeball the output**

```bash
uv run python -c "
from test.test_qc_calibration_summary import SETTINGS, dataset_values
from align_genotype.qc_calibration import summary, render
built = summary.build(
    [dataset_values('ds-a', [30,32,34,36,38,10], [10,10.5,11,11.5,12,40]),
     dataset_values('ds-b', [40,42,44,46,48,50], [7,7.5,8,8.5,9,9.5])],
    SETTINGS, current={'MEDIAN_COVERAGE': {'fail': 15, 'warn': 25}},
    skipped_datasets=[], generated='2026-08-12T00:00:00', ar_guid='local')
open('/tmp/calibration.html','w').write(render.render(built))
print('wrote /tmp/calibration.html')
"
open /tmp/calibration.html
```

Check the tables line up, the verdict pill is coloured, and the collapsible sections open.

- [ ] **Step 7: Commit**

```bash
git add src/align_genotype/qc_calibration/render.py src/align_genotype/templates/qc_calibration_report.html.jinja test/test_qc_calibration_render.py
git commit -m "feat(qc_calibration): render the calibration dashboard

Headline thresholds and MAD verdicts first, evidence in collapsible sections,
and the pasteable config block last. Autoescaped: dataset names and skip
reasons reach the page unmodified."
```

---

## Task 12: Slim `discovery.py` to a per-dataset Metamist lookup

Two changes. The dataset set now comes from the multicohort, so the `myProjects` query
and the `is_seqr` plus name-substring eligibility filter are gone entirely — dataset
selection is the operator's, expressed through `input_cohorts`. And a real bug is fixed:
`CramMultiQC` and `GvcfMultiQC` are both registered `analysis_type='qc'`, each producing
two analyses (`json` and `html`), so the current newest-wins selection can return the
GVCF report or an HTML file. CPG Flow merges `get_job_attrs()` into analysis meta, which
carries the stage name, so the filter can be server-side.

**Files:**
- Modify: `src/align_genotype/qc_calibration/discovery.py` (full rewrite)
- Replace: `test/test_qc_calibration_discovery.py` (full rewrite)

- [ ] **Step 1: Write the failing test**

Replace the whole of `test/test_qc_calibration_discovery.py` with:

```python
"""Unit tests for per-dataset CramMultiQC report discovery, using a fake query function."""

import pytest

from align_genotype.qc_calibration import discovery as discovery_mod


def analysis(analysis_id, output, timestamp) -> dict:
    return {'id': analysis_id, 'output': output, 'timestampCompleted': timestamp}


def fake_query(analyses):
    """A query function returning `analyses`, recording the variables it was given."""
    calls = []

    def _query(query_text, variables=None):  # noqa: ANN202, ARG001
        calls.append(variables)
        return {'project': {'analyses': analyses}}

    _query.calls = calls
    return _query


def test_returns_the_newest_completed_report():
    query = fake_query(
        [
            analysis(1, 'gs://a/old/multiqc_data.json', '2026-01-01T00:00:00'),
            analysis(2, 'gs://a/new/multiqc_data.json', '2026-06-01T00:00:00'),
        ],
    )
    found = discovery_mod.latest_cram_multiqc('ds-a', 'genome', query_fn=query)
    assert found.analysis_id == 2
    assert found.uri == 'gs://a/new/multiqc_data.json'
    assert found.timestamp == '2026-06-01T00:00:00'


def test_selects_on_timestamp_not_analysis_id():
    """A higher id with an older timestamp must lose."""
    query = fake_query(
        [
            analysis(99, 'gs://a/old/multiqc_data.json', '2026-01-01T00:00:00'),
            analysis(2, 'gs://a/new/multiqc_data.json', '2026-06-01T00:00:00'),
        ],
    )
    assert discovery_mod.latest_cram_multiqc('ds-a', 'genome', query_fn=query).analysis_id == 2


def test_filters_server_side_on_stage_and_sequencing_type():
    """CramMultiQC and GvcfMultiQC are both analysis_type='qc'."""
    query = fake_query([analysis(1, 'gs://a/multiqc_data.json', '2026-01-01T00:00:00')])
    discovery_mod.latest_cram_multiqc('ds-a', 'genome', query_fn=query)
    assert query.calls[0] == {
        'dataset': 'ds-a',
        'analysisType': 'qc',
        'metaFilter': {'stage': {'eq': 'CramMultiQC'}, 'sequencing_type': {'eq': 'genome'}},
    }


def test_ignores_the_html_analysis_of_the_same_stage():
    """analysis_keys=['json','html'] registers two analyses per run."""
    query = fake_query(
        [
            analysis(1, 'gs://a/multiqc.html', '2026-06-02T00:00:00'),
            analysis(2, 'gs://a/multiqc_data.json', '2026-06-01T00:00:00'),
        ],
    )
    found = discovery_mod.latest_cram_multiqc('ds-a', 'genome', query_fn=query)
    assert found.uri == 'gs://a/multiqc_data.json'


def test_no_matching_analysis_returns_none():
    assert discovery_mod.latest_cram_multiqc('ds-a', 'genome', query_fn=fake_query([])) is None


def test_only_html_analyses_returns_none():
    query = fake_query([analysis(1, 'gs://a/multiqc.html', '2026-06-01T00:00:00')])
    assert discovery_mod.latest_cram_multiqc('ds-a', 'genome', query_fn=query) is None


def test_analysis_without_an_output_is_skipped():
    query = fake_query(
        [
            analysis(1, None, '2026-06-02T00:00:00'),
            analysis(2, 'gs://a/multiqc_data.json', '2026-06-01T00:00:00'),
        ],
    )
    assert discovery_mod.latest_cram_multiqc('ds-a', 'genome', query_fn=query).analysis_id == 2


def test_unrankable_timestamp_is_skipped_with_a_warning(caplog):
    """A row we cannot order must be dropped, not compared - and never silently."""
    query = fake_query(
        [
            analysis(1, 'gs://a/broken/multiqc_data.json', None),
            analysis(2, 'gs://a/good/multiqc_data.json', '2026-01-01T00:00:00'),
        ],
    )
    with caplog.at_level('WARNING'):
        found = discovery_mod.latest_cram_multiqc('ds-a', 'genome', query_fn=query)
    assert found.analysis_id == 2
    assert 'timestampCompleted' in caplog.text


def test_a_null_project_returns_none():
    """`project` can be present-but-null, e.g. with no read access."""

    def _query(query_text, variables=None):  # noqa: ANN202, ARG001
        return {'project': None}

    assert discovery_mod.latest_cram_multiqc('ds-a', 'genome', query_fn=_query) is None
```

- [ ] **Step 2: Run the tests to verify they fail**

Run: `uv run pytest test/test_qc_calibration_discovery.py -v`

Expected: FAIL — `AttributeError: module ... has no attribute 'latest_cram_multiqc'`.

- [ ] **Step 3: Write the implementation**

Replace the whole of `src/align_genotype/qc_calibration/discovery.py` with:

```python
"""Find each dataset's latest CramMultiQC report in Metamist.

`CramMultiQC` and `GvcfMultiQC` both register `analysis_type='qc'`, and each registers
one analysis per entry in `analysis_keys` - so a plain type-and-sequencing-type query
returns GVCF reports and HTML files alongside the CRAM JSON we want. CPG Flow merges the
stage name into analysis meta, so the stage filter runs server-side and the output suffix
picks the JSON of the two.

The query function is injected so the selection logic is testable with a fake, leaving
only the four-line adapter untested.
"""

import functools
import logging
from collections.abc import Callable
from dataclasses import dataclass
from typing import Any

from metamist.graphql import gql, query

# The CramMultiQC output we want; the stage also registers its HTML under the same type.
JSON_SUFFIX = 'multiqc_data.json'
CRAM_MULTIQC_STAGE = 'CramMultiQC'

QueryFn = Callable[[str, dict[str, Any] | None], dict[str, Any]]

ANALYSES_QUERY = gql(
    """
    query CramMultiqc($dataset: String!, $analysisType: String!, $metaFilter: JSON!) {
        project(name: $dataset) {
            analyses(status: {eq: COMPLETED}, type: {eq: $analysisType}, meta: $metaFilter) {
                id
                output
                timestampCompleted
            }
        }
    }
    """,
)


@dataclass(frozen=True)
class MultiqcReport:
    """One dataset's MultiQC report, and the analysis it was registered under."""

    dataset: str
    uri: str
    analysis_id: int
    timestamp: str


def default_query(query_text: str, variables: dict[str, Any] | None = None) -> dict[str, Any]:
    return query(query_text, variables=variables or {})


def latest_cram_multiqc(
    dataset: str,
    seq_type: str,
    query_fn: QueryFn = default_query,
) -> MultiqcReport | None:
    """The newest completed CramMultiQC JSON for `dataset`, or None if there is none."""
    result = query_fn(
        ANALYSES_QUERY,
        {
            'dataset': dataset,
            'analysisType': 'qc',
            'metaFilter': {'stage': {'eq': CRAM_MULTIQC_STAGE}, 'sequencing_type': {'eq': seq_type}},
        },
    )
    # `project` can be present-but-null - no read access, say - so `.get('project', {})`
    # is not enough: the null wins over the default and the next `.get` crashes on None.
    analyses = (result.get('project') or {}).get('analyses') or []

    candidates = []
    for analysis in analyses:
        output = analysis.get('output')
        if not output or not str(output).endswith(JSON_SUFFIX):
            continue
        timestamp = analysis.get('timestampCompleted')
        if not isinstance(timestamp, str):
            # Cannot be ordered against the others. Dropping a candidate we cannot rank
            # beats picking the wrong one, and the cost is at worst a slightly older
            # report - inspectable, because the report records what it used.
            logging.warning(
                f'{dataset}: analysis {analysis.get("id")!r} has no usable timestampCompleted '
                f'({timestamp!r}); excluding it from latest-report selection',
            )
            continue
        candidates.append(analysis)

    if not candidates:
        return None
    newest = max(candidates, key=lambda a: a['timestampCompleted'])
    return MultiqcReport(
        dataset=dataset,
        uri=str(newest['output']),
        analysis_id=int(newest['id']),
        timestamp=newest['timestampCompleted'],
    )


@functools.lru_cache(maxsize=None)
def cached_latest_cram_multiqc(dataset: str, seq_type: str) -> MultiqcReport | None:
    """`latest_cram_multiqc`, memoised for the driver.

    `expected_outputs` keys the per-dataset output path on the analysis ID and is called
    repeatedly during DAG assembly, so without this each call would be a GraphQL round
    trip.
    """
    return latest_cram_multiqc(dataset, seq_type)
```

- [ ] **Step 4: Run the tests to verify they pass**

Run: `uv run pytest test/test_qc_calibration_discovery.py -v`

Expected: PASS, 9 tests.

- [ ] **Step 5: Commit**

```bash
git add src/align_genotype/qc_calibration/discovery.py test/test_qc_calibration_discovery.py
git commit -m "fix(qc_calibration): select the CRAM MultiQC report, not any qc analysis

CramMultiQC and GvcfMultiQC both register analysis_type='qc', each with a json
and an html analysis, so newest-wins could return the GVCF report or an HTML
file - neither carrying any CRAM metric. Filters server-side on the stage name
cpg-flow writes into analysis meta. Dataset eligibility filtering goes: the
dataset set now comes from the multicohort."
```

---

## Task 13: `scripts/qc_calibration_extract.py` — the per-dataset job entrypoint

**Files:**
- Create: `src/align_genotype/scripts/qc_calibration_extract.py`
- Test: `test/test_qc_calibration_scripts.py`

- [ ] **Step 1: Write the failing test**

Create `test/test_qc_calibration_scripts.py`:

```python
"""Unit tests for the two job entrypoint scripts, run in-process."""

import json

import pytest
from click.testing import CliRunner

from align_genotype.qc_calibration import settings as settings_mod
from align_genotype.qc_calibration import values as values_mod
from align_genotype.scripts import qc_calibration_extract

CALIBRATION_CONFIG = {
    'workflow': {'sequencing_type': 'genome'},
    'qc_calibration': {
        'genome': {
            'metrics': {
                'MEDIAN_COVERAGE': {'direction': 'min', 'unit': 'x'},
                'dup_pct': {'direction': 'max', 'unit': '%', 'relative': True},
            },
        },
    },
    'qc_thresholds': {'genome': {'fail': {'min': {'MEDIAN_COVERAGE': 15}}}},
}

REPORT = {
    'config_version': '1.33',
    'report_general_stats_data': {
        'picard_1': {f'CPG{i}': {'MEDIAN_COVERAGE': 30.0 + i} for i in range(6)},
        'samtools': {f'CPG{i}': {'dup_pct': 10.0 + i} for i in range(6)},
    },
}


@pytest.fixture
def patch_config(monkeypatch):
    """Point every config_retrieve used by the scripts at CALIBRATION_CONFIG."""

    def config_retrieve(keys, default=None):  # noqa: ANN202
        node = CALIBRATION_CONFIG
        for key in keys:
            if not isinstance(node, dict) or key not in node:
                return default
            node = node[key]
        return node

    from align_genotype.scripts import check_multiqc  # noqa: PLC0415

    monkeypatch.setattr(settings_mod.config, 'config_retrieve', config_retrieve)
    monkeypatch.setattr(check_multiqc.config, 'config_retrieve', config_retrieve)


def test_extract_writes_a_values_file(tmp_path, patch_config):  # noqa: ARG001
    report_path = tmp_path / 'multiqc_data.json'
    report_path.write_text(json.dumps(REPORT))
    output = tmp_path / 'values.json'

    result = CliRunner().invoke(
        qc_calibration_extract.main,
        [
            '--dataset', 'ds-a',
            '--multiqc-json', str(report_path),
            '--analysis-id', '42',
            '--timestamp', '2026-06-01T00:00:00',
            '--uri', 'gs://ds-a/multiqc_data.json',
            '--output', str(output),
        ],
    )
    assert result.exit_code == 0, result.output

    loaded = values_mod.load(output)
    assert loaded.dataset == 'ds-a'
    assert loaded.analysis_id == 42
    assert loaded.uri == 'gs://ds-a/multiqc_data.json'
    assert loaded.n_sequencing_groups == 6
    assert loaded.metric('MEDIAN_COVERAGE').n_values == 6


def test_extract_fails_loudly_on_an_unreadable_report(tmp_path, patch_config):  # noqa: ARG001
    report_path = tmp_path / 'multiqc_data.json'
    report_path.write_text(json.dumps({'config_version': '1.33'}))

    result = CliRunner().invoke(
        qc_calibration_extract.main,
        [
            '--dataset', 'ds-a',
            '--multiqc-json', str(report_path),
            '--analysis-id', '42',
            '--timestamp', '2026-06-01T00:00:00',
            '--uri', 'gs://ds-a/multiqc_data.json',
            '--output', str(tmp_path / 'values.json'),
        ],
    )
    assert result.exit_code != 0
    assert 'no usable report_general_stats_data' in str(result.exception)
```

- [ ] **Step 2: Run the tests to verify they fail**

Run: `uv run pytest test/test_qc_calibration_scripts.py -v`

Expected: FAIL — `ImportError: cannot import name 'qc_calibration_extract'`.

- [ ] **Step 3: Write the implementation**

Create `src/align_genotype/scripts/qc_calibration_extract.py`:

```python
"""Distil one dataset's MultiQC report into a values file.

Run as a Hail Batch job by `QcCalibrationDatasetMetrics`, with the report already
localised by `batch.read_input`. Also runnable off a checkout against a downloaded
report, which is the quickest way to check what a new MultiQC version surfaces.
"""

import json
import logging
from datetime import datetime, timezone

import click

from cpg_utils import to_path

from align_genotype.qc_calibration import extract, settings, values

logging.basicConfig(level=logging.INFO, format='%(levelname)s: %(message)s')


@click.command()
@click.option('--dataset', required=True, help='Dataset name, used to label the values file.')
@click.option('--multiqc-json', 'multiqc_json_path', required=True, help='Path to the MultiQC JSON.')
@click.option('--analysis-id', type=int, required=True, help='Metamist analysis ID the report came from.')
@click.option('--timestamp', required=True, help='timestampCompleted of that analysis.')
@click.option('--uri', required=True, help='The report URI, recorded for provenance.')
@click.option('--output', 'output_path', required=True, help='Where to write the values file.')
def main(
    dataset: str,
    multiqc_json_path: str,
    analysis_id: int,
    timestamp: str,
    uri: str,
    output_path: str,
) -> None:
    """Extract every configured calibration metric from one MultiQC report."""
    loaded = settings.load()
    logging.info(f'{dataset}: extracting {len(loaded.metrics)} {loaded.seq_type} metrics from {multiqc_json_path}')

    with to_path(multiqc_json_path).open() as f:
        document = json.load(f)

    result = extract.extract(
        document,
        loaded,
        dataset=dataset,
        analysis_id=analysis_id,
        timestamp=timestamp,
        uri=uri,
        generated=datetime.now(tz=timezone.utc).isoformat(timespec='seconds'),
    )
    values.save(result, output_path)

    for key, metric_values in result.metrics.items():
        if not metric_values.entries:
            logging.warning(f'{dataset}: {key} is absent from this report')
        elif metric_values.n_dropped:
            logging.warning(f'{dataset}: {key} lost {metric_values.n_dropped} unusable value(s)')
    logging.info(f'{dataset}: wrote {output_path} ({result.n_sequencing_groups} sequencing groups)')


if __name__ == '__main__':
    main()
```

- [ ] **Step 4: Run the tests to verify they pass**

Run: `uv run pytest test/test_qc_calibration_scripts.py -v`

Expected: PASS, 2 tests.

- [ ] **Step 5: Commit**

```bash
git add src/align_genotype/scripts/qc_calibration_extract.py test/test_qc_calibration_scripts.py
git commit -m "feat(qc_calibration): add the per-dataset extract job entrypoint"
```

---

## Task 14: `scripts/qc_calibration_report.py` — the report job entrypoint

**Files:**
- Create: `src/align_genotype/scripts/qc_calibration_report.py`
- Test: `test/test_qc_calibration_scripts.py` (append)

- [ ] **Step 1: Write the failing test**

Append to `test/test_qc_calibration_scripts.py`:

```python
def _write_values(tmp_path, name, coverage, dup):
    from align_genotype.qc_calibration import values as v  # noqa: PLC0415

    path = tmp_path / f'{name}.json'
    v.save(
        v.DatasetValues(
            dataset=name,
            seq_type='genome',
            analysis_id=1,
            timestamp='2026-06-01T00:00:00',
            uri=f'gs://{name}/multiqc_data.json',
            multiqc_version='1.33',
            generated='2026-08-12T00:00:00',
            n_sequencing_groups=len(coverage),
            section_sizes={'picard_1': len(coverage)},
            metrics={
                'MEDIAN_COVERAGE': v.MetricValues(
                    entries=tuple(('picard_1', f'{name}-CPG{i}', float(x)) for i, x in enumerate(coverage)),
                    n_dropped=0,
                ),
                'dup_pct': v.MetricValues(
                    entries=tuple(('samtools', f'{name}-CPG{i}', float(x)) for i, x in enumerate(dup)),
                    n_dropped=0,
                ),
            },
        ),
        path,
    )
    return path


def test_report_writes_json_and_html(tmp_path, patch_config, monkeypatch):  # noqa: ARG001
    from align_genotype.scripts import qc_calibration_report  # noqa: PLC0415

    monkeypatch.setattr(qc_calibration_report.config, 'try_get_ar_guid', lambda: 'test-ar-guid')
    a = _write_values(tmp_path, 'ds-a', [30, 32, 34, 36, 38, 10], [10, 10.5, 11, 11.5, 12, 40])
    b = _write_values(tmp_path, 'ds-b', [40, 42, 44, 46, 48, 50], [7, 7.5, 8, 8.5, 9, 9.5])
    out_json = tmp_path / 'calibration.json'
    out_html = tmp_path / 'calibration.html'

    result = CliRunner().invoke(
        qc_calibration_report.main,
        [
            '--values', str(a), '--values', str(b),
            '--skipped-dataset', 'ds-c',
            '--output-json', str(out_json),
            '--output-html', str(out_html),
        ],
    )
    assert result.exit_code == 0, result.output

    payload = json.loads(out_json.read_text())
    assert payload['sequencing_type'] == 'genome'
    assert [d['dataset'] for d in payload['datasets']] == ['ds-a', 'ds-b']
    assert payload['skipped_datasets'] == [
        {'dataset': 'ds-c', 'reason': 'no completed CramMultiQC qc analysis for genome'},
    ]
    assert payload['metrics']['MEDIAN_COVERAGE']['current'] == {'fail': 15}
    assert out_html.read_text().startswith('<!DOCTYPE html>')


def test_report_succeeds_even_when_a_metric_is_missing_everywhere(tmp_path, patch_config, monkeypatch, caplog):  # noqa: ARG001
    """Failing would destroy the HTML that explains the problem - Hail only copies
    write_output targets on job success."""
    from align_genotype.scripts import qc_calibration_report  # noqa: PLC0415

    monkeypatch.setattr(qc_calibration_report.config, 'try_get_ar_guid', lambda: 'x')
    path = tmp_path / 'ds-a.json'
    from align_genotype.qc_calibration import values as v  # noqa: PLC0415

    v.save(
        v.DatasetValues(
            dataset='ds-a', seq_type='genome', analysis_id=1, timestamp='t',
            uri='gs://a/multiqc_data.json', multiqc_version='1.33', generated='g',
            n_sequencing_groups=0, section_sizes={}, metrics={},
        ),
        path,
    )
    out_html = tmp_path / 'calibration.html'
    with caplog.at_level('ERROR'):
        result = CliRunner().invoke(
            qc_calibration_report.main,
            [
                '--values', str(path),
                '--output-json', str(tmp_path / 'calibration.json'),
                '--output-html', str(out_html),
            ],
        )
    assert result.exit_code == 0, result.output
    assert out_html.exists()
    assert 'checks nothing' in caplog.text
```

- [ ] **Step 2: Run the tests to verify they fail**

Run: `uv run pytest test/test_qc_calibration_scripts.py -v`

Expected: FAIL — `ImportError: cannot import name 'qc_calibration_report'`.

- [ ] **Step 3: Write the implementation**

Create `src/align_genotype/scripts/qc_calibration_report.py`:

```python
"""Turn every dataset's values file into the calibration report.

Run as a Hail Batch job by `QcCalibrationReport`, with the values files already localised
by `batch.read_input`. Also runnable off a checkout against downloaded values files,
which is the quickest way to try a different `k` or metric list without re-parsing any
MultiQC report:

    python -m align_genotype.scripts.qc_calibration_report \\
        --values a.json --values b.json --output-json out.json --output-html out.html

A metric absent from every dataset is logged at ERROR and banners the report, but does
not fail the job: Hail only copies `write_output` targets on success, so failing here
would destroy the page that explains the problem.
"""

import json
import logging
from datetime import datetime, timezone

import click

from cpg_utils import config, to_path

from align_genotype.qc_calibration import render, settings, summary, values

logging.basicConfig(level=logging.INFO, format='%(levelname)s: %(message)s')

SKIP_REASON = 'no completed CramMultiQC qc analysis for {seq_type}'


@click.command()
@click.option('--values', 'value_paths', multiple=True, required=True, help='A dataset values file. Repeatable.')
@click.option('--skipped-dataset', 'skipped', multiple=True, help='A dataset with no report. Repeatable.')
@click.option('--output-json', 'output_json', required=True, help='Where to write calibration.json.')
@click.option('--output-html', 'output_html', required=True, help='Where to write calibration.html.')
def main(value_paths: tuple[str, ...], skipped: tuple[str, ...], output_json: str, output_html: str) -> None:
    """Assemble the cross-dataset calibration report."""
    loaded = settings.load()
    datasets = [values.load(path) for path in value_paths]
    logging.info(f'Read {len(datasets)} dataset values file(s) for {loaded.seq_type}')

    built = summary.build(
        datasets,
        loaded,
        current=settings.current_thresholds(loaded.seq_type),
        skipped_datasets=[
            {'dataset': name, 'reason': SKIP_REASON.format(seq_type=loaded.seq_type)} for name in skipped
        ],
        generated=datetime.now(tz=timezone.utc).isoformat(timespec='seconds'),
        ar_guid=config.try_get_ar_guid() or 'unknown',
    )

    with to_path(output_json).open('w') as f:
        json.dump(built, f, indent=2, allow_nan=False)
    with to_path(output_html).open('w') as f:
        f.write(render.render(built))

    for warning in built['warnings']:
        # 'checks nothing' marks the PCT_PF_READS_ALIGNED class: a key that gates nothing
        # at all, as opposed to one merely absent from some datasets.
        logging.error(warning) if 'checks nothing' in warning else logging.warning(warning)
    for key, evaluation in built['relative'].items():
        logging.info(f'{key}: {evaluation["verdict"]} {evaluation["reason"]}')
    logging.info(f'Wrote {output_json} and {output_html}')


if __name__ == '__main__':
    main()
```

- [ ] **Step 4: Run the tests to verify they pass**

Run: `uv run pytest test/test_qc_calibration_scripts.py -v`

Expected: PASS, 4 tests.

- [ ] **Step 5: Commit**

```bash
git add src/align_genotype/scripts/qc_calibration_report.py test/test_qc_calibration_scripts.py
git commit -m "feat(qc_calibration): add the report job entrypoint

A metric missing everywhere is logged at ERROR and banners the report but does
not fail the job - Hail only copies write_output targets on success, so failing
would destroy the page that explains the problem."
```

---

## Task 15: `jobs/qc_calibration.py` — the two job builders

**Files:**
- Create: `src/align_genotype/jobs/qc_calibration.py`

No unit test: like every other module under `jobs/`, this only assembles Hail Batch
objects, and asserting on a `Job`'s command string tests the f-string rather than the
behaviour. It is exercised by the dry run in Task 18.

- [ ] **Step 1: Write the implementation**

Create `src/align_genotype/jobs/qc_calibration.py`:

```python
"""Batch jobs for QC threshold calibration."""

from hailtop.batch.job import Job

from cpg_utils import Path, config, hail_batch

from align_genotype.qc_calibration.discovery import MultiqcReport


def extract_dataset_metrics(report: MultiqcReport, output: Path, job_attrs: dict) -> Job:
    """Localise one dataset's MultiQC report and distil it to a values file.

    `read_input` rather than a `gcloud storage cp` in the command, matching every other
    job here. Reports run to hundreds of megabytes and parse into several GB of Python
    objects, so this asks for highmem and enough disk for one localised report.
    """
    batch = hail_batch.get_batch()

    job = batch.new_bash_job(f'QC calibration extract: {report.dataset}', job_attrs | {'tool': 'python'})
    job.image(config.config_retrieve(['workflow', 'driver_image']))
    job.cpu(2).memory('highmem').storage(config.config_retrieve(['qc_calibration', 'extract_storage'], '20Gi'))

    localised = batch.read_input(report.uri)

    job.command(
        f"""\
    python3 -m align_genotype.scripts.qc_calibration_extract \\
        --dataset {report.dataset} \\
        --multiqc-json {localised} \\
        --analysis-id {report.analysis_id} \\
        --timestamp {report.timestamp} \\
        --uri {report.uri} \\
        --output {job.values}
    """
    )

    batch.write_output(job.values, output)
    return job


def calibration_report(
    values_paths: dict[str, Path],
    skipped_datasets: list[str],
    outputs: dict[str, Path],
    job_attrs: dict,
) -> Job:
    """Read every dataset's values file and render the calibration report.

    Small by construction: the values files are a few hundred KB each, because the
    expensive parse already happened one dataset at a time in `extract_dataset_metrics`.
    """
    batch = hail_batch.get_batch()

    job = batch.new_bash_job('QC calibration report', job_attrs | {'tool': 'python'})
    job.image(config.config_retrieve(['workflow', 'driver_image']))
    job.cpu(2).memory('standard')

    localised = [batch.read_input(str(path)) for path in values_paths.values()]
    values_args = ' '.join(f'--values {resource}' for resource in localised)
    skipped_args = ' '.join(f'--skipped-dataset {name}' for name in skipped_datasets)

    job.command(
        f"""\
    python3 -m align_genotype.scripts.qc_calibration_report \\
        {values_args} {skipped_args} \\
        --output-json {job.report_json} \\
        --output-html {job.report_html}
    """
    )

    batch.write_output(job.report_json, outputs['json'])
    batch.write_output(job.report_html, outputs['html'])
    return job
```

- [ ] **Step 2: Verify it imports and lints**

Run: `uv run python -c "from align_genotype.jobs import qc_calibration; print('ok')"`

Expected: `ok`

Run: `uv run ruff check src/align_genotype/jobs/qc_calibration.py`

Expected: `All checks passed!`

- [ ] **Step 3: Commit**

```bash
git add src/align_genotype/jobs/qc_calibration.py
git commit -m "feat(qc_calibration): add the extract and report job builders

Inputs via read_input and outputs via write_output rather than gcloud storage
cp, matching every other job module here."
```

---

## Task 16: The two stages, wired into `run_workflow.py`

**Files:**
- Create: `src/align_genotype/qc_calibration_stages.py`
- Modify: `src/align_genotype/run_workflow.py`
- Test: `test/test_qc_calibration_stages.py`

- [ ] **Step 1: Write the failing test**

`Stage.__init__` calls `get_workflow()`, which raises unless a `Workflow` singleton
exists, so a stage class cannot be constructed outside a real run. Every decision the
stages make therefore lives in free functions, and the classes stay thin shells — the
same shape as every other stage module here.

Create `test/test_qc_calibration_stages.py`:

```python
"""Unit tests for the calibration stages' output paths and inert-by-default behaviour.

The stage classes are not instantiated here: `Stage.__init__` calls `get_workflow()`,
which raises unless a Workflow singleton exists. The `forced` flags and the DAG wiring
are verified by the dry run in Task 18, which prints `[forced]` beside a forced stage.
"""

from types import SimpleNamespace

import pytest

from cpg_utils import to_path

from align_genotype import qc_calibration_stages as stages_mod
from align_genotype.qc_calibration.discovery import MultiqcReport

REPORT = MultiqcReport(dataset='ds-a', uri='gs://x/multiqc_data.json', analysis_id=42, timestamp='t')


def fake_dataset(name: str) -> SimpleNamespace:
    return SimpleNamespace(
        name=name,
        prefix=lambda: to_path(f'gs://cpg-{name}-main'),
        web_prefix=lambda: to_path(f'gs://cpg-{name}-web'),
    )


@pytest.fixture
def patch_lookup(monkeypatch):
    """Control what discovery returns, the run's sequencing type, and the enabled flag."""

    def _apply(report, *, seq_type='genome', enabled=True):  # noqa: ANN202
        monkeypatch.setattr(stages_mod.settings, 'enabled', lambda: enabled)
        monkeypatch.setattr(stages_mod, 'sequencing_type', lambda: seq_type)
        monkeypatch.setattr(stages_mod, 'dataset_report', lambda _name: report)

    return _apply


def test_values_path_is_keyed_on_the_analysis_id():
    """A new MultiQC report means a new path, so extraction never reads a stale file."""
    path = stages_mod.dataset_values_path(fake_dataset('ds-a'), 'genome', REPORT)
    assert str(path) == 'gs://cpg-ds-a-main/qc_calibration/genome/values.42.json'


def test_values_path_is_prefixed_by_the_literal_sequencing_type():
    """Not sequencing_subdir(), which returns '' for genome and would look asymmetric."""
    path = stages_mod.dataset_values_path(fake_dataset('ds-a'), 'exome', REPORT)
    assert '/qc_calibration/exome/' in str(path)


def test_dataset_outputs_are_empty_when_calibration_is_disabled(patch_lookup):
    patch_lookup(REPORT, enabled=False)
    assert stages_mod.dataset_outputs(fake_dataset('ds-a')) == {}


def test_dataset_outputs_are_empty_without_a_report(patch_lookup):
    """cpg-flow treats a falsy expected output as reusable, so no job is queued."""
    patch_lookup(None)
    assert stages_mod.dataset_outputs(fake_dataset('ds-a')) == {}


def test_dataset_outputs_name_the_values_file_when_enabled(patch_lookup):
    patch_lookup(REPORT)
    outputs = stages_mod.dataset_outputs(fake_dataset('ds-a'))
    assert set(outputs) == {'values'}
    assert str(outputs['values']).endswith('values.42.json')


def test_report_outputs_are_empty_when_calibration_is_disabled(patch_lookup):
    patch_lookup(REPORT, enabled=False)
    assert stages_mod.report_outputs(fake_dataset('analysis')) == {}


def test_report_writes_json_to_main_and_html_to_web(patch_lookup):
    patch_lookup(REPORT)
    outputs = stages_mod.report_outputs(fake_dataset('analysis'))
    assert str(outputs['json']) == 'gs://cpg-analysis-main/qc_calibration/genome/calibration.json'
    assert str(outputs['html']) == 'gs://cpg-analysis-web/qc_calibration/genome/calibration.html'


def test_report_outputs_follow_the_sequencing_type(patch_lookup):
    patch_lookup(REPORT, seq_type='exome')
    assert '/qc_calibration/exome/' in str(stages_mod.report_outputs(fake_dataset('analysis'))['json'])


def test_collect_values_paths_separates_datasets_with_and_without_reports(monkeypatch):
    monkeypatch.setattr(stages_mod, 'dataset_report', lambda name: REPORT if name == 'ds-a' else None)
    found, skipped = stages_mod.collect_values_paths([fake_dataset('ds-a'), fake_dataset('ds-b')], 'genome')
    assert set(found) == {'ds-a'}
    assert skipped == ['ds-b']


def test_calibration_stages_are_wired_into_the_entrypoint():
    from align_genotype import run_workflow  # noqa: PLC0415

    assert stages_mod.QcCalibrationDatasetMetrics in run_workflow.STAGES
    assert stages_mod.QcCalibrationReport in run_workflow.STAGES
```

- [ ] **Step 2: Run the tests to verify they fail**

Run: `uv run pytest test/test_qc_calibration_stages.py -v`

Expected: FAIL — `ModuleNotFoundError: No module named 'align_genotype.qc_calibration_stages'`.

- [ ] **Step 3: Write the stages**

Create `src/align_genotype/qc_calibration_stages.py`:

```python
"""QC threshold calibration stages.

An isolated branch of the DAG: neither stage declares a `required_stages` dependency on
a production stage, so adding them to `run_workflow` never pulls alignment or genotyping
work into a calibration run, and never adds calibration work to a production one.

Both stages return no jobs unless `qc_calibration.enabled` is set. That flag ships false,
because a dependency-free stage whose outputs do not exist would otherwise queue a job
per dataset and register Metamist analyses on every production invocation. A calibration
run sets it true and passes
`only_stages = ['QcCalibrationDatasetMetrics', 'QcCalibrationReport']` - both names, since
`only_stages` marks every unlisted stage skipped, so naming just the report would make
the extract stage check for outputs rather than produce them.

Every decision lives in the free functions below rather than in the stage methods:
`Stage.__init__` calls `get_workflow()`, so a stage cannot be constructed outside a real
run and its methods cannot be unit tested.
"""

from cpg_flow import stage, targets
from cpg_utils import Path, config

from align_genotype.jobs import qc_calibration
from align_genotype.qc_calibration import discovery, settings
from align_genotype.qc_calibration.discovery import MultiqcReport


def sequencing_type() -> str:
    return config.config_retrieve(['workflow', 'sequencing_type'])


def dataset_report(dataset_name: str) -> MultiqcReport | None:
    """That dataset's latest CramMultiQC report, memoised across DAG assembly.

    The Metamist project name carries a `-test` suffix at test access level, and this
    name goes into a GraphQL query, so it must be resolved rather than used raw.
    """
    return discovery.cached_latest_cram_multiqc(
        config.dataset_for_access_level(dataset_name),
        sequencing_type(),
    )


def dataset_values_path(dataset: targets.Dataset, seq_type: str, report: MultiqcReport) -> Path:
    """Where one dataset's values file lives.

    Keyed on the analysis ID, so a new MultiQC report produces a new path and extraction
    re-runs exactly when the underlying data changed - never on a stale file.

    The literal sequencing type, not `targets.sequencing_subdir()`, which returns '' for
    genome and would leave the two runs' outputs asymmetric.
    """
    return dataset.prefix() / 'qc_calibration' / seq_type / f'values.{report.analysis_id}.json'


def dataset_outputs(dataset: targets.Dataset) -> dict[str, Path]:
    """The extract stage's expected outputs, empty when there is nothing to do.

    cpg-flow treats a falsy expected output as reusable, so the action becomes REUSE,
    `queue_jobs` is never called, and no Metamist analysis is created.
    """
    if not settings.enabled():
        return {}
    report = dataset_report(dataset.name)
    if report is None:
        return {}
    return {'values': dataset_values_path(dataset, sequencing_type(), report)}


def report_outputs(analysis_dataset: targets.Dataset) -> dict[str, Path]:
    """The report stage's expected outputs: JSON to main, HTML to web."""
    if not settings.enabled():
        return {}
    seq_type = sequencing_type()
    return {
        'json': analysis_dataset.prefix() / 'qc_calibration' / seq_type / 'calibration.json',
        'html': analysis_dataset.web_prefix() / 'qc_calibration' / seq_type / 'calibration.html',
    }


def collect_values_paths(
    datasets: list[targets.Dataset],
    seq_type: str,
) -> tuple[dict[str, Path], list[str]]:
    """Split the multicohort's datasets into those with a values file and those without.

    Paths are recomputed from the same memoised discovery the extract stage used rather
    than read out of `StageInput`: `StageInput._each` raises when *no* dataset produced an
    output, which is a legitimate state here. `required_stages` still supplies the job
    ordering, which is the part `inputs` is actually needed for.
    """
    found: dict[str, Path] = {}
    skipped: list[str] = []
    for dataset in datasets:
        report = dataset_report(dataset.name)
        if report is None:
            skipped.append(dataset.name)
        else:
            found[dataset.name] = dataset_values_path(dataset, seq_type, report)
    return found, skipped


@stage.stage
class QcCalibrationDatasetMetrics(stage.DatasetStage):
    """Distil one dataset's latest CramMultiQC report to a small values file."""

    def expected_outputs(self, dataset: targets.Dataset) -> dict[str, Path]:
        return dataset_outputs(dataset)

    def queue_jobs(self, dataset: targets.Dataset, inputs: stage.StageInput) -> stage.StageOutput:  # noqa: ARG002
        outputs = self.expected_outputs(dataset)
        if not outputs:
            return self.make_outputs(dataset, data=None, jobs=None)

        job = qc_calibration.extract_dataset_metrics(
            report=dataset_report(dataset.name),
            output=outputs['values'],
            job_attrs=self.get_job_attrs(dataset),
        )
        return self.make_outputs(dataset, data=outputs, jobs=job)


@stage.stage(
    required_stages=[QcCalibrationDatasetMetrics],
    analysis_type='web',
    analysis_keys=['html'],
    forced=True,
)
class QcCalibrationReport(stage.MultiCohortStage):
    """Cross-dataset analysis and the HTML dashboard the teams review."""

    def expected_outputs(self, multicohort: targets.MultiCohort) -> dict[str, Path]:
        return report_outputs(multicohort.analysis_dataset)

    def queue_jobs(self, multicohort: targets.MultiCohort, inputs: stage.StageInput) -> stage.StageOutput:  # noqa: ARG002
        outputs = self.expected_outputs(multicohort)
        if not outputs:
            return self.make_outputs(multicohort, data=None, jobs=None)

        values_paths, skipped = collect_values_paths(multicohort.get_datasets(), sequencing_type())
        if not values_paths:
            raise ValueError(
                f'No dataset in this multicohort has a completed CramMultiQC {sequencing_type()} analysis, '
                f'so there is nothing to calibrate from. Check that CramMultiQC has run for these '
                f'datasets at this sequencing type.',
            )

        job = qc_calibration.calibration_report(
            values_paths=values_paths,
            skipped_datasets=skipped,
            outputs=outputs,
            job_attrs=self.get_job_attrs(multicohort),
        )
        return self.make_outputs(multicohort, data=outputs, jobs=job)
```

- [ ] **Step 4: Wire the stages into the entrypoint**

Rewrite `src/align_genotype/run_workflow.py` so the stage list is a module-level constant
the test can import:

```python
#!/usr/bin/env python3

from argparse import ArgumentParser

from cpg_flow.workflow import run_workflow

from align_genotype.qc_calibration_stages import QcCalibrationDatasetMetrics, QcCalibrationReport
from align_genotype.stages import (
    CramQcPicardCollectMetrics,
    CramQcPicardMultiMetrics,
    CramQcSamtoolsStats,
    CramQcSomalier,
    CramQcVerifyBamId,
    GenotypeWithGatk,
    RunGvcfQc,
    VntyperIndexPage,
)

STAGES = [
    GenotypeWithGatk,
    CramQcPicardMultiMetrics,
    CramQcPicardCollectMetrics,
    CramQcSomalier,
    CramQcSamtoolsStats,
    CramQcVerifyBamId,
    RunGvcfQc,
    VntyperIndexPage,
    # An isolated branch: no dependency on the stages above, and inert unless
    # `qc_calibration.enabled` is set. See qc_calibration/README.md.
    QcCalibrationDatasetMetrics,
    QcCalibrationReport,
]


def cli_main():
    """
    CLI entrypoint - starts up the workflow
    """
    parser = ArgumentParser()
    parser.add_argument('--dry_run', action='store_true', help='Dry run')
    args = parser.parse_args()

    run_workflow(name='align_genotype', stages=STAGES, dry_run=args.dry_run)


if __name__ == '__main__':
    cli_main()
```

- [ ] **Step 5: Run the tests to verify they pass**

Run: `uv run pytest test/test_qc_calibration_stages.py -v`

Expected: PASS, 11 tests.

- [ ] **Step 6: Lint**

Run: `uv run ruff check src/align_genotype/qc_calibration_stages.py src/align_genotype/run_workflow.py`

Expected: `All checks passed!`

- [ ] **Step 7: Commit**

```bash
git add src/align_genotype/qc_calibration_stages.py src/align_genotype/run_workflow.py test/test_qc_calibration_stages.py
git commit -m "feat(qc_calibration): add the calibration stages to the workflow

A DatasetStage fans out one extract job per dataset, keyed on the Metamist
analysis ID so extraction re-runs only when the underlying report changed. A
forced MultiCohortStage does the cross-dataset analysis. Neither depends on a
production stage, and both are inert unless qc_calibration.enabled is set.

Stage decisions live in free functions because Stage.__init__ calls
get_workflow(), so a stage cannot be constructed outside a real run."
```

---

## Task 17: Delete the remaining superseded modules

The second half of the deletion begun in Task 6b. `spec.py`, `cache.py`, `tomlio.py` and
`manifest.py` survived that far because the un-rewritten `discovery.py` imported
`tomlio` and `manifest`; Task 12 removed those imports, so they can go now.

**Files:**
- Delete: `src/align_genotype/qc_calibration/{spec,cache,tomlio,manifest}.py`
- Delete: `test/test_qc_calibration_{spec,cache,tomlio,manifest}.py`
- Modify: `src/align_genotype/qc_calibration/__init__.py`

- [ ] **Step 1: Delete the modules and their tests**

```bash
git rm src/align_genotype/qc_calibration/{spec,cache,tomlio,manifest}.py
git rm test/test_qc_calibration_{spec,cache,tomlio,manifest}.py
```

- [ ] **Step 2: Rewrite the package docstring**

Replace the contents of `src/align_genotype/qc_calibration/__init__.py` with:

```python
"""Derive candidate QC thresholds from every dataset's latest MultiQC report.

Run as two CPG Flow stages - see `qc_calibration_stages.py` and README.md in this
package. Every module here is pure: no Hail Batch, no stage machinery, so the analysis is
testable without either.
"""
```

- [ ] **Step 3: Verify nothing still imports the deleted modules**

Run: `grep -rn 'qc_calibration import \(spec\|cache\|tomlio\|manifest\)\|qc_calibration\.\(spec\|cache\|tomlio\|manifest\)' src/ test/`

Expected: no output.

- [ ] **Step 4: Confirm the final module set**

Run: `ls src/align_genotype/qc_calibration/`

Expected exactly: `README.md`, `__init__.py`, `discovery.py`, `extract.py`, `relative.py`,
`render.py`, `settings.py`, `snippet.py`, `stats.py`, `summary.py`, `thresholds.py`.

Run: `ls test/test_qc_calibration_*.py`

Expected exactly twelve files: `discovery`, `extract`, `relative`, `render`, `scripts`,
`settings`, `snippet`, `stages`, `stats`, `summary`, `thresholds`, `values`.

- [ ] **Step 5: Run the whole suite**

Run: `uv run pytest -q`

Expected: PASS, across those twelve plus `test_check_multiqc.py` and
`test_sg_qc_report.py`.

- [ ] **Step 6: Lint and format**

```bash
uv run ruff check src/ test/
uv run ruff format --check src/ test/
```

Expected: both pass.

- [ ] **Step 7: Verify the package still installs and imports**

Run: `uv sync --reinstall-package align_genotype && uv run python -c "import align_genotype.qc_calibration_stages; print('ok')"`

Expected: `ok`

- [ ] **Step 8: Commit**

```bash
git add -A src/align_genotype/qc_calibration test/
git commit -m "refactor(qc_calibration): remove the spec, cache, tomlio and manifest modules

Superseded by settings, values and extract. Nothing reads or writes TOML any
more, and provenance lives in the report JSON rather than a manifest artifact."
```

---

## Task 18: README, and a dry run against a real multicohort

**Files:**
- Rewrite: `src/align_genotype/qc_calibration/README.md`

- [ ] **Step 1: Rewrite the README**

Replace the whole of `src/align_genotype/qc_calibration/README.md` with a document
covering, in this order:

1. **What this is.** Two stages that read every dataset's latest CramMultiQC report and
   produce candidate warn/fail thresholds plus dataset-relative (MAD) tier evidence, as
   an HTML dashboard in the analysis dataset's web bucket. Run it when onboarding a new
   capture kit or sequencing protocol, or when refreshing thresholds.
2. **How to run it.** The analysis-runner invocation, and the config it needs:
   `qc_calibration.enabled = true`, `workflow.sequencing_type`, the `input_cohorts` that
   determine which datasets are covered, and
   `only_stages = ['QcCalibrationDatasetMetrics', 'QcCalibrationReport']`. State that
   `enabled` belongs only in the calibration run's config — set in a shared config, every
   production run gains a job per dataset. State that one run covers one sequencing type,
   so calibrating both means two invocations.
3. **What it writes.** The three paths, and what each is for.
4. **How to read the report.** The judgement, which is the part no number supplies:
   - Aim for a healthy dataset flagging roughly 0% fail and single-digit % warn. `fail`
     means "do not analyse without a decision"; `warn` means "a human should look, and
     usually proceeds with a note".
   - Preserve the lab's intent for hard gates unless the data clearly contradicts it.
     Tighten `warn` freely — it costs a human a look. Moving `fail` changes what gets
     analysed at all.
   - The bar for overriding the lab is a specific measured number of good sequencing
     groups their line would have discarded.
   - Don't hard-fail on metrics that track ancestry, biology or chemistry — `error_rate`
     and `HET_SNP_SENSITIVITY` vary for reasons that are not sample quality. A `fail` on
     them fires on populations, not problems.
   - Candidates are percentile tails. Read them against the flag-rate tables and against
     what the lab asked for, not on their own.
5. **When to adopt a dataset-relative tier.** Two conditions, both required: the spread
   is genuinely dataset- or protocol-dependent (look at the per-dataset medians, not at
   intuition), and the flag set stays stable as the dataset grows. Explain that a relative
   tier ships only `direction`, `k` and `min_samples` — production recomputes the
   threshold every run — so the per-dataset thresholds in the report are illustrative.
   Explain growth churn as a forecast and merge churn as a stress test, and name the
   tension: a metric earns a relative tier because its level shifts between datasets, and
   merge punishes exactly that.
6. **Expect REJECT on both shipped tiers.** Over the full dataset set, exome
   `ZERO_CVG_TARGETS_PCT` measures ~1.0% growth / ~51.1% merge and genome
   `reads_duplicated_percent` ~8.6% / ~64.7%; both were adopted on a hand-picked subset.
   That is a live QC question for the team to settle, not a tool failure.
7. **Metric keys differ by sequencing type.** Genome uses Picard `CollectWgsMetrics`
   (`MEDIAN_COVERAGE`, `MEAN_COVERAGE`, `PCT_1X`…`PCT_100X`, `HET_SNP_SENSITIVITY`);
   exome uses `CollectHsMetrics` (`MEAN_TARGET_COVERAGE`, `PCT_TARGET_BASES_20X/50X`,
   `FOLD_80_BASE_PENALTY`, `ZERO_CVG_TARGETS_PCT`, `PCT_SELECTED_BASES`, `PCT_OFF_BAIT`,
   `AT_DROPOUT`/`GC_DROPOUT`). Shared via samtools: `reads_mapped_percent`,
   `reads_duplicated_percent`, `reads_properly_paired_percent`, `reads_MQ0_percent`,
   `error_rate`. Contamination via verifybamid: `FREEMIX`. Do not copy an exome metric
   list onto a genome run.
   One sentence on the trap, with no incident narrative: **`PCT_PF_READS_ALIGNED` is not
   in `report_general_stats_data`** — MultiQC writes it only to `report_saved_raw_data`,
   which the production check never reads, so it can never be gated. Use samtools
   `reads_mapped_percent`.
8. **Local iteration.** Both job scripts are `python -m` runnable, so downloaded values
   files can be re-reported with a different `k` or metric list without re-parsing any
   MultiQC report:
   `python -m align_genotype.scripts.qc_calibration_report --values a.json --values b.json --output-json out.json --output-html out.html`
9. **Troubleshooting.** A metric shows MISSING in the presence matrix (usually a MultiQC
   key rename; check the presence matrix against what the report actually carries). A
   dataset appears under "not included" (no completed CramMultiQC analysis at this
   sequencing type — check CramMultiQC has run for it). A metric is present but its value
   list is empty (the key exists and every value was unusable, typically Picard's `'?'`
   placeholder; the dropped count shows it) — a different problem from an absent key, and
   a different fix.

Do not reference `testing_scripts/`, "the old workflow", "the manual process this
replaces", or any local scratch script. This document stands alone.

- [ ] **Step 2: Verify no scratch references survive anywhere**

Run: `grep -rn 'testing_scripts\|qc_calibrate \|old workflow\|manual process\|fetch_multiqc_json_paths\|mad_relative_prototype\|verify_relative_6cohorts' src/`

Expected: no output.

- [ ] **Step 3: Check the terminology sweep is complete**

Run: `grep -rni 'cohort' src/align_genotype/qc_calibration/ src/align_genotype/qc_calibration_stages.py src/align_genotype/jobs/qc_calibration.py src/align_genotype/scripts/qc_calibration_*.py`

Expected: only hits referring to CPG Flow's own `MultiCohort`/`input_cohorts`/`CohortStage`
concepts. Any hit meaning "the set of samples in one MultiQC report" is a miss — fix it.

- [ ] **Step 4: Run everything**

```bash
uv run pytest -q
uv run ruff check src/ test/
uv run ruff format --check src/ test/
```

Expected: all pass.

- [ ] **Step 5: Commit**

```bash
git add src/align_genotype/qc_calibration/README.md
git commit -m "docs(qc_calibration): rewrite the README for the workflow

Stands alone: how to run the stages, what they write, and how to read the
report. Keeps the judgement prose - healthy-dataset flag rates, preserving the
lab's intent, not hard-failing on ancestry-driven metrics, the two conditions
for a relative tier - and drops every reference to the local scratch scripts."
```

- [ ] **Step 6: Dry run against a real multicohort**

Build a calibration run config containing `qc_calibration.enabled = true`, the
`input_cohorts` to cover, `workflow.sequencing_type = "genome"`, and
`only_stages = ['QcCalibrationDatasetMetrics', 'QcCalibrationReport']`. Then:

```bash
analysis-runner \
  --dataset <analysis-dataset> \
  --access-level test \
  --description "QC calibration dry run" \
  --output-dir qc_calibration_dryrun \
  --config <your-calibration-config>.toml \
  run_workflow --dry_run
```

Check in the printed DAG that:
- only the two calibration stages are queued, and no alignment or genotyping stage is;
- one extract job appears per dataset that has a CramMultiQC report;
- datasets without one are absent rather than erroring;
- the report job's expected outputs are under `qc_calibration/genome/`.

Then re-run without `--dry_run` at `--access-level test` and open the HTML from the web
bucket.

- [ ] **Step 7: Confirm production runs stay inert**

Run the same `analysis-runner` command with your ordinary production config (which does
not set `qc_calibration.enabled`) and `--dry_run`, and confirm neither calibration stage
queues a job.

---

## Notes for the implementer

- **Deletion is split across Tasks 6b and 17**, because the old modules import each
  other and `discovery.py` keeps `tomlio`/`manifest` alive until Task 12 rewrites it.
  Doing it in one pass at the end leaves the suite red from Task 7 onward. The full suite
  should be green at the end of every task.
- **Run order matters for Task 6.** `stats.py` is imported by `relative.py` and
  `summary.py`; do the wording pass before rewriting `relative.py` so the two commits
  stay separable.
- **`test/` is not a package.** `test_qc_calibration_render.py` imports fixtures from
  `test.test_qc_calibration_summary`. If that import fails, add an empty
  `test/__init__.py` in the same commit, or inline the two helpers.
- **Validate config values, never coerce them.** Every setting read from TOML goes
  through `settings._require_bool` / `_require_number` / `_require_int`. `bool('false')`
  is `True` and `int(True)` is `1`, and a bare `float()` raises a `ValueError` that never
  names the offending key. This applies to any new setting added in a later task.
- **The test lists in this plan are not exhaustive against their own code blocks.** Task 2
  shipped two branches (`parse_metric`'s unknown-key check, `CalibrationSettings.metric`'s
  `KeyError`) that the plan's own test list never exercised. When implementing, check each
  branch of the code you are given and add the missing test rather than transcribing the
  list as-is.
- **MAD fixtures need at least three points.** With two, the modified z-score is always
  exactly `±0.6745` (MAD equals half the range), so nothing can ever be flagged at any
  `k` this codebase uses, and a test built on two points passes whatever the code does.
  Task 1's original test had this defect. Check any small relative/churn fixture against
  it.
- **Ruff is strict here.** `ANN`, `BLE`, `TCH`, `PLC0415` (lazy imports) and `DTZ` are
  all enabled. Use `datetime.now(tz=timezone.utc)`, mark deliberate function-local
  imports `# noqa: PLC0415`, and annotate every signature.
- **Line length is 120**, single quotes, `uv run ruff format` to normalise.
