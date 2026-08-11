# QC Calibration Workflow Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Replace the five ad-hoc `testing_scripts/` calibration scripts with a committed, tested `qc_calibrate` CLI that derives warn/fail and cohort-relative QC thresholds from a set of MultiQC datasets and emits the `config_template.toml :: qc_thresholds.<seq_type>` block.

**Architecture:** A subpackage at `src/align_genotype/qc_calibration/` built around three durable artifacts — a **manifest** (which cohorts, from Metamist), a **value cache** (small per-metric numbers, produced by exactly one parse of each huge MultiQC JSON), and a **calibration spec** (metrics, directions, candidate thresholds). Eight subcommands read and write those artifacts. Section normalisation, value extraction and the MAD threshold all come from `check_multiqc` by import, so calibration and production enforcement can never disagree.

**Tech Stack:** Python 3.10–3.11, `click` (CLI), `numpy` (percentiles/masks), `tomllib`/`tomli` (reading TOML), `cpg_utils.to_path` (gs:// and local URIs), `cpg_utils.config.set_config_paths` (driving the production checker), `pytest`. No new declared dependencies.

**Design spec:** `docs/superpowers/specs/2026-08-11-qc-calibration-workflow-design.md`

**Branch:** `qc-calibration-workflow`, branched off `mad-relative-flagging`. PR #75 targets `mad-relative-flagging` (stacked under #74 → #73 → main).

---

## Conventions for every task

- **Style:** ruff, line-length 120, **single quotes**, type annotations on function signatures. Run `uvx ruff check <files>` and `uvx ruff format <files>` before each commit. The two pre-existing `PLR0917` warnings in `check_multiqc.py` are baseline — ignore them, don't "fix" them.
- **Tests:** `uv run --extra test pytest test/ -q`. Test files live in `test/`, flat, named `test_<module>.py`.
- **Never commit real CPG dataset names.** Use `dataset-a`, `dataset-b`, `cohort-1` in fixtures, docstrings and examples. Grep staged content before committing.
- **TDD is mandatory:** write the test, run it, watch it fail for the *right reason*, then implement.
- Do not create the `calibration/` working directory in git — it is added to `.gitignore` in Task 11.

---

## File structure

**Modified (production code — Task 0 only):**

| File | Change |
|---|---|
| `src/align_genotype/scripts/check_multiqc.py` | Add `normalise_sections()`; make `_gather_metric_values` public and return a drop count; call the normaliser in `run()` |
| `test/test_check_multiqc.py` | Tests for both section shapes and the new return signature |

**Created (the new subpackage):**

| File | Responsibility |
|---|---|
| `src/align_genotype/qc_calibration/__init__.py` | Package marker; re-export the error types |
| `src/align_genotype/qc_calibration/tomlio.py` | `tomllib` import shim (3.10/3.11) + scalar TOML writing helpers |
| `src/align_genotype/qc_calibration/spec.py` | `MetricSpec`/`RelativeSpec`/`CalibrationSpec`; load, validate, dump |
| `src/align_genotype/qc_calibration/manifest.py` | `Cohort`/`Manifest`; load, validate, dump |
| `src/align_genotype/qc_calibration/cache.py` | `CohortValues`/`ValueCache`; JSON load/save; staleness check |
| `src/align_genotype/qc_calibration/collect.py` | Parse one MultiQC JSON → survey row + extracted values |
| `src/align_genotype/qc_calibration/stats.py` | Percentiles, flag rates, cohort-growth churn |
| `src/align_genotype/qc_calibration/relative.py` | MAD evaluation driven through production `relative_flags()` |
| `src/align_genotype/qc_calibration/suggest.py` | Seed candidate thresholds from percentiles |
| `src/align_genotype/qc_calibration/emit.py` | Generate the `[qc_thresholds...]` TOML block + rationale comments |
| `src/align_genotype/qc_calibration/dryrun.py` | Temp config + real `check_multiqc.run()` |
| `src/align_genotype/qc_calibration/discovery.py` | Metamist query → manifest (lazy import; the only Metamist-aware module) |
| `src/align_genotype/qc_calibration/report.py` | Shared stdout table formatting |
| `src/align_genotype/qc_calibration/cli.py` | Click group and options; no logic |
| `src/align_genotype/qc_calibration/README.md` | Operator guide |

**Created (tests):** `test/test_qc_calibration_tomlio.py`, `_spec.py`, `_manifest.py`, `_cache.py`, `_collect.py`, `_stats.py`, `_relative.py`, `_suggest.py`, `_emit.py`, `_cli.py`, and `test/conftest.py` for shared fixtures.

> **Note vs the design spec:** the spec listed 12 modules; this plan adds `tomlio.py` (a 13th) so the 3.10 `tomli` shim and the scalar-writing helpers live in exactly one place instead of being duplicated across `spec.py` and `manifest.py`.

---

## Task 0: Shared section normaliser in `check_multiqc`

This is the enabling fix. `check_multiqc.run()` reads `d['report_general_stats_data']` and immediately calls `.items()` on it (`check_multiqc.py:280,285`), but MultiQC v1.14 stores that field as a positional **list**, not a dict — so a v1.14 report crashes production. The calibration tool must extract values through the same code path as production, so the normaliser has to live in `check_multiqc` and be imported, not re-implemented.

Also makes `_gather_metric_values` public (so the tool can import it without tripping ruff's `SLF001`) and has it report how many values were dropped as non-numeric — the calibration survey needs that count.

**Files:**
- Modify: `src/align_genotype/scripts/check_multiqc.py`
- Test: `test/test_check_multiqc.py`

- [ ] **Step 1: Write the failing tests**

Append to `test/test_check_multiqc.py`:

```python
# --- section shape normalisation ---------------------------------------------


def test_normalise_sections_dict_shape_passes_through():
    raw = {'picard': {'S1': {'MEDIAN_COVERAGE': 30}}, 'samtools': {'S1': {'error_rate': 0.01}}}
    assert check_multiqc.normalise_sections(raw) == raw


def test_normalise_sections_list_shape_gets_positional_names():
    raw = [{'S1': {'FREEMIX': 0.01}}, {'S1': {'MEDIAN_COVERAGE': 30}}]
    assert check_multiqc.normalise_sections(raw) == {
        'section_0': {'S1': {'FREEMIX': 0.01}},
        'section_1': {'S1': {'MEDIAN_COVERAGE': 30}},
    }


def test_normalise_sections_drops_non_dict_members():
    assert check_multiqc.normalise_sections([{'S1': {'a': 1}}, None, 'junk']) == {'section_0': {'S1': {'a': 1}}}
    assert check_multiqc.normalise_sections({'picard': {'S1': {'a': 1}}, 'broken': None}) == {
        'picard': {'S1': {'a': 1}},
    }


def test_normalise_sections_unexpected_type_is_empty():
    assert check_multiqc.normalise_sections(None) == {}
    assert check_multiqc.normalise_sections('nonsense') == {}


def test_run_handles_list_shaped_general_stats(tmp_path, patch_config):
    """A MultiQC v1.14 report stores general stats as a list; it must not crash."""
    patch_config('genome', GENOME_THRESHOLDS)
    path = _write_json(tmp_path, [{'CPG1|S1': {'MEDIAN_COVERAGE': 5}}])
    result = _run(path, tmp_path / 'out.json')
    assert _flags_by_metric(result, 'CPG1')['MEDIAN_COVERAGE']['severity'] == 'fail'


def test_run_raises_when_general_stats_absent(tmp_path, patch_config):
    """A report with no general stats must fail loudly, not silently check nothing."""
    patch_config('genome', GENOME_THRESHOLDS)
    path = tmp_path / 'multiqc_data.json'
    path.write_text(json.dumps({'report_saved_raw_data': {}}))
    with pytest.raises(ValueError, match='report_general_stats_data'):
        _run(str(path), tmp_path / 'out.json')


# --- gather_metric_values ------------------------------------------------------


def test_gather_metric_values_returns_entries_and_drop_count():
    sections = {
        'picard': {'S1': {'MEDIAN_COVERAGE': 30}, 'S2': {'MEDIAN_COVERAGE': '?'}},
        'samtools': {'S1': {'MEDIAN_COVERAGE': '28.5'}, 'S3': {'other': 1}},
    }
    entries, n_dropped = check_multiqc.gather_metric_values(sections, 'MEDIAN_COVERAGE')
    assert sorted(entries) == [('picard', 'S1', 30.0), ('samtools', 'S1', 28.5)]
    assert n_dropped == 1
```

- [ ] **Step 2: Run the tests to verify they fail**

Run: `uv run --extra test pytest test/test_check_multiqc.py -q`

Expected: 7 failures. `AttributeError: module 'align_genotype.scripts.check_multiqc' has no attribute 'normalise_sections'` for the normaliser tests, `has no attribute 'gather_metric_values'` for the last one, and `AttributeError: 'list' object has no attribute 'items'` for `test_run_handles_list_shaped_general_stats` — that last one is the production bug reproducing.

- [ ] **Step 3: Add `normalise_sections`**

In `src/align_genotype/scripts/check_multiqc.py`, insert immediately after the `MODIFIED_Z_CONST` constant (around line 88):

```python
def normalise_sections(raw: Any) -> dict[str, dict[str, Any]]:
    """Normalise ``report_general_stats_data`` to ``{section: {sample: {metric: value}}}``.

    MultiQC >=1.33 keys general-stats sections by tool name (a dict); v1.14 stores a
    positional list instead, and both shapes turn up in real archived reports. Every
    consumer - this check and the qc_calibration tooling - goes through here, so
    calibration and enforcement can never disagree about what a report contains.
    Members that aren't dicts are dropped.
    """
    if isinstance(raw, dict):
        return {str(name): section for name, section in raw.items() if isinstance(section, dict)}
    if isinstance(raw, list):
        return {f'section_{i}': section for i, section in enumerate(raw) if isinstance(section, dict)}
    return {}
```

- [ ] **Step 4: Make `gather_metric_values` public and return a drop count**

Replace the whole `_gather_metric_values` function (lines 159–172) with:

```python
def gather_metric_values(sections: dict[str, Any], metric: str) -> tuple[list[tuple[str, str, float]], int]:
    """``([(section, sample, value)], n_dropped)`` for every sample carrying `metric`.

    Non-numeric placeholders (Picard writes ``'?'`` when coverage is ~0) are skipped and
    counted rather than raising, so one bad cell can't sink the whole check.
    """
    entries: list[tuple[str, str, float]] = []
    n_dropped = 0
    for section_name, section in sections.items():
        for sample, val_by_metric in section.items():
            if metric not in val_by_metric:
                continue
            try:
                entries.append((section_name, sample, float(val_by_metric[metric])))
            except (TypeError, ValueError):
                n_dropped += 1
                logging.warning(f'{sample}: metric {metric!r} non-numeric {val_by_metric[metric]!r}; skipping.')
    return entries, n_dropped
```

- [ ] **Step 5: Update the one caller**

In `relative_flags` (line 237), change:

```python
        entries = _gather_metric_values(sections, metric)
```

to:

```python
        entries, _ = gather_metric_values(sections, metric)
```

- [ ] **Step 6: Call the normaliser in `run()`**

Replace lines 278–280:

```python
    with to_path(multiqc_json_path).open() as f:
        d = json.load(f)
        sections = d['report_general_stats_data']
```

with:

```python
    with to_path(multiqc_json_path).open() as f:
        d = json.load(f)
    sections = normalise_sections(d.get('report_general_stats_data'))
    if not sections:
        raise ValueError(
            f'No usable report_general_stats_data in {multiqc_json_path}; refusing to report a clean QC check '
            f'on a report we could not read.',
        )
```

- [ ] **Step 7: Run the tests to verify they pass**

Run: `uv run --extra test pytest test/test_check_multiqc.py -q`

Expected: all pass, including the pre-existing tests (the normaliser is a strict widening — dict input behaves identically).

- [ ] **Step 8: Lint and commit**

```bash
uvx ruff format src/align_genotype/scripts/check_multiqc.py test/test_check_multiqc.py
uvx ruff check src/align_genotype/scripts/check_multiqc.py test/test_check_multiqc.py
uv run --extra test pytest test/ -q
git add src/align_genotype/scripts/check_multiqc.py test/test_check_multiqc.py
git commit -m "fix: normalise MultiQC v1.14 list-shaped general stats in check_multiqc

report_general_stats_data is a dict keyed by section name in MultiQC >=1.33
but a positional list in v1.14. run() called .items() on it directly, so a
v1.14 report raised AttributeError and the QC check never ran.

Adds normalise_sections() as the single shared normaliser, and raises rather
than silently reporting a clean check when a report has no readable general
stats. Also makes gather_metric_values public and returns its non-numeric
drop count, so the calibration tooling can import the same extraction path
instead of re-implementing it."
```

---

## Task 1: `tomlio` — TOML read shim and scalar writers

One place for the 3.10/3.11 `tomllib` difference and for turning Python scalars into TOML. The manifest and spec writers both need it; the artifacts are flat tables of scalars, so a full TOML-writing dependency isn't warranted.

**Files:**
- Create: `src/align_genotype/qc_calibration/__init__.py`
- Create: `src/align_genotype/qc_calibration/tomlio.py`
- Test: `test/test_qc_calibration_tomlio.py`

- [ ] **Step 1: Write the failing test**

Create `test/test_qc_calibration_tomlio.py`:

```python
"""Unit tests for the calibration TOML helpers."""

import pytest

from align_genotype.qc_calibration import tomlio


@pytest.mark.parametrize(
    ('value', 'expected'),
    [
        (True, 'true'),
        (False, 'false'),
        (15, '15'),
        (3.5, '3.5'),
        (0.75, '0.75'),
        ('genome', '"genome"'),
    ],
)
def test_fmt_value(value, expected):
    assert tomlio.fmt_value(value) == expected


def test_fmt_value_bool_wins_over_int():
    """bool is a subclass of int; it must not render as 1/0."""
    assert tomlio.fmt_value(True) == 'true'


def test_fmt_value_escapes_quotes_and_backslashes():
    assert tomlio.fmt_value('a "quoted" c:\\path') == '"a \\"quoted\\" c:\\\\path"'


def test_fmt_kv():
    assert tomlio.fmt_kv('k', 3.5) == 'k = 3.5'


def test_fmt_kv_quoted_key():
    assert tomlio.fmt_kv('MEDIAN_COVERAGE', 15, quote_key=True) == '"MEDIAN_COVERAGE" = 15'


@pytest.mark.parametrize('key', ['MEDIAN_COVERAGE', 'reads_mapped_percent', 'dataset-a', 'PCT_20X'])
def test_require_bare_key_accepts_identifier_like_keys(key):
    tomlio.require_bare_key(key, 'metric')  # does not raise


@pytest.mark.parametrize('key', ['odd.key', 'has space', 'quote"d', '', 'sl/ash'])
def test_require_bare_key_rejects_anything_needing_quoting(key):
    with pytest.raises(ValueError, match='not a bare TOML key'):
        tomlio.require_bare_key(key, 'metric')


def test_loads_round_trips_written_scalars():
    text = '\n'.join([
        '[table]',
        tomlio.fmt_kv('name', 'dataset-a'),
        tomlio.fmt_kv('k', 3.5),
        tomlio.fmt_kv('gated', True),
    ])
    assert tomlio.loads(text) == {'table': {'name': 'dataset-a', 'k': 3.5, 'gated': True}}
```

- [ ] **Step 2: Run the test to verify it fails**

Run: `uv run --extra test pytest test/test_qc_calibration_tomlio.py -q`

Expected: FAIL — `ModuleNotFoundError: No module named 'align_genotype.qc_calibration'`

- [ ] **Step 3: Create the package marker**

Create `src/align_genotype/qc_calibration/__init__.py`:

```python
"""Derive warn/fail and cohort-relative QC thresholds from a set of MultiQC datasets.

See README.md in this package for the operator workflow. The entry point is the
``qc_calibrate`` console script (``qc_calibration.cli``).
"""
```

- [ ] **Step 4: Write the implementation**

Create `src/align_genotype/qc_calibration/tomlio.py`:

```python
"""TOML reading and scalar writing for the calibration artifacts.

Reading uses ``tomllib`` (stdlib from 3.11, ``tomli`` on 3.10 - already in the resolved
dependency graph under a version marker, so nothing new is declared). Writing is
hand-rolled: manifests and specs are flat tables of scalars, which doesn't justify a
style-preserving TOML writer dependency.
"""

import re
import sys
from pathlib import Path
from typing import Any

if sys.version_info >= (3, 11):
    import tomllib
else:  # pragma: no cover - 3.10 only
    import tomli as tomllib

# Keys we can write unquoted. Metric keys come from Picard/samtools and cohort labels
# from Metamist dataset names, so both are already identifier-like; anything else is
# rejected at load time rather than silently mangled when we write the file back out.
BARE_KEY_RE = re.compile(r'^[A-Za-z0-9_-]+$')


def require_bare_key(key: str, kind: str) -> None:
    """Raise unless `key` can be written as an unquoted TOML key."""
    if not BARE_KEY_RE.match(key):
        raise ValueError(f'{kind} {key!r} is not a bare TOML key; expected only letters, digits, underscore or hyphen')


def loads(text: str) -> dict[str, Any]:
    """Parse a TOML document from a string."""
    return tomllib.loads(text)


def load_path(path: str | Path) -> dict[str, Any]:
    """Parse a TOML document from a local or cloud path."""
    from cpg_utils import to_path

    with to_path(path).open('rb') as f:
        return tomllib.load(f)


def fmt_value(value: Any) -> str:
    """Render a Python scalar as a TOML value."""
    if isinstance(value, bool):  # bool subclasses int - must be checked first
        return 'true' if value else 'false'
    if isinstance(value, (int, float)):
        return repr(value)
    escaped = str(value).replace('\\', '\\\\').replace('"', '\\"')
    return f'"{escaped}"'


def fmt_kv(key: str, value: Any, quote_key: bool = False) -> str:
    """Render a ``key = value`` line. `quote_key` matches config_template.toml's style."""
    rendered_key = f'"{key}"' if quote_key else key
    return f'{rendered_key} = {fmt_value(value)}'
```

Note the `from cpg_utils import to_path` is inside `load_path` deliberately: it keeps `tomlio` importable in a bare interpreter, which makes the tests fast and independent of cpg-utils config state.

- [ ] **Step 5: Run the test to verify it passes**

Run: `uv run --extra test pytest test/test_qc_calibration_tomlio.py -q`

Expected: PASS (10 tests)

- [ ] **Step 6: Lint and commit**

```bash
uvx ruff format src/align_genotype/qc_calibration/ test/test_qc_calibration_tomlio.py
uvx ruff check src/align_genotype/qc_calibration/ test/test_qc_calibration_tomlio.py
git add src/align_genotype/qc_calibration/__init__.py src/align_genotype/qc_calibration/tomlio.py test/test_qc_calibration_tomlio.py
git commit -m "feat(qc_calibration): TOML read shim and scalar writers"
```

---

## Task 2: `spec` — the calibration spec

The spec is the single source of truth for a calibration run: which metrics are gated, in which direction, at what candidate values, and which take a cohort-relative warn tier. It replaces editing constants inside analysis scripts, which is what made the exome and genome calibrations irreproducible.

Four validation rules encode threshold-design decisions structurally rather than in prose:
1. A gated metric must define at least one of `fail` / `warn` / `relative` — otherwise it gates nothing.
2. `relative` requires `fail`: a cohort-relative warn tier always keeps an absolute hard gate behind it.
3. `relative` forbids `warn`: the warn tier *is* the relative one; an absolute warn alongside it would double-flag.
4. A non-gated metric carries no thresholds — it is surveyed and profiled for the record, deliberately not enforced.

**Files:**
- Create: `src/align_genotype/qc_calibration/spec.py`
- Test: `test/test_qc_calibration_spec.py`

- [ ] **Step 1: Write the failing test**

Create `test/test_qc_calibration_spec.py`:

```python
"""Unit tests for the calibration spec: load, validate, dump, round-trip."""

import pytest

from align_genotype.qc_calibration import spec as spec_mod
from align_genotype.qc_calibration.spec import CalibrationSpec, MetricSpec, RelativeSpec, SpecError

VALID = '''
seq_type = "genome"
cache = "calibration/genome_values.json"

[metrics.MEDIAN_COVERAGE]
direction = "min"
unit = "x"
fail = 15
warn = 25
reviewed = true
rationale = "Primary depth gate."

[metrics.reads_duplicated_percent]
direction = "max"
unit = "%"
fail = 40
reviewed = true
[metrics.reads_duplicated_percent.relative]
k = 3.5
min_cohort = 50

[metrics.error_rate]
direction = "max"
unit = "frac"
gated = false
rationale = "Varies with ancestry; not a defensible hard gate."
'''


def test_loads_valid_spec():
    spec = spec_mod.loads(VALID)
    assert spec.seq_type == 'genome'
    assert spec.cache == 'calibration/genome_values.json'
    assert [m.key for m in spec.metrics] == ['MEDIAN_COVERAGE', 'reads_duplicated_percent', 'error_rate']


def test_metric_defaults():
    depth = spec_mod.loads(VALID).metric('MEDIAN_COVERAGE')
    assert depth == MetricSpec(
        key='MEDIAN_COVERAGE',
        direction='min',
        unit='x',
        gated=True,
        fail=15,
        warn=25,
        relative=None,
        reviewed=True,
        rationale='Primary depth gate.',
    )


def test_relative_block_parsed():
    dup = spec_mod.loads(VALID).metric('reads_duplicated_percent')
    assert dup.relative == RelativeSpec(k=3.5, min_cohort=50)
    assert dup.warn is None


def test_gated_property_excludes_ungated():
    assert [m.key for m in spec_mod.loads(VALID).gated] == ['MEDIAN_COVERAGE', 'reads_duplicated_percent']


def test_metric_keys_property():
    assert spec_mod.loads(VALID).metric_keys == ('MEDIAN_COVERAGE', 'reads_duplicated_percent', 'error_rate')


def test_unknown_metric_lookup_raises():
    with pytest.raises(KeyError, match='NOPE'):
        spec_mod.loads(VALID).metric('NOPE')


# --- validation ---------------------------------------------------------------


@pytest.mark.parametrize(
    ('body', 'match'),
    [
        ('[metrics.M]\ndirection = "sideways"\nfail = 1\n', "direction must be 'min' or 'max'"),
        ('[metrics.M]\nfail = 1\n', 'missing required key: direction'),
        ('[metrics.M]\ndirection = "min"\nunit = "furlongs"\nfail = 1\n', 'unit must be one of'),
        ('[metrics.M]\ndirection = "min"\n', 'gated metric defines no fail, warn or relative tier'),
        ('[metrics.M]\ndirection = "max"\nwarn = 1\n[metrics.M.relative]\nk = 3.5\n', 'relative tier requires an absolute fail'),
        ('[metrics.M]\ndirection = "max"\ngated = false\nfail = 1\n', 'non-gated metric must not define'),
    ],
)
def test_validation_errors(body, match):
    text = f'seq_type = "genome"\ncache = "c.json"\n\n{body}'
    with pytest.raises(SpecError, match=match):
        spec_mod.loads(text)


def test_relative_with_absolute_warn_rejected():
    text = (
        'seq_type = "genome"\ncache = "c.json"\n\n'
        '[metrics.M]\ndirection = "max"\nfail = 40\nwarn = 30\n'
        '[metrics.M.relative]\nk = 3.5\n'
    )
    with pytest.raises(SpecError, match='relative tier replaces the absolute warn'):
        spec_mod.loads(text)


def test_missing_seq_type_rejected():
    with pytest.raises(SpecError, match='missing required key: seq_type'):
        spec_mod.loads('cache = "c.json"\n[metrics.M]\ndirection = "min"\nfail = 1\n')


def test_no_metrics_rejected():
    with pytest.raises(SpecError, match='defines no metrics'):
        spec_mod.loads('seq_type = "genome"\ncache = "c.json"\n')


def test_metric_key_needing_quoting_rejected():
    """dumps() writes metric keys unquoted, so a key that needs quoting can't round-trip."""
    text = 'seq_type = "genome"\ncache = "c.json"\n\n[metrics."odd.key"]\ndirection = "min"\nfail = 1\n'
    with pytest.raises(SpecError, match='not a bare TOML key'):
        spec_mod.loads(text)


# --- dump / round-trip ----------------------------------------------------------


def test_dump_round_trips():
    original = spec_mod.loads(VALID)
    assert spec_mod.loads(spec_mod.dumps(original)) == original


def test_dump_places_relative_subtable_after_its_parent_keys():
    text = spec_mod.dumps(spec_mod.loads(VALID))
    parent = text.index('[metrics.reads_duplicated_percent]')
    sub = text.index('[metrics.reads_duplicated_percent.relative]')
    nxt = text.index('[metrics.error_rate]')
    assert parent < sub < nxt


def test_with_thresholds_returns_a_new_spec():
    """Specs are immutable; suggest() builds an updated copy."""
    original = spec_mod.loads(VALID)
    updated = original.with_metric(
        replace_key='MEDIAN_COVERAGE',
        fail=20,
        warn=30,
        reviewed=False,
        rationale='seeded',
    )
    assert original.metric('MEDIAN_COVERAGE').fail == 15
    assert updated.metric('MEDIAN_COVERAGE').fail == 20
    assert updated.metric('MEDIAN_COVERAGE').reviewed is False
    assert [m.key for m in updated.metrics] == [m.key for m in original.metrics]


def test_save_and_load_path_round_trip(tmp_path):
    path = tmp_path / 'spec.toml'
    original = spec_mod.loads(VALID)
    spec_mod.save(original, path)
    assert spec_mod.load(path) == original
```

- [ ] **Step 2: Run the test to verify it fails**

Run: `uv run --extra test pytest test/test_qc_calibration_spec.py -q`

Expected: FAIL — `ModuleNotFoundError: No module named 'align_genotype.qc_calibration.spec'`

- [ ] **Step 3: Write the implementation**

Create `src/align_genotype/qc_calibration/spec.py`:

```python
"""The calibration spec - what to gate, in which direction, at what value.

A spec plus a manifest plus a value cache is everything needed to reproduce a
calibration, which is why thresholds live here rather than in editable Python
constants (the pain point in the original testing_scripts/ workflow).

Validation encodes the threshold-design rules structurally: a cohort-relative warn
tier always keeps an absolute fail gate behind it, and never coexists with an
absolute warn (the relative tier *is* the warn tier).
"""

from dataclasses import dataclass, replace
from pathlib import Path
from typing import Any

from align_genotype.qc_calibration import tomlio

DIRECTIONS = ('min', 'max')
UNITS = ('x', 'frac', '%')

_METRIC_KEYS = {'direction', 'unit', 'gated', 'fail', 'warn', 'relative', 'reviewed', 'rationale'}
_RELATIVE_KEYS = {'k', 'min_cohort'}


class SpecError(ValueError):
    """The calibration spec is missing something or is internally inconsistent."""


@dataclass(frozen=True)
class RelativeSpec:
    """Cohort-relative (MAD) warn tier settings; mirrors [qc_thresholds.*.relative.*]."""

    k: float = 3.5
    min_cohort: int = 50


@dataclass(frozen=True)
class MetricSpec:
    """One metric's direction, candidate thresholds and review state.

    `direction`: 'min' = higher is better (flag below the threshold); 'max' = lower is
    better (flag above it). `gated` metrics are enforced and are fatal if a cohort is
    missing them; un-gated metrics are surveyed and profiled for the record only.
    `reviewed` is False on anything `suggest` seeded and not yet signed off; emit-config
    refuses to emit an unreviewed gated metric.
    """

    key: str
    direction: str
    unit: str = 'frac'
    gated: bool = True
    fail: float | None = None
    warn: float | None = None
    relative: RelativeSpec | None = None
    reviewed: bool = False
    rationale: str = ''


@dataclass(frozen=True)
class CalibrationSpec:
    seq_type: str
    cache: str
    metrics: tuple[MetricSpec, ...]

    @property
    def gated(self) -> tuple[MetricSpec, ...]:
        return tuple(m for m in self.metrics if m.gated)

    @property
    def metric_keys(self) -> tuple[str, ...]:
        return tuple(m.key for m in self.metrics)

    def metric(self, key: str) -> MetricSpec:
        for m in self.metrics:
            if m.key == key:
                return m
        raise KeyError(f'{key!r} is not in this calibration spec')

    def with_metric(self, replace_key: str, **changes: Any) -> 'CalibrationSpec':
        """Return a copy with one metric's fields replaced; order preserved."""
        updated = tuple(replace(m, **changes) if m.key == replace_key else m for m in self.metrics)
        return replace(self, metrics=updated)


def _parse_relative(key: str, raw: Any) -> RelativeSpec:
    if not isinstance(raw, dict):
        raise SpecError(f'metric {key!r}: [metrics.{key}.relative] must be a table')
    if unknown := sorted(set(raw) - _RELATIVE_KEYS):
        raise SpecError(f'metric {key!r}: unknown relative key(s) {unknown}; expected {sorted(_RELATIVE_KEYS)}')
    return RelativeSpec(k=float(raw.get('k', 3.5)), min_cohort=int(raw.get('min_cohort', 50)))


def _parse_metric(key: str, raw: Any) -> MetricSpec:
    try:
        tomlio.require_bare_key(key, 'metric key')
    except ValueError as exc:
        raise SpecError(str(exc)) from exc
    if not isinstance(raw, dict):
        raise SpecError(f'metric {key!r}: [metrics.{key}] must be a table')
    if unknown := sorted(set(raw) - _METRIC_KEYS):
        raise SpecError(f'metric {key!r}: unknown key(s) {unknown}; expected {sorted(_METRIC_KEYS)}')
    if 'direction' not in raw:
        raise SpecError(f'metric {key!r}: missing required key: direction')
    direction = raw['direction']
    if direction not in DIRECTIONS:
        raise SpecError(f'metric {key!r}: direction must be \'min\' or \'max\', got {direction!r}')
    unit = raw.get('unit', 'frac')
    if unit not in UNITS:
        raise SpecError(f'metric {key!r}: unit must be one of {list(UNITS)}, got {unit!r}')

    metric = MetricSpec(
        key=key,
        direction=direction,
        unit=unit,
        gated=bool(raw.get('gated', True)),
        fail=raw.get('fail'),
        warn=raw.get('warn'),
        relative=_parse_relative(key, raw['relative']) if 'relative' in raw else None,
        reviewed=bool(raw.get('reviewed', False)),
        rationale=str(raw.get('rationale', '')),
    )
    _validate_metric(metric)
    return metric


def _validate_metric(metric: MetricSpec) -> None:
    key = metric.key
    has_tier = metric.fail is not None or metric.warn is not None or metric.relative is not None
    if not metric.gated:
        if has_tier:
            raise SpecError(
                f'metric {key!r}: non-gated metric must not define fail, warn or relative - '
                f'set gated = true to enforce it, or remove the thresholds.',
            )
        return
    if not has_tier:
        raise SpecError(f'metric {key!r}: gated metric defines no fail, warn or relative tier')
    if metric.relative is not None:
        if metric.fail is None:
            raise SpecError(
                f'metric {key!r}: a relative tier requires an absolute fail gate behind it - '
                f'cohort-relative flagging is warn-only.',
            )
        if metric.warn is not None:
            raise SpecError(
                f'metric {key!r}: the relative tier replaces the absolute warn; remove `warn` '
                f'or remove the [metrics.{key}.relative] block.',
            )


def loads(text: str) -> CalibrationSpec:
    """Parse a calibration spec from a TOML string."""
    return _from_dict(tomlio.loads(text))


def load(path: str | Path) -> CalibrationSpec:
    """Load a calibration spec from a local or cloud path."""
    return _from_dict(tomlio.load_path(path))


def _from_dict(raw: dict[str, Any]) -> CalibrationSpec:
    for required in ('seq_type', 'cache'):
        if required not in raw:
            raise SpecError(f'calibration spec: missing required key: {required}')
    metrics_raw = raw.get('metrics', {})
    if not metrics_raw:
        raise SpecError('calibration spec defines no metrics; add at least one [metrics.<KEY>] table')
    return CalibrationSpec(
        seq_type=str(raw['seq_type']),
        cache=str(raw['cache']),
        metrics=tuple(_parse_metric(key, value) for key, value in metrics_raw.items()),
    )


def dumps(spec: CalibrationSpec) -> str:
    """Render a spec back to TOML, preserving metric order."""
    lines = [tomlio.fmt_kv('seq_type', spec.seq_type), tomlio.fmt_kv('cache', spec.cache)]
    for metric in spec.metrics:
        lines += ['', f'[metrics.{metric.key}]', tomlio.fmt_kv('direction', metric.direction)]
        lines.append(tomlio.fmt_kv('unit', metric.unit))
        lines.append(tomlio.fmt_kv('gated', metric.gated))
        if metric.fail is not None:
            lines.append(tomlio.fmt_kv('fail', metric.fail))
        if metric.warn is not None:
            lines.append(tomlio.fmt_kv('warn', metric.warn))
        lines.append(tomlio.fmt_kv('reviewed', metric.reviewed))
        lines.append(tomlio.fmt_kv('rationale', metric.rationale))
        # The relative sub-table must follow all of its parent's scalar keys.
        if metric.relative is not None:
            lines += [
                f'[metrics.{metric.key}.relative]',
                tomlio.fmt_kv('k', metric.relative.k),
                tomlio.fmt_kv('min_cohort', metric.relative.min_cohort),
            ]
    return '\n'.join(lines) + '\n'


def save(spec: CalibrationSpec, path: str | Path) -> None:
    """Write a spec to disk."""
    Path(path).write_text(dumps(spec))
```

- [ ] **Step 4: Run the test to verify it passes**

Run: `uv run --extra test pytest test/test_qc_calibration_spec.py -q`

Expected: PASS (about 20 tests including the parametrised cases)

- [ ] **Step 5: Lint and commit**

```bash
uvx ruff format src/align_genotype/qc_calibration/spec.py test/test_qc_calibration_spec.py
uvx ruff check src/align_genotype/qc_calibration/spec.py test/test_qc_calibration_spec.py
git add src/align_genotype/qc_calibration/spec.py test/test_qc_calibration_spec.py
git commit -m "feat(qc_calibration): calibration spec with structural threshold validation

A relative tier requires an absolute fail behind it and forbids an absolute
warn, so the design rules for cohort-relative flagging are enforced by the
loader rather than documented in prose."
```

---

## Task 3: `manifest` — which cohorts a calibration used

The manifest is the reproducibility record: it names every cohort and the exact MultiQC JSON URI it came from. It is generated by `discover` (Task 11) but is deliberately hand-editable — dropping a bad cohort, pinning an older analysis, or pointing at a local file for testing are all normal operations.

`uri` is resolved with `cpg_utils.to_path`, which already handles `gs://` and local paths, so GCS and local files work with no extra code. Cohort labels come from the retrieval layer, never from parsing filenames.

**Files:**
- Create: `src/align_genotype/qc_calibration/manifest.py`
- Test: `test/test_qc_calibration_manifest.py`

- [ ] **Step 1: Write the failing test**

Create `test/test_qc_calibration_manifest.py`:

```python
"""Unit tests for the cohort manifest."""

import pytest

from align_genotype.qc_calibration import manifest as manifest_mod
from align_genotype.qc_calibration.manifest import Cohort, ManifestError

VALID = '''
seq_type = "genome"
generated = "2026-08-11T13:40:00"

[cohorts.dataset-a]
uri = "gs://cpg-dataset-a-main/qc/multiqc_data.json"
analysis_id = 84213
timestamp = "2026-06-02T04:11:09"

[cohorts.dataset-b]
uri = "file:///tmp/dataset-b_multiqc_data.json"
'''


def test_loads_valid_manifest():
    m = manifest_mod.loads(VALID)
    assert m.seq_type == 'genome'
    assert m.generated == '2026-08-11T13:40:00'
    assert m.labels == ('dataset-a', 'dataset-b')


def test_cohort_fields_and_defaults():
    m = manifest_mod.loads(VALID)
    assert m.cohort('dataset-a') == Cohort(
        label='dataset-a',
        uri='gs://cpg-dataset-a-main/qc/multiqc_data.json',
        analysis_id=84213,
        timestamp='2026-06-02T04:11:09',
    )
    assert m.cohort('dataset-b') == Cohort(label='dataset-b', uri='file:///tmp/dataset-b_multiqc_data.json')


def test_unknown_cohort_raises():
    with pytest.raises(KeyError, match='dataset-z'):
        manifest_mod.loads(VALID).cohort('dataset-z')


@pytest.mark.parametrize(
    ('text', 'match'),
    [
        ('generated = "x"\n[cohorts.a]\nuri = "u"\n', 'missing required key: seq_type'),
        ('seq_type = "genome"\ngenerated = "x"\n', 'lists no cohorts'),
        ('seq_type = "genome"\ngenerated = "x"\n[cohorts.a]\nanalysis_id = 1\n', "cohort 'a': missing required key: uri"),
        ('seq_type = "genome"\ngenerated = "x"\n[cohorts."a.b"]\nuri = "u"\n', 'not a bare TOML key'),
        ('seq_type = "genome"\ngenerated = "x"\n[cohorts.a]\nuri = "u"\nnope = 1\n', "unknown key(s) \\['nope'\\]"),
    ],
)
def test_validation_errors(text, match):
    with pytest.raises(ManifestError, match=match):
        manifest_mod.loads(text)


def test_dump_round_trips():
    original = manifest_mod.loads(VALID)
    assert manifest_mod.loads(manifest_mod.dumps(original)) == original


def test_dump_omits_absent_optional_fields():
    text = manifest_mod.dumps(manifest_mod.loads(VALID))
    dataset_b_block = text[text.index('[cohorts.dataset-b]') :]
    assert 'analysis_id' not in dataset_b_block
    assert 'timestamp' not in dataset_b_block


def test_save_and_load_path_round_trip(tmp_path):
    path = tmp_path / 'manifest.toml'
    original = manifest_mod.loads(VALID)
    manifest_mod.save(original, path)
    assert manifest_mod.load(path) == original
```

- [ ] **Step 2: Run the test to verify it fails**

Run: `uv run --extra test pytest test/test_qc_calibration_manifest.py -q`

Expected: FAIL — `ModuleNotFoundError: No module named 'align_genotype.qc_calibration.manifest'`

- [ ] **Step 3: Write the implementation**

Create `src/align_genotype/qc_calibration/manifest.py`:

```python
"""The cohort manifest - which datasets a calibration was derived from.

Generated by ``qc_calibrate discover`` but intended to be hand-edited: dropping a
cohort, pinning an older analysis, or substituting a local file are all normal. Because
it is a committed-alongside artifact rather than a live query, a calibration stays
reproducible after Metamist moves on.
"""

from dataclasses import dataclass
from pathlib import Path
from typing import Any

from align_genotype.qc_calibration import tomlio

_COHORT_KEYS = {'uri', 'analysis_id', 'timestamp'}


class ManifestError(ValueError):
    """The manifest is missing something or is malformed."""


@dataclass(frozen=True)
class Cohort:
    """One dataset's MultiQC report. `uri` is anything cpg_utils.to_path accepts."""

    label: str
    uri: str
    analysis_id: int | None = None
    timestamp: str | None = None


@dataclass(frozen=True)
class Manifest:
    seq_type: str
    generated: str
    cohorts: tuple[Cohort, ...]

    @property
    def labels(self) -> tuple[str, ...]:
        return tuple(c.label for c in self.cohorts)

    def cohort(self, label: str) -> Cohort:
        for c in self.cohorts:
            if c.label == label:
                return c
        raise KeyError(f'{label!r} is not in this manifest; known cohorts: {list(self.labels)}')


def _parse_cohort(label: str, raw: Any) -> Cohort:
    try:
        tomlio.require_bare_key(label, 'cohort label')
    except ValueError as exc:
        raise ManifestError(str(exc)) from exc
    if not isinstance(raw, dict):
        raise ManifestError(f'cohort {label!r}: [cohorts.{label}] must be a table')
    if unknown := sorted(set(raw) - _COHORT_KEYS):
        raise ManifestError(f'cohort {label!r}: unknown key(s) {unknown}; expected {sorted(_COHORT_KEYS)}')
    if 'uri' not in raw:
        raise ManifestError(f'cohort {label!r}: missing required key: uri')
    analysis_id = raw.get('analysis_id')
    timestamp = raw.get('timestamp')
    return Cohort(
        label=label,
        uri=str(raw['uri']),
        analysis_id=int(analysis_id) if analysis_id is not None else None,
        timestamp=str(timestamp) if timestamp is not None else None,
    )


def loads(text: str) -> Manifest:
    """Parse a manifest from a TOML string."""
    return _from_dict(tomlio.loads(text))


def load(path: str | Path) -> Manifest:
    """Load a manifest from a local or cloud path."""
    return _from_dict(tomlio.load_path(path))


def _from_dict(raw: dict[str, Any]) -> Manifest:
    if 'seq_type' not in raw:
        raise ManifestError('manifest: missing required key: seq_type')
    cohorts_raw = raw.get('cohorts', {})
    if not cohorts_raw:
        raise ManifestError('manifest lists no cohorts; add at least one [cohorts.<LABEL>] table')
    return Manifest(
        seq_type=str(raw['seq_type']),
        generated=str(raw.get('generated', '')),
        cohorts=tuple(_parse_cohort(label, value) for label, value in cohorts_raw.items()),
    )


def dumps(manifest: Manifest) -> str:
    """Render a manifest back to TOML, preserving cohort order."""
    lines = [tomlio.fmt_kv('seq_type', manifest.seq_type), tomlio.fmt_kv('generated', manifest.generated)]
    for cohort in manifest.cohorts:
        lines += ['', f'[cohorts.{cohort.label}]', tomlio.fmt_kv('uri', cohort.uri)]
        if cohort.analysis_id is not None:
            lines.append(tomlio.fmt_kv('analysis_id', cohort.analysis_id))
        if cohort.timestamp is not None:
            lines.append(tomlio.fmt_kv('timestamp', cohort.timestamp))
    return '\n'.join(lines) + '\n'


def save(manifest: Manifest, path: str | Path) -> None:
    """Write a manifest to disk."""
    Path(path).write_text(dumps(manifest))
```

- [ ] **Step 4: Run the test to verify it passes**

Run: `uv run --extra test pytest test/test_qc_calibration_manifest.py -q`

Expected: PASS (12 tests including parametrised cases)

- [ ] **Step 5: Lint and commit**

```bash
uvx ruff format src/align_genotype/qc_calibration/manifest.py test/test_qc_calibration_manifest.py
uvx ruff check src/align_genotype/qc_calibration/manifest.py test/test_qc_calibration_manifest.py
git add src/align_genotype/qc_calibration/manifest.py test/test_qc_calibration_manifest.py
git commit -m "feat(qc_calibration): cohort manifest as the reproducibility record"
```

---

## Task 4: `cache` — the small per-metric value cache

The whole workflow rests on parsing each huge MultiQC JSON exactly once. `collect` (Task 5) writes this cache; every tuning command reads it and completes in under a second, which is what makes threshold iteration practical.

Two guards live here. `complete = false` means a gated metric was missing during collect, so downstream commands must refuse to run. And the staleness rule: the spec's metric set must be a **subset** of the cache's, so narrowing the metric list is free but adding one forces a re-collect with an error that says so.

**Files:**
- Create: `src/align_genotype/qc_calibration/cache.py`
- Test: `test/test_qc_calibration_cache.py`

- [ ] **Step 1: Write the failing test**

Create `test/test_qc_calibration_cache.py`:

```python
"""Unit tests for the calibration value cache."""

import json

import numpy as np
import pytest

from align_genotype.qc_calibration import cache as cache_mod
from align_genotype.qc_calibration import spec as spec_mod
from align_genotype.qc_calibration.cache import CacheError, CohortValues, ValueCache

SPEC = spec_mod.loads(
    'seq_type = "genome"\ncache = "c.json"\n'
    '\n[metrics.MEDIAN_COVERAGE]\ndirection = "min"\nunit = "x"\nfail = 15\n'
    '\n[metrics.FREEMIX]\ndirection = "max"\nfail = 0.04\n',
)


def _cache(complete: bool = True, metrics: tuple[str, ...] = ('MEDIAN_COVERAGE', 'FREEMIX')) -> ValueCache:
    return ValueCache(
        seq_type='genome',
        generated='2026-08-11T13:52:00',
        complete=complete,
        metrics=metrics,
        cohorts=(
            CohortValues(
                label='dataset-a',
                n_samples=3,
                multiqc_version='1.33',
                shape='dict',
                n_dropped=1,
                values={'MEDIAN_COVERAGE': [30.0, 28.0, 12.0], 'FREEMIX': [0.001, 0.002, 0.05]},
            ),
        ),
    )


def test_save_and_load_round_trip(tmp_path):
    path = tmp_path / 'values.json'
    original = _cache()
    cache_mod.save(original, path)
    assert cache_mod.load(path) == original


def test_saved_json_is_readable_and_shaped_as_documented(tmp_path):
    path = tmp_path / 'values.json'
    cache_mod.save(_cache(), path)
    raw = json.loads(path.read_text())
    assert raw['seq_type'] == 'genome'
    assert raw['complete'] is True
    assert raw['metrics'] == ['MEDIAN_COVERAGE', 'FREEMIX']
    assert raw['cohorts']['dataset-a']['n_samples'] == 3
    assert raw['cohorts']['dataset-a']['values']['MEDIAN_COVERAGE'] == [30.0, 28.0, 12.0]


def test_save_refuses_to_write_nan(tmp_path):
    """NaN is filtered during collect; if one reaches here the cache must not go silently invalid."""
    bad = ValueCache(
        seq_type='genome',
        generated='x',
        complete=True,
        metrics=('M',),
        cohorts=(CohortValues('dataset-a', 1, '1.33', 'dict', 0, {'M': [float('nan')]}),),
    )
    with pytest.raises(ValueError, match='Out of range float values|NaN'):
        cache_mod.save(bad, tmp_path / 'values.json')


def test_labels_and_cohort_lookup():
    c = _cache()
    assert c.labels == ('dataset-a',)
    assert c.cohort('dataset-a').n_samples == 3
    with pytest.raises(KeyError, match='dataset-z'):
        c.cohort('dataset-z')


def test_series_returns_float_array():
    series = _cache().series('dataset-a', 'MEDIAN_COVERAGE')
    assert series.dtype == np.float64
    np.testing.assert_array_equal(series, np.array([30.0, 28.0, 12.0]))


def test_series_filters_nan_defensively():
    c = ValueCache('genome', 'x', True, ('M',), (CohortValues('dataset-a', 2, '1.33', 'dict', 0, {'M': [1.0, None]}),))
    np.testing.assert_array_equal(c.series('dataset-a', 'M'), np.array([1.0]))


def test_series_for_absent_metric_is_empty():
    assert _cache().series('dataset-a', 'NOT_COLLECTED').size == 0


def test_require_usable_accepts_a_superset_cache():
    cache_mod.require_usable(_cache(metrics=('MEDIAN_COVERAGE', 'FREEMIX', 'PCT_20X')), SPEC)  # does not raise


def test_require_usable_rejects_missing_metrics():
    with pytest.raises(CacheError, match=r"missing.*\['FREEMIX'\].*qc_calibrate collect"):
        cache_mod.require_usable(_cache(metrics=('MEDIAN_COVERAGE',)), SPEC)


def test_require_usable_rejects_incomplete_cache():
    with pytest.raises(CacheError, match='incomplete'):
        cache_mod.require_usable(_cache(complete=False), SPEC)
```

- [ ] **Step 2: Run the test to verify it fails**

Run: `uv run --extra test pytest test/test_qc_calibration_cache.py -q`

Expected: FAIL — `ModuleNotFoundError: No module named 'align_genotype.qc_calibration.cache'`

- [ ] **Step 3: Write the implementation**

Create `src/align_genotype/qc_calibration/cache.py`:

```python
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

from cpg_utils import to_path

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
        """This cohort's finite values for `metric`, as a float array (empty if absent)."""
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
                'values': c.values,
            }
            for c in cache.cohorts
        },
    }
    with to_path(path).open('w') as f:
        json.dump(payload, f, indent=2, allow_nan=False)


def load(path: str | Path) -> ValueCache:
    """Read a cache written by `save`."""
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
    the reports have to be parsed again.
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
```

- [ ] **Step 4: Run the test to verify it passes**

Run: `uv run --extra test pytest test/test_qc_calibration_cache.py -q`

Expected: PASS (11 tests)

- [ ] **Step 5: Lint and commit**

```bash
uvx ruff format src/align_genotype/qc_calibration/cache.py test/test_qc_calibration_cache.py
uvx ruff check src/align_genotype/qc_calibration/cache.py test/test_qc_calibration_cache.py
git add src/align_genotype/qc_calibration/cache.py test/test_qc_calibration_cache.py
git commit -m "feat(qc_calibration): per-metric value cache with completeness and staleness guards"
```

---

## Task 5: `collect` — survey and extraction in one parse

Surveying (is every gated metric actually present in general stats?) and extracting both require a full parse of every report, so doing them as separate commands would double the only expensive operation in the workflow. They are one pass here.

The survey is the guard against gating a key MultiQC renamed or never surfaces — the bug that made the old genome `PCT_PF_READS_ALIGNED` gate silently inert for who knows how long. Making it a hard error, with the cache marked `complete: false`, means the mandatory step is structurally unskippable rather than just documented as mandatory.

Reports are parsed one at a time and released (`gc.collect()` between cohorts); nothing ever holds two raw reports at once.

**Files:**
- Create: `src/align_genotype/qc_calibration/collect.py`
- Test: `test/test_qc_calibration_collect.py`

- [ ] **Step 1: Write the failing test**

Create `test/test_qc_calibration_collect.py`:

```python
"""Unit tests for the collect step: survey + extraction in a single parse."""

import json

import pytest

from align_genotype.qc_calibration import collect as collect_mod
from align_genotype.qc_calibration import manifest as manifest_mod
from align_genotype.qc_calibration import spec as spec_mod
from align_genotype.qc_calibration.collect import CollectError

SPEC = spec_mod.loads(
    'seq_type = "genome"\ncache = "c.json"\n'
    '\n[metrics.MEDIAN_COVERAGE]\ndirection = "min"\nunit = "x"\nfail = 15\n'
    '\n[metrics.FREEMIX]\ndirection = "max"\nfail = 0.04\n'
    '\n[metrics.error_rate]\ndirection = "max"\ngated = false\n',
)

DICT_SECTIONS = {
    'picard': {'CPG1|S1': {'MEDIAN_COVERAGE': 30}, 'CPG2|S2': {'MEDIAN_COVERAGE': '?'}},
    'verifybamid': {'CPG1|S1': {'FREEMIX': 0.001}, 'CPG2|S2': {'FREEMIX': 0.05}},
}


def _write_report(tmp_path, name, sections, version='1.33'):
    path = tmp_path / name
    path.write_text(json.dumps({'config_version': version, 'report_general_stats_data': sections}))
    return str(path)


def _manifest(**cohorts):
    body = ''.join(f'[cohorts.{label}]\nuri = "{uri}"\n\n' for label, uri in cohorts.items())
    return manifest_mod.loads(f'seq_type = "genome"\ngenerated = "x"\n\n{body}')


# --- single cohort ---------------------------------------------------------------


def test_collect_cohort_extracts_values(tmp_path):
    uri = _write_report(tmp_path, 'a.json', DICT_SECTIONS)
    values, row = collect_mod.collect_cohort(manifest_mod.Cohort('dataset-a', uri), SPEC)
    assert values.values['MEDIAN_COVERAGE'] == [30.0]
    assert sorted(values.values['FREEMIX']) == [0.001, 0.05]
    assert row.multiqc_version == '1.33'
    assert row.shape == 'dict'


def test_collect_cohort_counts_non_numeric_drops(tmp_path):
    uri = _write_report(tmp_path, 'a.json', DICT_SECTIONS)
    values, row = collect_mod.collect_cohort(manifest_mod.Cohort('dataset-a', uri), SPEC)
    assert values.n_dropped == 1  # the Picard '?' placeholder
    assert row.n_dropped == 1


def test_collect_cohort_filters_nan(tmp_path):
    uri = _write_report(tmp_path, 'a.json', {'picard': {'S1': {'MEDIAN_COVERAGE': 30}, 'S2': {'MEDIAN_COVERAGE': 'nan'}}})
    spec = spec_mod.loads('seq_type = "genome"\ncache = "c.json"\n[metrics.MEDIAN_COVERAGE]\ndirection = "min"\nfail = 1\n')
    values, _ = collect_mod.collect_cohort(manifest_mod.Cohort('dataset-a', uri), spec)
    assert values.values['MEDIAN_COVERAGE'] == [30.0]  # float('nan') coerces but isn't finite
    assert values.n_dropped == 1


def test_collect_cohort_dedupes_samples_across_sections(tmp_path):
    uri = _write_report(tmp_path, 'a.json', DICT_SECTIONS)
    values, _ = collect_mod.collect_cohort(manifest_mod.Cohort('dataset-a', uri), SPEC)
    assert values.n_samples == 2  # CPG1|S1 and CPG2|S2 appear in both sections


def test_collect_cohort_records_which_sections_carry_each_metric(tmp_path):
    uri = _write_report(tmp_path, 'a.json', DICT_SECTIONS)
    _, row = collect_mod.collect_cohort(manifest_mod.Cohort('dataset-a', uri), SPEC)
    assert row.where['MEDIAN_COVERAGE'] == ('picard',)
    assert row.where['FREEMIX'] == ('verifybamid',)
    assert row.where['error_rate'] == ()


def test_collect_cohort_handles_list_shaped_report(tmp_path):
    uri = _write_report(tmp_path, 'a.json', [{'S1': {'MEDIAN_COVERAGE': 30, 'FREEMIX': 0.01}}], version='1.14')
    values, row = collect_mod.collect_cohort(manifest_mod.Cohort('dataset-a', uri), SPEC)
    assert row.shape == 'list'
    assert row.where['MEDIAN_COVERAGE'] == ('section_0',)
    assert values.values['MEDIAN_COVERAGE'] == [30.0]


def test_collect_cohort_records_section_keys_for_rename_mapping(tmp_path):
    uri = _write_report(tmp_path, 'a.json', DICT_SECTIONS)
    _, row = collect_mod.collect_cohort(manifest_mod.Cohort('dataset-a', uri), SPEC)
    assert row.section_keys['picard'] == ('MEDIAN_COVERAGE',)


def test_collect_cohort_raises_without_general_stats(tmp_path):
    path = tmp_path / 'a.json'
    path.write_text(json.dumps({'config_version': '1.33', 'report_saved_raw_data': {}}))
    with pytest.raises(CollectError, match='no usable report_general_stats_data'):
        collect_mod.collect_cohort(manifest_mod.Cohort('dataset-a', str(path)), SPEC)


# --- all cohorts -----------------------------------------------------------------


def test_collect_all_builds_a_complete_cache(tmp_path):
    a = _write_report(tmp_path, 'a.json', DICT_SECTIONS)
    b = _write_report(tmp_path, 'b.json', DICT_SECTIONS)
    result = collect_mod.collect_all(_manifest(**{'dataset-a': a, 'dataset-b': b}), SPEC, generated='2026-08-11T00:00:00')
    assert result.ok
    assert result.cache.complete is True
    assert result.cache.labels == ('dataset-a', 'dataset-b')
    assert result.cache.metrics == ('MEDIAN_COVERAGE', 'FREEMIX', 'error_rate')


def test_collect_all_flags_missing_gated_metric_and_marks_cache_incomplete(tmp_path):
    good = _write_report(tmp_path, 'a.json', DICT_SECTIONS)
    without_freemix = _write_report(tmp_path, 'b.json', {'picard': {'S1': {'MEDIAN_COVERAGE': 30}}})
    result = collect_mod.collect_all(_manifest(**{'dataset-a': good, 'dataset-b': without_freemix}), SPEC)
    assert not result.ok
    assert result.missing_gated == {'dataset-b': ('FREEMIX',)}
    assert result.cache.complete is False


def test_collect_all_ignores_missing_ungated_metric(tmp_path):
    """error_rate is absent from every fixture but is gated = false, so it isn't fatal."""
    a = _write_report(tmp_path, 'a.json', DICT_SECTIONS)
    result = collect_mod.collect_all(_manifest(**{'dataset-a': a}), SPEC)
    assert result.ok


def test_collect_all_records_unreadable_cohort_but_continues(tmp_path):
    good = _write_report(tmp_path, 'a.json', DICT_SECTIONS)
    result = collect_mod.collect_all(
        _manifest(**{'dataset-a': good, 'dataset-b': str(tmp_path / 'nope.json')}),
        SPEC,
    )
    assert not result.ok
    assert result.cache.labels == ('dataset-a',)
    assert [label for label, _ in result.failures] == ['dataset-b']


def test_collect_all_rejects_seq_type_mismatch(tmp_path):
    a = _write_report(tmp_path, 'a.json', DICT_SECTIONS)
    exome_manifest = manifest_mod.loads(f'seq_type = "exome"\ngenerated = "x"\n\n[cohorts.dataset-a]\nuri = "{a}"\n')
    with pytest.raises(CollectError, match="manifest is for 'exome' but the spec is for 'genome'"):
        collect_mod.collect_all(exome_manifest, SPEC)
```

- [ ] **Step 2: Run the test to verify it fails**

Run: `uv run --extra test pytest test/test_qc_calibration_collect.py -q`

Expected: FAIL — `ModuleNotFoundError: No module named 'align_genotype.qc_calibration.collect'`

- [ ] **Step 3: Write the implementation**

Create `src/align_genotype/qc_calibration/collect.py`:

```python
"""Survey and extract in a single parse of each MultiQC report.

The survey answers "is every metric I intend to gate actually present in general
stats, in every cohort?". Getting that wrong is how a gate goes silently inert - the
genome reads-mapped gate was configured on ``PCT_PF_READS_ALIGNED``, which MultiQC only
writes to ``report_saved_raw_data``, so it checked nothing at all. A gated metric
missing from any cohort is therefore fatal here, and the cache it writes is marked
incomplete so nothing downstream will run on it.

Extraction goes through ``check_multiqc.normalise_sections`` and
``check_multiqc.gather_metric_values`` so calibration sees precisely what enforcement
will see.
"""

import gc
import json
import math
from dataclasses import dataclass
from datetime import datetime, timezone
from typing import Any

from cpg_utils import to_path

from align_genotype.qc_calibration.cache import CohortValues, ValueCache
from align_genotype.qc_calibration.manifest import Cohort, Manifest
from align_genotype.qc_calibration.spec import CalibrationSpec
from align_genotype.scripts import check_multiqc


class CollectError(RuntimeError):
    """A report could not be read, or is unusable for calibration."""


@dataclass(frozen=True)
class SurveyRow:
    """What one report looked like: provenance, shape, and where each metric lives."""

    label: str
    uri: str
    multiqc_version: str
    shape: str
    n_samples: int
    section_sizes: dict[str, int]
    where: dict[str, tuple[str, ...]]
    section_keys: dict[str, tuple[str, ...]]
    n_dropped: int


@dataclass(frozen=True)
class CollectResult:
    cache: ValueCache
    rows: tuple[SurveyRow, ...]
    missing_gated: dict[str, tuple[str, ...]]
    failures: tuple[tuple[str, str], ...]

    @property
    def ok(self) -> bool:
        return not self.missing_gated and not self.failures


def sections_carrying(sections: dict[str, Any], metric: str) -> tuple[str, ...]:
    """Section names in which at least one sample carries `metric`.

    Presence is judged on the key, not on whether its value parses as a number, so a
    metric that is present but entirely non-numeric still reads as present - the drop
    count is what surfaces that.
    """
    return tuple(
        name
        for name, section in sections.items()
        if any(metric in values for values in section.values() if isinstance(values, dict))
    )


def _all_keys(section: dict[str, Any]) -> tuple[str, ...]:
    keys: set[str] = set()
    for values in section.values():
        if isinstance(values, dict):
            keys |= set(values)
    return tuple(sorted(keys))


def _extract(sections: dict[str, Any], spec: CalibrationSpec) -> tuple[dict[str, list[float]], int]:
    values: dict[str, list[float]] = {}
    n_dropped = 0
    for metric in spec.metrics:
        entries, dropped = check_multiqc.gather_metric_values(sections, metric.key)
        finite = [value for _, _, value in entries if math.isfinite(value)]
        n_dropped += dropped + (len(entries) - len(finite))
        values[metric.key] = finite
    return values, n_dropped


def collect_cohort(cohort: Cohort, spec: CalibrationSpec) -> tuple[CohortValues, SurveyRow]:
    """Parse one report and distil it to values plus a survey row."""
    with to_path(cohort.uri).open() as f:
        document = json.load(f)

    version = str(document.get('config_version', 'unknown'))
    raw = document.get('report_general_stats_data')
    shape = 'dict' if isinstance(raw, dict) else 'list' if isinstance(raw, list) else type(raw).__name__
    sections = check_multiqc.normalise_sections(raw)
    if not sections:
        raise CollectError(
            f'{cohort.label}: no usable report_general_stats_data (multiqc {version}, shape {shape}) in {cohort.uri}',
        )

    n_samples = len({sample for section in sections.values() for sample in section})
    values, n_dropped = _extract(sections, spec)
    row = SurveyRow(
        label=cohort.label,
        uri=cohort.uri,
        multiqc_version=version,
        shape=shape,
        n_samples=n_samples,
        section_sizes={name: len(section) for name, section in sections.items()},
        where={metric.key: sections_carrying(sections, metric.key) for metric in spec.metrics},
        section_keys={name: _all_keys(section) for name, section in sections.items()},
        n_dropped=n_dropped,
    )
    cohort_values = CohortValues(
        label=cohort.label,
        n_samples=n_samples,
        multiqc_version=version,
        shape=shape,
        n_dropped=n_dropped,
        values=values,
    )
    return cohort_values, row


def collect_all(manifest: Manifest, spec: CalibrationSpec, generated: str | None = None) -> CollectResult:
    """Collect every cohort in `manifest`, one report at a time.

    A cohort that can't be read is recorded and the rest continue, so one bad URI
    doesn't waste a long run - but the result is not `ok`, and the caller must exit
    non-zero.
    """
    if manifest.seq_type != spec.seq_type:
        raise CollectError(
            f'manifest is for {manifest.seq_type!r} but the spec is for {spec.seq_type!r}',
        )

    collected: list[CohortValues] = []
    rows: list[SurveyRow] = []
    missing_gated: dict[str, tuple[str, ...]] = {}
    failures: list[tuple[str, str]] = []

    for cohort in manifest.cohorts:
        try:
            values, row = collect_cohort(cohort, spec)
        except (OSError, ValueError, CollectError) as exc:
            failures.append((cohort.label, str(exc)))
            continue
        finally:
            gc.collect()  # release the parsed report before the next one is read

        collected.append(values)
        rows.append(row)
        if gaps := tuple(metric.key for metric in spec.gated if not row.where.get(metric.key)):
            missing_gated[cohort.label] = gaps

    cache = ValueCache(
        seq_type=spec.seq_type,
        generated=generated or datetime.now(tz=timezone.utc).isoformat(timespec='seconds'),
        complete=not missing_gated and not failures,
        metrics=spec.metric_keys,
        cohorts=tuple(collected),
    )
    return CollectResult(cache=cache, rows=tuple(rows), missing_gated=missing_gated, failures=tuple(failures))
```

- [ ] **Step 4: Run the test to verify it passes**

Run: `uv run --extra test pytest test/test_qc_calibration_collect.py -q`

Expected: PASS (14 tests)

- [ ] **Step 5: Lint and commit**

```bash
uvx ruff format src/align_genotype/qc_calibration/collect.py test/test_qc_calibration_collect.py
uvx ruff check src/align_genotype/qc_calibration/collect.py test/test_qc_calibration_collect.py
git add src/align_genotype/qc_calibration/collect.py test/test_qc_calibration_collect.py
git commit -m "feat(qc_calibration): fold the mandatory key survey into extraction

Surveying and extracting both need a full parse of every report, so they run
as one pass. A gated metric missing from any cohort is fatal and marks the
cache incomplete - the guard against configuring a gate on a key MultiQC
never surfaces, which is how the old PCT_PF_READS_ALIGNED gate was inert."
```

---

## Task 6: `stats` — percentiles, flag rates, churn

Pure numeric functions over cached values, with no I/O. These back the `distributions` and `flagrates` commands and supply the churn calculation that `mad` uses in Task 7.

Two details that matter for correctness:
- **Warn excludes fail.** Production evaluates fail before warn and records one flag at the worst tier, so a warn rate that counted failing samples would overstate what operators actually see.
- **Churn uses production rounding.** `check_multiqc` rounds derived relative thresholds to 4 dp specifically to damp sub-0.0001 jitter that would otherwise churn "updated" flags in the database. The churn simulation has to round identically or it measures the wrong thing.

**Files:**
- Create: `src/align_genotype/qc_calibration/stats.py`
- Test: `test/test_qc_calibration_stats.py`

- [ ] **Step 1: Write the failing test**

Create `test/test_qc_calibration_stats.py`:

```python
"""Unit tests for the calibration statistics."""

import numpy as np
import pytest

from align_genotype.qc_calibration import stats


def test_percentiles_of_a_known_range():
    values = np.arange(1, 101, dtype=float)
    result = dict(zip(stats.PERCENTILES, stats.percentiles(values), strict=True))
    assert result[1] == pytest.approx(1.99)
    assert result[25] == pytest.approx(25.75)
    assert result[50] == pytest.approx(50.5)
    assert result[99] == pytest.approx(99.01)


def test_percentiles_of_empty_is_empty():
    assert stats.percentiles(np.array([])) == ()


def test_breach_min_flags_below():
    np.testing.assert_array_equal(
        stats.breach(np.array([1.0, 5.0, 10.0]), 5.0, 'min'),
        np.array([True, False, False]),
    )


def test_breach_max_flags_above():
    np.testing.assert_array_equal(
        stats.breach(np.array([1.0, 5.0, 10.0]), 5.0, 'max'),
        np.array([False, False, True]),
    )


def test_flag_rates_warn_excludes_fail():
    """Production records one flag at the worst tier, so warn must not double-count fails."""
    values = np.array([1.0, 2.0, 3.0, 100.0])  # max metric: fail>50, warn>2
    fail_rate, warn_rate = stats.flag_rates(values, 'max', fail=50, warn=2)
    assert fail_rate == pytest.approx(0.25)  # just the 100
    assert warn_rate == pytest.approx(0.25)  # just the 3; the 100 is already a fail


def test_flag_rates_with_only_one_tier():
    values = np.array([1.0, 2.0, 3.0, 100.0])
    assert stats.flag_rates(values, 'max', fail=50, warn=None) == (pytest.approx(0.25), 0.0)
    assert stats.flag_rates(values, 'max', fail=None, warn=2) == (0.0, pytest.approx(0.5))


def test_flag_rates_of_empty_is_nan():
    fail_rate, warn_rate = stats.flag_rates(np.array([]), 'max', fail=1, warn=2)
    assert np.isnan(fail_rate) and np.isnan(warn_rate)


def test_needs_review_uses_the_healthy_cohort_guardrail():
    assert stats.needs_review(0.0, 0.04) is False
    assert stats.needs_review(0.05, 0.04) is True  # too many fails
    assert stats.needs_review(0.0, 0.30) is True  # too many warns


# --- churn ---------------------------------------------------------------------

INITIAL = np.array([1.0, 2.0, 3.0, 4.0, 20.0])
GROWN = np.concatenate([INITIAL, np.array([18.0, 19.0, 20.0, 21.0, 22.0])])


def test_churn_measures_flips_caused_purely_by_cohort_growth():
    result = stats.churn(INITIAL, GROWN, 'max', 3.5)
    assert result.threshold_before == pytest.approx(8.189, abs=0.001)
    assert result.threshold_after == pytest.approx(34.067, abs=0.001)
    assert result.flagged_before == 1  # 20 is an outlier in the small cohort
    assert result.flagged_after == 0  # once the cohort includes similar values, it isn't
    assert result.flips == 1
    assert result.flip_rate == pytest.approx(0.2)


def test_churn_rounds_to_four_dp_like_production():
    result = stats.churn(INITIAL, GROWN, 'max', 3.5)
    assert result.threshold_before == round(result.threshold_before, 4)
    assert result.threshold_after == round(result.threshold_after, 4)


def test_churn_returns_none_on_degenerate_mad():
    identical = np.array([5.0] * 10)
    assert stats.churn(identical, identical, 'max', 3.5) is None


def test_churn_flip_rate_of_empty_initial_is_zero():
    assert stats.churn(np.array([]), GROWN, 'max', 3.5) is None
```

- [ ] **Step 2: Run the test to verify it fails**

Run: `uv run --extra test pytest test/test_qc_calibration_stats.py -q`

Expected: FAIL — `ModuleNotFoundError: No module named 'align_genotype.qc_calibration.stats'`

- [ ] **Step 3: Write the implementation**

Create `src/align_genotype/qc_calibration/stats.py`:

```python
"""Numeric analysis over cached values - percentiles, flag rates, cohort-growth churn.

No I/O and no config: everything here is a pure function of arrays already in memory,
which is what lets the flagrates/mad tuning loop run in under a second.
"""

from dataclasses import dataclass

import numpy as np

from align_genotype.scripts import check_multiqc

PERCENTILES: tuple[int, ...] = (1, 5, 10, 25, 50, 75, 90, 95, 99)

# A healthy cohort should sit near 0% fail and single-digit % warn. Beyond these, the
# candidate threshold gets flagged for a second look - it is advice, not a rejection:
# "healthy cohort" is the operator's judgement, not a computable property.
FAIL_RATE_LIMIT = 0.02
WARN_RATE_LIMIT = 0.10


def percentiles(values: np.ndarray, pcts: tuple[int, ...] = PERCENTILES) -> tuple[float, ...]:
    """Percentiles of `values`, in the order given. Empty input yields an empty tuple."""
    if values.size == 0:
        return ()
    return tuple(float(p) for p in np.percentile(values, pcts))


def breach(values: np.ndarray, threshold: float, direction: str) -> np.ndarray:
    """Boolean mask of samples on the bad side of `threshold`.

    'min' = higher is better, so low values breach; 'max' = lower is better.
    """
    return values < threshold if direction == 'min' else values > threshold


def flag_rates(
    values: np.ndarray,
    direction: str,
    fail: float | None = None,
    warn: float | None = None,
) -> tuple[float, float]:
    """``(fail_rate, warn_rate)`` as fractions, warn excluding samples already failing.

    Mirrors production, which evaluates fail before warn and records one flag per
    metric at the worst tier. Empty input yields ``(nan, nan)``.
    """
    if values.size == 0:
        return (float('nan'), float('nan'))
    is_fail = breach(values, fail, direction) if fail is not None else np.zeros(values.size, dtype=bool)
    is_warn = breach(values, warn, direction) & ~is_fail if warn is not None else np.zeros(values.size, dtype=bool)
    return float(is_fail.mean()), float(is_warn.mean())


def needs_review(fail_rate: float, warn_rate: float) -> bool:
    """Whether a candidate threshold's flag rate is outside the healthy-cohort guardrail."""
    if np.isnan(fail_rate) or np.isnan(warn_rate):
        return False
    return fail_rate > FAIL_RATE_LIMIT or warn_rate > WARN_RATE_LIMIT


@dataclass(frozen=True)
class ChurnResult:
    """How a cohort-relative threshold moved, and who changed status because of it."""

    threshold_before: float
    threshold_after: float
    n_initial: int
    flagged_before: int
    flagged_after: int
    flips: int

    @property
    def flip_rate(self) -> float:
        return self.flips / self.n_initial if self.n_initial else 0.0


def churn(initial: np.ndarray, grown: np.ndarray, direction: str, k: float) -> ChurnResult | None:
    """Re-score the *initial* samples against the *grown* cohort's threshold.

    `flips` counts samples whose flag status changes purely because the cohort grew -
    each one would be a spurious "updated" flag in the database. Returns None when
    either cohort has a degenerate (zero) MAD.

    Thresholds are rounded to 4 dp exactly as production does, so the simulation
    measures the churn operators would actually see rather than sub-0.0001 jitter.
    """
    before_threshold = check_multiqc.robust_threshold(list(initial), direction, k)
    after_threshold = check_multiqc.robust_threshold(list(grown), direction, k)
    if before_threshold is None or after_threshold is None:
        return None
    before_threshold = round(before_threshold, 4)
    after_threshold = round(after_threshold, 4)
    before = breach(initial, before_threshold, direction)
    after = breach(initial, after_threshold, direction)
    return ChurnResult(
        threshold_before=before_threshold,
        threshold_after=after_threshold,
        n_initial=int(initial.size),
        flagged_before=int(before.sum()),
        flagged_after=int(after.sum()),
        flips=int((before != after).sum()),
    )
```

- [ ] **Step 4: Run the test to verify it passes**

Run: `uv run --extra test pytest test/test_qc_calibration_stats.py -q`

Expected: PASS (13 tests)

- [ ] **Step 5: Lint and commit**

```bash
uvx ruff format src/align_genotype/qc_calibration/stats.py test/test_qc_calibration_stats.py
uvx ruff check src/align_genotype/qc_calibration/stats.py test/test_qc_calibration_stats.py
git add src/align_genotype/qc_calibration/stats.py test/test_qc_calibration_stats.py
git commit -m "feat(qc_calibration): percentiles, flag rates and cohort-growth churn

Warn rates exclude samples already failing, and churn rounds thresholds to
4 dp, both matching what production actually does."
```

---

## Task 7: `relative` — MAD evaluation through the production code path

The manual process had a MAD prototype *and* a second script that verified the prototype against production. Keeping both invites exactly the divergence the verification was meant to catch. Here there is one implementation: warn rates come from calling `check_multiqc.relative_flags()` for real.

Driving it needs config, and rather than monkeypatching `config_retrieve` (what the old scripts did), this writes a throwaway TOML and points `cpg_utils.config` at it. `set_config_paths` only validates that files exist, end in `.toml`, and parse — no other keys are required. So the numbers come from the genuine config → `load_thresholds` → `relative_flags` path, which also means this step incidentally proves the config shape is loadable.

**Files:**
- Create: `src/align_genotype/qc_calibration/relative.py`
- Test: `test/test_qc_calibration_relative.py`

- [ ] **Step 1: Write the failing test**

Create `test/test_qc_calibration_relative.py`:

```python
"""Unit tests for cohort-relative (MAD) evaluation."""

import os

import pytest

from cpg_utils import config

from align_genotype.qc_calibration import relative as relative_mod
from align_genotype.qc_calibration import spec as spec_mod
from align_genotype.qc_calibration.cache import CohortValues, ValueCache
from align_genotype.qc_calibration.relative import CohortMad, MadEvaluation
from align_genotype.qc_calibration.stats import ChurnResult

METRIC = spec_mod.loads(
    'seq_type = "genome"\ncache = "c.json"\n'
    '[metrics.DUP]\ndirection = "max"\nunit = "%"\nfail = 1000.0\n'
    '[metrics.DUP.relative]\nk = 3.5\nmin_cohort = 5\n',
).metric('DUP')

# 49 evenly spread values plus one extreme: median 25.5, raw MAD 12.5,
# so the modified-z threshold lands at 25.5 + 3.5*12.5/0.6745 = 90.3629.
OUTLIER_COHORT = [*range(1, 50), 500.0]


@pytest.fixture(autouse=True)
def _restore_config_paths():
    """Keep the global cpg-utils config state from leaking between tests."""
    previous = os.environ.get('CPG_CONFIG_PATH', '')
    yield
    config.set_config_paths([p for p in previous.split(',') if p])


def _cache(**cohorts: list[float]) -> ValueCache:
    return ValueCache(
        seq_type='genome',
        generated='x',
        complete=True,
        metrics=('DUP',),
        cohorts=tuple(
            CohortValues(label, len(values), '1.33', 'dict', 0, {'DUP': [float(v) for v in values]})
            for label, values in cohorts.items()
        ),
    )


# --- per-cohort numbers -----------------------------------------------------------


def test_evaluate_reports_median_mad_and_threshold():
    evaluation = relative_mod.evaluate(_cache(**{'dataset-a': OUTLIER_COHORT}), METRIC, 'genome')
    cohort = evaluation.cohorts[0]
    assert cohort.n == 50
    assert cohort.median == pytest.approx(25.5)
    assert cohort.mad_raw == pytest.approx(12.5)
    assert cohort.threshold == pytest.approx(90.363, abs=0.001)


def test_warn_count_comes_from_the_production_relative_flags_path():
    evaluation = relative_mod.evaluate(_cache(**{'dataset-a': OUTLIER_COHORT}), METRIC, 'genome')
    assert evaluation.cohorts[0].n_warn == 1
    assert evaluation.cohorts[0].warn_rate == pytest.approx(0.02)


def test_cohort_below_min_cohort_is_skipped():
    evaluation = relative_mod.evaluate(_cache(**{'dataset-a': [1.0, 2.0, 3.0]}), METRIC, 'genome')
    cohort = evaluation.cohorts[0]
    assert cohort.skipped is not None
    assert 'min_cohort' in cohort.skipped
    assert cohort.n_warn == 0


def test_degenerate_mad_cohort_is_skipped():
    evaluation = relative_mod.evaluate(_cache(**{'dataset-a': [5.0] * 10}), METRIC, 'genome')
    cohort = evaluation.cohorts[0]
    assert cohort.threshold is None
    assert 'zero MAD' in cohort.skipped
    assert cohort.n_warn == 0


def test_evaluate_restores_previous_config_paths(tmp_path):
    existing = tmp_path / 'existing.toml'
    existing.write_text('[workflow]\nsequencing_type = "genome"\n')
    config.set_config_paths([str(existing)])
    relative_mod.evaluate(_cache(**{'dataset-a': OUTLIER_COHORT}), METRIC, 'genome')
    assert config.get_config_paths() == [str(existing)]


# --- churn --------------------------------------------------------------------


def test_evaluate_simulates_homogeneous_and_heterogeneous_growth():
    evaluation = relative_mod.evaluate(
        _cache(**{'dataset-a': OUTLIER_COHORT, 'dataset-b': [float(v) for v in range(200, 250)]}),
        METRIC,
        'genome',
    )
    assert {label for label, _ in evaluation.homogeneous} == {'dataset-a', 'dataset-b'}
    assert {(a, b) for a, b, _ in evaluation.heterogeneous} == {('dataset-a', 'dataset-b')}


def test_cohorts_too_small_are_excluded_from_churn():
    evaluation = relative_mod.evaluate(_cache(**{'dataset-a': [1.0, 2.0, 3.0]}), METRIC, 'genome')
    assert evaluation.homogeneous == ()
    assert evaluation.heterogeneous == ()


# --- verdict ------------------------------------------------------------------


def _evaluation(warn_rates: list[float], flip_rates: list[float]) -> MadEvaluation:
    cohorts = tuple(
        CohortMad(f'cohort-{i}', n=100, median=1.0, mad_raw=1.0, threshold=5.0, n_warn=int(rate * 100), skipped=None)
        for i, rate in enumerate(warn_rates)
    )
    churns = tuple(
        (f'cohort-{i}', ChurnResult(1.0, 1.0, 100, 0, 0, int(rate * 100))) for i, rate in enumerate(flip_rates)
    )
    return MadEvaluation('DUP', 'max', cohorts, churns, ())


def test_verdict_recommends_when_warn_and_churn_are_low():
    evaluation = _evaluation(warn_rates=[0.0, 0.042], flip_rates=[0.01, 0.021])
    assert evaluation.max_warn_rate == pytest.approx(0.042)
    assert evaluation.max_churn == pytest.approx(0.021)
    assert evaluation.verdict == 'RECOMMEND'


def test_verdict_rejects_on_high_churn():
    evaluation = _evaluation(warn_rates=[0.02], flip_rates=[0.245])
    assert evaluation.verdict == 'REJECT'
    assert 'churn' in evaluation.verdict_reason


def test_verdict_rejects_on_high_warn_rate():
    evaluation = _evaluation(warn_rates=[0.4], flip_rates=[0.0])
    assert evaluation.verdict == 'REJECT'
    assert 'warn' in evaluation.verdict_reason
```

- [ ] **Step 2: Run the test to verify it fails**

Run: `uv run --extra test pytest test/test_qc_calibration_relative.py -q`

Expected: FAIL — `ModuleNotFoundError: No module named 'align_genotype.qc_calibration.relative'`

- [ ] **Step 3: Write the implementation**

Create `src/align_genotype/qc_calibration/relative.py`:

```python
"""Evaluate a cohort-relative (MAD) warn tier - through the production code path.

The manual workflow had a prototype and a separate script to check the prototype
against production. There is one implementation here instead: warn rates come from
calling ``check_multiqc.relative_flags()`` for real, against a throwaway config file,
so the numbers reported are by construction the numbers that will ship.

A relative tier is only worth adopting where the metric's normal level genuinely shifts
by cohort or protocol, *and* the flag set stays stable as the cohort grows. Churn is the
second half of that test: every flip is a spurious "updated" flag in the database.
"""

import itertools
import tempfile
from contextlib import contextmanager
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Iterator

import numpy as np

from cpg_utils import config

from align_genotype.qc_calibration.cache import ValueCache
from align_genotype.qc_calibration.spec import MetricSpec, RelativeSpec
from align_genotype.qc_calibration.stats import ChurnResult, churn
from align_genotype.scripts import check_multiqc

# Adoption bar, from the exome and genome calibrations: ZERO_CVG_TARGETS_PCT and
# reads_duplicated_percent cleared it; PCT_SELECTED_BASES / PCT_OFF_BAIT did not
# (churn up to 24.5%).
MAX_WARN_RATE = 0.10
MAX_CHURN = 0.02

# The flag date is irrelevant to a rate calculation, but relative_flags() needs one.
_FIXED_DATE = datetime(2000, 1, 1, tzinfo=timezone.utc)

# Fraction of a cohort used as the "before" set when simulating homogeneous growth.
_GROWTH_FRACTION = 0.6


@dataclass(frozen=True)
class CohortMad:
    """One cohort's relative threshold and how many samples it warns."""

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
    metric: str
    direction: str
    cohorts: tuple[CohortMad, ...]
    homogeneous: tuple[tuple[str, ChurnResult], ...]
    heterogeneous: tuple[tuple[str, str, ChurnResult], ...]

    @property
    def max_warn_rate(self) -> float:
        rates = [c.warn_rate for c in self.cohorts if c.skipped is None]
        return max(rates) if rates else 0.0

    @property
    def max_churn(self) -> float:
        rates = [r.flip_rate for _, r in self.homogeneous] + [r.flip_rate for _, _, r in self.heterogeneous]
        return max(rates) if rates else 0.0

    @property
    def verdict(self) -> str:
        return 'REJECT' if self.verdict_reason else 'RECOMMEND'

    @property
    def verdict_reason(self) -> str:
        """Why this metric was rejected, or '' if it clears the bar."""
        reasons = []
        if self.max_warn_rate > MAX_WARN_RATE:
            reasons.append(f'peak warn rate {self.max_warn_rate:.1%} exceeds {MAX_WARN_RATE:.0%}')
        if self.max_churn > MAX_CHURN:
            reasons.append(f'peak churn {self.max_churn:.1%} exceeds {MAX_CHURN:.0%}')
        return '; '.join(reasons)


@contextmanager
def _production_config(seq_type: str, metric: MetricSpec, relative: RelativeSpec) -> Iterator[None]:
    """Point cpg-utils config at a throwaway TOML holding just this relative tier.

    Using the real config mechanism rather than monkeypatching config_retrieve means
    relative_flags() resolves its spec exactly as it will in production - and proves
    the config shape parses.
    """
    previous = config.get_config_paths()
    with tempfile.TemporaryDirectory() as tmpdir:
        path = Path(tmpdir) / 'calibration.toml'
        path.write_text(
            f'[workflow]\nsequencing_type = "{seq_type}"\n\n'
            f'[qc_thresholds.{seq_type}.relative.{metric.key}]\n'
            f'direction = "{metric.direction}"\n'
            f'k = {relative.k}\n'
            f'min_cohort = {relative.min_cohort}\n',
        )
        config.set_config_paths([str(path)])
        try:
            yield
        finally:
            config.set_config_paths(list(previous))


def _warn_count(values: np.ndarray, metric: MetricSpec, seq_type: str) -> int:
    """Warn flags the production relative pass raises for this cohort."""
    sections = {'calibration': {f'S{i}': {metric.key: float(v)} for i, v in enumerate(values)}}
    flags = check_multiqc.relative_flags(sections, seq_type, _FIXED_DATE, already_flagged={})
    return len(flags)


def _evaluate_cohort(
    label: str,
    values: np.ndarray,
    metric: MetricSpec,
    relative: RelativeSpec,
    seq_type: str,
) -> CohortMad:
    n = int(values.size)
    median = float(np.median(values)) if n else 0.0
    mad_raw = float(np.median(np.abs(values - median))) if n else 0.0

    if n < relative.min_cohort:
        skipped = f'cohort {n} < min_cohort {relative.min_cohort}; MAD is too noisy, use the absolute gate'
        return CohortMad(label, n, median, mad_raw, None, 0, skipped)

    threshold = check_multiqc.robust_threshold(list(values), metric.direction, relative.k)
    if threshold is None:
        return CohortMad(label, n, median, mad_raw, None, 0, 'zero MAD (degenerate cohort)')

    return CohortMad(
        label=label,
        n=n,
        median=median,
        mad_raw=mad_raw,
        threshold=round(threshold, 4),
        n_warn=_warn_count(values, metric, seq_type),
        skipped=None,
    )


def evaluate(cache: ValueCache, metric: MetricSpec, seq_type: str) -> MadEvaluation:
    """Per-cohort warn rates plus cohort-growth churn for one relative metric."""
    relative = metric.relative
    if relative is None:
        raise ValueError(f'metric {metric.key!r} has no [metrics.{metric.key}.relative] block to evaluate')

    series = {label: cache.series(label, metric.key) for label in cache.labels}

    with _production_config(seq_type, metric, relative):
        cohorts = tuple(
            _evaluate_cohort(label, values, metric, relative, seq_type) for label, values in series.items()
        )

    usable = {label: values for label, values in series.items() if values.size >= relative.min_cohort}

    homogeneous = tuple(
        (label, result)
        for label, values in usable.items()
        for result in [churn(values[: int(values.size * _GROWTH_FRACTION)], values, metric.direction, relative.k)]
        if result is not None
    )
    heterogeneous = tuple(
        (a, b, result)
        for a, b in itertools.combinations(usable, 2)
        for result in [churn(usable[a], np.concatenate([usable[a], usable[b]]), metric.direction, relative.k)]
        if result is not None
    )
    return MadEvaluation(metric.key, metric.direction, cohorts, homogeneous, heterogeneous)
```

- [ ] **Step 4: Run the test to verify it passes**

Run: `uv run --extra test pytest test/test_qc_calibration_relative.py -q`

Expected: PASS (11 tests)

- [ ] **Step 5: Lint and commit**

```bash
uvx ruff format src/align_genotype/qc_calibration/relative.py test/test_qc_calibration_relative.py
uvx ruff check src/align_genotype/qc_calibration/relative.py test/test_qc_calibration_relative.py
git add src/align_genotype/qc_calibration/relative.py test/test_qc_calibration_relative.py
git commit -m "feat(qc_calibration): MAD evaluation driven through production relative_flags

Replaces the prototype-plus-verification script pair with a single path: warn
rates come from calling relative_flags() against a throwaway config file, so
there is no second implementation of the modified z-score to diverge."
```

---

## Task 8: `suggest` — seed a first draft from the distributions

Saves the operator starting from a blank spec, without pretending the numbers are decisions. Everything seeded is written `reviewed = false`, and `emit-config` (Task 9) refuses to emit an unreviewed gated metric — so a seed physically cannot reach production config without a human changing that flag.

Seeds come from the tail the metric's direction cares about: for `min` metrics the low tail (`fail` ← lowest p1 across cohorts, `warn` ← lowest p5); for `max` metrics the high tail (`fail` ← highest p99, `warn` ← highest p95). A metric with a `relative` block gets only a `fail` seed, because its warn tier is cohort-relative by construction.

**Files:**
- Create: `src/align_genotype/qc_calibration/suggest.py`
- Test: `test/test_qc_calibration_suggest.py`

- [ ] **Step 1: Write the failing test**

Create `test/test_qc_calibration_suggest.py`:

```python
"""Unit tests for threshold seeding."""

import pytest

from align_genotype.qc_calibration import spec as spec_mod
from align_genotype.qc_calibration import suggest as suggest_mod
from align_genotype.qc_calibration.cache import CohortValues, ValueCache

SPEC = spec_mod.loads(
    'seq_type = "genome"\ncache = "c.json"\n'
    # min metric, no thresholds yet
    '\n[metrics.MEDIAN_COVERAGE]\ndirection = "min"\nunit = "x"\n gated = true\nfail = 1\n'
    # max metric, no thresholds yet
    '\n[metrics.FREEMIX]\ndirection = "max"\nunit = "frac"\nfail = 1\n'
    # relative metric - warn tier is cohort-relative, so only fail should be seeded
    '\n[metrics.DUP]\ndirection = "max"\nunit = "%"\nfail = 1\n[metrics.DUP.relative]\nk = 3.5\nmin_cohort = 5\n'
    # already signed off - must not be touched
    '\n[metrics.PCT_20X]\ndirection = "min"\nunit = "frac"\nfail = 0.75\nwarn = 0.85\nreviewed = true\n'
    # not gated - must not be touched
    '\n[metrics.error_rate]\ndirection = "max"\nunit = "frac"\ngated = false\n',
)

# Cohort A: 10..109, Cohort B: 20..119. Lowest p1 = 10.99, lowest p5 = 14.95.
# FREEMIX 0.000..0.099 in both: highest p99 = 0.09801, highest p95 = 0.09405.
CACHE = ValueCache(
    seq_type='genome',
    generated='x',
    complete=True,
    metrics=('MEDIAN_COVERAGE', 'FREEMIX', 'DUP', 'PCT_20X', 'error_rate'),
    cohorts=(
        CohortValues(
            'dataset-a',
            100,
            '1.33',
            'dict',
            0,
            {
                'MEDIAN_COVERAGE': [float(v) for v in range(10, 110)],
                'FREEMIX': [v / 1000 for v in range(100)],
                'DUP': [float(v) for v in range(100)],
                'PCT_20X': [0.9] * 100,
                'error_rate': [0.01] * 100,
            },
        ),
        CohortValues(
            'dataset-b',
            100,
            '1.33',
            'dict',
            0,
            {
                'MEDIAN_COVERAGE': [float(v) for v in range(20, 120)],
                'FREEMIX': [v / 1000 for v in range(100)],
                'DUP': [float(v) for v in range(100)],
                'PCT_20X': [0.9] * 100,
                'error_rate': [0.01] * 100,
            },
        ),
    ),
)


def test_min_metric_seeded_from_the_low_tail():
    updated, _ = suggest_mod.seed(CACHE, SPEC)
    metric = updated.metric('MEDIAN_COVERAGE')
    assert metric.fail == 11  # lowest p1 across cohorts (10.99), rounded for unit 'x'
    assert metric.warn == 15  # lowest p5 across cohorts (14.95)


def test_max_metric_seeded_from_the_high_tail():
    updated, _ = suggest_mod.seed(CACHE, SPEC)
    metric = updated.metric('FREEMIX')
    assert metric.fail == pytest.approx(0.1)  # highest p99 (0.09801), 2 dp for unit 'frac'
    assert metric.warn == pytest.approx(0.09)  # highest p95 (0.09405)


def test_relative_metric_seeds_fail_only():
    """Its warn tier is cohort-relative; an absolute warn alongside would double-flag."""
    updated, _ = suggest_mod.seed(CACHE, SPEC)
    assert updated.metric('DUP').fail is not None
    assert updated.metric('DUP').warn is None


def test_seeded_metrics_are_marked_unreviewed():
    updated, _ = suggest_mod.seed(CACHE, SPEC)
    assert updated.metric('MEDIAN_COVERAGE').reviewed is False


def test_seeded_metrics_record_their_evidence():
    updated, seeded = suggest_mod.seed(CACHE, SPEC)
    assert 'p1' in updated.metric('MEDIAN_COVERAGE').rationale
    assert {s.key for s in seeded} == {'MEDIAN_COVERAGE', 'FREEMIX', 'DUP'}


def test_reviewed_metric_is_never_overwritten():
    updated, seeded = suggest_mod.seed(CACHE, SPEC)
    assert updated.metric('PCT_20X') == SPEC.metric('PCT_20X')
    assert 'PCT_20X' not in {s.key for s in seeded}


def test_ungated_metric_is_never_seeded():
    updated, seeded = suggest_mod.seed(CACHE, SPEC)
    assert updated.metric('error_rate') == SPEC.metric('error_rate')
    assert 'error_rate' not in {s.key for s in seeded}


def test_seed_does_not_mutate_the_input_spec():
    suggest_mod.seed(CACHE, SPEC)
    assert SPEC.metric('MEDIAN_COVERAGE').warn is None


def test_metric_with_no_data_is_reported_not_seeded():
    empty = ValueCache(
        seq_type='genome',
        generated='x',
        complete=True,
        metrics=('MEDIAN_COVERAGE',),
        cohorts=(CohortValues('dataset-a', 0, '1.33', 'dict', 0, {'MEDIAN_COVERAGE': []}),),
    )
    spec = spec_mod.loads('seq_type = "genome"\ncache = "c.json"\n[metrics.MEDIAN_COVERAGE]\ndirection = "min"\nfail = 1\n')
    updated, seeded = suggest_mod.seed(empty, spec)
    assert seeded == ()
    assert updated == spec
```

- [ ] **Step 2: Run the test to verify it fails**

Run: `uv run --extra test pytest test/test_qc_calibration_suggest.py -q`

Expected: FAIL — `ModuleNotFoundError: No module named 'align_genotype.qc_calibration.suggest'`

- [ ] **Step 3: Write the implementation**

Create `src/align_genotype/qc_calibration/suggest.py`:

```python
"""Seed candidate thresholds from the cohort distributions.

These are starting points, not answers. The judgement the manual calibrations turned
on - preserving the lab's intent for a hard gate, refusing to hard-fail on metrics that
track ancestry or chemistry, relaxing genome duplication from 25 to 40 once the data
showed 25 would fail a third of legitimate preps - is not derivable from percentiles.
So everything seeded here is marked ``reviewed = false``, and emit-config will not
write an unreviewed gated metric into config.
"""

from dataclasses import dataclass

import numpy as np

from align_genotype.qc_calibration.cache import ValueCache
from align_genotype.qc_calibration.spec import CalibrationSpec, MetricSpec

# Which percentile seeds which tier, by direction. 'min' metrics go bad in the low
# tail, 'max' metrics in the high tail.
_SEED_PERCENTILES = {'min': {'fail': 1, 'warn': 5}, 'max': {'fail': 99, 'warn': 95}}


@dataclass(frozen=True)
class Seeded:
    """One metric's seeded values and the evidence behind them."""

    key: str
    fail: float | None
    warn: float | None
    evidence: str


def _round_for_unit(value: float, unit: str) -> float:
    """Round to a precision that reads sensibly for the unit."""
    if unit in ('x', '%'):
        return int(round(value))
    return round(value, 2)


def _tail(cache: ValueCache, metric: MetricSpec, percentile: int) -> float | None:
    """The worst per-cohort value of `percentile` - lowest for 'min', highest for 'max'.

    Taking the worst rather than the average means a seeded threshold doesn't already
    flag a large slice of the most marginal cohort in the set.
    """
    per_cohort = [
        float(np.percentile(series, percentile))
        for label in cache.labels
        if (series := cache.series(label, metric.key)).size
    ]
    if not per_cohort:
        return None
    return min(per_cohort) if metric.direction == 'min' else max(per_cohort)


def _seed_metric(cache: ValueCache, metric: MetricSpec) -> Seeded | None:
    percentiles = _SEED_PERCENTILES[metric.direction]
    fail_raw = _tail(cache, metric, percentiles['fail'])
    if fail_raw is None:
        return None
    fail = _round_for_unit(fail_raw, metric.unit)

    # A relative metric's warn tier is cohort-derived, so seeding an absolute warn
    # would create the double-flagging the spec validation forbids.
    warn = None
    warn_note = ''
    if metric.relative is None:
        warn_raw = _tail(cache, metric, percentiles['warn'])
        if warn_raw is not None:
            warn = _round_for_unit(warn_raw, metric.unit)
            warn_note = f', warn from p{percentiles["warn"]} ({warn_raw:.4g})'

    evidence = (
        f'Seeded from cohort percentiles across {len(cache.labels)} cohorts: '
        f'fail from p{percentiles["fail"]} ({fail_raw:.4g}){warn_note}. Not yet reviewed.'
    )
    return Seeded(key=metric.key, fail=fail, warn=warn, evidence=evidence)


def seed(cache: ValueCache, spec: CalibrationSpec) -> tuple[CalibrationSpec, tuple[Seeded, ...]]:
    """Return a copy of `spec` with unreviewed gated metrics seeded, plus what changed.

    Metrics that are un-gated, already ``reviewed = true``, or absent from every cohort
    are left exactly as they are.
    """
    updated = spec
    seeded: list[Seeded] = []
    for metric in spec.gated:
        if metric.reviewed:
            continue
        if (result := _seed_metric(cache, metric)) is None:
            continue
        updated = updated.with_metric(
            replace_key=metric.key,
            fail=result.fail,
            warn=result.warn,
            reviewed=False,
            rationale=result.evidence,
        )
        seeded.append(result)
    return updated, tuple(seeded)
```

- [ ] **Step 4: Run the test to verify it passes**

Run: `uv run --extra test pytest test/test_qc_calibration_suggest.py -q`

Expected: PASS (9 tests)

- [ ] **Step 5: Lint and commit**

```bash
uvx ruff format src/align_genotype/qc_calibration/suggest.py test/test_qc_calibration_suggest.py
uvx ruff check src/align_genotype/qc_calibration/suggest.py test/test_qc_calibration_suggest.py
git add src/align_genotype/qc_calibration/suggest.py test/test_qc_calibration_suggest.py
git commit -m "feat(qc_calibration): seed candidate thresholds from cohort percentiles

Seeds are written reviewed = false so they cannot reach config without a
human signing off; the judgement calls a calibration turns on aren't
derivable from a distribution."
```

---

## Task 9: `emit` — generate the config block

Prints the `[qc_thresholds.<seq_type>...]` block to stdout for the operator to paste. It never writes to `config_template.toml`; the tool doesn't get to clobber hand-maintained config.

> **Correction to the design spec.** The spec said the header should record "the cohort count and labels". It must record the **count only**. This block is pasted into `config_template.toml`, which is committed and pushed — and cohort labels are real CPG dataset names. The manifest (in gitignored `calibration/`) is where the labels live; the header cites its path instead.

Each metric's comment is its spec `rationale` followed by generated evidence — the cohort median range and the observed flag-rate range at that threshold — matching the hand-written style already in the config.

**Files:**
- Create: `src/align_genotype/qc_calibration/emit.py`
- Create: `src/align_genotype/qc_calibration/report.py` (just the measurement formatter; Task 12 adds the tables)
- Test: `test/test_qc_calibration_emit.py`

- [ ] **Step 1: Write the failing test**

Create `test/test_qc_calibration_emit.py`:

```python
"""Golden test for the emitted config block."""

import os

import pytest

from cpg_utils import config

from align_genotype.qc_calibration import emit as emit_mod
from align_genotype.qc_calibration import spec as spec_mod
from align_genotype.qc_calibration import tomlio
from align_genotype.qc_calibration.cache import CohortValues, ValueCache
from align_genotype.qc_calibration.emit import EmitError
from align_genotype.scripts import check_multiqc

SPEC = spec_mod.loads(
    'seq_type = "genome"\ncache = "c.json"\n'
    '\n[metrics.MEDIAN_COVERAGE]\ndirection = "min"\nunit = "x"\nfail = 15\nwarn = 25\nreviewed = true\n'
    'rationale = "Primary depth gate."\n'
    '\n[metrics.FREEMIX]\ndirection = "max"\nunit = "frac"\nfail = 0.04\nwarn = 0.02\nreviewed = true\n'
    'rationale = "Contamination safety gate."\n'
    '\n[metrics.DUP]\ndirection = "max"\nunit = "%"\nfail = 40\nreviewed = true\n'
    'rationale = "Library-prep dependent."\n[metrics.DUP.relative]\nk = 3.5\nmin_cohort = 50\n'
    '\n[metrics.error_rate]\ndirection = "max"\nunit = "frac"\ngated = false\n',
)

CACHE = ValueCache(
    seq_type='genome',
    generated='2026-08-11T00:00:00',
    complete=True,
    metrics=('MEDIAN_COVERAGE', 'FREEMIX', 'DUP', 'error_rate'),
    cohorts=(
        CohortValues(
            'dataset-a',
            5,
            '1.33',
            'dict',
            0,
            {
                'MEDIAN_COVERAGE': [10.0, 20.0, 30.0, 40.0, 50.0],
                'FREEMIX': [0.001] * 5,
                'DUP': [5.0, 6.0, 7.0, 8.0, 9.0],
                'error_rate': [0.01] * 5,
            },
        ),
    ),
)


@pytest.fixture(autouse=True)
def _restore_config_paths():
    previous = os.environ.get('CPG_CONFIG_PATH', '')
    yield
    config.set_config_paths([p for p in previous.split(',') if p])


def _rendered() -> str:
    return emit_mod.render(SPEC, CACHE, spec_path='calibration/spec.genome.toml', generated='2026-08-11')


# --- the golden assertion: the parsed block, exactly --------------------------


def test_emitted_block_parses_to_exactly_the_expected_config():
    assert tomlio.loads(_rendered()) == {
        'qc_thresholds': {
            'genome': {
                'fail': {'min': {'MEDIAN_COVERAGE': 15}, 'max': {'FREEMIX': 0.04, 'DUP': 40}},
                'warn': {'min': {'MEDIAN_COVERAGE': 25}, 'max': {'FREEMIX': 0.02}},
                'relative': {'DUP': {'direction': 'max', 'k': 3.5, 'min_cohort': 50}},
            },
        },
    }


def test_ungated_metric_is_not_emitted():
    assert 'error_rate' not in _rendered()


def test_relative_metric_has_no_absolute_warn_entry():
    warn_max = tomlio.loads(_rendered())['qc_thresholds']['genome']['warn']['max']
    assert 'DUP' not in warn_max


def test_round_trips_through_production_load_thresholds(tmp_path):
    """The block must be loadable by the code that will actually enforce it."""
    path = tmp_path / 'emitted.toml'
    path.write_text(_rendered())
    config.set_config_paths([str(path)])
    assert check_multiqc.load_thresholds('genome') == {
        'min': {'MEDIAN_COVERAGE': {'fail': 15, 'warn': 25}},
        'max': {'FREEMIX': {'fail': 0.04, 'warn': 0.02}, 'DUP': {'fail': 40}},
    }


# --- structure and provenance --------------------------------------------------


def test_sections_appear_in_config_template_order():
    text = _rendered()
    order = [
        '[qc_thresholds.genome.fail.min]',
        '[qc_thresholds.genome.fail.max]',
        '[qc_thresholds.genome.warn.min]',
        '[qc_thresholds.genome.warn.max]',
        '[qc_thresholds.genome.relative.DUP]',
    ]
    positions = [text.index(section) for section in order]
    assert positions == sorted(positions)


def test_header_records_provenance_without_naming_datasets():
    """Cohort labels are real dataset names; this block gets committed."""
    text = _rendered()
    assert '1 cohorts' in text or '1 cohort' in text
    assert 'calibration/spec.genome.toml' in text
    assert 'dataset-a' not in text


def test_metric_comment_carries_rationale_and_generated_evidence():
    text = _rendered()
    assert '# Primary depth gate.' in text
    assert 'Cohort medians 30.0-30.0' in text
    assert '20%' in text  # 1 of 5 samples below the fail line of 15


def test_empty_section_is_omitted():
    spec = spec_mod.loads(
        'seq_type = "genome"\ncache = "c.json"\n'
        '[metrics.MEDIAN_COVERAGE]\ndirection = "min"\nunit = "x"\nfail = 15\nreviewed = true\n',
    )
    text = emit_mod.render(spec, CACHE, generated='2026-08-11')
    assert '[qc_thresholds.genome.fail.min]' in text
    assert 'warn' not in text.split('[qc_thresholds.genome.fail.min]')[1]


# --- the review gate -----------------------------------------------------------


def test_refuses_to_emit_an_unreviewed_gated_metric():
    unreviewed = SPEC.with_metric(replace_key='MEDIAN_COVERAGE', reviewed=False)
    with pytest.raises(EmitError, match=r"unreviewed.*\['MEDIAN_COVERAGE'\]"):
        emit_mod.render(unreviewed, CACHE, generated='2026-08-11')
```

- [ ] **Step 2: Run the test to verify it fails**

Run: `uv run --extra test pytest test/test_qc_calibration_emit.py -q`

Expected: FAIL — `ModuleNotFoundError: No module named 'align_genotype.qc_calibration.emit'`

- [ ] **Step 3: Write the measurement formatter**

Create `src/align_genotype/qc_calibration/report.py`:

```python
"""Human-readable formatting for the calibration commands."""


def fmt_measure(value: float, unit: str) -> str:
    """Render a metric value at a precision that reads sensibly for its unit."""
    if unit == 'frac':
        return f'{value:.3f}'
    if unit in ('x', '%'):
        return f'{value:.1f}'
    return f'{value:.2f}'
```

- [ ] **Step 4: Write the emit implementation**

Create `src/align_genotype/qc_calibration/emit.py`:

```python
"""Generate the [qc_thresholds.<seq_type>...] block for config_template.toml.

Prints; never writes. The tool doesn't get to clobber hand-maintained config, and an
operator eyeballing a diff before pasting is a cheap last line of defence.

Note that what this emits is committed and pushed, so it must never contain cohort
labels - those are real dataset names. The header cites the manifest path instead; the
manifest itself stays in the gitignored working directory.
"""

import textwrap
from dataclasses import dataclass

import numpy as np

from align_genotype.qc_calibration import stats, tomlio
from align_genotype.qc_calibration.cache import ValueCache
from align_genotype.qc_calibration.report import fmt_measure
from align_genotype.qc_calibration.spec import CalibrationSpec, MetricSpec

_WIDTH = 118

_RELATIVE_PREAMBLE = (
    'Cohort-relative (MAD / modified z-score) warn flagging - warn only. For metrics whose absolute level is '
    'strongly protocol- or kit-dependent, a fixed warn line either floods or never fires; instead we warn on '
    'samples that are statistical outliers within the current run. A sample is warned when its modified z-score '
    'exceeds `k` on the `direction` side: mz = 0.6745 * (value - median) / MAD; threshold = median +/- k*MAD/0.6745. '
    'Relative flags never fail, and are skipped when the cohort has fewer than `min_cohort` samples or zero MAD. '
    'The absolute fail gate above still applies; a sample caught there is not double-flagged here.'
)


class EmitError(RuntimeError):
    """The spec isn't ready to be turned into config."""


@dataclass(frozen=True)
class _Tier:
    severity: str
    direction: str


_TIERS = (_Tier('fail', 'min'), _Tier('fail', 'max'), _Tier('warn', 'min'), _Tier('warn', 'max'))


def _comment(text: str) -> list[str]:
    """Wrap `text` as TOML comment lines within the project's line length."""
    return [f'# {line}' for line in textwrap.wrap(text, width=_WIDTH - 2)]


def _evidence(cache: ValueCache, metric: MetricSpec, severity: str) -> str:
    """Cohort median range and observed flag rate for this metric at this tier."""
    medians: list[float] = []
    rates: list[float] = []
    for label in cache.labels:
        series = cache.series(label, metric.key)
        if not series.size:
            continue
        medians.append(float(np.median(series)))
        fail_rate, warn_rate = stats.flag_rates(series, metric.direction, fail=metric.fail, warn=metric.warn)
        rates.append(fail_rate if severity == 'fail' else warn_rate)
    if not medians:
        return ''
    return (
        f'Cohort medians {fmt_measure(min(medians), metric.unit)}-{fmt_measure(max(medians), metric.unit)}; '
        f'flags {min(rates):.0%}-{max(rates):.0%} per cohort at this tier.'
    )


def _require_reviewed(spec: CalibrationSpec) -> None:
    if unreviewed := sorted(m.key for m in spec.gated if not m.reviewed):
        raise EmitError(
            f'Refusing to emit: unreviewed gated metric(s) {unreviewed}. These were seeded by `suggest` and have '
            f'not been signed off - check them against the flagrates output, then set reviewed = true in the spec.',
        )


def render(
    spec: CalibrationSpec,
    cache: ValueCache,
    spec_path: str = '',
    manifest_path: str = '',
    generated: str = '',
) -> str:
    """Render the config block for `spec`, with evidence drawn from `cache`."""
    _require_reviewed(spec)
    seq_type = spec.seq_type
    n_cohorts = len(cache.labels)

    lines = [
        f'# {seq_type} QC thresholds, calibrated with `qc_calibrate` on {generated}',
        f'# against {n_cohorts} cohorts. Re-run the calibration when onboarding a new capture kit or protocol.',
    ]
    if spec_path:
        lines.append(f'# Spec: {spec_path}')
    if manifest_path:
        lines.append(f'# Manifest (cohort list): {manifest_path}')

    for tier in _TIERS:
        metrics = [
            m for m in spec.gated if m.direction == tier.direction and getattr(m, tier.severity) is not None
        ]
        if not metrics:
            continue
        lines += ['', f'[qc_thresholds.{seq_type}.{tier.severity}.{tier.direction}]']
        for metric in metrics:
            lines += _comment(f'{metric.rationale} {_evidence(cache, metric, tier.severity)}'.strip())
            lines.append(tomlio.fmt_kv(metric.key, getattr(metric, tier.severity), quote_key=True))

    if relative_metrics := [m for m in spec.gated if m.relative is not None]:
        lines += ['', *_comment(_RELATIVE_PREAMBLE)]
        for metric in relative_metrics:
            lines += _comment(f'{metric.rationale} {_evidence(cache, metric, "warn")}'.strip())
            lines += [
                f'[qc_thresholds.{seq_type}.relative.{metric.key}]',
                tomlio.fmt_kv('direction', metric.direction),
                tomlio.fmt_kv('k', metric.relative.k),
                tomlio.fmt_kv('min_cohort', metric.relative.min_cohort),
            ]

    return '\n'.join(lines) + '\n'
```

- [ ] **Step 5: Run the test to verify it passes**

Run: `uv run --extra test pytest test/test_qc_calibration_emit.py -q`

Expected: PASS (9 tests). If `test_metric_comment_carries_rationale_and_generated_evidence` fails on the exact `20%` substring, print the rendered block and check the arithmetic before changing the assertion — 1 of 5 samples below 15 is 20%, so a different number means `flag_rates` is being called wrongly, not that the test is too strict.

- [ ] **Step 6: Lint and commit**

```bash
uvx ruff format src/align_genotype/qc_calibration/emit.py src/align_genotype/qc_calibration/report.py test/test_qc_calibration_emit.py
uvx ruff check src/align_genotype/qc_calibration/emit.py src/align_genotype/qc_calibration/report.py test/test_qc_calibration_emit.py
git add src/align_genotype/qc_calibration/emit.py src/align_genotype/qc_calibration/report.py test/test_qc_calibration_emit.py
git commit -m "feat(qc_calibration): emit the qc_thresholds config block with generated rationale

Golden-tested by parsing the block and round-tripping it through the
production load_thresholds(). Refuses to emit an unreviewed gated metric,
and never names a cohort - the emitted block is committed config."
```

---

## Task 10: `dryrun` — end-to-end on one real report

Runs the genuine `check_multiqc.run()` against one cohort's report using the config `emit` would produce, and reports flag counts by `(metric, severity, method)`. This is what tells the operator "here is what will actually happen when this ships".

Two improvements over the original `dryrun_exome_check.py` / `dryrun_genome_check.py`:
- **One command, not two.** The seq_type comes from the spec, so there's no exome/genome copy to keep in sync.
- **No monkeypatching.** The thresholds are written to a temp TOML via `emit.render` and loaded through `set_config_paths`, so the dry run also proves the emitted block is what production reads. The old scripts patched `config_retrieve`, which meant they could pass against config that wouldn't actually load.

The report is parsed exactly once (by `check_multiqc.run` itself); the old scripts pre-parsed and injected the result to avoid a double parse, which is unnecessary once we're not inspecting the structure separately.

**Files:**
- Create: `src/align_genotype/qc_calibration/dryrun.py`
- Test: `test/test_qc_calibration_dryrun.py`

- [ ] **Step 1: Write the failing test**

Create `test/test_qc_calibration_dryrun.py`:

```python
"""End-to-end dry run of the production check against a small report."""

import json
import os

import pytest

from cpg_utils import config

from align_genotype.qc_calibration import dryrun as dryrun_mod
from align_genotype.qc_calibration import spec as spec_mod
from align_genotype.qc_calibration.cache import CohortValues, ValueCache
from align_genotype.qc_calibration.manifest import Cohort

SPEC = spec_mod.loads(
    'seq_type = "genome"\ncache = "c.json"\n'
    '\n[metrics.MEDIAN_COVERAGE]\ndirection = "min"\nunit = "x"\nfail = 15\nwarn = 25\nreviewed = true\n'
    'rationale = "Depth gate."\n'
    '\n[metrics.DUP]\ndirection = "max"\nunit = "%"\nfail = 40\nreviewed = true\nrationale = "Dup gate."\n'
    '[metrics.DUP.relative]\nk = 3.5\nmin_cohort = 5\n',
)

CACHE = ValueCache(
    seq_type='genome',
    generated='x',
    complete=True,
    metrics=('MEDIAN_COVERAGE', 'DUP'),
    cohorts=(
        CohortValues('dataset-a', 6, '1.33', 'dict', 0, {'MEDIAN_COVERAGE': [30.0] * 6, 'DUP': [5.0] * 6}),
    ),
)

# One sample below the fail line, one in the warn band, four healthy; DUP has one
# clear within-cohort outlier so the relative pass has something to find.
SECTIONS = {
    'picard': {
        'CPG1|S1': {'MEDIAN_COVERAGE': 10.0},
        'CPG2|S2': {'MEDIAN_COVERAGE': 20.0},
        'CPG3|S3': {'MEDIAN_COVERAGE': 30.0},
        'CPG4|S4': {'MEDIAN_COVERAGE': 31.0},
        'CPG5|S5': {'MEDIAN_COVERAGE': 32.0},
        'CPG6|S6': {'MEDIAN_COVERAGE': 33.0},
    },
    'samtools': {
        'CPG1|S1': {'DUP': 5.0},
        'CPG2|S2': {'DUP': 5.5},
        'CPG3|S3': {'DUP': 6.0},
        'CPG4|S4': {'DUP': 6.5},
        'CPG5|S5': {'DUP': 7.0},
        'CPG6|S6': {'DUP': 35.0},
    },
}


@pytest.fixture(autouse=True)
def _restore_config_paths():
    previous = os.environ.get('CPG_CONFIG_PATH', '')
    yield
    config.set_config_paths([p for p in previous.split(',') if p])


@pytest.fixture
def report(tmp_path):
    path = tmp_path / 'multiqc_data.json'
    path.write_text(json.dumps({'config_version': '1.33', 'report_general_stats_data': SECTIONS}))
    return Cohort('dataset-a', str(path))


def test_build_config_includes_sequencing_type_and_thresholds():
    text = dryrun_mod.build_config(SPEC, CACHE)
    assert 'sequencing_type = "genome"' in text
    assert '[qc_thresholds.genome.fail.min]' in text


def test_dryrun_reports_absolute_fail_and_warn(report, tmp_path):
    result = dryrun_mod.execute(SPEC, CACHE, report, output_dir=tmp_path)
    assert result.counts[('MEDIAN_COVERAGE', 'fail', 'absolute')] == 1  # the 10.0
    assert result.counts[('MEDIAN_COVERAGE', 'warn', 'absolute')] == 1  # the 20.0


def test_dryrun_exercises_the_relative_pass(report, tmp_path):
    result = dryrun_mod.execute(SPEC, CACHE, report, output_dir=tmp_path)
    assert result.counts[('DUP', 'warn', 'relative')] == 1  # the 35.0 outlier


def test_dryrun_counts_flagged_samples_and_writes_structured_output(report, tmp_path):
    result = dryrun_mod.execute(SPEC, CACHE, report, output_dir=tmp_path)
    assert result.n_samples_flagged == 3  # CPG1, CPG2, CPG6
    written = json.loads((tmp_path / 'dryrun_dataset-a.json').read_text())
    assert written['sequencing_type'] == 'genome'


def test_dryrun_records_timing_and_memory(report, tmp_path):
    result = dryrun_mod.execute(SPEC, CACHE, report, output_dir=tmp_path)
    assert result.seconds >= 0
    assert result.peak_rss_gb > 0


def test_dryrun_restores_config_paths(report, tmp_path):
    existing = tmp_path / 'existing.toml'
    existing.write_text('[workflow]\nsequencing_type = "exome"\n')
    config.set_config_paths([str(existing)])
    dryrun_mod.execute(SPEC, CACHE, report, output_dir=tmp_path)
    assert config.get_config_paths() == [str(existing)]
```

- [ ] **Step 2: Run the test to verify it fails**

Run: `uv run --extra test pytest test/test_qc_calibration_dryrun.py -q`

Expected: FAIL — `ModuleNotFoundError: No module named 'align_genotype.qc_calibration.dryrun'`

- [ ] **Step 3: Write the implementation**

Create `src/align_genotype/qc_calibration/dryrun.py`:

```python
"""Run the production QC check against one real report, using the emitted config.

This is the last check before pasting thresholds into config_template.toml: it answers
"what will this actually flag?" using the genuine code path, including the
cohort-relative pass.

Because the config comes from ``emit.render`` and is loaded through
``set_config_paths``, a dry run that passes also proves the emitted block parses and is
readable by ``load_thresholds`` - which the old monkeypatching scripts could not tell
you.
"""

import resource
import sys
import tempfile
import time
from collections import Counter
from contextlib import contextmanager
from dataclasses import dataclass
from pathlib import Path
from typing import Iterator

from cpg_utils import config

from align_genotype.qc_calibration import emit
from align_genotype.qc_calibration.cache import ValueCache
from align_genotype.qc_calibration.manifest import Cohort
from align_genotype.qc_calibration.spec import CalibrationSpec
from align_genotype.scripts import check_multiqc


@dataclass(frozen=True)
class DryRunResult:
    cohort: str
    n_samples_flagged: int
    counts: Counter
    seconds: float
    peak_rss_gb: float
    output_path: str


def peak_rss_gb() -> float:
    """Peak resident set size in GB (ru_maxrss is bytes on macOS, KiB on Linux)."""
    maxrss = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    return maxrss / (1024**3) if sys.platform == 'darwin' else maxrss / (1024**2)


def build_config(spec: CalibrationSpec, cache: ValueCache) -> str:
    """The full config a production run would need: seq type plus the emitted thresholds."""
    thresholds = emit.render(spec, cache, generated='dryrun')
    return f'[workflow]\nsequencing_type = "{spec.seq_type}"\n\n{thresholds}'


@contextmanager
def _config_from(text: str) -> Iterator[None]:
    previous = config.get_config_paths()
    with tempfile.TemporaryDirectory() as tmpdir:
        path = Path(tmpdir) / 'dryrun.toml'
        path.write_text(text)
        config.set_config_paths([str(path)])
        try:
            yield
        finally:
            config.set_config_paths(list(previous))


def execute(spec: CalibrationSpec, cache: ValueCache, cohort: Cohort, output_dir: str | Path) -> DryRunResult:
    """Run the real check against `cohort`'s report and summarise what it flagged."""
    output_path = Path(output_dir) / f'dryrun_{cohort.label}.json'

    with _config_from(build_config(spec, cache)):
        started = time.perf_counter()
        result = check_multiqc.run(
            multiqc_json_path=cohort.uri,
            html_url=None,
            dataset=cohort.label,
            title=f'{spec.seq_type} calibration dry run',
            send_to_slack=False,
            output_json_path=str(output_path),
        )
        elapsed = time.perf_counter() - started

    counts: Counter = Counter(
        (flag['flag'], flag['severity'], flag.get('method', 'absolute'))
        for flags in result['qc_flags'].values()
        for flag in flags
    )
    return DryRunResult(
        cohort=cohort.label,
        n_samples_flagged=result['n_samples_flagged'],
        counts=counts,
        seconds=elapsed,
        peak_rss_gb=peak_rss_gb(),
        output_path=str(output_path),
    )
```

- [ ] **Step 4: Run the test to verify it passes**

Run: `uv run --extra test pytest test/test_qc_calibration_dryrun.py -q`

Expected: PASS (6 tests)

- [ ] **Step 5: Lint and commit**

```bash
uvx ruff format src/align_genotype/qc_calibration/dryrun.py test/test_qc_calibration_dryrun.py
uvx ruff check src/align_genotype/qc_calibration/dryrun.py test/test_qc_calibration_dryrun.py
git add src/align_genotype/qc_calibration/dryrun.py test/test_qc_calibration_dryrun.py
git commit -m "feat(qc_calibration): seq_type-agnostic dry run through the real config path

One command replaces the exome and genome dry-run scripts. Loading the
emitted block via set_config_paths instead of monkeypatching means a passing
dry run also proves the block is readable by load_thresholds."
```

---

## Task 11: `discovery` — Metamist to manifest

Formalises `testing_scripts/fetch_multiqc_json_paths.py`: find the eligible datasets, take the most recent completed QC analysis for the requested sequencing type in each, and write a manifest.

> **Improvement on the design spec.** The spec said discovery would go untested because it's "a thin Metamist query". Passing the query function in as an argument makes the selection logic — dataset filtering, picking the latest analysis, skipping datasets with no report — fully unit-testable with a fake, and leaves only the ~4-line GraphQL adapter untested. Worth it: "pick the latest by timestamp" is exactly the kind of logic that silently picks the wrong one.

`metamist` is imported inside the adapter function, so nothing else in the tool — and no test — needs it installed.

**Files:**
- Create: `src/align_genotype/qc_calibration/discovery.py`
- Test: `test/test_qc_calibration_discovery.py`

- [ ] **Step 1: Write the failing test**

Create `test/test_qc_calibration_discovery.py`:

```python
"""Unit tests for Metamist cohort discovery, using a fake query function."""

import pytest

from align_genotype.qc_calibration import discovery as discovery_mod

PROJECTS = {
    'myProjects': [
        {'name': 'dataset-a', 'dataset': 'dataset-a', 'meta': {'is_seqr': True}},
        {'name': 'dataset-b', 'dataset': 'dataset-b', 'meta': {'is_seqr': True}},
        {'name': 'dataset-a-test', 'dataset': 'dataset-a-test', 'meta': {'is_seqr': True}},
        {'name': 'dataset-c-training', 'dataset': 'dataset-c-training', 'meta': {'is_seqr': True}},
        {'name': 'dataset-d', 'dataset': 'dataset-d', 'meta': {'is_seqr': False}},
    ],
}

ANALYSES = {
    'dataset-a': [
        {'id': 1, 'meta': {'sequencing_type': 'genome'}, 'output': 'gs://a/old.json', 'timestampCompleted': '2026-01-01T00:00:00'},
        {'id': 2, 'meta': {'sequencing_type': 'genome'}, 'output': 'gs://a/new.json', 'timestampCompleted': '2026-06-01T00:00:00'},
        {'id': 3, 'meta': {'sequencing_type': 'exome'}, 'output': 'gs://a/exome.json', 'timestampCompleted': '2026-07-01T00:00:00'},
    ],
    'dataset-b': [
        {'id': 4, 'meta': {'sequencing_type': 'genome'}, 'output': 'gs://b/g.json', 'timestampCompleted': '2026-03-01T00:00:00'},
    ],
}


def _fake_query(query_text: str, variables: dict | None = None) -> dict:
    if 'myProjects' in query_text:
        return PROJECTS
    return {'project': {'analyses': ANALYSES.get(variables['datasetName'], [])}}


def test_selects_eligible_datasets_only():
    manifest = discovery_mod.build_manifest('genome', query_fn=_fake_query, generated='2026-08-11')
    assert manifest.labels == ('dataset-a', 'dataset-b')


@pytest.mark.parametrize(
    ('project', 'eligible'),
    [
        ({'name': 'dataset-a', 'meta': {'is_seqr': True}}, True),
        ({'name': 'dataset-a-test', 'meta': {'is_seqr': True}}, False),
        ({'name': 'dataset-c-training', 'meta': {'is_seqr': True}}, False),
        ({'name': 'dataset-seqr-x', 'meta': {'is_seqr': True}}, False),
        ({'name': 'dataset-d', 'meta': {'is_seqr': False}}, False),
        ({'name': 'dataset-e', 'meta': {}}, False),
    ],
)
def test_dataset_eligibility(project, eligible):
    assert discovery_mod.is_eligible(project) is eligible


def test_picks_the_latest_analysis_for_the_requested_seq_type():
    manifest = discovery_mod.build_manifest('genome', query_fn=_fake_query, generated='2026-08-11')
    cohort = manifest.cohort('dataset-a')
    assert cohort.uri == 'gs://a/new.json'
    assert cohort.analysis_id == 2
    assert cohort.timestamp == '2026-06-01T00:00:00'


def test_other_seq_types_are_ignored():
    manifest = discovery_mod.build_manifest('exome', query_fn=_fake_query, generated='2026-08-11')
    assert manifest.labels == ('dataset-a',)
    assert manifest.cohort('dataset-a').uri == 'gs://a/exome.json'


def test_dataset_with_no_matching_analysis_is_omitted():
    manifest = discovery_mod.build_manifest('exome', query_fn=_fake_query, generated='2026-08-11')
    assert 'dataset-b' not in manifest.labels


def test_analysis_without_an_output_path_is_ignored():
    def query_fn(query_text, variables=None):
        if 'myProjects' in query_text:
            return {'myProjects': [{'name': 'dataset-a', 'dataset': 'dataset-a', 'meta': {'is_seqr': True}}]}
        return {
            'project': {
                'analyses': [
                    {'id': 1, 'meta': {'sequencing_type': 'genome'}, 'output': None, 'timestampCompleted': '2026-09-01T00:00:00'},
                    {'id': 2, 'meta': {'sequencing_type': 'genome'}, 'output': 'gs://a/g.json', 'timestampCompleted': '2026-01-01T00:00:00'},
                ],
            },
        }

    manifest = discovery_mod.build_manifest('genome', query_fn=query_fn, generated='x')
    assert manifest.cohort('dataset-a').analysis_id == 2


def test_dataset_label_that_is_not_a_bare_toml_key_is_skipped(caplog):
    def query_fn(query_text, variables=None):
        if 'myProjects' in query_text:
            return {'myProjects': [{'name': 'odd.name', 'dataset': 'odd.name', 'meta': {'is_seqr': True}}]}
        return {'project': {'analyses': [
            {'id': 1, 'meta': {'sequencing_type': 'genome'}, 'output': 'gs://x/g.json', 'timestampCompleted': '2026-01-01T00:00:00'},
        ]}}

    with pytest.raises(discovery_mod.DiscoveryError, match='no cohorts'):
        discovery_mod.build_manifest('genome', query_fn=query_fn, generated='x')
    assert 'odd.name' in caplog.text


def test_no_eligible_datasets_raises():
    def query_fn(query_text, variables=None):
        return {'myProjects': []}

    with pytest.raises(discovery_mod.DiscoveryError, match='no cohorts'):
        discovery_mod.build_manifest('genome', query_fn=query_fn, generated='x')


def test_manifest_round_trips_through_toml():
    from align_genotype.qc_calibration import manifest as manifest_mod

    original = discovery_mod.build_manifest('genome', query_fn=_fake_query, generated='2026-08-11')
    assert manifest_mod.loads(manifest_mod.dumps(original)) == original
```

- [ ] **Step 2: Run the test to verify it fails**

Run: `uv run --extra test pytest test/test_qc_calibration_discovery.py -q`

Expected: FAIL — `ModuleNotFoundError: No module named 'align_genotype.qc_calibration.discovery'`

- [ ] **Step 3: Write the implementation**

Create `src/align_genotype/qc_calibration/discovery.py`:

```python
"""Find the MultiQC reports to calibrate against, and write them to a manifest.

The only Metamist-aware module, and the import is inside the adapter so nothing else
in the tool depends on it. The selection logic takes the query function as an argument
so it can be tested with a fake - "pick the latest analysis" is exactly the sort of
rule that silently picks the wrong one.

Discovery is a separate step, not something the other commands do implicitly: the
manifest it writes is the durable record of which data a calibration used, and it is
meant to be reviewed and edited before anything is parsed.
"""

import logging
from datetime import datetime, timezone
from typing import Any, Callable

from align_genotype.qc_calibration import tomlio
from align_genotype.qc_calibration.manifest import Cohort, Manifest

QueryFn = Callable[[str, dict[str, Any] | None], dict[str, Any]]

DATASETS_QUERY = """
    query Datasets {
        myProjects {
            name
            dataset
            meta
        }
    }
"""

ANALYSES_QUERY = """
    query DatasetData($datasetName: String!) {
        project(name: $datasetName) {
            analyses(status: {eq: COMPLETED}, type: {eq: "qc"}) {
                id
                meta
                output
                timestampCompleted
            }
        }
    }
"""

# Substrings that mark a project as not a real production cohort.
EXCLUDED_NAME_TOKENS = ('test', 'training', 'seqr')


class DiscoveryError(RuntimeError):
    """Discovery found nothing usable."""


def default_query(query_text: str, variables: dict[str, Any] | None = None) -> dict[str, Any]:
    """Run a GraphQL query against Metamist. Imported lazily - only this path needs it."""
    from metamist.graphql import gql, query

    return query(gql(query_text), variables=variables or {})


def is_eligible(project: dict[str, Any]) -> bool:
    """Whether a project is a real production cohort worth calibrating against."""
    name = project.get('name', '')
    if not project.get('meta', {}).get('is_seqr', False):
        return False
    return not any(token in name for token in EXCLUDED_NAME_TOKENS)


def latest_analysis(analyses: list[dict[str, Any]], seq_type: str) -> dict[str, Any] | None:
    """The most recently completed QC analysis of `seq_type` that produced an output."""
    candidates = [
        analysis
        for analysis in analyses
        if analysis.get('meta', {}).get('sequencing_type') == seq_type and analysis.get('output')
    ]
    if not candidates:
        return None
    return max(candidates, key=lambda analysis: analysis['timestampCompleted'])


def build_manifest(seq_type: str, query_fn: QueryFn = default_query, generated: str | None = None) -> Manifest:
    """Build a manifest of the latest `seq_type` MultiQC report per eligible dataset."""
    projects = query_fn(DATASETS_QUERY, None)['myProjects']
    datasets = [project['dataset'] for project in projects if is_eligible(project)]

    cohorts: list[Cohort] = []
    for dataset in datasets:
        try:
            tomlio.require_bare_key(dataset, 'cohort label')
        except ValueError as exc:
            logging.warning(f'Skipping dataset {dataset!r}: {exc}')
            continue
        analyses = query_fn(ANALYSES_QUERY, {'datasetName': dataset})['project']['analyses']
        if (analysis := latest_analysis(analyses, seq_type)) is None:
            logging.info(f'{dataset}: no completed {seq_type} QC analysis with an output; skipping')
            continue
        cohorts.append(
            Cohort(
                label=dataset,
                uri=analysis['output'],
                analysis_id=analysis['id'],
                timestamp=analysis['timestampCompleted'],
            ),
        )

    if not cohorts:
        raise DiscoveryError(
            f'Discovery found no cohorts with a completed {seq_type} QC analysis. Check your Metamist access '
            f'and that the eligibility filter ({EXCLUDED_NAME_TOKENS} excluded, is_seqr required) suits your intent.',
        )
    return Manifest(
        seq_type=seq_type,
        generated=generated or datetime.now(tz=timezone.utc).isoformat(timespec='seconds'),
        cohorts=tuple(cohorts),
    )
```

- [ ] **Step 4: Run the test to verify it passes**

Run: `uv run --extra test pytest test/test_qc_calibration_discovery.py -q`

Expected: PASS (14 tests including parametrised cases)

- [ ] **Step 5: Lint and commit**

```bash
uvx ruff format src/align_genotype/qc_calibration/discovery.py test/test_qc_calibration_discovery.py
uvx ruff check src/align_genotype/qc_calibration/discovery.py test/test_qc_calibration_discovery.py
git add src/align_genotype/qc_calibration/discovery.py test/test_qc_calibration_discovery.py
git commit -m "feat(qc_calibration): Metamist cohort discovery behind an injectable query fn

Selection logic (eligibility, latest-analysis-wins) is unit-tested with a
fake; only the GraphQL adapter needs metamist, and it imports it lazily."
```

---

## Task 12: `report` — the stdout tables

Every reporting function returns a string rather than printing, so the output is testable and the CLI stays a thin shell around it.

**Files:**
- Modify: `src/align_genotype/qc_calibration/report.py` (currently just `fmt_measure`)
- Test: `test/test_qc_calibration_report.py`

- [ ] **Step 1: Write the failing test**

Create `test/test_qc_calibration_report.py`:

```python
"""Unit tests for the calibration report formatting."""

from align_genotype.qc_calibration import report, spec as spec_mod
from align_genotype.qc_calibration.cache import CohortValues, ValueCache
from align_genotype.qc_calibration.collect import CollectResult, SurveyRow
from align_genotype.qc_calibration.dryrun import DryRunResult
from align_genotype.qc_calibration.relative import CohortMad, MadEvaluation
from align_genotype.qc_calibration.stats import ChurnResult
from align_genotype.qc_calibration.suggest import Seeded

from collections import Counter

SPEC = spec_mod.loads(
    'seq_type = "genome"\ncache = "c.json"\n'
    '\n[metrics.MEDIAN_COVERAGE]\ndirection = "min"\nunit = "x"\nfail = 15\nwarn = 25\nreviewed = true\n'
    '\n[metrics.FREEMIX]\ndirection = "max"\nunit = "frac"\nfail = 0.04\nreviewed = true\n',
)

CACHE = ValueCache(
    seq_type='genome',
    generated='x',
    complete=True,
    metrics=('MEDIAN_COVERAGE', 'FREEMIX'),
    cohorts=(
        CohortValues(
            'dataset-a',
            5,
            '1.33',
            'dict',
            1,
            {'MEDIAN_COVERAGE': [10.0, 20.0, 30.0, 40.0, 50.0], 'FREEMIX': [0.001] * 5},
        ),
    ),
)

ROW = SurveyRow(
    label='dataset-a',
    uri='gs://x/multiqc_data.json',
    multiqc_version='1.33',
    shape='dict',
    n_samples=5,
    section_sizes={'picard': 5, 'verifybamid': 5},
    where={'MEDIAN_COVERAGE': ('picard',), 'FREEMIX': ('verifybamid',)},
    section_keys={'picard': ('MEDIAN_COVERAGE',), 'verifybamid': ('FREEMIX',)},
    n_dropped=1,
)


def test_table_aligns_columns():
    text = report.table(['a', 'longheader'], [['1', '2'], ['333', '4']])
    lines = text.splitlines()
    assert lines[0].startswith('a    longheader')
    assert len(lines) == 4  # header, rule, two rows


def test_survey_report_shows_version_shape_and_drops():
    text = report.survey_report(CollectResult(CACHE, (ROW,), {}, ()), SPEC)
    assert '1.33' in text
    assert 'dict' in text
    assert 'dataset-a' in text
    assert '1' in text  # the dropped non-numeric value


def test_survey_report_marks_missing_gated_metric_and_dumps_candidate_keys():
    result = CollectResult(CACHE, (ROW,), {'dataset-a': ('FREEMIX',)}, ())
    text = report.survey_report(result, SPEC)
    assert 'MISSING' in text
    assert 'candidate keys' in text.lower()


def test_survey_report_lists_unreadable_cohorts():
    result = CollectResult(CACHE, (ROW,), {}, (('dataset-b', 'no such file'),))
    text = report.survey_report(result, SPEC)
    assert 'dataset-b' in text
    assert 'no such file' in text


def test_distributions_report_has_a_row_per_cohort_and_a_percentile_per_column():
    text = report.distributions_report(CACHE, SPEC)
    assert 'MEDIAN_COVERAGE' in text
    assert 'p50' in text
    assert '30.0' in text  # the median of 10..50


def test_flagrates_report_shows_fail_and_warn_per_cohort():
    text = report.flagrates_report(CACHE, SPEC)
    assert 'MEDIAN_COVERAGE' in text
    assert '20% / 20%' in text  # 1 of 5 fails (<15), 1 of 5 warns (<25 but not failing)


def test_flagrates_report_marks_metrics_outside_the_guardrail():
    text = report.flagrates_report(CACHE, SPEC)
    assert 'REVIEW' in text  # a 20% fail rate is well over the healthy-cohort limit


def test_mad_report_shows_the_verdict_and_its_reason():
    evaluation = MadEvaluation(
        metric='DUP',
        direction='max',
        cohorts=(CohortMad('dataset-a', 100, 7.0, 1.0, 12.2, 4, None),),
        homogeneous=(('dataset-a', ChurnResult(12.2, 12.4, 60, 3, 3, 0)),),
        heterogeneous=(),
    )
    text = report.mad_report(evaluation)
    assert 'RECOMMEND' in text
    assert 'dataset-a' in text
    assert '12.2' in text


def test_mad_report_explains_a_rejection():
    evaluation = MadEvaluation(
        metric='PCT_OFF_BAIT',
        direction='max',
        cohorts=(CohortMad('dataset-a', 100, 7.0, 1.0, 12.2, 4, None),),
        homogeneous=(('dataset-a', ChurnResult(12.2, 30.0, 100, 5, 30, 25)),),
        heterogeneous=(),
    )
    text = report.mad_report(evaluation)
    assert 'REJECT' in text
    assert 'churn' in text


def test_mad_report_notes_skipped_cohorts():
    evaluation = MadEvaluation(
        metric='DUP',
        direction='max',
        cohorts=(CohortMad('dataset-a', 3, 7.0, 1.0, None, 0, 'cohort 3 < min_cohort 50'),),
        homogeneous=(),
        heterogeneous=(),
    )
    assert 'min_cohort' in report.mad_report(evaluation)


def test_suggest_summary_flags_everything_as_unreviewed():
    text = report.suggest_summary((Seeded('MEDIAN_COVERAGE', 11, 15, 'Seeded from p1 (10.99).'),))
    assert 'MEDIAN_COVERAGE' in text
    assert 'reviewed = false' in text


def test_dryrun_summary_groups_by_metric_severity_and_method():
    result = DryRunResult(
        cohort='dataset-a',
        n_samples_flagged=3,
        counts=Counter({('MEDIAN_COVERAGE', 'fail', 'absolute'): 1, ('DUP', 'warn', 'relative'): 2}),
        seconds=1.5,
        peak_rss_gb=0.4,
        output_path='/tmp/out.json',
    )
    text = report.dryrun_summary(result)
    assert 'MEDIAN_COVERAGE' in text
    assert 'relative' in text
    assert '3' in text
```

- [ ] **Step 2: Run the test to verify it fails**

Run: `uv run --extra test pytest test/test_qc_calibration_report.py -q`

Expected: FAIL — `AttributeError: module 'align_genotype.qc_calibration.report' has no attribute 'table'`

- [ ] **Step 3: Write the implementation**

Replace the whole of `src/align_genotype/qc_calibration/report.py` with:

```python
"""Human-readable formatting for the calibration commands.

Every function returns a string rather than printing, so output is testable and the
CLI stays a thin shell.

The result types are imported under TYPE_CHECKING only. `emit` needs `fmt_measure` from
here and `dryrun` imports `emit`, so a runtime import of `DryRunResult` would close the
cycle report -> dryrun -> emit -> report.
"""

from __future__ import annotations

from collections.abc import Sequence
from typing import TYPE_CHECKING

import numpy as np

from align_genotype.qc_calibration import stats

if TYPE_CHECKING:
    from align_genotype.qc_calibration.cache import ValueCache
    from align_genotype.qc_calibration.collect import CollectResult
    from align_genotype.qc_calibration.dryrun import DryRunResult
    from align_genotype.qc_calibration.relative import MadEvaluation
    from align_genotype.qc_calibration.spec import CalibrationSpec
    from align_genotype.qc_calibration.suggest import Seeded


def fmt_measure(value: float, unit: str) -> str:
    """Render a metric value at a precision that reads sensibly for its unit."""
    if unit == 'frac':
        return f'{value:.3f}'
    if unit in ('x', '%'):
        return f'{value:.1f}'
    return f'{value:.2f}'


def table(headers: Sequence[str], rows: Sequence[Sequence[str]]) -> str:
    """Render a fixed-width table with a rule under the header."""
    cells = [[str(c) for c in row] for row in [headers, *rows]]
    widths = [max(len(row[i]) for row in cells) for i in range(len(headers))]

    def line(row: Sequence[str]) -> str:
        return '  '.join(str(c).ljust(w) for c, w in zip(row, widths, strict=True)).rstrip()

    return '\n'.join([line(headers), line(['-' * w for w in widths]), *(line(row) for row in cells[1:])])


def survey_report(result: CollectResult, spec: CalibrationSpec) -> str:
    """Provenance per cohort, where each metric lives, and anything that went wrong."""
    parts = [
        '=== survey: report provenance ===',
        table(
            ['cohort', 'multiqc', 'shape', 'samples', 'sections', 'non-numeric dropped'],
            [
                [
                    row.label,
                    row.multiqc_version,
                    row.shape,
                    str(row.n_samples),
                    ', '.join(f'{name}={n}' for name, n in row.section_sizes.items()),
                    str(row.n_dropped),
                ]
                for row in result.rows
            ],
        ),
        '',
        '=== survey: metric presence (gated metrics marked *) ===',
        table(
            ['metric', *[row.label for row in result.rows]],
            [
                [
                    f'{metric.key}*' if metric.gated else metric.key,
                    *[', '.join(row.where.get(metric.key, ())) or 'MISSING' for row in result.rows],
                ]
                for metric in spec.metrics
            ],
        ),
    ]

    if result.missing_gated:
        parts += ['', '!! GATED METRICS MISSING - these would silently check nothing:']
        for label, metrics in result.missing_gated.items():
            parts.append(f'   {label}: {list(metrics)}')
        parts.append('   Candidate keys present in each section, for mapping a MultiQC rename:')
        for row in result.rows:
            if row.label not in result.missing_gated:
                continue
            for section, keys in row.section_keys.items():
                parts.append(f'   [{row.label}/{section}] candidate keys: {list(keys)}')

    if result.failures:
        parts += ['', '!! COHORTS THAT COULD NOT BE READ:']
        parts += [f'   {label}: {message}' for label, message in result.failures]

    return '\n'.join(parts)


def distributions_report(cache: ValueCache, spec: CalibrationSpec) -> str:
    """A percentile table per metric, one row per cohort."""
    parts = []
    for metric in spec.metrics:
        sense = 'higher=better' if metric.direction == 'min' else 'lower=better'
        parts.append(f'\n=== {metric.key}  ({sense}, unit={metric.unit}) ===')
        rows = []
        for label in cache.labels:
            series = cache.series(label, metric.key)
            if not series.size:
                rows.append([label, '0', *['-'] * len(stats.PERCENTILES)])
                continue
            values = stats.percentiles(series)
            rows.append([label, str(series.size), *[fmt_measure(v, metric.unit) for v in values]])
        parts.append(table(['cohort', 'n', *[f'p{p}' for p in stats.PERCENTILES]], rows))
    return '\n'.join(parts).lstrip()


def flagrates_report(cache: ValueCache, spec: CalibrationSpec) -> str:
    """What each candidate threshold would flag, per cohort, as fail% / warn%."""
    rows = []
    for metric in spec.gated:
        tiers = f'{metric.fail} / {metric.warn}'
        cells, review = [], False
        for label in cache.labels:
            series = cache.series(label, metric.key)
            fail_rate, warn_rate = stats.flag_rates(series, metric.direction, metric.fail, metric.warn)
            if np.isnan(fail_rate):
                cells.append('n/a')
                continue
            review = review or stats.needs_review(fail_rate, warn_rate)
            cells.append(f'{fail_rate:.0%} / {warn_rate:.0%}')
        rows.append([metric.key, tiers, *cells, 'REVIEW' if review else ''])
    header = ['metric', 'fail / warn', *cache.labels, '']
    return '\n'.join([
        'Flag rates per cohort (fail% / warn%; warn excludes samples already failing).',
        f'REVIEW marks a metric exceeding {stats.FAIL_RATE_LIMIT:.0%} fail or {stats.WARN_RATE_LIMIT:.0%} warn',
        'in some cohort - a prompt to look, not a rejection.',
        '',
        table(header, rows),
    ])


def mad_report(evaluation: MadEvaluation) -> str:
    """Per-cohort relative thresholds and warn rates, churn, and the adoption verdict."""
    parts = [f'=== {evaluation.metric}  (cohort-relative, direction={evaluation.direction}) ===']

    rows = []
    for cohort in evaluation.cohorts:
        if cohort.skipped:
            rows.append([cohort.label, str(cohort.n), '-', '-', '-', '-', cohort.skipped])
            continue
        rows.append([
            cohort.label,
            str(cohort.n),
            f'{cohort.median:.4f}',
            f'{cohort.mad_raw:.4f}',
            f'{cohort.threshold:.4f}',
            f'{cohort.warn_rate:.1%}',
            '',
        ])
    parts.append(table(['cohort', 'n', 'median', 'MAD', 'threshold', 'warn', 'note'], rows))

    if evaluation.homogeneous:
        parts += ['', 'Churn - homogeneous growth (first 60% of a cohort -> all of it):']
        parts.append(
            table(
                ['cohort', 'threshold', 'flagged', 'flips'],
                [
                    [
                        label,
                        f'{result.threshold_before:.4f} -> {result.threshold_after:.4f}',
                        f'{result.flagged_before} -> {result.flagged_after}',
                        f'{result.flips} ({result.flip_rate:.1%})',
                    ]
                    for label, result in evaluation.homogeneous
                ],
            ),
        )

    if evaluation.heterogeneous:
        worst = sorted(evaluation.heterogeneous, key=lambda item: -item[2].flip_rate)[:5]
        parts += ['', 'Churn - heterogeneous growth (cohort A joined by cohort B), 5 worst pairs:']
        parts.append(
            table(
                ['cohort A', 'joined by', 'threshold', 'flips'],
                [
                    [
                        a,
                        b,
                        f'{result.threshold_before:.4f} -> {result.threshold_after:.4f}',
                        f'{result.flips} ({result.flip_rate:.1%})',
                    ]
                    for a, b, result in worst
                ],
            ),
        )

    reason = evaluation.verdict_reason or (
        f'peak warn {evaluation.max_warn_rate:.1%}, peak churn {evaluation.max_churn:.1%} - both within the bar'
    )
    parts += ['', f'VERDICT: {evaluation.verdict} - {reason}']
    parts.append('(Adopting this is your call: set [metrics.<KEY>.relative] in the spec.)')
    return '\n'.join(parts)


def suggest_summary(seeded: Sequence[Seeded]) -> str:
    """What `suggest` wrote, and the reminder that none of it is signed off."""
    if not seeded:
        return 'Nothing seeded: every gated metric is already reviewed, or has no data in any cohort.'
    rows = [[s.key, str(s.fail), '-' if s.warn is None else str(s.warn), s.evidence] for s in seeded]
    return '\n'.join([
        table(['metric', 'fail', 'warn', 'evidence'], rows),
        '',
        f'{len(seeded)} metric(s) seeded with reviewed = false. These are starting points drawn from the',
        'distributions, not decisions - check them against `flagrates`, then set reviewed = true per metric.',
        '`emit-config` will refuse to run until you do.',
    ])


def dryrun_summary(result: DryRunResult) -> str:
    """Flag counts from a real run of the production check."""
    rows = [
        [metric, severity, method, str(count)]
        for (metric, severity, method), count in sorted(result.counts.items(), key=lambda kv: -kv[1])
    ]
    return '\n'.join([
        f'Dry run on {result.cohort}: {result.n_samples_flagged} samples flagged '
        f'in {result.seconds:.1f}s (peak RSS {result.peak_rss_gb:.2f} GB)',
        '',
        table(['metric', 'severity', 'method', 'flags'], rows),
        '',
        f'Structured output: {result.output_path}',
    ])
```

- [ ] **Step 4: Run the test to verify it passes**

Run: `uv run --extra test pytest test/test_qc_calibration_report.py -q`

Expected: PASS (13 tests)

- [ ] **Step 5: Lint and commit**

```bash
uvx ruff format src/align_genotype/qc_calibration/report.py test/test_qc_calibration_report.py
uvx ruff check src/align_genotype/qc_calibration/report.py test/test_qc_calibration_report.py
git add src/align_genotype/qc_calibration/report.py test/test_qc_calibration_report.py
git commit -m "feat(qc_calibration): stdout tables for every calibration command"
```

---

## Task 13: `cli` — wire the eight subcommands together

Thin shell: parse options, call the module, print the string, choose an exit code. All the logic already exists and is tested.

**Files:**
- Create: `src/align_genotype/qc_calibration/cli.py`
- Modify: `pyproject.toml` (add the console script)
- Modify: `.gitignore` (ignore the operator's working directory)
- Test: `test/test_qc_calibration_cli.py`

- [ ] **Step 1: Write the failing test**

Create `test/test_qc_calibration_cli.py`:

```python
"""CLI behaviour, especially exit codes on the failure paths."""

import json
import os

import pytest
from click.testing import CliRunner

from cpg_utils import config

from align_genotype.qc_calibration import cache as cache_mod
from align_genotype.qc_calibration import spec as spec_mod
from align_genotype.qc_calibration.cli import main

SECTIONS = {
    'picard': {f'CPG{i}|S{i}': {'MEDIAN_COVERAGE': float(10 + i)} for i in range(8)},
    'verifybamid': {f'CPG{i}|S{i}': {'FREEMIX': 0.001} for i in range(8)},
}

SPEC_TEXT = (
    'seq_type = "genome"\ncache = "{cache}"\n'
    '\n[metrics.MEDIAN_COVERAGE]\ndirection = "min"\nunit = "x"\nfail = 15\nwarn = 25\nreviewed = true\n'
    'rationale = "Depth gate."\n'
    '\n[metrics.FREEMIX]\ndirection = "max"\nunit = "frac"\nfail = 0.04\nreviewed = true\n'
    'rationale = "Contamination gate."\n'
)


@pytest.fixture(autouse=True)
def _restore_config_paths():
    previous = os.environ.get('CPG_CONFIG_PATH', '')
    yield
    config.set_config_paths([p for p in previous.split(',') if p])


@pytest.fixture
def workspace(tmp_path):
    """A spec, a manifest and one small report, all under tmp_path."""
    report = tmp_path / 'dataset-a.json'
    report.write_text(json.dumps({'config_version': '1.33', 'report_general_stats_data': SECTIONS}))

    manifest = tmp_path / 'manifest.toml'
    manifest.write_text(f'seq_type = "genome"\ngenerated = "x"\n\n[cohorts.dataset-a]\nuri = "{report}"\n')

    spec = tmp_path / 'spec.toml'
    spec.write_text(SPEC_TEXT.format(cache=tmp_path / 'values.json'))
    return {'spec': str(spec), 'manifest': str(manifest), 'cache': str(tmp_path / 'values.json'), 'dir': tmp_path}


def _run(*args):
    return CliRunner().invoke(main, list(args))


def test_help_lists_every_subcommand():
    result = _run('--help')
    assert result.exit_code == 0
    for command in ('discover', 'collect', 'distributions', 'suggest', 'flagrates', 'mad', 'emit-config', 'dryrun'):
        assert command in result.output


def test_collect_writes_a_complete_cache(workspace):
    result = _run('collect', '--spec', workspace['spec'], '--manifest', workspace['manifest'])
    assert result.exit_code == 0, result.output
    assert cache_mod.load(workspace['cache']).complete is True


def test_collect_exits_nonzero_when_a_gated_metric_is_missing(workspace):
    spec = workspace['dir'] / 'spec.toml'
    spec.write_text(spec.read_text() + '\n[metrics.NOT_A_REAL_KEY]\ndirection = "min"\nfail = 1\n')
    result = _run('collect', '--spec', str(spec), '--manifest', workspace['manifest'])
    assert result.exit_code == 1
    assert 'NOT_A_REAL_KEY' in result.output
    assert cache_mod.load(workspace['cache']).complete is False  # written, but unusable


def test_downstream_command_refuses_an_incomplete_cache(workspace):
    spec = workspace['dir'] / 'spec.toml'
    spec.write_text(spec.read_text() + '\n[metrics.NOT_A_REAL_KEY]\ndirection = "min"\nfail = 1\n')
    _run('collect', '--spec', str(spec), '--manifest', workspace['manifest'])
    result = _run('flagrates', '--spec', str(spec))
    assert result.exit_code == 1
    assert 'incomplete' in result.output


def test_flagrates_reports_rates(workspace):
    _run('collect', '--spec', workspace['spec'], '--manifest', workspace['manifest'])
    result = _run('flagrates', '--spec', workspace['spec'])
    assert result.exit_code == 0, result.output
    assert 'MEDIAN_COVERAGE' in result.output


def test_distributions_reports_percentiles(workspace):
    _run('collect', '--spec', workspace['spec'], '--manifest', workspace['manifest'])
    result = _run('distributions', '--spec', workspace['spec'])
    assert result.exit_code == 0, result.output
    assert 'p50' in result.output


def test_suggest_rewrites_the_spec_in_place(workspace):
    _run('collect', '--spec', workspace['spec'], '--manifest', workspace['manifest'])
    spec_path = workspace['dir'] / 'unreviewed.toml'
    spec_path.write_text(SPEC_TEXT.format(cache=workspace['cache']).replace('reviewed = true', 'reviewed = false'))
    result = _run('suggest', '--spec', str(spec_path))
    assert result.exit_code == 0, result.output
    assert spec_mod.load(spec_path).metric('MEDIAN_COVERAGE').reviewed is False
    assert 'reviewed = false' in result.output


def test_emit_config_prints_the_block(workspace):
    _run('collect', '--spec', workspace['spec'], '--manifest', workspace['manifest'])
    result = _run('emit-config', '--spec', workspace['spec'])
    assert result.exit_code == 0, result.output
    assert '[qc_thresholds.genome.fail.min]' in result.output


def test_emit_config_writes_to_a_file_when_asked(workspace):
    _run('collect', '--spec', workspace['spec'], '--manifest', workspace['manifest'])
    out = workspace['dir'] / 'block.toml'
    result = _run('emit-config', '--spec', workspace['spec'], '--output', str(out))
    assert result.exit_code == 0, result.output
    assert '[qc_thresholds.genome.fail.min]' in out.read_text()


def test_emit_config_refuses_unreviewed_metrics(workspace):
    _run('collect', '--spec', workspace['spec'], '--manifest', workspace['manifest'])
    spec_path = workspace['dir'] / 'unreviewed.toml'
    spec_path.write_text(SPEC_TEXT.format(cache=workspace['cache']).replace('reviewed = true', 'reviewed = false'))
    result = _run('emit-config', '--spec', str(spec_path))
    assert result.exit_code == 1
    assert 'unreviewed' in result.output


def test_mad_says_so_when_no_metric_has_a_relative_block(workspace):
    _run('collect', '--spec', workspace['spec'], '--manifest', workspace['manifest'])
    result = _run('mad', '--spec', workspace['spec'])
    assert result.exit_code == 0
    assert 'no metrics' in result.output.lower()


def test_dryrun_runs_the_real_check(workspace):
    _run('collect', '--spec', workspace['spec'], '--manifest', workspace['manifest'])
    result = _run(
        'dryrun',
        '--spec', workspace['spec'],
        '--manifest', workspace['manifest'],
        '--cohort', 'dataset-a',
        '--output-dir', str(workspace['dir']),
    )
    assert result.exit_code == 0, result.output
    assert 'MEDIAN_COVERAGE' in result.output


def test_a_malformed_spec_is_a_clean_error_not_a_traceback(workspace):
    bad = workspace['dir'] / 'bad.toml'
    bad.write_text('seq_type = "genome"\ncache = "c.json"\n[metrics.M]\ndirection = "sideways"\nfail = 1\n')
    result = _run('flagrates', '--spec', str(bad))
    assert result.exit_code == 1
    assert 'direction' in result.output
    assert 'Traceback' not in result.output
```

- [ ] **Step 2: Run the test to verify it fails**

Run: `uv run --extra test pytest test/test_qc_calibration_cli.py -q`

Expected: FAIL — `ModuleNotFoundError: No module named 'align_genotype.qc_calibration.cli'`

- [ ] **Step 3: Write the implementation**

Create `src/align_genotype/qc_calibration/cli.py`:

```python
"""`qc_calibrate` - derive QC thresholds from a set of MultiQC reports.

Typical run:

    qc_calibrate discover     --seq-type genome --output calibration/manifest.genome.toml
    qc_calibrate collect      --spec calibration/spec.genome.toml --manifest calibration/manifest.genome.toml
    qc_calibrate distributions --spec calibration/spec.genome.toml
    qc_calibrate suggest      --spec calibration/spec.genome.toml
    qc_calibrate flagrates    --spec calibration/spec.genome.toml      # iterate here
    qc_calibrate mad          --spec calibration/spec.genome.toml
    qc_calibrate emit-config  --spec calibration/spec.genome.toml
    qc_calibrate dryrun       --spec ... --manifest ... --cohort <label>

See README.md in this package for the full workflow.
"""

import functools
import logging
from datetime import datetime, timezone
from typing import Any, Callable

import click

from align_genotype.qc_calibration import cache as cache_mod
from align_genotype.qc_calibration import collect as collect_mod
from align_genotype.qc_calibration import discovery, dryrun, emit, relative, report, suggest
from align_genotype.qc_calibration import manifest as manifest_mod
from align_genotype.qc_calibration import spec as spec_mod

logging.basicConfig(level=logging.INFO, format='%(levelname)s: %(message)s')

CALIBRATION_ERRORS = (
    spec_mod.SpecError,
    manifest_mod.ManifestError,
    cache_mod.CacheError,
    collect_mod.CollectError,
    emit.EmitError,
    discovery.DiscoveryError,
)

_spec_option = click.option('--spec', 'spec_path', required=True, help='Path to the calibration spec TOML.')
_manifest_option = click.option('--manifest', 'manifest_path', required=True, help='Path to the cohort manifest TOML.')


def clean_errors(fn: Callable[..., Any]) -> Callable[..., Any]:
    """Turn expected calibration failures into a one-line error, not a traceback."""

    @functools.wraps(fn)
    def wrapper(*args: Any, **kwargs: Any) -> Any:
        try:
            return fn(*args, **kwargs)
        except CALIBRATION_ERRORS as exc:
            raise click.ClickException(str(exc)) from exc

    return wrapper


def _load_spec_and_cache(spec_path: str) -> tuple[spec_mod.CalibrationSpec, cache_mod.ValueCache]:
    spec = spec_mod.load(spec_path)
    values = cache_mod.load(spec.cache)
    cache_mod.require_usable(values, spec)
    return spec, values


def _now() -> str:
    return datetime.now(tz=timezone.utc).isoformat(timespec='seconds')


@click.group()
def main() -> None:
    """Derive warn/fail and cohort-relative QC thresholds from MultiQC reports."""


@main.command()
@click.option('--seq-type', required=True, help="Sequencing type, e.g. 'genome' or 'exome'.")
@click.option('--output', 'output_path', required=True, help='Where to write the manifest TOML.')
@clean_errors
def discover(seq_type: str, output_path: str) -> None:
    """Query Metamist for the latest MultiQC report per eligible dataset."""
    built = discovery.build_manifest(seq_type, generated=_now())
    manifest_mod.save(built, output_path)
    click.echo(f'Wrote {len(built.cohorts)} cohorts to {output_path}. Review it before collecting.')


@main.command()
@_spec_option
@_manifest_option
@clean_errors
def collect(spec_path: str, manifest_path: str) -> None:
    """Survey and extract every cohort - the one pass over the large reports."""
    spec = spec_mod.load(spec_path)
    result = collect_mod.collect_all(manifest_mod.load(manifest_path), spec, generated=_now())
    cache_mod.save(result.cache, spec.cache)
    click.echo(report.survey_report(result, spec))
    click.echo(f'\nCache written to {spec.cache} (complete={result.cache.complete}).')
    if not result.ok:
        raise click.ClickException('Collection incomplete - see the survey output above. Nothing downstream will run.')


@main.command()
@_spec_option
@clean_errors
def distributions(spec_path: str) -> None:
    """Percentile tables per metric, per cohort."""
    spec, values = _load_spec_and_cache(spec_path)
    click.echo(report.distributions_report(values, spec))


@main.command('suggest')
@_spec_option
@click.option('--output', 'output_path', default=None, help='Where to write the updated spec (default: in place).')
@clean_errors
def suggest_cmd(spec_path: str, output_path: str | None) -> None:
    """Seed candidate thresholds from the distributions, marked unreviewed."""
    spec, values = _load_spec_and_cache(spec_path)
    updated, seeded = suggest.seed(values, spec)
    spec_mod.save(updated, output_path or spec_path)
    click.echo(report.suggest_summary(seeded))


@main.command()
@_spec_option
@clean_errors
def flagrates(spec_path: str) -> None:
    """What the spec's candidate thresholds would flag, per cohort."""
    spec, values = _load_spec_and_cache(spec_path)
    click.echo(report.flagrates_report(values, spec))


@main.command()
@_spec_option
@click.option('--metric', 'metric_key', default=None, help='Evaluate one metric only.')
@clean_errors
def mad(spec_path: str, metric_key: str | None) -> None:
    """Cohort-relative (MAD) spread, warn rates and cohort-growth churn."""
    spec, values = _load_spec_and_cache(spec_path)
    candidates = [m for m in spec.gated if m.relative is not None and metric_key in (None, m.key)]
    if not candidates:
        click.echo(
            'There are no metrics with a [metrics.<KEY>.relative] block to evaluate. To assess a candidate, add '
            'one (it needs an absolute `fail` too) and re-run - the verdict tells you whether to keep it.',
        )
        return
    for metric in candidates:
        click.echo(report.mad_report(relative.evaluate(values, metric, spec.seq_type)))
        click.echo('')


@main.command('emit-config')
@_spec_option
@click.option('--manifest', 'manifest_path', default='', help='Manifest path, cited in the generated header.')
@click.option('--output', 'output_path', default=None, help='Write to a file instead of stdout.')
@clean_errors
def emit_config(spec_path: str, manifest_path: str, output_path: str | None) -> None:
    """Print the [qc_thresholds.<seq_type>...] block to paste into config_template.toml."""
    spec, values = _load_spec_and_cache(spec_path)
    block = emit.render(
        spec,
        values,
        spec_path=spec_path,
        manifest_path=manifest_path,
        generated=datetime.now(tz=timezone.utc).date().isoformat(),
    )
    if output_path:
        with open(output_path, 'w') as f:
            f.write(block)
        click.echo(f'Wrote the config block to {output_path}.')
    else:
        click.echo(block)


@main.command('dryrun')
@_spec_option
@_manifest_option
@click.option('--cohort', 'cohort_label', required=True, help='Which cohort in the manifest to run against.')
@click.option('--output-dir', default='.', help='Where to write the structured QC flags JSON.')
@clean_errors
def dryrun_cmd(spec_path: str, manifest_path: str, cohort_label: str, output_dir: str) -> None:
    """Run the real check_multiqc against one cohort using the emitted config."""
    spec, values = _load_spec_and_cache(spec_path)
    cohort = manifest_mod.load(manifest_path).cohort(cohort_label)
    click.echo(report.dryrun_summary(dryrun.execute(spec, values, cohort, output_dir)))


if __name__ == '__main__':
    main()
```

Note the explicit command names on `suggest`, `dryrun` and `emit-config`: Click derives a command's name from the function name, and `suggest_cmd` / `dryrun_cmd` would otherwise register as `suggest-cmd` / `dryrun-cmd`. The function names carry the `_cmd` suffix to avoid shadowing the imported `suggest` and `dryrun` modules.

- [ ] **Step 4: Run the test to verify it passes**

Run: `uv run --extra test pytest test/test_qc_calibration_cli.py -q`

Expected: PASS (13 tests)

- [ ] **Step 5: Register the console script**

In `pyproject.toml`, under `[project.scripts]`, add below the existing two entries:

```toml
# the QC threshold calibration tool - operator-facing, not used by the pipeline
qc_calibrate = 'align_genotype.qc_calibration.cli:main'
```

- [ ] **Step 6: Ignore the operator's working directory**

Append to `.gitignore`:

```gitignore
# qc_calibrate working artifacts: manifests and specs name real datasets, and the
# value caches are derived data. None of it belongs in git.
calibration/
```

- [ ] **Step 7: Verify the installed entry point works**

Run: `uv run --extra test qc_calibrate --help`

Expected: the group help, listing `collect`, `discover`, `distributions`, `dryrun`, `emit-config`, `flagrates`, `mad`, `suggest`.

- [ ] **Step 8: Lint and commit**

```bash
uvx ruff format src/align_genotype/qc_calibration/cli.py test/test_qc_calibration_cli.py
uvx ruff check src/align_genotype/qc_calibration/cli.py test/test_qc_calibration_cli.py
uv run --extra test pytest test/ -q
git add src/align_genotype/qc_calibration/cli.py test/test_qc_calibration_cli.py pyproject.toml .gitignore
git commit -m "feat(qc_calibration): qc_calibrate CLI with the eight workflow subcommands"
```

---

## Task 14: Operator README

The guide someone reads when onboarding a new capture kit eighteen months from now. It has to carry the judgement that isn't in the code.

**Files:**
- Create: `src/align_genotype/qc_calibration/README.md`

- [ ] **Step 1: Write the README**

Create `src/align_genotype/qc_calibration/README.md`:

````markdown
# `qc_calibrate` — QC threshold calibration

Derives warn/fail and cohort-relative QC thresholds for a sequencing type from a set of
real MultiQC reports, and emits the `[qc_thresholds.<seq_type>...]` block for
`src/align_genotype/config_template.toml`.

Run it when onboarding a new capture kit or sequencing protocol, or when refreshing
thresholds against a newer set of cohorts.

## The three artifacts

Everything lives in a gitignored `calibration/` directory:

| File | What it is |
|---|---|
| `manifest.<seq_type>.toml` | which cohorts, and the exact report URI for each. The reproducibility record. |
| `<seq_type>_values.json` | the value cache — small per-metric numbers, produced by one parse of each report. |
| `spec.<seq_type>.toml` | which metrics you gate, in which direction, at what value. The thing you edit. |

**Never commit these.** Manifests and specs name real datasets.

## The workflow

```bash
# 1. Find the cohorts. Review the manifest before going further - drop anything
#    that isn't a representative production run.
qc_calibrate discover --seq-type genome --output calibration/manifest.genome.toml

# 2. Write calibration/spec.genome.toml listing every candidate metric with its
#    direction and unit (see "Writing a spec" below). Then parse the reports - this
#    is the slow step, and the only one that touches the large files.
qc_calibrate collect --spec calibration/spec.genome.toml --manifest calibration/manifest.genome.toml

# 3. Look at the distributions. For 'min' metrics the bad samples are the low tail;
#    for 'max' metrics the high tail.
qc_calibrate distributions --spec calibration/spec.genome.toml

# 4. Optional: seed a first draft of the thresholds.
qc_calibrate suggest --spec calibration/spec.genome.toml

# 5. Iterate. Edit thresholds in the spec, re-run - it reads the cache, so it's instant.
qc_calibrate flagrates --spec calibration/spec.genome.toml

# 6. For any metric whose normal level shifts by cohort, assess a relative warn tier.
qc_calibrate mad --spec calibration/spec.genome.toml

# 7. Set reviewed = true per metric once you're satisfied, then emit.
qc_calibrate emit-config --spec calibration/spec.genome.toml --manifest calibration/manifest.genome.toml

# 8. Paste the block into config_template.toml, then sanity-check end to end.
qc_calibrate dryrun --spec calibration/spec.genome.toml \
    --manifest calibration/manifest.genome.toml --cohort <label> --output-dir calibration/
```

## Writing a spec

```toml
seq_type = "genome"
cache = "calibration/genome_values.json"

[metrics.MEDIAN_COVERAGE]
direction = "min"    # min = higher is better (flag below); max = lower is better (flag above)
unit = "x"           # x | frac | %  - display precision only
gated = true         # enforced, and fatal if a cohort is missing it
fail = 15
warn = 25
reviewed = true      # emit-config refuses while this is false
rationale = "Primary depth gate. Cohort medians ~32-37x; p1 ~15-28x."
```

Rules the loader enforces:

- A gated metric needs at least one of `fail`, `warn`, `relative`.
- A `relative` block needs an absolute `fail` behind it — cohort-relative flagging is
  warn-only, and there must always be a hard stop.
- A `relative` block forbids an absolute `warn` — the relative tier *is* the warn tier.
- An un-gated metric carries no thresholds. Keep rejected candidates in the spec with
  `gated = false` and a `rationale` saying why, so the next person doesn't re-litigate it.

## Metric keys differ by sequencing type

- **Exome** uses Picard `CollectHsMetrics`: `MEAN_TARGET_COVERAGE`,
  `PCT_TARGET_BASES_20X/50X`, `FOLD_80_BASE_PENALTY`, `ZERO_CVG_TARGETS_PCT`,
  `PCT_SELECTED_BASES`, `PCT_OFF_BAIT`, `AT/GC_DROPOUT`.
- **Genome** uses Picard `CollectWgsMetrics`: `MEDIAN_COVERAGE`, `MEAN_COVERAGE`,
  `SD/MAD_COVERAGE`, `PCT_1X`…`PCT_100X`, `PCT_EXC_*`, `HET_SNP_SENSITIVITY`,
  `GENOME_TERRITORY`. There are **no** `*_TARGET_*` / `FOLD_80` / `ZERO_CVG` keys.
- **Shared** (samtools): `reads_mapped_percent`, `reads_duplicated_percent`,
  `reads_properly_paired_percent`, `reads_MQ0_percent`, `error_rate`. Contamination
  (verifybamid): `FREEMIX`.

`PCT_PF_READS_ALIGNED` is **not** in `report_general_stats_data` — only in
`report_saved_raw_data` — so it can never be gated. Use `reads_mapped_percent`. The old
genome reads-mapped gate was configured on it and checked nothing for as long as it
existed; the `collect` survey is what catches that class of bug, which is why a missing
gated metric is fatal rather than a warning.

## Choosing thresholds

- **A healthy cohort should flag ≈ 0% fail and single-digit % warn.** `fail` means "do
  not analyse without a decision"; `warn` means "a human should look".
- **Preserve the lab's intent for hard gates** unless the data clearly contradicts it.
  You can tighten `warn` freely. When genome `reads_duplicated_percent` fail was relaxed
  from 25 to 40, it was because 25 would have hard-failed a quarter to a third of
  legitimately higher-duplication preps — that's the bar for overriding lab intent.
- **Don't hard-fail on metrics that track ancestry, biology or chemistry**
  (`error_rate`, `HET_SNP_SENSITIVITY`). Warn, or go cohort-relative.
- **`suggest` produces starting points, not answers.** Everything it writes is
  `reviewed = false` and `emit-config` won't emit it. That gate is only worth anything
  if you actually check the number before flipping it.

## When to adopt a cohort-relative tier

Only where a metric's normal level genuinely shifts by cohort or protocol — bimodal or
wide-ranging medians across your cohort set — **and** the flag set stays stable as the
cohort grows. `mad` gives a RECOMMEND/REJECT verdict against peak warn ≤ 10% and peak
churn ≤ 2%, but the decision is yours; the spec's `relative` block is what adopts it.

Churn is the important half. Every flip is a spurious "updated" flag in the database
caused by nothing but the cohort composition changing.

Adopted so far: exome `ZERO_CVG_TARGETS_PCT`, genome `reads_duplicated_percent`.
Rejected: exome `PCT_SELECTED_BASES` and `PCT_OFF_BAIT` (churn up to 24.5%).

## Memory

The reports are large — tens to ~500 MB each. `collect` parses one at a time and
releases it, and it is the only command that touches them. If you find yourself wanting
to load several at once, use the cache instead; that's what it's for.
````

- [ ] **Step 2: Check for dataset names before committing**

Derive the names to search for from your local (gitignored) manifest, so no real dataset
name is ever written into a tracked file — not even into a grep pattern:

```bash
LABELS=$(grep -oE '^\[cohorts\.[^]]+\]' calibration/manifest.genome.toml | sed 's/\[cohorts\.//;s/\]//' | paste -sd'|' -)
grep -nEi "$LABELS" src/align_genotype/qc_calibration/README.md || echo 'clean'
```

Expected: `clean`.

- [ ] **Step 3: Commit**

```bash
git add src/align_genotype/qc_calibration/README.md
git commit -m "docs(qc_calibration): operator guide for the calibration workflow"
```

---

## Task 15: Reproduce the known results, then finish the branch

The correctness check that matters: the formalised tool must reproduce the analysis it replaces. The exome and genome numbers in `docs/qc_thresholds_tiers_plan.md` came from the old scripts; the new tool should agree with them.

The two existing caches use the old flat format (`{label: {metric: [values]}}`), so they need converting to the new one first. This is a one-off throwaway — put it in gitignored `testing_scripts/`, don't commit it.

**Files:**
- Create (gitignored, not committed): `testing_scripts/convert_old_cache.py`
- Create (gitignored, not committed): `calibration/spec.genome.toml`, `calibration/spec.exome.toml`

- [ ] **Step 1: Convert the old genome cache**

Create `testing_scripts/convert_old_cache.py`:

```python
"""One-off: convert an old flat metric cache to the qc_calibration ValueCache format.

    python testing_scripts/convert_old_cache.py \
        testing_scripts/testing_data/wgs_metric_values.json calibration/genome_values.json genome
"""

import json
import sys

from align_genotype.qc_calibration.cache import CohortValues, ValueCache, save

old_path, new_path, seq_type = sys.argv[1], sys.argv[2], sys.argv[3]
with open(old_path) as f:
    old = json.load(f)

metrics = sorted({metric for by_metric in old.values() for metric in by_metric})
cache = ValueCache(
    seq_type=seq_type,
    generated='converted-from-legacy-cache',
    complete=True,
    metrics=tuple(metrics),
    cohorts=tuple(
        CohortValues(
            label=label,
            n_samples=max((len(values) for values in by_metric.values()), default=0),
            multiqc_version='unknown',
            shape='unknown',
            n_dropped=0,
            values={m: [v for v in vals if v is not None and v == v] for m, vals in by_metric.items()},
        )
        for label, by_metric in old.items()
    ),
)
save(cache, new_path)
print(f'Wrote {new_path}: {len(cache.cohorts)} cohorts, {len(metrics)} metrics')
```

Run:

```bash
mkdir -p calibration
uv run --extra test python3 testing_scripts/convert_old_cache.py \
    testing_scripts/testing_data/wgs_metric_values.json calibration/genome_values.json genome
uv run --extra test python3 testing_scripts/convert_old_cache.py \
    testing_scripts/testing_data/wes_metric_values.json calibration/exome_values.json exome
```

- [ ] **Step 2: Write specs matching the shipped config**

Create `calibration/spec.genome.toml` transcribing the current
`[qc_thresholds.genome.*]` from `src/align_genotype/config_template.toml`:
`MEDIAN_COVERAGE` (min, x, fail 15, warn 25), `PCT_20X` (min, frac, fail 0.75, warn 0.85),
`reads_mapped_percent` (min, %, fail 80, warn 97), `reads_properly_paired_percent`
(min, %, warn 92), `FREEMIX` (max, frac, fail 0.04, warn 0.02), and
`reads_duplicated_percent` (max, %, fail 40, relative k 3.5 / min_cohort 50). Set
`reviewed = true` on all of them and `cache = "calibration/genome_values.json"`.

Do the same for `calibration/spec.exome.toml` from `[qc_thresholds.exome.*]`:
`MEAN_TARGET_COVERAGE` (min, x, fail 30, warn 50), `PCT_TARGET_BASES_20X` (min, frac,
fail 0.80, warn 0.90), `reads_mapped_percent` (min, %, fail 90, warn 97), `FREEMIX`
(max, frac, fail 0.04, warn 0.02), `FOLD_80_BASE_PENALTY` (max, ratio → use `frac`,
fail 3.5, warn 2.5), `reads_duplicated_percent` (max, %, fail 50, warn 35), and
`ZERO_CVG_TARGETS_PCT` (max, frac, fail 0.10, relative k 3.5 / min_cohort 50).

- [ ] **Step 3: Check the distributions and flag rates reproduce the documented numbers**

```bash
uv run --extra test qc_calibrate distributions --spec calibration/spec.genome.toml
uv run --extra test qc_calibrate flagrates --spec calibration/spec.genome.toml
uv run --extra test qc_calibrate distributions --spec calibration/spec.exome.toml
uv run --extra test qc_calibrate flagrates --spec calibration/spec.exome.toml
```

Check against `docs/qc_thresholds_tiers_plan.md` and the rationale comments in
`config_template.toml`. Specifically, the genome run should show:

- `MEDIAN_COVERAGE` cohort medians in the 32–37x range, fail ~0–1% per cohort.
- `PCT_20X` cohort medians ~0.94–0.95, warn tail ~0–8%.
- `reads_mapped_percent` median ~99%, with one visibly messier cohort failing ~5%.
- `FREEMIX` p99 ≤ 0.008 across all cohorts, 0% fail.
- `reads_duplicated_percent` cohort medians spanning ~7.2–18.0%, 0% fail at 40.

If a number disagrees, **investigate before adjusting anything** — a mismatch means the
tool extracts or aggregates differently from the scripts it replaces, which is exactly
what this check exists to catch. Note which cohorts these numbers refer to only in your
local notes, never in a committed file.

- [ ] **Step 4: Check the MAD verdict reproduces the adoption decisions**

```bash
uv run --extra test qc_calibrate mad --spec calibration/spec.genome.toml
uv run --extra test qc_calibrate mad --spec calibration/spec.exome.toml
```

Expected: `reads_duplicated_percent` RECOMMEND with per-cohort warn 0–4.2% and peak churn
around 1% homogeneous / 2.1% heterogeneous; `ZERO_CVG_TARGETS_PCT` RECOMMEND with warn
0–8% and churn under 1%.

- [ ] **Step 5: Check the emitted block matches what already ships**

```bash
uv run --extra test qc_calibrate emit-config --spec calibration/spec.genome.toml
```

Compare against `[qc_thresholds.genome.*]` in `config_template.toml`. The threshold
*values* must match exactly. The comments will differ — the shipped ones are
hand-written and the generated ones cite the same evidence in a uniform format. Don't
update `config_template.toml` in this PR; the point here is only to confirm agreement.

- [ ] **Step 6: Full verification**

```bash
uv run --extra test pytest test/ -q
uvx ruff check src/align_genotype/ test/
uvx ruff format --check src/align_genotype/ test/
uv run --extra test python3 -c "import tomllib; tomllib.load(open('src/align_genotype/config_template.toml','rb')); print('config OK')"
```

Expected: all tests pass; ruff clean apart from the two baseline `PLR0917` warnings in
`check_multiqc.py`.

- [ ] **Step 7: Confirm nothing committed names a dataset**

```bash
git diff mad-relative-flagging...HEAD --name-only

# Search the diff for every cohort label in the local manifests, without hardcoding
# any of them into a tracked file.
LABELS=$(cat calibration/manifest.*.toml | grep -oE '^\[cohorts\.[^]]+\]' | sed 's/\[cohorts\.//;s/\]//' | sort -u | paste -sd'|' -)
git diff mad-relative-flagging...HEAD | grep -nEi "$LABELS" || echo 'clean'

git status --short  # calibration/ and testing_scripts/ must not appear as tracked
```

Expected: `clean`, and no `calibration/` entries.

- [ ] **Step 8: Push and open the stacked PR**

```bash
git push -u origin qc-calibration-workflow
gh pr create --base mad-relative-flagging --title "feat: formalise QC threshold calibration into a reusable tool" --body "$(cat <<'EOF'
Replaces the five ad-hoc `testing_scripts/` calibration scripts with a committed,
tested `qc_calibrate` CLI. Stacked on #74 (which is stacked on #73).

## What this adds
`src/align_genotype/qc_calibration/` — eight subcommands over three durable artifacts
(a Metamist-derived cohort manifest, a small value cache, and a calibration spec):
`discover`, `collect`, `distributions`, `suggest`, `flagrates`, `mad`, `emit-config`,
`dryrun`. The output is the ready-to-paste `[qc_thresholds.<seq_type>...]` block.

## Production change
`check_multiqc.run()` read `report_general_stats_data` and called `.items()` on it, but
MultiQC v1.14 stores that field as a positional list — so a v1.14 report raised
`AttributeError` and the QC check never ran. Adds a shared `normalise_sections()`, used
by both the checker and the calibration tool, and raises rather than silently reporting
a clean check on an unreadable report.

## Design decisions
- The mandatory key survey is folded into `collect`, so it costs no extra parse and a
  missing gated metric is fatal rather than skippable. This is the guard against the
  `PCT_PF_READS_ALIGNED` class of bug.
- `mad` calls the production `relative_flags()` directly, so there is no prototype
  implementation of the modified z-score to diverge from what ships.
- `dryrun` loads the emitted config through `set_config_paths` instead of
  monkeypatching, so it also proves the emitted block is readable by `load_thresholds`.
- `suggest` seeds thresholds from percentiles but marks them `reviewed = false`, and
  `emit-config` refuses to emit an unreviewed metric.

## Test plan
- [x] Unit tests per module, including both MultiQC section shapes
- [x] Golden test: emitted block parses and round-trips through `load_thresholds`
- [x] Full suite green, ruff clean
- [ ] Reproduces the exome and genome numbers in `docs/qc_thresholds_tiers_plan.md`
      from the existing caches (run locally — caches are gitignored)

🤖 Generated with [Claude Code](https://claude.com/claude-code)
EOF
)"
```

---

## Self-review checklist

Before declaring the plan complete, the implementing engineer should confirm:

- [ ] Every module listed in the File structure table exists and is imported by something.
- [ ] `uv run --extra test pytest test/ -q` is green, and coverage of `src/align_genotype/qc_calibration/` is at least 80% (`uv run --extra test pytest test/ --cov=src/align_genotype/qc_calibration --cov-report=term-missing`).
- [ ] No committed file contains a real dataset name.
- [ ] `calibration/` is gitignored and untracked.
- [ ] `qc_calibrate --help` lists all eight subcommands with the names in Task 13's test.
- [ ] `config_template.toml` is **unchanged** by this branch — the tool emits config, it doesn't apply it.

