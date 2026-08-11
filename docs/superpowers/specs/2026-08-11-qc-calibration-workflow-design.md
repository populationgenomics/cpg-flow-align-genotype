# Design: `qc_calibration` — a reusable QC threshold calibration workflow

**Date:** 2026-08-11
**Status:** approved, ready for implementation planning
**Supersedes the ad-hoc process in:** `docs/qc_calibration_workflow_handover.md`

## Problem

Deriving warn/fail (and cohort-relative MAD) QC thresholds for a sequencing type is
currently a manual process run through five throwaway scripts in `testing_scripts/`
(gitignored). Each script is edited in place between runs — the exome calibration and
the genome calibration mutated the same files, so neither run is reproducible and the
next operator inherits nothing but a handover document.

This design replaces those scripts with a committed, tested tool. The tool is run when
onboarding a new capture kit or sequencing protocol, or when refreshing thresholds for
an existing one. Its ultimate output is the `[qc_thresholds.<seq_type>...]` block of
`src/align_genotype/config_template.toml`.

**Out of scope:** changing any threshold values. This is tooling only; the exome and
genome numbers already shipped on branches `QC_flags_report` (PR #73) and
`mad-relative-flagging` (PR #74) stand as-is and serve as the regression target.

## Enabling fix: `check_multiqc` does not normalise the v1.14 section shape

The core requirement is *calibration ≡ enforcement* — the tool must extract metric
values exactly the way the production checker does, so the two can never disagree.
Realising that surfaced a real defect.

`check_multiqc.run()` reads `sections = d['report_general_stats_data']` and immediately
calls `sections.items()` (`check_multiqc.py:280,285`); `warn_unmatched_metrics` and
`_gather_metric_values` likewise assume a dict. MultiQC v1.33 stores general stats as a
dict keyed by section name, but v1.14 stores a positional **list** — the shape present in
the real WES dataset set. A v1.14 JSON raises `AttributeError` in production. There is no
list-shape test in `test_check_multiqc.py`.

The two calibration scripts each carried their own private normaliser (`iter_sections`,
`_section_dicts`) precisely because of this, which is how the gap went unnoticed.

**Fix (step 0 of implementation):** extract a shared normaliser into `check_multiqc` and
call it from `run()`:

```python
def normalise_sections(raw: Any) -> dict[str, dict[str, dict[str, Any]]]:
    """Normalise report_general_stats_data to {section: {sample: {metric: value}}}.

    MultiQC >=1.33 keys sections by name (dict); v1.14 uses a positional list.
    """
```

The calibration tool imports it rather than re-implementing it. This is the only change
to shipped production behaviour in this work, and it is a strict widening — dict inputs
behave identically.

## Architecture

A subpackage at `src/align_genotype/qc_calibration/`, exposed as the console script
`qc_calibrate`. Three durable, human-readable, diffable artifacts:

| Artifact | Produced by | Purpose |
|---|---|---|
| `manifest.<seq_type>.toml` | `discover`, hand-editable | which cohorts, which URIs — the reproducibility record |
| `<seq_type>_values.json` | `collect` | small per-metric value cache |
| `spec.<seq_type>.toml` | human, seeded by `suggest` | metrics, directions, candidate thresholds, relative settings |

The large MultiQC JSONs (tens to ~500 MB each) are parsed **exactly once**, by `collect`,
one file at a time. Every tuning step reads only the cache and completes in under a
second, so threshold iteration stays fast — the two-phase property that made the manual
process workable.

### Data flow

```
metamist ──discover──► manifest.toml
                            │
                            ▼
                        collect ──► survey report (fatal on MISSING gated metric)
                            │   └──► values cache
                            ▼
        ┌───────────┬───────────┬──────────┬──────────┐
    distributions  suggest   flagrates    mad      emit-config ──► TOML block (stdout)
                      │          │         │            │
                      └──────────┴─────────┴────────────┘
                              spec.toml (read/written)

    manifest + spec ──dryrun──► real check_multiqc.run() on one cohort
```

## Subcommands

Eight commands. Two deliberate consolidations relative to the five original scripts.

| Command | Input | Output |
|---|---|---|
| `discover` | Metamist | manifest.toml |
| `collect` | manifest | survey table + values cache |
| `distributions` | cache | percentile tables (p1…p99) per cohort |
| `suggest` | cache + spec | candidate thresholds seeded into spec (`reviewed = false`) |
| `flagrates` | cache + spec | per-cohort fail% / warn% per metric |
| `mad` | cache + spec | spread, warn rate, cohort-growth churn, verdict |
| `emit-config` | cache + spec | `[qc_thresholds.<seq_type>...]` block on stdout |
| `dryrun` | manifest + spec | real `check_multiqc.run()` flag counts on one cohort |

### Consolidation 1 — `survey` folds into `collect`

The survey step is mandatory: it is the guard against gating a key MultiQC renamed or
never surfaces (the `PCT_PF_READS_ALIGNED` bug, which made the old genome reads-mapped
gate silently inert). Surveying and extracting both require a full parse of every JSON,
so running them as separate commands doubles the only expensive operation in the
workflow.

`collect` therefore emits the survey table *and* the value cache in one pass. A metric
marked `gated = true` in the spec that is absent from any cohort is a **hard error**:
non-zero exit, and the cache is written with `"complete": false` so downstream commands
refuse to run. This makes the mandatory step structurally unskippable rather than merely
documented as such.

### Consolidation 2 — `verify_relative_6cohorts.py` folds into `mad`

The manual process had a MAD prototype and a separate script that verified the prototype
against the production code path. Keeping both invites divergence. In this tool, `mad`
**is** the production path:

- warn rates come from `check_multiqc.relative_flags()`, driven with a MultiQC-shaped
  cohort constructed from the cache;
- churn uses `check_multiqc.robust_threshold()`, including its 4-dp rounding.

There is no second implementation of the modified z-score anywhere in the tool.

### `dryrun` drops the monkeypatching

The original `dryrun_*` scripts monkeypatched `config.config_retrieve`. Instead, `dryrun`
writes the spec's thresholds to a temporary TOML, points `CPG_CONFIG_PATH` at it, and
calls `check_multiqc.run()` unmodified. This exercises `load_thresholds` for real, so the
dry run additionally proves the emitted config block is loadable — something the
monkeypatched version never checked.

It keeps the two memory optimisations from the originals: the JSON is loaded once and
injected (avoiding a double parse), and Slack is disabled.

The two seq_type-specific scripts collapse into one command; seq_type comes from the spec.

## Module layout

Each module has one job and stays small (target ≤200 lines, hard ceiling 400).

```
src/align_genotype/qc_calibration/
  __init__.py
  cli.py           click group and options only; no logic
  manifest.py      Cohort dataclass; manifest read/write
  discovery.py     Metamist query → cohorts (lazy import; only metamist-aware module)
  spec.py          MetricSpec / CalibrationSpec; load, validate
  cache.py         value cache I/O; staleness check against the spec's metric list
  collect.py       parse one JSON → survey row + extracted values
  stats.py         percentiles, flag rates, churn simulation
  relative.py      MAD evaluation via check_multiqc.relative_flags / robust_threshold
  emit.py          TOML block + generated rationale comments
  dryrun.py        temp config + check_multiqc.run()
  report.py        table formatting shared by the reporting commands
  README.md        operator guide: the workflow, step by step
```

`collect.py` imports `normalise_sections` and `_gather_metric_values` from
`check_multiqc`. No section normalisation or non-numeric coercion is re-implemented
anywhere in the subpackage.

`discovery.py` is the only module that imports `metamist`, and it does so inside the
function. Nothing else in the tool — and no test — needs Metamist available.

## Data formats

### Manifest

Generated by `discover`, then hand-editable: drop a bad cohort, pin an older analysis,
or point at a local file for testing.

```toml
seq_type = "genome"
generated = "2026-08-11T13:40:00"

[cohorts.dataset-a]
uri = "gs://cpg-dataset-a-main/qc/.../multiqc_data.json"
analysis_id = 84213
timestamp = "2026-06-02T04:11:09"

[cohorts.local-test]
uri = "file:///path/to/dataset-b_multiqc_data.json"
```

`uri` is resolved with `cpg_utils.to_path`, which already handles `gs://` and local
paths — so GCS, local files and anything else cpg-utils supports work without new code.
Cohort labels come from the retrieval layer (Metamist dataset name), never from parsing
filenames.

### Value cache

```json
{
  "seq_type": "genome",
  "generated": "2026-08-11T13:52:00",
  "complete": true,
  "metrics": ["MEDIAN_COVERAGE", "PCT_20X", "..."],
  "cohorts": {
    "dataset-a": {
      "n_samples": 841,
      "multiqc_version": "1.33",
      "shape": "dict",
      "n_dropped": 3,
      "values": {"MEDIAN_COVERAGE": [34.1, 29.7], "PCT_20X": [0.94, 0.91]}
    }
  }
}
```

**Staleness rule:** the cache is usable if the spec's metric set is a *subset* of
`cache["metrics"]`. Narrowing the metric list does not force a re-collect; adding a
metric does, and the error names the missing metrics and the command to run.

### Calibration spec

The single source of truth for what gets gated and at what value. Replaces editing
Python constants in place — the central pain point of the manual process.

```toml
seq_type = "genome"
cache = "calibration/genome_values.json"

[metrics.MEDIAN_COVERAGE]
direction = "min"          # min = higher is better (flag below); max = lower is better
unit = "x"                 # x | frac | % — controls display and rounding only
gated = true               # survey treats absence as fatal; included in emitted config
fail = 15
warn = 25
reviewed = true            # `suggest` writes false; emit-config refuses to emit false
rationale = "Primary depth gate. Cohort medians ~32-37x; p1 ~15-28x."

[metrics.reads_duplicated_percent]
direction = "max"
unit = "%"
gated = true
fail = 40
reviewed = true
rationale = "Strongly library-prep dependent; a fixed warn line over-flags high-dup preps."
# Warn tier is cohort-relative, so no absolute `warn` key.
[metrics.reads_duplicated_percent.relative]
k = 3.5
min_cohort = 50

[metrics.error_rate]
direction = "max"
unit = "frac"
gated = false              # surveyed and profiled, deliberately not gated
rationale = "Varies with ancestry and chemistry; not a defensible hard gate."
```

Three properties this buys:

- **`reviewed`** is the gate that stops a `suggest`-seeded number reaching production
  config unexamined.
- **Un-gated metrics stay in the spec**, so the next operator sees what was considered
  and rejected, and why — the "log what you rejected" guardrail, made structural.
- **`relative` without `warn`** encodes the rule that a MAD metric keeps an absolute
  `fail` hard gate and takes its warn tier cohort-relatively. Validated in `spec.py`.

## Command behaviour

### `suggest`

Seeds a first draft so the operator tunes rather than starts blank. For each gated metric
without a threshold:

- `min` metrics: `fail` ← min across cohorts of p1; `warn` ← min across cohorts of p5.
- `max` metrics: `fail` ← max across cohorts of p99; `warn` ← max across cohorts of p95.

Rounded by unit (`x` → integer, `%` → integer, `frac` → 2 dp). Every seeded metric is
written with `reviewed = false` and a `rationale` recording the percentiles it came from.
`suggest` never overwrites a threshold that is already `reviewed = true`.

These are starting points, not answers. The judgement calls the handover documents —
preserving the lab's intent for hard gates, not hard-failing on metrics that vary with
ancestry or chemistry, relaxing genome `reads_duplicated_percent` from 25 to 40 — are the
operator's, and the `reviewed` flag is where they sign off.

### `flagrates`

Per metric, per cohort: fail% and warn%, with warn computed excluding samples already
failing (matching production's fail-before-warn evaluation). Flags any metric for review
where a cohort exceeds ~2% fail or ~10% warn, against the guardrail that a healthy cohort
should sit near 0% fail and single-digit % warn. It flags for attention; it does not
reject — "healthy cohort" is the operator's judgement, not a computable property.

### `mad`

For each metric with a `relative` block, per cohort: n, median, raw MAD, derived
threshold, and warn rate from the production `relative_flags()` path. Then churn:

- **homogeneous growth** — first 60% of a cohort → 100%, per cohort;
- **heterogeneous growth** — cohort A + cohort B for all pairs (45 pairs at 10 cohorts;
  trivial on cached values), reporting the worst case.

Churn is the fraction of the initial cohort's samples whose flag status flips purely
because the cohort grew — the thing that would generate spurious "updated" flags in the
database.

Verdict per metric: **RECOMMEND** when max warn ≤ 10% and max churn ≤ 2%, **REJECT**
otherwise, always printed with the driving number. The operator decides; the spec's
`relative` block is what actually adopts it. This is the bar that admitted exome
`ZERO_CVG_TARGETS_PCT` and genome `reads_duplicated_percent` and rejected
`PCT_SELECTED_BASES` / `PCT_OFF_BAIT` (churn up to 24.5%).

### `emit-config`

Prints the block to stdout (or `--output` to a file) for the operator to paste. It never
writes to `config_template.toml` — the tool does not get to clobber hand-maintained
config.

Section order matches the existing file: `fail.min`, `fail.max`, `warn.min`, `warn.max`,
`relative.<METRIC>`. Each metric carries a comment built from its spec `rationale` plus
an auto-generated evidence line citing the cohort median range, the relevant percentile
range and the observed flag-rate range — matching the hand-written style already in the
config. A header records the date, the cohort count and labels, and the manifest and spec
paths the block was derived from.

Refuses to emit (non-zero exit) if any gated metric has `reviewed = false`.

## Error handling

| Condition | Behaviour |
|---|---|
| Gated metric absent from any cohort | Fatal. Names the cohort; dumps candidate section keys so a rename can be mapped. Cache marked `complete: false`. |
| `report_general_stats_data` absent or unparseable | Loud per-cohort error, never a silent skip. |
| Non-numeric value (Picard `'?'`) | Dropped and counted, per production behaviour; count surfaced in the survey. |
| NaN | Filtered before percentiles and MAD. |
| Spec metric not in cache | Refuse to report; name the metrics and tell the operator to re-`collect`. |
| Unreachable manifest URI | Record the failure, continue other cohorts, exit non-zero at the end with the full list. |
| `reviewed = false` on a gated metric | `emit-config` exits non-zero. |
| Cohort below `min_cohort` | Skipped for relative flagging with an explicit note (production already guards this). |

Memory discipline is preserved throughout: `collect` parses one JSON at a time, extracts,
and releases it (`gc.collect()` between files). No command ever holds two cohorts' raw
JSON simultaneously.

## Testing

Fixtures are tiny and synthetic. The large JSONs and the existing caches are **not**
committed.

**Production fix**
- `normalise_sections`: v1.33 dict, v1.14 positional list, malformed input, non-dict
  members. Added to `test_check_multiqc.py`.
- `check_multiqc.run()` end-to-end against a list-shaped JSON — the regression that
  currently crashes.

**Tool**
- `spec`: unknown direction, missing direction, `relative` without `fail`, `reviewed`
  round-trip, un-gated metric excluded from emitted config.
- `manifest`: round-trip; hand-edited file with a `file://` URI.
- `cache`: subset metric list accepted; superset rejected with the missing names.
- `collect`: `'?'` dropped and counted, NaN filtered, gated-metric-missing raises and
  marks the cache incomplete, both section shapes handled.
- `stats`: percentiles and flag rates against hand-computed values; warn excludes fails;
  churn simulation on a fixed array.
- `relative`: min_cohort skip and zero-MAD skip at the tool level (the underlying
  `robust_threshold` / `relative_flags` behaviour is already covered).
- `emit`: **golden test** — a 3-cohort synthetic cache plus a fixed spec produces an exact
  expected TOML block; the block parses with `tomllib` and round-trips through
  `check_multiqc.load_thresholds`.
- `cli`: each subcommand's exit code on the failure paths above.

`discovery.py` is not unit-tested (it is a thin Metamist query); it is exercised manually
and isolated so nothing else depends on it.

**Manual regression, not committed:** run `distributions` and `flagrates` against the
existing `testing_scripts/testing_data/wgs_metric_values.json` (10 WGS cohorts) and
`wes_metric_values.json` (6 WES cohorts) and confirm they reproduce the numbers recorded
in `docs/qc_thresholds_tiers_plan.md`. This is the correctness check that the formalised
tool matches the manual analysis it replaces.

Target: 80%+ coverage on the subpackage.

## Dependencies and housekeeping

- **No new declared dependency.** `numpy` is already imported by
  `src/align_genotype/utils.py`; `click` and `cpg_utils` are direct dependencies.
  `tomllib` is stdlib on 3.11, and the project supports `>=3.10,<3.12` — but `tomli` is
  already in the resolved graph under a `python_full_version < '3.11'` marker, so a
  `try: import tomllib / except ImportError: import tomli as tomllib` shim covers 3.10
  without adding anything to `pyproject.toml`.
- **`metamist`** is imported lazily inside `discovery.py`.
- Writing TOML (manifest, spec) is done with a small hand-rolled writer in `manifest.py` /
  `spec.py` rather than adding `tomli-w` — the structures are flat and fully controlled.
- Operator artifacts (manifests, caches, specs) live under a gitignored `calibration/`
  directory. The tool and its tests are committed.
- Style: ruff line-length 120, single quotes, existing lint selection. Baseline has two
  pre-existing `PLR0917` in `check_multiqc` — unrelated, leave them.
- The existing `docs/*.md` QC notes stay untracked (operator's working notes). The
  operator-facing guide ships as `src/align_genotype/qc_calibration/README.md`. This
  design document is committed under `docs/superpowers/specs/` as the record of the
  decisions above; it is the only tracked file under `docs/`.

## Implementation order

TDD throughout — test first, watch it fail, then implement.

0. `normalise_sections` in `check_multiqc` + list-shape tests (the enabling fix).
1. `spec.py`, `manifest.py`, `cache.py` — the data layer.
2. `collect.py` — survey + extraction in one pass.
3. `stats.py` — `distributions`, `flagrates`.
4. `relative.py` — `mad`, driven through the production path.
5. `suggest`.
6. `emit.py` + the golden test.
7. `dryrun.py`.
8. `discovery.py` — Metamist.
9. `cli.py` wiring, `[project.scripts]` entry, README.
10. Manual regression against the two existing caches.

Steps 1–7 need no Metamist access and no large JSONs, so the bulk of the tool is
buildable and testable offline.

## Risks

- **`dryrun` config injection.** `CPG_CONFIG_PATH` may require config keys beyond
  `workflow.sequencing_type` for `check_multiqc.run()` to reach the threshold logic. If
  so, the temp TOML gains the minimum extra keys; monkeypatching `config_retrieve` (as the
  original scripts did) is the fallback. Resolved during step 7.
- **`normalise_sections` behaviour change.** Strictly widening — dict inputs behave
  identically — but it touches shipped code, so it lands as its own commit with its own
  tests, reviewable in isolation.
- **`suggest` seeds getting rubber-stamped.** Mitigated by `reviewed = false` plus the
  `emit-config` refusal, but it is a social risk the README must call out explicitly:
  seeds are evidence, not decisions.
