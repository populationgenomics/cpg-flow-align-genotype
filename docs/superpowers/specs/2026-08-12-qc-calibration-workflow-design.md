# Design: QC threshold calibration as a CPG Flow workflow

**Date:** 2026-08-12
**Status:** approved, ready for implementation planning

## Problem

`src/align_genotype/qc_calibration/` is currently an eight-command operator CLI
(`qc_calibrate`) driven by three hand-edited local artifacts: a manifest TOML, a value
cache JSON and a calibration spec TOML. It is run off a checkout, against files the
operator downloads or reads directly from GCS, and its final output is a
`[qc_thresholds.<seq_type>...]` block printed to stdout for pasting into
`config_template.toml`.

That shape is wrong for what this actually is. The analysis is a batch job over data the
pipeline already produces, its inputs are Metamist-registered `qc` analyses, and its
output is a report to review with the genomics analysis and bioinformatics teams. It
should be CPG Flow stages, run on Hail Batch, writing to the buckets every other stage
writes to.

This design converts it. Nothing in `qc_calibration/` has reached `main` — the branch is
68 commits ahead and the whole package, `scripts/sg_qc_report.py`, the
`sg_qc_overview.html.jinja` template and the `check_multiqc.py` changes are all new here.
So the restructure costs nothing in eventual PR diff; deleting the CLI surface makes the
diff against `main` smaller, not larger.

**Out of scope:** changing any shipped threshold *value*. The exome and genome numbers in
`config_template.toml` stand as-is and are the comparison baseline the report renders
against.

## Terminology

The current package calls "the set of samples in one MultiQC report" a **cohort**. Under
CPG Flow that is a **dataset**, and CPG Flow's own `Cohort` is a different thing
altogether — a collection of datasets and sequencing groups, assembled from
`input_cohorts` into a `MultiCohort`. Likewise "sample" is a **sequencing group**.

Throughout the new code: cohort → dataset, sample → sequencing group. `CohortValues`
becomes `DatasetValues`, `CohortMad` becomes `DatasetMad`, `min_cohort` becomes
`min_samples`, and so on. The one place the old word survives is CPG Flow's own API.

This includes the shipped production config key. `[qc_thresholds.<seq>.relative.<METRIC>]`
currently carries `min_cohort`, read by `check_multiqc.relative_flags`. In CPG Flow terms
that value means "the minimum number of sequencing groups in this MultiQC run", and since
`CramMultiQC` is a `DatasetStage` that run is exactly one dataset — the precise conflation
this rename exists to remove. It is renamed to `min_samples`: one line in
`check_multiqc.py` (`cfg.get('min_cohort', 0)`), two keys in `config_template.toml`. No
other consumer reads it.

## Corrections to the statistical picture

Recorded here because they change what the output can honestly claim.

**1. There is no "best" threshold, and none is computable.** `fail` is seeded from the
worst per-dataset p1 (`min` metrics) or p99 (`max` metrics); `warn` from p5 / p95. Those
are descriptions of where the tails sit in the datasets supplied, not decisions. The
decisions are explicitly non-computable — preserve the lab's intent for hard gates unless
the data contradicts it, don't hard-fail on metrics that track ancestry, biology or
chemistry, and override the lab only with a measured count of good samples their line
would have discarded. The report therefore presents **candidates**, and its most
decision-relevant content is not the candidate number but the per-dataset flag rates a
threshold would produce.

**2. A dataset-relative (MAD) tier has no threshold to ship.** `check_multiqc.relative_flags`
recomputes median and MAD from whatever sequencing groups are in that run's MultiQC report,
every run, and derives `median ± k·MAD/0.6745` on the fly. What reaches config is only
`direction`, `k` and `min_samples`. The per-dataset thresholds this tool reports are
**illustrative** — what today's data would give — and they move as a dataset grows. That
movement is what churn measures. Relative tiers are warn-only by construction and always
sit behind an absolute `fail` gate: within-dataset outlier detection finds nothing at all
in a uniformly poor dataset.

**3. Both currently-shipped relative tiers are expected to REJECT.** The advisory bar is
peak warn rate ≤ 10%, growth churn ≤ 2%, merge churn ≤ 5%. Over the full dataset set,
exome `ZERO_CVG_TARGETS_PCT` measures 1.0% growth / 51.1% merge and genome
`reads_duplicated_percent` measures 8.6% / 64.7%. Both were originally adopted on a
hand-picked subset (three datasets, leading slice only, two ordered pairs). There is a
structural tension the report must state plainly: a metric earns a relative tier
*because* its normal level shifts between datasets, and the merge simulation punishes
exactly that property. This is a live QC question for the team, not a defect, and the
report should frame it as one rather than let it read as a tool failure.

**4. Value counts and sequencing-group counts differ, and the fix is cheap.** MultiQC 1.33
can split one tool across general-stats sections (`picard_1` and `picard_4` both in the
Picard namespace), so one sequencing group contributes one value per section. The current
cache stores values without identity, making per-sequencing-group rates uncomputable — the
existing code documents a case reporting 4/28 = 14.3% where the true rate is 2/26 = 7.7%,
enough to flip a verdict. The new per-dataset values file **keeps the sequencing group ID
alongside each value**. Production keeps the duplication when deriving a threshold, so
calibration must too, but both counts can now be reported exactly instead of
approximated-with-a-caveat.

**5. A latent bug in discovery.** It selects the newest `qc` analysis whose
`meta.sequencing_type` matches. But `CramMultiQC` and `GvcfMultiQC` are both registered
`analysis_type='qc'`, and each registers two analyses (`analysis_keys=['json', 'html']`).
So it can return the GVCF MultiQC JSON, or an HTML file — neither carries any CRAM metric.
CPG Flow merges `get_job_attrs()` into analysis meta (`cpg_flow/stage.py:555`), which
includes `stage=<StageName>`, so the fix is a server-side meta filter on
`{stage: "CramMultiQC", sequencing_type: <type>}` plus a requirement that the output ends
in `multiqc_data.json`.

Incidentally noted: `dataset_stages.py:184` `_update_meta` is dead code — defined, never
referenced by any stage decorator.

## Architecture

Two stages in a new `qc_calibration_stages.py`, added to `run_workflow.py`'s stage list.
Neither declares a `required_stages` dependency on any production stage, so they form an
isolated branch of the DAG and never pull alignment or genotyping work in.

```
              driver: discovery.latest_cram_multiqc(dataset, seq_type)
                      metamist, lru_cached, one query per dataset
                                    │
QcCalibrationDatasetMetrics  (DatasetStage)
    ├── dataset-a ── job: read_input(multiqc_data.json) → extract → values.<analysis_id>.json
    ├── dataset-b ── job: …                              in that dataset's own main bucket
    └── dataset-n ── job: …           parallel; one report per job
                                    │
QcCalibrationReport  (MultiCohortStage, forced=True,
                      analysis_type='web', analysis_keys=['html'])
    └── one small job: read all N values files → cross-dataset analysis
          → <analysis-dataset>-main/qc_calibration/<seq_type>/calibration.json
          → <analysis-dataset>-web/qc_calibration/<seq_type>/calibration.html
```

### Why the split

Parsing a MultiQC report (tens to ~500 MB) is the only expensive operation, and it is
per-dataset and independent. Fanning it out means each job holds exactly one report,
which retires both the `gc.collect()`-between-datasets loop and the broad
`except Exception` containment block in `collect_all` — Hail already isolates per-job
failure. The cross-dataset work (percentile tails, flag rates, merge churn) then reads N
files of a few hundred KB each, so the report job is small and fast, preserving the
iterate-on-config loop the CLI's two-phase design existed to provide.

### Per-dataset output path

`dataset.prefix() / 'qc_calibration' / <seq_type> / f'values.{analysis_id}.json'`

Keying on the Metamist analysis ID means a new MultiQC report produces a new path, so
extraction re-runs exactly when the underlying data changed and never reads a stale
cache. Discovery therefore has to run inside `expected_outputs`, which is why it is
`lru_cache`d on `(dataset_name, seq_type)`.

`QcCalibrationDatasetMetrics` is deliberately **not** `forced` — that would defeat the
keying. `QcCalibrationReport` **is** `forced=True`, like `GenerateSgQcReport`, so the
report always regenerates from whatever values files exist.

### Datasets with no MultiQC report

`expected_outputs` returns `{}`. CPG Flow's `_is_reusable` treats a falsy expected output
as reusable (`cpg_flow/stage.py:723`), so the action is REUSE, `queue_jobs` is never
called, and no jobs and no Metamist analyses are created. The report stage filters these
out and names them in the provenance section as `skipped_datasets` with a reason.

### Staying inert during production runs

Being in `run_workflow`'s stage list with no dependencies means the DAG includes these
stages on every invocation. Their outputs won't exist, so they would queue a job per
dataset and register analyses. Both stages therefore return
`make_outputs(target, data=None, jobs=None)` unless `qc_calibration.enabled` is true, and
`config_template.toml` ships it `false`. `skip_stages` would also work but depends on
someone remembering; a default-off flag fails safe.

A calibration run sets `enabled = true` and
`only_stages = ['QcCalibrationDatasetMetrics', 'QcCalibrationReport']`. Both must be
listed: `_process_only_stages` marks every unlisted stage `skipped=True` with
`assume_outputs_exist=True`, and grants `assume_outputs_exist=False` to immediate
predecessors — so naming only the report would make the extract stage check for outputs
rather than produce them. Note `only_stages` is mutually exclusive with `skip_stages`,
`first_stages` and `last_stages`.

### Hail Batch conventions

No `gcloud storage cp`. Inputs via `batch.read_input(uri)`, outputs via
`batch.write_output(j.values, path)`, matching `jobs/multiqc.py` and `jobs/sg_qc_report.py`.

| Job | Resources | Rationale |
| --- | --- | --- |
| extract (per dataset) | `HIGHMEM`, `storage` configurable, default `20Gi` | a 500 MB JSON parses into several GB of Python objects; one report localised per job |
| report (per multicohort) | `STANDARD`, 2 CPU | reads N small JSON files |

Both jobs run the `workflow.driver_image` and invoke `python3 -m align_genotype.scripts.…`,
as `check_report_job` and `sg_qc_report_job` already do.

### Output paths use the literal sequencing type

`.../qc_calibration/genome/…` and `.../qc_calibration/exome/…`. This deliberately differs
from `CramMultiQC`, which uses `targets.sequencing_subdir()` — that helper returns `''`
for genome, so genome outputs would land unprefixed and the two runs would be asymmetric.

### Sequencing type scope

One sequencing type per run, from `workflow.sequencing_type`, as every other stage in the
repo does. Calibrating both means two invocations with two configs. A multicohort's
`input_cohorts` are chosen per sequencing type in practice, so a dual-type run would be
querying for data the multicohort wasn't assembled to cover.

## Configuration

```toml
[qc_calibration]
enabled = false          # default off, so production runs never queue these stages
k = 3.5                  # Iglewicz-Hoaglin outlier constant
min_samples = 50         # below this, MAD is too noisy to derive a relative tier
max_warn_rate = 0.10     # advisory bars behind the RECOMMEND / REJECT verdict
max_growth_churn = 0.02
max_merge_churn = 0.05
extract_storage = "20Gi"

[qc_calibration.genome.metrics.MEDIAN_COVERAGE]
direction = "min"        # min = higher is better, so a low value is flagged
unit = "x"               # x | % | frac — affects display and rounding only

[qc_calibration.genome.metrics.reads_duplicated_percent]
direction = "max"
unit = "%"
relative = true          # additionally evaluate a dataset-relative warn tier
```

`config_template.toml` ships both metric lists as the documented default. The candidate
set is not shared between sequencing types, because the Picard module differs — genome
uses `CollectWgsMetrics` (`MEDIAN_COVERAGE`, `PCT_20X`, …), exome uses `CollectHsMetrics`
(`MEAN_TARGET_COVERAGE`, `PCT_TARGET_BASES_20X`, `FOLD_80_BASE_PENALTY`,
`ZERO_CVG_TARGETS_PCT`, …), with samtools (`reads_mapped_percent`,
`reads_duplicated_percent`, `reads_properly_paired_percent`) and verifybamid (`FREEMIX`)
shared. The metric list stays a config parameter rather than a Python constant for
exactly that reason. Scope is limited to the metrics already found empirically relevant;
the mechanism generalises to any general-stats key, but no effort is spent on metrics
nobody has asked to gate.

Dropped relative to the current spec format: `gated`, `fail`, `warn`, `reviewed` and
`rationale`. They existed because an operator hand-authored candidate thresholds into a
file and iterated on them, with `reviewed` as the interlock stopping a raw percentile
reaching production config. The tool now derives the candidates and writes them to a
report, not to config, so there is nothing to sign off inside the tool. The judgement
those fields encoded moves to the README and to the report's own framing text.

## Module layout

| Current | Fate | Target size | Notes |
| --- | --- | --- | --- |
| `cli.py` (241) | delete | — | the stage is the entrypoint; `qc_calibrate` leaves `pyproject.toml` |
| `tomlio.py` (77) | delete | — | nothing reads or writes TOML; the emitted block is f-strings |
| `manifest.py` (131) | delete | — | no manifest artifact; provenance lives in the report JSON |
| `dryrun.py` (179) | delete | — | proved the emitted block loads; production proves that every run |
| `spec.py` (246) | → `settings.py` | ~90 | `MetricSpec` + direction/unit validation, read via `config_retrieve` |
| `cache.py` (227) | → `values.py` | ~90 | drops the `mkstemp`/`replace` dance, the `CloudPath` branch and `require_usable`; gains SG IDs |
| `collect.py` (220) | → `extract.py` | ~110 | one dataset per call; drops `collect_all`, `gc.collect()`, broad `except` |
| `stats.py` (113) | keep | 113 | percentiles, flag rates, churn — wording only |
| `relative.py` (499) | keep, cut | ~230 | see below |
| `suggest.py` (173) | → `thresholds.py` | ~110 | keeps tail selection and unit rounding; drops spec mutation and `reviewed` |
| `emit.py` (239) | → `snippet.py` | ~70 | keeps section ordering for a clean diff; drops both guards and `tomlio` |
| `report.py` (587) | → `summary.py` + `render.py` + template | ~190 | the ~250 lines of ASCII table machinery go; pure data assembly survives, plus a thin Jinja wrapper |
| `discovery.py` (196) | keep, cut | ~80 | drops `myProjects` and the eligibility filter; `gql()`, meta filter, non-lazy import |

New:

```
src/align_genotype/
  qc_calibration_stages.py                      ~110  two stages
  jobs/qc_calibration.py                         ~90  two job builders
  scripts/qc_calibration_extract.py              ~60  per-dataset job entry
  scripts/qc_calibration_report.py               ~90  report job entry
  templates/qc_calibration_report.html.jinja    ~350  mostly CSS, mirrors sg_qc_overview
  qc_calibration/render.py                       ~40  Jinja environment and template render
  qc_calibration/README.md                             rewritten
```

The package goes from 3,133 lines to roughly 1,050 of library plus ~350 of
stages/jobs/scripts.

### The one substantive cut in `relative.py`

Delete `_production_config`, `_restore_config_paths` and `_warn_count` — about 80 lines
that write a throwaway TOML, install it with `config.set_config_paths()`, call
`check_multiqc.relative_flags` to count warns, then attempt to restore the previous paths.
Defensible in a local CLI; inside a Hail Batch job, where CPG Flow has already installed
the run's config and every later `config_retrieve` depends on it, mutating global config
paths for a warn count is a hazard.

Replacement: call `check_multiqc.robust_threshold` and `stats.breach` directly. The "no
second MAD implementation anywhere" invariant is fully preserved — `robust_threshold` is
production's own function, and `stats.churn` already calls it that way. What routing
through `relative_flags` additionally bought was exercising `load_thresholds`, which was
`dryrun`'s purpose and is now redundant.

### Discovery

```python
ANALYSES_QUERY = gql("""
    query CramMultiqc($dataset: String!, $analysisType: String!, $metaFilter: JSON!) {
        project(name: $dataset) {
            analyses(status: {eq: COMPLETED}, type: {eq: $analysisType}, meta: $metaFilter) {
                id
                output
                timestampCompleted
            }
        }
    }
""")
```

Meta filtering takes a per-key operator dict, not a flat mapping, and both it and the
analysis type are passed as query variables:

```python
meta_filter = {
    'stage': {'eq': 'CramMultiQC'},
    'sequencing_type': {'eq': seq_type},
}
result = query(
    ANALYSES_QUERY,
    variables={
        'dataset': dataset_name,
        'analysisType': 'qc',
        'metaFilter': meta_filter,
    },
)
```

Then keep only outputs ending `multiqc_data.json` and take the newest
`timestampCompleted`. Datasets come
from `multicohort.get_datasets()`, so the `is_seqr` check and the
`('test', 'training', 'seqr')` name-substring exclusions are gone — dataset selection is
now the operator's, expressed through `input_cohorts`. `metamist` is imported normally
(it ships with cpg-flow). An analysis with a non-string `timestampCompleted` can't be
ranked and is skipped with a warning.

## Outputs

### `calibration.json` — analysis dataset main bucket

```text
sequencing_type, generated, ar_guid
settings            k, min_samples, max_warn_rate, max_growth_churn, max_merge_churn
datasets[]          dataset, analysis_id, timestamp, uri, n_sequencing_groups,
                    multiqc_version
skipped_datasets[]  dataset, reason
metrics{KEY}        direction, unit, n_values, n_groups_with_values, n_datasets,
                    n_dropped, present_in[], missing_from[],
                    current{fail, warn}, candidate{fail, warn, basis},
                    flag_rates{current{dataset: [fail, warn]}, candidate{…}},
                    percentiles{dataset: {p1 … p99}}
relative{KEY}       verdict, reason, max_warn_rate, max_growth_churn, max_merge_churn,
                    datasets[]{n_values, n_groups_with_values, median, mad, threshold,
                               n_warn, warn_rate, skipped},
                    growth[]{ordered, shuffled, worse, ordering_sensitive},
                    merge_worst[]
warnings[]
config_snippet      the qc_thresholds TOML block, as a string
```

`n_values` and `n_groups_with_values` both appear per metric — correction 4 made visible
rather than caveated. The two are named differently on purpose: `n_sequencing_groups` on a
dataset is that dataset's total, while `n_groups_with_values` counts only the groups
carrying that metric, and they diverge whenever a metric is partially missing. A rate over
a metric divides by `n_values` or `n_groups_with_values`, never by the dataset total. `current` is read from the live
`[qc_thresholds.<seq_type>...]` config, so every candidate is presented against what is
shipped today.

### `calibration.html` — analysis dataset web bucket

Rendered from a Jinja template, following `sg_qc_overview.html.jinja`. Sections in order:

1. **Header** — sequencing type, sequencing group count, dataset count, metric count,
   `k`, `min_samples`, date, ar-guid.
2. **Banner** — red, if any configured metric is missing from any dataset.
3. **Recommended fixed thresholds** — metric, direction, current fail/warn, candidate
   fail/warn, flag-rate delta. Carries the caveat that candidates are percentile tails,
   not decisions.
4. **Dataset-relative (MAD) tiers** — per metric: verdict, and the three measured numbers
   against their bars; then per dataset: n, median, MAD, illustrative threshold, warn
   count and rate, or the skip reason. Carries the merge-churn framing paragraph.
5. **Collapsible evidence** — metric presence matrix; percentile distributions per
   dataset; flag rates per dataset for current and candidate; churn detail (ordered vs
   shuffled growth with the `ORDERING-SENSITIVE` marker, worst merge pairs); provenance
   and skipped datasets.
6. **Copy-pasteable `qc_thresholds` block** in a `<pre>` with a copy button.

No separate `.toml` output file — the snippet exists in the HTML and in the JSON.

### Failure policy

| Condition | Behaviour |
| --- | --- |
| Configured metric missing from *some* datasets | red banner, `warnings[]` entry, WARNING log. Not a job failure. |
| Configured metric missing from *every* dataset | loudest banner, ERROR log. Still not a job failure. |
| Report has no usable `report_general_stats_data` | extract job fails, mirroring `check_multiqc.load_sections` |
| Non-numeric value (Picard `'?'`) | dropped and counted per metric, as production does; surfaced in the report |
| NaN / inf | filtered before percentiles and MAD |
| Dataset below `min_samples` | no relative tier derived; skip reason recorded and displayed |

The report job never fails on a missing metric because Hail only copies `write_output`
targets on job success — failing would destroy the HTML that explains the problem. This
is a deliberate departure from the CLI, where `collect` exited non-zero and marked the
cache `complete = false` precisely so nothing downstream ran on it. The banner replaces
that interlock: there is no longer a downstream consumer to protect, only a human reader
to inform, and the reader is guaranteed to open the page since it is the whole reason the
run happened.

## Testing

Twelve files, ~1,320 lines, down from thirteen and 3,692. Almost all of the reduction is
a consequence of deleted module surface rather than thinned coverage. The file *count*
barely moves because two new suites appear — `render` (template smoke tests) and
`scripts` (the two job entrypoints) — that the old CLI had no equivalent of; the line
count is what falls.

| | now | after | why |
| --- | --- | --- | --- |
| `tomlio`, `manifest`, `cli`, `dryrun` | 714 | 0 | modules deleted |
| `spec` → `settings` | 284 | ~80 | tested validation rules and TOML round-trips that no longer exist |
| `report` → `summary` | 703 | ~150 | see below |
| `relative` | 538 | ~250 | `set_config_paths` restore-path tests go with `_production_config` |
| `discovery` | 324 | ~120 | eligibility-filter tests go; gains the CramMultiQC-vs-GvcfMultiQC test |
| `collect` → `extract` | 288 | ~140 | gains sequencing-group-ID retention |
| `cache` → `values` | 289 | ~70 | atomic-write and CloudPath-branch tests go with the code |
| `emit` → `snippet` | 226 | ~60 | one golden-block test |
| `suggest` → `thresholds` | 184 | ~90 | |
| `stats` | 142 | 142 | unchanged |
| `stages` | — | ~120 | new |
| `render` | — | ~90 | new: template smoke tests |
| `scripts` | — | ~110 | new: the two job entrypoints |

What must be covered:

- **`settings`** — bad direction, bad unit, missing metrics table, `relative` defaulting.
- **`values`** — round-trip with sequencing group IDs; non-numeric value rejected by name.
- **`extract`** — `'?'` dropped and counted per metric; NaN filtered; both MultiQC section
  shapes (v1.33 dict, v1.14 positional list); presence recorded per section; one
  sequencing group appearing in two sections yields two values and one SG count.
- **`stats`** — percentiles and flag rates against hand-computed values; warn excludes
  already-failing; churn on a fixed array.
- **`relative`** — `min_samples` skip; zero-MAD skip; growth taking the worse of ordered
  and shuffled; `ORDERING-SENSITIVE` detection; merge simulating both directions of every
  pair; each verdict reason naming its own bar.
- **`thresholds`** — tail selection per direction; unit rounding; a metric with no data
  skipped rather than fabricated.
- **`snippet`** — golden block, parses with `tomllib`, round-trips through
  `check_multiqc.load_thresholds`.
- **`summary`** — numeric assertions against a `CalibrationSummary` built from a synthetic
  three-dataset fixture; one smoke test that the template renders and contains the
  expected section headings.
- **`discovery`** — newest-by-timestamp wins over higher ID; unrankable timestamp skipped;
  a GvcfMultiQC analysis and an HTML output are both excluded.
- **`stages`** — output paths per sequencing type; `{}` from `expected_outputs` when no
  analysis exists; no jobs when `enabled = false`.
- **`check_multiqc`** — existing coverage stands; add the `min_samples` rename.

`test_qc_calibration_report.py` is cut on principle regardless of the refactor: 703 lines
asserting exact ASCII column alignment tests a renderer, not a calculation, and the
renderer is being deleted.

Fixtures stay tiny and synthetic. Real MultiQC reports and real dataset names are never
committed.

## Documentation

`qc_calibration/README.md` is rewritten to stand alone: what the workflow is for, the
analysis-runner invocation and the config it needs, what each output is, and how to read
the report. The judgement prose is the most valuable content in the current README and it
survives — aim for ~0% fail and single-digit % warn on a healthy dataset; preserve the
lab's intent for hard gates; the bar for overriding the lab is a measured count of good
samples their line would discard; don't hard-fail on metrics tracking ancestry, biology or
chemistry; the two conditions for adopting a relative tier are genuinely
dataset-dependent spread *and* a flag set that stays stable as the dataset grows.

Removed everywhere: references to `testing_scripts/`, "the old workflow", "the manual
process this replaces", "ported verbatim from `fetch_multiqc_json_paths.py`", "the old
scripts monkeypatched `config_retrieve`", and the `PCT_PF_READS_ALIGNED` incident
narrative. The operative fact stays as one sentence — that metric exists only in
`report_saved_raw_data`, so it can never be gated; use samtools `reads_mapped_percent` —
because that is still a live trap. How it was discovered is not.

Comment policy: a one-line summary plus only what a reader would otherwise get wrong.
`discovery.latest_cram_multiqc` gets about three lines, not twenty about malformed GraphQL
rows. Long comments survive only where they encode a decision someone would otherwise
undo — why merge churn is judged against a looser bar than growth, why growth takes the
worse of the two slice orderings, why duplicated cross-section values are kept rather
than de-duplicated.

Also in scope:

- `config_template.toml` gains the `[qc_calibration]` block with both metric lists, and
  its `qc_thresholds` comments lose the `testing_scripts/mad_relative_prototype.py` /
  `verify_relative_6cohorts.py` citations (lines 146–147) and the similar reference at
  157–161. Both `min_cohort` keys become `min_samples`.
- `docs/superpowers/specs/2026-08-11-qc-calibration-workflow-design.md` and
  `docs/superpowers/plans/2026-08-11-qc-calibration-workflow.md` were deleted alongside
  this document's first commit. They describe the CLI being removed and reference
  `testing_scripts/` throughout; this document replaces them.
- `.gitignore` keeps both `testing_scripts/` and `calibration/`. Nothing writes to them
  any more, but local copies must stay ignored.

**Not a CLI, but locally runnable.** Both job entry scripts are ordinary `python -m`
modules, so
`python -m align_genotype.scripts.qc_calibration_report --values a.json b.json --output out.html`
works off a checkout against downloaded per-dataset values files. That keeps the fast
iterate-on-thresholds loop without a `click` group to maintain.

## Implementation order

TDD throughout — test first, watch it fail, then implement.

1. `check_multiqc.py`: `min_cohort` → `min_samples`, plus `config_template.toml`. Its own
   commit: it touches shipped behaviour and should be reviewable in isolation.
2. `settings.py` and `values.py` — the data layer, plus the `[qc_calibration]` config
   block in `config_template.toml`.
3. `extract.py` — one dataset, one report, with sequencing group IDs retained.
4. `stats.py` and `relative.py` — the rename, and the `_production_config` excision.
5. `thresholds.py` and `snippet.py`.
6. `summary.py` and the Jinja template.
7. `scripts/qc_calibration_extract.py`, `scripts/qc_calibration_report.py`,
   `jobs/qc_calibration.py`.
8. `qc_calibration_stages.py`, `run_workflow.py` wiring.
9. Delete `cli.py`, `tomlio.py`, `manifest.py`, `dryrun.py`, the `qc_calibrate` entry
   point, and the four obsolete test files.
10. README, and a dry-run of the workflow against a real multicohort.

Steps 1–6 need no Metamist access and no large reports, so the bulk is buildable and
testable offline.

## Risks

- **Report size at 16 datasets.** The percentile and flag-rate tables are metric ×
  dataset. At ~14 metrics and ~16 datasets that is legible in HTML where it was not in a
  terminal, but the presence matrix should still be checked for width once real labels
  are in it.
- **`enabled` flag drift.** If someone sets `enabled = true` in a shared config, every
  production run gains a job per dataset. It belongs only in the calibration run's config,
  and the README must say so.
- **`min_samples` rename.** Any run config that overrides `min_cohort` silently stops
  taking effect — the value then comes from `config_template.toml`'s `min_samples = 50`
  rather than the operator's intended override, and `cfg.get('min_samples', 0)` would fall
  back to 0 for any relative metric the template doesn't cover. Grep for the key across
  run configs before merging; there are no in-repo occurrences beyond
  `config_template.toml` and `check_multiqc.py`.
