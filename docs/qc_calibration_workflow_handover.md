# Handover: formalise the QC threshold calibration into a reusable workflow

**Goal:** turn the ad-hoc, machine-local calibration scripts in `testing_scripts/`
(gitignored) into a **committed, executable, reusable workflow** that derives
warn/fail (and cohort-relative MAD) QC thresholds from a set of MultiQC datasets and
emits the `config_template.toml :: qc_thresholds.<seq_type>` block. The analysis has
already been done twice by hand (exome + genome); this task is to consolidate the
throwaway scripts into a single tool a future operator can run when onboarding a new
capture kit / sequencing protocol or refreshing thresholds.

**Explicitly out of scope for you (the human operator will provide this):** how to
**retrieve the MultiQC JSONs from the database / Metamist / GCS**. Today all inputs are
hand-downloaded into local folders. The operator will give you the retrieval mechanism
separately — your job is to design the tool so the data source is pluggable (see
["Data source abstraction"](#data-source-abstraction)). Do **not** hardcode a local
folder as the only input path. (operator note: see @testing_scripts/fetch_multiqc_json_paths.py for this)

---

## 1. What already exists (read these first)

- **`docs/qc_thresholds_tiers_plan.md`** — the design rationale for the warn/fail tier
  system and the calibrated exome + genome results. Authoritative on *what* the
  thresholds are and *why*.
- **`docs/qc_flagging_update.md`** — the original lab guidance the tiers derive from.
- **`docs/wgs_threshold_calibration_handover.md`** — the handover that drove the genome
  calibration; a good example of the manual per-run process you are formalising.
- **Shipped config + code** (branches, both have open PRs):
  - `QC_flags_report` (**PR #73**) — warn/fail tier system, `QcFlag.severity`, the
    report, **calibrated absolute exome + genome thresholds**.
  - `mad-relative-flagging` (**PR #74**, base = `QC_flags_report`) — the cohort-relative
    (MAD) code and the two relative tiers: exome `ZERO_CVG_TARGETS_PCT` + genome
    `reads_duplicated_percent`. These two branches are the *output* of the manual
    calibration you are automating.
- **Production consumer:** `src/align_genotype/scripts/check_multiqc.py`. Key functions
  you must stay compatible with:
  - `load_thresholds(seq_type)` → reads `['qc_thresholds', seq_type, severity, direction]`.
  - `robust_threshold(values, direction, k)` → the Iglewicz–Hoaglin modified-z line
    (`median ± k·MAD/0.6745`); returns `None` on zero MAD.
  - `relative_flags(sections, seq_type, today, already_flagged)` →
    reads `['qc_thresholds', seq_type, 'relative']` **generically** (any seq_type),
    warn-only, defers to absolute flags. `apply_relative_flags(...)` wires it in.
  - `_gather_metric_values` / `warn_unmatched_metrics` — metric extraction + missing-key
    logging. **The calibration tool should extract values the same way the production
    checker does, so calibration and enforcement never disagree.**

---

## 2. The scripts to consolidate (all in `testing_scripts/`, gitignored)

These are the raw material. Each is a standalone script edited in place between the
exome and genome runs (that mutation-in-place is exactly the smell to design out).

| Script | Role | Key I/O |
|---|---|---|
| `survey_wes_keys.py` | **Step 1 — survey.** Confirms each candidate metric key is present in `report_general_stats_data` across every dataset; reports MultiQC version + section shape; warns on `MISSING` (possible key rename). | in: folder of JSONs; out: stdout table |
| `qc_metric_distributions.py` | **Step 2 — distributions.** Parses each JSON once (memory-safe), extracts curated metrics, caches values to a small JSON, prints per-cohort percentile tables (p1…p99). | in: folder; out: `*_metric_values.json` cache + stdout |
| `qc_threshold_flagrates.py` | **Step 3 — flag rates.** Reads the cache, reports the fraction of samples each candidate `fail`/`warn` threshold would flag, per cohort. Fast — iterate here. | in: cache; out: stdout |
| `mad_relative_prototype.py` | **Step 4 — MAD spread + churn.** For protocol-dependent metrics, computes per-cohort MAD thresholds/flag-rates and simulates cohort-growth churn (homogeneous + heterogeneous). | in: cache; out: stdout |
| `verify_relative_6cohorts.py` | **Step 4b — production verification.** Imports `check_multiqc.relative_flags`/`robust_threshold` and drives them with the cache (config monkeypatched) so the numbers reflect exactly what ships. | in: cache; out: stdout |
| `dryrun_exome_check.py` / `dryrun_genome_check.py` | **Step 5 — end-to-end sanity.** Runs the real `check_multiqc.run()` against one big JSON (single-load + Slack off + `json.load` monkeypatched to avoid a double parse), prints flag counts by `(metric, severity, method)`. | in: one JSON; out: stdout + report html |

Caches currently on disk: `testing_scripts/testing_data/wes_metric_values.json`
(6 WES cohorts) and `wgs_metric_values.json` (10 WGS cohorts). Treat these as
**disposable artifacts**, not committed fixtures — but they're handy golden data to
test the formalised tool against (its distribution/flag-rate output should reproduce
the numbers in `qc_thresholds_tiers_plan.md`).

---

## 3. The calibration methodology (the process to encode)

This is the reusable pipeline, independent of exome/genome:

1. **Survey** — for the candidate metric list, confirm every key is present in general
   stats across *all* datasets. Never trust a metric that's `MISSING` on some.
2. **Distributions** — parse each JSON once, extract values, cache them, read the
   percentile tables. For `min` metrics (higher=better) the bad samples are the *left*
   tail; for `max` metrics the *right* tail.
3. **Flag rates** — from the cache, tabulate what each candidate threshold flags per
   cohort. Tune so a **healthy cohort ≈ 0% fail / single-digit % warn** while genuinely
   poor batches still surface.
4. **MAD (only where justified)** — if a metric's normal level shifts by cohort/protocol
   (bimodal or wide-ranging medians), evaluate cohort-relative warn flagging: per-cohort
   warn rate should land in single digits and cohort-growth churn should be low (aim
   <1%; we accepted ~2% worst-case for a contrived cross-project merge). Keep it
   **warn-only** and always keep an absolute `fail` hard gate. Verify via the **production**
   `relative_flags()` path, not just the prototype.
5. **Apply + verify** — write the `config_template.toml` block with a per-metric
   rationale comment citing the cohort percentiles / flag rates; validate the TOML parses;
   run the test suite; optional end-to-end dry run on one real JSON.

---

## 4. Domain knowledge / gotchas learned this session (don't re-derive these)

These are the things that cost time to discover. Bake them into the tool.

### MultiQC file shapes differ by version
- **v1.33** stores `report_general_stats_data` as a **dict** keyed by section name
  (`verifybamid`, `picard`, `picard_3`, `picard_4`, `samtools`).
- **v1.14** stores it as a **list** (sections are positional: seen as
  `section_0=verifybamid`, `section_2=coverage`, `section_3=samtools`).
- Both must be normalised. The existing helpers do this: `iter_sections` (survey) and
  `_section_dicts` (distributions). The production `check_multiqc` has its own
  normaliser. **Pick one normaliser and reuse it everywhere** — ideally import the
  production one so calibration and enforcement agree.
- Always print/inspect the MultiQC version + shape per file; a file that yields *no*
  `report_general_stats_data` should be surfaced loudly, not silently skipped.

### Metric keys differ by sequencing type (Picard module)
- **Exome** = Picard `CollectHsMetrics`: `MEAN_TARGET_COVERAGE`,
  `PCT_TARGET_BASES_20X/50X`, `FOLD_80_BASE_PENALTY`, `ZERO_CVG_TARGETS_PCT`,
  `PCT_SELECTED_BASES`, `PCT_OFF_BAIT`, `AT/GC_DROPOUT`.
- **Genome** = Picard `CollectWgsMetrics`: `MEDIAN_COVERAGE`, `MEAN_COVERAGE`,
  `SD/MAD_COVERAGE`, `PCT_1X…PCT_100X` (breadth), `PCT_EXC_*`, `HET_SNP_SENSITIVITY`,
  `GENOME_TERRITORY`. **No** `*_TARGET_*` / `FOLD_80` / `ZERO_CVG` keys exist.
- Shared across both (samtools): `reads_mapped_percent`, `reads_duplicated_percent`,
  `reads_properly_paired_percent`, `reads_MQ0_percent`, `error_rate`. Contamination
  (verifybamid): `FREEMIX`.
- **Trap:** `PCT_PF_READS_ALIGNED` is **not** in `report_general_stats_data` (only in
  `report_saved_raw_data`), so it can never be gated — use samtools `reads_mapped_percent`
  instead. This is the bug that made the old genome reads-mapped gate inert. **The tool's
  survey step is what catches this class of bug — keep it mandatory.**
- So the candidate metric list is **seq_type-specific**. Make it a config/parameter,
  not a hardcoded constant.

### Memory
- These JSONs are **large** (tens to ~400 MB; one WES file was ~500 MB). **Parse one at
  a time** and release it (`gc.collect()` between files helped). Never load all of a
  cohort set at once.
- The two-phase design (parse-once → small value cache → iterate on the cache) is
  essential: Step 3/4 tuning must be fast and must not re-parse hundreds of MB. Preserve
  this in the formalised tool (cache keyed by cohort label + metric).

### Data hygiene
- Picard writes non-numeric placeholders (e.g. `'?'`); coerce with `float()` and drop
  failures (count them). NaNs must be filtered before percentiles/MAD.
- `min_cohort` matters for MAD: below ~50 samples MAD is too noisy — skip relative
  flagging and fall back to absolute only. The production code already guards this.

### Direction convention
- `min` = higher is better (flag when **below** threshold; bad = low tail).
- `max` = lower is better (flag when **above** threshold; bad = high tail).
  Every metric needs an explicit direction; the tool should carry `(key, direction, unit)`.

### MAD specifics (matches production)
- Iglewicz–Hoaglin modified z-score, **k=3.5** standard outlier line,
  `threshold = median ± k·MAD_raw/0.6745`, `MAD_raw = median(|x − median|)`.
- Round the derived threshold (production rounds to 4 dp) to damp sub-0.0001 jitter that
  would otherwise churn "updated" flags.
- Only adopt MAD where spread is genuinely protocol-dependent **and** churn is low.
  Rejected examples: exome `PCT_SELECTED_BASES` / `PCT_OFF_BAIT` (churn up to 24.5%).
  Adopted: exome `ZERO_CVG_TARGETS_PCT`, genome `reads_duplicated_percent`.

---

## 5. Config schema + production integration

The tool's ultimate output is TOML under `[qc_thresholds.<seq_type>....]`:

```toml
[qc_thresholds.<seq_type>.fail.min]   # flag when value < threshold
[qc_thresholds.<seq_type>.fail.max]   # flag when value > threshold
[qc_thresholds.<seq_type>.warn.min]
[qc_thresholds.<seq_type>.warn.max]
[qc_thresholds.<seq_type>.relative.<METRIC>]   # warn-only MAD tier
direction = "max"   # or "min"
k = 3.5
min_cohort = 50
```

- A metric may appear in only the tiers that make sense (e.g. warn-only metrics live
  only in `warn`; a MAD metric has an absolute `fail` **and** a `relative` entry, no
  absolute `warn`).
- `check_multiqc` evaluates **fail before warn** and records one flag per
  `(section, metric)` at the worst tier; the relative pass is warn-only and won't
  double-flag a sample already caught absolutely.
- **No production code change is needed to add a new `relative` metric or a new
  seq_type** — `relative_flags` and `load_thresholds` are seq_type-generic. The tool
  only needs to emit config.

---

## 6. Design guidance for the formalised workflow

### Data source abstraction
Model the input as an iterable of **`(cohort_label, json_source)`**, where `json_source`
resolves to a parseable MultiQC JSON (local path, GCS blob, DB blob, stream — the
operator will provide the concrete retrieval). Requirements:
- **Lazy + one-at-a-time**: yield/parse a single cohort at a time; never materialise all
  JSONs. A generator that downloads → parses → extracts → discards is ideal.
- Keep a **local value cache** (the small per-metric JSON) so re-tuning thresholds never
  re-downloads/re-parses. Cache key = cohort label; invalidate when the metric list
  changes.
- Don't assume filenames; derive `cohort_label` from whatever the retrieval layer gives
  (project/run id). Today labels came from stripping `_multiqc_data.json`.

### Suggested shape
- A single module/CLI (e.g. `src/align_genotype/scripts/qc_calibration.py`, or a small
  package) with subcommands mirroring the 5 steps: `survey`, `distributions`,
  `flagrates`, `mad`, `emit-config` (+ optional `dryrun`). Wire it into
  `[project.scripts]` in `pyproject.toml` if it should be a first-class entrypoint.
- Parameterise: `--seq-type`, a metric spec (list of `(key, direction, unit)` +
  candidate `fail`/`warn` values and optional `relative` flag), `--cache`, and the
  data-source selector. Consider a small TOML/YAML "calibration spec" file rather than
  editing Python constants (that in-place editing is the current pain point).
- **Reuse production extraction/normalisation** by importing from `check_multiqc`
  instead of re-implementing `_section_dicts`/`_gather_metric_values`. Fold the
  standalone `verify_*` step into the tool so "what the flag rate will actually be in
  production" is a first-class output, not a separate script.
- Outputs: the percentile tables, flag-rate tables, MAD spread/churn, and a
  ready-to-paste `[qc_thresholds.<seq_type>...]` block **with generated rationale
  comments** citing the numbers (match the hand-written style already in the config).

### Tests
- Unit-test the normaliser against both a **v1.33 dict** and a **v1.14 list** fixture
  (tiny synthetic ones — do *not* commit the huge JSONs).
- Test direction/threshold flag logic, MAD `robust_threshold` (incl. zero-MAD → skip and
  `min_cohort` skip), and the non-numeric/NaN dropping.
- A golden test: feed a small synthetic multi-cohort set and assert the emitted config
  block + flag-rate numbers. The existing caches can seed a regression fixture.

### Housekeeping
- The current scripts are gitignored under `testing_scripts/`. The formalised tool
  should be **committed** under `src/` (or a `tools/`), with its heavy inputs/caches kept
  out of git. Keep `dryrun_*` capability but generalise the two into one seq_type-aware
  command.
- `docs/` is currently **untracked** in this repo — decide with the operator whether the
  calibration docs get committed alongside the tool.

---

## 7. Guardrails / principles (carry over from the manual process)

- **Healthy cohort ≈ 0% fail; fail = "do not analyse"; warn = "a human should look".**
- Preserve the lab's intent for hard gates unless the data clearly says otherwise; you
  can tighten `warn` freely. (Example: genome `reads_duplicated_percent` fail was
  *relaxed* 25 → 40 because the data showed 25 would hard-fail ~25–37% of legitimately
  high-dup preps.)
- Don't hard-fail on metrics that vary with ancestry/biology/chemistry (`error_rate`,
  `HET_SNP_SENSITIVITY`) — prefer warn or cohort-relative.
- Don't shoehorn a MAD tier onto every metric — only where spread is clearly
  protocol-dependent *and* churn stays low. Log what you rejected and why.
- Parse one big JSON at a time; never load them all.
- The survey step is mandatory and non-negotiable — it's the guard against silently
  gating a key that MultiQC renamed or never surfaces.

---

## 8. First steps checklist

1. Read `docs/qc_thresholds_tiers_plan.md` + skim PRs #73/#74 to see the target output.
2. Get the DB/GCS retrieval mechanism from the operator; wrap it behind the
   `(cohort_label, json_source)` iterator.
3. Import extraction/normalisation from `check_multiqc` (or refactor a shared helper out
   of it) so calibration ≡ enforcement.
4. Port the 5 steps into subcommands; replace in-place constant editing with a
   calibration-spec input.
5. Reproduce the known results (exome + genome numbers in the plan doc) from the existing
   caches as a correctness check, then add tests with tiny synthetic fixtures.
6. Emit the config block with generated rationale comments; validate TOML; run
   `uv run --extra test pytest test/ -q`; optional `dryrun` on one real JSON.

Run scripts with `uv run --extra test python3 <script>` (how the current scripts are
invoked). Baseline lint: `uvx ruff check <files>` (two pre-existing `PLR0917` in
`check_multiqc` are expected — ignore).

---

## 9. Deliverables
- A **committed, reusable, documented tool** that implements the 5-step calibration
  workflow, with a pluggable data source and a small value cache.
- A **test suite** that covers the normaliser, direction/threshold logic, MAD, and a golden regression test.
- A **calibration spec** (TOML/YAML) that captures the thresholds, directions, and any cohort-relative settings, which will ultimately be used to inform the thresholds in the config_template.toml.