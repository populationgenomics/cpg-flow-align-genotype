# Handover: calibrate genome (WGS) QC thresholds from real cohorts

**Goal:** derive sensible warn/fail QC thresholds for **genome** sequencing by running
the existing analysis scripts over 10 real WGS MultiQC datasets, then update
`config_template.toml :: qc_thresholds.genome.*`. This mirrors the WES calibration
already completed — read that first as the template.

## Context (what's already done)

- Warn/fail severity tiers, exome thresholds, and the report template are on branch
  `QC_flags_report` (merged/base).
- Exome thresholds were **calibrated** against 6 real WES cohorts; MAD cohort-relative
  warn flagging for the kit-dependent `ZERO_CVG_TARGETS_PCT` is in **PR #74**
  (`mad-relative-flagging` → `QC_flags_report`).
- The full method, rationale, and results are in `docs/qc_thresholds_tiers_plan.md`.
- **Genome thresholds are NOT yet calibrated.** Current state in
  `config_template.toml`:
  - `[qc_thresholds.genome.fail.*]` = the lab's original single-tier values
    (`MEDIAN_COVERAGE` 10, `reads_mapped_percent` 80, `FREEMIX` 0.04,
    `reads_duplicated_percent` 25).
  - `[qc_thresholds.genome.warn.*]` = generic values from `docs/qc_flagging_update.md`
    (`MEDIAN_COVERAGE` 25, `reads_mapped_percent` 98, `FREEMIX` 0.01,
    `reads_duplicated_percent` 10).
  Your job is to replace both tiers with cohort-calibrated values. Additionally, consider any other genome metrics which may be sensible to include. Also, if any of the genome metrics show strong cohort-dependent spread, consider adding a cohort-relative MAD-based warn tier (like exome `ZERO_CVG_TARGETS_PCT`) - but only if the spread is stable across cohorts and the MAD-based flagging is not too sensitive to cohort growth (see Step 4). Don't shoehorn a MAD-based tier for every metric — only if the data clearly supports it and it would be useful to the lab.

## The data

10 WGS MultiQC JSONs live in `testing_scripts/testing_data/WGS/` (gitignored; large,
tens–hundreds of MB each):

```
acute-care  afhcs  chop-gliadx  circa  ghfm-kidgen
heartkids   ibmdx  perth-neuro  rdnow  schr-neuro
```

(Some names overlap the WES set — these are the **WGS** runs of those projects.)

## Key differences from the WES work — READ THIS

1. **Genome uses Picard `CollectWgsMetrics`, not `CollectHsMetrics`.** So the
   general-stats metric keys are different. There are **no** `*_TARGET_*` / `FOLD_80` /
   `ZERO_CVG` / on-bait keys. The WGS keys (seen in a real genome JSON) live in the
   Picard section as: `MEDIAN_COVERAGE`, `MEAN_COVERAGE`, `SD_COVERAGE`,
   `MAD_COVERAGE`, `PCT_1X`…`PCT_100X` (coverage breadth), `PCT_EXC_TOTAL`,
   `PCT_EXC_DUPE`, `PCT_EXC_MAPQ`, `HET_SNP_SENSITIVITY`, `GENOME_TERRITORY`; and in
   samtools: `reads_mapped_percent`, `reads_duplicated_percent`,
   `reads_properly_paired_percent`, `error_rate`, `reads_MQ0_percent`; and in
   verifybamid: `FREEMIX`.
2. **`PCT_20X` (breadth) is the recommended primary depth gate for WGS** (per
   `docs/qc_flagging_update.md`), alongside `MEDIAN_COVERAGE`.
3. **Duplication is low on PCR-free WGS** (~1–8%), so `reads_duplicated_percent`
   thresholds should be much tighter than exome's.
4. **`PCT_PF_READS_ALIGNED` is NOT in `report_general_stats_data`** (it's only in
   `report_saved_raw_data`) — this is the bug that made the old genome gate inert.
   Use samtools `reads_mapped_percent` for the reads-mapped gate. **Always run the
   survey step first to confirm each metric key is actually present in general stats.**
5. **Older MultiQC versions** store `report_general_stats_data` as a *list*, not a
   dict (seen in the WES set at v1.14). The scripts already normalise both shapes
   (`_section_dicts` / `iter_sections`), so key detection works regardless — but check
   the survey output for any `MISSING (possible key rename)` warnings.

## Step-by-step

All scripts are in `testing_scripts/` (gitignored). Run with `uv run --extra test
python3 <script>`. They parse **one file at a time** (memory-safe) — keep it that way.

### Step 1 — survey keys, versions, and shapes
Edit `testing_scripts/survey_wes_keys.py`:
- set `WES_DIR = 'testing_scripts/testing_data/WGS'`
- replace the `METRICS` / `CONFIG_METRICS` lists with the WGS candidates (Step 2 list).

Run it. Confirm every metric you intend to gate is present in general stats across all
10 datasets; note any older-version list-form files. Do **not** trust a metric that
shows `MISSING` on some datasets.

### Step 2 — per-metric distributions
Edit `testing_scripts/qc_metric_distributions.py`:
- `WES_DIR = 'testing_scripts/testing_data/WGS'`
- `CACHE = 'testing_scripts/testing_data/wgs_metric_values.json'` (don't clobber the WES cache)
- set `METRICS` to the WGS candidate list with (key, direction, unit), e.g.:
  ```
  ('MEDIAN_COVERAGE', 'min', 'x'), ('MEAN_COVERAGE', 'min', 'x'),
  ('PCT_20X', 'min', 'frac'), ('PCT_30X', 'min', 'frac'), ('PCT_10X', 'min', 'frac'),
  ('reads_mapped_percent', 'min', '%'), ('reads_properly_paired_percent', 'min', '%'),
  ('HET_SNP_SENSITIVITY', 'min', 'frac'),
  ('FREEMIX', 'max', 'frac'), ('reads_duplicated_percent', 'max', '%'),
  ('PCT_EXC_TOTAL', 'max', 'frac'), ('error_rate', 'max', 'frac'),
  ('reads_MQ0_percent', 'max', '%'),
  ```
  (Adjust to whatever the survey confirmed is present.)

Run it. Reads distributions across all 10 cohorts and writes the cache. Interpret the
percentile tables: for `min` metrics the left tail (low percentiles) are the bad
samples; for `max` metrics the right tail.

### Step 3 — flag rates for candidate thresholds
Edit `testing_scripts/qc_threshold_flagrates.py`:
- `CACHE = 'testing_scripts/testing_data/wgs_metric_values.json'`
- set `DIRECTION` to your WGS metrics
- set `CONFIGS['CURRENT']` to the existing genome config values, and add a
  `CONFIGS['PROPOSED']` with your data-driven candidates.

Run it. Choose thresholds so a **healthy cohort flags at ~0% fail / single-digit % warn**
while genuinely poor batches still surface. Iterate on the PROPOSED numbers and re-run
(fast — reads the cache, no re-parse).

### Step 4 (optional) — MAD for any protocol-dependent genome metric
If any genome metric shows strong cohort-dependent spread (bimodal medians across the
10 cohorts, like exome `ZERO_CVG` did), consider cohort-relative flagging:
- reuse `testing_scripts/mad_relative_prototype.py` (point it at the WGS cache, set the
  metric list) to check spread stability and cohort-growth churn (<1% is the bar);
- if it qualifies, add `[qc_thresholds.genome.relative.<METRIC>]` (same schema as the
  exome one) — the production `check_multiqc.relative_flags` already reads
  `['qc_thresholds', <seq_type>, 'relative']` for any seq_type, so **no code change is
  needed**, only config. Keep it warn-only, and keep an absolute `fail` hard gate.

### Step 5 — apply and verify
- Update `[qc_thresholds.genome.fail.*]` and `[qc_thresholds.genome.warn.*]` in
  `src/align_genotype/config_template.toml`, with an inline rationale comment per
  metric (cite the cohort percentiles / flag rates), matching the exome section's style.
- Validate the TOML parses: `python3 -c "import tomllib; tomllib.load(open('src/align_genotype/config_template.toml','rb'))"`.
- Run the suite: `uv run --extra test pytest test/ -q` (tests don't assert config
  values, so they should stay green; add tests only if you change code).
- Lint touched files: `uvx ruff check <files>` (baseline has 2 pre-existing `PLR0917`
  in check_multiqc — ignore those).
- Optional end-to-end sanity: adapt `testing_scripts/dryrun_exome_check.py` (set
  `sequencing_type='genome'`, mirror the genome config) to confirm flag counts on one
  real WGS JSON look sane.

## Guardrails / principles
- Healthy cohort ≈ 0% fail; fail = "do not analyse", warn = "a human should look".
- Preserve the lab's intent for hard gates unless the data clearly says otherwise; you
  can tighten `warn` freely.
- Don't over-flag on metrics that vary with ancestry/biology (error_rate, het
  sensitivity) — prefer warn, or cohort-relative, over hard fail.
- Parse one big JSON at a time; never load all 10 at once.

## Deliverable
Do this on a new branch (e.g. `genome-threshold-calibration`) off `QC_flags_report`,
commit the config change (+ any new tests), and open a PR — same pattern as PR #74.
Update `docs/qc_thresholds_tiers_plan.md`'s "genome tier not yet calibrated" note once
done.
