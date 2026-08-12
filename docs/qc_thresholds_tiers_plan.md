# Plan: warn/fail QC threshold tiers

Status: **implemented**. Config (`config_template.toml :: qc_thresholds`, nested
`<seq_type>.<warn|fail>.<min|max>`), `QcFlag.severity`, the fail-before-warn
evaluation in `scripts/check_multiqc.py`, severity round-trip in
`scripts/record_qc_flags.py`, and the report display (`scripts/sg_qc_report.py` +
`templates/sg_qc_overview.html.jinja`) are all in place, with tests in
`test/test_check_multiqc.py` and `test/test_sg_qc_report.py`. The notes below are
retained as the design rationale.

A 647-sample WES dry run produced 1006 fail + 1021 warn flags — confirming the
tiers work end-to-end but also that the exome thresholds still need kit-specific
calibration (see the caveat in the config).

## Motivation

Today a metric has a single threshold, so every gate must sit at the
"catastrophic failure" line — we're blind to samples that are merely mediocre.
We want two tiers per metric:

- **fail** — do not analyse without a decision.
- **warn** — a human looks, usually proceeds with a note.

This mirrors the guidance in `docs/qc_flagging_update.md`. (Cohort-relative
outlier detection and identity/relatedness checks are explicitly **out of scope** —
relatedness lives in the separate relatedness pipeline.)

## Config schema (decided: nested severity)

Add a severity level between `<seq_type>` and `<min|max>`:

```toml
[qc_thresholds.exome.fail.min]
"MEAN_TARGET_COVERAGE" = 50
[qc_thresholds.exome.warn.min]
"MEAN_TARGET_COVERAGE" = 80

[qc_thresholds.exome.fail.max]
"FOLD_80_BASE_PENALTY" = 3.0
[qc_thresholds.exome.warn.max]
"FOLD_80_BASE_PENALTY" = 2.0
```

Rationale: minimal change to the existing `min`/`max` structure (adds one
dimension), keeps warn/fail visually grouped, and lets a metric define only the
tiers that make sense (e.g. warn-only metrics like `PCT_TARGET_BASES_50X`).

The live thresholds are in `config_template.toml :: qc_thresholds` with per-metric
rationale. The **exome** tier has since been calibrated against two real WES
cohorts (647- and 522-sample MultiQC runs) using
`testing_scripts/qc_metric_distributions.py` and `qc_threshold_flagrates.py`:
thresholds were chosen so a healthy cohort flags at ~0% fail / single-digit % warn
while a genuinely poor batch still surfaces (poor cohort: 242 fail / 542 warn;
healthy cohort: 4 fail / 116 warn). `ZERO_CVG_TARGETS_PCT` and `PCT_SELECTED_BASES`
proved strongly capture-kit dependent. Re-run those scripts on a known-good batch
when onboarding a new kit.

**Cohort-relative (MAD) flagging — implemented for `ZERO_CVG_TARGETS_PCT`** (PR #74,
branch `mad-relative-flagging` → `QC_flags_report`). Because its absolute level is
kit-dependent, its warn tier is now cohort-relative (median/MAD modified z-score,
warn-only) while the absolute `fail > 0.10` gate is kept as a hard stop. Validated
across 6 WES cohorts: warn rate 0–8%, <1% flag churn on cohort growth. See
`[qc_thresholds.exome.relative]` in the config and `check_multiqc.relative_flags`.
`PCT_SELECTED_BASES` / `PCT_OFF_BAIT` were evaluated and rejected (cohort-dependent
spread, up to 24.5% churn) — left ungated. Prototype/verification tooling:
`testing_scripts/mad_relative_prototype.py`, `verify_relative_6cohorts.py`.

The **genome** tier has now been cohort-calibrated against 10 real WGS MultiQC
cohorts (110-841 samples each), the same way as exome (branch
`genome-threshold-calibration`, based on `mad-relative-flagging` so the relative
code is present). WGS uses Picard `CollectWgsMetrics`, so the keys differ from
exome (no `*_TARGET_*` / `FOLD_80` / `ZERO_CVG`). Calibrated gates:

- `MEDIAN_COVERAGE` fail 15 / warn 25 (cohort medians ~32-37x; primary depth gate).
- `PCT_20X` fail 0.75 / warn 0.85 (primary breadth gate for WGS; medians ~0.94-0.95).
- `reads_mapped_percent` fail 80 (lab intent) / warn 97 (~99% median; only the messy
  chop-gliadx cohort fails ~5%).
- `reads_properly_paired_percent` warn 92, **warn-only** (correlated with mapped%,
  added as an extra human-look signal rather than a second hard gate) - a new metric.
- `FREEMIX` fail 0.04 / warn 0.02 (contamination ~0 across all cohorts; safety gate).
- `reads_duplicated_percent` absolute fail 40 **+ cohort-relative (MAD) warn** - dup
  is strongly cohort/library-prep dependent on WGS (medians 7.2-18.0%), so like exome
  `ZERO_CVG` its warn is cohort-relative (`[qc_thresholds.genome.relative]`, warn 0-4.2%
  per cohort, ~1% homogeneous / 2.1% worst-case churn). The lab's original fixed fail
  of 25 would have hard-failed ~25-37% of legitimately higher-dup preps, so it was
  relaxed to 40 (0% across all cohorts). Validated via the production `relative_flags()`
  path on all 10 cohorts (`testing_scripts/mad_relative_prototype.py`,
  `verify_relative_6cohorts.py`; end-to-end `dryrun_genome_check.py`).

## Code changes

1. **`utils.QcFlag`** — add `severity: str = 'fail'`. The default keeps existing
   DB-stored flags (which have no `severity`) loadable via `QcFlag(**flag)`, and
   means a plain re-run treats old flags as `fail` (their historical meaning).

2. **`scripts/check_multiqc.py`** — the threshold loop gains a severity dimension.
   Evaluate **fail before warn** and record **one** flag per `(section, metric)`
   at the *highest* tier tripped (a fail value also breaches warn; we don't want
   two flags for one metric). Sketch:

   ```python
   for section_name, section in sections.items():
       for sample, val_by_metric in section.items():
           for metric in metrics_seen(section):
               sev = worst_breach(metric, val_by_metric, seq_type)  # 'fail' | 'warn' | None
               if sev:
                   emit QcFlag(..., severity=sev, threshold=<that tier's threshold>)
   ```

   `worst_breach` checks the `fail` table first, then `warn`, for both `min` and
   `max`. Keep the `warn_unmatched_metrics` check (extend it to scan both tiers).

3. **`scripts/record_qc_flags.py`** — largely unchanged. Flags are keyed by
   `(section, flag)`, so a metric that escalates warn→fail (or de-escalates) keeps
   the same key; because the `threshold` value differs between tiers,
   `compare_qc_flag` returns False and the existing **"updated"** branch overwrites
   it with the new severity. Add `severity` to the fields carried through
   `reconcile_sg_qc_flags` (it's already `asdict`/`QcFlag(**flag)` round-tripped).
   Consider logging escalations explicitly.

4. **`scripts/sg_qc_report.py` + template** —
   - `_flag_to_dict`: pass `severity` through.
   - Template: amber `badge-warn` alongside red `badge-fail`; drive the
     `flag-line` left-border colour by severity.
   - `summarise_flags`: split active counts into `active_warn` / `active_fail`;
     the "Active flags" summary card shows both.
   - Add a severity filter chip next to the source/metric chips.

5. **`scripts/check_multiqc.py` Slack message** — prefix warn lines with ⚠ and
   fail lines with ❗, and summarise counts (e.g. "3 failing, 5 warnings"). Low
   priority — Slack is being retired in favour of the rendered report.

## Tests

- `test_check_multiqc.py`: metric at fail value → one `fail` flag; at warn value →
  one `warn` flag; warn-only metric never emits `fail`; a value breaching both
  tiers emits a single `fail` flag (dedup by severity).
- `test_sg_qc_report.py`: `summarise_flags` warn/fail split; template renders
  `badge-warn`; severity chip present.
- `record_qc_flags`: warn→fail escalation lands in the "updated" branch.

## Migration / backfill

Existing stored flags have no `severity`; the `= 'fail'` default makes them read
as fail until the next pipeline run re-evaluates and rewrites them with an explicit
tier. No manual DB migration required.
