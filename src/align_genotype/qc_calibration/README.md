# QC threshold calibration

## What this is

Two CPG Flow stages that turn every dataset's latest CramMultiQC report into candidate
`fail`/`warn` thresholds and dataset-relative (MAD) tier evidence, rendered as an HTML
dashboard in the analysis dataset's web bucket. Run it when onboarding a new capture kit
or sequencing protocol, or when refreshing thresholds after the data underlying them has
moved.

- `QcCalibrationDatasetMetrics` (`qc_calibration_stages.py`) - a `DatasetStage`. One job
  per dataset: localises that dataset's latest CramMultiQC report and distils every
  configured metric to a small values file
  (`jobs/qc_calibration.py::extract_dataset_metrics`,
  `scripts/qc_calibration_extract.py`).
- `QcCalibrationReport` (`qc_calibration_stages.py`) - a `MultiCohortStage`, `forced=True`.
  One job: reads every dataset's values file and writes `calibration.json` and
  `calibration.html` (`jobs/qc_calibration.py::calibration_report`,
  `scripts/qc_calibration_report.py`).

Both stages are dependency-free - neither declares `required_stages` on a production
stage - so adding them to `run_workflow.py`'s `STAGES` list never pulls alignment or
genotyping work into a calibration run, and never adds a calibration job to a production
one.

## How to run it

```bash
analysis-runner --skip-repo-checkout \
    --image australia-southeast1-docker.pkg.dev/cpg-common/images/cpg-flow-align-genotype:<version> \
    --config src/align_genotype/config_template.toml \
    --config <calibration-run-config.toml> \
    --dataset <analysis-dataset> \
    --description 'QC threshold calibration' \
    --access-level full \
    --output-dir OUTPUT_DIR \
    run_workflow
```

The calibration-run config needs:

- `qc_calibration.enabled = true`
- `workflow.sequencing_type` - `'genome'` or `'exome'`
- `workflow.input_cohorts` - which datasets the run covers
- `workflow.only_stages = ['QcCalibrationDatasetMetrics', 'QcCalibrationReport']` - both
  stage names are required. `only_stages` marks every unlisted stage skipped, so naming
  just the report stage would make the extract stage check for outputs rather than
  produce them.

`enabled` belongs **only** in a calibration run's config, never in a config a production
run also loads. The two stages sit in the production `STAGES` list with no dependency on
anything else, so a truthy flag in a shared config would give every ordinary production
invocation a job per dataset (and register a Metamist analysis) for no reason.

One run covers one sequencing type - both stages read a single
`workflow.sequencing_type` - so calibrating genome and exome thresholds means two
separate invocations.

### Dry-run first

Add `--dry_run` to the `run_workflow` invocation and read the printed DAG before
committing real compute. Four things are worth checking, because none of them can be
verified without a live Metamist and a real multicohort:

1. **Discovery returns rows at all.** The Metamist `meta` filter is a flat mapping
   (`{'stage': 'CramMultiQC', 'sequencing_type': <type>}`), matching
   `scripts/build_vntyper_index.py`, the one other place in this repo that filters
   analyses on `meta`. `meta` is an opaque `JSON` scalar in the schema with no typed
   filter object, so nothing validates the convention offline - the unit tests use a fake
   query function. If the filter is wrong it returns nothing, and that reads as "no
   dataset has a report" rather than as a query bug.
2. **Only the two calibration stages queue**, and no alignment or genotyping stage is
   pulled in alongside them.
3. **Datasets without a CramMultiQC report are simply absent** from the DAG rather than
   erroring.
4. **A production config queues neither stage.** Run the same command with your ordinary
   config, which does not set `qc_calibration.enabled`, and confirm both stages are inert.

Run at `--access-level test` first. The invocation above uses `full` because the report
is written to the analysis dataset's main and web buckets.

## What it writes

| Path | Stage | What it is |
|---|---|---|
| `dataset.prefix()/qc_calibration/<seq_type>/values.<analysis_id>.json` | `QcCalibrationDatasetMetrics` | One dataset's extracted metric values, keyed on the CramMultiQC analysis ID so a new report produces a new path and re-extraction happens exactly when the underlying data changed. |
| `analysis_dataset.prefix()/qc_calibration/<seq_type>/calibration.json` | `QcCalibrationReport` | The full cross-dataset result: percentiles, flag rates, candidate thresholds, relative-tier evaluation, warnings. Everything the HTML shows, in one structure. |
| `analysis_dataset.web_prefix()/qc_calibration/<seq_type>/calibration.html` | `QcCalibrationReport` | The dashboard, registered as a `web` Metamist analysis. |

The `qc_thresholds` block this feeds lives in `config_template.toml`, one directory up.
Treat what the report proposes as **merge-carefully, not paste-over-the-top**. The
candidates carry per-threshold statistical evidence, but they cannot reproduce the
committed file's domain prose - which Picard tool a metric comes from
(`CollectWgsMetrics` vs `CollectHsMetrics`), the lab's framing of an "ideal ask" for a
given coverage or breadth figure. The committed thresholds' comments also record what
they were calibrated against - genome against 10 real WGS datasets, exome against two
real WES datasets of 647 and 522 sequencing groups - and a wholesale paste erases that
provenance along with the domain prose. Pasting a generated block over the top trades one
kind of documentation for another, not an upgrade. Diff the candidate against what's
already there, then merge by hand.

## How to read the report

The judgement is the part no number in the report supplies.

- Aim for a healthy dataset flagging roughly 0% fail and single-digit % warn. `fail`
  means "do not analyse without a decision"; `warn` means "a human should look, and
  usually proceeds with a note".
- Preserve the lab's intent for hard gates unless the data clearly contradicts it.
  Tighten `warn` freely - it only costs a human a look. Moving `fail` changes what gets
  analysed at all.
- The bar for overriding the lab is a specific measured count of good sequencing groups
  their line would have discarded, not a general sense that the line looks tight.
- Don't hard-fail on metrics that track ancestry, biology or chemistry rather than
  sample quality. `error_rate` and `HET_SNP_SENSITIVITY` vary for reasons that have
  nothing to do with sample quality; a `fail` on a metric like that fires on
  populations, not problems.
- Candidates are percentile tails, nothing more. Read them against the flag-rate tables
  and against what the lab actually asked for - never on their own.
- **Peak figures carry their own denominators.** The report states how many datasets
  each peak (warn rate, growth churn, merge churn) was computed over, because the
  filters differ per simulation: skipped datasets drop out of the warn-rate peak, growth
  additionally drops any dataset whose before-slice would itself be too small to
  threshold, and merge needs at least two usable datasets to produce a pair at all. A
  reader who sees only "peak merge churn 64.7%" cannot tell whether that describes one
  bad pair among ninety or the general case, so always read the count next to the
  percentage.

## When to adopt a dataset-relative tier

Two conditions, both required:

1. The spread is genuinely dataset- or protocol-dependent - look at the per-dataset
   medians the report shows, not at intuition about what "should" vary.
2. The flag set stays stable as the dataset grows.

Such a tier ships only `direction`, `k` and `min_samples` into
`qc_thresholds.<seq_type>.relative.<metric>` - production recomputes the actual
threshold from that run's own median and MAD every time it runs. The per-dataset
thresholds shown in the report are illustrative of the data you calibrated against; they
are not what gets shipped, and they will move.

The report simulates two different things and judges each against its own bar.
**Growth** - a dataset's leading slice re-scored against the threshold the whole dataset
produces - is a forecast: datasets accrete sequencing groups over time, and that's what
a shipped tier actually faces in production. **Merge** - two datasets pooled into one
run - is a stress test: nothing schedules two projects into one run, but it asks the
harsher question. Name the tension plainly: a metric earns a relative tier *because* its
normal level shifts between datasets, and the merge simulation punishes exactly that
property. A high merge figure may be intrinsic to the whole class of metric relative
tiers exist for, not evidence that one metric's tier is broken.

## Expect REJECT on both shipped tiers

Run over the full dataset set, this tool currently reports:

| tier | growth churn | merge churn |
|---|---|---|
| exome `ZERO_CVG_TARGETS_PCT` | ~1.0% | ~51.1% |
| genome `reads_duplicated_percent` | ~8.6% | ~64.7% |

Both tiers were adopted on a hand-picked subset, not the full set the figures above come
from. That gap is a **live QC question for the team to settle** - are these tiers
churnier than intended, or is gating on a cross-project merge the wrong test for this
class of metric - and not a tool failure or a sign the bars themselves are miscalibrated.

## Metric keys differ by sequencing type

Genome gates Picard `CollectWgsMetrics` output; exome gates Picard `CollectHsMetrics`
output, and the two lists are disjoint where the capture matters. Shared via samtools:
`reads_mapped_percent` (both sequencing types) and `reads_duplicated_percent` (both;
dataset-relative on genome). `reads_properly_paired_percent` is samtools-derived too but
is currently configured for genome only - exome simply has no properly-paired gate
today, not for a capture-kit reason. Contamination via verifybamid: `FREEMIX` (both).
Read the shipped `[qc_calibration.<seq_type>.metrics]` tables in `config_template.toml`
for the exact, current list rather than trusting a list to stay in sync with it - don't
copy an exome metric list onto a genome run, or vice versa.

One trap worth naming: `PCT_PF_READS_ALIGNED` is not in `report_general_stats_data` -
MultiQC writes it only to `report_saved_raw_data`, which the production check never
reads, so a metric configured under that key can never be gated. Use samtools
`reads_mapped_percent` instead.

## What the report will tell you loudly

Two checks worth knowing about, both catching the same class of defect - a gate that
looks configured and checks nothing:

- A metric absent from every dataset in the run.
- A warn tier that `fail` can never let fire - `warn <= fail` for a `min` metric,
  mirrored for a `max` one, because production evaluates `fail` before `warn` and
  records only the worst tier breached.

Both checks run against the *shipped* `config_template.toml` thresholds as well as
against this run's candidates, and catching either defect in the shipped thresholds
already running in production is reported as a live production defect, not filed
alongside routine candidate-review notes.

## Local iteration

Both job scripts are `python -m` runnable directly, so a downloaded set of values files
can be re-reported with a different `k` or metric list without re-parsing any MultiQC
report:

```bash
python -m align_genotype.scripts.qc_calibration_report \
    --values a.json --values b.json \
    --output-json out.json --output-html out.html
```

`--values` and `--skipped-dataset` are repeatable. The report reads
`workflow.sequencing_type` and `qc_calibration.<seq_type>.metrics`/`k`/`min_samples` the
same way the Hail Batch job does, through `cpg_utils.config` - point `CPG_CONFIG_PATH`
at a local copy of `config_template.toml` (plus any override) before running.

`scripts/qc_calibration_extract.py` is runnable the same way, against a downloaded
MultiQC report - see its `--help` for the full option list - which is the quickest way
to check what a new MultiQC version actually surfaces.

## Troubleshooting

**A metric shows `MISSING` in the presence matrix**
That metric key is absent from that dataset's general-stats sections entirely - usually
a MultiQC key rename between versions, or the key belongs to the other sequencing type.
Check the exact key against the shipped `[qc_calibration.<seq_type>.metrics]` list.

**A dataset appears under "Datasets not included"**
`QcCalibrationDatasetMetrics` found no completed CramMultiQC report for that dataset at
this sequencing type. Check that CramMultiQC has actually run for that dataset at that
sequencing type.

**A metric is present but its value list is empty**
A different problem from `MISSING`: the key exists in the report, but every value for it
was unusable. Picard writes `'?'` when coverage is effectively zero, and extraction
drops and counts those rather than raising. Check the dropped-value count next to the
metric - a key that's structurally absent needs a different fix (a rename, or the wrong
sequencing type) than one that's present but full of placeholders.
