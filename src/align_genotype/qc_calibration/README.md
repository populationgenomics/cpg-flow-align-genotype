# `qc_calibrate` - deriving QC thresholds from real cohorts

This package turns a set of dataset-level MultiQC reports into the
`[qc_thresholds.<seq_type>...]` block in
[`config_template.toml`](../config_template.toml). Run it when you are onboarding a new
capture kit or sequencing protocol, or refreshing thresholds after the data has moved.

It replaces a set of hand-edited local scripts. The point of formalising it was not
speed - it was that the *judgement* behind a threshold used to live only in the head of
whoever ran the scripts. The commands produce numbers; this README carries the judgement.
Read the whole "Choosing thresholds" and "When to adopt a cohort-relative tier" sections
before you sign anything off.

Every command is `uv run qc_calibrate <cmd>`; each one's options are in `--help` and are
not repeated here.

## The three artifacts

Everything lives in a `calibration/` working directory, which is **gitignored**
(see the entry at the end of [`.gitignore`](../../../.gitignore)). Keep it that way:

| Artifact | Written by | What it is |
|---|---|---|
| `manifest.<seq_type>.toml` | `discover` (hand-edited after) | Which datasets, and which MultiQC report per dataset (`uri`, `analysis_id`, `timestamp`). The durable record of what a calibration was derived from, so it stays reproducible after Metamist moves on. |
| `values.<seq_type>.json` | `collect` | The value cache: every metric's finite values per cohort, plus provenance (sample count, MultiQC version, general-stats shape, dropped-value count). Small. Everything after `collect` reads only this. |
| `spec.<seq_type>.toml` | you (seeded by `suggest`) | The calibration spec: what to gate, in which direction, at what value, and why. Thresholds live here rather than in Python constants - editing constants in place was the pain point in the old workflow. |

**None of these may ever be committed.** Manifest cohort labels and spec `rationale` text
name real CPG datasets. `emit.render` checks the rationale independently before printing
anything, and refuses to emit a block whose rationale names a cohort in the cache - but
that guard only covers the emitted block, not the artifacts themselves. The `docs/`
directory in this repo is untracked for the same reason.

## The workflow

```bash
# 1. Which datasets and reports. Metamist query; writes a manifest you then review.
uv run qc_calibrate discover --seq-type genome --output calibration/manifest.genome.toml

# 2. The only slow step. Parses every report once, surveys it, distils it to the cache.
uv run qc_calibrate collect --spec calibration/spec.genome.toml \
                            --manifest calibration/manifest.genome.toml

# 3. Where the data actually sits: p1..p99 per metric, per cohort.
uv run qc_calibrate distributions --spec calibration/spec.genome.toml

# 4. A first draft of fail/warn from the distribution tails, written back as unreviewed.
uv run qc_calibrate suggest --spec calibration/spec.genome.toml

# 5. What those candidates would flag, per cohort. Iterate here: edit the spec, re-run.
uv run qc_calibrate flagrates --spec calibration/spec.genome.toml

# 6. Only for a metric with a [metrics.<KEY>.relative] block: warn rate and churn.
uv run qc_calibrate mad --spec calibration/spec.genome.toml

# 7. The block to diff and paste into config_template.toml. Refuses if anything is unreviewed.
uv run qc_calibrate emit-config --spec calibration/spec.genome.toml \
                                --manifest calibration/manifest.genome.toml

# 8. Optional end-to-end proof: the real check_multiqc, on one real report, under this spec.
uv run qc_calibrate dryrun --spec calibration/spec.genome.toml \
                           --manifest calibration/manifest.genome.toml \
                           --cohort <label> --output-dir calibration/dryrun
```

Review the manifest before step 2. Dropping a bad cohort, pinning an older analysis or
substituting a local file are all normal, and `collect` is the step that costs you ten
minutes.

**`collect` is the only slow command.** It parses reports of tens to ~500 MB, one at a
time. Everything after it reads the cache and returns in under a second - that is what
makes the `flagrates` loop practical: change a number in the spec, re-run, look, change
it again. If you find yourself waiting, you are re-collecting when you didn't need to.

`emit-config` prints; it never edits `config_template.toml`. That file carries comments,
ordering and judgement no generator reproduces, and the paste-after-diff step is the
cheapest guard against a calibration run quietly rewriting a production gate. The emitted
block is ordered to match the committed one so the diff reads as a change, not a rewrite.

## Writing a spec

```toml
seq_type = "genome"
cache = "calibration/values.genome.json"

# A plain absolute gate. direction = 'min' means higher is better, so a value *below*
# the threshold is flagged; 'max' is the mirror. unit is 'x', '%' or 'frac' and only
# affects how values are rounded and displayed.
[metrics.MEDIAN_COVERAGE]
direction = "min"
unit = "x"
gated = true
fail = 15
warn = 25
reviewed = true
rationale = "Primary depth gate. 30x is the usual ask; fail below 15x catches the genuinely under-sequenced tail."

# Warn-only is legitimate: a metric can define just one tier.
[metrics.reads_properly_paired_percent]
direction = "min"
unit = "%"
gated = true
warn = 92
reviewed = true
rationale = "Correlated with reads_mapped_percent, so an extra human-look signal rather than a second hard gate."

# A cohort-relative warn tier. The absolute fail gate stays; there is no absolute warn.
[metrics.reads_duplicated_percent]
direction = "max"
unit = "%"
gated = true
fail = 40
reviewed = true
rationale = "Duplication is strongly library-prep dependent on WGS, so the warn tier is cohort-relative."
[metrics.reads_duplicated_percent.relative]
k = 3.5          # Iglewicz-Hoaglin outlier line; 3.5 is the standard
min_cohort = 50  # below this MAD is too noisy - production skips relative flagging

# A rejected candidate, kept on the record. Surveyed and profiled, never enforced.
[metrics.PCT_SELECTED_BASES]
direction = "min"
unit = "frac"
gated = false
reviewed = true
rationale = "Evaluated for a cohort-relative warn tier and rejected: cohort-dependent spread, churn up to 24.5%."
```

Metric keys and cohort labels are written as unquoted TOML table headers, so both must be
bare keys - letters, digits, underscore, hyphen. Anything else is rejected on load rather
than silently mangled on write.

### The four rules the loader enforces

1. **A gated metric needs at least one tier** (`fail`, `warn` or `relative`). A gated
   metric with no threshold is a gate that checks nothing - the exact failure mode this
   whole package exists to prevent.
2. **A `relative` block requires an absolute `fail` behind it.** Cohort-relative flagging
   is warn-only by construction: it finds outliers *within* a cohort, so a uniformly
   terrible cohort produces no flags at all. Without a hard absolute floor there is
   nothing to catch that.
3. **A `relative` block forbids an absolute `warn`.** The relative tier *is* the warn
   tier. Both would double-flag the same samples and disagree about which threshold was
   breached.
4. **An un-gated metric carries no thresholds.** `gated = false` means surveyed and
   profiled for the record, not enforced. A threshold sitting on an un-gated metric reads
   as active and isn't - so the loader makes you say which you meant.

### Keep your rejects

Leave a rejected candidate in the spec as `gated = false` with a rationale saying what you
measured and why you said no. It costs nothing (un-gated metrics are still surveyed and
still appear in `distributions`), and it stops the next operator re-litigating a decision
from scratch eighteen months from now. The three worth having on the record for exome are
`PCT_SELECTED_BASES`, `PCT_OFF_BAIT` and - the reason the survey is fatal -
`PCT_PF_READS_ALIGNED`.

## Metric keys differ by sequencing type

The candidate metric list is *not* shared between exome and genome, because the Picard
module differs. This is a spec parameter, not a constant, for exactly that reason.

- **Exome** - Picard `CollectHsMetrics`: `MEAN_TARGET_COVERAGE`,
  `PCT_TARGET_BASES_20X` / `PCT_TARGET_BASES_50X`, `FOLD_80_BASE_PENALTY`,
  `ZERO_CVG_TARGETS_PCT`, `PCT_SELECTED_BASES`, `PCT_OFF_BAIT`, `AT_DROPOUT` /
  `GC_DROPOUT`.
- **Genome** - Picard `CollectWgsMetrics`: `MEDIAN_COVERAGE`, `MEAN_COVERAGE`,
  `SD_COVERAGE` / `MAD_COVERAGE`, `PCT_1X` ... `PCT_100X` (breadth), `PCT_EXC_*`,
  `HET_SNP_SENSITIVITY`, `GENOME_TERRITORY`. There are **no** `*_TARGET_*`, `FOLD_80` or
  `ZERO_CVG` keys - don't copy an exome spec across.
- **Shared, via samtools**: `reads_mapped_percent`, `reads_duplicated_percent`,
  `reads_properly_paired_percent`, `reads_MQ0_percent`, `error_rate`.
- **Contamination, via verifybamid**: `FREEMIX`.

### The `PCT_PF_READS_ALIGNED` incident

`PCT_PF_READS_ALIGNED` is **not** in `report_general_stats_data`. MultiQC only writes it
to `report_saved_raw_data`, which the production check never reads, so it can never be
gated. The old genome reads-mapped gate was configured on it and therefore checked
nothing at all, for as long as it existed. Use samtools `reads_mapped_percent` instead.

That is why a gated metric missing from *any* cohort is a hard error in `collect` rather
than a warning: the cache is written but marked `complete = false`, `require_usable`
rejects it, and nothing downstream will run. A warning here would scroll past, and the
result of it scrolling past is a gate that silently protects nothing. The survey is not
optional.

## Choosing thresholds

`suggest` seeds `fail` from the worst per-cohort p1 (`min` metrics) or p99 (`max`
metrics), and `warn` from p5 / p95. Those are starting points, not answers. What follows
is what a percentile cannot supply.

- **Aim for a healthy cohort flagging roughly 0% fail and single-digit % warn.** `fail`
  means "do not analyse without a decision"; `warn` means "a human should look, and
  usually proceeds with a note". `flagrates` marks a cohort with `*` and the metric with
  `[REVIEW]` when the fail rate exceeds `FAIL_RATE_LIMIT` (2%) or the warn rate exceeds
  `WARN_RATE_LIMIT` (10%) - a prompt to look, not a rejection. Whether a cohort is
  healthy is your call; it is not a computable property.
- **Preserve the lab's intent for hard gates unless the data clearly contradicts it.** You
  can tighten `warn` freely - it costs a human a look. Moving a `fail` line changes what
  gets analysed at all.
- **The bar for overriding the lab is real evidence of harm.** Genome
  `reads_duplicated_percent` `fail` was *relaxed* from 25 to 40 because 25 would have
  hard-failed roughly a quarter to a third of legitimately higher-duplication library
  preps. That is the standard: a specific measured number of good samples the lab's line
  would have thrown away.
- **Don't hard-fail on metrics that track ancestry, biology or chemistry.** `error_rate`
  and `HET_SNP_SENSITIVITY` vary for reasons that are not sample quality. Warn on them,
  or go cohort-relative, or leave them un-gated - but a `fail` on them will fire on
  populations, not problems.
- **`suggest` cannot decide anything, and the interlock says so.** Everything it writes
  lands `reviewed = false`, and `emit-config` refuses to emit while any gated metric is
  unreviewed, naming all of them at once. **That gate is only worth something if you
  actually check the number before flipping the flag.** Read it against `flagrates`,
  against the distribution, and against what the lab asked for. Flipping `reviewed = true`
  to make the tool stop complaining converts the one structural guard between a raw
  percentile and production config into a formality.

`suggest` refreshes `rationale` only when it is empty or still carries its own
`Seeded: ` marker, so prose you write yourself survives a re-run. Write prose - it is
pasted verbatim into the committed config as the per-metric comment. Never put a cohort
label in it.

## When to adopt a cohort-relative tier

A cohort-relative (MAD / modified z-score) warn tier is right for a metric whose *normal
level* genuinely shifts by cohort or protocol - bimodal or wide-ranging medians across
your cohort set. Genome duplication rate is the canonical case: cohort medians span
roughly 7% to 18% depending on library prep, so any fixed warn line either floods the
high-duplication cohorts or never fires on the low-duplication ones. Don't reach for it
just because a fixed line is awkward to pick.

Two conditions, both required:

1. **The spread is genuinely cohort- or protocol-dependent** - look at the per-cohort
   medians in `distributions`, not at your intuition.
2. **The flag set stays stable as the cohort grows.** This is the important half. A
   relative threshold moves when the cohort's median moves, and every sample that changes
   flag status because of that is a spurious "updated" flag in the database, caused by
   nothing but cohort composition. A tier that churns is worse than no tier at all.

The advisory bar in `relative.py` is three constants: `MAX_WARN_RATE = 0.10` (peak
per-cohort warn rate), `MAX_GROWTH_CHURN = 0.02` and `MAX_MERGE_CHURN = 0.05`. Clearing
all three prints `RECOMMEND`; missing any prints `REJECT` naming which bar and by how
much. It is advice for you to sign off, not an automatic gate, and `mad` never changes
the spec.

`mad` runs two simulations, and they are judged against **separate bars** because they
model different things. **Cohort growth** — a 60% before-slice of each cohort, re-scored
against the threshold the full cohort produces — is a *forecast*: samples get added to a
project over time, and that is what a shipped tier actually faces. **Cohort merge** —
each cohort re-scored against the threshold it gets once a second cohort joins it — is a
*stress test*: nothing schedules two projects into one run. So growth is held to the
strict 2% and merge to a looser 5%.

The merge bar is looser, not absent. A metric churning a quarter of its flag set on a
merge is unstable however you frame the scenario, and dropping merge from the verdict
entirely would have admitted the two exome metrics recorded as rejected below.

Both bars are **conservative defaults you may revisit with evidence**, not values derived
from a shipped decision. Read the record below before treating either as settled.

The merge simulation runs every *ordered* pair — `n*(n-1)`, so 90 at ten cohorts, not 45
— because a merge disturbs both projects' flag sets and the two directions measure
differently. Note the structural tension: a metric earns a relative tier precisely
because its normal level shifts between cohorts, and the merge simulation punishes
exactly that property. A high merge figure may be intrinsic to the whole class of metric
relative tiers exist for, rather than evidence that one metric is broken.

The growth simulation is run twice per cohort, on the leading 60% and on a
seeded-shuffled 60%, and the verdict uses the **worse** of the two. That is not paranoia,
and the ordered reading is not noise: across ten real WGS cohorts, four show the leading
60% and trailing 40% of the report differing in median duplication by 2.6 to 9.5
percentage points — one goes 14.8% to 24.4%. Samples later in a MultiQC report have
systematically higher duplication, which is what batch-ordered sequencing looks like. So
the leading slice is a genuine forecast of batch growth, and an `ORDERING-SENSITIVE`
marker means the ordered and shuffled readings disagree about whether that cohort clears
the bar — treat it as the batch-growth reading being the one that matters, not as an
artifact to discount.

### The record so far

- **Adopted**: exome `ZERO_CVG_TARGETS_PCT` (kit-dependent zero-coverage rate) and genome
  `reads_duplicated_percent` (library-prep dependent, cohort medians 7.2-18.0%). Both
  keep their absolute `fail` gate. Warn rates are comfortable: exome 0-8.2% per cohort,
  genome 0-4.2%.
- **Rejected**: exome `PCT_SELECTED_BASES` and `PCT_OFF_BAIT` - cohort-dependent spread
  and churn up to 24.5% of the flag set. Left un-gated, with the rejection recorded in
  the spec.

**Read this before trusting either adoption.** Both tiers were adopted on a hand-picked
subset — growth measured on 3 of 10 cohorts (leading slice only, no shuffle), merge on 2
ordered pairs — which gave ~1% growth and ~2.1% merge. Run over the full cohort set, this
tool reports:

| tier | growth churn | merge churn | verdict |
|---|---|---|---|
| exome `ZERO_CVG_TARGETS_PCT` | 1.0% | 51.1% | `REJECT` (merge only) |
| genome `reads_duplicated_percent` | 8.6% | 64.7% | `REJECT` (both bars) |

On the same three cohorts the original analysis used, this tool reproduces its 1.1% /
0.0% / 0.0% exactly — so the arithmetic agrees and the gap is sampling, not a defect. The
8.6% comes from a cohort the manual analysis never examined, and it is the batch-ordering
effect described above rather than an artifact.

So `mad` will print `REJECT` for both currently-shipped tiers. That is a **live QC
question** — are these tiers churnier than intended, or is gating on a cross-project
merge the wrong test for this class of metric? — and not evidence that the bars are
miscalibrated. Whichever way it resolves, record the reasoning in the metric's
`rationale`.

## Memory

`collect` is the only command that touches the large reports. It parses one at a time and
releases it (with an explicit `gc.collect()` between cohorts) precisely so a ten-cohort
set of ~500 MB files fits in an ordinary machine. Everything else in the package works
from the cache.

If you find yourself wanting several reports in memory at once, you want the cache
instead - that is what it is for. `dryrun` also loads one full report, by design, since
its whole purpose is to run the production check against real data; it reports peak RSS
so you can answer "will this fit in a 4 GB job?".

## Troubleshooting

**"Cache is incomplete: a gated metric was missing from at least one cohort"**
`collect` found a gated metric absent from a cohort's general stats, so it marked the
cache `complete = false` and nothing will run on it. Read the survey's `MISSING GATED
METRICS` block: either the key is wrong for this `seq_type` (see "Metric keys differ"),
or that cohort genuinely lacks the metric and belongs out of the manifest, or the metric
should be `gated = false`. Fix one of those, then re-run `collect`.

**"Cache is missing [...] - it was collected for a different metric list"**
You added a metric to the spec after collecting. Narrowing the spec's metric list is free
(a subset of the cache is fine); adding one means the reports have to be parsed again.
Re-run `collect`.

**`emit-config` refuses: "unreviewed gated metric(s) [...]"**
Working as intended. Check each named metric's numbers against `flagrates` (and `mad`, if
it has a relative tier), then set `reviewed = true` on it in the spec. Don't flip them all
to silence the error - see "Choosing thresholds".

**`emit-config` refuses: "rationale names cohort ..."**
A `rationale` contains a cohort label, and the emitted block is pasted into a file that is
committed and pushed publicly. Rewrite it in aggregate terms - "cohort medians 7-18% across
10 cohorts" - and re-run.

**A cohort you expected is missing from the manifest after `discover`**
`discover` keeps only projects whose Metamist `meta.is_seqr` is true and whose name
contains none of `test`, `training`, `seqr` (substring, not token - so a name containing
any of those anywhere is dropped). Then it needs a *completed* `qc` analysis whose
`meta.sequencing_type` matches `--seq-type`, which has an `output`, and whose
`timestampCompleted` is a usable string; the newest such analysis wins, and skips are
logged at WARNING/INFO. Labels that aren't bare TOML keys are dropped too. If none of
that explains it, add the cohort to the manifest by hand - the manifest is meant to be
edited, and `uri` is anything `cpg_utils.to_path` accepts, including a local file.

**A metric shows as `MISSING` in the survey's presence matrix**
Usually a MultiQC key rename between versions. The survey's missing-metrics block dumps
every key present in each general-stats section of the affected cohorts for exactly this
- find the new name there and update the spec. Note the distinction: an empty presence
cell means the key is absent or renamed, while a metric that *is* present but whose value
list is empty means the key exists and every value was unusable (Picard's `'?'`
placeholder, typically - the survey's per-metric dropped counts will show it). Those need
different fixes.

**A cohort in the `UNREADABLE COHORTS` block**
Different problem: the report didn't parse at all, so that cohort is absent from every
table above. Check the URI, the file and your credentials, then re-run. One bad cohort
doesn't discard the reports already parsed, but the run still exits non-zero.

**`mad` prints "There are no metrics with a `[metrics.<KEY>.relative]` block"**
Not an error. A relative tier is a candidate you propose, not something the tool finds.
Add a `relative` block (with an absolute `fail` behind it) and re-run.
