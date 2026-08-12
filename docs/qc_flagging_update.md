User:
What are some very general and sensible QC flag thresholds for human genomics data? We typically work with WGS or WES data, aligned from fastq, bam, or cram data, typically given to us via sequencing providers. We re-align to hg38 and perform genotyping and various analysis through our pipelines. Before the analysis, we usually do some "QC Checks" - basically, collect heaps of QC metrics from samtools stats, picard, VerifyBamID; then we aggregate these using MultiQC, then we parse the raw multiqc data and check a few metrics against some thresholds. This has mostly been "vibes based" so far, I don't have strong justifications for the given thresholds. 

So far we just have a couple designed to flag issues in WGS samples:

```
[qc_thresholds.min]
# Picard - Median coverage must be >10x (obviously 30x is ideal)
"MEDIAN_COVERAGE" = 10
# samtools - Proportion of reads mapped to the reference must be >80%
"reads_mapped_percent" = 80

[qc_thresholds.max]
# VerifyBamID - Contamination must be <4%
"FREEMIX" = 0.04
# samtools - proportion of duplicate reads must be <25% 
"reads_duplicated_percent" = 25
```

A few questions

1. Are these sensible thresholds? They were determined by my boss who has a lot of experience both in the lab and working with data in silico, so I think his intuition is probably correct, but might need refining?
2. Can we apply these as is to WES data as well?
3. Are there any other core QC metrics that should be surfaced and flagged based on some threshold?

Agent:
This is a good set of questions, and the honest framing ("vibes based") is the right starting point — most published thresholds are also vibes, just older vibes. Your boss's numbers aren't wrong, but they're almost all set at the "catastrophic failure" level rather than the "something is off, look at this" level. That's a defensible design choice, but it means you're currently blind to samples that are merely mediocre.

## 1. Your four thresholds

**`MEDIAN_COVERAGE >= 10` — too low to be useful as your only coverage gate.**

At 10x, germline het sensitivity is genuinely bad: by binomial sampling alone a meaningful fraction of het sites get 0–2 alt reads, and callers then either miss them or emit them at low quality. Het sensitivity is roughly ~90% at 15x, ~95–97% at 20x, and plateaus around 30x. So 10x isn't "suboptimal but usable" — it's a different sensitivity regime.

Better: two tiers, plus switch your primary gate to coverage *breadth*, which captures depth and uniformity in one number:

- `PCT_20X >= 0.90` (warn), `>= 0.80` (fail)
- `MEDIAN_COVERAGE`: warn `< 25`, fail `< 15`

Two practical gotchas with Picard `CollectWgsMetrics`:
- If your hg38 has alt/decoy contigs, `GENOME_TERRITORY` includes them and dilutes your coverage numbers. Pass an `INTERVALS` list restricted to chr1–22 + chrX/chrY.
- It excludes duplicates, MAPQ<20 and BQ<20 by default, so your `MEDIAN_COVERAGE` is already effective post-filter coverage. Worth also surfacing `PCT_EXC_TOTAL` — if >20–25% of raw bases are being discarded, that's a signal in itself even when median coverage passes.

**`reads_mapped_percent >= 80` — far too lenient.**

Human WGS to hg38 with BWA-MEM/DRAGMAP routinely lands at 98–99.5%. WES similar. An 80% floor only catches wrong species, wrong reference build, or badly untrimmed adapters. Move to warn `< 98`, fail `< 95`. Note this metric is nearly useless for detecting human-on-human contamination — mapping rate stays high — so it's really a "did the right data reach the right reference" check.

**`FREEMIX <= 0.04` — permissive, and depends on your downstream use.**

For germline diagnostic work, contamination above ~1–2% starts producing spurious low-VAF calls that look like mosaicism, and inflates het counts. Common practice is 0.01–0.02 warn, 0.03–0.05 fail. gnomAD-style large cohorts used ~5% as the exclusion line; TOPMed and similar used tighter. I'd suggest warn `> 0.01`, fail `> 0.03`.

Two caveats: use VerifyBamID2 (the maintained version) with an ancestry-appropriate reference panel — misspecified panels inflate FREEMIX for non-European samples. And FREEMIX variance grows at low coverage, so don't over-interpret it on a 12x sample.

**`reads_duplicated_percent <= 25` — too lenient for WGS, too strict for WES.**

PCR-free WGS on patterned flowcells typically comes in at 1–8%. 25% would mean something went badly wrong with library input or amplification. Warn `> 10`, fail `> 20` is more informative for WGS.

Also: make sure MarkDuplicates ran with `OPTICAL_DUPLICATE_PIXEL_DISTANCE=2500` for patterned flowcells (NovaSeq etc.), otherwise ExAmp cluster duplicates get miscounted. And conceptually, duplicate rate is a *library complexity* metric, not a data quality metric — if post-dedup coverage is fine, high duplication cost you money but didn't hurt your calls. Consider treating it as warn-only and letting coverage breadth be the fail gate. `ESTIMATED_LIBRARY_SIZE` from Picard is the more direct measure of the underlying problem.

## 2. Can you apply these to WES?

Mostly no, and one of them is meaningless.

`MEDIAN_COVERAGE` from `CollectWgsMetrics` is computed over the whole genome, so on an exome it will be some small number that reflects off-target spill, not your actual data quality. You need `CollectHsMetrics` with bait and target interval lists:

| Metric | Warn | Fail | Note |
|---|---|---|---|
| `MEAN_TARGET_COVERAGE` | < 80 | < 50 | 100x is the usual ask |
| `PCT_TARGET_BASES_20X` | < 0.95 | < 0.90 | clinical exomes often demand ≥0.98 |
| `PCT_TARGET_BASES_50X` | < 0.85 | — | |
| `FOLD_80_BASE_PENALTY` | > 2.0 | > 3.0 | capture uniformity |
| `PCT_SELECTED_BASES` (on-bait) | < 0.60 | < 0.40 | strongly kit-dependent — calibrate to yours |
| `ZERO_CVG_TARGETS_PCT` | > 0.01 | > 0.05 | |
| `AT_DROPOUT` / `GC_DROPOUT` | > 5 | > 10 | |

`reads_mapped_percent` and `FREEMIX` thresholds carry over (use an exome-appropriate VerifyBamID2 marker set). Duplication does **not** carry over — 10–30% is normal for exomes given the small target space and higher depth, so a 25% max would flag healthy samples. Raise to warn `> 35`, fail `> 50`, or drop it as a gate and rely on post-dedup target coverage.

Practical implication: your threshold config wants to be keyed by assay type, not global. Something like `[wgs.min]`, `[wes.min]`, and later `[panel.min]`.

## 3. Metrics worth adding

I'd rank these by how often they catch real problems.

**Sample identity — the highest-value addition by far.** Coverage problems waste time; sample swaps produce confidently wrong answers.
- **Inferred vs reported sex**: normalised chrX and chrY coverage plus chrX heterozygosity. Catches swaps and mislabelled manifests constantly.
- **Fingerprint concordance** against an array, a prior sequencing run, or other samples from the same individual (`somalier`, Picard `CrosscheckFingerprints`, or `NGSCheckMate`). `somalier` is cheap, runs on BAM/CRAM directly, and also gives you relatedness and rough ancestry.
- **Within-cohort relatedness** (`somalier relate` or KING): catches accidental duplicates and pedigree errors.
- **Inferred ancestry**: use it as a cross-check on FREEMIX panel choice and to avoid flagging ancestry-driven variation as failure.

**Reference and provenance sanity** — cheap and catches provider mistakes:
- Verify `@SQ` M5 checksums in the BAM/CRAM header match your hg38. Providers do occasionally send hg19, or a different hg38 flavour (alts/no-alts, `chr` prefixes, HLA contigs). This is a five-second check that saves days.
- Verify `@RG SM` matches the expected sample ID.

**Paired-end and alignment integrity:**
- `reads_properly_paired_percent`: warn `< 95`, fail `< 90`
- `PCT_CHIMERAS` (Picard `AlignmentSummaryMetrics`): warn `> 0.03`, fail `> 0.05`
- `reads_MQ0_percent`: warn `> 5` — elevated MAPQ0 suggests contamination or reference mismatch
- `PCT_ADAPTER`: warn `> 0.005`

**Insert size** (`CollectInsertSizeMetrics`) — matters a lot if you do any SV/CNV calling:
- `MEDIAN_INSERT_SIZE`: warn outside ~200–600 (calibrate to your prep)
- High `MEDIAN_ABSOLUTE_DEVIATION` relative to median, or a visibly multimodal distribution, indicates library problems. Short inserts mean read-through and heavy overlap, which inflates `PCT_EXC_OVERLAP`.

**Error and base quality:**
- `error_rate` from samtools stats: typically 0.2–0.6% for good Illumina→hg38. Warn `> 0.01`. Note it's confounded by genuine divergence from the reference, so it's ancestry-sensitive — best used as a cohort outlier detector, not an absolute gate.
- Fraction of bases ≥Q30: warn `< 0.85`.

**Post-genotyping metrics** (if your QC step runs after joint calling, these are excellent global sanity checks):
- Variant count per sample: WGS roughly 4.0–5.5M total, heavily ancestry-dependent (African-ancestry genomes sit at the top of that range). WES roughly 20–25k coding variants.
- **Ti/Tv**: ~2.0–2.1 for WGS, ~3.0–3.3 for WES. This is a wonderfully sensitive indicator — a drop toward 1.5 usually means false positives from contamination or a systematic caller/filter problem.
- Het/hom ratio (~1.5–2.0) and inbreeding coefficient F. Strongly negative F (excess het) is a classic contamination signature and an independent check on FREEMIX.
- Insertion/deletion ratio (~0.8–1.0), novel (non-dbSNP) variant fraction, per-sample call rate/missingness.

## The structural fix for the "vibes" problem

Three changes matter more than any individual number:

**Two tiers, not one.** `warn` (human looks at it, usually proceed with a note) and `fail` (do not analyse without a decision). Your current single tier forces every threshold to sit at the catastrophe line, which is why they're all lenient.

**Add cohort-relative outlier detection alongside absolute thresholds.** For each batch, flag samples more than ~4–5 MADs from the batch median on each continuous metric. This catches failure modes no fixed threshold anticipates, and it self-calibrates to your assays and providers. Fixed thresholds catch known problems; robust z-scores catch unknown ones. You want both.

**Store every metric in a database, not just pass/fail.** After 100–200 samples per assay type you can set your `warn` thresholds empirically at the 1st/99th percentile of your own historically-good samples. At that point the thresholds stop being vibes and become "this sample is unlike the samples we know worked" — which is a much easier claim to defend, and much more useful than any number I or a paper can give you, because it's specific to your providers, kits, and pipelines.

One caution on the cohort-relative approach: don't normalise away real biology. Het rate, novel variant fraction, error rate and Ti/Tv all vary with ancestry, so a genuinely fine African-ancestry sample in a mostly-European batch will look like an outlier on several metrics. Either stratify by inferred ancestry or exempt those metrics from purely statistical flagging.

Finally, worth separating two things that often get conflated in one config: **contract compliance** (did the provider deliver the 30x / Q30 / on-target rate you paid for — useful for chasing re-sequencing) versus **analytical fitness** (can this sample support the analysis you're about to do). Same metrics, different thresholds, different actions.

I can write this up as a full stratified TOML config with WGS/WES sections, warn/fail tiers, and inline justifications for each number if that'd be useful as a starting point to edit.