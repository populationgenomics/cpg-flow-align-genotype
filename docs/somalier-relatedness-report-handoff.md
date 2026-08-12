# Handoff: building a Somalier relatedness report (from the SG QC report)

**Audience:** an agent building a Somalier/relatedness flag report, modelled on the
SG QC report in this repo. You'll likely work in `cpg-flow-relatedness`, but the
report pattern to copy lives here in `cpg-flow-align-genotype`.

**What this doc is:** everything we learned building the QC report — the pattern,
the reusable pieces, the gotchas that cost us time, and (most importantly) the
places where the relatedness report will **need to diverge** because Somalier
flags have a fundamentally different shape.

---

## 1. The pattern in one paragraph

A metric-checking script (`check_multiqc.py`) compares results against config
thresholds and emits structured flag objects. A recording script
(`record_qc_flags.py`) reconciles those flags into each sequencing group's
`meta` blob in Metamist, tracking a resolved/unresolved lifecycle. A report
script (`sg_qc_report.py`) queries every SG in a dataset, keeps the ones with
flags, enriches them with sample/participant/family/reads metadata, and renders
a single self-contained HTML dashboard (`templates/sg_qc_overview.html.jinja`).
The relatedness report is the **third** step for a different flag type — the
first two already exist in `cpg-flow-relatedness`
(`src/rd_qc/scripts/record_somalier_flags.py`).

---

## 2. Reference files (read these first)

| File | Role |
|---|---|
| `src/align_genotype/utils.py` | `QcFlag` frozen dataclass — the flag schema |
| `src/align_genotype/scripts/check_multiqc.py` | Produces flags from MultiQC JSON vs `qc_thresholds` config |
| `src/align_genotype/scripts/record_qc_flags.py` | Reconciles flags into `SG.meta` (the lifecycle) |
| `src/align_genotype/scripts/sg_qc_report.py` | **The report — your main template** |
| `src/align_genotype/templates/sg_qc_overview.html.jinja` | The HTML/CSS/JS |
| `test/test_sg_qc_report.py` | Offline tests — the model for how to test a report |
| `testing_scripts/test_qc_flags_e2e.py` | End-to-end smoke test against a live SG |

In `cpg-flow-relatedness`: `src/rd_qc/scripts/record_somalier_flags.py` and
`src/rd_qc/utils.py` (the three `Somalier*Flag` dataclasses).

---

## 3. The flag lifecycle (shared by both report types)

Flags are stored **in the SG's `meta`** as a JSON list and reconciled on each run.
`record_qc_flags.py :: reconcile_sg_qc_flags` and
`record_somalier_flags.py :: reconcile_sg_somalier_flags` are the same algorithm:

- **Flag identity** = a subset of fields *excluding the measured value* (values
  drift between runs). QC uses `(section, flag, comparison, threshold)`. Somalier
  uses category-specific keys (see §7).
- Each existing flag is: **resolved** (was present, now absent → set
  `resolved=True` + `resolution_date`), **retained** (still present, refresh the
  measured value only), **updated** (present but identity changed, or a resolved
  flag reappeared → overwrite, `resolved=False`), or left alone (already resolved
  and still absent).
- New flags with no matching identity are **added**.
- A flag therefore carries `resolved: bool` and `resolution_date: str | None`.
  **The report's entire active-vs-resolved split keys off `resolved`.**

**Two meta caveats that matter for the report:**

1. The mutation sends only the flag key (`{qc_flags_key: [...]}` /
   `{'somalier_flags': [...]}`). Metamist **merges at the top-level meta key**, so
   other meta keys survive (the e2e test asserts `test_field` is preserved). This
   is why `cram_qc_flags`, `gvcf_qc_flags`, and `somalier_flags` coexist happily.
2. **Within** the `somalier_flags` list the *whole list is overwritten* each run.
   `record_somalier_flags.py` reconciles all three categories unconditionally for
   this reason — if you skip a category, its flags vanish instead of resolving.
   (Not your problem for the report, but explains why all categories always
   appear together in one list.)

---

## 4. Report architecture (what to copy)

`sg_qc_report.py` is ~430 lines and cleanly separated. The shape to reuse:

1. **`collect_qc_flags(sequencing_groups)`** — pull the flag list(s) out of each
   SG's `meta`, instantiate dataclasses. *For Somalier: read the single
   `meta['somalier_flags']` list; you'll partition by `category` yourself.*
2. **`get_sg_infos(sg_ids) -> dict[str, SGInfo]`** — one Metamist query for
   sample/participant/family/reads metadata, **keyed by sg_id** (see §6, ordering
   bug). Reusable almost verbatim.
3. **`_flag_to_dict(flag, ...)`** — flatten a dataclass into a template-ready dict
   with display-ready fields (`metric_label`, `value_display`, `date_short`, …).
   *This is where per-category display logic goes for Somalier.*
4. **`build_sections(reports) -> (unresolved, resolved)`** — split each SG's flags
   into active/resolved, build row dicts, sort. Two-section page structure.
5. **`summarise_flags` / `metric_histogram` / `source_histogram`** — header-card
   counts and filter-chip counts.
6. **`render_report(...)`** — build context, render Jinja. Keep it thin.
7. **`main(dataset, output)`** — query all SGs → collect → filter to flagged →
   enrich → render → write. Phase-timed logging (see §6, transient hangs).

The template (`sg_qc_overview.html.jinja`) is one self-contained file (inline CSS
+ JS, no external assets — important, these get shared as static files / emailed).
Key pieces: a Jinja **macro** `flag_table(rows, resolved)` renders both sections;
a **filter bar** (search + source chips + metric chips) driven by `data-*`
attributes and ~40 lines of vanilla JS.

---

## 5. Design decisions we made (and why) — likely to carry over

These were deliberate choices validated with the user; start from them:

- **Two top-level sections, not a status column:** `⚠ Unresolved flags` (the
  headline) and `✓ Resolved — past incidents` (greyed backlog at 72% opacity). An
  entity with both appears in both, showing only the relevant flags in each. This
  fixed an incoherence where the rendered set (SGs with *any* flag) didn't match
  the "flagged" count (SGs with *active* flags).
- **Flag-centric header cards, not SG-cleanliness:** `SGs scanned` / `Active
  flags` (+ breakdown subtext) / `SGs affected` / `Resolved`. **Colour semantics:
  red only when active > 0; green only for the all-clear state; grey for resolved
  and neutral counts. Never green a flag count** (a flag isn't "good"). When zero
  active, show a green "✓ All clear" banner instead of an empty table.
- **At-a-glance inline flags:** show up to 2 flags *in the top-level row* (badge +
  bold label + result), with "+N more — click to expand". Users should not have to
  click every row to see what's wrong. Full detail (all flags + metadata) stays
  behind the row expand. Threshold: ≤2 inline, >2 collapses.
- **Collaborator-facing identifiers:** lead with **family** (bold/dark) then
  **participant** (lighter); sample external ID secondary; **CPG SG id demoted to
  small muted monospace** (kept for internal reference, not prominent). Do *not*
  prefix with the words "Family"/"Sample" — rely on styling/position. Sort rows by
  family → participant → id.
- **Client-side filter bar:** metric chips (with per-metric SG counts), a
  source/category toggle, and a free-text search over
  family/participant/sample/id, plus a live "Showing X of Y" counter. All done
  with `data-*` attributes — no re-query, works in the static file. Only render
  the bar when there's something to filter (>1 metric, >1 source, or >5 rows).
- **Number formatting** (`_fmt_num`): integers stay integers; values ≥1 → 2
  decimal places; values <1 → 2 significant figures; strip trailing zeros. e.g.
  `64.579124 → 64.58`, `0.0616722 → 0.062`. Raw floats from the tools are ugly.
- **Value-vs-threshold in plain language** using the comparison direction:
  `22× (below minimum 30×)`, `0.12 (above maximum 0.04)`.
- **Dates truncated** to `YYYY-MM-DD` for display, full ISO timestamp kept in a
  hover `title`.
- **Human labels via a hand-maintained map** (`METRIC_LABELS`, `SECTION_LABELS`).
  The raw tool keys (`reads_mapped_percent`, `picard_4`) aren't user-friendly and
  the friendly names aren't in the source JSON. Small, stable set — just maintain
  a dict, with a fallback to the raw key.

---

## 6. Metamist gotchas (these cost us time)

- **External IDs are keyed under the empty string.** `sample['externalIds']` is
  `{'': 'HG003_NA24149'}`. Use a helper: `ext.get('') or next(iter(ext.values()),
  '')`. Same for participant and family. Null-guard the whole chain
  (`sample → participant → families[0]`) — any of them can be missing.
- **`sequencingGroups(id: {in_: [...]})` (plural) not `sequencingGroup(id:)`.**
  The plural form with an `in_` filter lets you fetch many SGs in one query. The
  response order is **not guaranteed** — build a `dict[sg_id -> info]` and look up
  by id; do **not** `zip()` the response against your request list (we hit this
  bug).
- **`meta` is a JSON scalar** — you get the whole blob, you can't sub-select
  fields inside it in GraphQL.
- **Assay `reads`**: `assay.meta['reads']` is a *list* of file dicts
  (`{basename, location, size, ...}`); `assay.meta['reads_type']` is
  `'fastq' | 'bam' | 'cram'`. A fastq assay is `[R1, R2]`. Group by type and pair
  R1/R2; one stress-test SG (`CPG276402`, validation-test) has ~76 fastq pairs, so
  **collapse long read lists** behind a scroll/expand.
- **Transient GraphQL hangs are real and recurring.** A query occasionally stalls
  for minutes. Root cause: `metamist.graphql.query`'s backoff only retries
  `HTTPError`/`JSONDecodeError`/`TransportServerError`; a stalled socket with no
  timeout raises none of those, so it blocks forever. Mitigations we applied:
  **phase + timing logs** around each query (so a hang is instantly attributable
  to the query, not your code) — log the SG count/ids before and elapsed after.
  If it becomes a real problem, set a socket/transport timeout so it fails fast.
  When it stalls during dev, just retry; it's not your code.
- **Auth / running locally:** `export SM_ENVIRONMENT=production` to hit the real
  API (default env otherwise points at localhost:8000 and refuses). The query
  itself is fast (~1s for 6 SGs; ~45KB for 76 assays' meta) — if it's slow it's
  the transient stall, not payload size.

---

## 7. How Somalier flags DIFFER — the important part

The QC report assumed: **one flag = one SG, one metric vs one threshold, split
across two meta keys (`cram_`/`gvcf_`).** Almost none of that holds for Somalier.
Plan for these differences up front:

### 7a. One list, three heterogeneous schemas, discriminated by `category`
Somalier flags all live in `meta['somalier_flags']` as a single list. Each item
has a `category` field:

| category | dataclass | scope | identity fields | measured fields |
|---|---|---|---|---|
| `sex_inference_mismatch` | `SomalierSexInferenceFlag` | **per-SG** | `provided`, `inferred` | `mean_depth`, `x_het_ratio`, `x_depth_ratio`, `y_depth_ratio`, `p_middling_ab` |
| `self_relatedness_mismatch` | `SomalierSelfRelatednessFlag` | **pairwise** (same participant, 2 samples) | `sg_id_1`, `sg_id_2`, `participant_external_id`, `threshold` | `relatedness`, `ibs0`, `ibs2` |
| `relatedness_mismatch` | `SomalierRelatednessFlag` | **pairwise** (across a family / dataset) | `sg_id_1`, `sg_id_2`, `family_external_id`, `expected_relationship`, `inferred_relationship` | `relatedness`, `ibs0`, `ibs2` |

So `collect_*` reads one key and partitions by `category`. `_flag_to_dict` needs a
**per-category branch** for labels and the "result" string — there is no single
`value/comparison/threshold` triple. e.g.:
- sex: `provided ♀ vs inferred ♂` (no numeric direction — the "value_display"
  concept is a provided-vs-inferred comparison, not below/above a threshold).
- self-relatedness: `relatedness 0.62 (expected ~1.0 for same individual)`.
- relatedness: `expected parent–child (~0.5), inferred unrelated (0.02)`.

### 7b. Pairwise flags break the per-SG row model
This is the biggest design decision. A relatedness flag is about a **pair** of
samples, and the recording script writes it into an SG's meta keyed by `sg_id`
(often on *both* members). Consequences:

- **The same logical flag can appear on two SGs' meta** → naive per-SG rendering
  **double-counts and duplicates** it. You need to dedup pairwise flags by their
  identity tuple before rendering.
- **What is a "row"?** The QC report's row = one SG. For relatedness that's
  awkward. Options to raise with the user (this is a genuine design fork worth an
  `AskUserQuestion`):
  1. **Per-pair rows** in the relatedness section (row = `sample_A ↔ sample_B`,
     showing both identities), and **per-SG rows** only for sex inference.
  2. **Per-family grouping** for the relatedness section (a family is the natural
     unit — "does inferred pedigree match expected?").
  3. Keep per-SG but list the pairwise partner inline.
  My recommendation: **three sections keyed by category** (Sex inference /
  Self-relatedness / Relatedness), each with the row unit that suits it
  (per-SG for sex, per-pair for the two relatedness kinds), each still split into
  active/resolved. The "source" facet (CRAM/GVCF) becomes the **category** facet
  for filters.

### 7c. Identifiers
Family-primary already fits relatedness well. But pairwise rows need **two**
identifiers — render both samples (and, for `relatedness_mismatch`, the
expected-vs-inferred relationship prominently, since that's the whole point).
`get_sg_infos` currently returns one SGInfo per id; for pairwise flags you'll want
to look up **both** `sg_id_1` and `sg_id_2` (make sure both are in the id list you
query).

### 7d. Filters / summary
- Filter facet: **category** (sex / self / relatedness) instead of source
  (CRAM/GVCF). Metric chips could become "flag type" chips.
- Header cards: count by category rather than CRAM/GVCF.
- Search should match both members of a pair.

---

## 8. Testing & tooling conventions

- **Test the report offline** — no DB. Build fixtures that mirror the exact DB
  `meta` shape, feed them through `collect_* → build_sections → render_report`,
  and assert on the builder output *and* on substrings in the rendered HTML
  (`test/test_sg_qc_report.py` is the template; 25 tests, all pure functions).
  Do the same for Somalier — a fixture per category, plus a pairwise-dedup test.
- **Run tests:** `uv run --with pytest python -m pytest test/ -q`
  (pytest isn't in the base env; `uv run --with pytest` injects it). `testpaths`
  is `test/` in `pyproject.toml`.
- **Lint:** `uv run --with ruff ruff check <files>`. Conventions we set:
  - Tests are exempted from `S101, PLR2004, RUF001, RUF059` via
    `[tool.ruff.lint.per-file-ignores]` `"test/**"` — add the same block if the
    relatedness repo doesn't have it.
  - `×` (U+00D7) trips `RUF001`; we used `# noqa: RUF001` on those dict lines.
- **Jinja gotcha:** a `{% macro %}` must be **defined before it's called** in the
  template (define it right after `<body>`). Autoescape is **on** — good, but
  means you build display strings in Python, not with `|safe`.
- **Verify visually with an offline render**, not the live DB — write a tiny
  script that feeds fixtures through `render_report` and writes to `/tmp/*.html`.
  Deterministic, and dodges the transient Metamist hangs.

---

## 9. Suggested order of work

1. Read `record_somalier_flags.py` + the three `Somalier*Flag` dataclasses; nail
   down exactly which fields each category carries and how pairwise flags are keyed
   into `meta` (on one SG or both).
2. Decide the **row/section model** with the user (§7b) before writing template
   code — it drives everything.
3. Fork `sg_qc_report.py`: keep `get_sg_infos`, the two-section split, header
   cards, filter bar, `_fmt_num`, date/label helpers. Rewrite `collect_*`,
   `_flag_to_dict` (per-category), and the row-building to the chosen model.
   Add **pairwise dedup**.
4. Fork the template: reuse the CSS/JS wholesale; adapt the macro to render
   per-category result strings and pairwise identities; swap the source facet for
   category.
5. Write offline tests mirroring each category + a dedup test. Render a fixture to
   `/tmp` and eyeball. Only then try a live run with `SM_ENVIRONMENT=production`.

---

## 10. Commands

```bash
# Run the existing QC report against a real dataset (needs metamist auth)
SM_ENVIRONMENT=production python -m align_genotype.scripts.sg_qc_report \
    --dataset validation-test --output /tmp/report.html

# Tests + lint
uv run --with pytest python -m pytest test/ -q
uv run --with ruff ruff check src/ test/
```

`validation-test` / SG `CPG276402` is the go-to real example (its QC flags are all
*resolved*, so it exercises the resolved section + all-clear banner; it also has
~76 fastq pairs for the read-collapse path). It already has a `somalier_flags` key
in its meta.
