# Gaps and friction — pathway enrichment B2

Live document — appended to at each step's decide phase. Sections per artifacts-guide §gaps_and_friction.md format. Pre-seeded items will be added by Task 12 Step 5 (from `docs/superpowers/specs/2026-04-18-research-methodology-v3-improvements-from-b2.md` §4); this file records the execution-time delta.

---

## KG data bugs

_(none observed yet)_

---

## KG gaps

_(none observed yet)_

---

## MCP friction

_(pre-seeded items from v3 meta-doc §4.5 A1–A4 will be added in Task 12 Step 5; execution-time items below)_

---

## Skill / methodology friction

### 2026-04-20 — Step 1a: publication attribution drift (anti-hallucination / Rule 7)

**What happened.** During MCP orientation in Step 1a, I presented the proposed NC1 experiment (`10.1038/ismej.2016.70_coculture_alteromonas_hot1a3_med4_rnaseq`) as "Hennon coculture MED4." The researcher corrected: the DOI belongs to **Aharonovich & Sher 2016** (ISME J, "Transcriptional response of Prochlorococcus to co-culture with a marine Alteromonas: differences between strains and the involvement of putative infochemicals"), not Hennon. Verification via `list_publications(organism="MED4", verbose=true)` confirmed the correct attribution.

**Root cause.** I attributed the DOI from intrinsic knowledge without running `list_publications` first. The `list_experiments` response carries `publication_doi` but not author/title — paper attribution requires a separate `list_publications` call. I skipped that call and fell back on training-era memory, which produced a silent mis-attribution.

**Scope of drift.** Only one item was checked when the researcher caught this (NC1). Running `list_publications` across all proposed experiments afterward revealed the other attributions were correct, but the check happened only because of the researcher's catch, not because I did it upfront.

**Rule violation.** Rule 7 (anti-hallucination) + Rule 1 (KG as sole data source for data). Paper-to-DOI mapping is a data claim, not interpretation, and must come from the KG.

**Proposed skill change (candidate for v3).**
- Add to `kg-rules.md` or `anti-hallucination.md`: **"Before presenting publication attributions for any experiment_id, run `list_publications(...)` or look up `publication_doi` against the publications table. Author names from training-era memory are not acceptable — DOIs are stable, author lists in the KG are canonical."**
- Concretely: every Step-1a-style experiment-classification flow should include a publication-resolution sub-step that joins `experiment_id → publication_doi → author list` via `list_publications`, committed alongside `experiments_classified.csv` (e.g., as an extra column or a sibling CSV). This makes the attribution check a hard data dependency, not a soft reminder.

**How to apply.** Before any future experiment-classification presentation, call `list_publications` with the relevant organism filter (or by DOI if the tool adds a `publication_doi` parameter — currently uses `author` / `organism` / `search_text` / DOI filter via experiment's `publication_doi`). Do not name authors from memory.

### 2026-04-20 — Step 1a: gene_count misreported as cumulative instead of per-timepoint

**What happened.** In the show-phase table, I reported "Tolonen 2006 R1 gene_count range 10,182 (min) to 10,182 (max)" for the R class — taking the top-level `gene_count` field from `list_experiments`. This is the *cumulative sum* across timepoints (6 TPs × 1697 genes = 10,182), not unique gene count. The researcher caught it: "gene count - report per timepoint not total - this is misleading - mean or range per timepoint."

**Root cause.** `list_experiments` returns both a top-level `gene_count` (summed across timepoints) and per-timepoint `tp_gene_count` (via `experiments_to_dataframe`). I used the top-level field without thinking about its aggregation. The KG schema is not wrong — the top-level number is legitimately the cumulative count — but presenting it as "genes per experiment" in a per-class summary table is misleading because it scales with #timepoints, not with detection power.

**Impact on reasoning.** A reader seeing "Tolonen 10,182 genes" assumes ~10k unique genes detected — 6× the actual MED4 ORFome. This distorts expectations about pathway-background size, per-cluster Fisher 2×2 dimensionality, and enrichment detection power. For Tolonen specifically it's obvious on second look (MED4 only has ~1,700 ORFs), but for larger organisms or deeper time courses it could silently mislead.

**Rule violation.** Not strictly anti-hallucination — the number is in the KG — but a presentation-layer drift that propagates wrong mental models. Related to Rule 7 in spirit: numbers must be presented in a form the researcher can interpret without needing to know the aggregation semantics.

**Proposed skill change (candidate for v3).** Add to `python-api-guide.md` or a new presentation-layer guidance section:
- **"When summarizing experiment coverage, always use per-timepoint gene counts (`tp_gene_count`), never the top-level `gene_count`.** The top-level field sums across timepoints and scales with #TPs. Report as single value if constant, or min/median/max range if it varies. For cross-class summary tables, aggregate `median(tp_gene_count)` per experiment, then min/median/max across experiments."
- Ideally the API itself should rename the top-level field to `total_row_count` or add a sibling `genes_per_timepoint` / `n_timepoints` pair so the aggregation semantics are self-documenting. Cross-link to A-series MCP friction (`api_coverage.md`).

**How to apply.** Every show-phase table, notebook entry, manifest row, or methods description that mentions gene count must use `tp_gene_count` from `experiments_to_dataframe` (which is the per-timepoint value). If reporting a single scalar per experiment, use the median (or single value for non-time-course). The top-level `gene_count` should never appear as "N genes" in a human-facing summary without the word "cumulative" adjacent.

Also saved as user-memory feedback (`feedback_gene_count_per_timepoint.md`) so future sessions start with the right convention.

---

## Process retrospective

_(populated at Task 14)_
