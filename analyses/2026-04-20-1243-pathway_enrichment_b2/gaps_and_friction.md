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

---

## Process retrospective

_(populated at Task 14)_
