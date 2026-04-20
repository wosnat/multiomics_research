# Gaps and friction — pathway enrichment B2

Live document — appended to at each step's decide phase. Sections per artifacts-guide §gaps_and_friction.md format. Pre-seeded items will be added by Task 12 Step 5 (from `docs/superpowers/specs/2026-04-18-research-methodology-v3-improvements-from-b2.md` §4); this file records the execution-time delta.

---

## KG data bugs

### 2026-04-20 — `de_enrichment_inputs` background collapsed to foreground (FIXED same day)

**What happened.** Step 2 enrichment run produced 0 significant pathways across 43 MED4 × `table_scope` bg clusters — biologically impossible for Tolonen R1 (6 TPs × 2 dir, strong N-deprivation time course). Debugging showed `gene_ratio == bg_ratio` in every row — Fisher 2×2 collapsed because `background` == `foreground`.

**Root cause.** `multiomics_explorer.de_enrichment_inputs` iterated over `differential_expression_by_gene` rows and applied `continue` when `row_direction not in ("up", "down")` BEFORE populating `background`. Since `_STATUS_TO_DIR` doesn't map `"not_significant"` → anything, every not-significant row was skipped. Result: per-cluster `background` contained only significant rows in that (TP, direction) — always equal to the foreground gene set, never the full quantified universe.

Compounded by: cluster key was `{exp_id}|{tp}|{direction}`, so even within significant rows the up-cluster bg got only `significant_up` genes, not all sig genes at that TP.

**Detection signal.** `bg_count == count` and `bg_ratio == gene_ratio` across all rows in the enrichment_all.csv. If this repeats in the future, the reader should suspect the background-construction pipeline immediately.

**Fix.** Two-pass background construction in `de_enrichment_inputs`: pass 1 builds per-(exp, tp) quantified universes from all rows (including not_significant); pass 2 builds per-(exp, tp, direction) foregrounds from significant rows. Up/down clusters at the same (exp, tp) share the same background (the full quantified set at that TP). Researcher merged the fix in sibling repo `multiomics_explorer` on 2026-04-20. Editable install picked it up automatically.

**Empirical verification.** Re-running `03_run_enrichment.py`:
- Tolonen R1 cluster backgrounds: all 12 clusters show bg=1697 (MED4 Affy array size, 100% match to expected quantified universe).
- Up/down clusters at the same TP: identical background sets (set-equality verified).
- 3h|down cluster (zero significant genes) now exists with full bg + empty gene_set (previously missing entirely).
- Significant hits: 225 across 11,239 tests total. R × key-pathway agreement with expected direction: 40/41 significant hits concordant (97.6%). Biological sanity restored.

**Impact on B2.** None on final results — bug found and fixed before committing any enrichment-derived outputs. The Step 2 do-phase commit captures only post-fix artifacts.

**Impact on B1** (prior analysis, `analyses/2026-04-09-1713-pathway_enrichment_b1/`): B1 used a custom `enrich_utils/` package that likely built backgrounds manually, so B1 may be unaffected. Worth confirming; if B1 used `pathway_enrichment` via MCP with default `table_scope` bg, its results are invalid and need re-running. Flagged for researcher review.

**Skill-relevance.** This was a collaboration loop: I (Claude) detected the signal (`gene_ratio == bg_ratio` everywhere), diagnosed the root cause from the `de_enrichment_inputs` source, and wrote a hand-off prompt for the explorer team. Researcher merged the fix in a separate session. The workflow pattern (detect → diagnose → hand off with reproducing case → resume downstream) is worth preserving as a process template. Candidate for v3 skill addition: "When an upstream primitive returns biologically implausible results, diagnose the primitive before adjusting the analysis — don't compensate in the analysis code for an upstream bug."

---

## KG gaps

_(none observed yet)_

---

## MCP friction

_(pre-seeded items from v3 meta-doc §4.5 A1–A4 will be added in Task 12 Step 5; execution-time items below)_

### 2026-04-20 — Step 1b: ~~`pathway_enrichment` has no `term_ids` filter~~ — RETRACTED (my error)

**RETRACTED.** I claimed `pathway_enrichment` only accepts `level`, forcing all-or-nothing at a single depth. **This was wrong.** Loading the tool schema directly shows:

```
"term_ids": {"anyOf": [{"items": {"type": "string"}, "type": "array"}, {"type": "null"}],
             "description": "Specific term IDs to test. Combines with level to scope rollup."}
```

`level` description confirms: "At least one of level or term_ids required." Three valid modes exist: `level`-only, `term_ids`-only, or `level + term_ids` (scoped rollup).

**What this means for the analysis.** go_bp is fully usable via a hand-curated `term_ids` panel across L3–L5. The N-specific terms we want (`go:0071941` L3, `go:0071705` L4, `go:0006995` L5 if MED4-annotated, etc.) can all be passed as a single list. The "DAG ontology friction" I described does not exist — the API already supports the workflow.

**Final decision unchanged.** KEGG L2 remains Pick 2, because at the same single level it bundles N-metabolism + photosynthesis + ribosome + AA-metabolism cleanly — no hand-curation required. But the reasoning that drove that choice was partly based on my incorrect friction claim. If curation overhead were the blocker, we'd have an easier alternative than I thought. Treating this as a pragmatic pick, not a forced one.

### 2026-04-20 — Step 1b: Meta-pattern — claimed API limitation without checking schema (anti-hallucination / Rule 7)

**What happened.** I asserted that `pathway_enrichment` accepts only `ontology` + `level`, no `term_ids`. This was based on memory of the tool signature from earlier orientation, not a direct schema re-read. The researcher pushed back ("are you sure?"); I loaded the schema via `ToolSearch` and saw the claim was wrong.

**Root cause.** Same pattern as the Hennon → Aharonovich error earlier: presenting a claim about KG/MCP state from intrinsic/stale memory instead of verifying against the current source. Tool signatures evolve; I hadn't re-checked the `pathway_enrichment` schema in this session. The earlier orientation pull of tool schemas may have been for older versions, or I may have confused tools in my summary.

**Rule violation.** Rule 7 (anti-hallucination) + Rule 1 (KG/tools as source of truth). API surface claims are data claims, not interpretation. Must come from current schema, not memory.

**Proposed skill change (v3 candidate, generalizes the earlier `list_publications` attribution entry).**
- Add to `kg-rules.md` or `anti-hallucination.md`: **"Before asserting what any MCP tool accepts or returns, re-load the current schema via `ToolSearch select:<tool_name>` and read the parameter definitions directly. Memory of tool signatures from earlier in the session can be stale (tools evolve), and from training-era knowledge is unreliable."**
- Concretely: whenever claiming a tool limitation ("X doesn't support Y", "X can only take Z"), the claim MUST be preceded by a fresh schema load in the same session turn. Otherwise treat the claim as unverified hypothesis.

**How to apply.** When deciding ontology/tool choices based on API capability: always verify the capability by schema inspection FIRST, not after the researcher catches the error. Two cases now caught in this step — pattern is real.

### 2026-04-20 — Step 1b: `pathway_enrichment` level-only mode for DAG ontologies is a UX refinement opportunity (minor)

**What happened.** Tool accepts `level`, `term_ids`, or both. However, the default / first-reach convention in the docs (`docs://tools/pathway_enrichment`) and existing example (`docs://examples/pathway_enrichment.py`) shows `level`-only usage. For shallow ontologies (cyanorak_role, kegg pathway) that's fine; for DAG ontologies like GO, a researcher unaware of `term_ids` will reach for `level` and silently drop biologically meaningful terms that sit off-level.

**Not a bug — a discoverability concern.** The parameter exists and works. But when the research question is narrow (specific N-biology), the level-only default guides a weaker choice.

**Proposed doc / example change (minor, candidate for API-docs revision).**
- Add an example to `docs://examples/pathway_enrichment.py` showing hand-curated `term_ids` for DAG ontologies (GO), paired with `search_ontology` to find the panel.
- In `docs://tools/pathway_enrichment` Response format or Common mistakes section: note that for DAG ontologies (go_*), level-only selection may miss biologically-meaningful terms at heterogeneous depths; recommend `term_ids` for narrow research questions.
- In `ontology_landscape` summary, flag DAG ontologies (go_*) explicitly so the researcher knows to consider `term_ids` selection upstream.

**Impact on B2.** None for the enrichment itself (KEGG L2 chosen for unrelated reasons). Captured as API UX feedback for the explorer team — goes to `api_coverage.md` at Step 5.

### 2026-04-20 — Step 1b: design notes for future GO-aware pathway enrichment (captures why we did NOT use GO BP for B2)

**Context.** Captured so that future work on "redefining `pathway_enrichment` for GO" has the statistical + MED4-data-availability context that led us to pick KEGG L2 over GO BP for B2.

**The structural issue for DAG ontologies.** GO is a DAG; genes annotated to a specific child term (e.g., `go:0006995` cellular response to nitrogen starvation, L5) are also rolled up to parents (`go:0043562` L4 cellular response to N levels → `go:1901698` L3 response to N compound → ... → root). Three consequences for Fisher-ORA enrichment:

1. **Single-level enrichment (`level=N` only)** either misses specific terms at L>N (too shallow) or misses breadth at L<N (too deep). For MED4 N-biology, the biologically-direct term (`go:0006995` L5) has <5 MED4 genes → filtered out regardless of level choice.
2. **Multi-level enrichment (`term_ids=[L3∪L4∪L5]`)** introduces parent-child redundancy: same genes counted in overlapping terms. BH correction assumes independence; FDR inflates. A clean 68-test L3 Fisher becomes a non-independent 343-test multi-level Fisher with distorted p-value distribution.
3. **Hand-curated `term_ids`** avoids redundancy but introduces curation bias: "enriched relative to my chosen panel" is not "enriched relative to the GO BP background."

**Concrete MED4-data constraint.** MED4's GO BP annotation is sparse for N-limitation-specific terms across ALL levels. Only 2 N-related GO BP terms pass min-gene-set-size=5 in MED4:

| Term | Level | MED4 genes |
|---|---|---|
| `go:0071941` nitrogen cycle metabolic process | L3 | 13 |
| `go:0071705` nitrogen compound transport | L4 | 18 |

Every other N-labeled GO BP term (nitrogen fixation, starvation response, utilization regulation, compound response) has <5 MED4 genes. Going deeper does not add N-specific signal for this organism.

**Why KEGG L2 wins for B2.** KEGG is a shallow orthology-pathway mapping (not a DAG), so level-based enrichment gives a clean Fisher null. At a single L2, KEGG bundles all canonical N-response anchors:

| KEGG pathway (L2) | MED4 genes |
|---|---|
| `ko03010` Ribosome | 54 |
| `ko00195` Photosynthesis | 51 |
| `ko00260` Glycine/serine/threonine metabolism | 25 |
| `ko00240` Pyrimidine metabolism | 23 |
| `ko00250` Alanine/aspartate/glutamate metabolism | ~10–15 |
| `ko00910` Nitrogen metabolism | ~10–12 |
| `ko00710` Carbon fixation in photosynth org. | ~5–8 |

One-level, non-DAG, clean interpretation of p_adjust.

**Considered alternative: GO BP `level=4` (132 terms)**, which has `nitrogen compound transport` (18 genes, the biggest usable N-GO term) at its native level. 50% genome coverage, median 13 genes/term — sweet spot. Would be a valid single-level GO BP choice. Rejected because (a) KEGG L2 still bundles all 4 canonical anchors, GO BP L4 bundles only transport and lacks clean photosynthesis/ribosome anchors at that level; (b) cross-ontology agreement is more informative between orthogonally-framed ontologies (functional role vs pathway map) than between two hierarchical ontologies.

**Candidate feature for a future "GO-aware" revision of `pathway_enrichment`.**

A DAG-aware enrichment mode would be a genuinely novel API surface (vs the existing `term_ids` workaround). Design principles to carry into that feature work:

1. **Parent-child deduplication.** For a given DE foreground, roll up to the most-specific term that's significant; don't score the parent separately if the signal is fully explained by the child. Prior art: topGO's "elim" and "weight" algorithms.
2. **Respect annotation sparsity.** For organisms where the DAG has thin annotation at deep levels (like MED4), dynamically choose the most-specific level with ≥min_gene_set_size coverage per subtree, rather than a fixed numeric level.
3. **Keep Fisher null interpretable.** Any multi-level mode must clearly document the "pathway background" definition (full genome? all annotated terms? subtree of a root term?). The current `level=N` has a clean answer ("all level-N terms that pass min_gene_set_size"); multi-level does not until the null is defined.
4. **Ontology-type tag in `ontology_landscape`.** Flag ontologies as flat (cyanorak_role, tigr_role, cog_category), hierarchical-tree (KEGG), or DAG (GO), so researchers choose enrichment mode appropriately. Current landscape just returns per-(ontology, level) stats — researcher has to know DAG-vs-flat from outside the tool.
5. **Cross-reference with `genes_by_ontology` DAG handling.** Parent-child rollup semantics should be consistent across `pathway_enrichment` and `genes_by_ontology` so drill-down (`result.explain`) works sensibly on multi-level results.

These are feature-scope thoughts, not urgent B2 blockers. Captured here so the design lineage is traceable when the tool gets revised.

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
