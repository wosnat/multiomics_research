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

### 2026-04-20 — Step 2 decide: signed_score cap needs explicit methodology-layer reasoning (candidate statistical-rigor.md addition)

**What happened.** Spec §5 Step 4 M2 set `SCORE_CAP = 10` for Layer A scoring with the rationale "corresponds to padj = 1e-10; beyond this, magnitude differences no longer carry useful information." That's correct reasoning but compressed into a one-liner; the researcher pushed back during Step 2 decide with "is 10 the right cutoff?" which surfaced the deeper question: *at what point does signed_score stop measuring biology and start measuring measurement precision?*

Empirical analysis of this dataset (225 significant rows, `explore_step2_score_distribution.py`):

| \|signed_score\| regime | padj range | interpretation |
|---|---|---|
| 0–3 | > 1e-3 | below sig threshold (excluded by padj<0.05) |
| 3–10 | 1e-3 to 1e-10 | **biologically meaningful** — cluster differences reflect real coverage/fold |
| 10–15 | 1e-10 to 1e-15 | semi-saturated — differences mix biology and detection precision |
| 15–25 | 1e-15 to 1e-25 | saturation — padj depends on integer overlap count; ±1 gene shifts padj by many orders of magnitude |
| 25+ | < 1e-25 | statistical ceiling — differences are detection-threshold artifacts |

~15% of this dataset's significant hits exceed \|s\|=10; most are ribosome-DOWN at R/CTX clusters at \|s\|=23–33.

**Root cause (methodology-layer gap).** The spec presented `SCORE_CAP = 10` as a single value with one purpose. In practice, capping has two purposes:
- **Scoring cap** (for mean-based Layer A): prevents single-pathway domination. Principled choice ≈ "end of biologically-meaningful regime" ≈ \|s\|=10.
- **Visualization cap** (for heatmap colormap): spreads color dynamic range over the biologically-interesting band. For this dataset, ±5 keeps the 0–5 region visually resolved; anything above is flagged with saturation stars.

Collapsing both into one cap loses information on the viz side. Using the same number for both is a coincidence of this dataset, not a principle.

**Proposed v3 skill change (candidate for statistical-rigor.md).**

> **When scoring with `-log10(padj)`-based metrics, distinguish scoring-layer caps from visualization-layer caps explicitly.**
>
> - **Scoring cap:** set at the boundary between the "biological dynamic range" and the "statistical saturation regime." For Fisher ORA on moderate-sized pathway sets (10s to 100s of genes), this is typically \|s\| = 10 (padj = 1e-10). Check empirically: plot `|signed_score|` quantiles on significant rows; the knee of the distribution (~p90–p95) is often where saturation kicks in.
> - **Visualization cap:** set lower than the scoring cap to preserve color resolution in the biologically-meaningful band. Combine with saturation markers (`*`/`**`/`***` in cells, or alternate visual annotation) to preserve the "this cell is beyond cap" info that readers need for interpretation.
> - **Document both caps in methods.md** with the empirical distribution that motivated them, not as inherited spec defaults.

**Impact on B2.** Scoring cap stays at spec's ±10 (principled). Visualization cap for Step 2 QC and Step 5 Fig 1 = ±5 + saturation stars. Documented in Step 2 decide notebook entry.

### 2026-04-20 — Step 2 decide: within-ontology pathway redundancy is not covered by LOO (candidate statistical-rigor.md + spec-revision)

**What happened.** Step 2 explore drill-down on J.1 ATP synthase / ko00190 Ox phos / ko00195 Photosynthesis revealed that all three share the same 9 MED4 atp genes (PMM1438–PMM1457, atpA-I operon) as their enriched subset. If all three land in the Step 3 signature, the same biological event (atp operon coordinately downregulated) gets counted three times in Layer A scoring.

Checking the spec: LOO stability (§5 Step 4 M3) detects **fragility** (single-pathway or single-experiment dominance). It does NOT detect **redundancy** (correlated gene overlap inflating confidence):
- Remove J.1 alone → ko00190 + ko00195 still contribute atp signal. Score barely moves.
- LOO declares the score "robust," but the robustness is illusory — three pathways voting for the same 9 genes.

LOO and redundancy are orthogonal failure modes.

**Root cause.** Fisher ORA tests each (cluster × term) independently. When pathway sets overlap in gene membership, their enrichment p-values are not statistically independent — but the `pathway_enrichment` output doesn't surface this, and the signature derivation rule (`n ≥ 3 same-direction clusters`) doesn't audit it either. Within-ontology overlap is more acute than cross-ontology (per spec, `score_A` is computed per-ontology, so cross-ontology overlap is limited to inflating the cross-ontology agreement stability check).

**Per-ontology concrete concerns from this dataset:**
- **Within-KEGG:** `ko00190 Oxphos (33 genes) ∩ ko00195 Photosynthesis (51 genes) = 9 atp genes.` If both enter signature, atp signal is 2× counted.
- **Within-cyanorak:** J.1/J.2/J.7/J.8/K.2 are hierarchically disjoint L1 categories — likely gene-disjoint, but unverified. Audit during Step 3 explore.

**Proposed v3 skill change (candidate for statistical-rigor.md, new stability-check spec section).**

> **Stability check M4 — Pathway-gene-overlap audit (within-ontology).**
>
> For each ontology's signature pathways, compute pairwise `Jaccard(pathway_i.gene_members, pathway_j.gene_members)` and detect strict subset relations. Flag pairs with Jaccard > 0.5 or subset relations as redundancy candidates. Report as heatmap + flag list.
>
> Resolution options (researcher-chosen per analysis):
> - **(A) Audit-only** — report the overlaps as an M4 stability observation, keep all pathways in the signature, document redundancy in caveats.md.
> - **(B) Post-filter in Step 3** — within each overlap cluster (terms sharing >threshold genes), keep only the most-enriched term in the signature.
> - **(C) De-weight during scoring** — scale each pathway's contribution to Layer A by `1 - max_jaccard_with_higher_enriched_signature_term`, so redundant contributions are dampened.
>
> Default recommendation: (A) unless redundancy is extreme (Jaccard > 0.8 or strict subset).
>
> **Why within-ontology:** per spec §5 Step 4, Layer A score is per-ontology; cross-ontology overlap inflates the cross-ontology-agreement stability check but not the score itself. Within-ontology overlap inflates the score directly.

**Impact on B2.** Added new Step 3 explore sub-task before Task 8 — compute pairwise Jaccard on signature pathways per ontology; decide A/B/C based on what the audit surfaces.

### 2026-04-20 — Step 2 decide: a priori anchor lists stay locked through QC; discovered-strong pathways enter via derivation (candidate research-notebook.md / step-protocol.md addition)

**What happened.** During Step 2 explore, discovered-strong pathways emerged in the R clusters (J.1 ATP synthase `n=6/12 down`; ko00190 Ox phos `n=6/12 down`; ko00710 Calvin cycle `n=5/12 down`; D.1 Adaptation `n=3/12 up`). My first instinct was to amend `key_pathways.csv` to include them, since they're clearly canonical N-limit signatures ("ATP synthase should have been on the list!"). The researcher correctly pushed back: adding pathways to the a priori list after observing them in the data is fitting the anchor to the observation — confirmation bias on the very validation step that's supposed to be independent of the outcome.

**Root cause (methodology-layer insight).** The spec distinguishes two lists with different epistemic roles:
- **key_pathways.csv** (locked at Step 1b): *a priori* prediction based on textbook biology. Used as Step 2 QC sanity check — "did the signal behave at known N-response anchors?" Failure here falsifies the ontology/signal choice.
- **reference_signature.csv** (derived at Step 3): *a posteriori* discovery from R cluster enrichment data. Used to score T clusters.

Mixing these is a category error. Adding J.1 to the a priori list after seeing it hit means "the signal is validated because the pathways I added after seeing them work also work." Tautology. Backporting would also weaken the 3 zero-hit AA-biosynthesis anchors' falsification value (decision #3) — those are the anchors that DIDN'T behave as predicted, and their presence in the key panel is informative *because* they falsify the prediction.

**Proposed v3 skill change (candidate for research-notebook.md or step-protocol.md Redo path).**

> **Once an a priori list (key-pathway panel, expected marker genes, hypothesized enrichment directions) is committed to disk in the pre-data step, it is locked through the analysis. Discoveries from the data enter via the designated a posteriori mechanism (signature derivation, discovered-strong panel), not by amending the a priori list.**
>
> Exceptions:
> - If the a priori list contains a **data error** (wrong term_id, wrong organism mapping, typo), correct it in a redo-path commit with an explicit "correction, not amendment" annotation in the notebook entry.
> - If the a priori list contains a **falsified prediction** (observed no signal in either direction, contrary to the expected direction), the prior stays but gets a caveats.md entry documenting the falsification. Option (c) in B2 decision #3.
>
> The derived signature, heatmap row selection, or any other a posteriori artifact is the right place for discovered-strong pathways to appear. The key-pathway panel is not.

**Impact on B2.** J.1, ko00190, ko00710, D.1, ko01200 stay out of key_pathways.csv. They enter Step 3 signature via the `n ≥ 3` core rule automatically. The v2 heatmap (decision #2) displays them in a second block ("discovered") alongside the a priori key rows, with bold font distinguishing the two epistemic sources — that's the right place for the union view. Key panel itself unchanged.

### 2026-04-20 — Step 2 decide: NC-class calibration assumes NC is noise, but NC clusters can carry real biology that contaminates the noise floor (candidate statistical-rigor.md + spec-revision)

**What happened.** Step 2 explore found that 4 of the 225 significant enrichment rows come from NC clusters — biologically interpretable signals, not noise:
- Steglich 2006 high-light 45min (NC) → photosynthesis-DOWN at padj 3e-3 (cyanorak J.7 PSI), 3e-2 (kegg ko00195).
- Weissberg 2025 coculture day 11 (NC, no explicit N-starvation) → N-metabolism-UP at padj 8.6e-5 (cyanorak E.4), 1.6e-5 (kegg ko00910). Real coculture-induced N-scavenging biology.

These "negative control" clusters carry real signal on N-limit signature markers. Per spec §5 Step 4 M3, the NC noise floor is computed as `nc_mean_{o,b} = mean(nc_scores_{o,b})` and `nc_std_{o,b}` — i.e., whatever NC clusters happen to be in each `(ontology, background_used)` group. When those clusters carry real biology overlapping with signature anchors, the mean and SD inflate; T-score classification thresholds (`nc_mean + 2σ`) rise accordingly; real-but-weak T signature responses get misclassified as "no signal."

**Root cause (methodology-layer gap).** The spec's NC definition is "MED4 experiments unrelated to N (light, phage, salt, dark, glucose)" — a class label applied at experiment-selection time (Step 1a). But the KG contains heterogeneous NC candidates, and "unrelated to N at the experimental-design level" does not guarantee "unrelated to N at the transcriptional/proteomic signature level." Several canonical stress responses (HL-stress shutting down photosynthesis; coculture-enabled N-scavenging) overlap with the signature axes the analysis is calibrating against.

The spec treats NC as a statistical fact ("NC is the noise floor") rather than a testable hypothesis ("NC *should be* the noise floor, verified by no significant key-pathway enrichment"). Step 2 QC surfaces the violation, but the spec has no built-in mechanism to adjust.

**Proposed v3 skill change (candidate for statistical-rigor.md, spec-language refinement).**

> **NC calibration sanity check + exclusion rule.** Before computing `nc_mean` and `nc_std` per `(ontology, background_used)` group:
>
> 1. Identify NC clusters with `padj < 1e-3` significant enrichment on any a priori key-pathway anchor for that ontology.
> 2. Exclude those clusters from NC calibration (log the exclusion with cluster ID, pathway, and padj).
> 3. If exclusion reduces a calibration group below 3 clusters, flag the group as "weakly calibrated" in caveats.md; threshold-based classifications should be read as narrative indicators rather than hypothesis-test conclusions.
>
> The padj<1e-3 threshold (3 orders of magnitude below standard significance) is deliberately stricter than 0.05. At that threshold, an NC cluster's enrichment is implausibly attributable to sampling noise — retaining it in the calibration biases the threshold. Borderline hits (padj 1e-3 to 0.05) stay in NC on the principle that true noise floors include noise-adjacent fluctuations. The threshold is a design choice, not a theorem — can be tuned per dataset, but should be committed at the step-boundary, not tuned post-hoc after seeing T scores.

**Complementary v3 suggestion (for recipe skill or step-protocol).**

> **NC is a hypothesis, not a fact.** When selecting NC candidates in Step 1a, list each proposed NC experiment with "expected signature behavior" — i.e., what biology might leak across the class boundary and inflate the calibration? (e.g., "Steglich high-light may show PS-down contamination of the PS-anchor calibration"; "Weissberg coculture-replete-N may show coculture-induced N-scavenging contamination of the N-metab anchor"). Treat NC class membership as provisional until Step 2 QC confirms or refutes via the key-pathway enrichment check.

**Impact on B2.** Decision #4 applied option (b): exclude Weissberg coculture d11 up cluster from `(cyanorak_role, table_scope)` and `(kegg, table_scope)` NC calibration groups (padj<1e-3 criterion triggered). Steglich high-light 45min down cluster retained (hits at padj 3e-3 and 3e-2, above the threshold). `(*, organism)` NC groups have only 2 clusters each — flagged as weakly calibrated in caveats.md C5. Implementation requirement surfaces as a Task 10 followup (see plan).

### 2026-04-20 — Step 2 decide: heatmap-visualization conventions for pathway-enrichment QC (candidate artifacts.md addition)

**What happened.** The original Step 2 QC heatmap ([step2_key_pathway_heatmap.png](step2_key_pathway_heatmap.png), committed earlier) was unreadable at full zoom — 11 rows × 70 columns crammed into a single axis with tiny labels. Iterated through v2 design in Step 2 decide:

1. **Row selection = a priori key panel ∪ a posteriori discovered-strong (`n_sig ≥ 3` in signature-eligible R clusters).** Bold font for key, regular for discovered — separates epistemic sources visually.
2. **Columns grouped by class** (R | PC | CTX | NC) with vertical dividers and class labels above. Within each class, ordered by experiment + timepoint.
3. **T panel separate (stacked below non-T).** Different epistemic role — T is scored *against* the signature derived from R/PC. Keeping T visually isolated makes the "predict from R, validate on T" structure readable.
4. **T panel split by biological contrast** (axenic vs coculture for B2) rather than by omics type. The primary scientific question drives the divider; the lab-instrument type is secondary.
5. **Display cap ±5 (not ±10) with saturation stars** (`*`/`**`/`***` for \|s\| ≥ 5/10/20) — spreads color resolution over biologically-meaningful 0–5 band; stars preserve "beyond cap by how much" info.
6. **Uniform cell size across panels** via explicit axis positioning in figure coordinates — short panels don't stretch, wide panels don't compress.
7. **Column labels = `firstauthor_6chars | tp_short`** (non-T) or `| omics | tp` (T). Author truncation and timepoint abbreviation (`day 14 → 14d`, `steady state → ss`) keeps 70 columns readable.

Collective artifact: [step2_heatmap_cyanorak_role.png](step2_heatmap_cyanorak_role.png), [step2_heatmap_kegg.png](step2_heatmap_kegg.png), via [explore_step2_heatmap_v2.py](../scripts/explore_step2_heatmap_v2.py).

**Proposed v3 skill change (candidate for artifacts.md, heatmap-conventions section).**

A recipe skill (`recipes/pathway-enrichment-qc-heatmap`?) or an artifacts.md subsection could capture these conventions as a reusable pattern. Low priority — practical enough to apply in future analyses without a formal skill, but consolidating conventions would reduce design-iteration overhead.

**Impact on B2.** Step 2 QC artifact is now publication-adjacent quality. Step 5 Fig 1 will reuse the same design with minor refinements (extra pathways from the final signature, maybe per-panel independent normalization).

---

## Mid-way Step 2 retrospective (2026-04-21)

Captured at end-of-session after completing Step 2 decide, before Step 3 begins. Narrative retrospective that motivates the distilled lessons in the next section. Comparable to a sprint-retro: honest reflection on what worked, what friction I introduced, what I'd change next session.

### What worked

- **Every analytical claim came from a file read, not memory.** The prior session (pre-restart) had drifted into memory-based chat tables with inaccuracies. This session held the line: CSVs + script output + schema probes before any numeric claim. Zero attribution bugs this time.
- **Researcher pushbacks were productive, not dismissed.** Three moments where corrective feedback improved the analysis: (a) "don't amend key_pathways" — I'd proposed adding discovered-strong pathways (J.1 ATP synthase, ko00190 Ox phos) and the researcher caught the confirmation bias; (b) "is 5→10 vs 10→20 meaningful?" — I'd been pushing ±20 display cap, the question forced genuine re-examination and I updated to keep spec ±10 for scoring, use ±5 for visualization; (c) "axenic vs coculture, not omics" — I'd defaulted to Prot/RNA as T divider and the researcher corrected to the real biological contrast.
- **Explore scripts stayed committed with outputs, linked from the notebook.** 6 scripts + 5 artifacts (`step2_cluster_totals.csv`, `step2_r_top_pathways.csv`, 2 heatmap PNGs, `step2_score_distribution.png`) all cross-referenced. No orphaned work.
- **Skill-friction capture was substantive.** 5 entries with concrete v3 skill additions, each tied to a specific incident rather than generic advice (SCORE_CAP scoring-vs-viz, within-ontology redundancy, a priori locked, heatmap conventions, NC calibration).

### What I did badly

- **Too much heatmap iteration.** Rendered the heatmap ~5 times before it landed. Could have been 2 if I'd front-loaded the design questions (author format, timepoint shortening, direction encoding, axenic-vs-omics divider). I started coding on partial spec and paid the cost in iterations. User-visible cost was modest, but researcher time and cache burn accumulated.
- **Missed the `omics_type=NaN` bug before rendering.** The T panel first render showed "nan" as a divider label — a dataset-assumption violation catchable by `df.groupby("omics_type").size()` before building. This is the same failure class as the prior session's issues: trusting merged columns were populated without verifying. One diagnostic line would have caught it upfront.
- **Let documentation lag.** `DATA_MANIFEST.md` was 5 artifacts stale at session close — a Gate 2 violation. `decisions.md`, `hypotheses.md`, `api_coverage.md` didn't exist until the researcher asked "anything else worth capturing?" The plan's step-boundary manifest updates invite exactly this drift. Fixed in follow-up commit `160f10b`, but cleaner discipline would have been: any new artifact in a commit gets its manifest row in the same commit.
- **Over-elaborated first proposals.** When asked for a heatmap, I jumped to full design structure before confirming shape. When discussing the score cap, I pushed a specific number (±20) before establishing the biological-vs-saturation regime framework. Better pattern: confirm framing, then propose concrete values.
- **Notebook Step 2 decide entry grew too long.** ~350 lines covering 7 Q&A blocks, 4 decisions, legend, caveats. Exhaustive but not scannable for a fresh reader who needs to orient quickly. A 10-line summary at the top would help — added in follow-up commit per researcher request.

### Risk flag for next session

The `>=` → `>` spec amendment for the temporal filter (decision #1) is captured in memory + plan + notebook + decisions.md, but it's a one-character change that's easy to overlook in the scaffolded code block. Next session's first Task 7 action must be a deliberate grep check on `scripts/signature.py` and `scripts/04_derive_signature.py` to confirm the filter uses `>`, not `>=`. A 30-second verification prevents a whole-step redo. Captured as lesson L5 below.

### Distilled lessons

Below in "Session-level process lessons" — L1–L5, each with rule, concrete application, and candidate v3 skill addition. These are the actionable form of the narrative above; the retrospective itself is the evidence base.

---

## Session-level process lessons

Captured at end-of-session 2026-04-21 after completing Step 2 decide. These are process patterns that surfaced during the session — distinct from the whole-analysis Process retrospective (Task 14). Each is a candidate v3 skill addition in its own right or reinforcement of existing rules.

### L1 — Verify dataset assumptions before visualization

**What happened.** First T-panel render of the v2 heatmap showed "nan" as a divider label. Cause: `omics_type` column in `enrichment_all.csv` is populated for non-Weissberg experiments but NaN for Weissberg T experiments; I had merged on classified for `first_author` but not for `omics_type`, trusting the column was uniformly populated. One `df.groupby("omics_type").size()` check before rendering would have caught it upfront.

**Rule.** Before building a visualization or analytical transform that depends on a merged/joined column, run one diagnostic line to verify the column is populated for the relevant subset. Cheap check, high-value catch.

**Concretely.** Before any plot, aggregation, or groupby on a derived dataframe, run `df[key].isna().sum()` or `df.groupby(key).size()` and eyeball it. Same class of check as the "check API schema with 2-line call" rule.

**Candidate v3 addition** (to `research-methodology` skill, likely `kg-rules.md` or a new `data-hygiene.md`):

> **Pre-render dataset-assumption check.** Before visualization, aggregation, or any transform that assumes a column is populated: run a diagnostic line (`isna().sum()`, `groupby(col).size()`, `head()`) and inspect the result. Merged / joined / derived columns are the high-risk category — do not trust that the column is uniformly populated across the classes, organisms, or omics types the downstream code assumes.

### L2 — Confirm framing before values in interactive design

**What happened.** When the researcher asked for a better heatmap, I proposed a complete design structure (row selection, column grouping, cap, stars, labels) in one go, then iterated 4–5 times through revisions (authors format, tp shortening, dividers, omics vs condition split, legend, newlines). Some revisions were real refinement, but several could have been avoided by confirming the rough shape first. Similarly for the SCORE_CAP discussion, I proposed ±20 before establishing the evidence-saturation framework; the researcher's "is 5→10 vs 10→20 meaningful?" forced a re-framing and value change.

**Rule.** In interactive design tasks, separate framing decisions ("what axes of choice exist? what's the goal?") from value decisions ("what specific number / format / layout?"). Get framing right first, then propose values; iterate on values only.

**Concretely.** When asked "build X," structure the response as:
1. Identify the framing questions ("by-ontology vs combined?" "biological-contrast divider vs lab-instrument divider?" "display-layer cap vs scoring-layer cap?")
2. Confirm framing with researcher (short exchange).
3. Propose specific values once framing is locked.

This is a version of the brainstorming skill applied to sub-decisions within a larger task.

**Candidate v3 addition** (to `research-notebook.md` or `step-protocol.md`):

> **Separate framing from values in interactive design.** For any design task where multiple axes of choice exist (what to plot, what to cap, what to group by), identify the axes first and confirm framing with the researcher before proposing specific values. Iterating on values under a wrong frame wastes cycles; iterating under a confirmed frame is cheap.

### L3 — Manifest currency is per-commit, not per-step

**What happened.** `DATA_MANIFEST.md` was 5 artifacts stale when closing the session. Gate 2 ("manifest currency") exists in the spec but was treated as a step-boundary check — i.e., update at each Step's commit 1 (do-phase). That's what the plan says. But `cc61de1` (Step 2 decide) landed 5 new artifacts without manifest rows because they were produced during the explore phase (not the do-phase), and the plan's Task 6 Step 6 commit list didn't include manifest update. Documentation drift by design.

**Rule.** Any commit that adds or modifies a data artifact in `data/`, `results/`, or `exploration/qc/` must update the relevant manifest file in the same commit.

**Concretely.** Before every commit that includes a new file in `data/`, `results/`, or `exploration/qc/`: check whether the manifest covers it. If not, add the manifest row in the same commit. No "update later" commits. Treat manifest updates as part of the commit content, not a separate task.

**Candidate v3 addition** (to `step-protocol.md` or `artifacts.md`):

> **Manifest currency is a per-commit hard gate.** Any commit that touches files in `data/`, `results/`, or `exploration/qc/` must update the corresponding manifest file (`data/DATA_MANIFEST.md` or `results/RESULTS_MANIFEST.md`) in the same commit. This is stricter than step-boundary checking: explore-phase and decide-phase commits often land artifacts, and each needs manifest coverage at commit time. The plan's per-step commit script should always include the manifest path in `git add`.

### L4 — Notebook entry summary-at-top

**What happened.** The Step 2 decide notebook entry grew to ~350 lines (7 Q&A blocks, 4 decision writeups, heatmap legend, caveats text). Comprehensive but not scannable — a fresh reader (e.g., next-session agent, or the researcher returning after a week) has to parse the whole block to orient.

**Rule.** Long notebook entries should open with a 10–15-line "what this step resolved" summary. Details below. The summary stands alone as "if you read nothing else."

**Concretely.** For any notebook entry covering a decide-phase with ≥3 decisions or ≥5 Q&A blocks, lead with:
- One-sentence entry scope ("Step N decide: X resolved")
- Bulleted list of decisions made (with numeric labels referencing `decisions.md`)
- One-line preview of key observations / hypotheses
- Pointer to "Full detail below" for the deeper content

**Candidate v3 addition** (to `research-notebook.md`):

> **Summary-at-top for long decide-phase notebook entries.** When a notebook decide-phase entry exceeds ~50 lines (roughly: ≥3 decisions or ≥5 Q&A blocks), open with a concise summary block: (1) one-sentence entry scope; (2) bulleted decisions with `decisions.md` cross-references; (3) one-line preview of key observations; (4) pointer to full detail below. Long-entry discipline — the researcher rereads the summary; the agent rereads details when needed.

### L5 — Verify spec amendments against scaffolded code before executing

**What happened.** Decision #1 amended the spec from `hours < 3` to `hours > 3`. The plan's Task 7 scaffold has the cutoff constant + the filter mask embedded in a code block. A fresh-session agent executing Task 7 by copy-pasting the scaffold would miss the `>=` → `>` change unless they remember it from memory + notebook + plan.

**Rule.** When a decision amends something that's already scaffolded in plan/code, the scaffold must be updated (or explicitly flagged for per-execution verification) — otherwise the amendment silently fails to apply.

**Concretely.** At decide-phase commits that amend scaffolded code, do one of:
- Update the scaffold in place to match the new decision (preferred).
- If the scaffold can't be updated (e.g., it's a template that serves multiple call sites), add a **deliberate verification step** at the execution task's Step 1: "grep/diff the live code against the latest decision before running." Document the specific check ("grep for `>=` and `>` in `signature.py`'s temporal filter; confirm `>` is used").

**Candidate v3 addition** (to `step-protocol.md` Redo path section):

> **Spec-amendment ↔ scaffold reconciliation.** When a decide-phase commit amends the spec or plan semantics, either (a) update the scaffolded code in the same commit, or (b) add a deliberate verification step to the next execution task's Step 1 ("before running, grep/diff for the amended semantic"). Do not rely on future-session agents to remember multi-file amendments from notebook + memory + plan alone — a one-character semantic change is too easy to overlook.

---

## Process retrospective

_(populated at Task 14 — whole-analysis retrospective. Session-level lessons above are the per-session layer.)_
