# Pathway enrichment B2 — research notebook

Analysis directory: `analyses/2026-04-20-1243-pathway_enrichment_b2/`
Spec: `superpowers/spec.md`
Plan: `superpowers/plan.md`

---

## 2026-04-20 12:43 — Spec walkthrough

### Context
Redo of B1 with new MCP/Python enrichment API. Goals: (1) N-limitation score
per Weissberg 2025 MED4 condition, (2) evaluate research-methodology skill,
(3) evaluate new pathway_enrichment / ontology_landscape tools.

### Spec sections reviewed
- §1 Research question, §2 Goals, §3 Scope, §4 Classification (T/R/PC/NC/CTX)
- §5 Phase structure (6 steps with do/show/explore/decide per spec §5.0)
- §7 Methodology details (M1 signature, M2 scoring, M3 calibration)
- §8 Risks and contingencies
- §9 Artifact plan

### Open items to resolve in Step 1a
- Specific experiments to place in each class (T/R/PC/NC/CTX)
- Whether restricted-`table_scope` experiments will be accepted
- Non-MED4 CTX selection (3–5 experiments, 1–3 organisms)

### Decision
Proceed to Step 1a.

---

## 2026-04-20 13:41 — Step 1a: experiment discovery + classification

### Command
```bash
uv run scripts/01_select_experiments.py
```

### Outputs
- [data/experiments_classified.csv](../data/experiments_classified.csv) — 40 (experiment × timepoint) rows; 16 unique experiments
- [logs/step1a.log](../logs/step1a.log)

### QC

**Class × organism × omics (unique experiments)** [KG]

| class | organism | omics | # |
|---|---|---|---|
| T | MED4 | RNASEQ | 2 |
| T | MED4 | PROTEOMICS | 2 |
| R | MED4 | MICROARRAY | 1 |
| R | MED4 | RNASEQ | 1 |
| PC | MED4 | MICROARRAY | 2 |
| NC | MED4 | RNASEQ | 2 |
| NC | MED4 | PROTEOMICS | 2 |
| NC | MED4 | MICROARRAY | 1 |
| CTX | HOT1A3 (Alteromonas) | RNASEQ | 1 |
| CTX | MIT9313 (Prochlorococcus) | MICROARRAY | 1 |
| CTX | SS120 (Prochlorococcus) | PROTEOMICS | 1 |

**`table_scope` → background split** [KG]

| class | all_detected → `table_scope` bg | restricted → `organism` bg |
|---|---|---|
| T | 4 | 0 |
| R | 1 (Tolonen 2006) | 1 (Read 2017, filtered_subset) |
| PC | 2 | 0 |
| NC | 4 | 1 (Steglich 2006, filtered_subset) |
| CTX | 2 | 1 (Domínguez-Martín 2017 SS120 azaserine, significant_only) |

**Per-timepoint gene counts per class (median genes/TP across experiments in class)** [KG]

| class | n_exps | min | median | max |
|---|---|---|---|---|
| T | 4 | 1424 | 1637 | 1849 |
| R | 2 | 841 | 1269 | 1697 |
| PC | 2 | 1697 | 1697 | 1697 |
| NC | 5 | 198 | 327 | 1849 |
| CTX | 3 | 389 | 2241 | 3950 |

**Attribution (from `list_publications` join)** [KG]

- T1–T4: Osnat Weissberg 2025 (`10.1101/2025.11.24.690089`)
- R1, PC1, PC2, CTX2: Andrew C Tolonen 2006 (`10.1038/msb4100087`)
- R2: Robert W Read 2017 (`10.1038/ismej.2017.88`)
- NC1: Dikla Aharonovich 2016 (`10.1038/ismej.2016.70`)
- NC2: Osnat Weissberg 2025 (same paper as T, coculture-vs-axenic contrast under N-replete)
- NC3, NC4: José Ángel Moreno-Cabezuelo 2023 (`10.1128/spectrum.03275-22`)
- NC5: Claudia Steglich 2006 (`10.1128/JB.01097-06`)
- CTX1: Osnat Weissberg 2025 (HOT1A3 axenic N-starvation)
- CTX3: Maria Agustina Domínguez-Martín 2017 (`10.1128/mSystems.00008-17`)

**R cluster-count arithmetic for Step 3** [KG + interpretation]

- Step 2 enrichment: Tolonen 6 TPs × 2 dir = 12 clusters; Read 3 TPs × 2 dir = 6 clusters. Total **18 R enrichment clusters.**
- Step 3 signature derivation (after `timepoint_hours < 3` filter): Tolonen 0h dropped, 5 late × 2 = 10 clusters (table_scope bg); Read all 3 retained × 2 = 6 clusters (organism bg). Total **16 signature-eligible clusters** — core rule `≥3` has comfortable margin.

### Chat exploration

**Q1: Is the NC1 attribution "Hennon"?**
Data: `list_publications(organism="MED4", verbose=true)` returned `doi=10.1038/ismej.2016.70, authors=["Dikla Aharonovich","Daniel Sher"], year=2016`.
Finding: [KG] NC1 is **Aharonovich 2016**, not Hennon. Paper title: "Transcriptional response of Prochlorococcus to co-culture with a marine Alteromonas: differences between strains and the involvement of putative infochemicals."
Impact: classification unchanged; attribution corrected in rationale and CSV. Logged as anti-hallucination skill-friction item in [gaps_and_friction.md](../gaps_and_friction.md).

**Q2: Is there a Weissberg MED4 coculture-vs-axenic contrast at the proteomics level (matched NC for T3/T4)?**
Data: `list_experiments(publication_doi=["10.1101/2025.11.24.690089"])` returned 10 experiments; the 2 `coculture` treatment_type experiments are both RNASEQ (MED4 and HOT1A3 sides).
Finding: [KG] no matching **proteomics** coculture-vs-axenic contrast in Weissberg. Proteomics NC coverage comes from Moreno-Cabezuelo 2023 glucose experiments (NC3, NC4).
Impact: accepted — proteomics NC is less matched-design than RNASEQ NC2, but Moreno-Cabezuelo proteomics provides reasonable noise floor for `table_scope` bg proteomics clusters.

**Q3: Tolonen table_scope verified?**
Data: `experiments_classified.csv` rows for `msb4100087...med4_med4_microarray` show `table_scope=all_detected_genes`, `tp_gene_count=1697` constant across all 6 TPs (0/3/6/12/24/48h). MED4 ORFome ≈ 1716 (Affymetrix array covers 1697 of them).
Finding: [KG] Tolonen R1 correctly tagged `all_detected_genes`; per-TP coverage is 1697 = Affymetrix array detected set.
Impact: R1 routes to Call A (`background="table_scope"`) in Step 2. No change.

**Q4: Should enrichment be run for all timepoints, or filter early TPs upstream?**
Data: spec §5 Step 3 — "Early clusters are *not dropped from the analysis* — they remain in `enrichment_all.csv`, are scored in Step 4, and appear in Fig 2." The temporal filter (`<3h`) applies ONLY at Step 3 signature derivation.
Finding: [interpretation, following spec] all timepoints run through `pathway_enrichment` in Step 2. 0h, 3h cluster divergence from late cluster behavior is a *feature* of Fig 2 trajectory, not noise to filter.
Impact: Step 2 script runs on all 16 experiments' full timepoint set. No upstream TP filter.

**Q5: Is per-timepoint gene_count the right presentation metric?**
Data: top-level `gene_count` = sum across TPs (Tolonen 6 TPs × 1697 = 10,182); per-TP `tp_gene_count` = detection coverage per TP.
Finding: [interpretation] top-level `gene_count` is misleading because it scales with #TPs, not detection power. Per-TP is what the researcher needs to gauge Fisher ORA background.
Impact: QC table above uses median `tp_gene_count` per experiment; top-level `gene_count` banned from human-facing summaries. Logged in [gaps_and_friction.md](../gaps_and_friction.md) Skill/methodology friction; saved as feedback memory.

### Exploration (agent-driven QC)

- Restricted-`table_scope` inclusions accepted with caveats (to be documented in `caveats.md` at Step 5):
  - R2 Read 2017 `filtered_subset` (top 50% by expression): signed_score remains comparable across bg groups; the "absent" genes bias is documented.
  - NC5 Steglich high-light `filtered_subset`: provides organism-bg NC calibration for R2 and CTX3 interpretation.
  - CTX3 Domínguez-Martín SS120 azaserine `significant_only`: expected inflated fold_enrichment (spec §8 risk 3); used for Fig 1 context only, not for MED4 signature anchoring.
- **T3 axenic proteomics enters `death` phase at day 31+89** (per `growth_phases`) — biology feature, not classification issue. Expected behavior: late axenic proteomics may diverge from mid-course transcriptional pattern. Relevant to [interpretation] in Step 4.
- PC size concern: 2 PC experiments × 1 TP = 4 total clusters (all `table_scope` bg). Adequate for a PC range per `(ontology, bg)` group, but tight — any PC outliers will dominate. Researcher noted cyclic fallback: if PC calibration is too weak, R score can stand as a "forced PC" upper reference.
- No CTX4 proposed: soft cap (3–5) respected with 3 experiments across 3 organisms for Fig 1 legibility.

### Decision

Classification locked. 16 experiments: 4 T / 2 R / 2 PC / 5 NC / 3 CTX. Proceed to Step 1b (ontology landscape + selection).

Skill-friction items captured:
- [publication attribution drift](../gaps_and_friction.md#2026-04-20--step-1a-publication-attribution-drift-anti-hallucination--rule-7) — mandate `list_publications` before naming authors.
- [cumulative gene_count misreporting](../gaps_and_friction.md#2026-04-20--step-1a-gene_count-misreported-as-cumulative-instead-of-per-timepoint) — use per-TP values; top-level field banned from summaries.

---

## 2026-04-20 16:XX — Step 1a REDO: drop CTX1 per Step 1b ontology-compatibility review

### Why the redo
During Step 1b MCP orientation (`ontology_landscape` per organism), discovered that **Alteromonas macleodii HOT1A3 has no `cyanorak_role` or `tigr_role` annotations** (cyanobacteria-specific ontologies; HOT1A3 is Gammaproteobacteria). With cyanorak_role L1 selected as Pick 1 (MED4-optimal), CTX1 (Weissberg HOT1A3 axenic N-starvation RNA-seq) would appear as a blank column in the cyanorak enrichment panel — no signal, no conservation signal to read from the figure.

### Decision
Drop CTX1 from the classification. New CTX set is a sister-Prochlorococcus conservation test (MIT9313 Tolonen companion + SS120 Domínguez-Martín azaserine). Both covered by all Pick 1/Pick 2 ontologies.

### Impact
- **Classification:** 15 experiments total (was 16). 4 T / 2 R / 2 PC / 5 NC / **2 CTX** (was 3). Spec §4 CTX soft cap (3–5) relaxed to 2, justified by ontology compatibility.
- **Biology narrative:** CTX becomes a clean within-Prochlorococcus conservation test — the right scope for scoring Weissberg T against a MED4-derived signature.
- **Artifact:** `data/experiments_classified.csv` regenerated (37 rows, 15 unique experiments). Old 16-experiment version remains traceable in commit `1b1878a`.

Redo follows spec §Redo path: new commit, new notebook entry appended (not overwriting).

---

## 2026-04-20 — Step 1b: ontology landscape + selection (do-phase)

### Command
```bash
uv run scripts/02_ontology_landscape.py
```

### Outputs
- [data/landscape_Prochlorococcus_MED4.csv](../data/landscape_Prochlorococcus_MED4.csv) — 62 rows
- [data/landscape_Prochlorococcus_MIT9313.csv](../data/landscape_Prochlorococcus_MIT9313.csv) — 64 rows
- [data/landscape_Prochlorococcus_marinus_subsp._marinus_CCMP1375_SS120.csv](../data/landscape_Prochlorococcus_marinus_subsp._marinus_CCMP1375_SS120.csv) — 62 rows
- [data/nitrogen_ontology_search.csv](../data/nitrogen_ontology_search.csv) — 56 N-term rows across 7 ontologies
- [ontology_selection.md](../ontology_selection.md) — design doc with ranked picks + rejected alternatives
- [exploration/key_pathways.csv](key_pathways.csv) — 11-row key-pathway panel with expected directions

### QC — coverage per organism

**MED4 (13 exps, T+R+PC+NC)** [KG] top 5 by `relevance_rank`:

| rank | ontology | level | n_terms | genome_cov | median genes/term |
|---|---|---|---|---|---|
| 1 | tigr_role | 0 | 77 | 87% | 13 |
| **2** | **cyanorak_role** | **1** | **69** | **73%** | **14** ← Pick 1 |
| 3 | go_mf | 2 | 35 | 58% | 25 |
| 5 | go_bp | 3 | 68 | 54% | 13 |
| 32 | **kegg** | **2** | **97** | **38%** | **12** ← Pick 2 |

**MIT9313** and **SS120**: both cyanobacteria, both covered by cyanorak_role and kegg. See [ontology_selection.md](../ontology_selection.md) for full per-organism tables.

### QC — nitrogen-term density across ontologies

| ontology | N-term matches | top scoring |
|---|---|---|
| cyanorak_role | 2 | `E.4` Nitrogen metabolism (L1) |
| tigr_role | 1 | `role:160` Nitrogen metabolism |
| go_bp | 14 | `go:0009399` N fixation, `go:0019740` N utilization |
| kegg | 12 | `kegg.pathway:ko00910` Nitrogen metabolism (L2) |
| go_mf | 12 | mostly enzyme-activity |
| ec | 13 | N-related EC classes |
| pfam | 2 | Carbon-nitrogen hydrolase, P-II regulator |

### QC — GO BP N-term coverage per level in MED4 [KG]

Investigated because go_bp was a Pick 2 candidate. **All N-labeled GO BP terms across L3-L5 and their MED4 gene counts** (from `genes_by_ontology`):

| Term | Level | MED4 n_genes | Passes ≥5 filter? |
|---|---|---|---|
| `go:0071941` nitrogen cycle metabolic process | L3 | 13 | ✓ |
| `go:0071705` nitrogen compound transport | L4 | 18 | ✓ |
| `go:0019740` nitrogen utilization | L3 | 2 | ✗ |
| `go:1901698` response to nitrogen compound | L3 | <5 | ✗ |
| `go:0009399` nitrogen fixation | L4 | 0 | ✗ |
| `go:0141067` intracellular nitrogen homeostasis | L4 | 2 | ✗ |
| `go:0043562` cellular response to nitrogen levels | L4 | <5 | ✗ |
| `go:1901699` cellular response to nitrogen compound | L4 | <5 | ✗ |
| `go:0006995` **cellular response to nitrogen starvation** | L5 | **<5** ← the ideal term, MED4-sparse | ✗ |
| `go:0006808` regulation of N utilization | L5 | 3 | ✗ |
| `go:0090293` N catabolite regulation of transcription | L5 | <5 | ✗ |
| `go:1902170` cellular response to reactive N species | L5 | <5 | ✗ |

**Finding** [interpretation]: MED4's GO BP annotation is annotation-sparse at N-limitation-specific terms across all levels. Only 2 N-specific GO BP terms pass min-gene-set-size=5 regardless of level. Going deeper doesn't add N-specific signal.

### QC — KEGG pathway L2 anchor coverage in MED4 [KG]

| KEGG pathway L2 | MED4 n_genes |
|---|---|
| `kegg.pathway:ko03010` Ribosome | 54 |
| `kegg.pathway:ko00195` Photosynthesis | 51 |
| `kegg.pathway:ko00260` Glycine/serine/threonine metabolism | 25 |
| `kegg.pathway:ko00250` Alanine/aspartate/**glutamate** metabolism | 15 |
| `kegg.pathway:ko00910` **Nitrogen metabolism** | 6 |
| `kegg.pathway:ko00710` Carbon fixation | ~5–8 |
| `kegg.pathway:ko01310` Nitrogen cycle | 0 (not MED4-annotated) |

All 4 canonical key-pathway anchors (N-metab + photosynthesis + ribosome + AA-metab) present at the same L2. Clean Fisher null over 97 pathways.

### Chat exploration

**Q1: Can MCP orientation run without running the landscape script?** [interpretation]
Data: plan §Task 3 Step 1 says "Run MCP orientation queries … in chat, not the script." Spec §5 Step 1b do names the interactive MCP sub-phase.
Finding: the MCP-orientation queries are done interactively before the script, to let the researcher see available ontologies and choose. The script (`02_ontology_landscape.py`) locks the choice by reproducing the per-organism landscape as CSVs; the researcher's decisions live in `ontology_selection.md` + `key_pathways.csv`.
Impact: script structure finalized — three loops (per-organism landscape), one search_ontology loop (nitrogen across 7 ontologies), no researcher-decision logic embedded in code.

**Q2: Does HOT1A3 (Alteromonas, Gammaproteobacteria) have any cyanobacteria-ontology annotations?**
Data: `ontology_landscape(organism="Alteromonas macleodii HOT1A3", ...)` summary — cyanorak_role and tigr_role both absent from the returned `by_ontology` dict; top ranks are pfam L0 (65%), go_mf L2 (50%), go_bp L3 (45%).
Finding: [KG] HOT1A3 has no cyanorak_role or tigr_role annotation — expected (these are cyanobacteria-specific ontologies).
Impact: with cyanorak_role L1 as MED4-optimal Pick 1, CTX1 would be a blank column. Dropped CTX1 in Step 1a redo (see previous notebook entry) — sister-Prochlorococcus conservation test (MIT9313 + SS120) is the right scope.

**Q3: Is there a "relevant BRITE tree" to use as Pick 3?**
Data: MED4 landscape — brite trees range from 3–34 terms, 4–32% genome coverage. `transporters` (most N-relevant): 4.3% coverage, 3 L0 terms, mixed P+N "Phosphate and amino acid transporters" — no pure N term.
Finding: [KG + interpretation] no BRITE tree is a clear Pick 3 for N-biology. `transporters` is closest but conflates P+N and has <5% coverage. `enzymes` L1 at 32% is usable but largely redundant with KEGG pathway membership.
Impact: rejected BRITE as Pick 3. BRITE drill-down deferred to `explore_*.py` scripts post-Step 2 if needed. Captured in `ontology_selection.md` Considered alternatives.

**Q4: Can we use GO BP with multi-level `term_ids` instead of level-only?**
Data: `pathway_enrichment` tool schema (`ToolSearch select:mcp__multiomics-kg__pathway_enrichment`) confirms `term_ids: list[str] | None = None` parameter exists. Description: "Specific term IDs to test. Combines with level to scope rollup. At least one of level or term_ids required."
Finding: [KG] API supports hand-curated `term_ids` — my earlier claim that it didn't was wrong. Logged as anti-hallucination meta-pattern friction in [gaps_and_friction.md](../gaps_and_friction.md).
Impact: go_bp with cross-level `term_ids` panel is mechanically feasible, but the N-signal doesn't grow (same 2 usable N-terms regardless of level), DAG redundancy inflates BH correction (343 non-independent tests), and curated `term_ids` = key-pathway panel by another name. Rejected for B2. Design thoughts for future "DAG-aware pathway_enrichment" captured in `gaps_and_friction.md`.

**Q5: Between GO BP and KEGG for Pick 2, what's the real trade-off?**
Data: see `ontology_selection.md` Head-to-head section and the QC tables above.
Finding: [KG + interpretation] KEGG L2 has all 4 canonical anchors (N-metab + photosynthesis + ribosome + AA-metab) at a single level with clean gene sets (6–54 genes). GO BP L3 has only 1 N-specific term passing the min-5 filter (`go:0071941` nitrogen cycle metabolic process, 13 genes) and requires indirect inference for photosynthesis/ribosome anchors via broader process terms (`carboxylic acid metabolic process`, `gene expression`). KEGG's orthology-based pathway-map framing is also genuinely complementary to cyanorak_role's functional-role framing, making cross-ontology agreement (spec §5 Step 4 M3) more informative.
Impact: Pick 2 = `kegg` L2.

### Decision (researcher-approved 2026-04-20)

Selected ontologies (2):
1. **`cyanorak_role` level 1** — MED4-optimal, 4 canonical anchors at this level, cyano-specific.
2. **`kegg` (pathway) level 2** — orthogonal pathway-map framing, bundles N-metab + photosynthesis + ribosome + AA-metab.

Key-pathway panel (11 terms): 6 cyanorak L1 + 5 kegg L2. All term IDs validated via `genes_by_ontology(term_ids=..., organism="MED4")`; canonical marker genes (`glnA`, `cynA-S`, `psbA`, `psaA`, `rbcL`, `rplA`) confirmed present. Panel serves as the biological anchor for Step 2 QC (R clusters must show expected directions). Direction assignments: 5 UP (N-metab × 2, Glu-family AA, Ala/Asp/Glu, AA biosynth family), 5 DOWN (PSI, PSII, CO2 fixation, ribosome × 2, photosynthesis combined), 1 ambiguous (Gly/Ser/Thr).

Step 1b locked. Gate 1 (step-boundary) satisfied — proceed to Step 2 do (Task 5).

Rejected alternatives (all documented in `ontology_selection.md`): tigr_role L0 (redundant), go_bp L3 (sparse N annotation), go_bp L3-5 via term_ids (DAG redundancy + no new N-signal), BRITE subtrees (low coverage), 3rd ontology (reduces cross-ontology agreement informativeness).

MCP / skill-friction items captured:
- [pathway_enrichment `term_ids` limitation claim — retracted](../gaps_and_friction.md#2026-04-20--step-1b-pathway_enrichment-has-no-term_ids-filter-dag-ontology-friction--retracted-my-error) — anti-hallucination meta-pattern (same class as the Hennon→Aharonovich correction earlier).
- [pathway_enrichment UX refinement](../gaps_and_friction.md#2026-04-20--step-1b-pathway_enrichment-level-only-mode-for-dag-ontologies-is-a-ux-refinement-opportunity-minor) — example + doc emphasize `level`-only pattern; for DAG ontologies `term_ids` would be the better first-reach.
- [GO-aware enrichment design notes](../gaps_and_friction.md#2026-04-20--step-1b-design-notes-for-future-go-aware-pathway-enrichment-captures-why-we-did-not-use-go-bp-for-b2) — feature-scope thoughts for a future DAG-aware `pathway_enrichment` mode.

Proceeds to Step 2 do (Task 5) after Task 4 formal decide gate.

---

## 2026-04-20 15:20 — Step 2 do: enrichment + explore phase (IN PROGRESS — handoff to fresh session)

### Commands run
```bash
uv run scripts/03_run_enrichment.py              # main Step 2 do
uv run scripts/explore_step2_rkey_agreement.py   # explore-phase R×keypath drill-down
```

### Outputs committed as 3511318 (do-phase) + pending commit (explore artifacts)
- [data/enrichment_all.csv](../data/enrichment_all.csv) — 11,239 rows, 70 clusters, **225 p_adjust<0.05**
- [data/enrichment_results.pkl](../data/enrichment_results.pkl) — 8-entry dict (organism, ontology, bg) → EnrichmentResult
- [exploration/qc/step2_key_pathway_heatmap.png](qc/step2_key_pathway_heatmap.png) — diagnostic heatmap (crammed; re-render pending)
- [exploration/qc/step2_rkey_matrix.csv](qc/step2_rkey_matrix.csv) + 4 summary CSVs (by_term, by_expTp, disagreements, NC)

### Significance + agreement (from CSVs, not memory)

Per (org, ontology, bg) — from `logs/step2.log`:

| org | ontology | bg | clusters | tests | significant |
|---|---|---|---|---|---|
| MED4 | cyanorak | organism | 8 | 552 | 17 |
| MED4 | cyanorak | table_scope | 51 | 3309 | 78 |
| MED4 | kegg | organism | 8 | 776 | 13 |
| MED4 | kegg | table_scope | 51 | 4695 | 72 |
| MIT9313 | cyanorak | table_scope | 9 | 702 | 23 |
| MIT9313 | kegg | table_scope | 9 | 873 | 17 |
| SS120 | cyanorak | organism | 2 | 144 | 2 |
| SS120 | kegg | organism | 2 | 188 | 3 |
| **total** | | | **70** | **11,239** | **225** |

R cluster × key-pathway agreement per term — from `step2_rkey_summary_by_term.csv`:

| term_id | expected | n_sig | n_agree | min_padj |
|---|---|---|---|---|
| cyanorak:E.4 N-metab | up | 6 | 6 | 4e-5 |
| cyanorak:J.2 CO2 fixation | down | 6 | 5 | 1.5e-6 |
| cyanorak:J.7 PSI | down | 5 | 5 | 2.1e-9 |
| cyanorak:J.8 PSII | down | 4 | 4 | 1.9e-3 |
| cyanorak:K.2 Ribosome | down | 5 | 5 | 1.6e-30 |
| ko00195 Photosynthesis | down | 5 | 5 | 2.0e-18 |
| ko00910 N-metab | up | 5 | 5 | 6.6e-4 |
| ko03010 Ribosome | down | 5 | 5 | 6.2e-34 |

**Total agreement: 40/41 significant R×key-pathway hits concordant with expected (97.6%).** [KG]

**3 key pathways had 0 significant R hits** (signature weakness): `cyanorak.role:A.3`, `ko00250`, `ko00260`. [interpretation] AA biosynthesis is likely suppressed rather than induced under N-limit; the UP expected direction I assigned in `key_pathways.csv` was biologically naive.

### One R disagreement

From `step2_rkey_disagreements.csv`: Read 2017 3h|up × `cyanorak.role:J.2` CO2 fixation, signed=+1.35, padj=0.045 (borderline), expected=down. Early-timepoint transient flip — exactly what the spec §5 Step 3 temporal filter (`timepoint_hours < 3`) is meant to catch; 3h is on the boundary. **Open decision (see below).**

### NC noise floor — interpretable biology, not pure noise

From `step2_rkey_nc_enrichment.csv` (4 significant NC hits):
- Steglich high-light 45min down × PSI/photosynthesis (padj 3e-3, 3e-2) — real high-light-stress response, not noise.
- Weissberg coculture (no N-starvation) day 11 up × N-metabolism (padj 9e-5, 1.6e-5) — real coculture-induced N-scavenging without explicit starvation.

[interpretation] NC mean+SD will be non-zero on photosynthesis (Steglich) and N-metabolism (Weissberg coculture) anchors. Signature score calibration will be biased — T scores must be interpreted against this elevated noise floor.

### OPEN DECISIONS — resolve in fresh session before Step 2 decide commit

1. **Temporal filter boundary.** Keep spec `< 3h` (retains Read 3h disagreement) or tighten to `≤ 3` (excludes). Recommendation: keep spec default.
2. **Heatmap quality.** Re-render with two-panel by ontology + larger fontsize, or accept current diagnostic for QC purpose and rely on Step 5 Fig 1 for publication quality.
3. **Signature weakness for AA-biosynthesis anchors.** Update `key_pathways.csv` expected directions? (current: A.3 up, ko00250 up → signal suggests neither direction is reliably enriched).
4. **NC noise-floor bias.** Document in `caveats.md` that NC floor is non-zero on photosynthesis + N-metabolism anchors. Decide whether to flag specific T scores differently as a result.

Fresh session picks up here: read this entry, read the CSVs listed above, resolve decisions with researcher, commit Step 2 decide gate, then proceed to Step 3 (signature derivation, Task 7).

---


