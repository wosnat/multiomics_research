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

## 2026-04-20 20:30 — Step 2 decide: 4 decisions resolved, extended explore phase in fresh session

### Summary (read this first)

**Step 2 decide: all 4 open decisions resolved; Step 2 QC passed; ready for Step 3.**

Decisions locked (full detail in "Decisions resolved" block below and in [`../decisions.md`](../decisions.md)):
- **D1** — Temporal filter tightened to `hours > 3` (excludes 3h; 12 signature-eligible R clusters remain).
- **D2** — v2 heatmap locked: two-panel per ontology, class dividers + axenic/coculture T split, ±5 display cap + saturation stars, inline legend.
- **D3** — AA-biosynthesis anchors (`A.3`, `ko00250`) kept as `expected_direction="up"` (falsified prior preserved; caveats.md at Step 5).
- **D4** — NC calibration option (b) with padj<1e-3 exclusion criterion; Weissberg coculture d11 up excluded from `(*, table_scope)`; Steglich retained.

New Step 3 sub-task added (plan Task 8 Step 2b): within-ontology pathway-gene-overlap audit — researcher decides A (audit-only) / B (post-filter) / C (de-weight).

Preliminary hypothesis H1 ([`../hypotheses.md`](../hypotheses.md)): coculture > axenic on N-limitation signature strength. **PRELIMINARY — requires Step 4 scoring with LOO + NC-exclusion-sensitivity check before it becomes a finding.**

6 explore scripts + 5 new `exploration/qc/` artifacts committed (`cc61de1` step 2 decide, `160f10b` follow-ups).

**Fresh-session agents resuming: start by reading this summary, then [`../decisions.md`](../decisions.md), then plan Task 7. The full Q&A and exploration detail below is the evidence base — drill in when needed.**

### Entry scope

Fresh session resumed from commit `02ef861`. Read handover memory, notebook through interim Step 2 entry, `gaps_and_friction.md`, `ontology_selection.md`, and all `exploration/qc/step2_rkey_*.csv` artifacts from files (not memory). Extended the explore phase before committing decide.

### Extended exploration — scripts produced during this session

All in `scripts/`, outputs in `exploration/qc/`:

| script | output | purpose |
|---|---|---|
| [explore_step2_cluster_totals.py](../scripts/explore_step2_cluster_totals.py) | [step2_cluster_totals.csv](qc/step2_cluster_totals.csv) | Total significant hits per cluster across full ontology (disambiguates "N sig" in by_expTp which was counted over 11 key pathways only) |
| [explore_step2_r_top_pathways.py](../scripts/explore_step2_r_top_pathways.py) | [step2_r_top_pathways.csv](qc/step2_r_top_pathways.csv) | Main enriched pathways in R clusters at `timepoint_hours > 3` — preview of Step 3 signature eligibility |
| [explore_step2_probe_schema.py](../scripts/explore_step2_probe_schema.py) | stdout | One-off: probe EnrichmentResult public API (confirmed `term2gene` is a DataFrame, `overlap_genes` returns GeneRef objects) |
| [explore_step2_atpsynth_nmetab.py](../scripts/explore_step2_atpsynth_nmetab.py) | stdout | `.explain()` + `overlap_genes` drill-down on J.1 ATP synthase, ko00190 Ox phos, E.4 N-metab; cross-pathway gene-set Jaccards |
| [explore_step2_score_distribution.py](../scripts/explore_step2_score_distribution.py) | [step2_score_distribution.png](qc/step2_score_distribution.png) | Distribution of signed_score across all 11,239 tests — informs SCORE_CAP choice |
| [explore_step2_heatmap_v2.py](../scripts/explore_step2_heatmap_v2.py) | [step2_heatmap_cyanorak_role.png](qc/step2_heatmap_cyanorak_role.png), [step2_heatmap_kegg.png](qc/step2_heatmap_kegg.png) | Publication-style two-panel heatmap per ontology, class-grouped columns, key ∪ discovered rows |

### Chat exploration (Q → data → finding)

**Q1: "Tolonen 3h up → 2 sig, 2 agree" — only 2 significant terms in all enrichment, or something else?**

Data: `step2_cluster_totals.csv` — Tolonen 3h up has 2 sig cyanorak + 1 sig kegg = **3 total significant** ontology terms across 166 tested (69 cyanorak + 97 kegg). The "2" in the rkey summary referred to 2 of the 11 key pathways specifically.
Finding: [KG] the by_expTp summary's `n_significant` column was scoped to the 11-pathway key panel, not whole-ontology significance. Both numbers are real; different denominators.
Impact: clarified semantics in the Step 2 decide narrative; no change to data.

**Q2: What are the main enriched pathways in R clusters at `timepoint_hours > 3`?**

Data: `step2_r_top_pathways.csv`. Top pathways by n_clusters in 12 signature-eligible R clusters (2 exps):

cyanorak_role:
| term | dir | n/12 | mean \|signed\| | exps | in key panel? |
|---|---|---:|---:|---|---|
| **J.1 ATP synthase** | down | **6** | 7.4 | both | **no** |
| K.2 Ribosomes | down | 5 | 24.4 | both | yes |
| J.7 PSI | down | 5 | 6.2 | both | yes |
| J.2 CO2 fixation | down | 5 | 4.6 | both | yes |
| E.4 N metabolism | up | 5 | 3.8 | both | yes |
| J.8 PSII | down | 4 | 2.2 | both | yes |
| D.1 Adaptation/acclimation | up | 3 | 3.9 | both | no (catch-all) |

kegg:
| term | dir | n/12 | mean \|signed\| | exps | in key panel? |
|---|---|---:|---:|---|---|
| **ko00190 Oxidative phosphorylation** | down | **6** | 2.3 | both | **no** |
| ko03010 Ribosome | down | 5 | 27.1 | both | yes |
| ko00195 Photosynthesis | down | 5 | 14.1 | both | yes |
| **ko00710 Calvin cycle** | down | 5 | 4.2 | both | **no** |
| ko00910 N metabolism | up | 4 | 2.4 | Tolonen only | yes |
| ko01200 Carbon metabolism | down | 4 | 1.9 | both | no |

Finding: [KG] 5 pathways beyond the a priori key panel reach `n ≥ 3` same-direction cluster support: J.1, D.1 (cyanorak); ko00190, ko00710, ko01200 (kegg). They will likely enter the Step 3 signature automatically. One Tolonen-only weakness: kegg `ko00910` N-metab (Read's 6-gene kegg N-metab doesn't survive padj<0.05 but its 28-gene cyanorak E.4 N-metab does).
Impact: [interpretation] discovered-strong pathways are captured by the Step 3 derivation mechanism, not by amending the a priori key panel — see decision #3 discussion.

**Q3: `.explain()` drill-down on J.1 ATP synthase, ko00190 Ox phos, E.4 N-metab. What genes drive them? Do they share genes?**

Data: `explore_step2_atpsynth_nmetab.py` stdout. Tolonen 12h-down cluster:

- **J.1 ATP synthase** (10/10 MED4 atp genes): atpA-I operon, log2fc −2.2 to −4.7, padj 1.68e-9, fold 10.48. **Every atp gene coordinately down.**
- **ko00190 Ox phos** (11/33 genes): same 9 atp genes + ppa (pyrophosphatase). MED4 ox phos's 11 ndh* (NADH dehydrogenase) and 4 cta* (cyt c oxidase) **don't hit at 12h down** — ko00190 significance is entirely driven by atp operon.
- **E.4 N-metab** (12/28 genes, 12h up cluster): cynA/B/D/S (cyanate), glnA (GS), ntcA (N-response TF), glnB (PII), ureA-G (urea catabolism), urtA-D (urea ABC transporter). Canonical N-scavenging cassette.

Cross-pathway gene-set overlap (MED4 membership, not cluster-specific):

| pair | overlap | meaning |
|---|---|---|
| J.1 ∩ E.4 | ∅ | N-metab disjoint from ATP synthase — independent evidence |
| J.1 ⊂ ko00190 | 9 atp genes | ATP synthase is the enriched subset of oxphos |
| J.1 ⊂ ko00195 | 9 atp genes | ATP synthase also in KEGG Photosynthesis map |
| ko00190 ∩ ko00195 | 9 atp genes (only) | oxphos and photosynthesis share ONLY the atp operon |
| E.4 ∩ ko00195 | ∅ | N-metab independent of photosynthesis machinery |

Finding: [KG] double-counting risk in the Step 3 signature if J.1, ko00190, ko00195 all enter — same 9 atp genes count three times in Layer A scoring. **Per-ontology scoring** (spec §5 Step 4) mitigates cross-ontology overlap, but within-kegg overlap (ko00190 ∩ ko00195) is real.
Impact: added new Step 3 explore task to audit within-ontology pairwise Jaccard on signature pathways. Decision between "M4 audit in Step 4" and "post-filter in Step 3" deferred to Step 3 explore (see Decisions below).

**Q4: Is redundancy covered by LOO, or separate?**

Data: spec §5 Step 4 M3 (LOO-pathway, LOO-R-experiment, cross-ontology agreement).
Finding: [interpretation] LOO catches **fragility** (single-pathway or single-experiment dominance). Redundancy catches **inflation** by correlated gene overlap. Orthogonal failure modes. If J.1, ko00190, ko00195 all vote "atp down" in a T cluster, LOO shows the score is robust but the robustness is illusory — three pathways agreeing because they share 9 genes, not because three independent biological events confirm each other.
Impact: within-ontology overlap is a separate check, added as Step 3 explore task (see skill-friction note at end).

**Q5: Should we amend the key_pathways.csv list based on discovered-strong pathways (J.1 ATP synthase, ko00190 Ox phos, ko00710 Calvin cycle)?**

Data: key_pathways.csv was locked at Step 1b as an *a priori* biological anchor list for Step 2 QC — "did the signal behave at canonical N-response markers?"
Finding: [interpretation] adding discovered-strong pathways after seeing them hit would fit the anchor to the data — confirmation bias on the validation step. The signature (Step 3 derivation) is the right mechanism to surface discovered pathways. Keep key_pathways.csv stable.
Impact: decision — do NOT add J.1, ko00190, ko00710. They'll enter via Step 3 `n_down ≥ 3` rule. Logged as skill-friction meta-pattern (see below).

**Q6: `signed_score = -log10(padj)` — is the difference between 5→10 and 10→20 meaningful, or just noise?**

Data: [step2_score_distribution.png](qc/step2_score_distribution.png) — |signed_score| quantiles on significant rows (N=225):

| quantile | \|s\| |
|---:|---:|
| p50 | 2.92 |
| p75 | 5.43 |
| p90 | 15.70 |
| p95 | 24.13 |
| p99 | 31.49 |
| max | 33.21 |

15% of significant hits exceed \|s\| = 10. Top hits all ribosome-DOWN (ko03010 / K.2) at R or CTX clusters, |s|=23–33.

Finding: [KG + interpretation] `signed_score` is log-scale evidence strength. Three regimes:
- **0–3**: below threshold (padj > 1e-3) — noise; excluded by padj<0.05 anyway.
- **3–10**: biologically meaningful (padj 1e-3 to 1e-10) — cluster-to-cluster differences reflect real coverage/fold differences.
- **10+**: statistical saturation. In Fisher ORA, at this scale padj depends on integer overlap count; one additional gene hit can drop padj by many orders of magnitude. Differences are measurement-precision artifacts, not biology.

Example: ribosome at Cluster A: 50/54 ribosomal proteins hit → padj 1e-25; Cluster B: 54/54 → padj 1e-33. Same biology ("ribosome shut down"), 8-point |s| difference = 4-gene detection noise.

Impact: [interpretation] spec's ±10 SCORE_CAP is principled. The tail |s|>10 is statistical saturation, not biology. Capping prevents single-pathway mean-domination of Layer A scores without losing interpretive signal.
Decision: **keep spec ±10 for Layer A scoring**. Use **±5 display cap** for visualization (full color range on biologically-meaningful 0–5 band) + saturation stars (`*` ≥5, `**` ≥10, `***` ≥20) to preserve "how much beyond cap" info.

**Q7: Axenic vs coculture in T panel — is this the biologically meaningful contrast?**

Data: [step2_heatmap_cyanorak_role.png](qc/step2_heatmap_cyanorak_role.png) + [step2_heatmap_kegg.png](qc/step2_heatmap_kegg.png), T panel split axenic | coculture.
Finding: [KG + interpretation] coculture T clusters show STRONGER N-limitation signature than axenic:
- Coculture Prot day 31 is the loudest T cluster — ko03010 Ribosome `***` (|s|=23.6, padj 2.3e-24), ko00910 N-metab `*` red.
- Coculture band on kegg N-metab red across days 14/18/31/60+89; axenic band faint/absent.
- Cyanorak E.4 N-metab strongest in coculture days 18+31+60+89 (Prot).

Against the memory-held biological expectation: "MED4 dies axenically in ~2 weeks, survives 90 days in coculture." Naively, axenic (dying under N-limit) should show stronger N-limit signature. But data shows the opposite: coculture has the stronger and cleaner signature.

Impact: [interpretation] **hypothesis worth pursuing through Step 4–5**: the N-limit response program in MED4 is *actively engaged* by coculture — axenic cells shut down generally without engaging specific N-scavenging circuitry before dying, whereas coculture enables sustained engagement of the N-response. The paper-worthy story, if it survives Step 3+4 scoring. Documented here for the methodology→interpretation trail.

### Decisions resolved (all 4)

**Decision #1 — Temporal filter boundary.** Tighten spec `timepoint_hours < 3` to `> 3` (exclude 3h from signature derivation). Researcher clarified original intent was to drop all early TPs; the strict `<` wording in spec §5 Step 3 was imprecise. Signature-eligible R clusters: 12 (down from 16 under `< 3`). Core rule `n ≥ 3` remains comfortable. Read 3h-up × J.2 CO2 fixation disagreement now excluded cleanly.
- Spec refinement to apply in Task 7: `TIMEPOINT_HOURS_CUTOFF = 3.0` with keep-mask `hours > cutoff` (not `>=`).

**Decision #2 — Heatmap re-render.** Build publication-quality v2 heatmap per ontology with:
- Rows: key panel ∪ discovered-strong pathways (`n_sig ≥ 3` in R clusters at `hours > 3`).
- Columns: all 70 clusters per ontology, grouped by class (R | PC | CTX | NC) with dividers.
- T on second row panel below, split axenic | coculture (biological contrast, not Prot/RNA).
- Display cap ±5 + saturation stars (`*`/`**`/`***`).
- Author names truncated to 6 chars (Weissb, Tolone, Read, Aharon, Moreno, Stegli).
- Uniform cell sizes across both panels via explicit axis positioning.
- Inline legend at top of figure spelling out the conventions (see legend block below).

Artifacts: [step2_heatmap_cyanorak_role.png](qc/step2_heatmap_cyanorak_role.png), [step2_heatmap_kegg.png](qc/step2_heatmap_kegg.png).

#### Heatmap legend (identical on both panels, also rendered atop each figure)

| glyph | meaning |
|---|---|
| **color** (cell) | `signed_score = sign × −log10(padj)` — red = up, blue = down. Display color-capped at ±5 (preserves biology-range resolution 0–5; cells above are all deepest shade). |
| **blank cell** | term was not tested in that cluster (either `min_gene_set_size` filter, or cluster failed at the `pathway_enrichment` call) |
| **`*`** in cell | uncapped \|signed_score\| ≥ 5 (padj ≤ 1e-5) — strong |
| **`**`** in cell | uncapped \|signed_score\| ≥ 10 (padj ≤ 1e-10) — saturated |
| **`***`** in cell | uncapped \|signed_score\| ≥ 20 (padj ≤ 1e-20) — statistical ceiling, differences at this level are detection-precision artifacts (see Q6) |
| `↑` / `↓` / `~` beside row label | a priori expected direction for that key pathway (up / down / ambiguous). Rows without a marker are discovered-strong pathways, not a priori anchors. |
| **bold row label** | member of a priori key panel (`exploration/key_pathways.csv`) |
| regular row label | discovered-strong pathway (n_sig ≥ 3 in R clusters at hours > 3 — signature-eligibility threshold) |
| **vertical dividers** (non-T panel) | class boundaries: R \| PC \| CTX \| NC |
| **vertical divider** (T panel) | axenic \| coculture (Weissberg 2025 biological contrast) |
| column label format | `author6|tp_short` (non-T) or `author6|omics|tp_short` (T). Authors truncated to 6 chars: Weissb, Tolone, Read, Aharon, Moreno, Stegli, Domíng. Timepoints: `day 14 → 14d`, `steady state → ss`, `darkness → dark`, etc. |

(The scoring cap for Layer A in Step 4 remains the spec's ±10. The ±5 visualization cap here is for diagnostic readability only — see gaps_and_friction.md "signed_score cap needs explicit methodology-layer reasoning" for the distinction.)

**Decision #3 — AA-biosynthesis anchors.** Keep `cyanorak.role:A.3` and `kegg.pathway:ko00250` as `expected_direction="up"` in `key_pathways.csv` (option c: preserve original prior as falsified hypothesis). Add caveats.md entry at Step 5 documenting the up-prediction was biologically naive — AA biosynthesis under N-limit is typically suppressed, not induced. `ko00260` stays `ambiguous` unchanged.
Rationale: changing expected direction to match observation after the fact is post-hoc. Leaving as `up` with caveat preserves the honest scientific record (prior → observation → falsification).

**Decision #4 — NC noise-floor bias.** Option (b) chosen: **exclude ontology-misclassified NC clusters from NC calibration, with an explicit programmatic criterion.**

Rationale: the spec's Step 4 M3 calibration computes `nc_mean` and `nc_std` across all NC clusters, assuming NC is true noise. When specific NC clusters carry real biology that overlaps with N-limit signature markers, they inflate the calibration threshold, silently biasing T-score classification toward "no signal." Option (a) (document-only) leaves the biased threshold in place; option (b) corrects it with a transparent rule; option (c) (side-by-side both calibrations) doubles the reporting burden without improving interpretability.

**Exclusion rule (programmatic, for Step 4 implementation):**

> An NC cluster is excluded from the `(ontology, background_used)` NC calibration group if it shows significant enrichment at padj < 1e-3 on any term in the a priori key-pathway panel for that ontology.
>
> The padj < 1e-3 threshold is deliberately stricter than the standard 0.05 significance bar — at 1e-3, the enrichment is 3 orders of magnitude below nominal significance and implausibly attributable to sampling noise alone. Borderline hits (1e-3 ≤ padj < 0.05) stay in NC calibration on the principle that true noise floors include noise-adjacent fluctuations.

Applied to Step 2 NC hits (from `step2_rkey_nc_enrichment.csv`), 4 significant NC × key-pathway rows fall across 2 distinct NC clusters:

| NC cluster | sig key-pathway hits (padj) | threshold check | calibration status |
|---|---|---|---|
| Steglich high-light 45min \| down | cyanorak J.7 PSI (3.2e-3), kegg ko00195 Photosynthesis (3.0e-2) | both ≥ 1e-3 | **keep in NC calibration** (both ontologies) — borderline noise-adjacent |
| Weissberg coculture d11 \| up | cyanorak E.4 N-metab (8.6e-5), kegg ko00910 N-metab (1.6e-5) | both < 1e-3 | **EXCLUDE from NC calibration** for (cyanorak_role, table_scope) AND (kegg, table_scope) — Weissberg NC is `all_detected_genes` → table_scope bg |

The exclusion is per-cluster, not per-experiment: the Weissberg d11 **up** cluster is excluded (strong N-metab signal); the Weissberg d11 **down** cluster is retained (no significant key-pathway hits).

Impact on NC calibration group sizes (verified from `step2_cluster_totals.csv` + `step2_rkey_nc_enrichment.csv`):

| (ontology, background_used) | NC clusters before | excluded | **NC clusters for calibration** |
|---|---:|---:|---:|
| cyanorak_role × table_scope | 8 | 1 (Weissberg d11 up) | **7** |
| cyanorak_role × organism | 2 | 0 | 2 |
| kegg × table_scope | 8 | 1 (Weissberg d11 up) | **7** |
| kegg × organism | 2 | 0 | 2 |

[caveat] The `(*, organism)` groups only have 2 NC clusters each (both from Steglich high-light). 2-cluster SD is extremely fragile — any single extreme value dominates. T clusters scored against `organism` background (Read R2, Domínguez-Martín CTX3, and potentially some Weissberg axenic-protein T clusters if they carry `significant_only` tables) have very weak NC-noise anchoring. Flag in T-score interpretation at Step 4 explore. No action here beyond noting; the underlying constraint is that N-related organism-bg experiments are scarce in the KG.

[interpretation] Steglich shows a real HL-stress response on photosynthesis but it's borderline — conservative to keep it in NC. Weissberg coculture day 11 with no explicit N-starvation showing strong N-metab upregulation is real biology ("coculture-enabled N-scavenging"), clearly misclassified as noise for N-metab calibration purposes.

**CRITICAL FOLLOWUP — Task 10 (Step 4 scoring) implementation requirement:**

The `05_compute_scores.py` script at Task 10 MUST implement this filter when computing `nc_mean_{o,b}` and `nc_std_{o,b}`:

```python
# Pseudocode for the NC calibration filter:
KEY_PATHWAY_EXCLUSION_PADJ = 1e-3

for (ontology, bg), nc_group in nc_clusters.groupby(["ontology", "background_used"]):
    # Drop NC clusters with any key-pathway hit at padj < 1e-3 in this (o, b) group.
    key_ids_this_ontology = set(key_panel.loc[key_panel["ontology"] == ontology, "term_id"])
    disqualifying_hits = enrichment_all[
        (enrichment_all["cluster"].isin(nc_group["cluster"])) &
        (enrichment_all["term_id"].isin(key_ids_this_ontology)) &
        (enrichment_all["p_adjust"] < KEY_PATHWAY_EXCLUSION_PADJ)
    ]
    disqualified_clusters = set(disqualifying_hits["cluster"])
    nc_calibration_group = nc_group[~nc_group["cluster"].isin(disqualified_clusters)]
    # Log the exclusion to step4.log and to decisions.md per the spec.
    # Compute nc_mean and nc_std on nc_calibration_group, not the original nc_group.
```

**Pre-committed caveats.md text** (collated into `caveats.md` at Task 12 / Step 5 do):

> **C5 — NC calibration exclusion for biology-contaminated clusters.** One NC cluster was excluded from NC calibration in Step 4 per the "ontology-misclassified NC" rule (padj<1e-3 criterion, see notebook Step 2 decide, decision #4): the Weissberg 2025 coculture day-11 **up** cluster was excluded from the `(cyanorak_role, table_scope)` and `(kegg, table_scope)` NC calibration groups due to strong N-metab enrichment at padj ~1e-5 (real coculture-induced N-scavenging biology, not noise). The matched Weissberg d11 **down** cluster shows no significant key-pathway hits and remains in calibration. The Steglich 2006 high-light 45min down cluster was **kept** in NC calibration for both `(cyanorak_role, organism)` and `(kegg, organism)` — its PS-down hits at padj 3e-3 / 3e-2 are near the noise-adjacent threshold (padj ≥ 1e-3) and conservatively retained. Resulting NC calibration group sizes: `(*, table_scope)` = 7 clusters (was 8); `(*, organism)` = 2 clusters (unchanged). T-score calibration thresholds (`nc_mean ± 2σ`) are computed against these filtered groups. Interpretive impact: T clusters scoring just above `nc_mean + 2σ` on N-metab anchors should be read with this exclusion in mind — the threshold is cleaner than it would have been with Weissberg coculture d11 up in the calibration, but the `(*, organism)` calibration with only 2 clusters is statistically weak regardless. Readers should treat threshold classifications as narrative indicators, not hypothesis-test conclusions.

### Gate checks

### New Step 3 task: within-ontology pathway redundancy audit

Before Task 8, Step 3 `explore` phase gains a new sub-task:
> For each signature pathway pair within the same ontology, compute Jaccard(pathway_i.term2gene, pathway_j.term2gene) and flag ⊂ relationships or Jaccard > 0.5. Present as heatmap / flag list. Researcher decides between: (A) audit-only reported as M4 stability check in Step 4, (B) post-filter in Step 3 (keep most-enriched of each overlap cluster), or (C) caveats-only.

Known upfront concern from Step 2 exploration: within-KEGG, `ko00190 Oxphos ∩ ko00195 Photosynthesis = 9 atp genes`, and J.1 ATP synthase (cyanorak, 10 atp genes) is effectively within-scope if kegg carries ko00190+ko00195 both. Within-cyanorak: J.1/J.2/J.7/J.8/K.2 are hierarchically-disjoint L1 categories — likely gene-disjoint, but verify during audit.

### Gate checks

- Gate 1 (step-boundary): Step 2 do-phase output (`enrichment_all.csv`, pickle, heatmap) committed `3511318`. Explore-phase interim artifacts committed `02ef861`. This decide commit is the boundary.
- Gate 2 (manifest currency): `data/DATA_MANIFEST.md` + `exploration/qc/` artifacts current (to be verified in commit).
- Gate 3 (chat-capture): all exploration Q→data→finding→impact captured above. No chat-only reasoning.

### Skill / methodology friction — captured for gaps_and_friction.md

Added in this session (see [../gaps_and_friction.md](../gaps_and_friction.md) Skill/methodology friction section):

1. **SCORE_CAP choice needs explicit methodology-layer reasoning**, not just a spec default. `signed_score = -log10(padj)` has saturation regime above ~|s|=10 where differences are measurement artifacts, not biology. Scoring caps (±10) and visualization caps (±5 here) serve different purposes and should be distinguished.
2. **Within-ontology pathway redundancy is not covered by LOO.** LOO catches fragility (single-pathway dominance); redundancy catches inflation (correlated gene overlap). Orthogonal failure modes. Needs a separate audit.
3. **A priori anchor lists stay locked through QC**. Discovered-strong pathways go into the a posteriori signature, not backported to the key panel. Mixing the two confuses validation with discovery.
4. **Author-label truncation convention** (6-char last names) as figure-convention for multi-cluster heatmaps — captured as practical artifact convention.

### Decision

Step 2 decide gate passed. Proceed to Task 7 (Step 3 do — signature derivation) with:
- Temporal filter `hours > 3` (decision #1).
- New within-ontology redundancy audit sub-task inserted before Task 8.
- caveats.md entries deferred to Step 5 (decisions #3, #4).

---

## 2026-04-21 10:XX — Step 3 do + show + explore + decide: MED4 reference signature

### Summary (read this first)

**Step 3 completed; signature locked.** 13 signature pathways across 2 ontologies, all under the core rule (`n_clusters_supporting ≥ 3`, no fallback) on 12 signature-eligible MED4 R clusters (Tolonen 6h/12h/24h/48h × 2dir + Read 12h/24h × 2dir, after D1 `hours > 3.0` filter).

- **cyanorak_role (7):** 5 down (J.1 ATP synthase, J.2 CO2 fix, J.7 PSI, J.8 PSII, K.2 Ribosome), 2 up (D.1 Adaptation/acclimation, E.4 N-metab).
- **kegg (6):** 5 down (ko00190 Oxphos, ko00195 Photosynthesis, ko00710 Calvin, ko01200 C-metab, ko03010 Ribosome), 1 up (ko00910 N-metab).
- **Key-panel coverage:** 8/11 a priori anchors in signature; 3 silently omitted (cyanorak A.3, ko00250, ko00260) — all D3-documented falsifications of the "AA biosynthesis up under N-limit" prior.
- **Discovered-strong new in signature (not in a priori panel):** 5 — cyanorak D.1, J.1; kegg ko00190, ko00710, ko01200.

**Decisions locked in this step (full detail in [`../decisions.md`](../decisions.md)):**
- **D5** — Watchlist curation: `kegg.pathway:ko05152` (Tuberculosis) excluded from `signature_dropped.csv` via `WATCHLIST_EXCLUDE_TERMS` in the runner. Curation affects only the "terms to watch" shortlist, never the signature. 5 retained watchlist terms (D.4, L.3, R.2, ko00061, ko01212).
- **D6** — Signature stability handling: option (A) audit-only for the 1 flagged within-kegg redundancy (`ko00710 Calvin ⊂ ko01200 C-metab` strict subset) + 1 soft note (`ko00190 ∩ ko00195 = 9 atp genes`). `ko00910 N-metab (up)` flagged as single-R-experiment dominated (Tolonen-only) — forward-flag to Step 4 LOO-R. Task 10 adds `results/kegg_redundancy_sensitivity.csv` artifact; Task 12 adds caveats C6 (redundancy) and C7 (ko00910 dominance).

4 explore scripts committed (dominance + redundancy) with CSV/PNG outputs in `exploration/qc/step3_*`.

**Fresh-session agents resuming Step 4: start from this summary + [`../decisions.md`](../decisions.md) D1–D6, then plan Task 9.**

### Do-phase — [signature.py](../scripts/signature.py) + [04_derive_signature.py](../scripts/04_derive_signature.py)

Primitive module `signature.py` holds pure `derive_for_ontology(df, support_threshold)` (no I/O). Runner `04_derive_signature.py` performs I/O, applies the temporal filter `hours > TIMEPOINT_HOURS_CUTOFF=3.0` (strict `>`, per D1) and then calls the primitive per ontology, with the `n_sig < 5 → fallback` rule wired but unused (both ontologies satisfied core).

**Run sanity check** (`logs/step3.log`):

- 2 MED4 R experiments × 17 total R clusters before filter.
- Temporal filter dropped 830 rows from 5 clusters: Tolonen 0h up+down, Tolonen 3h up, Read 3h up+down. Expected minimum was 3 cluster rows (Tolonen 0h up + 3h up + Read 3h up); actual was 5 (also dropped Tolonen 0h down + Read 3h down as they had enrichment rows). Clean pass.
- 12 signature-eligible clusters remaining = Tolonen 6h/12h/24h/48h × 2dir (8) + Read 12h/24h × 2dir (4). Matches D1 budget.
- cyanorak signature size = 7 core; kegg signature size = 6 core; no fallback triggered.

### Show — composition tables

**Signature distribution per ontology (from `reference_signature.csv`, 13 rows):**

| ontology | up | down | total | rule |
|---|---:|---:|---:|---|
| cyanorak_role | 2 | 5 | 7 | core (≥3) |
| kegg | 1 | 5 | 6 | core (≥3) |

**Per-R-experiment support (from [step3_per_exp_support.csv](qc/step3_per_exp_support.csv)):**

| term | dir | n_exps | Tolonen | Read | share_max_exp |
|---|:---:|:---:|---:|---:|---:|
| cyanorak J.1 ATP synthase | ↓ | 2 | 4 | 2 | 0.67 |
| cyanorak J.2 CO2 fixation | ↓ | 2 | 4 | 1 | 0.80 |
| cyanorak J.7 PSI | ↓ | 2 | 3 | 2 | 0.60 |
| cyanorak J.8 PSII | ↓ | 2 | 3 | 1 | 0.75 |
| cyanorak K.2 Ribosome | ↓ | 2 | 3 | 2 | 0.60 |
| cyanorak D.1 Adaptation | ↑ | 2 | 1 | 2 | 0.67 |
| cyanorak E.4 N-metabolism | ↑ | 2 | 3 | 2 | 0.60 |
| kegg ko00190 Oxphos | ↓ | 2 | 4 | 2 | 0.67 |
| kegg ko00195 Photosynthesis | ↓ | 2 | 3 | 2 | 0.60 |
| kegg ko00710 Calvin cycle | ↓ | 2 | 4 | 1 | 0.80 |
| kegg ko01200 C-metabolism | ↓ | 2 | 3 | 1 | 0.75 |
| kegg ko03010 Ribosome | ↓ | 2 | 3 | 2 | 0.60 |
| **kegg ko00910 N-metabolism** | **↑** | **1** | **4** | **0** | **1.00 ← flag** |

**Dropped-notable (post-D5 curation) = 5 terms** — all n=2, max |s| below 3.3:

| ontology | term | dir | n | max \|s\| |
|---|---|:---:|:---:|:---:|
| cyanorak | D.4 Chaperones | ↓ | 2 | 1.57 |
| cyanorak | L.3 Protein folding and stabilization | ↓ | 2 | 1.80 |
| cyanorak | R.2 Conserved hypothetical proteins | ↑ | 2 | 3.28 |
| kegg | ko00061 Fatty acid biosynthesis | ↓ | 2 | 2.97 |
| kegg | ko01212 Fatty acid metabolism | ↓ | 2 | 2.87 |

(ko05152 Tuberculosis excluded per D5.)

**Key-panel coverage (11 a priori anchors):**

| status | count | terms |
|---|:---:|---|
| In signature | 8/11 | cyanorak E.4↑, J.7↓, J.8↓, J.2↓, K.2↓; kegg ko00910↑, ko00195↓, ko03010↓ |
| Silently omitted (0 R hits, D3 falsified) | 3/11 | cyanorak A.3 (up expected, 0/12 R sig); kegg ko00250 (up expected, 0/12), ko00260 (ambiguous, 0/12) |
| New discovered-strong (not in key panel) | +5 | cyanorak D.1, J.1; kegg ko00190, ko00710, ko01200 |

### Explore — dominance + redundancy

#### [explore_step3_single_exp_dominance.py](../scripts/explore_step3_single_exp_dominance.py) → [step3_single_exp_dominated.csv](qc/step3_single_exp_dominated.csv)

**1 flag:** `kegg.pathway:ko00910 Nitrogen metabolism (up)` — 4 Tolonen 2006 clusters (6h, 12h, 24h, 48h), 0 Read 2017. `share_max_exp = 1.00`.

[interpretation] Confirms Step 2 preview ("Tolonen-only weakness on kegg ko00910"). Read 2017 does hit `cyanorak.role:E.4 Nitrogen metabolism` (3 clusters, padj<0.05), so the broader N-metabolism signal has cross-experiment support in cyanorak. But within kegg, the narrower 6-gene `ko00910` pathway needs more detection power than Read's RNA-seq provides post-BH, so Read's N-metab support is only readable via cyanorak.

[impact] Task 10 M4 LOO-R: expected drop of ko00910 from kegg signature when Tolonen removed; expected retention when Read removed. Recorded as forward-flag in D6 Part 2; caveats C7 at Step 5.

#### [explore_step3_redundancy_audit.py](../scripts/explore_step3_redundancy_audit.py) → [step3_signature_redundancy_{cyanorak_role,kegg}.{csv,png}](qc/)

**cyanorak_role:** 0 flagged pairs. All pairwise Jaccard < 0.08. Hierarchical disjointness at L1 confirmed.

| top-5 cyanorak pairs | n_a | n_b | ∩ | Jaccard |
|---|:---:|:---:|:---:|:---:|
| D.1 Adaptation (213) ∩ E.4 N-metab (28) | 213 | 28 | 17 | 0.076 |
| D.1 ∩ J.8 PSII (31) | 213 | 31 | 4 | 0.017 |
| D.1 ∩ K.2 Ribosome (61) | 213 | 61 | 4 | 0.015 |
| D.1 ∩ J.7 PSI (17) | 213 | 17 | 3 | 0.013 |
| D.1 ∩ J.2 CO2 fix (22) | 213 | 22 | 0 | 0.000 |

D.1 "Adaptation/acclimation to atypical conditions" is a 213-gene cyanorak catch-all — it contains many genes but its overlap with focused pathways is small relative to its own size. D.1 itself carries N-adaptation biology (includes stress-response and regulatory genes) and is biologically distinct from the focal pathways even with small shared gene counts. No filter warranted.

**kegg:** **1 flagged pair + 1 soft note:**

| kegg pair | n_a | n_b | ∩ | Jaccard | relation | flag |
|---|:---:|:---:|:---:|:---:|:---|:---:|
| **ko00710 Calvin cycle ⊂ ko01200 Carbon metabolism** | 16 | 58 | 16 | 0.276 | **a_in_b strict subset** | ✓ |
| ko00190 Oxphos ∩ ko00195 Photosynthesis | 33 | 51 | 9 | 0.120 | none (9 atp genes) | — |

[interpretation]
- **ko00710 ⊂ ko01200:** every one of the 16 MED4 Calvin-cycle genes is also in the Carbon-metabolism umbrella. KEGG annotates ko01200 as covering TCA + glycolysis + gluconeogenesis + Calvin cycle + pentose phosphate; Calvin is one branch of it. Keeping both in the signature means the 16 Calvin genes' evidence contributes to the kegg Layer A aggregate twice.
- **ko00190 ∩ ko00195 = 9 atp genes:** the atpA-I operon. Step 2 Q3 `.explain()` showed this is the ONLY overlap between oxphos and photosynthesis in MED4. Jaccard 0.12 is below the 0.5 flag threshold but biologically consequential — coordinate atp-down ≠ independent confirmation of "oxphos down" AND "photosynthesis down."
- **cyanorak J.1 ATP synthase (10 atp genes)** does not overlap any cyanorak signature term. It's a distinct category at L1 per cyanorak's hierarchical schema.

**Decide — option (A) audit-only** for both kegg overlap findings (see D6 Part 1 for reasoning). Task 10 adds an **M4 redundancy sensitivity check** that computes kegg `score_A(T)` with full signature AND with `ko00710` removed AND with `{ko00710, ko00195}` removed, records deltas in `results/kegg_redundancy_sensitivity.csv`, and flags T-cluster classification flips.

### Decisions resolved (2)

- **D5** — Watchlist curation: ko05152 Tuberculosis excluded from signature_dropped.csv via runner-level constant. 5 watchlist terms retained.
- **D6** — Signature stability handling: option (A) audit-only for within-kegg redundancy (ko00710⊂ko01200 subset + ko00190∩ko00195 atp-operon soft overlap); ko00910 N-metab single-R-exp dominance flagged for Task 10 LOO-R expected-drop check and caveats C7.

### Gate checks

- **Gate 1 (step-boundary):** Step 3 do committed `36aee00`. This decide commit bundles explore scripts + signature_dropped.csv re-generation after D5 + decisions.md D5/D6 append + notebook entry + manifest + qc artifacts.
- **Gate 2 (manifest currency):** DATA_MANIFEST.md updated for reference_signature.csv + signature_dropped.csv (do-phase commit) + the 4 new Step 3 explore qc files (this commit).
- **Gate 3 (chat-capture):** all exploration Q→data→finding→impact captured above. Recommendation framing captured inline; researcher approved option (A) audit-only in single-turn "ok".

### Decision

Step 3 decide gate passed. Signature locked at 13 terms (7 cyanorak + 6 kegg, all core rule). Proceed to Task 9 (Step 4 pre-registration of T outcomes) then Task 10 (scoring + stability, with D4 NC exclusion + D6 M4 redundancy sensitivity + D6 LOO-R ko00910 expected-drop check all wired).

---


