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

