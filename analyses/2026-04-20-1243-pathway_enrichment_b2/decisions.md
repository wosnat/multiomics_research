# Decision Log — Pathway enrichment B2

Formal record of decisions taken during the analysis. Each entry captures what was decided, when, on what evidence, and what it affects downstream. Companion to `exploration/notebook.md` (full reasoning) and `gaps_and_friction.md` (skill-level lessons).

Numbering continues across steps; Task 9 (pre-registration) will add subsequent entries.

---

## D1 — Temporal filter semantics: `hours > 3`, not `hours < 3`

**Decided:** 2026-04-20, Step 2 decide (commit `cc61de1`)
**Supersedes:** spec §5 Step 3 wording `timepoint_hours < 3`

Spec §5 Step 3 uses `timepoint_hours < 3` for the signature-derivation temporal filter. Under that wording, 3h is retained; under the researcher's original intent (drop all early timepoints including 3h), it should be excluded. Reviewed during Step 2 decide when Read 2017 3h × cyanorak.role:J.2 CO2 fixation appeared as the sole borderline R disagreement (padj=0.045, transient early-TP flip); confirmed the original intent was strict early-TP exclusion.

**What changes.** `TIMEPOINT_HOURS_CUTOFF = 3.0` with keep-mask `hours > cutoff` (not `>=`). Signature-eligible R clusters drop from 16 → 12 (Tolonen 6h+12h+24h+48h × 2 dir = 8; Read 12h+24h × 2 dir = 4). Core rule `n_same_direction ≥ 3` remains comfortable (up to 6 up-clusters + 6 down-clusters per term).

**Affects:** Task 7 (`scripts/signature.py` / `04_derive_signature.py` runner) — keep-mask must be `>`, not `>=`.

---

## D2 — Heatmap design v2: two-panel per ontology, display-cap ±5, saturation stars, axenic/coculture T split

**Decided:** 2026-04-20, Step 2 decide (commit `cc61de1`)
**Supersedes:** initial single-axis heatmap `step2_key_pathway_heatmap.png` (crammed, illegible)

Step 2 QC heatmap needs to separate concerns: (i) a priori anchors vs a posteriori discovered; (ii) R/PC/CTX/NC reference/calibration clusters vs T target clusters; (iii) biological dynamic range vs statistical-saturation tail.

**Design locked (see `scripts/explore_step2_heatmap_v2.py`):**

- **Rows** = a priori key panel (bold) ∪ discovered-strong pathways (`n_sig ≥ 3` in signature-eligible R clusters, regular font).
- **Columns** per ontology, two stacked panels:
  - Top: non-T clusters grouped `R | PC | CTX | NC` with vertical dividers and class labels.
  - Bottom: T clusters grouped `axenic | coculture` with divider and condition labels (biological contrast; omics type is in the cell label as `Prot` / `RNA`).
- **Color scale**: diverging RdBu_r, display cap ±5 (preserves resolution on biologically-meaningful 0–5 band). Scoring cap stays ±10 per spec §5 Step 4 M2 (separate concern).
- **Saturation stars** on cells: `*` ≥ \|s\|=5, `**` ≥ \|s\|=10, `***` ≥ \|s\|=20. Preserves "beyond cap by how much" without visual domination by saturated cells.
- **Labels**: authors truncated to 6 chars (Weissb, Tolone, Read, Aharon, Moreno, Stegli, Domíng); timepoints abbreviated (`day 14 → 14d`, `steady state → ss`, `45 min → 45m`, etc.).
- **Cell sizing** uniform across both panels via explicit figure-coord axis positioning.
- **Inline 3-line legend** above the title spelling out cell color, blank-cell meaning, star tiers with padj thresholds, row-label direction markers.

**Affects:** Step 5 Fig 1 will reuse this design with minor refinements (final signature rows, maybe per-panel normalization).

---

## D3 — AA-biosynthesis anchors (A.3, ko00250) kept as `expected_direction="up"` despite falsified prior

**Decided:** 2026-04-20, Step 2 decide (commit `cc61de1`)
**Alternative rejected:** Option (a) change to `ambiguous`; option (d) remove from panel.

Three key-pathway anchors had 0 significant R hits in the enrichment data:
- `cyanorak.role:A.3` AA biosynthesis Glu family (16 genes, expected `up`) — 0/12 eligible R clusters hit at padj<0.05.
- `kegg.pathway:ko00250` Ala/Asp/Glu metabolism (15 genes, expected `up`) — 0/12 R clusters hit.
- `kegg.pathway:ko00260` Gly/Ser/Thr metabolism (25 genes, expected `ambiguous`) — stayed ambiguous as predicted.

The first two were wrong-direction predictions: under N-limit, AA biosynthesis is typically suppressed (protein synthesis throttled), not induced. The up-direction prediction was biologically naive.

**Option (c) chosen: leave expected_direction as `up`**, add a caveats.md entry at Step 5 documenting the falsification. Changing the prediction post-hoc to match the observation would erase the honest record (prior → observation → falsification). Keeping it falsified-visible is informative to readers.

**Affects:** `exploration/key_pathways.csv` unchanged. Task 12 Step 5 caveats.md must include an entry noting the falsification — see notebook Step 2 decide for pre-committed text.

---

## D4 — NC calibration excludes ontology-misclassified NC clusters, padj<1e-3 criterion

**Decided:** 2026-04-20, Step 2 decide (commit `cc61de1`)
**Alternative rejected:** Option (a) document-only; option (c) side-by-side both calibrations.

Spec §5 Step 4 M3 defines NC noise floor as `nc_mean_{o,b} = mean(nc_scores_{o,b})` across all NC clusters in each `(ontology, background_used)` group, assuming NC is true noise. Step 2 exploration found 4 significant NC × key-pathway enrichment rows coming from 2 distinct NC clusters that carry real biology (not noise):

- **Weissberg 2025 coculture day 11 up** cluster — cyanorak E.4 N-metab (padj 8.6e-5) and kegg ko00910 N-metab (padj 1.6e-5). Real coculture-enabled N-scavenging under nutrient-replete conditions.
- **Steglich 2006 high-light 45min down** cluster — cyanorak J.7 PSI (padj 3.2e-3) and kegg ko00195 Photosynthesis (padj 3.0e-2). Real high-light PS-stress response.

Without filtering, these inflate NC mean+SD on PS and N-metab anchors and bias T-score classification toward false "no signal."

**Option (b) chosen — programmatic exclusion with padj<1e-3 criterion.**

> An NC cluster is excluded from the `(ontology, background_used)` NC calibration group if it shows significant enrichment at padj < 1e-3 on any term in the a priori key-pathway panel for that ontology. padj<1e-3 is deliberately stricter than 0.05 — 3 orders of magnitude below nominal significance, implausibly noise. Borderline hits (1e-3 ≤ padj < 0.05) stay in NC calibration on the principle that true noise floors include noise-adjacent fluctuations. Per-cluster granularity (up and down clusters at same exp+tp evaluated separately).

**Applied outcome:**
- Weissberg d11 **up** cluster → excluded from `(cyanorak_role, table_scope)` and `(kegg, table_scope)` NC calibration groups.
- Weissberg d11 **down** cluster → retained (no significant key-pathway hits).
- Steglich high-light 45m → retained in `(cyanorak_role, organism)` and `(kegg, organism)` NC calibration (hits above 1e-3 threshold).

Resulting NC calibration group sizes: `(*, table_scope)` = 7 (was 8); `(*, organism)` = 2 (unchanged). The `organism`-bg groups are flagged weakly-calibrated in caveats (2-cluster SD is fragile).

**Affects:** Task 10 (`scripts/05_compute_scores.py`) — `disqualified_nc_clusters()` helper + `compute_calibration(scores, nc_exclusions=...)` parameter wired per plan. `results/nc_calibration_exclusions.csv` artifact. LOO-R also applies the exclusions (disqualification is a property of the NC cluster, not the signature).

---

---

## D5 — Watchlist curation: exclude `kegg.pathway:ko05152` (Tuberculosis) from `signature_dropped.csv`

**Decided:** 2026-04-21, Step 3 decide
**Scope affected:** `signature_dropped.csv` only. Reference signature is never touched by this curation.

The `below_threshold_notable` tie-breaker rule (|max_signed_score| ≥ 3.0 OR `n_up + n_down ≥ 2`) in `signature.py` surfaced 6 dropped-notable terms from the Step 3 run. One of them — `kegg.pathway:ko05152 Tuberculosis` — is a KEGG disease-pathway map built around *Mycobacterium tuberculosis* host–pathogen biology. Its 2 down-direction hits in MED4 R clusters (Tolonen 24h + 48h, `|s|` ≤ 2.48) arise from generic-enzyme cross-annotations, not from biology relevant to *Prochlorococcus* N-limitation. Keeping it on the "terms to look at" shortlist clutters downstream review.

**Implementation.** `scripts/04_derive_signature.py` carries an explicit `WATCHLIST_EXCLUDE_TERMS: set[str]` constant; terms listed there are filtered out of `signature_dropped.csv` after derivation (not from the signature). Re-runs reproduce the exclusion; step3.log records a line `Watchlist curation (D5): excluded N term(s) from signature_dropped.csv: [...]` so the curation is auditable. The constant currently contains only `kegg.pathway:ko05152`; future entries require an amendment to this log.

**Retained watchlist terms (5):** cyanorak `D.4 Chaperones` ↓, `L.3 Protein folding` ↓, `R.2 Conserved hypothetical` ↑; kegg `ko00061 Fatty acid biosynthesis` ↓, `ko01212 Fatty acid metabolism` ↓. All biologically plausible under N-limit (chaperone/proteostasis throttling; fatty-acid metabolism reorganization). To be tracked through Step 4 inspection — if any show T-cluster signal above the noise floor, document as post-signature observations, not pre-registered tests.

**Affects:** `signature_dropped.csv` = 5 rows (down from 6). Runner prints `Watchlist curation (D5): excluded 1 term(s)` at each run.

---

## D6 — Signature stability handling: option (A) audit-only for within-ontology redundancy; ko00910 flagged for Step 4 LOO-R

**Decided:** 2026-04-21, Step 3 decide (after `explore_step3_redundancy_audit.py` + `explore_step3_single_exp_dominance.py`)
**Alternatives considered:** option (B) post-filter in Step 3, option (C) de-weight during Step 4 scoring

### Part 1 — Within-ontology redundancy

**Findings from `exploration/qc/step3_signature_redundancy_*.csv`:**

- **cyanorak_role:** 0 flagged pairs. Max pairwise Jaccard = 0.076 (D.1 ∩ E.4, diluted by D.1's 213-gene catch-all). L1 hierarchical disjointness holds across the 7 signature terms; J.1 ATP synthase does NOT overlap J.7/J.8/K.2 within cyanorak.
- **kegg:** 1 flagged pair — strict subset `ko00710 Calvin cycle (16 genes) ⊂ ko01200 Carbon metabolism (58 genes)` — every Calvin gene KEGG annotates is also in the umbrella Carbon-metabolism map. Jaccard = 0.276. 1 known-biology soft overlap — `ko00190 Oxidative phosphorylation ∩ ko00195 Photosynthesis = 9 atp genes` (the atpA-I operon, per Step 2 Q3). Jaccard = 0.120, below the 0.5 flag threshold but substantively meaningful because coordinate atp-down double-counts as evidence.

**Option (A) chosen — audit-only.** Both kegg overlap findings reflect KEGG's annotation structure (hierarchical umbrella; shared core-metabolic operon), not analytical error. Dropping either member of an overlap pair loses distinct biology:
- Dropping `ko00710`: loses Calvin-specific evidence. `ko01200` is broader (TCA + glycolysis + gluconeogenesis + Calvin) and the down-direction in it could be driven by any subset.
- Dropping `ko01200`: loses the TCA/glycolysis/gluconeogenesis genes not contained in `ko00710` or `ko00190`.
- Dropping `ko00190` or `ko00195`: both speak to distinct N-limit biology (energy throttling vs light-reaction throttling) despite sharing the atp operon.

Option (C) (per-pathway de-weighting in Step 4) was rejected because the Jaccard-based scaling rule would also dampen biologically independent pairs and inject a methodology-level tuning parameter for marginal gain. Option (B) (post-filter now) was rejected because no single drop is defensible.

**What "audit-only" means operationally:**

1. No change to `reference_signature.csv` at Step 3.
2. Task 10 (`05_compute_scores.py`) adds an **M4 redundancy sensitivity sub-check** alongside LOO-pathway and LOO-R:
   - Compute kegg `score_A(T_cluster)` with the full signature AND with `ko00710` removed (largest within-kegg redundancy).
   - Compute kegg `score_A(T_cluster)` with the full signature AND with BOTH `ko00710` and `ko00195` removed (aggressive check — drops both of the atp-operon's home maps except ko00190).
   - Record per-T-cluster deltas in `results/kegg_redundancy_sensitivity.csv`.
   - Flag any T-cluster whose classification flips (significant ↔ not-significant) under these subsets.
3. Task 12 (Step 5 caveats.md): add caveat **C6** documenting both overlaps (Calvin⊂C-metab subset; atp-operon ko00190∩ko00195) with the sensitivity-check results summarized.
4. **Cross-ontology atp reference:** Step 2 Q3 noted `cyanorak.role:J.1 ATP synthase (10 atp genes) ⊂ kegg.pathway:ko00190 ⊂ kegg.pathway:ko00195`. This is NOT a scoring concern because Step 4 `score_A` aggregates per-ontology — cyanorak and kegg scores are independent. But the cross-ontology agreement test (§5 Step 4 M3) should flag "atp operon agreement" as distinct from independent cross-ontology confirmation of oxphos/photosynthesis. Note to Task 10: cross-ontology agreement comparison must disclose the atp-operon shared-evidence nature when coordinating J.1↓ with ko00190↓ and ko00195↓.

### Part 2 — Single-R-experiment dominance

**Finding from `exploration/qc/step3_single_exp_dominated.csv`:**

- `kegg.pathway:ko00910 Nitrogen metabolism (up)` — all 4 supporting clusters from Tolonen 2006 microarray; Read 2017 contributes 0 up-direction significant clusters on this term. `share_max_exp = 1.00`. Confirms the Step 2 preview: Read's 6-gene kegg N-metab pathway does not survive padj<0.05 despite Read's 28-gene cyanorak E.4 N-metab pathway reaching significance. Read's detection floor on the narrower kegg pathway is plausibly the explanation (fewer genes in the term, smaller effective statistical power).
- All 12 other signature terms have 2-experiment support (both Tolonen and Read contribute ≥1 same-direction significant cluster).

**No Step 3 action; flag forward to Step 4 LOO-R (M4).** Expected LOO-R behavior: removing Tolonen → ko00910 drops entirely from the kegg signature (0 R support); removing Read → ko00910 retained (4 Tolonen clusters remain, still ≥3).

Task 10 M4 LOO-R must record this expected drop and verify it; Task 12 caveats adds entry **C7** noting ko00910 is a Tolonen-only signature member, with the interpretation-layer note that cyanorak E.4 N-metab retains both-experiment support and therefore carries the cross-ontology N-metab agreement.

**Affects:** Task 10 — adds `results/kegg_redundancy_sensitivity.csv` artifact + LOO-R ko00910 expected-drop check. Task 12 — caveats C6 (redundancy), C7 (ko00910 dominance) to be written at Step 5.

---

_Subsequent decisions (D7+) will be added by Task 9 (pre-registration of T outcomes) and Task 11 (Step 4 decide)._
