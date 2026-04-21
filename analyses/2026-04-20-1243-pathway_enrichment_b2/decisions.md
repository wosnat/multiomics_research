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

_Subsequent decisions (D5+) will be added by Task 9 (pre-registration of T outcomes) and Task 11 (Step 4 decide)._
