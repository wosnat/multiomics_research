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

---

## D7 — Step 4 pre-registration: 3-category framework for Weissberg T conditions

**Written:** 2026-04-21, Step 4 do (Task 9, pre-registration) — BEFORE any Step 4 scoring runs.
**In tension with:** hypotheses.md H1 (preliminary Step 2 QC observation "coculture > axenic"). D7 pre-registers the opposite ordering. Data arbitrates at Step 4 decide (Task 11) → D8.
**Uncommitted note:** this entry is saved to `decisions.md` but will NOT be committed until Task 10 Step 5 (per plan), bundled with the scoring-script commit so pre-registration lands in the same commit as the code that tests it.

### Framework — 3 biologically-motivated categories

The 28 Weissberg T clusters split into 3 categories based on cell physiological state, not omics type or timepoint directly:

| cat | biology | T clusters | n |
|---|---|---|---:|
| **A. Coculture** | MED4 rescued by *Alteromonas* partner; cells surviving through 90 days; N-stress attenuated by cross-feeding | T2 RNA coculture (d18/d31/d60/d89/d60+89) × 2 dir + T4 Prot coculture (same 5 TPs) × 2 dir | 20 |
| **B. Axenic alive-but-dying** | Acute N-stress just before death; cells actively engaging canonical N-response program (the R-cluster analog in the T set) | T1 RNA axenic d14 × 2 dir + T3 Prot axenic d14 × 2 dir | 4 |
| **C. Axenic decomposing** | Post-death cells; no active gene regulation, only protein/RNA decay dynamics; signature signal may reflect passive decay, not response | T3 Prot axenic d31 × 2 dir + T3 Prot axenic d89 × 2 dir | 4 |

Rationale for B being the "strongest N-response" category: per biology context memory, MED4 dies axenically under N-limit in ~14 days. Day 14 is the acute stress / peri-mortem window — cells have engaged the full N-limit response program for days and are still transcriptionally/translationally active. Tolonen 48h (R cluster) is the closest axenic-kinetics analog on the R side: "acute stress, still responding." The prediction is that axenic d14 clusters are essentially T-set analogs of the Tolonen 48h R clusters that defined the signature.

### Operationalization — per-direction, per-ontology scoring

**Note added 2026-04-22 (post-walkthrough amendment):** the operational formula below describes the pre-registration *as a claim*. The actual scoring mechanics for evaluating these claims are defined in decisions.md D8 (scoring methodology), which introduces (i) per-term max-abs cluster combination, (ii) significance threshold at |contrib| ≥ 1.301, (iii) omics tagging, and (iv) per-(exp × tp × ontology) aggregation. The D7 formal check (P1-P4, T1-T4) uses category means computed with omics='ALL' — i.e., RNA and Prot rows averaged into one category mean per (category × ontology × direction). This omics-mixed aggregation is the D7-specific summary only; the **primary scientific interpretation at Task 11 relies on the per-omics breakdown** (RNA and Prot as separate columns), because RNA-protein correlation is known to be low and mixing them in a mean obscures mechanism.

### Original pre-registration formula (frozen as written 2026-04-21)

For each T cluster × ontology × direction (up or down), compute:

```
dir_score(cluster, ontology, d) = mean over signature terms with expected_direction=d of
    (+1 if expected_direction=up else -1) * signed_score(cluster, term, ontology)
```

Positive `dir_score` = cluster matches the signature's expectation in that direction (up-expected pathways up, or down-expected pathways down). Negative = opposes. Magnitude = strength of canonical engagement.

**Direction × ontology coverage per signature (from reference_signature.csv):**

| ontology | up-expected terms | down-expected terms |
|---|---|---|
| cyanorak_role | D.1 Adaptation, E.4 N-metab (n=2) | J.1 ATP synthase, J.2 CO2 fix, J.7 PSI, J.8 PSII, K.2 Ribosome (n=5) |
| kegg | **ko00910 N-metab (n=1)** — single term, Tolonen-only per D6 Part 2 (fragile evidence base) | ko00190 Oxphos, ko00195 Photosynth, ko00710 Calvin, ko01200 C-metab, ko03010 Ribosome (n=5) |

⚠ **kegg up-direction fragility.** The kegg up-expected direction has only 1 signature term (ko00910), which is itself single-R-experiment dominated (Tolonen-only per D6 Part 2). `kegg_up_score` is effectively the signed contribution of ko00910 alone. The cross-ontology agreement check (M3) with `cyanorak_up_score` (which has 2 terms and is supported by both R experiments) is the more robust test of the up-direction claim.

**Category-level aggregation:**

```
cat_dir_score(category, ontology, direction) = mean over clusters in category of dir_score
```

### Pre-registered predictions (4 ordering claims + 4 threshold claims + 2 no-prediction claims)

**Primary ordering claims (test H1 reframing):**

| # | claim | interpretation if TRUE | interpretation if FALSE |
|---|---|---|---|
| P1 | `cat_dir_score(B, cyanorak, up) > cat_dir_score(A, cyanorak, up)` | axenic-alive engages N-scavenging (E.4 + D.1) more than coculture | coculture N-scavenging engaged more — H1's surface observation reinstated on cyanorak up |
| P2 | `cat_dir_score(B, cyanorak, down) > cat_dir_score(A, cyanorak, down)` | axenic-alive shuts down PS/ribosome more than coculture (stronger canonical response) | coculture shuts down more — rescue reframing weakened |
| P3 | `cat_dir_score(B, kegg, up) > cat_dir_score(A, kegg, up)` | single-term ko00910 — fragile, robustness depends on LOO-R. Cross-ontology agreement with P1 required. | kegg ko00910 fragile; interpret only in conjunction with P1 |
| P4 | `cat_dir_score(B, kegg, down) > cat_dir_score(A, kegg, down)` | axenic-alive throttles oxphos/C-metab/photosynth/ribosome more than coculture | coculture throttles more — reframing falsified on kegg down |

**Threshold claims (coculture category A still engages the signature above noise):**

| # | claim | threshold source |
|---|---|---|
| T1 | `cat_dir_score(A, cyanorak, up) > nc_mean_{cyanorak, table_scope} + 2σ_{cyanorak, table_scope}` | spec §5 Step 4 M3, with D4 NC exclusions applied |
| T2 | `cat_dir_score(A, cyanorak, down) > nc_mean + 2σ` | " |
| T3 | `cat_dir_score(A, kegg, up) > nc_mean + 2σ` | " (fragile — see P3 note) |
| T4 | `cat_dir_score(A, kegg, down) > nc_mean + 2σ` | " |

**No-prediction claims (category C):**

- `cat_dir_score(C, *, *)` — NO pre-registered prediction on any direction or ontology. Rationale: post-mortem cells have no active gene regulation; signature signal reflects protein/RNA decay kinetics, not canonical response engagement. Whatever C shows is a Step 4 explore-phase observation to document, not a pre-registered test.

### Falsification summary

- **D7 reframing (B > A) is CONFIRMED** if all 4 of P1–P4 hold (strong), OR if P1–P2 hold and P4 holds (kegg up P3 may fail due to ko00910 fragility; cyanorak cross-ontology test still required).
- **D7 reframing is FALSIFIED** if P1 AND P2 AND P4 all fail (A > B on cyanorak up+down and kegg down). In that case H1's surface observation is reinstated and promoted from preliminary to confirmed in D8.
- **Partial outcomes** (e.g., P2 holds but P1 fails — axenic throttles more but doesn't scavenge more) require biological reinterpretation: e.g., axenic shows passive-decay asymmetry rather than active-response symmetry. Document at Task 11 (Step 4 decide → D8).
- **Threshold claims T1–T4** test a secondary point: coculture is not FULLY rescued to baseline. If T1–T4 fail, coculture cells show no signature engagement above noise — implying full rescue, not attenuated response. That's a third reframing (not pre-registered, documented if it emerges).

### Stability requirements — LOO-R applies to ordering claims

Per D6 Part 2, Task 10 M4 LOO-R already recomputes signature (and hence `dir_score`) with each R experiment removed. D7 claims P1–P4 must also survive LOO-R:
- If P3 flips when Tolonen removed (expected, since ko00910 drops from signature entirely), flag explicitly — claim is only meaningful with Tolonen in.
- If P1, P2, P4 flip on ANY LOO-R removal, claim is fragile and must be caveated at Task 12 (C8).

### Affects

- **Task 10 (`05_compute_scores.py`):** must compute `dir_score` per (cluster × ontology × direction) and category-means per (category × ontology × direction) as new artifacts `results/t_cluster_dir_scores.csv` and `results/category_mean_dir_scores.csv`. Category column `category ∈ {A, B, C}` derived from (experiment_id, timepoint) per the membership table above.
- **Task 11 (Step 4 decide → D8):** must evaluate P1–P4 + T1–T4 outcomes; promote H1 to CONFIRMED/FALSIFIED; document LOO-R behavior per D7 stability requirements.
- **Task 12 caveats:** C8 for D7 outcome narrative + LOO-R behavior; C6 (redundancy) and C7 (ko00910 dominance) also get refreshed with the Step 4 observed values.

### Biology memo (for Task 11 interpretation)

Three possible biological stories from D7 outcomes:

1. **Reframing confirmed (D7 predictions hold):** Coculture buffers MED4 from acute N-stress, allowing survival by reducing the need for full canonical N-response engagement. Axenic alive-but-dying cells show the strongest canonical response because they have the most acute stress. The paper's mechanistic claim: "coculture rescues survival by attenuating the N-response itself, not by boosting it."

2. **H1 reinstated (D7 falsified):** Coculture engages N-response MORE than axenic — the observation from Step 2 QC was not heatmap artifact. The mechanism: coculture enables *sustained* engagement of N-response (axenic cells shut down generally before fully engaging specific circuitry; coculture keeps cells alive long enough to mount a full response). Paper claim: "coculture rescues survival by enabling sustained N-response engagement."

3. **Asymmetric outcome (partial falsification):** Some directions/ontologies pass, others fail. Most likely scenario: P1/P2 (cyanorak up/down) hold because cyanorak has more robust support, but P4 (kegg down) fails because kegg response is pathway-specific in a way that differs between axenic and coculture. Interpretation is pathway-level biology and more granular narrative needed.

Story choice is NOT pre-registered. The data arbitrates at Step 4 decide.

---

---

## D8 — Step 4 scoring methodology (compute_signature_score primitive + max-abs + significance threshold + omics split)

**Decided:** 2026-04-22, Task 10 (Step 4 do) scoring walkthrough — before any Step 4 commit.
**Status:** methodology decision. Scoring code in `scripts/05_compute_scores.py` implements this; all Step 4 results (scores, category means, pre-reg check, LOO, kegg redundancy sensitivity) are produced under D8 mechanics.

The plain-English version: instead of scoring UP and DOWN gene-list clusters separately and averaging them into a category later, we combine the two clusters into a single per-(experiment × timepoint × ontology) signature-engagement row by picking whichever cluster has the stronger evidence for each signature term. We also zero out subthreshold noise (anything with padj ≥ 0.05) so weak whispers don't drag or inflate the score. RNA and proteomics are tagged and never mixed except at the D7 formal-check summary level. The biology-facing category labels replace the letter codes.

### D8-a — compute_signature_score primitive: max-abs cluster combination per term

For each (experiment × timepoint × ontology), we have two clusters: the UP gene set and the DOWN gene set. For each signature term:

- Compute the per-cluster contribution = `sign_ref × capped_signed_score` where `sign_ref = +1 if signature expects "up", -1 if signature expects "down"`, and `capped_signed_score = sign(signed_score) × min(|signed_score|, 10)`. Positive contribution = cluster matches the signature's expectation in this term's direction; negative = opposes.
- Take the max-absolute-value of the two cluster contributions for this term.
- Result: one evidence number per signature term, with its sign preserved.

Under this rule, up-expected terms naturally draw their evidence from the UP cluster (because UP genes enrich up-expected pathways, not down-expected ones) and down-expected terms from the DOWN cluster — no explicit direction-matching rule needed. Edge cases (both clusters enrich same term, rare anti-signature scenario) are handled by letting whichever carries more magnitude dominate; the dropped contribution is diagnostic for anti-signature anomaly reporting.

### D8-b — Significance threshold at |contrib| ≥ 1.301 (≈ −log₁₀ 0.05)

Subthreshold contributions are zeroed out. The threshold matches the signature derivation bar from Step 3 (padj < 0.05 was the cutoff for counting a cluster as supporting a term's direction). Rationale: scoring should use the same evidence bar as the signature derivation; subthreshold values introduce directional whispers that confuse interpretation without adding statistical justification. Implementation: in `scripts/05_compute_scores.py::thresholded_capped_contrib`, `return 0.0 if abs(capped) < SIGNIFICANCE_THRESHOLD else sign_ref * capped`.

### D8-c — Per-direction aggregation with omics tagging

After computing per-term max-abs contributions for a (exp × tp × ontology), split them by expected direction:

- `score_up = mean of max-abs contributions over up-expected signature terms`
- `score_down = mean of max-abs contributions over down-expected signature terms`
- `final_score = (n_up_total × score_up + n_down_total × score_down) / (n_up_total + n_down_total)` — n-weighted mean

Each row carries `omics` as a separate column (RNA / Prot / Microarray) derived from the experiment_id. At no point in the core aggregation are RNA and Prot rows averaged into a single number; omics-agnostic summaries appear only as the explicit D7 formal-check aggregation and are labeled `omics='ALL'` in `results/category_mean_scores.csv`.

### D8-d — Category names

T experiments are categorized by cell physiological state (per D7), with biology-facing labels replacing letter codes:

| label | cells | D7 category |
|---|---|---|
| `coculture` | MED4 rescued by Alteromonas partner, surviving 90 days | A |
| `axenic_dying` | Axenic MED4 at day 14, alive but about to die (acute N-stress window) | B |
| `axenic_dead` | Axenic MED4 day 31/89 proteomics, post-mortem decomposition | C |

Applied in `scripts/05_compute_scores.py::assign_category`.

### D8-e — Classification labels are non-informative; scores + ordering are the primary reporting

Spec §5 Step 4 M3 defines a classification scheme (`detectable / no_signal / pc_like / below_nc / insufficient_nc`) based on NC noise floor. Under D8's thresholded scoring, NC cluster scores collapse to 0 in most (ontology, background) groups (the signature is so specific that diverse non-N-limit biology doesn't trigger it at padj<0.05). `nc_std → 0`, so the classification function always returns `insufficient_nc`. This is itself informative — it says the signature is highly N-specific — but classification labels become non-informative as an analytical tool.

**Decision:** accept classification labels as non-informative under D8; treat `insufficient_nc` as an expected byproduct of signature specificity; rely on raw `final_score`, `score_up`, `score_down` + ordering comparisons (D7 P1-P4) as the primary reporting at Task 11 decide. Classification column is still written to `score_summary.csv` for audit, but downstream interpretation ignores it. Add **caveats.md C9** entry documenting this.

### D8-f — D7 formal-check uses omics-agnostic aggregation; interpretation uses per-omics breakdown

D7 was written with 4 ordering claims (P1-P4) and 4 threshold claims (T1-T4) at the `category × ontology × direction` level — no omics split. Evaluating these under D8 requires collapsing RNA and Prot into a single category mean, which violates the "can't mix RNA and Prot" principle.

**Resolution:** the D7 formal check (P1-P4, T1-T4 holds/doesn't hold) uses `omics='ALL'` category means — an explicit omics-mixed summary labeled as such in `results/category_mean_scores.csv` and `results/preregistration_check.csv`. At Task 11 decide, the **primary interpretation** uses the per-omics (RNA vs Prot) breakdown — which reveals mechanism-level asymmetry that the omics-mixed mean obscures. The D7 formal-check outcome is a summary signal; the per-omics pattern is the scientific conclusion. Both get recorded; Task 11's D-level decision records the asymmetric biology (not a binary "D7 confirmed/falsified").

### D8-g — LOO-R stability: raw-score change + sign flip (replaces classification flip as primary flag)

Spec §5 M3 specifies LOO-R as rederiving the signature with each R experiment removed, rescoring T clusters, and flagging classification flips. Under D8-e, classification is non-informative, so classification flips are trivially zero (everything stays `insufficient_nc`). New LOO-R stability flags on raw scores:

- `flag_sign_flip_raw` — `np.sign(orig_final_score) != np.sign(new_final_score)` (both non-zero)
- `flag_large_change_raw` — `|new - orig| / |orig| > 0.5` (when `|orig| > 1e-9`); if `|orig|` is ~0, flag if `|new| > 0.5`; if either is NaN, flag as large change

Threshold choice (0.5 = 50% relative change) is a judgment call: 30% feels too sensitive (common noise level); 100% too lax (misses material shifts). 50% is a middle-ground convention and can be revisited if Task 11 interpretation shows it's mis-tuned. Rationale for flagging: LOO-R stability failures at the raw-score level are the right signal for "is this finding robust to which R experiment anchors the signature?"

### D8-h — LOO-signature stability flags

Per-T × ontology × signature-term LOO: for each signature term removed, recompute `final_score` on the same (exp × tp × ontology) row and flag:

- `flag_sign_flip` — sign of new_score ≠ sign of orig_score
- `flag_large_drop` — `|new| < 0.5 × |orig|` (>50% magnitude drop)

Unchanged from previous scaffold; criteria stay at 50% for consistency with D8-g.

### D8-i — Whole-ontology QA columns (n_total_sig, n_total_terms)

Each scoring row includes two QA columns orthogonal to the signature:

- `n_total_sig` — count of ALL ontology terms (not just signature terms) significant at padj<0.05 in this (exp × tp × ontology), summed across UP and DOWN clusters.
- `n_total_terms` — total signature level term count for the ontology (69 for cyanorak L1, 97 for kegg L2).

Purpose: QA on the enrichment step. If `n_total_sig` is 0 for a cluster, the enrichment call returned nothing — something upstream is wrong (empty gene list, missing ontology, etc.). Doesn't affect scoring; lives alongside for auditability.

### Interactions with other decisions

- **D1 (temporal filter):** D8 inherits `hours > 3.0` via `signature.py::TIMEPOINT_HOURS_CUTOFF`; `rederive_signature_loo` in `05_compute_scores.py` uses `>` per D1 (verified via grep-check at Task 10 start).
- **D4 (NC calibration exclusion):** D4 was defined at per-cluster granularity; under D8's per-(exp × tp) aggregation, D4 exclusions collapse via `disqualified_nc_exp_tp`. Same excluded clusters (Weissberg coculture day 11 up) translate to the same per-(exp × tp) exclusions in new calibration.
- **D6 (redundancy sensitivity):** the D6 kegg signature variants (`full / no_ko00710 / no_ko00710_ko00195`) are applied through the D8 primitive — variants pass through `compute_signature_score` just like the main signature. Material-shift flag at >50% relative change (same threshold as D8-g for consistency).
- **D7 (pre-registration):** see D8-f. D7's formal claims are evaluated on omics-mixed category means; the biology story emerges from per-omics breakdown.

### Affects

- **Task 12 (Step 5 caveats):** add C9 documenting classification-label degeneracy under D8-e; amend C6 (redundancy) and C7 (ko00910 dominance) with D8-computed values.
- **Task 11 (Step 4 decide):** primary interpretation uses per-omics breakdown from `results/category_mean_scores.csv`; D7 omics='ALL' check is a summary, not the conclusion.

---

_Subsequent decisions (D9+) will be added by Task 11 (Step 4 decide)._
