# Results Manifest

All files produced by Step 4 scoring + stability. See [decisions.md](../decisions.md) D7 (pre-registration), D8 (scoring methodology), D4 (NC calibration exclusion), D6 (signature stability handling).

## Scoring (D8 methodology — max-abs + 1.301 threshold + per-(exp×tp×ontology) aggregation)

| File | Rows | Produced by | Description |
|------|------|-------------|-------------|
| `scores_all.csv` | 72 | `05_compute_scores.py::score_all_exp_tp` | Per (experiment × timepoint × ontology × bg) signature-score rows. Columns: category (coculture/axenic_dying/axenic_dead/empty), omics (RNA/Prot/Microarray), timepoint, ontology, experiment_id, class, background_used, organism, score_up, score_down, final_score, n_up_sig, n_down_sig, n_up_total, n_down_total, n_total_sig, n_total_terms. Covers T + R + PC + NC + CTX classes (36 unique exp×tp × 2 ontologies). |
| `score_summary.csv` | 28 | `05_compute_scores.py::main` | T-only subset of scores_all.csv with classification label appended. 14 T (exp × tp) × 2 ontologies. Note: classification labels are all `insufficient_nc` under D8 scoring because NC noise variance collapses (see D8-e); column retained for audit only, not used for interpretation. |
| `category_mean_scores.csv` | 36 | `05_compute_scores.py::compute_category_means` | Category × omics × ontology × direction means. `direction ∈ {up, down, final}`, `omics ∈ {RNA, Prot, ALL}`. `omics=ALL` rows are the D7 formal-check aggregation (RNA+Prot mixed); `omics=RNA` and `omics=Prot` are the primary interpretation per D8-f. |
| `preregistration_check.csv` | 8 | `05_compute_scores.py::check_preregistration` | D7 P1-P4 ordering claims + T1-T4 threshold claims evaluated on `omics=ALL` category means. Outcomes: P1 FALSE (cyanorak×up: B=1.33 < A=2.37, H1 surface reading holds); P2, P3, P4 TRUE; T1-T4 all vacuously TRUE (threshold≈0 per D8-e classification degeneracy). |

## Calibration + NC exclusions

| File | Rows | Produced by | Description |
|------|------|-------------|-------------|
| `nc_calibration.csv` | 12 | `05_compute_scores.py::compute_calibration` | NC noise + PC peak calibration per (ontology × background_used × value_col). `value_col ∈ {final_score, score_up, score_down}`. Under D8, most groups have nc_std≈0 because the signature almost never triggers in non-N-limit NC biology (signature-specificity signal). cyanorak table_scope nc_std: final=0.02, up=0.02, down=0.01; kegg table_scope nc_std all 0.0. |
| `nc_calibration_exclusions.csv` | 2 | `05_compute_scores.py::disqualified_nc_exp_tp` | Per-(experiment × timepoint × ontology × background) NC exclusions applied per D4 (padj<1e-3 on key-pathway anchors). Weissberg 2025 coculture RNA day 11 UP cluster excluded from (cyanorak_role × table_scope) and (kegg × table_scope) — real coculture-induced N-scavenging biology, not noise. |

## D7 pre-registration evaluation

See `preregistration_check.csv` above. Summary on omics-ALL level:

- P1 cyanorak_role × up (B > A): **FALSE** (B=1.33, A=2.37). Coculture engages N-scavenging more than axenic-dying on cyanorak up.
- P2 cyanorak_role × down (B > A): **TRUE** (B=2.69, A=0.30; ~9× ratio).
- P3 kegg × up (B > A): **TRUE** marginal (B=1.57, A=1.47). Fragile — kegg up is single-term ko00910, Tolonen-only per D6.
- P4 kegg × down (B > A): **TRUE** (B=2.38, A=0.65; ~3.7× ratio).
- T1-T4 threshold claims: all TRUE but vacuously (thresholds ≈ 0 under D8-e).

**Per-omics breakdown (primary scientific interpretation per D8-f)** in `category_mean_scores.csv` reveals asymmetric pattern: coculture's up-signature lives in Prot (up=4.74 cyanorak), not RNA (up=0.00); axenic-dying's down-signature lives in RNA (down=4.99 cyanorak), only partially in Prot (down=0.38 — only ribosome at protein level). Task 11 interpretation uses this, not the omics-mixed summary.

## Stability (D6 redundancy + LOO-signature + LOO-R)

| File | Rows | Produced by | Description |
|------|------|-------------|-------------|
| `kegg_redundancy_sensitivity.csv` | 14 | `05_compute_scores.py::kegg_redundancy_sensitivity` | Per T × kegg (exp × tp): final_score under 3 kegg signature variants (full / no_ko00710 / no_ko00710_ko00195) + relative change + material-shift flag (>50%). Removing ko00710 alone: 0 material shifts. Removing ko00710 + ko00195: 4 of 14 shift >50% (score magnitudes ~2× due to denominator effect when signature shrinks from 6→4 terms). Subset redundancy per D6 Part 1 is robust. |
| `loo_signature.csv` | 182 | `05_compute_scores.py::main` | Per T × ontology × signature-term-removed: final_score before/after, flag_sign_flip (sign changes), flag_large_drop (\|new\| < 0.5×\|orig\|). 6 sign flips (3%) + 14 large drops (8%) out of 182 — modest sensitivity, no single-term dominance. Flags per D8-h. |
| `loo_r_experiments.csv` | 56 | `05_compute_scores.py::main` | Per T (exp × tp × ontology × bg) × R-experiment-removed: orig/new final_score, rel_change, flag_sign_flip_raw, flag_large_change_raw (\|rel_change\|>50%), orig/new classification + flip, new_signature_size. D8-g raw-score flags replace classification-flip as primary stability signal because D8-e makes classification vacuous. Dropping Read: 0 raw flags (Read is redundant). Dropping Tolonen: 8 raw large changes + 1 sign flip; cyanorak signature shrinks 7→5, kegg shrinks 6→3. Category-mean ordering preserved under both. |

## Interpretation for Task 11 decide

The asymmetric-outcome framing (B > A on down-direction, A > B on cyanorak up-direction, B > A on kegg up-direction marginal) is the Step 4 finding. Per-omics breakdown (cyanorak Prot vs RNA; kegg Prot vs RNA) shows the mechanism: coculture engages N-scavenging at the protein level without shutting down machinery; axenic-dying engages classical stress shutdown strongly at the transcript level but only partially at the protein level (proteins haven't turned over in 14 days). See [exploration/notebook.md](../exploration/notebook.md) 2026-04-22 entry for full data-vs-interpretation audit and the prose summary.
