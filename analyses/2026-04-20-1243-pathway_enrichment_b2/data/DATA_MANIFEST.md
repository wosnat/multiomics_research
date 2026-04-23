# Data Manifest

All files produced by extraction scripts from the multiomics KG.
Run scripts from `multiomics_research` root with `uv run`.

## Experiment landscape

| File | Rows | Genes | Timepoints | Produced by | Description |
|------|------|-------|------------|-------------|-------------|
| `experiments_classified.csv` | 37 (15 unique exps) | n/a | varies (1–6 per exp) | `01_select_experiments.py` | T/R/PC/NC/CTX-classified experiments joined with `list_experiments` metadata + `list_publications` first_author / title attribution. 4 T (Weissberg), 2 R (Tolonen 2006, Read 2017), 2 PC (Tolonen alt-N), 5 NC (Aharonovich, Weissberg, Moreno-Cabezuelo ×2, Steglich), 2 CTX (Tolonen MIT9313, Domínguez-Martín SS120). HOT1A3 CTX1 dropped 2026-04-20 during Step 1b ontology-compatibility review — see notebook Step 1a redo entry. |
| `landscape_Prochlorococcus_MED4.csv` | 62 | n/a | n/a | `02_ontology_landscape.py` | `ontology_landscape(organism="MED4", experiment_ids=[13 MED4 exps])` — per-(ontology, level) stats: relevance_rank, n_terms_with_genes, genome_coverage, median_genes_per_term. `example_terms` column dropped by `to_dataframe` (nested list). |
| `landscape_Prochlorococcus_MIT9313.csv` | 64 | n/a | n/a | `02_ontology_landscape.py` | As above for MIT9313 CTX2 experiment. |
| `landscape_Prochlorococcus_marinus_subsp._marinus_CCMP1375_SS120.csv` | 62 | n/a | n/a | `02_ontology_landscape.py` | As above for SS120 CTX3 experiment. |
| `nitrogen_ontology_search.csv` | 56 | n/a | n/a | `02_ontology_landscape.py` | `search_ontology("nitrogen", ontology=X)` across 7 ontologies (cyanorak_role, tigr_role, go_bp, go_mf, kegg, ec, pfam; BRITE excluded — needs tree filter). Fields: id, name, score, level, ontology. |

## Enrichment outputs

| File | Rows | Genes | Timepoints | Produced by | Description |
|------|------|-------|------------|-------------|-------------|
| `enrichment_all.csv` | 11,239 | n/a | per-cluster | `03_run_enrichment.py` | Fisher ORA results concatenated across 8 (organism, ontology, background) calls. 70 unique clusters. Columns include cluster, term_id, term_name, gene_ratio, bg_ratio, fold_enrichment, pvalue, p_adjust, signed_score, timepoint, direction, experiment_id, omics_type, table_scope, organism, ontology, background_used. 225 rows with p_adjust<0.05. Requires upstream fix in multiomics_explorer `de_enrichment_inputs` (merged 2026-04-20) — prior runs had bg collapsed to foreground. |
| `enrichment_results.pkl` | 8 dict entries | n/a | n/a | `03_run_enrichment.py` | Dict keyed by (organism, ontology, background_used) → EnrichmentResult instance. Preserves `.explain()` / `.overlap_genes()` / `.background_genes()` accessors. Verified via stage-1 + stage-2 pickle round-trip. |

## Exploration artifacts (step 2 QC)

| File | Rows | Produced by | Description |
|------|------|-------------|-------------|
| `exploration/qc/step2_rkey_matrix.csv` | 187 | `explore_step2_rkey_agreement.py` | R clusters × 11 key pathways, one row per combination. Columns include class, first_author, experiment_short, timepoint, direction, term_id, expected_direction, signed_score, p_adjust, is_significant, agrees_with_expected. |
| `exploration/qc/step2_rkey_nc_enrichment.csv` | 98 | `explore_step2_rkey_agreement.py` | Same structure for NC clusters × key pathways. Used for noise-floor inspection. |
| `exploration/qc/step2_rkey_summary_by_term.csv` | 8 | `explore_step2_rkey_agreement.py` | Per-term summary: n_significant, n_agree, mean_abs_signed_score, min_padj across 17 R clusters (Tolonen 12 + Read 5 after temporal filter; Read 3h up is the 1 disagreement). 3 key pathways (A.3, ko00250, ko00260) had 0 significant R hits. |
| `exploration/qc/step2_rkey_summary_by_expTp.csv` | 17 | `explore_step2_rkey_agreement.py` | Per-(experiment, timepoint, direction) count: n_keypath_tested=11, n_significant, n_agree_expected, n_disagree_expected. |
| `exploration/qc/step2_rkey_disagreements.csv` | 1 | `explore_step2_rkey_agreement.py` | Significant R cluster × key pathway rows where direction disagrees with expected. Single entry: Read 2017 3h up × cyanorak.role:J.2 CO2 fixation, padj=0.045 (borderline; early-TP transient flip per B1 observation and spec §5 Step 3 temporal filter rationale). |
| `exploration/qc/step2_cluster_totals.csv` | 140 | `explore_step2_cluster_totals.py` | Total significant hits per (cluster × ontology × bg) across full ontology — disambiguates rkey-summary's "N significant" (which was scoped to 11 key pathways only). Columns: organism, ontology, background_used, experiment_id, timepoint, direction, n_tested, n_sig, experiment_short. |
| `exploration/qc/step2_r_top_pathways.csv` | 26 | `explore_step2_r_top_pathways.py` | Main enriched pathways in signature-eligible R clusters (`timepoint_hours > 3`, decision #1) — per (ontology, term, direction): n_clusters (padj<0.05), mean_abs_signed, min_padj, contributing experiments. Preview of Step 3 signature eligibility. |
| `exploration/qc/step2_heatmap_cyanorak_role.png` | — | `explore_step2_heatmap_v2.py` | Two-panel QC heatmap for cyanorak_role: top = non-T (R\|PC\|CTX\|NC, class dividers), bottom = T (axenic\|coculture divider). Display cap ±5 + saturation stars `*`/`**`/`***` at \|s\|≥5/10/20. Rows = key panel (bold) + discovered-strong (n_sig≥3 in R at hours>3). Authors truncated 6-char, timepoints abbreviated. Inline 3-line legend. |
| `exploration/qc/step2_heatmap_kegg.png` | — | `explore_step2_heatmap_v2.py` | As above for kegg. |
| `exploration/qc/step2_score_distribution.png` | — | `explore_step2_score_distribution.py` | 4-panel analysis of `signed_score` distribution over all 11,239 enrichment rows: histogram, per-class |s| distribution, CDF, top-50 by |s|. Informed decision to keep spec's ±10 scoring cap (saturation regime beyond ~|s|=10) but use ±5 display cap for heatmaps. |

## Signature

| File | Rows | Genes | Timepoints | Produced by | Description |
|------|------|-------|------------|-------------|-------------|
| `reference_signature.csv` | 13 | n/a | n/a | `04_derive_signature.py` | MED4 reference N-limitation signature. Per-ontology derivation on 12 signature-eligible R clusters (Tolonen 6h/12h/24h/48h × 2dir = 8 + Read 12h/24h × 2dir = 4, after `timepoint_hours > 3.0` filter per D1). Core rule `n_clusters_supporting ≥ 3` applied to both ontologies (no fallback). cyanorak_role: 7 terms (5 down: J.1 ATP synthase, J.2 CO2 fixation, J.7 PSI, J.8 PSII, K.2 Ribosome; 2 up: D.1 Adaptation, E.4 N-metabolism). kegg: 6 terms (5 down: ko00190 Oxphos, ko00195 Photosynthesis, ko00710 Calvin, ko01200 C-metab, ko03010 Ribosome; 1 up: ko00910 N-metab). Columns: term_id, term_name, direction, n_clusters_supporting, n_up, n_down, contributing_clusters, per_experiment_breakdown (JSON), ontology, rule_applied. |
| `signature_dropped.csv` | 5 | n/a | n/a | `04_derive_signature.py` | Terms below the `≥3` threshold but notable (`|signed_score| ≥ 3.0` or `n_up+n_down ≥ 2`), curated via `WATCHLIST_EXCLUDE_TERMS` (D5: ko05152 Tuberculosis excluded as KEGG disease map, biologically N/A). cyanorak: D.4 Chaperones down, L.3 Protein folding down, R.2 Conserved hypothetical up. kegg: ko00061 Fatty acid biosynth down, ko01212 Fatty acid metab down. All `drop_reason=below_threshold_notable` (no bidirectional drops). Columns: term_id, term_name, n_up, n_down, drop_reason, max_signed_score, contributing_clusters_up, contributing_clusters_down, per_experiment_breakdown, ontology. |

## Exploration artifacts (step 3 QC)

| File | Rows | Produced by | Description |
|------|------|-------------|-------------|
| `exploration/qc/step3_per_exp_support.csv` | 26 | `explore_step3_single_exp_dominance.py` | Long-form per-(signature term × R experiment) support. Columns: ontology, term_id, term_name, direction, experiment_id, n_up, n_down, n_in_direction. 13 signature terms × 2 R experiments = 26 rows. |
| `exploration/qc/step3_single_exp_dominated.csv` | 1 | `explore_step3_single_exp_dominance.py` | Signature terms where only 1 R experiment contributes same-direction significant clusters. Single flag: `kegg.pathway:ko00910 Nitrogen metabolism` (up) — 4 Tolonen 2006 clusters, 0 Read 2017. `share_max_exp=1.00`. Forward-flagged to Task 10 LOO-R (M4). |
| `exploration/qc/step3_signature_redundancy_cyanorak_role.csv` | 21 | `explore_step3_redundancy_audit.py` | Pairwise Jaccard on MED4 gene membership, 7×6/2=21 signature pathway pairs. 0 flagged (all Jaccard < 0.08, no subsets). Hierarchical L1 disjointness confirmed. |
| `exploration/qc/step3_signature_redundancy_kegg.csv` | 15 | `explore_step3_redundancy_audit.py` | Pairwise Jaccard on MED4 gene membership, 6×5/2=15 signature pathway pairs. 1 flag: `ko00710 Calvin ⊂ ko01200 C-metab` (16/16, Jaccard 0.276). Soft note: `ko00190 ∩ ko00195 = 9 atp genes` (Jaccard 0.120, below flag but biologically meaningful — per Step 2 Q3 atpA-I operon cross-reference). |
| `exploration/qc/step3_signature_redundancy_cyanorak_role.png` | — | `explore_step3_redundancy_audit.py` | Jaccard heatmap, OrRd colormap, symmetric matrix. |
| `exploration/qc/step3_signature_redundancy_kegg.png` | — | `explore_step3_redundancy_audit.py` | Jaccard heatmap. The Calvin/C-metab strict-subset cell is the darkest off-diagonal. |

## Exploration artifacts (step 4 QC)

| File | Rows | Produced by | Description |
|------|------|-------------|-------------|
| `exploration/qc/step4_t_contribution_decomposition.csv` | 182 | `explore_step4_contribution_decomposition.py` | Long-form per-(T cluster × signature term): up_signed_score, down_signed_score, up_contribution, down_contribution, max_abs_contribution, contributing_cluster (UP/DOWN/-). Source data for Fig B. 14 T rows × 13 signature terms = 182. |
| `exploration/qc/step4_t_top_contributors.csv` | 42 | `explore_step4_contribution_decomposition.py` | Top-3 contributors per (T × ontology) by \|max_abs_contribution\|, zeros excluded. 14 T rows × 3 = up to 42 rows. |
| `exploration/qc/step4_cross_ontology_agreement.csv` | 14 | `explore_step4_contribution_decomposition.py` | Per T (exp × tp × bg), cyanorak vs kegg agreement on final_score, score_up, score_down. Categories: same_sign, opposite_sign, both_weak, missing. 0 opposite-sign disagreements on all three. |
| `exploration/qc/step4_figA_trajectory.png` | — | `explore_step4_figures.py` | 6-panel per-TP trajectory: 2 ontologies × (RNA + Prot), 3 categories overlaid. Up-direction solid, down-direction dashed. Dotted grey line at display cap ±5. |
| `exploration/qc/step4_figB_contributions.png` | — | `explore_step4_figures.py` | 13-row × 28-column heatmap of max-abs contributions. Rows: signature terms grouped by ontology + expected direction. Columns: T clusters grouped by category \| omics \| tp \| ontology. Diverging blue-red, display cap ±5, saturation stars at \|5\| and \|10\|. Anti-signature (negative) annotated with leading '-'. Category + ontology dividers as solid black lines; direction dividers as dotted grey. |
| `exploration/qc/step4_figC_category_means.png` | — | `explore_step4_figures.py` | 2-panel (cyanorak, kegg) bar chart. Per category on x-axis, 4 bars: RNA ↑ (blue solid), RNA ↓ (blue hatched), Prot ↑ (orange solid), Prot ↓ (orange hatched). Summarizes asymmetric-outcome pattern. |
