# Results Manifest

All files produced by scoring, plotting, and annotation scripts.

## Score tables

| File | Rows | Produced by | Description |
|------|------|-------------|-------------|
| `signature_scores_core.csv` | 21 | `score_signature.py` | Core signature (198 genes) scored against all experiments: 14 Weissberg condition x timepoints + 7 reference timepoints. Columns: hit_rate_concordant, hit_rate_reversed, mean_signed_log2fc, rank_score, concordant_hits, reversed_hits, not_significant, study, platform, condition, timepoint, etc. |
| `signature_scores_extended.csv` | 21 | `score_signature.py` | Same as above for extended signature (367 genes). |
| `core_vs_extended_comparison.csv` | 21 | `score_signature.py` | Side-by-side: core vs extended hit_rate, mean_signed_log2fc, rank_score for each condition x timepoint. |
| `permutation_tests.csv` | 14 | `score_signature.py` | Permutation test results (1000 permutations) for mean_signed_log2fc of core signature against Weissberg experiments only. Columns: observed, empirical_p, n_permutations, platform, condition, timepoint. |

## Annotation tables

| File | Rows | Produced by | Description |
|------|------|-------------|-------------|
| `pathway_summary.csv` | 28 | inline script (Task 11) | Core signature genes grouped by gene_category x direction. Columns: gene_category, direction, count, genes (locus_tags), gene_names. |
| `alteromonas_nrecycling.csv` | 210 | `explore_alteromonas.py` | Significant HOT1A3 N-recycling genes per timepoint x omics_type. 65 unique genes across 5 timepoints. |
| `unnamed_genes_annotated.csv` | 47 | inline (gene_overview + ontology lookup) | Annotation of the 47 core signature genes lacking gene_name. Columns: locus_tag, direction, rank, product, category, ontology_notes. |

## Figures

| File | Produced by | Description |
|------|-------------|-------------|
| `trajectory_rnaseq.png` | `plot_trajectories.py` | 3-panel trajectory (hit_rate, mean_signed_log2fc, rank_score) for RNA-seq: axenic (dashed) vs coculture (solid), core + extended, with reference baseline bands. |
| `trajectory_proteomics.png` | `plot_trajectories.py` | Same for proteomics. |
| `score_overview.png` | `plot_score_overview.py` | 6-panel strip plot (3 metrics x 2 signature types). All timepoints as jittered dots, grouped by condition, colored by platform. |
