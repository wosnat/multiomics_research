# Results Manifest

All files produced by scoring, plotting, and analysis scripts.

## Score tables

| File | Rows | Produced by | Description |
|------|------|-------------|-------------|
| `scores_all.csv` | 72 | `05_score_experiments.py` | All experiments × timepoints × tiers (top/core/extended) scored. Columns: rank_score, empirical_p, n_matched, n_absent, n_concordant, n_reversed, n_not_significant, hit_rate, tier, label, role, timepoint, timepoint_hours, total_genes_in_experiment. |

## Figures

| File | Produced by | Description |
|------|-------------|-------------|
| `trajectory_rnaseq.png` | `06_plot_results.py` | 2-panel: rank score and hit rate over time for RNA-seq axenic (star) vs coculture (line), with reference range band. Core tier. |
| `trajectory_proteomics.png` | `06_plot_results.py` | Same for proteomics. Shows coculture scoring higher than axenic. |
| `control_separation.png` | `06_plot_results.py` | 2-panel: rank score and hit rate by experiment role (reference/negative/target). Core tier. Annotates outliers. |
| `tier_comparison.png` | `06_plot_results.py` | Rank score for all target experiments, top vs core vs extended tiers side by side. |
