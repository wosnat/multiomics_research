# Data Manifest

All files produced by extraction scripts from the multiomics KG.
Run scripts from `multiomics_research` root with `uv run`.

## Experiment scoping

| File | Rows | Produced by | Description |
|------|------|-------------|-------------|
| `experiment_scoping.csv` | 10 | Step 1 (interactive) | All selected experiments with classification, gene counts, table_scope, treatment_type, background_factors. |

## Reference DE extracts

| File | Rows | Genes | Timepoints | Produced by | Description |
|------|------|-------|------------|-------------|-------------|
| `de_ref_tolonen_ndep.csv` | 10,182 | 1,697 | 0h, 3h, 6h, 12h, 24h, 48h | `02_extract_reference_de.py` | Tolonen 2006 N-deprivation. all_detected_genes. |
| `de_ref_read_ndep.csv` | 2,534 | 857 | 3h, 12h, 24h | `02_extract_reference_de.py` | Read 2017 N-depleted. filtered_subset (top 50%). |

## Negative control DE extracts

| File | Rows | Genes | Timepoints | Produced by | Description |
|------|------|-------|------------|-------------|-------------|
| `de_ctrl_tolonen_cyanate.csv` | 1,697 | 1,697 | single | `02_extract_reference_de.py` | Tolonen 2006 cyanate. all_detected_genes. |
| `de_ctrl_tolonen_urea.csv` | 1,697 | 1,697 | single | `02_extract_reference_de.py` | Tolonen 2006 urea. all_detected_genes. |
| `de_ctrl_aharonovich_coculture.csv` | 1,714 | 1,714 | 20h | `02_extract_reference_de.py` | Aharonovich 2016 coculture. all_detected_genes. |
| `de_ctrl_steglich_high_white_light.csv` | 198 | 198 | single | `02_extract_reference_de.py` | Steglich 2006 high white light. filtered_subset. |

## Target DE extracts

| File | Rows | Genes | Timepoints | Produced by | Description |
|------|------|-------|------------|-------------|-------------|
| `de_weissberg_rnaseq_axenic.csv` | 1,849 | 1,849 | day 14 | `04_extract_target_de.py` | Weissberg RNA-seq axenic. all_detected_genes. |
| `de_weissberg_rnaseq_coculture.csv` | 9,245 | 1,849 | day 18, 31, 60, 89, days 60+89 | `04_extract_target_de.py` | Weissberg RNA-seq coculture. |
| `de_weissberg_proteomics_axenic.csv` | 4,272 | 1,424 | day 14, 31, 89 | `04_extract_target_de.py` | Weissberg proteomics axenic. |
| `de_weissberg_proteomics_coculture.csv` | 7,120 | 1,424 | day 18, 31, 60, 89, days 60+89 | `04_extract_target_de.py` | Weissberg proteomics coculture. |

## Signature files

| File | Rows | Produced by | Description |
|------|------|-------------|-------------|
| `core_signature.csv` | 189 | `03_build_signature.py` | Core N-limitation signature: concordant in Tolonen + Read. Columns: locus_tag, gene_name, direction, signature_type, per-study ranks, cross_study_best_dir_rank, product, gene_category. |
| `extended_signature.csv` | 342 | `03_build_signature.py` | One-study-only genes with symmetric classification (_ns/_absent). |
| `discordant_genes.csv` | 4 | `03_build_signature.py` | Genes significant in both studies with opposite direction. |

## Applied subsets

| File | Produced by | Description |
|------|-------------|-------------|
| `applied_de_*.csv` | `05_score_experiments.py` | Signature genes matched to each experiment's DE data (core tier). Shows concordance, dir_rank per gene per timepoint. |

## Common columns (all DE files)

Extracted with `verbose=True`: locus_tag, gene_name, experiment_id, treatment_type, timepoint, timepoint_hours, timepoint_order, log2fc, padj, rank, rank_up, rank_down, expression_status, product, experiment_name, treatment, gene_category, omics_type, coculture_partner, table_scope, table_scope_detail, background_factors.
