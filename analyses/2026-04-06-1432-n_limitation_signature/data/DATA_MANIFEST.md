# Data Manifest

All files produced by extraction scripts from the multiomics KG.
Run scripts from `multiomics_research` root with `uv run`.

## Reference DE extracts

| File | Rows | Genes | Timepoints | Produced by | Description |
|------|------|-------|------------|-------------|-------------|
| `de_reference_tolonen_ndep.csv` | 10,182 | 1,697 | 0h, 3h, 6h, 12h, 24h, 48h | `extract_reference_de.py` | Tolonen 2006 N-deprivation time course (microarray). All genes, all timepoints including non-significant. |
| `de_reference_read_ndep.csv` | 2,534 | 857 | 3h, 12h, 24h | `extract_reference_de.py` | Read 2017 N-depleted vs N-replete (RNA-seq). **Filtered to top 50% by expression** (table_scope: filtered_subset). |
| `de_reference_tolonen_cyanate.csv` | 1,697 | 1,697 | single | `extract_reference_de.py` | Tolonen 2006 cyanate as sole N source (microarray). Supplementary, not used for signature building. |
| `de_reference_tolonen_urea.csv` | 1,697 | 1,697 | single | `extract_reference_de.py` | Tolonen 2006 urea as sole N source (microarray). Supplementary, not used for signature building. |

## Weissberg 2025 DE extracts

| File | Rows | Genes | Timepoints | Produced by | Description |
|------|------|-------|------------|-------------|-------------|
| `de_weissberg_med4_rnaseq_axenic.csv` | 1,849 | 1,849 | single (day ~11) | `extract_weissberg_de.py` | MED4 RNA-seq axenic. Single contrast: starvation vs exponential. Labeled `timepoint=single` in data. |
| `de_weissberg_med4_rnaseq_coculture.csv` | 9,245 | 1,849 | day 18, 31, 60, 89, days 60+89 | `extract_weissberg_de.py` | MED4 RNA-seq coculture with HOT1A3. Time course. |
| `de_weissberg_med4_proteomics_axenic.csv` | 4,272 | 1,424 | day 14, 31, 89 | `extract_weissberg_de.py` | MED4 proteomics axenic. Fewer genes than RNA-seq (1,424 vs 1,849). |
| `de_weissberg_med4_proteomics_coculture.csv` | 7,120 | 1,424 | day 18, 31, 60, 89, days 60+89 | `extract_weissberg_de.py` | MED4 proteomics coculture with HOT1A3. Time course. |
| `de_weissberg_med4_all.csv` | 22,486 | — | all | `extract_weissberg_de.py` | Combined file with `source_file` column. Convenience for cross-experiment queries. |
| `de_weissberg_hot1a3_nrecycling.csv` | 1,460 | 171 | day 18, 31, 60, 89, days 60+89 | `explore_alteromonas.py` | HOT1A3 N-recycling candidate genes in coculture (RNA-seq + proteomics). |

## Signature files

| File | Rows | Produced by | Description |
|------|------|-------------|-------------|
| `core_signature_genes.csv` | 198 | `build_signature.py` | Core N-limitation signature: genes significant in BOTH Tolonen and Read, concordant direction. Columns: locus_tag, gene_name, direction, signature_type, per-study ranks, cross_study_best_dir_rank, product, gene_category. |
| `extended_signature_genes.csv` | 367 | `build_signature.py` | Extended signature: genes significant in only one study. Tagged: `tolonen_only_read_absent` (72), `tolonen_only_read_ns` (122), `read_only` (173). Same columns as core. |

## Common columns (all DE files)

Extracted with `verbose=True`:

`locus_tag, gene_name, experiment_id, treatment_type, timepoint, timepoint_hours, timepoint_order, log2fc, padj, rank, rank_up, rank_down, expression_status, product, experiment_name, treatment, gene_category, omics_type, coculture_partner, table_scope, table_scope_detail, background_factors`

- `rank`: by |log2FC| within experiment x timepoint; 1 = strongest effect
- `rank_up` / `rank_down`: directional rank among significant_up / significant_down genes. Null if not significant in that direction.
- `expression_status`: `significant_up`, `significant_down`, or `not_significant` (publication-specific thresholds)
