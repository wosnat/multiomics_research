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

## Signature

| File | Rows | Genes | Timepoints | Produced by | Description |
|------|------|-------|------------|-------------|-------------|
