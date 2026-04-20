# Data Manifest

All files produced by extraction scripts from the multiomics KG.
Run scripts from `multiomics_research` root with `uv run`.

## Experiment landscape

| File | Rows | Genes | Timepoints | Produced by | Description |
|------|------|-------|------------|-------------|-------------|
| `experiments_classified.csv` | 37 (15 unique exps) | n/a | varies (1–6 per exp) | `01_select_experiments.py` | T/R/PC/NC/CTX-classified experiments joined with `list_experiments` metadata + `list_publications` first_author / title attribution. 4 T (Weissberg), 2 R (Tolonen 2006, Read 2017), 2 PC (Tolonen alt-N), 5 NC (Aharonovich, Weissberg, Moreno-Cabezuelo ×2, Steglich), 2 CTX (Tolonen MIT9313, Domínguez-Martín SS120). HOT1A3 CTX1 dropped 2026-04-20 during Step 1b ontology-compatibility review — see notebook Step 1a redo entry. |

## Enrichment outputs

| File | Rows | Genes | Timepoints | Produced by | Description |
|------|------|-------|------------|-------------|-------------|

## Signature

| File | Rows | Genes | Timepoints | Produced by | Description |
|------|------|-------|------------|-------------|-------------|
