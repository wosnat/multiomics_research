# Data Manifest

All files produced by extraction scripts from the multiomics KG.
Run scripts from `multiomics_research` root with `uv run`.

## Experiment scoping

| File | Rows | Produced by | Description |
|------|------|-------------|-------------|
| `experiment_scoping.csv` | 10 | Step 1 (interactive KG discovery) | All selected experiments with classification (reference/negative_control/target), gene counts, table_scope, treatment_type, background_factors. |
