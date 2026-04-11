# References and Citations

## Publications

- **Weissberg et al. 2025** — The primary experimental data source. MED4 N-deprivation under axenic and coculture (with Alteromonas macleodii) conditions. RNA-seq and proteomics. bioRxiv preprint.
- **Tolonen et al.** — Reference N-deprivation RNA-seq dataset for Prochlorococcus MED4. Used as positive control for N-limitation enrichment patterns.
- **Read et al.** — Reference N-deprivation RNA-seq dataset. Used as second positive control.
- **Aharonovich et al.** — Prochlorococcus coculture dataset (different organism pair). Used as negative control.
- **Steglich et al.** — Prochlorococcus high light dataset. Used as negative control.

## Databases

- **CyanoRak** — Functional role ontology curated for cyanobacteria. Provides the hierarchical annotation used for pathway enrichment. Three-level hierarchy (root → subcategory → leaf). Version as indexed in the multiomics KG (see KG version below). URL: http://cyanorak.sb-roscoff.fr/

## Software versions

- **scipy** — Fisher's exact test (`scipy.stats.fisher_exact`, one-sided). Version as specified in `multiomics_research` project dependencies.
- **statsmodels** — Benjamini-Hochberg FDR correction (`statsmodels.stats.multitest.multipletests`, method='fdr_bh'). Version as specified in project dependencies.
- **pandas** — Data manipulation and test result aggregation.
- **matplotlib / seaborn** — Figure generation.
- **enrich_utils** — Analysis utility package developed in this analysis (`analyses/2026-04-09-1713-pathway_enrichment_b1/enrich_utils/`). Provides annotation roll-up, enrichment runner, and scoring functions.

## KG version

- **multiomics-kg MCP server** — All gene annotations extracted via the multiomics KG Python API and MCP tools from the `multiomics_explorer` sibling repo. KG built from the `biocypher_kg` repo. Extraction date: 2026-04-09.
