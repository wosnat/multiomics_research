# Results Manifest

All files produced by scoring, plotting, and analysis scripts.

---

## Ontology ranking

| File | Script | Description |
|------|--------|-------------|
| `ontology_ranking.csv` | `02_survey_landscape.py` | Ranked comparison of all 9 ontologies on genome coverage, term count, median term size, max term size at the best hierarchy level. CyanoRak level 1 ranked first. |

---

## Enrichment tests

| File | Script | Description |
|------|--------|-------------|
| `enrichment_all.csv` | `04_run_enrichment.py` | All 5,589 Fisher's exact test results. Columns: `pathway_id`, `pathway_name`, `experiment`, `timepoint`, `direction` (up/down/combined), `n_de`, `n_pathway`, `n_de_in_pathway`, `n_background`, `odds_ratio`, `pval`, `padj`. |
| `enrichment_significant.csv` | `04_run_enrichment.py` | Subset of `enrichment_all.csv` where padj < 0.05. 133 rows. |

---

## Figures

### Primary figure

| File | Description |
|------|-------------|
| `heatmap_enrichment_signed.png` | Signed enrichment score heatmap: all 18 enriched pathways × all experiments × timepoints. Color encodes signed -log10(padj). Positive = up-enriched, negative = down-enriched. Primary result figure. |
| `heatmap_enrichment_signed.svg` | Vector format of above. |

### Supporting figures

| File | Description |
|------|-------------|
| `heatmap_enrichment_up.png` | Up-enriched pathways only. Supplementary. |
| `heatmap_enrichment_up.svg` | Vector format of above. |
| `heatmap_enrichment_down.png` | Down-enriched pathways only. Supplementary. |
| `heatmap_enrichment_down.svg` | Vector format of above. |
| `discordance_scatter.png` | RNA vs protein enrichment scatter (coculture conditions). X-axis: RNA-seq coculture signed score (mean across timepoints). Y-axis: proteomics coculture signed score. Labels: E.4 (N-metabolism), Q.1 (amino acid transport), K.2 (ribosomal proteins). Illustrates discordance finding. |
| `discordance_scatter.svg` | Vector format of above. |
| `trajectory_proteomics_coculture.png` | E.4, Q.1, K.2 signed enrichment score across proteomics coculture timepoints (d5, d18, d31, d60, d89). Shows persistence of N-stress and ribosomal shutdown signal. |
| `trajectory_proteomics_coculture.svg` | Vector format of above. |
| `control_separation.png` | Negative control vs target experiment enrichment comparison. Confirms clean separation for Aharonovich, Steglich, and Tolonen urea; shows partial signal in Tolonen cyanate (photosynthesis only). |
| `control_separation.svg` | Vector format of above. |
