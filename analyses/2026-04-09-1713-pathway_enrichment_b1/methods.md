# Methods

## Research question

Which functional pathways are enriched among differentially expressed genes in Weissberg 2025 conditions (MED4 axenic N-deprivation, MED4/Alteromonas coculture N-deprivation), compared to reference N-deprivation experiments and negative controls?

## Data scope

10 experiments from Weissberg 2025 v2 differential expression data:

**Reference N-deprivation (positive controls):**
- Tolonen N-dep (RNA-seq, multiple timepoints)
- Read N-dep (RNA-seq, multiple timepoints)

**Target conditions (Weissberg 2025):**
- Weissberg RNA-seq axenic N-dep
- Weissberg RNA-seq coculture N-dep
- Weissberg proteomics axenic N-dep
- Weissberg proteomics coculture N-dep

**Negative controls:**
- Tolonen cyanate (poor N-source, partial N-stress)
- Tolonen urea (non-N-limiting N-source)
- Aharonovich coculture (different organism pair, no N-stress)
- Steglich high light (light stress, no N perturbation)

Ontology: CyanoRak level 1, 110 terms. 69 terms retained after minimum-size filter (≥5 genes). Full MED4 genome: 1,976 genes annotated.

## Gene selection

The full MED4 genome (1,976 genes) was retrieved from the multiomics KG and annotated with all available ontologies (9 ontologies extracted, CyanoRak selected; see decisions.md). The background set for each Fisher's exact test is the set of genes present in the differential expression table for that experiment × timepoint — i.e., genes that were quantified and passed QC for that dataset. This varies by experiment and timepoint (typically 1,500–1,900 genes for RNA-seq; 600–1,200 genes for proteomics).

## Statistical methods

**Test:** Fisher's exact test (one-sided, testing over-representation), applied per pathway × experiment × timepoint × direction (up-regulated DE genes, down-regulated DE genes, and combined DE genes regardless of direction). Differentially expressed genes defined as the DE set in each v2 `de_*.csv` file.

**Multiple testing correction:** Benjamini-Hochberg FDR applied within each experiment × timepoint × direction group. This controls within-group FDR but does not correct across experiments.

**Significance threshold:** padj < 0.05.

**Signed enrichment score:** For visualization, each pathway × experiment × timepoint is summarized by a signed score: sign × -log10(padj), where sign is +1 if the dominant enrichment is up-regulated and −1 if down. When both up and down pass the threshold, the direction with the smaller padj determines the sign. Score = 0 when neither direction is significant.

## Results summary

5,589 tests total. 133 significant at padj < 0.05 (2.4%). 18 unique pathways enriched in at least one experiment × timepoint.

Key finding: RNA/protein discordance observed in Approach A (reference signature scoring) resolves at the pathway level. N-metabolism (E.4) and amino acid transport (Q.1) are strongly enriched among up-regulated proteins in proteomics coculture (padj 1e-15 and 5e-5, respectively) but are not enriched in RNA-seq coculture. Ribosomal proteins (K.2) are enriched down in proteomics coculture (padj 3e-23 at d31) but not in RNA-seq coculture.

Photosynthesis pathways (J.1 ATP synthase, J.2 CO2 fixation, J.7 PSI, J.8 PSII) are enriched down in references and in RNA-seq axenic, consistent with canonical N-deprivation physiology.

Negative controls are largely clean: Aharonovich coculture and Steglich high light show 0 significant enrichments. Tolonen urea shows 2 (the catch-all D.1 adaptation). Tolonen cyanate shows 8 — photosynthesis pathways enriched down, consistent with partial N-stress from a poor nitrogen source, but N-metabolism (E.4) is NOT enriched, distinguishing cyanate from true N-limitation.

## Limitations

- **Single ontology:** Results are specific to CyanoRak. Different ontologies (KEGG, GO BP) may reveal additional or different pathways, though CyanoRak was selected for its superior genome coverage and biological specificity for cyanobacteria (see decisions.md).
- **20% genome unannotated:** ~470 MED4 genes have no CyanoRak annotation and are excluded from enrichment tests. These are predominantly uncharacterized hypothetical proteins.
- **Catch-all categories:** R.2 (Conserved hypothetical proteins, 340 genes) and D.1 (Adaptation/acclimation, 213 genes) are large heterogeneous categories with high statistical power. Enrichment in these terms should not be over-interpreted.
- **No cross-experiment FDR correction:** BH-FDR is applied within each experiment × timepoint × direction group. The 133 significant results are not corrected for testing across 10 experiments, which inflates the family-wise type I error rate.
- **Steglich low power:** Steglich high light has only 198 genes in its DE universe (49/69 pathways represented, mean coverage 10%). Absence of enrichment in this experiment reflects low power, not biological absence of signal.
- **Signed score information loss:** When both up and down directions are significant for a pathway in the same experiment × timepoint, the signed score retains only the dominant direction. This edge case was rare in these data but should be noted.
