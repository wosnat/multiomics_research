# Methods

## Research question

Can we distinguish nitrogen-specific stress from general nutrient stress, general stress, and cell death in the transcriptional response of Prochlorococcus MED4 to nitrogen deprivation?

## Data scope

**Organism:** Prochlorococcus MED4

**Knowledge graph:** multiomics-kg (Neo4j), queried via MCP tools.

**Experiments used for cross-stress comparison:**

All differential expression experiments for MED4 in the KG (30 experiments, 11 publications, 8 treatment types):
- nitrogen_stress (8 experiments): Tolonen 2006, Read 2017, Weissberg 2025
- iron_stress (3): Thompson 2011
- phosphorus_stress (1): Martiny 2006
- carbon_stress (4): Hopkinson 2015
- light_stress (8): Steglich 2006, Vislova 2022
- coculture (3): Aharonovich & Sher 2016, Weissberg 2025
- viral (3): Vislova 2022
- salt_stress (1): Fang 2022

## Gene identification

### Nitrogen first responders
Extracted from Tolonen 2006 microarray (10.1038/msb4100087), N-deprivation time course.
Criterion: significantly differentially expressed at 3h (the earliest non-baseline timepoint).
Result: 6 genes (PMM0370/cynA, PMM0970/urtA, PMM0920/glnA, PMM0958, PMM1462, PMM0687), all upregulated.

### N-metabolism gene set
Identified via functional text search in KG: "nitrogen AND transport", "ammonium OR ammonia OR urea OR cyanate OR nitrite".
Curated to 17 genes: amt1, cynABDS, glnA, glnB, glsF, ntcA, ureABC, urtABCDE.

## Cross-stress analysis

For each gene set, queried differential_expression_by_gene across ALL MED4 experiments (no treatment filter) to identify which stresses each gene responds to.

Classification:
- **N-specific**: significantly DE only under nitrogen_stress (and potentially coculture, which involves N-provision)
- **Nutrient-general**: responds to multiple nutrient stresses (N, P, Fe, C)
- **General stress**: responds to non-nutrient stresses (light, viral, salt)

## Caveats

- Microarray (Tolonen) and RNA-seq (Read, Weissberg) platforms have different sensitivity and dynamic range. Fold changes are not directly comparable across platforms.
- Statistical tests differ: Goldenspike (Tolonen), Rockhopper (Read), DESeq2 (Weissberg). P-value thresholds and significance calls are per-study.
- Table scopes differ: all_detected_genes (Tolonen, Weissberg), filtered_subset/top 50% (Read). Absence from Read's results could mean filtered out, not absent.
- The Tolonen 48h timepoint may represent dying cells — growth and Fv/Fm data (not in KG) would help interpret this.
