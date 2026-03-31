# Methods

## Research question

How do we know when Prochlorococcus is nitrogen stressed? What are the reliable transcriptional markers of N-stress, and how specific are they to nitrogen vs other stresses?

## Data scope

**Organism:** Prochlorococcus MED4 (HLI clade, NCBI taxon 59919)

**Knowledge graph:** multiomics-kg (Neo4j), queried via MCP tools and multiomics_explorer Python API.

**Nitrogen stress experiments (8):**

| Publication | Year | Platform | Design | Table scope |
|---|---|---|---|---|
| Tolonen et al. (10.1038/msb4100087) | 2006 | Microarray | N-deprivation time course (0, 3, 6, 12, 24, 48h); urea as sole N; cyanate as sole N | all_detected_genes |
| Read et al. (10.1038/ismej.2017.88) | 2017 | RNA-seq | N-depleted vs N-replete (3, 12, 24h) | filtered_subset (top 50%) |
| Weissberg et al. (10.1101/2025.11.24.690089) | 2025 | RNA-seq + Proteomics | Long-term N-starvation (14-89 days), axenic vs coculture | all_detected_genes |

**Cross-stress comparison (30 MED4 experiments total, 8 treatment types):**
nitrogen_stress (8), light_stress (8), carbon_stress (4), coculture (3), iron_stress (3), viral (3), phosphorus_stress (1), salt_stress (1)

## Gene identification

### Candidate N-stress markers (21 genes)

Identified via functional text search in KG: "nitrogen AND transport", "ammonium OR ammonia OR urea OR cyanate OR nitrite", "glutamine synthetase OR ntcA OR glnB OR PII".

Curated to 21 genes across 5 functional groups:
- **N-assimilation:** glnA (PMM0920), glsF (PMM1512)
- **N-regulation:** ntcA (PMM0246), glnB (PMM1463), pipX (PMM0393)
- **N-transport (ammonium):** amt1 (PMM0263)
- **N-transport (urea):** urtABCDE (PMM0970-0974), ureABC (PMM0963-0965)
- **N-transport (cyanate):** cynABDS (PMM0370-0373)
- **Hypotheticals (3h first responders from Tolonen):** PMM0958 (DUF1830), PMM1462 (adj. glnB), PMM0687
