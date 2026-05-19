# Nitrogen-related gene inventory across cyanobacterial strains in the KG

## Question

What nitrogen-related genes — anchored by Cyanorak functional roles `E.4` (Central intermediary metabolism > Nitrogen metabolism) and `D.1.3` (Cellular processes > Adaptation > Nitrogen) across the 13 Cyanorak-annotated cyanobacterial strains in the KG, and expanded via their eggnog Cyanobacteria-level ortholog groups — exist across the 19 cyanobacterial genome_strains in the KG, and how does this inventory differ across *Prochlorococcus* clades (HL vs LL) and across *Synechococcus* lineages?

The question is scoped to **annotation inventory and cross-strain comparison**, with the gene-set expressed as ortholog groups rather than per-strain locus tags. Expression behavior under N perturbation or other conditions is out of scope for this analysis.

## Background

This analysis draws on **19 cyanobacterial genome strains** in the multi-omics KG: 11 *Prochlorococcus* (clades HLI, HLII, LLI, LLII, LLIV; including the strains AS9601, MED4, MIT0801, MIT9301, MIT9303, MIT9312, MIT9313, NATL1A, NATL2A, RSP50, SS120), 4 marine *Synechococcus* (WH8102, WH7803, CC9311, BL107 — Cyanorak SubCluster 5.1, Clades I/III/IV/V), and 4 cyano outside the Cyanorak picocyano scope (*Synechococcus* PCC 7002, *S. elongatus* PCC 7942 and UTEX 2973, *Thermosynechococcus vestitus* BP-1). Clade assignments use the Cyanorak source CSV; the KG `Organism.clade` field is mis-populated for *Synechococcus* (null) and for NATL1A/NATL2A (LLII instead of LLI), per `gaps_and_friction.md`.

The **anchor gene set** is the union of Cyanorak functional roles `E.4` (Central intermediary metabolism > Nitrogen metabolism) and `D.1.3` (Cellular processes > Adaptation > Nitrogen) across the 13 Cyanorak-annotated strains in scope — 375 distinct locus_tags total, ranging from 12 (SS120) to 43 (WH8102, CC9311) per strain. The anchor set is then expanded via **eggnog Cyanobacteria-level ortholog groups** (`taxonomic_level='Cyanobacteria'`, specificity_rank=2) to reach the 6 non-Cyanorak strains; 9 anchor locus_tags lacking a Cyanobacteria-level eggnog group are recovered as **synthetic singleton ortholog groups** (one per strain × gene_name combination) so no anchor is dropped. The resulting inventory is a 19 strain × 61 ortholog group matrix with copy-count cells, sourced from 54 eggnog Cyanobacteria-level groups + 7 anchor singletons.

Field-relevant N-perturbation literature in the KG (cited as Background; experiments not analyzed here):
- Tolonen et al. 2006 [@tolonen2006] — MED4 + MIT9313 microarray time-course under N depletion.
- Read et al. 2017 [@read2017] — MED4 RNA-seq under N deprivation with transcription start site mapping.
- Domínguez-Martín et al. 2017 [@dominguezmartin2017] — SS120 proteomics under N limitation.
- Weissberg et al. 2025 [@weissberg2025] — MED4 + *Alteromonas* HOT1A3 coculture multi-omics, framed around N recycling.

Two strains warrant a Background note. **SS120 lacks the entire urea pathway** (urtABCDE transporter + ureABCDEFG urease) — verified directly from the KG anchor set (12 anchor genes vs 24–28 in other *Prochlorococcus*); the 14 missing genes are all real biological absences, not annotation gaps. **WH7803** shows a similar pattern in marine *Synechococcus*: lacks urease and most urea transport while keeping the nitrate/nitrite pathway (nrtP, narB, narM, nirA). Both losses are consistent with the established field literature on these strains.

## Methods

**Gene-set construction.** The inventory was assembled by querying the multi-omics KG via the `multiomics_explorer` Python API (mirror of the MCP tools). Per-strain Cyanorak `E.4 ∪ D.1.3` membership was pulled with `genes_by_ontology`. Each anchor locus_tag was then mapped to its eggnog Cyanobacteria-level ortholog group (`gene_homologs(source='eggnog', taxonomic_level='Cyanobacteria')`); 366 / 375 (97.6%) mapped to one of 54 unique eggnog groups. The remaining 9 orphan anchor locus_tags (no Cyanobacteria-level group) were recovered as **synthetic singleton groups** (one per (strain, gene_name) pair) so no anchor was dropped; this yielded 7 additional groups for a final inventory dimension of 19 strains × 61 ortholog groups. Group members across all 19 cyano strains were pulled with `genes_by_homolog_group`. Strain clade assignments were taken from the Cyanorak source CSVs (`prochlorococcus.csv`, `synechococcus.csv`) rather than the KG `Organism.clade` field, which is mis-populated for *Synechococcus* (null) and NATL1A/NATL2A (LLII instead of LLI).

**Cross-strain comparison.** The 19 × 61 copy-count matrix was clustered hierarchically on both axes. The clustering used **Jaccard distance** computed on the binarized presence/absence matrix (cell > 0 → 1), with **UPGMA (average linkage)** for both rows and columns. The binarized clustering pairs with copy-count cell encoding in the heatmap, which keeps paralog signal visible without committing to a paralog-distance metric. The colormap is discrete (0, 1, 2, ≥3) and cells with copy_count ≥ 2 carry an integer overlay. Row annotations: canonical clade (Cyanorak), genus, Cyanorak Pigment type (where applicable), genome gene count.

**Reproducibility.** All extraction, clustering, and rendering live in `analyses/2026-05-19-1037-cyano_n_gene_inventory/` with one script per pipeline stage (`2_kg_selection/scripts/01..06`, `4_methods/heatmap_clustering.py`, `5_analyze/scripts/01_cluster_and_plot.py`). The clustering module was toy-data-verified against hand-computed Jaccard distances before running on the real inventory.

## Results

The clustered heatmap is presented in `5_analyze/figures/01_clustered_heatmap.png` (see Figure 1). Across 19 cyanobacterial genome strains and 61 N-gene ortholog groups, 539 of 1159 (46.5%) strain × group cells are present, with 18 cells carrying paralog copies (copy_count ≥ 2).

**Strain clustering recovers taxonomy.** The data-driven dendrogram splits cleanly into two top-level clades: all 11 *Prochlorococcus* in one branch, all 8 non-*Prochlorococcus* (4 marine *Synechococcus* + 1 *Synechococcus* PCC 7002 + 2 *S. elongatus* + 1 *Thermosynechococcus*) in the other. Within *Prochlorococcus*, sub-clades pair as expected: LLIV (MIT9303 + MIT9313), LLI (NATL2A + MIT0801 + NATL1A), HLII (MIT9301 + AS9601 + MIT9312), HLI (MED4 + RSP50). The two systematic deviations are biologically explainable: SS120 (LLII) separates from the LL group because of complete urea-pathway loss (no urtABCDE + no ureABCDEFG); WH7803 (marine *Synechococcus* clade V) separates from the other marine *Synechococcus* because of urease loss with partial urea-transport loss.

**Ortholog-group structure reflects functional units.** The column clustering recovers operon-level co-occurrence: the urea pathway (urtABCDE + ureABCDEFG) cluster tightly — strains either have the entire 12-gene module or none of it — consistent with horizontal acquisition / loss as a unit. The nitrate / nitrite pathway (nrtP, narB, narM, nirA, focA) cluster as a separate unit, present in marine *Synechococcus* + PCC 7002 + *Thermosynechococcus* + LL Prochlorococcus NATL2A and absent from HL *Prochlorococcus*. The cyanate operon (cynABDSH) splits into Pro-lineage and Syn-lineage paralog sub-clusters — the orthology bridge correctly preserves this lineage asymmetry.

**Universally-conserved N machinery is the dense right block of the heatmap.** Seven ortholog groups are present in all 19 strains: NtcA (global N regulator), GlnB (P-II regulator), GlnA (glutamine synthetase), GlsF (ferredoxin-dependent glutamate synthase), CarA (carbamoyl-P synthase small subunit), MetC (cystathionine β-lyase), and MoeB (molybdopterin biosynthesis). These define the cyano-wide core N regulatory + assimilation machinery.

**Notable paralog signal.** Of the 18 multi-copy cells, the highest is *Synechococcus* PCC 7002 amt1 × 3. *Prochlorococcus* MIT9303 carries a duplicated ntcA — a Crp/FNR-family paralog (`ptrA`, P9303_11071) in the same Cyanobacteria-level ortholog group as the canonical ntcA. The cysK ortholog group has paralog copies in all 4 marine *Synechococcus*. WH8102 carries a duplicated urtA.

### Pathway-level cross-strain comparison

Collapsing the 61 ortholog groups into 10 N-pathway categories (regulation, assimilation core, ammonium uptake, urea uptake, urease, cyanate, nitrate/nitrite, Mo cofactor biosynthesis, Met biosynthesis, other) yields three complementary views (Figures 2–4).

**Figure 2 (Strain × Pathway).** With strains ordered by Figure 1's clustering, the strain × pathway heatmap shows pathway-level coverage at a glance. SS120 is the visible outlier: zero coverage in urea uptake, urease, cyanate, nitrate/nitrite, and Mo cofactor. The two *Synechococcus elongatus* strains (PCC 7942 and UTEX 2973) have identical pathway profiles. PCC 7002 is the broadest non-marine cyano.

**Figure 3 (Clade-group × Pathway summary).** Aggregating strains into 8 clade groups (HLI, HLII, LLI, LLII, LLIV, marine *Synechococcus*, non-marine *Synechococcus*, *Thermosynechococcus*), the highest-level cross-clade comparison surfaces several patterns:
- All *Prochlorococcus* clades except LLII (=SS120) retain ~70% of urea uptake and ~70% of urease.
- Cyanate utilization is HLI *Prochlorococcus*-specific (38% in HLI; 0% in HLII / LLII / LLIV; 12% in LLI from MIT0801 alone) — the canonical MED4-lineage cyanate operon (cynABDS).
- Nitrate / nitrite is a *Synechococcus* / *Thermosynechococcus* signature: marine *Syn* 68%, non-marine *Syn* 52%, *Thermosynechococcus* 43%; absent in HL *Prochlorococcus*; partial in LL *Prochlorococcus* (29%).
- Mo cofactor biosynthesis is broadly conserved in *Synechococcus* (~60-90%) but rare in *Prochlorococcus* (~11%).
- Met biosynthesis is universal in marine *Synechococcus* (100%), high in *Prochlorococcus* (~80%), reduced in non-marine *Synechococcus* (47%), and minimal in *Thermosynechococcus* (20%) — *Thermosynechococcus* uses a different Met biosynthesis route.

**Figure 4 (Per-strain composition).** Stacked horizontal bars of pathway-group counts per strain, sorted descending by total. CC9311 (43 groups) and WH8102 (41) are the broadest; SS120 (12) is the most reduced. Strains with higher total N-gene counts diversify across more pathways (not just more copies of the same pathway).

## Discussion

*[Populated in step 6.]*

## References

All references resolved via `list_publications`; DOIs verified against the KG, not from intrinsic knowledge.

1. **Tolonen et al. 2006.** Global gene expression of *Prochlorococcus* ecotypes in response to changes in nitrogen availability. *Molecular Systems Biology*. doi:10.1038/msb4100087. KG: 6 microarray experiments, time-course, *Prochlorococcus* MED4 and MIT9313.
2. **Read et al. 2017.** Nitrogen cost minimization is promoted by structural changes in the transcriptome of N-deprived *Prochlorococcus* cells. *The ISME Journal*. doi:10.1038/ismej.2017.88. KG: 1 RNA-seq experiment with transcription start site mapping, MED4.
3. **Domínguez-Martín et al. 2017.** Quantitative Proteomics Shows Extensive Remodeling Induced by Nitrogen Limitation in *Prochlorococcus marinus* SS120. *mSystems*. doi:10.1128/mSystems.00008-17. KG: 1 proteomics experiment, SS120.
4. **Weissberg et al. 2025.** Transcriptomic and Proteomic Analysis Reveals Nitrogen Recycling as a Core Mechanism for *Prochlorococcus* Prolonged Survival. *bioRxiv*. doi:10.1101/2025.11.24.690089. KG: 10 multi-omics (RNA-seq + proteomics) experiments, MED4 axenic + *Alteromonas* HOT1A3 coculture, N recycling.
