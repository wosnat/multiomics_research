# Nitrogen-related gene inventory across cyanobacterial strains in the KG

## Question

What nitrogen-related genes — anchored by Cyanorak functional roles `E.4` (Central intermediary metabolism > Nitrogen metabolism) and `D.1.3` (Cellular processes > Adaptation > Nitrogen) across the 13 Cyanorak-annotated cyanobacterial strains in the KG, and expanded via their eggnog Cyanobacteria-level ortholog groups — exist across the 19 cyanobacterial genome_strains in the KG, and how does this inventory differ across *Prochlorococcus* clades (HL vs LL) and across *Synechococcus* lineages?

The question is scoped to **annotation inventory and cross-strain comparison**, with the gene-set expressed as ortholog groups rather than per-strain locus tags. Expression behavior under N perturbation or other conditions is out of scope for this analysis.

## Background

This analysis draws on **19 cyanobacterial genome strains** in the multi-omics KG: 11 *Prochlorococcus* (clades HLI, HLII, LLI, LLII, LLIV; including the strains AS9601, MED4, MIT0801, MIT9301, MIT9303, MIT9312, MIT9313, NATL1A, NATL2A, RSP50, SS120), 4 marine *Synechococcus* (WH8102, WH7803, CC9311, BL107 — Cyanorak SubCluster 5.1, Clades I/III/IV/V), and 4 cyano outside the Cyanorak picocyano scope (*Synechococcus* PCC 7002, *S. elongatus* PCC 7942 and UTEX 2973, *Thermosynechococcus vestitus* BP-1). Clade assignments use the Cyanorak source CSV; the KG `Organism.clade` field is mis-populated for *Synechococcus* (null) and for NATL1A/NATL2A (LLII instead of LLI), per `gaps_and_friction.md`.

The **anchor gene set** is the union of Cyanorak functional roles `E.4` (Central intermediary metabolism > Nitrogen metabolism) and `D.1.3` (Cellular processes > Adaptation > Nitrogen) across the 13 Cyanorak-annotated strains in scope — 375 distinct locus_tags total, ranging from 12 (SS120) to 43 (WH8102, CC9311) per strain. The anchor set is then expanded via **eggnog Cyanobacteria-level ortholog groups** (`taxonomic_level='Cyanobacteria'`, specificity_rank=2) to reach the 6 non-Cyanorak strains. The resulting inventory is expressed as ortholog groups, not per-strain locus tags.

Field-relevant N-perturbation literature in the KG (cited as Background; experiments not analyzed here):
- Tolonen et al. 2006 [@tolonen2006] — MED4 + MIT9313 microarray time-course under N depletion.
- Read et al. 2017 [@read2017] — MED4 RNA-seq under N deprivation with transcription start site mapping.
- Domínguez-Martín et al. 2017 [@dominguezmartin2017] — SS120 proteomics under N limitation.
- Weissberg et al. 2025 [@weissberg2025] — MED4 + *Alteromonas* HOT1A3 coculture multi-omics, framed around N recycling.

Two strains warrant a Background note. **SS120 lacks the entire urea pathway** (urtABCDE transporter + ureABCDEFG urease) — verified directly from the KG anchor set (12 anchor genes vs 24–28 in other *Prochlorococcus*); the 14 missing genes are all real biological absences, not annotation gaps. **WH7803** shows a similar pattern in marine *Synechococcus*: lacks urease and most urea transport while keeping the nitrate/nitrite pathway (nrtP, narB, narM, nirA). Both losses are consistent with the established field literature on these strains.

## Methods

*[Populated in steps 3–4: framing decisions, comparison axes, methods module.]*

## Results

*[Populated in step 5.]*

## Discussion

*[Populated in step 6.]*

## References

All references resolved via `list_publications`; DOIs verified against the KG, not from intrinsic knowledge.

1. **Tolonen et al. 2006.** Global gene expression of *Prochlorococcus* ecotypes in response to changes in nitrogen availability. *Molecular Systems Biology*. doi:10.1038/msb4100087. KG: 6 microarray experiments, time-course, *Prochlorococcus* MED4 and MIT9313.
2. **Read et al. 2017.** Nitrogen cost minimization is promoted by structural changes in the transcriptome of N-deprived *Prochlorococcus* cells. *The ISME Journal*. doi:10.1038/ismej.2017.88. KG: 1 RNA-seq experiment with transcription start site mapping, MED4.
3. **Domínguez-Martín et al. 2017.** Quantitative Proteomics Shows Extensive Remodeling Induced by Nitrogen Limitation in *Prochlorococcus marinus* SS120. *mSystems*. doi:10.1128/mSystems.00008-17. KG: 1 proteomics experiment, SS120.
4. **Weissberg et al. 2025.** Transcriptomic and Proteomic Analysis Reveals Nitrogen Recycling as a Core Mechanism for *Prochlorococcus* Prolonged Survival. *bioRxiv*. doi:10.1101/2025.11.24.690089. KG: 10 multi-omics (RNA-seq + proteomics) experiments, MED4 axenic + *Alteromonas* HOT1A3 coculture, N recycling.
