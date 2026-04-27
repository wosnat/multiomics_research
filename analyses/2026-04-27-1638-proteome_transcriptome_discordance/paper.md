# Proteome–transcriptome discordance under nitrogen starvation in *Prochlorococcus* MED4 and *Alteromonas macleodii* HOT1A3 (Weissberg 2025)

## Question

In *Prochlorococcus* MED4 and *Alteromonas macleodii* HOT1A3 under PRO99-lowN nitrogen starvation (axenic and in coculture, Weissberg et al. 2025 [1]), where mRNA and protein responses disagree per gene at matched timepoints, what gene-level features explain the disagreement?

The analysis is scoped to a single publication (Weissberg 2025, DOI `10.1101/2025.11.24.690089`). For each organism the paired RNA-seq + proteomics fold-change comparisons that share matched timepoints are: MED4 coculture (d18, 31, 60, 89), MED4 axenic (d14 only — RNA was extractable from MED4 axenic only at d14; later proteome timepoints d31 and d89 have no matching RNA-seq and are excluded from this analysis), HOT1A3 coculture (d18, 31, 60, 89), and HOT1A3 axenic (d18, 31). Per-organism analysis (separate locus-tag spaces); cross-organism comparison runs at the pathway / ortholog-group level.

The discordance question is decomposed across multiple gene-level classification axes — pathway / function, protein size, protein architecture (cellular localization), direction-of-discordance category, ortholog conservation breadth, cross-condition consistency, cross-organism conservation, temporal pattern, operon / genomic context, hydrophobicity, and annotation quality. Step 3 framing will prioritize among these axes; later steps may surface additional caveats. The diel paired RNA-seq / proteomics dataset of Waldbauer et al. 2012 [2] is held in reserve as a step-6 cross-condition replication slot.

## Background

The data underlying this analysis comes from Weissberg et al. 2025 [1], a multi-omics time-course of *Prochlorococcus* MED4 and *Alteromonas macleodii* HOT1A3 in PRO99-lowN medium under continuous light at 24 °C, in axenic culture and in coculture with each other. The published narrative is that axenic *Prochlorococcus* perishes under extreme N-deprivation while coculture cultures persist 90+ days, supported by N-recycling activity attributed to *Alteromonas*.

The KG indexes 10 differential-expression experiments for this publication. After dropping the two single-timepoint coculture-vs-axenic RNA-seq contrasts (no matching proteomics, out of scope per step 1), the analysis uses the remaining **8 N-starvation experiments** — for each (organism × omics × condition) cell, a fold-change comparison between PRO99-lowN nutrient-starvation timepoints and PRO99-lowN exponential-phase controls. All 8 experiments use `table_scope == "all_detected_genes"`, supporting fair cross-experiment comparison.

Joining RNA-seq and proteomics on (organism, condition, timepoint_hours) produces **11 paired observations**, all in `nutrient_limited` growth phase: MED4 coculture days 18 / 31 / 60 / 89, MED4 axenic day 14, HOT1A3 coculture days 18 / 31 / 60 / 89, HOT1A3 axenic days 18 / 31. The MED4 axenic RNA-seq experiment has no per-timepoint metadata in the KG and is asserted to map to day 14 in the analysis pipeline (the only timepoint where RNA was extractable from axenic MED4 cells; later proteome timepoints have no matching RNA-seq — see `gaps_and_friction.md` F1, F2). The pooled "days 60+89" rows present in four experiments are dropped — they are statistical-power summaries of the individual late timepoints, which we use directly.

The paired-discordance gene pool — proteomics-detected ∩ RNA-seq-detected — is **1424 genes for MED4** (100% of the MED4 proteome) and **2221 genes for HOT1A3** (99.8% of the HOT1A3 proteome; the 5 missing genes are all `istB` IS21-element transposase paralogs, a known paralog-resolution mismatch between MS peptide assignment and RNA-seq multi-mapper filtering, see `gaps_and_friction.md` F3). Within each organism the RNA-seq and proteomics detected sets are identical between axenic and coculture, so axis-6 (cross-condition consistency) compares the same genes' fold-changes across conditions, not different gene sets.

The 11 classification axes for discordant genes (locked at step 1 — pathway / function, protein size, protein architecture / localization including GO Cellular Component and KEGG BRITE, direction-of-discordance, ortholog conservation breadth, cross-condition consistency, cross-organism conservation, temporal pattern, operon / genomic context, hydrophobicity, annotation quality) will be prioritized at step 3 framing.

A descriptive overview of the paired data was produced at step 2 close. Per-gene paired log2FC correlations across the paired-detection pool are: HOT1A3 Pearson r ∈ [−0.09, +0.03] at every observation; MED4 Pearson r ∈ [−0.06, +0.23], strongest at axenic d14 (+0.23), coculture d60 (+0.18), coculture d89 (+0.22). Categorizing each paired-gene observation by RNA-significance × protein-significance × direction yields, pooled across all 20 446 paired-gene-observations: concordant-significant (both omics significant in the same direction) 2.3 %, asymmetric-significant (exactly one omic significant) 29.3 %, opposite-significant (both significant, opposite directions) 1.6 %, neither significant 66.8 %. A canonical N-stress marker panel (ntcA, glnA, glnB, urtA, amt1) shows concordant_up at MED4 axenic d14 (textbook response); at MED4 coculture timepoints d18 / d31 / d60 / d89 the same five markers show RNA log2FC ≤ 0 with protein log2FC > 0, producing `prot_only_up` or `opposite_rna_down` categories. rpoB is a clean housekeeping reference (|log2FC| ≤ 0.6 in both omics across all 5 MED4 paired observations); atpA shows RNA log2FC = −4.73 at MED4 axenic d14 and is therefore not a uniformly-flat reference.

These observations [interpretation] are compatible with multiple non-exclusive mechanisms — post-transcriptional persistence, mRNA degradation / sampling artefact at coculture, or proteomics-side baseline drift — and are not yet formally tested. Step 3 framing must (a) redesign the control panel for coculture observations because the canonical N-regulon concordance positive controls do not behave concordantly there, and (b) propose competing alternative explanations as a basis for step 5 testing.

## Methods

### Framing (step 3)

*To be populated at step 3 decide-close.*

### Implementation (step 4)

*To be populated at step 4 decide-close.*

## Results

*To be populated at step 5 decide-close.*

## Discussion

*To be populated at step 6 decide-close.*

## References

[1] Weissberg O, Aharonovich D, Sher D. *Transcriptomic and Proteomic Analysis Reveals Nitrogen Recycling as a Core Mechanism for Prochlorococcus Prolonged Survival.* bioRxiv (2025). DOI: `10.1101/2025.11.24.690089`.

[2] Waldbauer JR, Rodrigue S, Coleman ML, Chisholm SW. *Transcriptome and Proteome Dynamics of a Light-Dark Synchronized Bacterial Cell Cycle.* PLoS ONE (2012). DOI: `10.1371/journal.pone.0043432`. — *Reserved for step-6 replication; not used in primary analysis.*
