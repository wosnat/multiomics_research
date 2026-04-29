# Per-gene KG dossier of upregulated hypothetical-function genes in *Prochlorococcus* MED4 axenic N-starvation (Weissberg 2025)

## Question

In the *Prochlorococcus* MED4 axenic RNA-seq experiment of Weissberg et al. 2025 [1] — PRO99-lowN nitrogen starvation vs PRO99-lowN exponential growth, single contrast, `nutrient_limited` growth phase — among genes called `significant_up` by the publication's DESeq2 analysis (padj < 0.05, log2FC > 0) that are annotated as hypothetical or poorly characterized (`annotation_quality ≤ 1`), what does the KG already tell us — per gene — that bears on each gene's potential role under N-starvation?

For every gene in the candidate set, the per-gene dossier assembles: annotation (product, description, gene_category, BRITE if any), ontology terms (GO, KEGG, COG, etc.), genomic neighborhood / operon-mates, cluster assignments (any clustering analyses the gene participates in), derived metrics (every `metric_type` × `value_kind` value the KG carries for the gene), and expression in other KG studies (other organisms, treatments, timepoints). The same facets are assembled for each ortholog of each gene via its homolog group.

A stretch question is held in reserve: if patterns emerge across the dossier — shared ortholog groups, neighborhood proximity to known N-regulon members, co-occurrence in clusters, shared derived-metric profiles, or co-expression with known N-stress genes — document candidate functional bins for the upregulated hypotheticals.

A meta sub-question runs alongside the analysis and is harvested at step 6: which questions about these genes does the KG fail to answer well, and what additions to the KG would help future analyses of this kind? Concrete leads already identified at step 1: amino-acid sequences (currently absent from the KG due to per-node size limits, needed for ad-hoc BLAST / family searches / structure prediction) and batch bioinformatics layers such as Pfam / InterPro domain hits, SignalP, TMHMM transmembrane predictions, and AlphaFold structure summaries. New leads accumulate in `gaps_and_friction.md` and are consolidated in step 6 evaluation.

The analysis is RNA-seq-only and single-organism — HOT1A3 axenic, the proteome dimension, and the coculture conditions are out of scope.

## Background

*To be populated at step 2 decide-close.*

## Methods

### Framing (step 3)

*To be populated at step 3 decide-close.*

### Implementation (step 4)

*To be populated at step 4 decide-close.*

## Results

*To be populated at step 5 decide-close.*

## Discussion

*To be populated at step 6 decide-close. Includes the consolidated KG-enhancement proposals from `gaps_and_friction.md`.*

## References

[1] Weissberg O, Aharonovich D, Sher D. *Transcriptomic and Proteomic Analysis Reveals Nitrogen Recycling as a Core Mechanism for Prochlorococcus Prolonged Survival.* bioRxiv (2025). DOI: `10.1101/2025.11.24.690089`. — KG experiment used: `10.1101/2025.11.24.690089_growth_state_pro99lown_nutrient_starvation_med4_rnaseq_axenic`.
