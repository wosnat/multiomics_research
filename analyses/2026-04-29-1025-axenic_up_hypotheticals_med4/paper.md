# Per-gene KG dossier of upregulated hypothetical-function genes in *Prochlorococcus* MED4 axenic N-starvation (Weissberg 2025)

## Question

In the *Prochlorococcus* MED4 axenic RNA-seq experiment of Weissberg et al. 2025 [1] — PRO99-lowN nitrogen starvation vs PRO99-lowN exponential growth, single contrast, `nutrient_limited` growth phase — among genes called `significant_up` by the publication's DESeq2 analysis (padj < 0.05, log2FC > 0) that are annotated as hypothetical or poorly characterized (`annotation_quality ≤ 1`), what does the KG already tell us — per gene — that bears on each gene's potential role under N-starvation?

For every gene in the candidate set, the per-gene dossier assembles: annotation (product, description, gene_category, BRITE if any), ontology terms (GO, KEGG, COG, etc.), genomic neighborhood / operon-mates, cluster assignments (any clustering analyses the gene participates in), derived metrics (every `metric_type` × `value_kind` value the KG carries for the gene), and expression in other KG studies (other organisms, treatments, timepoints). The same facets are assembled for each ortholog of each gene via its homolog group.

A stretch question is held in reserve: if patterns emerge across the dossier — shared ortholog groups, neighborhood proximity to known N-regulon members, co-occurrence in clusters, shared derived-metric profiles, or co-expression with known N-stress genes — document candidate functional bins for the upregulated hypotheticals.

A meta sub-question runs alongside the analysis and is harvested at step 6: which questions about these genes does the KG fail to answer well, and what additions to the KG would help future analyses of this kind? Concrete leads already identified at step 1: amino-acid sequences (currently absent from the KG due to per-node size limits, needed for ad-hoc BLAST / family searches / structure prediction) and batch bioinformatics layers such as Pfam / InterPro domain hits, SignalP, TMHMM transmembrane predictions, and AlphaFold structure summaries. New leads accumulate in `gaps_and_friction.md` and are consolidated in step 6 evaluation.

The analysis is RNA-seq-only and single-organism — HOT1A3 axenic, the proteome dimension, and the coculture conditions are out of scope.

## Background

The data underlying this analysis comes from Weissberg et al. 2025 [1], a multi-omics time-course of *Prochlorococcus* MED4 and *Alteromonas macleodii* HOT1A3 in PRO99-lowN medium under continuous light at 24 °C, in axenic culture and in coculture. The KG indexes 10 differential-expression experiments under this DOI; this analysis uses exactly one — the MED4 axenic RNA-seq experiment, whose sole contrast is PRO99-lowN nitrogen starvation versus PRO99-lowN exponential-phase growth (DESeq2, table_scope = `all_detected_genes`, `is_time_course = false`). The DE edges carry per-row timepoint metadata of `day 14` / `336 h`, even though the experiment node's `is_time_course` and `timepoints` fields are null.

The experiment reports 1849 detected genes; 405 are called significant_up by the publication's DESeq2 thresholds. Their gene-category breakdown (from the KG `gene_category` field) is dominated by "Unknown" (133 of 405), followed by "Stress response and adaptation" (53), "Coenzyme metabolism" (37), "Carbohydrate metabolism" (25), and "Replication and repair" (21). Restricting to genes with `annotation_quality ≤ 1` — the analysis's "hypothetical or poorly characterized" definition — yields a candidate set of **116 genes**, split 28 with `annotation_quality = 0` (no product, no description) and 88 with `annotation_quality = 1` (has description, no named product). Within the candidate set, 103 of 116 carry the "Unknown" gene_category.

Within the candidate set the log2fc distribution spans 1.02 to 6.72 (median 2.06, top-quartile cut-off 3.12); the absolute max log2fc among all 405 sig_up genes is 9.58 (PMM1404, hli paralog, AQ=2 — outside the candidate set), so the candidate set has the very highest-FC tail truncated by the AQ filter but a substantial top-quartile remains. Two top-10 candidates carry an explicit Pfam DUF assignment (PMM0958 = DUF1830, PMM0684 = DUF1651); two top-10 candidates have no ontology data at all (PMM1898, PMM1939, both AQ=0).

The candidate set's gene_category distribution is dominated by the catch-all category "Unknown" (103 / 116, 89 %); the remaining 13 candidates are spread across 6 functional categories (6 Stress response and adaptation, 2 Coenzyme metabolism, 2 Cell wall and membrane, and 1 each in Post-translational modification, Translation, Transcription). Eleven other gene_categories present in the broader 405 sig_up set (Carbohydrate, Replication, Amino acid, Transport, Nucleotide, Photosynthesis, Lipid, Regulatory, Cellular processes, Energy production, Cell cycle, Secondary metabolites, Signal transduction) contribute zero candidates — those categories are universally well-annotated. Conversely, 77 % of all "Unknown"-category sig_up genes (103 / 133) land in the candidate set.

Per-gene KG facet availability across the candidate set was snapshotted via `gene_overview`. Cross-study expression coverage is universal (116 / 116 with `expression_edge_count > 0`); ortholog-group memberships are near-universal (114 / 116); cluster memberships are near-universal (113 / 116); derived-metric coverage is sparse at 9 / 116 (consistent with the KG carrying few DerivedMetric annotations for MED4 N-stress experiments). Ontology-term presence is 99 / 116 (85 %) at the row level. The candidate-set ontology landscape was characterized at the term level via the `gene_ontology_terms` envelope (`by_term`, `by_ontology` summaries): 309 gene × term rows across 69 distinct terms, median 3 terms per gene. The three highest-frequency terms — TIGR "Hypothetical proteins / Conserved" (84 genes), COG "Function unknown" (83 genes), Cyanorak "Other > Conserved hypothetical proteins" (57 genes) — are catch-all "this is a hypothetical" terms. The same three sources, however, also carry a small tail of informative terms (e.g. COG "Coenzyme transport and metabolism", Cyanorak "Cellular processes > Adaptation/acclimation > Phosphorus / Iron / Other", TIGR "Transcription / Other"); so a per-gene "informative vs uninformative" classification cannot be done at the source level. Step 2 publishes the term-level landscape and defers per-gene informativeness assignment to step 4/5 dossier construction; an initial source-level per-gene classification was retracted (`gaps_and_friction.md` F2). Pfam terms in the candidate set are dominated by DUF / Domain-of-unknown-function assignments (35 candidates with at least one Pfam term, mostly DUFs); GO_BP / GO_MF / GO_CC and KEGG are sparse (≤ 6 each). The 17-gene subset with no ontology rows at all (F1) defines the analysis's most under-annotated slice and is the primary concrete motivation for the KG-enhancement sub-question.

A smoke test of the seven per-gene endpoints on the top-3 candidates by log2fc (PMM1828, PMM1813, PMM0958) confirmed each endpoint returns informative content (term names, cluster names, group ids, response profiles), with one exception: `gene_details` returned no sparse fields beyond what `gene_overview` already exposes for hypothetical-class genes — the planned dossier-construction surface drops `gene_details` accordingly. The three smoke-tested candidates all share the MED4 cluster "VEG (variably expressed)", and PMM0958 carries a single derived-metric annotation flagging it among the top-50 ORFs in MED4 vesicle DNA payload (Biller 2014 Table S4) — both are step-5 leads, not yet formally tested.

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
