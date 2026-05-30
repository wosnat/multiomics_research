# PaperBLAST-discovered candidate publications for KG ingestion

## Question

Build a deduplicated, prioritized shortlist of **~10–20+ candidate publications to add to the multi-omics KG**, discovered by running PaperBLAST on a breadth-first set of high-response genes spanning all organism types in the KG. Each candidate is triaged for what it offers (a high-throughput omics dataset vs single-protein functional characterization) and filtered against the publications already in the KG so the shortlist is genuinely new.

The end goal is **KG enrichment**: surfacing real, ingestible literature (primarily new omics datasets, secondarily functional-characterization references) that the builder repo can later load. This analysis produces the candidate list and a handoff note; it does not perform the ingestion itself.

The workflow is **hybrid**, because PaperBLAST's web interface sits behind Cloudflare bot protection (plain HTTP returns HTTP 403): (1) Claude builds a prioritized gene table from the KG; (2) the researcher runs PaperBLAST offline in a browser and saves the result pages; (3) Claude parses, dedups, and triages the downloads into the candidate-publication shortlist.

## Background

The KG holds **43 publications** across **197 experiments** and ~44 organisms
(resolved via `list_publications` / `list_experiments`), spanning *Prochlorococcus*,
marine and freshwater *Synechococcus*, and heterotroph partners (*Alteromonas*,
*Marinobacter*, *Ruegeria*, *Pseudomonas*, *Shewanella*, others). These 43
publications are the **dedup target** for the eventual candidate shortlist — a
PaperBLAST hit already represented here is not new.

The **seed-gene set** for PaperBLAST was built from differential expression across
all experiments. A gene is a "top responder" of an experiment if its KG-supplied
rank among significant rows (`rank_up`/`rank_down`, 1 = strongest |log2FC|) is in
the top 10 in any timepoint, either direction; for time-course experiments the
gene must qualify in ≥2 timepoints (transient-spike filter), for single-contrast
experiments 1 suffices. Pooled across the 155 experiments that returned gene-level
DE (42 errored inside the package's diagnostic query — see `gaps_and_friction.md`),
this gave 2290 (experiment×gene) rows → **1653 distinct responder genes**.

Genes were partitioned into **four pools** by genus — `pro` (*Prochlorococcus*),
`syn` (marine *Synechococcus* incl. GTDB *Para/Picosynechococcus* renames), `alt`
(*Alteromonas*), and `other_hetero` (all other heterotrophs) — and deduplicated
**within pool** at the ortholog-group level (finest available group: cyanorak
curated for Pro/Syn, eggNOG family/phylum/COG for heterotrophs), collapsing
paralog/operon redundancy (e.g. the seven *hli* paralogs → one) to **1430 seeds**
(pro 575 / syn 308 / alt 282 / other_hetero 265). Pools are kept separate so a
*Prochlorococcus* homolog and its *Alteromonas* counterpart can each surface their
own literature.

Field-relevant prior work in the KG cited as background:
- Tolonen et al. 2006 [@tolonen2006] — *Prochlorococcus* MED4 + MIT9313 microarray, N availability (the seed-method validation case).
- Read et al. 2017 [@read2017] — MED4 RNA-seq under N deprivation.
- Domínguez-Martín et al. 2017 [@dominguezmartin2017] — SS120 proteomics under N limitation.
- Weissberg et al. 2025 [@weissberg2025] — MED4 + *Alteromonas* HOT1A3 coculture multi-omics.

## Methods

_(populated in steps 3–4 — prioritization / dedup / triage criteria, and the gene-list builder + PaperBLAST-result parser modules)_

## Results

_(populated in step 5 — the ranked candidate-publication shortlist)_

## Discussion

_(populated in step 6 — yield assessment, what to ingest, handoff to the builder repo)_

## References

Resolved via `list_publications` (DOIs, not from memory):
- [@tolonen2006] Tolonen et al. 2006, *Mol. Syst. Biol.* — doi:10.1038/msb4100087
- [@read2017] Read et al. 2017, *ISME J.* — doi:10.1038/ismej.2017.88
- [@dominguezmartin2017] Domínguez-Martín et al. 2017, *mSystems* — doi:10.1128/mSystems.00008-17
- [@weissberg2025] Weissberg et al. 2025, *bioRxiv* — doi:10.1101/2025.11.24.690089
