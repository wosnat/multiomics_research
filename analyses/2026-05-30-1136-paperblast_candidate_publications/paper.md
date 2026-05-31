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

**Seed selection** is described in Background. **PaperBLAST run:** the 32 seed
protein IDs were submitted to PaperBLAST (papers.genomics.lbl.gov) manually in a
browser — the interface is Cloudflare-gated, so this step cannot be scripted — and
each result saved as one "HTML Only" page per seed.

**Triage pipeline** (`3_paperblast/scripts/01–03`, typed parsing of the saved
HTML, no live PaperBLAST calls):
1. *Inventory* — each saved page is mapped to its seed by the in-page query line
   (filenames are unreliable), with missing/duplicate detection.
2. *Parse + filter* — each page is parsed into (homolog hit → paper) rows
   carrying organism, %identity, and the linked papers (title, PMCID/PMID/DOI,
   year). A **pool-aware relevance filter** keeps only homologs in the seed's
   lineage: cyanobacteria for *Prochlorococcus*/*Synechococcus* (`pro`/`syn`)
   seeds, Proteobacteria (the KG heterotrophs' phylum) for `alt`/`other_hetero`
   seeds; MAG/SAG/metagenome sources are dropped (isolate genomes only). Each hit
   is tagged `in_cyanorak` by matching its strain against the Cyanorak organism
   tables — distinguishing a KG marine-picocyano organism from any other
   cyanobacterium.
3. *Dedup + triage* — papers are collapsed to distinct works; PMCID/PMID are
   resolved to DOI via the NCBI id-converter; papers already among the KG's 43
   publications are dropped; each survivor is scored for **KG-fit** (omics = 3,
   comparative-genomics/bioinformatics = 2, other = 1, single-protein
   characterization = 0) and assigned a **priority tier** (1 = in-Cyanorak
   homolog, 2 = other cyanobacterium from a `pro`/`syn` seed, 3 = heterotroph).
   The ranked shortlist sorts priority-first, then KG-fit, then cross-seed
   breadth, then identity.

KG version, the `multiomics_explorer` package, and Cyanorak source tables
(`multiomics_biocypher_kg/data`) are the data sources; the NCBI id-converter
provides PMCID→DOI resolution.

## Results

_Preliminary — pro pool only (8 of 32 seeds); the other three pools are pending
their PaperBLAST runs._

The 8 *Prochlorococcus* seeds returned 624 homolog hits; the pool-aware filter
kept 88 (cyanobacterial, isolate-genome), yielding 189 paper rows → 102 distinct
papers. Eight were already in the KG (correctly identified — Tolonen N,
the salinity transcriptome, the glucose proteomics, etc.) and dropped, leaving
**94 new candidate papers**. Of these, 12 have an in-Cyanorak homolog (priority 1)
and 9 score KG-fit ≥2 (5 omics, 4 comparative-genomics).

The two ranking axes pull apart for *Prochlorococcus*: the in-scope (in-Cyanorak)
papers are mostly not ingestible datasets — the MED4 genome paper, ncRNA surveys,
ProPortal, regulon inference — while the genuine omics/genomics datasets are about
other cyanobacteria (*Synechocystis*, *Crocosphaera*, *Anabaena*, *Leptolyngbya*),
surfaced because the conserved Pro seed genes hit their homologs. This is
consistent with *Prochlorococcus* being heavily studied and already well
represented among the KG's 43 publications; the heterotroph pools are expected to
be the richer source of novel ingestible papers.

## Discussion

_(populated in step 6 — yield assessment, what to ingest, handoff to the builder repo)_

## References

Resolved via `list_publications` (DOIs, not from memory):
- [@tolonen2006] Tolonen et al. 2006, *Mol. Syst. Biol.* — doi:10.1038/msb4100087
- [@read2017] Read et al. 2017, *ISME J.* — doi:10.1038/ismej.2017.88
- [@dominguezmartin2017] Domínguez-Martín et al. 2017, *mSystems* — doi:10.1128/mSystems.00008-17
- [@weissberg2025] Weissberg et al. 2025, *bioRxiv* — doi:10.1101/2025.11.24.690089
