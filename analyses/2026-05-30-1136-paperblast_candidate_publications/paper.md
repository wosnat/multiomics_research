# PaperBLAST-discovered candidate publications for KG ingestion

## Question

Build a deduplicated, prioritized shortlist of **~10–20+ candidate publications to add to the multi-omics KG**, discovered by running PaperBLAST on a breadth-first set of high-response genes spanning all organism types in the KG. Each candidate is triaged for what it offers (a high-throughput omics dataset vs single-protein functional characterization) and filtered against the publications already in the KG so the shortlist is genuinely new.

The end goal is **KG enrichment**: surfacing real, ingestible literature (primarily new omics datasets, secondarily functional-characterization references) that the builder repo can later load. This analysis produces the candidate list and a handoff note; it does not perform the ingestion itself.

The workflow is **hybrid**, because PaperBLAST's web interface sits behind Cloudflare bot protection (plain HTTP returns HTTP 403): (1) Claude builds a prioritized gene table from the KG; (2) the researcher runs PaperBLAST offline in a browser and saves the result pages; (3) Claude parses, dedups, and triages the downloads into the candidate-publication shortlist.

## Background

_(populated in step 2 — KG entries: organisms, experiments, the existing 42 publications, and the prioritized gene set)_

## Methods

_(populated in steps 3–4 — prioritization / dedup / triage criteria, and the gene-list builder + PaperBLAST-result parser modules)_

## Results

_(populated in step 5 — the ranked candidate-publication shortlist)_

## Discussion

_(populated in step 6 — yield assessment, what to ingest, handoff to the builder repo)_

## References

_(accumulates as KG publications are cited — resolved via `list_publications`, never from memory)_
