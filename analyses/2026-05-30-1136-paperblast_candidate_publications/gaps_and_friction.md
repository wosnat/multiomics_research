# Gaps and friction

Append-only log of methodology / KG / tooling friction encountered during this analysis. Distinct from decisions (which live in each step's `notebook.md`).

---

## 2026-05-30 — PaperBLAST is behind Cloudflare bot protection

**What happened.** Probing PaperBLAST (`papers.genomics.lbl.gov/cgi-bin/litSearch.cgi`) programmatically — via WebFetch and via `curl` with a browser user-agent — returned HTTP 403 with a Cloudflare "Just a moment…" JS-challenge page (`challenges.cloudflare.com`). Simple HTTP cannot reach the live PaperBLAST interface.

**Impact on methodology.** Rules out a fully-scripted end-to-end pipeline that calls PaperBLAST directly. Resolved by a **hybrid workflow**: Claude builds the gene list and parses results (scripted); the researcher runs PaperBLAST in a real browser and downloads the result pages (manual — the one non-scriptable link, documented per Rule 5). A local PaperBLAST install would restore full scripting if the volume ever justifies it.

## 2026-05-30 — Stored DE rank fields can be null

**What happened.** `Changes_expression_of` edges carry `rank_up` and `rank_by_effect`, but for the inspected experiment (`experiment_Tolonen2006_ntcA_depletion_MED4`) both were null. Only one experiment checked so far.

**Impact on methodology.** "High-rank per experiment" may not be able to rely on the stored rank fields uniformly. Likely fallback: compute rank per experiment from significant up-regulation + |log2FC|. To be confirmed across experiments in step 2.

## 2026-05-30 — Some KG genes have null protein sequences

**What happened.** `gene_aa_sequence` returned `not_matched` for PMM0838 — the gene exists but its `sequence` is null (consistent with the ~3% expression-only genes noted in `docs://guide/concepts`).

**Impact on methodology.** Seed genes must be filtered to those with a usable identifier for PaperBLAST. Sequence-null genes are flagged/dropped from the gene list (or rely on `protein_id` lookup if available). (Step 2 update: sequences dropped entirely per researcher — locus_tag / RefSeq `protein_id` is the PaperBLAST input, so KG sequence nulls no longer gate seed eligibility.)

## 2026-05-30 — `differential_expression_by_gene` raises IndexError on 42/197 experiments

**What happened.** Calling `differential_expression_by_gene(experiment_ids=[eid], ...)` raises `IndexError: list index out of range` inside the package (`functions.py` `_validate_organism_inputs` → `conn.execute_query(diag_cypher)[0]`) for 42 of 197 experiments. The failing set is stable across runs and spans EXOPROTEOMICS (3), PROTEOMICS (11), VESICLE_PROTEOMICS (10), METABOLOMICS (12), RNASEQ (3), plus a few others — i.e. mostly non-whole-cell-RNA experiments, several with no gene-level DE at all.

**Impact on methodology.** The all-experiments seed builder wraps each experiment in try/except, skips+logs the failures, and proceeds on the 155 that work. Acceptable for this analysis (the skipped experiments carry little/no gene DE), but it is an upstream package bug worth reporting to `multiomics_explorer`: the diagnostic-query indexing should handle experiments whose DE table is empty/absent rather than raising.

## 2026-05-30 — list_experiments field is `experiment_id`, not `id`; DE fields are `log2fc`/`rank_up`

**What happened.** Early seed-builder attempts used `e["id"]` and `r["log2_fold_change"]`/`r["rank_by_effect"]` (from the `kg_schema` edge property names). The typed-tool results actually expose `experiment_id`, `log2fc`, `rank_up`, `rank_down`, `expression_status` (per `docs://tools/differential_expression_by_gene`). The schema property names ≠ the tool's result field names.

**Impact on methodology.** Wasted a couple of iterations producing empty/zero pools. Lesson logged for the methodology: read the per-tool doc's "Per-result fields" before scripting against result dicts; don't assume `kg_schema` property names carry through to tool output.

## 2026-05-31 — PaperBLAST save: filename collision + one missave

**What happened.** Browser "Save as" defaults to a generic `PaperBLAST.html`, so saving multiple queries overwrites. One pro seed (PMT1742, rplN) was saved as a byte-identical copy of another (PMT9312_0549, csoS1A) — the inventory script caught it via md5/query-line, and it was re-saved correctly.

**Impact on methodology.** The save step needs one file per seed named with the locus_tag (or protein_id). `01_inventory.py` is the guard — it matches each file to its seed by in-page query line (not filename) and flags missing/duplicate/unmatched, so a missave can't silently corrupt the results.

## 2026-05-31 — relevance filter: "same phylum" too loose; Cyanorak ≠ cyanobacteria

**What happened.** First relevance pass kept any Proteobacteria hit (the KG heterotrophs' phylum), which readmitted clinical/gut bacteria (P. aeruginosa, E. coli, Salmonella) for cyano seeds. Also, a genus-token match would have wrongly flagged freshwater Synechococcus elongatus / PCC 7002 as "in Cyanorak".

**Impact on methodology.** Two fixes: (1) the relevance filter is **pool-aware** — cyano-only for pro/syn, Proteobacteria for the heterotroph pools; (2) `in_cyanorak` is computed from the **Cyanorak strain tables** (strain-code match), not a genus match, because Cyanorak is marine-picocyano only (excludes Synechocystis, Anabaena, freshwater elongatus). The Cyanorak CSVs live in `multiomics_biocypher_kg/data`.

## 2026-05-31 — PaperBLAST links papers by PMCID; NCBI id-converter for DOI dedup

**What happened.** PaperBLAST result pages link papers almost entirely by PMC ID (102 PMCIDs vs 3 explicit DOIs for the pro pool). The KG dedup target is keyed by DOI, so PMCID→DOI resolution is required.

**Impact on methodology.** Added an NCBI id-converter step (cached on disk, dedup-before-call so each paper sends one id, exponential backoff on HTTP 429). NCBI throttles to 3 req/s anonymously, 10 with an API key (read from `.env` via python-dotenv). Without resolution the KG dedup silently passes everything as "new" — so this is load-bearing, not cosmetic.

## 2026-06-03 — heterotroph relevance: KG-genus + identity floor; Pseudomonas off-domain

**What happened.** For alt/other_hetero seeds the "keep all Proteobacteria (same phylum)" rule readmitted clinical/plant-pathogen Proteobacteria, ballooning alt to 919 candidates. Even restricting to KG genera, P. putida (a KG genus) hits are dominated by clinical P. aeruginosa / agricultural P. syringae literature — right genus, wrong ecology. A genus filter alone can't separate these.

**Impact on methodology.** The heterotroph keep-rule is now **KG genera only + identity floor** (`--min-identity-hetero`, default 50%): a close homolog is more likely the actual marine relative; distant clinical paralogs fall below the floor. Cut alt 919→36. But the residual is Pseudomonas-specific and inherent to seeding from P. putida genes — flagged for the researcher rather than auto-filtered (a topic/ecology classifier, not just taxonomy, would be needed to catch it). Genuinely marine alt yield is ~6 (Shewanella oneidensis, Ruegeria pomeroyi).

## 2026-06-03 — protein_id version-suffix mismatch in file→seed matching

**What happened.** Some saves were named by protein_id without the version suffix (`WP_014948722` vs seed `WP_014948722.1`), and PaperBLAST's in-page query echoes a UniProt/UniParc accession, not the locus_tag — so 4 alt files came back unmatched in the inventory.

**Impact on methodology.** Both matchers (`01_inventory.py`, `02_parse_hits.py` `query_id`) now also try the version-stripped protein_id. A side observation: PaperBLAST may report "No hits to characterized proteins" for a query yet still list papers attached to *uncharacterized* homologs — those are kept (they can still be relevant literature).
