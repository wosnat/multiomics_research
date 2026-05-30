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
