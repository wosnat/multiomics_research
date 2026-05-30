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

**Impact on methodology.** Seed genes must be filtered to those with a usable identifier for PaperBLAST. Sequence-null genes are flagged/dropped from the gene list (or rely on `protein_id` lookup if available).
