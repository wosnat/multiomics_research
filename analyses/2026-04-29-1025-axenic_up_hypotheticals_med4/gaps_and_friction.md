# Gaps and friction — Axenic-upregulated hypotheticals, MED4 (Weissberg 2025)

Append-only log of methodology / KG / tooling friction encountered during this analysis.
Decisions live in each step's `notebook.md`; this file captures friction (gaps, schema mismatches, anti-hallucination corrections, process slowdowns).

This analysis carries an **active KG-enhancement sub-question**: as the per-gene dossier is built, log every place where the KG fails to answer a question we want to ask, and what addition to the KG would have helped. Step 6 evaluation consolidates the entries into a prioritized list of KG-enhancement proposals. Concrete leads already on the table at step 1 (researcher-supplied, not yet hit in the analysis): amino-acid sequences absent due to per-node size limits; batch-bioinformatics layers such as Pfam / InterPro domain hits, SignalP, TMHMM, AlphaFold structure summaries. These get F-entries when the analysis actually attempts the corresponding question and is blocked.

---

## F1 — 17 / 116 candidate genes have no ontology annotations in the KG

**Date:** 2026-04-29 (encountered in step 2, candidate-set extraction)

**What happened.** Of the 116 candidate genes (sig_up × annotation_quality ≤ 1 in MED4 axenic RNA-seq), `gene_overview` returns an empty/NaN `annotation_types` for 17 (15 %). These genes have no GO, no KEGG, no EC, no COG, no Pfam, no Cyanorak, no TIGR, no BRITE term in the KG. Two of the 17 are in the top-10 by log2fc (PMM1898 log2fc=4.68 AQ=0, PMM1939 log2fc=4.65 AQ=0). They have orthologs (mostly) and clusters but the per-gene "what IS this gene" question is unanswerable from the KG content alone.

**Workaround.** None at step 2 — the candidate set keeps these 17 as full participants, and the dossier at step 4/5 will record them with explicit "no ontology" rows. They cluster in their orthologs and on the response-profile axis; that's all the KG can offer them.

**Downstream impact.**
- *Step 5 (dossier):* a "no ontology" subset will surface as a distinct slice. Its inferential ceiling is much lower than the rest of the candidate set.
- *KG-enhancement proposal (step-1 sub-question lead now activated).* The same set of 17 genes is exactly where the researcher-supplied step-1 leads would help most: (a) **AA sequences in the KG** would let us run ad-hoc Pfam / InterPro / structure searches and recover something about these genes from sequence alone; (b) **batch-bioinformatics layers** — Pfam, InterPro, SignalP, TMHMM, AlphaFold structure summaries — pre-computed and stored as ontology rows on Gene nodes — would close most of this annotation gap without us needing to run scans per analysis. Adding either or both of these to the KG would change "17 candidate genes have nothing" to "17 candidate genes have at least domain / topology / structural priors". This is the first concrete instance of the step-1 KG-enhancement sub-question being hit by real analysis content; it is consolidated at step 6 alongside any further entries.

---

## F2 — Source-level "hypclass-only" classification was over-reaching; retracted

**Date:** 2026-04-29 (encountered in step 2; first instance of an anti-hallucination correction during this analysis)

**What happened.** While reporting ontology coverage on the candidate set I introduced a 3-bucket per-gene classification: "informative" = at least one term from `pfam` / `kegg` / `go_*` / `brite` / `ec`; "hypclass-only" = terms ONLY from `cog_category` / `cyanorak_role` / `tigr_role`; "none" = no ontology data. This produced the headline "61 / 116 are hypclass-only". Osnat flagged that the assumption "all cog/cyanorak/tigr terms are hypothetical-class" was unverified.

Empirical check via `gene_ontology_terms(locus_tags=[…116 candidates…], organism="MED4")`: of the 309 candidate-set ontology rows, the cog/cyanorak/tigr sources do contribute mostly catch-all terms (TIGR "Hypothetical proteins / Conserved" 84/90, COG "Function unknown" 83/88, Cyanorak "Other > Conserved hypothetical proteins" 57/71), but each source also carries a small tail of informative terms in this candidate set — e.g., COG "Coenzyme transport and metabolism" (2 genes), "Cell wall/membrane/envelope biogenesis" (1), "Translation, ribosomal structure and biogenesis" (1), "Post-translational modification, protein turnover, chaperones" (1); Cyanorak "Cellular processes > Adaptation/acclimation > Phosphorus / Iron / Other" (~7 genes), "Cell envelope > Surface structures" (1), "Transcription > Other" (1); TIGR "Transcription / Other" (1), "Cell envelope / Surface structures" (1). About 16 candidate genes carry an informative cog/cyanorak/tigr term, so "61 hypclass-only" was wrong as a per-gene classification.

**Workaround / retraction.** The per-gene 3-bucket classification has been removed from the analysis. Step 2 now publishes (a) descriptive per-source presence flags (no claim about content), and (b) the term-level landscape from the `gene_ontology_terms` envelope (`by_term` and `by_ontology`). Per-gene informativeness is deferred to step 4/5 dossier construction, where each candidate's actual term names are inspected. `qc_ontology_per_source_flags.csv` (descriptive flags) replaces the retracted `qc_ontology_buckets.csv`; the latter file's classifications should not be trusted.

**Generalizable lesson.** Ontology informativeness is a property of the *term*, not of the *source*. Future hypothetical-gene analyses should classify at the term level — either by inspecting term names or by maintaining an explicit list of catch-all terms (e.g. "Function unknown", "Conserved hypothetical proteins", "Hypothetical proteins / Conserved", "uncharacterized protein", "Domain of unknown function (DUF*)") and treating those as uninformative regardless of source. Source-level proxies are tempting because they are one envelope field away, but they encode an assumption that doesn't hold uniformly.

**Downstream impact.**
- *Step 4/5 (dossier):* term-level evaluation is required. The `data/qc_candidate_ontology_terms_by_term.csv` artifact gives the full term list ranked by gene-frequency for the candidate set — a starting point for deciding which terms count as informative.
- *KG-enhancement proposal (still relevant for F1).* The original motivation — that AA sequences and batch bioinformatics layers (Pfam / InterPro / SignalP / TMHMM / AlphaFold) would close the annotation gap on the most under-annotated candidates — still holds for the 17 no-ontology genes (F1) and for the broader subset of candidates carrying only catch-all terms in the term-level analysis. The remediation is the same; only the per-gene count is now uncertain until step 4/5 confirms it.

---

## F3 — 2 RefSeq-only locus tags carry no cluster, no homolog, no ontology

**Date:** 2026-04-29 (encountered in step 3, dossier-axis discovery)

**What happened.** The 116-gene candidate set contains two locus tags with the RefSeq-style prefix `TX50_RS` rather than the expected MED4 `PMM` prefix: **`TX50_RS09500`** (log2fc 2.70, AQ=0, "hypothetical protein") and **`TX50_RS09520`** (log2fc 1.93, AQ=0, "hypothetical protein"). Per-axis discovery returned: no homolog groups (`gene_homologs` reports both in `no_groups`), no cluster memberships (`gene_clusters_by_gene` reports both in `not_matched`), no ontology terms (already counted in the F1 17-gene no-ontology subset). The dossier card on these two genes reduces to identity + DE evidence + cross-study response profile only; cluster, homolog, and ontology axes all surface as "no data" rows.

A third RefSeq-only candidate **`TX50_RS09860`** (log2fc 1.75, AQ=0, "conserved hypothetical protein") has no clusters and no ontology but DOES have homolog groups — partially less stripped than the two above.

**Workaround.** No workaround at step 3 — the dossier card honestly reports "no data" on the affected axes. The two genes pass through to step 4/5 dossier construction with most-stripped floor-case behavior. Step 4 method-design will exercise the empty-axis branches against PMM1898 (driving example B, F1 floor with singleton homolog group) and may also test `TX50_RS09500` as an even-more-stripped check (no cluster + no homolog + no ontology + minimal locus identity).

**Interpretation.** [interpretation] The `TX50_RS` prefix is the RefSeq automated-locus-tagging system. PMM-prefix locus tags are the curated MED4 genome annotation. The two `TX50_RS09500/09520` genes are probably late additions to the genome (RefSeq-only entries that never received PMM curation, possibly small ORFs or re-annotated regions) and consequently were never folded into the cyanorak / eggNOG ortholog grouping pipelines or any of the four clustering analyses. The TX50 RefSeq tags don't show up in the legacy MED4 annotation pipelines because those pipelines were built on the PMM identifier set.

**Downstream impact.**
- *Step 4/5 (dossier):* the two RefSeq-only genes define the dossier's true floor. The card layout's empty-row branches must handle them gracefully — three "no data" rows in addition to the F1 ontology-empty row.
- *KG-enhancement proposal (step 6):* the lack of ortholog grouping and cluster membership for these two genes traces back to the **upstream pipeline scope** (cyanorak curation, eggNOG inference, microarray-era clustering analyses) — these were all run on the PMM-identifier set without reprocessing for late-added RefSeq entries. The remediation is at the pipeline boundary, not the KG schema: re-run the ortholog-inference and cluster-membership pipelines with the current full locus-tag set (PMM ∪ TX50_RS) so late-added genes get assigned to groups and clusters where appropriate. This is a different remediation from the F1 / earlier KG-enhancement leads (AA sequences + batch bioinformatics layers), which addressed annotation-content gaps; F3 addresses identifier-coverage gaps in the precomputed analyses.

---

*Append further entries as encountered.*
