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

## F4 — `"N/A"` stored as literal string in cluster description fields (not NULL)

**Date:** 2026-04-29 (encountered in step 4 dossier rendering and confirmed in step 5 sweep)

**What happened.** `gene_clusters_by_gene(verbose=True)` returns the literal string `"N/A"` (case-sensitive, with no whitespace variation observed) for `cluster_functional_description`, `cluster_expression_dynamics`, and `cluster_temporal_pattern` when those fields are not curated for a given cluster. The KG does not return `null` / `None` / `""` / SQL `NULL` for these — it returns the three-character string `"N/A"`. Consumers that filter by `pd.isna(...)` or `field is None` will treat `"N/A"` as populated content and surface uninformative cells in downstream tables.

Concrete example: `cluster:1471-2180-14-11:med4_expression_level:VEG` (the cluster touching 72 / 116 candidates) returns `cluster_functional_description="N/A"`, `cluster_expression_dynamics="N/A"`, `cluster_temporal_pattern="N/A"`. The dossier rendering code in `4_methods/dossier.py::_fmt_na` strips the literal `"N/A"` (and `None`, `""`) to em-dash for display, but the JSON dump preserves the original string.

**Workaround.** `dossier.py::_fmt_na` handles rendering; downstream aggregations on `agg_clusters_pivot.csv` use the raw string. Future scripts working off `dossiers.json` need to apply the same convention — `value is None or value == "" or str(value).strip().upper() == "N/A"` is the test for "no content".

**Downstream impact.**
- *Step 5 aggregations (already produced):* `agg_clusters_pivot.csv`'s `cluster_functional_description` column carries the literal `"N/A"` for un-curated clusters; the `agg_top_fc_signals.csv` "first informative cluster" column uses the strip-N/A rule and correctly returns empty for VEG-only candidates.
- *KG-enhancement proposal (step 6 consolidation).* Either (a) normalize at the KG layer — replace `"N/A"` with `NULL` / empty in the underlying graph database, or (b) document the convention in the KG schema docs so consumers default to the strip-N/A rule. Option (a) is the cleaner remediation; option (b) is a quicker workaround that requires every downstream tool to know the convention.

---

## F5 — `function_description` carries placeholder text `"Alternative locus ID"`

**Date:** 2026-04-29 (encountered in step 4 anchor cards and confirmed in step 5 sweep)

**What happened.** The `function_description` field on Gene nodes (returned by `gene_overview(verbose=True)`) sometimes contains real curated text (e.g., for ntcA: "Required for full expression of proteins subject to ammonium repression. Transcriptional activator of genes subject to nitrogen control") and sometimes contains placeholder text `"Alternative locus ID"` that is meaningless as a function description. Examples in the candidate set: PMM0958 → `"Alternative locus ID"`; PMM0689 (hli) → `"Alternative locus ID"`; PMM1813 → `"Alternative locus ID"`; PMM1404 (hli) → `null`.

Reading the dossier card top-to-bottom, a reader sees `function_description: Alternative locus ID` on PMM0958's identity section and momentarily mistakes it for a real description before realizing it's metadata stub. The card is honest (it shows the field as the KG returned it) but the field itself is polluted at the source.

**Workaround.** None at the dossier level. Step 5 cards display the field as-returned. The dossier card reader must mentally filter `"Alternative locus ID"` as a non-description.

**Downstream impact.**
- *Future analyses:* same risk for any per-gene tool that surfaces `function_description`.
- *KG-enhancement proposal (step 6 consolidation).* Clean up the `function_description` content at the source: either populate with real curated text where available or replace placeholder strings with `NULL`. The placeholder appears to be an artifact of an upstream import that filled the field with a row-identifier-like string when no real description was available — should be undone at the import layer.

---

## F6 — Curated cluster description coverage is uneven across clustering analyses

**Date:** 2026-04-29 (encountered in step 5 cluster aggregation)

**What happened.** From `3_analysis_framing/data/discover_clusters_aggregate.csv` (re-confirmed in step 5):

| analysis | candidate rows touching it | rows with curated `cluster_functional_description` | rate |
|---|---:|---:|---:|
| MED4 K-means N-starvation clusters (Tolonen 2006) | 37 | 37 | 100 % |
| MED4 phage-upregulated transcription groups (Lindell 2007) | 8 | 8 | 100 % |
| MED4 diel cycling expression clusters (Zinser 2009) | 80 | 58 | 72 % |
| MED4 gene expression level classification (Wang 2014) | 113 | 32 | 28 % |

The expression-level analysis (VEG / HEG / MEG / LEG / NEG / "--") is touched by every candidate that has any cluster membership (113 of 116) but only 28 % of those rows carry a curated functional description. The K-means N-starvation analysis is the smallest in candidate-coverage (37 rows) but has 100 % curated descriptions. So the analysis with the highest dossier impact (most candidate cards touched) is also the worst-curated; the analysis with the most informative content per row (K-means N-starvation, with curated descriptions like "Contains nitrogen transport genes such as urtA and cynA") touches the fewest candidates.

For 5 of the top 10 candidates by log2fc (PMM1828, PMM1813, PMM1966, PMM1898, PMM1939) the only cluster membership is in the under-curated expression-level analysis — they get an N/A description on their only cluster axis. This is the dossier's worst-served subset and is captured in step 5's S4 surprise.

**Workaround.** None within the analysis. The dossier honestly reports the curated description as `"N/A"` (after F4 normalization in rendering); the candidate's role-suggestion content from the cluster axis is just absent for those 81 rows in the expression-level analysis (113 − 32 = 81).

**Downstream impact.**
- *Step 5 conclusions:* the candidate set is *cluster-axis-informative for ~30 % of cluster rows* (37 + 8 + 58 + 32 = 135 of 238 total cluster rows = 57 %, but the share is much lower for VEG / HEG / MEG / LEG / NEG — see breakdown above). The dossier's "potential role" anchor on the cluster axis is most reliable for K-means N-starvation cluster members.
- *KG-enhancement proposal (step 6 consolidation).* The expression-level (VEG / HEG / MEG / LEG / NEG) clustering is a per-gene "expression-level bin" classification, not a functional grouping per se — so the lack of curated functional descriptions may be by design rather than oversight. If so, the KG could (a) flag this analysis as "non-functional clustering" so consumers know not to expect role-bearing descriptions, or (b) add a different metadata field that describes what each VEG / HEG / MEG / LEG / NEG bin SHOULD mean (e.g. "VEG = top RPKM quartile, expected to include constitutive housekeeping + highly-expressed regulators") rather than leaving the per-cluster description as `"N/A"`.

---

*Append further entries as encountered.*
