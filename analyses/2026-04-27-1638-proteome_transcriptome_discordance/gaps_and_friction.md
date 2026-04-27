# Gaps and friction — Proteome–transcriptome discordance, Weissberg 2025

Append-only log of methodology / KG / tooling friction encountered during this analysis.
Decisions live in each step's `notebook.md`; this file captures friction (gaps, schema mismatches, anti-hallucination corrections, process slowdowns).

---

## F1 — MED4 axenic RNA-seq carries no timepoint label in the KG

**Date:** 2026-04-27 (encountered in step 1)

**What happened.** The KG experiment node `10.1101/2025.11.24.690089_growth_state_pro99lown_nutrient_starvation_med4_rnaseq_axenic` has `is_time_course: false` and `timepoints: null`. Its `experimental_context` is "Axenic Prochlorococcus MED4 in PRO99-lowN at 24C under continuous light" with no day mentioned. The reported fold-change is treatment-vs-control without a per-timepoint breakdown — but the corresponding axenic proteomics experiment has three timepoint contrasts (d14, d31, d89). Without timepoint metadata on the RNA-seq side, automated pairing of RNA-seq and proteomics fold-changes by timepoint is not possible from the KG alone.

**Resolution at step 1 (from researcher).** Osnat clarified that the axenic RNA-seq corresponds to **day 14** — the first proteome timepoint. The later axenic proteome timepoints (d31, d89) had no matching RNA-seq because RNA could not be extracted from axenic MED4 cells at those points (cell death, RNA degradation, insufficient biomass — see F2). This timepoint pairing came from outside-the-KG knowledge.

**Workaround (closed in step 2, 2026-04-27).** The mapping is hardcoded in `2_kg_selection/scripts/02_build_paired_observations.py` as constants `F1_TIMEPOINT_HOURS = 336.0`, `F1_TIMEPOINT_LABEL = "day 14"`, `F1_GROWTH_PHASE = "nutrient_limited"`, `F1_TP_GENE_COUNT = 1849`. The script log records the row state before and after the patch for auditability. All downstream pair-construction goes through this fixup.

**Downstream impact.**
- Anyone re-running this analysis from the KG alone will not be able to pair MED4 axenic RNA-seq to MED4 axenic proteomics by timepoint — they need this gap-and-friction entry to know the pairing is d14.
- *Suggested KG fix:* add a `time_point_hours` / `time_point_label` value to the MED4 axenic RNA-seq experiment node so the pairing becomes machine-derivable.

---

## F2 — "No extractable RNA" ≠ "mRNA biologically absent"

**Date:** 2026-04-27 (encountered in step 1; anti-hallucination correction)

**What happened.** During step-1 dialogue I sketched an analysis sub-question framed as "protein persists while mRNA is gone" for the MED4 axenic d31 / d89 proteome-only data. Osnat corrected: the reason RNA was not extractable (cell death, RNA degradation, insufficient biomass for extraction protocol, technical extraction failure) is **not equivalent to** mRNA being biologically absent in the cell. Treating the missing RNA-seq data as "mRNA absent" would have produced a false post-transcriptional-persistence claim.

**Workaround.** The axenic MED4 d31 and d89 proteome-only data is excluded from the locked discordance analysis. It is treated as missing data (no information at all), not as "below biological detection." Any future use of this data must explicitly handle the missing-not-at-random structure.

**Downstream impact.**
- Step 1 scope locked to matched-timepoint paired RNA-seq + proteomics only (option (a) in the Q5 dialogue). The "extreme discordance" sub-narrative was abandoned.
- General methodological rule retained in user-level memory (`feedback_measurement_failure_vs_biology.md`): treat measurement failure as missing data, not as zero.
- *Suggested methodology improvement:* the research-methodology skill could explicitly call out this pattern (measurement failure vs biological absence) under [anti-hallucination.md](../../skills/research-methodology/references/anti-hallucination.md) or [statistical-rigor.md](../../skills/research-methodology/references/statistical-rigor.md), since it is a generic discordance / multi-omics failure mode.

---

## F3 — HOT1A3 has 5 `istB` paralogs detected by proteomics but not RNA-seq

**Date:** 2026-04-27 (encountered in step 2 gene-coverage QC)

**What happened.** The `qc_gene_coverage.py` per-(organism × condition) overlap report flagged 5 HOT1A3 locus tags present in proteomics but absent from RNA-seq: `ACZ81_17130`, `ACZ81_09590`, `ACZ81_07515`, `ACZ81_04425`, `ACZ81_01340`. A follow-up Cypher query showed all five are annotated as **IS21-like element helper ATPase IstB** (`gene_name = istB`). The pattern repeated identically in both axenic and coculture conditions.

**Interpretation.** This is a known methodological artefact from paralog-resolution differences between the two omics platforms. MS-based proteomics assigns shared peptides to multiple genome-redundant copies — when several IS-element copies have identical or near-identical IstB protein sequences, the peptide evidence cannot distinguish between them and the pipeline reports all 5. RNA-seq with short-read multi-mapper filtering (the default in DESeq2 pipelines) discards reads that map ambiguously to multiple loci, so genome-redundant transposase paralogs are effectively dropped from the count matrix.

**Workaround.** The 5 paralogs are excluded from the paired-discordance pool naturally by the inner-join in `02_build_paired_observations.py` (no protein FC × RNA FC pair can be formed without both). HOT1A3 paired pool is therefore 2221, not the 2226 distinct proteomics genes. Documented in `2_kg_selection/notebook.md` Surprise S1.

**Downstream impact.**
- *Step 4 / 5:* the analysis pool for HOT1A3 is 2221 paired genes; this is the universe over which discordance scores can be computed. Repeated explicitly in `paper.md` Background and Methods.
- *General lesson:* any paired-omics analysis on bacterial genomes with expanded transposase / IS-element families will encounter this pattern. It is not a KG bug — the KG correctly reports what each upstream pipeline detected — but it is worth noting in framing.
- *Suggested KG addition:* a "paralog-shared-peptide" flag on Polypeptide nodes, derivable from sequence-identity clustering of the proteome, would let downstream consumers detect these cases without needing to inspect annotations.

---

## F4 — Initial `qc_gene_coverage.py` reached for `run_cypher` when the API covered the case

**Date:** 2026-04-27 (encountered in step 2)

**What happened.** The first version of `2_kg_selection/scripts/qc_gene_coverage.py` used `run_cypher` to fetch the distinct locus_tag set per experiment via `MATCH (e:Experiment)-[:Changes_expression_of]->(g:Gene) RETURN DISTINCT g.locus_tag`. Osnat asked whether this was a gap in the API. It was not — `differential_expression_by_gene(experiment_ids=[exp_id], limit=None)` accepts `experiment_ids` without `locus_tags` (both are `None`-default optional filters) and returns every DE row for the experiment, from which the locus_tag set is one set-comprehension away. The envelope additionally exposes `matching_genes`, which serves as a built-in cross-check against the comprehension count.

**Workaround / fix.** Script switched to `differential_expression_by_gene(experiment_ids=[exp_id], limit=None)`. Re-running produced identical numbers (1849 / 1424 / 1424 for MED4 in both conditions; 3950 / 2226 / 2221 for HOT1A3 in both conditions). The comprehension count agreed with `matching_genes` for all 8 experiments; no warnings emitted.

**Generalizable lesson.** Function names in `multiomics_explorer` reflect the most common use case, not the only one. `differential_expression_by_gene` looks gene-centric but doubles as an experiment-scoped DE-row dump when `experiment_ids` is supplied alone. Before reaching for `run_cypher`, inspect the signatures of the existing functions (`inspect.signature`, or the import list in [python-api-guide.md](../../skills/research-methodology/references/python-api-guide.md)) — broader applicability is hidden behind `None`-default filter parameters.

**Downstream impact.**
- Methodology-side: a user-level feedback memory was added (`feedback_check_api_before_cypher.md`) so future analyses default to "check the API surface first."
- *Suggested methodology improvement:* [python-api-guide.md](../../skills/research-methodology/references/python-api-guide.md) could note explicitly that several functions named for one entity (e.g., `_by_gene`, `_by_ortholog`, `_by_function`) accept filters for adjacent entities with `None` defaults and double as broader queries — the current "Writing raw Neo4j or `requests` calls → Use the API functions" rule covers the spirit but doesn't surface this specific pattern.

---

## F5 — Editorializing data observations at step 2 before the formal analysis runs

**Date:** 2026-04-27 (encountered in step 2 close)

**What happened.** While presenting the QC #4 canonical-marker heatmap (5 N-regulon genes showing mRNA-down / protein-up in MED4 coculture), I described it as "biologically loaded", "massive finding", "biologically explosive", and claimed the data "reframes the analysis" with a "strongly motivated default hypothesis." Osnat asked me to flag every place where my interpretation read as hand-wavy. Two failures: (a) using emotive vocabulary for data observed but not yet formally tested — discordance metric is step-4 work, axis-stratified testing is step-5; (b) pre-committing to a post-transcriptional-regulation interpretation when the same data is also compatible with mRNA-degradation artefacts at coculture sampling, normalization issues, paralog mapping, etc. The 6-step methodology puts hypothesis testing at step 5 and caveat harvesting at step 6 to prevent step-2 enthusiasm from biasing later analysis. Also flagged: my use of "confirms the prior-analysis flag" for atpA, which leaned on a prior conclusion about the same single data point as authority rather than describing the data here.

**Workaround.** Notebook + paper text describing the marker heatmap and the agreement-overview Pearson values has been re-written in factual terms ("5 N-regulon markers show RNA log2FC < 0 and protein log2FC > 0 at MED4 coculture days 31–89", "HOT1A3 paired-pool Pearson r is in [−0.09, +0.03] at every observation"). Interpretive framings tagged `[interpretation]` and listed alongside competing alternatives where applicable.

**Generalizable lesson.** Vocabulary itself signals premature commitment. Emotive words like "striking", "massive", "rich", "central", "explosive", "reframes" should not appear in step-2 / step-3 outputs. Reserve them for step-6 evaluation, where preregistered predictions have been tested and the framing has earned them.

**Downstream impact.**
- A user-level feedback memory was added (`feedback_no_pre_celebration.md`) so future analyses default to factual-at-early-steps phrasing.
- *Suggested methodology improvement:* the [research-methodology/references/anti-hallucination.md] could add a category covering "interpretive vocabulary as premature commitment" — distinct from the existing categories that focus on hallucinated facts. The failure mode here was real-data + over-claimed interpretation, not invented data.
