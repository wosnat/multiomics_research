# Gaps and friction — Proteome–transcriptome discordance, Weissberg 2025

Append-only log of methodology / KG / tooling friction encountered during this analysis.
Decisions live in each step's `notebook.md`; this file captures friction (gaps, schema mismatches, anti-hallucination corrections, process slowdowns).

---

## F1 — MED4 axenic RNA-seq carries no timepoint label in the KG

**Date:** 2026-04-27 (encountered in step 1)

**What happened.** The KG experiment node `10.1101/2025.11.24.690089_growth_state_pro99lown_nutrient_starvation_med4_rnaseq_axenic` has `is_time_course: false` and `timepoints: null`. Its `experimental_context` is "Axenic Prochlorococcus MED4 in PRO99-lowN at 24C under continuous light" with no day mentioned. The reported fold-change is treatment-vs-control without a per-timepoint breakdown — but the corresponding axenic proteomics experiment has three timepoint contrasts (d14, d31, d89). Without timepoint metadata on the RNA-seq side, automated pairing of RNA-seq and proteomics fold-changes by timepoint is not possible from the KG alone.

**Resolution at step 1 (from researcher).** Osnat clarified that the axenic RNA-seq corresponds to **day 14** — the first proteome timepoint. The later axenic proteome timepoints (d31, d89) had no matching RNA-seq because RNA could not be extracted from axenic MED4 cells at those points (cell death, RNA degradation, insufficient biomass — see F2). This timepoint pairing came from outside-the-KG knowledge.

**Workaround.** For this analysis, the timepoint label "d14" is asserted in step 2 selection scripts as a hard-coded mapping for that experiment. This is documented in the step-2 notebook and propagated to all downstream pairing.

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
