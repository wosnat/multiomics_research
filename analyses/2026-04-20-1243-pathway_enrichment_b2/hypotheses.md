# Hypotheses log — pathway enrichment B2

Preliminary observations generated during the analysis that require downstream confirmation. Each entry is tagged with its confirmation state. **Do not cite as findings until state = CONFIRMED.** This file exists to prevent confirmation bias during write-up: observations surfaced at early-step QC can drift into conclusions if not bounded explicitly.

Entries are ordered by when the observation was made, not by importance.

---

## H1 — Coculture > axenic on N-limitation signature strength in Weissberg 2025 MED4

**Status:** PRELIMINARY (generated at Step 2 QC, 2026-04-20)
**Confirmation requires:** Step 4 scoring with proper NC calibration + LOO stability (Task 10). Specifically:
- Layer A `score_A(T_coculture, ontology)` ≥ Layer A `score_A(T_axenic, ontology)` across both ontologies, after matched-background calibration.
- Classification thresholds hold when LOO-R removes any single R experiment.
- Not attributable to the single T cluster that carries the observation (Weissberg coculture-Prot day 31) — the trend must hold across multiple T timepoints and both omics types.

**Observation (from Step 2 QC heatmap `step2_heatmap_{cyanorak_role,kegg}.png`):**

Coculture T clusters appear to show stronger N-limitation signature than axenic T clusters on multiple anchors:
- `cyanorak.role:E.4 N-metabolism`: red cells stronger in coculture block (days 18/31/60+89) than axenic block.
- `kegg.pathway:ko00910 N-metabolism`: same pattern — visible red across coculture timepoints, faint in axenic.
- `kegg.pathway:ko03010 Ribosome`: coculture-Prot day 31 shows the single loudest T hit (`***`, |s|=23.6, padj=2.3e-24).
- `cyanorak.role:K.2 Ribosomes`: coculture-Prot day 31 and day 60+89 strongest T blue cells.

**Why it's notable (would-be biology, if confirmed):**

The naive biological expectation is the opposite: MED4 dies axenically under N-limit in ~2 weeks; MED4 survives 90 days in coculture (per biology-context memory). One would predict **axenic** to show the strongest N-stress response (cells actively dying under N-deprivation). Instead the data appears to suggest:
- Axenic cells shut down generally without fully engaging the specific N-response program before dying.
- Coculture enables sustained transcriptional/proteomic engagement of N-scavenging pathways — MED4 *actively responds* to N-limit in coculture, which may be part of why it survives.

If this holds up through Step 4 scoring + calibration + LOO, it's a paper-worthy mechanistic observation about why coculture rescues N-limited survival.

**Risks (things that could falsify or weaken H1):**

1. **Timepoint mismatch.** Axenic and coculture T clusters aren't at matched timepoints. Weissberg axenic proteomics runs to day 14 (cells dying); coculture proteomics runs to days 31/60/89 (cells surviving). Comparing "dying axenic at day 14" to "surviving coculture at day 31" may compare different stages of a dynamic response rather than different responses per se.
2. **NC calibration artifact.** Decision #4 excludes Weissberg coculture day 11 up (NC) from N-metab calibration. If the same biology (coculture-enabled N-scavenging) also drives the coculture T signal, we may be calibrating against a noise floor that's artificially low BECAUSE we removed the biologically-relevant NC. Need to check whether coculture-T signal is above `nc_mean + 2σ` in both with-exclusion and without-exclusion calibrations.
3. **Omics bias.** Weissberg coculture has day 31+60+89 PROTEOMICS but not RNA-seq at those late TPs. The strongest signals sit in proteomics. Is the coculture-vs-axenic contrast actually an early-RNA-vs-late-Prot contrast, or a biology-driven contrast?
4. **Small-n.** 4 T experiments × ~7 timepoints ≈ 28 T clusters. Drawing cross-condition conclusions from visual pattern alone is weak; formal scoring is required.

**Where to verify in Step 4:**

When reading `results/score_summary.csv` for T clusters, specifically check:
- Per-ontology axenic-Prot vs coculture-Prot score comparison (matched omics, paired by TP where possible).
- Per-ontology axenic-RNA vs coculture-RNA score comparison at shared TPs (day 14, 18, 31 if RNA runs to those).
- LOO-R results: does the coculture > axenic pattern hold when Tolonen is removed? When Read is removed?
- Manual inspection of `loo_r_experiments.csv` `classification_flip` column for coculture T clusters specifically.

**Promotion to CONFIRMED / FALSIFIED happens at Step 4 decide (Task 11)**, with commit entry in `decisions.md` (D6 or later).

---
