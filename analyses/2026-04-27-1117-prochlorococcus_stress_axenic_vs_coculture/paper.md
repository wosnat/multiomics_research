# Is Prochlorococcus stressed in axenic vs co-culture conditions? (Weissberg 2025)

## Question

Within each of axenic and coculture *Prochlorococcus* MED4 cultures from Weissberg et al. 2025 [1], how do five stress axes — N-starvation, oxidative, proteotoxic, photosynthetic, and cell-death/late-stationary — evolve from an exponential-phase baseline through nutrient-limited and death-phase timepoints, scored independently in transcriptome and proteome?

The question is structured as two parallel within-condition trajectory analyses (axenic, coculture). Cross-condition statistical contrast at matched timepoints is not included in this analysis (data limitation; see `gaps_and_friction.md`). Cross-condition reading is performed by inspection of the curve shapes only.

The analysis is scoped to a single publication (Weissberg 2025) and a single organism (*Prochlorococcus* MED4 — the only *Prochlorococcus* strain in the study). *Alteromonas macleodii* HOT1A3 data, although available in the KG for this study, is out of scope: the question concerns *Prochlorococcus* alone.

## Background

The data underlying this analysis comes from Weissberg et al. 2025 [1], a 90-day multi-omics time-course comparing *Prochlorococcus* MED4 in axenic culture vs co-culture with *Alteromonas macleodii* HOT1A3, both grown in PRO99-lowN medium under continuous light at 24 °C. The published narrative is that axenic *Prochlorococcus* perishes under extreme N-deprivation while coculture cultures persist for 90+ days, supported by N-recycling activity attributed to *Alteromonas*.

The KG indexes 10 differential-expression experiments for this publication: 5 with *Prochlorococcus* MED4 as the profiled organism, 5 with *Alteromonas* HOT1A3. This analysis uses the 5 MED4 experiments. They cover both omics platforms (RNA-seq, proteomics) and both conditions (axenic, coculture), and all use `table_scope == "all_detected_genes"`, supporting fair cross-experiment comparison. All 4 nutrient-starvation experiments are structured as `PRO99-lowN nutrient starvation` (treatment) vs `PRO99-lowN exponential growth` (control) within the same condition — i.e., each timepoint's reported log2FC is already against that condition's own exponential baseline, which matches the within-condition trajectory framing locked at step 1.

The TP / phase structure has two important asymmetries that constrain what this analysis can claim:

1. **Axenic RNA-seq is not a time-course.** It is a single contrast (nutrient_limited vs exponential), with no per-timepoint breakdown. The corresponding axenic RNA "trajectory" is therefore a single point.
2. **Calendar-shared proteomics timepoints are in different growth phases across conditions.** Day 31 and day 89 exist in both axenic and coculture proteomics, but axenic is in *death* phase by day 31 while coculture remains in *nutrient_limited* at all sampled days. Cross-condition statistical contrast at calendar-shared days would conflate "stress mitigation" with "different physiological states at the same calendar day."

Together these confirm and characterize the data limitation flagged at step 1 (`gaps_and_friction.md` F1): cross-condition statistical contrast is not supported. The analysis produces three multi-point within-condition trajectories (axenic proteomics — 3 timepoints; coculture proteomics — 4 timepoints; coculture RNA-seq — 4 timepoints) plus one single-point summary (axenic RNA-seq). Trajectory plots align by `timepoint_hours` and annotate growth phase per point.

The five stress axes (N-starvation, oxidative, proteotoxic, photosynthetic, cell-death/late-stationary) and the within-condition trajectory framing follow from step 1.

## Methods

### Framing (step 3)

**Hypothesis.** If *Prochlorococcus* MED4 is stressed in a given axis (N-starvation, oxidative, proteotoxic, photosynthetic, or late-stationary) within a given condition (axenic, coculture), then across the within-condition trajectory from exponential phase through nutrient-limited and (where present) death-phase timepoints, the panel of canonical positive-control genes for that axis will show coordinated, statistically significant differential expression in the direction empirically calibrated against this dataset (rather than a textbook-presumed direction).

**Control panel.** We curated a 26-gene control panel — 23 positive controls across the 5 stress axes plus 3 negative-control housekeeping candidates — anchored in the cyanorak hand-curated functional ontology where available (`D.1.3` Adaptation/Nitrogen, `D.1.4` Oxidative stress, `L.3` Protein folding and stabilization, `J.8` Photosystem II), supplemented by GO Biological Process terms and canonical literature markers (psbA, lrtA) where the ontologies were sparse for *Prochlorococcus*. The cell-death axis was reframed as **late-stationary / starvation response** because formal cell-death annotations are absent in both GO BP and cyanorak for MED4 (`gaps_and_friction.md` F4).

**Direction calibration via validation.** For every panel gene we pulled all available DE rows from the 4 trajectory + 1 single-point experiments (298 rows) and inspected each gene's log2-fold-change and significance pattern across (omics × condition × timepoint) cells. The validation revealed that **direction is axis-specific, not always upregulated**:

- N-stress markers up (ntcA, glnA, glnB, urtA, amt1)
- Photosystem proteins down (psbA, psbD, ftsH2 — PSII disassembly), but high-light-inducible proteins up (HLI; PMM1404 reaches log2FC = +9.6)
- Proteotoxic markers weakly up (htpG, groES) or flat (dnaK1)
- Oxidative markers mostly flat (sodN, ahpC, gor) — see explicit limitation below
- Late-stationary marker lrtA strongly down; spoT mixed; isiB down (wrong-handed for N-stress, retained as a transparency anchor)

**Sensitivity limits made explicit.** Three limits were calibrated by the validation:

1. *Oxidative axis* — the canonical *intracellular* oxidative-defense responders (sodN, ahpC, gor) are unresponsive in this experiment. *Prochlorococcus* MED4 lacks catalase and depends on heterotroph / extracellular peroxiredoxin activity for ROS detoxification (consistent with the Black Queen Hypothesis); the intracellular transcript / protein readout is therefore an insensitive assay. A "no oxidative-axis signal" finding from this analysis is read as *insensitivity to the relevant biology*, not as "no oxidative stress."
2. *HLI proteins are RNA-only* in the photo axis — bottom-up proteomics does not detect them due to small size. Cross-omics concordance for photo is computed only over psbA, psbD, ftsH2; HLI is RNA-side only.
3. *Negative controls are imperfect under prolonged starvation* — atpA crashes (max|log2FC| = 4.7) under axenic-RNA d18; rpoB is the only near-flat reference. We retain all three and report the imperfect flatness rather than swap to candidates that may also fail.

**Operational success criterion.** For each (axis × condition × omics) cell the step-5 output is: (a) a stress-score trajectory aggregating signed positive-control contributions; (b) a direction-aware verdict per cell; (c) a cross-omics concordance flag with the photo-axis caveat; and (d) the explicit oxidative-axis insensitivity caveat.

### Implementation (step 4)

The per-axis stress score is a direction-aware signed-Z statistic against the genome-wide log2-fold-change distribution at each (experiment × timepoint). For axis gene set G with calibrated per-gene direction d_g ∈ {+1, -1}:

$$\text{axis\_score} = \frac{\overline{d_g \cdot \log_2 FC_g}_{g \in G} - \overline{\log_2 FC_b}_{b \in \text{background}}}{\text{SD}\!\left(\log_2 FC_b\right)_{b \in \text{background}}}$$

where the background is all non-axis genes quantified at the same (experiment, timepoint), with NaN values excluded. Population SD (ddof=0) is used. The score is reported alongside its components — the raw axis-mean signed log2FC, the background mean and SD, and the count of axis-genes-with-data — because the score and the raw axis response can disagree under heavy-tailed background distributions and that disagreement is itself biological information (see Results / Discussion).

Direction calibration (the d_g values) comes from the step-3 validation — for n_stress all five positive controls are +1 (UP), for the photo axis HLI proteins are +1 and psbA/psbD/ftsH2 are -1 (PSII disassembly), and so on. The implementation requires direction be specified per gene (`dict[locus_tag, ±1]`) and raises `KeyError` on missing entries to prevent silent default-to-+1 bugs.

Verification: the formula was hand-computed for six toy cases (basic, all-negative-with-symmetry, bidirectional, missing-direction error, NaN handling, single-gene axis) and the implementation reproduces every case to floating-point precision. The single-line docstring example was corrected during this step from `sqrt(0.13/3) ≈ 0.208` (sample variance) to `sqrt(0.13/4) ≈ 0.180` (population variance, matching the implementation).

Driving-example application: N-stress axis × axenic-proteomics across 3 timepoints (day 14 nutrient_limited, day 31 death, day 89 death) using the 5 N-stress positives (ntcA, glnA, amt1, glnB, urtA) all with direction +1. The axis score is +1.71 at day 14, drops to +0.93 at day 31, and recovers slightly to +1.01 at day 89 — driven not by the axis genes (whose raw signed log2FC stays near +1.27) but by the background SD widening at death phase as the global proteome response intensifies. This is the methodology insight that motivates reporting axis_score and raw axis-mean signed log2FC as two complementary panels in step 5 trajectories.

## Results

We scored five stress axes — N-starvation, photosynthetic, proteotoxic, oxidative, late-stationary/cell-death — across all (axis × condition × omics × timepoint) cells defined by the four trajectory experiments and one single-point summary. Each cell was scored with the validated 3–5-gene **positive panel** (per-gene direction calibrated against the data) and, where available, with the broader **cyanorak hand-curated panel** (~21–31 genes per axis) as a sensitivity control. The two panels agree qualitatively across all engaged cells; the positive panel produces larger absolute scores because it concentrates strongly-responding genes, while the cyanorak panel is more conservative.

### N-starvation: cross-omics divergence in coculture, condition-specific drop in axenic

The N-starvation axis is the strongest signal in the analysis and reveals a clear cross-omics divergence in the coculture condition. At the protein level, the five validated N-stress positives (ntcA, glnA, amt1, glnB, urtA) are upregulated and stable across all 4 sampled TPs in coculture (axis_score +2.79 → +2.66 → +2.76 → +2.49 from day 18 to day 89). At the RNA level in the same condition, the same five genes are *downregulated* across all 4 TPs (axis_score −0.89 → −1.07 → −1.76 → −1.25; raw axis-mean signed log2FC −0.57 to −1.82). This is consistent with a "transcription off, protein still on" pattern: coculture cells reduce N-scavenging transcription — likely because Alteromonas-supplied N relieves the need to make more N-scavenging machinery — while the existing N-scavenging proteome remains engaged.

In axenic, the same five genes are upregulated at both omics (axenic-RNA single point at +1.45; axenic-proteomics +1.71 → +0.93 → +1.01 across nutrient_limited day 14 → death day 31 → death day 89). The drop in axenic-proteomics axis_score from day 14 to day 31 reflects not a fading axis response (raw axis-mean signed log2FC stays around +1.16-1.35 throughout) but a widening genome-wide background SD (0.77 → 1.32) at death phase: the N-stress axis no longer stands out because the entire proteome is in flux.

### Photosynthesis: PSII disassembly accelerates at axenic death phase

The photo axis is bidirectional by construction: PSII reaction-center proteins (psbA, psbD, ftsH2) carry direction = -1 (DOWN under stress = engaged); high-light-inducible proteins (HLI; PMM1404, PMM0064) carry direction = +1 (UP under stress = engaged). The largest single-cell axis_score in the analysis is photo × axenic-proteomics × day 31 = **+3.15** (death phase), up from +1.17 at day 14 (nutrient_limited) — a near-tripling of the distinctiveness score concomitant with the physiological transition to death phase. Coculture-proteomics also rises (+1.07 → +2.14 at day 31) but plateaus at lower magnitude (+1.07 to +1.15 at later TPs). At the RNA level, axenic shows +0.86 at the single timepoint (day 14); coculture-RNA stays around +0.3 throughout and dips slightly negative at day 89 (−0.76), the only cell where the photo axis is in the un-engaged direction — consistent with coculture cells preserving photosynthetic capacity even at 90 days.

### Cell-death/late-stationary: RNA-side engagement, protein-side fade

The cell-death axis (spoT, lrtA, isiB; lrtA and isiB calibrated direction = -1 because they go DOWN under stress in this dataset) is engaged at the RNA level — axenic-RNA single-point score +1.51 is the strongest RNA-side signal of any axis × condition cell in the analysis. Coculture-RNA shows a milder but consistent positive signal across TPs (+0.83 to +1.24). In contrast, both proteomes show declining cell-death axis scores over time (axenic +0.87 → +0.55 → +0.48; coculture +1.10 → +0.48 across the 4 TPs). This is the only axis where RNA and protein read in opposite directions, plausibly because cell-death-marker proteins themselves degrade (or fall below MS detection) in cells that are dying.

### Oxidative and proteotoxic: weakly engaged across all cells

Both axes are weakly engaged across all 4 trajectory experiments × all TPs. Oxidative axis_scores stay between +0.19 and +0.63; no cell crosses any reasonable engagement threshold. Proteotoxic axis_scores stay between -0.13 and +1.25; the maximum is in axenic-RNA (single point), and even there the score is below z = +2. This is the calibrated sensitivity-limit prediction from step 3: *Prochlorococcus* MED4 lacks catalase and depends on heterotroph / extracellular ROS handling for oxidative stress relief (Black Queen Hypothesis), so the *intracellular* transcriptional / translational oxidative-defense response is a poor assay for "is *Prochlorococcus* oxidatively stressed in this dataset" regardless of condition. The proteotoxic axis is similarly muted — N-starvation does not appear to produce strong canonical chaperone-and-protease induction in this organism under these conditions.

### Cross-omics concordance summary

| Axis | Axenic (RNA vs Prot) | Coculture (RNA vs Prot) |
|---|---|---|
| N-starvation | concordant UP | **divergent** — RNA DOWN, Prot UP |
| Photosynthesis | concordant UP at single-point and rising-prot | RNA flat / dipping negative late, Prot UP |
| Cell-death | RNA strong UP, Prot declining | RNA modest UP, Prot declining |
| Proteotoxic | concordant low | concordant low |
| Oxidative | concordant null | concordant null |

The most informative cross-omics signal is the N-starvation × coculture divergence, which is precisely the cell where the multi-omics framing of the underlying paper is most relevant.

### Headline figure

The headline figure is [5_analyze/figures/trajectories_positive_panel.png](5_analyze/figures/trajectories_positive_panel.png) — a 5-axis × 2-omics grid showing axenic vs coculture trajectories per (axis, omics) cell, with axis_score on the y-axis and timepoint (days) on the x-axis. Companion figures: [trajectories_positive_panel_raw.png](5_analyze/figures/trajectories_positive_panel_raw.png) (axis-mean signed log2FC view), [trajectories_cyanorak_panel.png](5_analyze/figures/trajectories_cyanorak_panel.png) (broader gene-set sensitivity check), [panel_comparison.png](5_analyze/figures/panel_comparison.png) (positive-panel vs cyanorak-panel scatter).

The step-6 trajectory figure with permutation-significance overlaid: [6_evaluate/figures/trajectories_positive_panel_sig.png](6_evaluate/figures/trajectories_positive_panel_sig.png) (positive panel) and [6_evaluate/figures/trajectories_cyanorak_panel_sig.png](6_evaluate/figures/trajectories_cyanorak_panel_sig.png) (cyanorak panel).

The step-3 control-panel validation heatmap that grounded the per-gene direction calibration: [3_analysis_framing/figures/control_validation_heatmap.png](3_analysis_framing/figures/control_validation_heatmap.png).

The step-4 driving-example trajectory (N-stress × axenic-proteomics, 2-panel raw + Z): [4_methods/figures/n_stress_axenic_prot_trajectory.png](4_methods/figures/n_stress_axenic_prot_trajectory.png).

## Discussion

### Answering the locked question

*Is Prochlorococcus MED4 stressed in axenic and in coculture conditions under prolonged nitrogen deprivation, as measured across five stress axes?*

**Yes for both conditions, but with strikingly different multi-omics shapes.**

Within axenic culture, the photosynthetic axis emerges as the dominant signal at the physiological transition from nutrient_limited to death phase: the photo axis_score in axenic-proteomics jumps from +1.17 at day 14 to +3.15 at day 31 (BH-FDR p_bh = 0.0007), the largest single-cell engagement in the analysis. This is driven by significant downregulation of PSII reaction-center proteins (psbA, psbD, ftsH2) — the molecular signature of photosystem disassembly. Concurrently the N-starvation axis is engaged at both omics throughout the sampled timepoints (axis_score +0.93 to +1.71 in proteomics, p_bh 0.003-0.041; +1.45 at the single RNA timepoint, p_bh = 0.009). The cell-death/late-stationary axis is engaged at the RNA level (axis_score +1.51 at day 14, p_bh = 0.029) but proteomic cell-death scores decline over time, plausibly because the cell-death markers themselves degrade in failing cells.

Within coculture, the N-starvation axis tells a more nuanced story than axenic does. The five canonical N-stress positives (ntcA, glnA, amt1, glnB, urtA) are *strongly upregulated at the protein level* across all 4 sampled timepoints (axis_score +2.49 to +2.79; all p_bh ≈ 0.0007), the highest sustained engagement of any axis in any condition in the analysis. But the *same five genes* are *significantly downregulated at the RNA level* across the same timepoints (axis_score −1.07 to −1.76; p_bh 0.015-0.047 in the negative direction at 3 of 4 TPs). This RNA/protein decoupling is the most informative finding of the analysis: coculture cells reduce *transcription* of the N-scavenging machinery while maintaining the N-scavenging *proteome*. The most parsimonious interpretation is **homeostatic relief at the transcriptional level** — Alteromonas-supplied bioavailable N (per the published nitrogen-recycling mechanism) reduces the demand for *more* N-scavenging proteins while the existing proteome, synthesized under the early N-starvation stimulus, remains stable and functional. The protein-side axis is engaged because the cells maintain N-scavenging capacity; the RNA-side axis is *anti*-engaged because the cells no longer need to make more.

The photosynthetic axis is also engaged in coculture-proteomics (peak axis_score +2.14 at day 31, p_bh = 0.003), though at lower magnitude than axenic, and the coculture-RNA photo axis stays approximately null with one cell trending slightly negative at day 89 — coculture cells appear to preserve photosynthetic capacity even at 90 days. The proteotoxic axis is essentially null across all 12 cells (1 of 12 transient-significant); the oxidative axis is null in all 12 cells.

### How this maps onto the published nitrogen-recycling mechanism

The published narrative — that *Alteromonas* recycles organic N to bioavailable NH₄⁺ for *Prochlorococcus*, enabling 90+ day survival — predicts that coculture cells should be *less* N-stressed than axenic cells over time. Our analysis adds quantitative texture to that picture: coculture cells are *not* N-stress-free; rather, their N-scavenging proteome remains highly engaged (more so than axenic, in fact, because axenic distinctiveness drops at death phase as the entire proteome enters flux). The novel observation is that coculture cells appear to have *exited the transcriptional N-stress response* even while maintaining the corresponding proteome — consistent with the heterotroph providing *enough* N to relieve the transcriptional drive without enough N to dismantle the existing scavenging machinery. This is exactly the kind of "transcription off, protein on" pattern that multi-omics analysis is uniquely positioned to detect.

The accelerating PSII disassembly in axenic at death phase (photo axis_score doubling between day 14 and day 31) is consistent with the published axenic-collapse narrative and adds a quantitative timing anchor: photosystem disassembly does not increase smoothly with N-deprivation duration; it accelerates abruptly at the nutrient_limited → death transition. This may be a biomarker for the irreversible commitment to cell death.

### What the analysis cannot conclude

The intracellular oxidative-defense response is null across all conditions, all TPs, both omics, both panel kinds. We do *not* read this as "no oxidative stress occurs." *Prochlorococcus* MED4 lacks catalase and depends on extracellular / heterotroph-mediated peroxiredoxin activity for ROS detoxification (Black Queen Hypothesis). The intracellular transcript / protein readout of canonical oxidative-defense genes is therefore an insensitive assay for the relevant biology, and a null result here is a calibrated null, not a biological denial. Direct measurement of extracellular H₂O₂ kinetics — which is not in the KG — would be needed to test the oxidative-stress component of the Black Queen Hypothesis.

The proteomic cell-death axis declines over time in both conditions, opposite to its RNA-side engagement. We tentatively attribute this to degradation of cell-death-marker proteins themselves in dying cells, but cannot confirm without orthogonal data (e.g., viability counts at each TP, absolute protein quantification with internal standards). The cell-death axis is the only axis where RNA and protein read in opposite directions in our data.

Calendar-shared timepoints across conditions (day 31, day 89 in proteomics) are not physiologically comparable: axenic is in death phase by day 31 while coculture is still nutrient_limited. Cross-condition statistical contrast at matched calendar days would conflate "stress mitigation" with "different physiological states at the same calendar day." The analysis is therefore strictly within-condition; cross-condition reading is by trajectory shape and inspection only.

The per-gene direction calibration (psbA = -1, lrtA = -1, isiB = -1, etc.) is dataset-specific. lrtA going DOWN under N-starvation here may not generalize to other studies; another study (different *Pro* strain, different stressor) might invalidate the per-gene signs. The signed-Z methodology generalizes; the calibration table does not.

### Future work

The most natural next step is a **methodology-generalization analysis**: compute these axis_scores across all MED4 experiments in the KG (~30-50 experiments after filtering to those with ≥100 quantified genes) using the same Weissberg-2025-calibrated direction map. This would (a) ground our specific axis_scores against the broader landscape of MED4 stress responses in the KG; (b) test whether the direction calibration we derived generalizes across other stress types; (c) identify which axes are detectors-of-which-stress in a hypothesis-free way (e.g., "do oxidative-stress experiments engage our oxidative axis? do phage-challenge experiments engage proteotoxic?"). The compute is trivial (~3 minutes per experiment) and the data retrieval is bounded; the dev work is moderate (~2 hours) and primarily concerns experiment-metadata aware visualization. We propose this as a separate follow-up analysis under its own slug; it is not required to answer the current locked question but would significantly strengthen any publication-level claim that "this scoring methodology is the right way to detect stress in MED4 multi-omics data."

A second line of follow-up is **direct extracellular-ROS measurements** to test the oxidative-stress component of the Black Queen Hypothesis. The current null result for the oxidative axis is calibrated rather than biological; closing that loop requires data not in the KG.

A third line is **a cross-condition contrast analysis at the *one* physiologically-matched timepoint that may be feasible** — day 14-18 in nutrient_limited phase before axenic enters death — restricted to proteomics where both conditions sampled there. This would address F1's "no cross-condition contrast" limitation in the most defensible way and could be designed as a single follow-up question.

### Summary

The analysis answers the locked question with a multi-axis, multi-omics, statistically-defensible profile per condition. The central new finding is the **transcription-off / protein-on pattern in coculture for the N-starvation axis** — a homeostatic relief signature that is consistent with the published nitrogen-recycling mechanism and would be invisible to a single-omics analysis. The methodology, calibrated and tested in this analysis, is positioned for further use across the MED4 KG.

## References

[1] Weissberg O, Aharonovich D, Sher D. *Transcriptomic and Proteomic Analysis Reveals Nitrogen Recycling as a Core Mechanism for Prochlorococcus Prolonged Survival.* bioRxiv, 2025. DOI `10.1101/2025.11.24.690089`. Resolved via `list_publications(author="Weissberg")`.
