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

_Populated across steps 3 (framing) and 4 (implementation)._

## Results

_Populated at end of step 5._

## Discussion

_Populated at end of step 6._

## References

[1] Weissberg O, Aharonovich D, Sher D. *Transcriptomic and Proteomic Analysis Reveals Nitrogen Recycling as a Core Mechanism for Prochlorococcus Prolonged Survival.* bioRxiv, 2025. DOI `10.1101/2025.11.24.690089`. Resolved via `list_publications(author="Weissberg")`.
