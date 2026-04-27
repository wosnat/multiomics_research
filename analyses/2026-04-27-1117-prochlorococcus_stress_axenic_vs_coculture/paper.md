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

_Populated at end of step 4._

## Results

_Populated at end of step 5._

## Discussion

_Populated at end of step 6._

## References

[1] Weissberg O, Aharonovich D, Sher D. *Transcriptomic and Proteomic Analysis Reveals Nitrogen Recycling as a Core Mechanism for Prochlorococcus Prolonged Survival.* bioRxiv, 2025. DOI `10.1101/2025.11.24.690089`. Resolved via `list_publications(author="Weissberg")`.
