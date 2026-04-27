# Step 3 — Analysis framing

## Context

Step 1 locked the question (within-condition stress trajectories per axis × omics × condition). Step 2 confirmed the data structure (5 *Pro* MED4 experiments; 3 viable multi-point trajectories + 1 single-point summary; calendar-shared TPs in different growth phases between axenic and coculture). This step's job is to (a) define **positive and negative controls** per stress axis from the KG and (b) **validate them against the actual DE data** so step 4's methodology has real ground-truth anchors.

The just-in-time-formalization principle applies: I'm not defining the full per-axis gene set here (that's step 4); I'm only defining the controls that step 4's gene-set choice has to satisfy.

## What I did

### Building the control panel

For each of the 5 stress axes, picked 3–5 canonical positive-control genes by combining:
- **Cyanorak hand-curated roles** (preferred; 73 % MED4 genome coverage, hand-curated for cyanobacteria) for N-stress (`D.1.3` Adaptation/Nitrogen + `E.4` N metabolism), oxidative (`D.1.4`), proteotoxic (`L.3` Protein folding and stabilization), and photosynthesis (`J.8` PSII)
- **GO BP** as a secondary check
- **Canonical literature markers** for genes not in either ontology (psbA, lrtA)

Negative controls picked as housekeeping candidates: `rpoB` (RNA polymerase β), `rpsL` (30S ribosomal protein S12), `atpA` (ATP synthase α).

Final panel: 26 genes (23 positives across 5 axes + 3 negatives), saved as [`data/control_panel.csv`](data/control_panel.csv).

### Annotation gaps surfaced (logged in `../gaps_and_friction.md` F4)

- `go:0006995` cellular response to nitrogen starvation: **0 MED4 genes**. Cyanorak roles fill the gap.
- `go:0008219` cell death + `go:0012501` programmed cell death + `cyanorak.role:D.3` Cell growth and death: **all 0 MED4 genes**. The cell-death axis was reframed as **late-stationary / starvation response** and anchored on the 4 genes that *do* appear in `go:0042594` (response to starvation): spoT, isiB, lytB, phoR (plus the literature canonical lrtA).
- `psbA` (PMM0223) is the canonical D1 protein but is missing from photo-stress GO BP terms; properly annotated only in `cyanorak.role:J.8`. Kept as a positive control on literature-canonical grounds.
- HLI proteins are not in any photosystem ontology term — they live only in `gene_category="Stress response and adaptation"`.
- MED4 biologically lacks catalase (Black Queen Hypothesis); not a KG gap, but worth recording so future axes don't expect it.

### Validating positive controls against the actual DE data

Wrote [`scripts/01_pull_control_de.py`](scripts/01_pull_control_de.py) and ran:

```bash
uv run analyses/2026-04-27-1117-prochlorococcus_stress_axenic_vs_coculture/3_analysis_framing/scripts/01_pull_control_de.py
```

Pulls log2FC, padj, and per-row growth phase for all 26 control genes from the 4 trajectory + 1 single-point experiments via `differential_expression_by_gene`. Outputs:

- [`data/control_de_long.csv`](data/control_de_long.csv) — 298 rows: one per (gene × experiment × TP)
- [`data/control_de_summary.csv`](data/control_de_summary.csv) — one row per gene with max |log2FC|, n significant TPs

Then [`scripts/02_plot_control_validation.py`](scripts/02_plot_control_validation.py) builds the heatmap [`figures/control_validation_heatmap.png`](figures/control_validation_heatmap.png) (26 genes × 12 cells: 4 cocult-Prot + 3 axenic-Prot + 4 cocult-RNA + 1 axenic-RNA TPs). Cells colored by log2FC (red = UP, blue = DOWN), `*` marks publication-significant rows, `†` marks axenic death-phase TPs.

### KG queries issued during this step

```
resolve_gene  (16 calls): glnA, ntcA, amt1, glnB, urtA, sodB(0), katG(0), ahpC,
              groEL(0), dnaK, clpP, psbA, psbD, rpoB, rpsL, atpA, relA(0)
genes_by_function  (7 calls): superoxide dismutase, catalase(0), chaperonin,
              high light inducible, stringent/hibernation, cyanophycin(0),
              chlorophyll degradation
search_ontology  (6 calls in go_bp + 6 in cyanorak_role): per-axis term discovery
genes_by_ontology (5 calls in go_bp + 5 in cyanorak_role): materialize per-axis sets
ontology_landscape (1 call): coverage check for our 4 trajectory experiments
differential_expression_by_gene  (4 calls): one per trajectory experiment, for the 26 controls
```

## Results

### Per-axis validation (positives)

Per-axis significance rate and magnitude across all (experiment × TP) cells:

| Axis | n positives | Genes with ≥1 sig TP | Mean max\|log2FC\| | Direction | Comment |
|---|---|---|---|---|---|
| n_stress | 5 | 5 / 5 | 3.36 | clearly UP | strongest, cleanest axis; all 5 markers respond |
| photo | 5 | 5 / 5 | 4.64 | **bidirectional** | psbA / psbD / ftsH2 DOWN (PSII disassembly); HLI UP (high-light-inducible response). PMM1404 hits +9.6 in axenic-RNA d14 — single largest log2FC in the panel. |
| proteotoxic | 5 | 4 / 5 | 1.75 | weakly UP | dnaK1 (canonical) essentially flat; htpG only responds at axenic-Prot d89†; the chaperone response is muted under N-starvation specifically |
| oxidative | 5 | 2 / 5 | 1.33 | weak / DOWN | sodN, ahpC, gor all flat; only prxQ has a real signal (and it's DOWN) |
| cell_death | 3 | 3 / 3 | 2.66 | mostly DOWN | lrtA strongly DOWN (consistent across cocult-RNA cells); spoT weakly UP/DOWN; isiB DOWN — wrong-handed for N-stress, see decision below |

### Per-gene summary table (positives — by axis)

| Axis | Locus | Gene | KG evidence | n_sig | max\|log2FC\| | Direction | Notes |
|---|---|---|---|---|---|---|---|
| n_stress | PMM0246 | ntcA | cyanorak D.1.3 + E.4 | 9 | 3.59 | UP | global N regulator; strongest in axenic Prot d31† |
| n_stress | PMM0263 | amt1 | cyanorak E.4 | 6 | 3.86 | UP | NH4 transporter |
| n_stress | PMM0920 | glnA | cyanorak E.4 | 9 | 3.36 | UP | glutamine synthetase |
| n_stress | PMM0970 | urtA | cyanorak D.1.3 + E.4 | 7 | 2.67 | UP | urea transporter |
| n_stress | PMM1463 | glnB | cyanorak D.1.3 + E.4 | 10 | 3.30 | UP | PII regulator |
| photo | PMM0223 | psbA | cyanorak J.8 | 7 | 4.06 | DOWN | D1 protein; PSII disassembly |
| photo | PMM1157 | psbD | cyanorak J.8 + go:0009769 | 6 | 4.63 | DOWN | D2 protein |
| photo | PMM0743 | ftsH2 | cyanorak L.3 + go:0010206 | 5 | 3.53 | DOWN | D1 repair protease |
| photo | PMM1404 | hli | gene_category=stress_response | 3 | **9.58** | UP | RNA-only (proteomics doesn't detect; too small) |
| photo | PMM0064 | hli | gene_category=stress_response | 1 | 1.40 | weak | RNA-only |
| proteotoxic | PMM0901 | htpG | cyanorak L.3 + go:0006457 | 1 | 1.98 | UP at d89† | only axenic-Prot late |
| proteotoxic | PMM1432 | dnaK1 | cyanorak L.3 + go:0006457/6986 | 0 | 0.80 | flat | canonical heat-shock, but unresponsive here |
| proteotoxic | PMM1436 | groL1 | cyanorak L.3 + go:0006457 | 1 | 1.71 | weak UP | |
| proteotoxic | PMM1437 | groES | cyanorak L.3 + go:0006457 | 4 | 2.10 | UP | more responsive than groL1 |
| proteotoxic | PMM1657 | clpX | cyanorak L.3 + go:0006457 | 6 | 2.15 | DOWN | protease ATPase, decreases |
| oxidative | PMM0079 | prxQ | cyanorak D.1.4 + go:0006979 | 5 | 3.02 | DOWN | the only real oxidative signal — decreasing |
| oxidative | PMM0567 | gor | cyanorak D.1.4 + go:0006979 | 0 | 0.66 | flat | glutathione reductase |
| oxidative | PMM0856 | ahpC | cyanorak D.1.4 + go:0006979 | 0 | 0.91 | flat | canonical 2-Cys peroxiredoxin |
| oxidative | PMM1061 | trxA | cyanorak D.1.4 + go:0006979 | 1 | 1.24 | weak | thioredoxin |
| oxidative | PMM1294 | sodN | cyanorak D.1.4 + go:0006979 | 0 | 0.82 | flat | Ni-SOD — canonical, but flat under N-stress |
| cell_death | PMM0191 | spoT | go:0015968 + go:0042594 | 2 | 2.09 | mixed | stringent response |
| cell_death | PMM0400 | lrtA | literature canonical | 10 | 3.30 | DOWN | ribosome hibernation factor; consistently DOWN under N-stress |
| cell_death | PMM1171 | isiB | go:0042594 | 4 | 2.58 | DOWN | flavodoxin — actually an iron-stress marker; wrong-handed for N-stress |

### Negative-control behavior

| Locus | Gene | n_sig | max\|log2FC\| | Comment |
|---|---|---|---|---|
| PMM1485 | rpoB | 1 | 1.46 | cleanest; near-flat across cells |
| PMM1511 | rpsL | 3 | 2.15 | mostly flat with one outlier TP |
| PMM1451 | atpA | 3 | **4.73** | dramatic DOWN at axenic-RNA d18 — **bad as a flat-housekeeping reference**, but informative as an "energy collapse" marker |

### Cross-omics concordance (qualitative)

For most positive-control genes, RNA and protein signals are direction-concordant where both omics quantify the gene:

- **ntcA, glnA, glnB, urtA, amt1**: UP in both RNA and Prot, with similar magnitudes
- **groES, htpG, clpX**: direction-concordant; magnitudes differ
- **psbA, psbD**: direction-concordant DOWN in both, but the **protein signal is much stronger in axenic late TPs than the RNA signal** — protein levels at death phase reflect *accumulated* PSII loss; RNA is shutting down in parallel
- **HLI proteins**: RNA-only (proteomics doesn't detect; too small for bottom-up MS) — so cross-omics concordance is structurally undefined for HLI

This means cross-omics concordance is a viable check for most axes, but the photo axis specifically requires interpreting "RNA falling slowly while protein has crashed" as expected behavior at death-phase rather than as omics disagreement.

## Surprises

**S1. Oxidative axis is structurally weakly engaged in this experiment.** The canonical oxidative-stress responders sodN, ahpC, gor all fail to respond significantly across all 12 cells, and only prxQ has a real signal (and it goes DOWN, not up). This is biologically sensible — *Prochlorococcus* MED4 lacks catalase and depends on heterotrophs / extracellular peroxiredoxin activity for ROS detox (Black Queen Hypothesis), so the *intracellular* oxidative-defense transcriptional response may simply not be the right place to look for "oxidative stress" in this organism. **Implication: any "no oxidative stress" finding from this analysis must be reported with the explicit caveat that the assay is insensitive to the most likely mechanism (extracellular ROS handling).**

**S2. Direction of the stress signature is axis-dependent, not always UP.** N-stress markers go up. PSII proteins go DOWN (disassembly). Ribosome hibernation factor lrtA goes DOWN (not the bacterial-textbook "ribosome dormancy under starvation" pattern). atpA crashes dramatically (energy collapse). The methodology cannot assume "stress = upregulation"; step 4 scoring must respect axis-specific direction.

**S3. HLI PMM1404 reaches log2FC = +9.6 at axenic-RNA d14.** The largest single log2FC magnitude in the panel, ~12 cells. This is one of the strongest stress signals in the entire experiment. The HLI family broadly behaves as a generalized stress responder in cyanobacteria — its photo-stress label here is conventional but the gene is responsive to many stresses simultaneously.

**S4. atpA is a bad flat-housekeeping reference but a great energy-collapse readout.** For the analysis it means our negative-control bar is not as flat as we'd like; rpoB is the only really-flat one. But the dramatic atpA crash (max|log2FC|=4.7 in axenic-RNA d18) is itself biological information — it signals energy/ATP-synthesis breakdown, which fits the death-phase narrative for axenic *Prochlorococcus*.

**S5. lrtA strongly DOWN under N-stress.** The classical assumption is "starvation → ribosome hibernation → lrtA UP." Here it goes DOWN consistently (max|log2FC|=3.30, 10 sig TPs out of 12). Possibilities: (i) lrtA in *Prochlorococcus* MED4 is *light*-regulated rather than starvation-regulated despite the cyanobacterial-stress literature, and the cells under continuous light keep it down; (ii) overall transcription decline includes lrtA. Either way: **lrtA is a real cell-death-axis marker for this experiment, but in the inverted direction.**

## Decisions

**2026-04-27 — Stress is defined as coordinated DE in the axis panel, agnostic to direction.** Per-axis "stress engagement" is detected as positive controls showing significant DE in coordinated direction (whatever the direction is per axis, as observed in the validation), not as "all genes go up." This is forced by the bidirectional photo and cell-death signatures (S2).

**2026-04-27 — `isiB` reframed but kept in panel.** isiB (PMM1171) goes DOWN under N-stress in this experiment, not up. Likely because flavodoxin is an Fe-stress marker, not N-stress. **Keeping isiB in the panel transparently** as a "wrong-handed marker" — informative because (a) it confirms the cell-death axis is direction-bidirectional; (b) it provides a reality check on the axis interpretation; (c) replacing it would require finding another late-stationary marker, and the cell-death axis is structurally thin in MED4 annotations regardless. Adding a note to step 4: when defining the cell-death gene set, isiB should not be included as an "increasing-with-stress" marker.

**2026-04-27 — Negative-control `atpA` kept but downgraded.** atpA's dramatic axenic-RNA d18 signal (-4.7) is real and informative. We keep it in the panel for transparency but treat it as **a known-non-flat reference** — a "second-tier" negative. rpoB is the primary flat reference. This is more honest than swapping to another candidate that might also turn out responsive in 90-day starvation.

**2026-04-27 — Oxidative axis has a known sensitivity limit.** The intracellular oxidative-defense transcriptional/translational response is weakly engaged in this experiment regardless of condition. Step 6 evaluation must read "no significant oxidative axis signal" as **"insensitive to extracellular peroxide / heterotroph-mediated ROS handling"**, not as "no oxidative stress." This is also the central biological framing of the Black Queen Hypothesis, so the limitation is methodologically appropriate.

**2026-04-27 — HLI proteins are RNA-only for the photo axis.** Proteomics does not detect HLI (too small for bottom-up MS). Step 4 must build the photo-stress axis with separate RNA and protein gene sets — RNA includes HLI; protein excludes them. Cross-omics concordance for the photo axis is computed only over the genes detectable in both omics (psbA, psbD, ftsH2).

## Hypothesis statement

**If *Prochlorococcus* MED4 is stressed in a given axis (N-starvation / oxidative / proteotoxic / photosynthetic / late-stationary) within a given condition (axenic, coculture), then across the within-condition trajectory from exponential phase through nutrient_limited and (where present) death-phase TPs, the canonical positive-control genes for that axis will show coordinated, statistically significant differential expression in the direction calibrated by this validation step (UP for n_stress; DOWN for psbA/psbD/ftsH2 in photo and UP for HLI; DOWN for lrtA / mixed for spoT in cell_death; muted to weakly UP for proteotoxic; mostly flat with some DOWN for oxidative).**

The within-condition baseline is exponential phase (log2FC ≈ 0 by experimental construction; the per-experiment DE is already vs the same condition's exponential phase). The negative controls (rpoB primary; rpsL, atpA secondary) anchor the "what flat looks like" baseline.

## What success looks like operationally (step-5 deliverable in framing terms)

For each (axis × condition × omics) cell, the step-5 output is:

1. A **stress score trajectory** over TPs (single point for axenic-RNA), aggregating the signed contributions of axis positive-control genes calibrated by direction (this validation step's findings)
2. A **direction-aware verdict** per (axis × condition × omics): "engaged at TP X (direction Y, magnitude Z); not engaged at TP W (magnitude near baseline)"
3. A **cross-omics concordance flag** per (axis × condition): RNA vs protein direction agreement, with the photo-axis caveat (HLI RNA-only) and the photo-axis time-lag caveat (death-phase protein loss exceeds RNA loss)
4. An **explicit caveat** for the oxidative axis: insensitive-to-extracellular-ROS-handling

## Decide-gate checklist

- **Outputs produced** —
  - `data/control_panel.csv` (26 genes: 23 positives across 5 axes + 3 negatives)
  - `scripts/01_pull_control_de.py` (run as `uv run analyses/2026-04-27-1117-prochlorococcus_stress_axenic_vs_coculture/3_analysis_framing/scripts/01_pull_control_de.py`)
  - `data/control_de_long.csv` (298 rows; 26 genes × 12 (experiment, TP) cells, less the single-point axenic-RNA collapsed and pooled days-60+89 dropped per step 2 decision)
  - `data/control_de_summary.csv` (per-gene summary)
  - `data/01_pull_control_de.log`
  - `scripts/02_plot_control_validation.py` (run as `uv run analyses/2026-04-27-1117-prochlorococcus_stress_axenic_vs_coculture/3_analysis_framing/scripts/02_plot_control_validation.py`)
  - `figures/control_validation_heatmap.png` + `.pdf`
- **Results presented** — per-axis validation table, per-gene summary table, negative-control behavior table, cross-omics concordance discussion, hypothesis statement, operational-success criterion all inline above; full data and figures linked.
- **QC gate** — controls resolved either via cyanorak hand-curated roles (preferred) or canonical literature; KG-confirmed locus tags; per-row growth phase from `differential_expression_by_gene` agrees with the F3 manual zip; positive controls validated by direct DE inspection (not assumed); 5/5 N-stress + 5/5 photo + 4/5 proteotoxic + 3/3 cell-death + 2/5 oxidative positives have at least one significant TP.
- **Decisions made this step** — five logged above (direction-agnostic stress definition; isiB kept as wrong-handed marker; atpA kept as second-tier negative; oxidative-axis sensitivity limit documented; HLI RNA-only for photo axis).
- **Friction logged this step** — F4 (annotation gaps for stress axes in MED4) added to `../gaps_and_friction.md` covering the GO BP nitrogen-starvation gap, the cell-death gap (GO BP and cyanorak both empty), the photo-axis psbA-not-in-GO gap, the HLI-not-in-PSII-ontology cross-reference gap, and the catalase-absent (biological, not KG) note.
- **Advance rationale** — the locked question's framing now has (a) a 26-gene control panel grounded in cyanorak/GO/literature evidence, (b) per-axis direction calibration validated against actual DE, (c) explicit per-axis sensitivity limits (oxidative weak; HLI RNA-only; cell-death has bidirectional markers), and (d) a hypothesis + operational-success criterion that step 4 can build a methodology against. Step 4 can now define per-axis gene sets that include these positives, choose a direction-aware scoring function, and verify with toy-data tests before applying to all genes.
