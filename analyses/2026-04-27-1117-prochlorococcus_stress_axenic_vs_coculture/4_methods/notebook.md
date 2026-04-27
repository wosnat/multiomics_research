# Step 4 — Methods

## Context

Step 3 locked the framing: per-axis stress = coordinated DE of validated positive-control panel in the direction calibrated against the dataset (not "always up"). This step builds the **methodology** — a minimal Python module that computes a per-axis stress score, verified against hand-computed toy data, and applied to a **driving example** to confirm it works end-to-end before step 5 expands across all axes / conditions / omics.

Driving example: **N-stress axis × axenic-proteomics trajectory** — chosen because it had the cleanest validation in step 3 (5/5 positives strong, mean max|log2FC| = 3.36) and spans nutrient_limited (day 14) → death (day 31, 89), so the trajectory crosses a physiological transition.

## What I did

### The methodology

[`stress_score.py`](stress_score.py) defines `axis_stress_score(de_df, axis_genes, direction) -> dict`. For a single (experiment × timepoint) DE table, it computes:

```
signed_axis_lfc[g] = direction[g] * log2fc[g]    for g in axis_genes
axis_mean          = mean over g of signed_axis_lfc[g]

background         = log2fc for all non-axis genes (NaN dropped)
background_mean    = mean(background)
background_sd      = SD(background)              # ddof=0, population SD

axis_score         = (axis_mean - background_mean) / background_sd
```

The score is a **direction-aware signed Z** — units of "background SD of log2FC at this TP". Direction calibration (the ±1 per gene) comes from step 3 validation, **not** from a textbook default. The function fails loudly (`KeyError`) if direction is missing for any axis gene; this prevents silent default-to-+1 bugs.

Module returns axis_score, axis_mean, background_mean, background_sd, n_axis, n_background, axis_genes_with_data, axis_genes_missing_data — every diagnostic step 5 / 6 might want to reference.

### Toy verification

[`scripts/01_toy_verify.py`](scripts/01_toy_verify.py) hand-constructs six cases:

```bash
uv run analyses/2026-04-27-1117-prochlorococcus_stress_axenic_vs_coculture/4_methods/scripts/01_toy_verify.py
```

| Case | What it tests | Expected | Got | Status |
|---|---|---|---|---|
| 1 — basic, all +1 | log2fc {A:2.0, B:1.5, C:1.0} (axis +1), {D:0.1, E:0.3, F:-0.2, G:0.0} (bg) | axis_score = 8.0432 | 8.0432 | ✓ |
| 2 — all -1 with neg log2fc | symmetry: negating both should match case 1 | 8.0432 | 8.0432 | ✓ |
| 3 — bidirectional axis | mixed-direction gene set (HLI +1, psbA/D -1) | 23.3357 | 23.3357 | ✓ |
| 4 — missing direction | omitting a key must raise KeyError | KeyError | KeyError | ✓ |
| 5 — NaN in log2fc | exclude from both axis and background | n_axis=2, n_bg=1, NaN score | NaN | ✓ |
| 6 — single-gene axis | n_axis=1 must work (no min-size assumption) | finite score | 40.40 | ✓ |

The first iteration of the module docstring had `sqrt(0.13/3) ≈ 0.208` (sample variance, ddof=1) — but the implementation uses `np.std(..., ddof=0)` so the actual value is `sqrt(0.13/4) ≈ 0.180`. Toy-verification caught this; docstring corrected to match implementation. Without the toy step, this would have been a silent inconsistency between documented and actual behavior.

### Driving example

[`scripts/02_apply_n_stress_axenic_prot.py`](scripts/02_apply_n_stress_axenic_prot.py) pulls genome-wide DE for the axenic-proteomics experiment (4272 rows: 1424 genes × 3 TPs), applies `axis_stress_score` per TP using the 5 N-stress positives, saves per-TP scores, and plots a two-panel trajectory. Run:

```bash
uv run analyses/2026-04-27-1117-prochlorococcus_stress_axenic_vs_coculture/4_methods/scripts/02_apply_n_stress_axenic_prot.py
```

## Results

### N-stress × axenic-proteomics: per-TP scores

| TP | Phase | hours | axis_mean signed log2FC | bg_mean log2FC | bg_sd log2FC | **axis_score** | n_axis | n_background |
|---|---|---|---|---|---|---|---|---|
| day 14 | nutrient_limited | 336 | **1.271** | -0.038 | 0.767 | **+1.71** | 5 | 1419 |
| day 31 | death | 744 | **1.155** | -0.065 | 1.317 | **+0.93** | 5 | 1419 |
| day 89 | death | 2136 | **1.354** | +0.090 | 1.249 | **+1.01** | 5 | 1419 |

Full data → [`data/n_stress_axenic_prot_scores.csv`](data/n_stress_axenic_prot_scores.csv).
Genome-wide DE used → [`data/n_stress_axenic_prot_de.csv`](data/n_stress_axenic_prot_de.csv).
Trajectory figure → [`figures/n_stress_axenic_prot_trajectory.png`](figures/n_stress_axenic_prot_trajectory.png).

### Worked end-to-end example (day 14, nutrient_limited)

Per-axis-gene log2FC at day 14:

| Locus | Gene | log2FC |
|---|---|---|
| PMM0246 | ntcA | +1.8650 |
| PMM0920 | glnA | +2.1205 |
| PMM0263 | amt1 | +0.2705 |
| PMM1463 | glnB | +1.2747 |
| PMM0970 | urtA | +0.8244 |

axis_mean (signed, all dir = +1) = (1.8650 + 2.1205 + 0.2705 + 1.2747 + 0.8244) / 5 = **1.2710**
background_mean over the other 1419 genes = **-0.0375**
background_sd = **0.7666**
axis_score = (1.2710 − (-0.0375)) / 0.7666 = **+1.7069**

This matches the table row exactly — the formula propagates cleanly from per-gene log2FC to the per-TP score.

## Surprises

**S1. The axis score is highest at day 14 (early NL) and drops at the two death-phase TPs — even though the axis genes themselves stay strongly upregulated.** Raw signed-mean log2FC is essentially flat across the trajectory (1.16, 1.16, 1.35), but the signed-Z score drops from +1.71 to +0.93 to +1.01. The cause is visible in the data: **background SD increases from 0.77 (day 14) to 1.32 (day 31) to 1.25 (day 89)**. At death phase, the entire proteome is responding strongly in many directions; the N-stress axis no longer stands out *specifically* against the genome-wide response. The two views — raw axis response vs distinctiveness — answer different questions:

- *"Are N-stress genes still strongly engaged?"* — **yes**, raw mean signed log2FC stays around +1.3 across all TPs.
- *"Is N-stress engagement specifically elevated relative to typical genes at this TP?"* — **only at day 14**; by death phase, N-stress is one of many strongly-responding axes, not a standout.

This is a real and biologically meaningful methodology finding: **z-score-against-background measures distinctiveness, not absolute response.** When the background itself is shifted (death-phase global stress), distinctiveness goes down even if absolute response is constant. The two-panel trajectory plot makes both views visible.

**S2. The conventional z = +2 threshold is barely reached at any TP.** Day 14's score is +1.71, just below z = +2. This means a strict "axis is statistically engaged at z > 2" threshold would call no TP engaged — even though biologically we *know* axenic *Prochlorococcus* is N-stressed (it dies of N-deprivation). The threshold is too strict for a small (5-gene) panel against a heavy-tailed background. Step 5 / 6 should report scores with their raw components, not threshold them.

**S3. All 5 axis genes have data at all 3 TPs.** No missing-data complications for the driving example. The `n_axis_genes_with_data = 5` field will not always equal 5 in step 5 (HLI proteins missing from proteomics, axis-specific gaps possible) — the diagnostics field exists to make this transparent.

## Decisions

**2026-04-27 — Report axis_score AND axis_mean signed log2FC together in step 5 trajectories.** The driving example showed they answer different questions and can disagree (S1). The two-panel layout becomes the default per-(condition × omics × axis) trajectory visualization in step 5. Single-number summaries should always show both.

**2026-04-27 — Do not threshold scores in step 5 (no z > 2 cut-off).** S2 shows a strict threshold would dismiss real biology. Step 5 reports continuous trajectories; step 6 evaluation discusses thresholds in light of the gene-set size and the per-TP background SD.

**2026-04-27 — Population SD (ddof=0) is the convention for this analysis.** ddof=0 was the implementation choice; ddof=1 would change scores by a factor of ~sqrt(N/(N-1)) ≈ 1.0004 for our N=1419 background — negligible. Locking ddof=0 to keep the formula simple and matching the toy verification.

**2026-04-27 — Direction calibration is per-gene, not per-axis.** The function requires direction as `dict[locus_tag, ±1]`, not a single per-axis sign. This is required for bidirectional axes (photo: HLI = +1, psbA/psbD/ftsH2 = -1) and was tested in toy case 3.

**2026-04-27 — Driving example is sufficient for step 4 close.** Just-in-time formalization: the methodology is verified on the cleanest axis × condition × omics cell. Other axes / conditions / omics get applied in step 5 *with the same formula*. No need to pre-validate every cell here.

## Decide-gate checklist

- **Outputs produced** —
  - `stress_score.py` — module with `axis_stress_score()` and worked-example docstring
  - `scripts/01_toy_verify.py` — 6 hand-computed cases, all passing (run as `uv run analyses/.../4_methods/scripts/01_toy_verify.py`)
  - `scripts/02_apply_n_stress_axenic_prot.py` — driving-example application (run as `uv run analyses/.../4_methods/scripts/02_apply_n_stress_axenic_prot.py`)
  - `data/n_stress_axenic_prot_scores.csv` — 3 rows, per-TP scores + diagnostics
  - `data/n_stress_axenic_prot_de.csv` — 4272 rows of genome-wide DE for the driving experiment
  - `data/02_apply_n_stress_axenic_prot.log` — per-TP listing + worked example
  - `figures/n_stress_axenic_prot_trajectory.png` + `.pdf` — two-panel trajectory (raw + Z)
- **Results presented** — per-TP score table and worked end-to-end example for day 14 inline above; full data and figures linked.
- **QC gate** — toy verification 6/6 cases pass; one docstring discrepancy caught (sqrt(0.13/3) → sqrt(0.13/4)) and fixed; driving example reproduces by hand to within rounding (axis_score = (1.2710 − (-0.0375)) / 0.7666 = 1.7069 ✓); all 5 axis genes have data at all 3 TPs (no axis-side gaps in this driving cell).
- **Decisions made this step** — five logged above (two-panel default; no thresholding; ddof=0; per-gene direction; just-in-time scope close).
- **Friction logged this step** — none new. F1–F4 from prior steps remain in `../gaps_and_friction.md`.
- **Advance rationale** — methodology is a single ~70-line module with hand-computed test coverage; the driving example shows the pipeline produces coherent, biologically-interpretable output (N-stress is engaged, but distinctiveness drops at death phase as global response widens); two-panel reporting addresses the raw-vs-normalized tension. Step 5 can apply the same formula to all 4 trajectory cells × 5 axes × (where applicable) cross-omics, using the cyanorak-grounded broader gene sets where step 3 supports it.
