# Step 5 — Analyze

## Context

Step 4 verified the methodology on a single driving cell (N-stress × axenic-proteomics). This step applies the same `axis_stress_score` formula to **all 5 axes × 4 trajectory experiments × variable TPs** — full per-cell coverage of the locked question — and produces the central trajectory figures and per-cell scores. We score with **two panel kinds in parallel** as a sensitivity check:

- **Positive panel** — the 3-5 step-3-validated genes per axis, with per-gene direction calibrated from the validation findings (e.g., prxQ = -1, lrtA = -1, psbA/D/ftsH2 = -1, HLI = +1)
- **Cyanorak panel** — the broader hand-curated cyanorak functional role set per axis (~21-31 genes), with a single textbook per-axis direction default (n_stress / oxidative / proteotoxic = +1; photo (J.8 PSII) = -1). No cyanorak panel for cell_death (`D.3` is empty for MED4).

If the two panels agree qualitatively, the result is robust to gene-set choice. If they disagree, that is itself a step-5 finding.

## What I did

### Score all axes

```bash
uv run analyses/2026-04-27-1117-prochlorococcus_stress_axenic_vs_coculture/5_analyze/scripts/01_score_all_axes.py
```

[`scripts/01_score_all_axes.py`](scripts/01_score_all_axes.py) — pulls genome-wide DE for the 4 trajectory experiments (cached as `data/genome_de_<short_id>.csv`), defines positive and cyanorak gene sets per axis (cyanorak retrieved fresh via `genes_by_ontology`), and applies `axis_stress_score` per (axis × panel × experiment × TP). Output: 108 score rows in [`data/all_axes_scores.csv`](data/all_axes_scores.csv); panel definitions (gene list + per-gene direction) in [`data/panel_definitions.csv`](data/panel_definitions.csv).

Coverage at first TP per cell:

| axis | panel | proteomics_axenic | proteomics_coculture | rnaseq_axenic | rnaseq_coculture |
|---|---|---|---|---|---|
| n_stress | positive | 5 / 5 | 5 / 5 | 5 / 5 | 5 / 5 |
| n_stress | cyanorak | 27 / 28 | 27 / 28 | 28 / 28 | 28 / 28 |
| oxidative | positive | 5 / 5 | 5 / 5 | 5 / 5 | 5 / 5 |
| oxidative | cyanorak | 21 / 25 | 21 / 25 | 25 / 25 | 25 / 25 |
| proteotoxic | positive | 5 / 5 | 5 / 5 | 5 / 5 | 5 / 5 |
| proteotoxic | cyanorak | 27 / 28 | 27 / 28 | 28 / 28 | 28 / 28 |
| photo | positive | 3 / 5 | 3 / 5 | 5 / 5 | 5 / 5 |
| photo | cyanorak | 22 / 31 | 22 / 31 | 30 / 31 | 30 / 31 |
| cell_death | positive | 3 / 3 | 3 / 3 | 3 / 3 | 3 / 3 |

Proteomics misses the 2 HLI proteins from the photo positive panel and 4-9 small-protein gaps in the cyanorak panels — same pattern flagged in step 3.

### Generate trajectory figures

```bash
uv run analyses/2026-04-27-1117-prochlorococcus_stress_axenic_vs_coculture/5_analyze/scripts/02_plot_trajectories.py
```

Four figures:

- [`figures/trajectories_positive_panel.png`](figures/trajectories_positive_panel.png) — 5 axes × 2 omics, axenic vs coculture overlaid; **headline figure**
- [`figures/trajectories_positive_panel_raw.png`](figures/trajectories_positive_panel_raw.png) — same layout but `axis_mean signed log2FC` (un-normalized — answers "are the genes engaged?" without the distinctiveness-vs-genome layer)
- [`figures/trajectories_cyanorak_panel.png`](figures/trajectories_cyanorak_panel.png) — broader gene-set sensitivity check (4 axes; no cell_death)
- [`figures/panel_comparison.png`](figures/panel_comparison.png) — scatter of positive vs cyanorak per (axis × cell × TP); y = x reference line

## Results

### Per-axis findings (positive panel, signed-Z)

**N-stress** — strongest, cleanest cell-by-cell story.

| omics | axenic | coculture |
|---|---|---|
| RNA-seq | +1.45 (single point, day 14 NL) | **−0.89 → −1.07 → −1.76 → −1.25** (DOWN throughout) |
| proteomics | +1.71 → 0.93 → 1.01 (drops at death) | **+2.79 → +2.66 → +2.76 → +2.49** (high, stable, all 4 TPs) |

The **central new finding from step 5** is the cross-omics divergence in coculture: N-stress genes are strongly UP at the protein level (axis_score consistently ~2.5-2.8 across 90 days) but DOWN at the RNA level (axis_score ~−1 across all 4 TPs; raw axis_mean signed log2FC = −1.2 to −1.8). Axenic shows the opposite: both omics UP. The coculture pattern is consistent with **homeostatic relief at the RNA level** — Alteromonas-supplied N (per the published narrative) reduces the *need* for further N-scavenging transcription, while the existing N-scavenging proteome (synthesized earlier) remains stable. Step 3's heatmap was misread as "all 5 N-stress positives UP across both omics in both conditions" — the score-based view here makes the divergence quantitative and unmistakable. Note this is a step-3 reading error, not new biology — the data was always this way; I missed the negative cocult-RNA pattern when describing the heatmap. Logged as a methodology lesson in `../gaps_and_friction.md` F5.

**Photo (axis_score)**

| omics | axenic | coculture |
|---|---|---|
| RNA-seq | +0.86 (single point) | +0.31 → +0.27 → +0.42 → −0.76 |
| proteomics | +1.17 → **+3.15** → +2.27 (peak at death day 31) | +1.07 → +2.14 → +1.07 → +1.15 |

PSII disassembly accelerates dramatically in axenic-Prot at day 31 (death phase) — axis_score jumps from 1.17 (NL) to 3.15 (death). This is the largest condition-specific axis_score shift in the analysis. Coculture-Prot also rises (1.07 → 2.14) at d31, then plateaus around 1.07-1.15; axenic continues at 2.27 at d89. **The signal is direction-aware**: the PSII proteins (psbA, psbD, ftsH2) are going *down*, with direction = -1 making this a positive contribution. The coculture-RNA dip to -0.76 at day 89 is small but worth noting — it is the only cell where the photo axis is in the *wrong* (un-engaged) direction.

**Cell-death (axis_score)**

| omics | axenic | coculture |
|---|---|---|
| RNA-seq | **+1.51** (single point, the strongest RNA-side signal anywhere) | +0.83 → +0.88 → +1.24 → +0.92 |
| proteomics | +0.87 → +0.55 → +0.48 (drops over time) | +1.10 → +1.03 → +0.48 → +0.54 |

In raw response (not normalized), axenic-RNA cell-death = +2.5 — the strongest cell-death signal anywhere. The lrtA-DOWN signature (calibrated direction = -1) drives this. Curiously, the proteomic cell-death signal *drops* over time in both conditions — likely because lrtA / spoT proteins themselves degrade or fall under the detection threshold at later TPs. The cell-death axis is the only one where RNA-side and protein-side disagree on direction (RNA shows engagement, protein shows weakening signal).

**Proteotoxic (axis_score)** — weak across the board, confirming step-3 prediction. Maxes:
- RNA: axenic d14 single point +1.25; coculture d18 +0.58
- Prot: axenic d31 +0.40, d89 +0.41; coculture d31 +0.47

No cell exceeds z = +2; raw axis_mean signed log2FC stays at 0-0.6 across all cells. The dnaK1 / groL1 / groES / htpG / clpX panel is simply not strongly engaged in this experiment.

**Oxidative (axis_score)** — weakest axis. All scores between +0.19 and +0.63; never crosses any reasonable threshold for "engaged." This is exactly the sensitivity-limit warning step 3 calibrated. The intracellular oxidative-defense response is not the right assay for *Prochlorococcus* under the Black Queen Hypothesis — extracellular peroxiredoxin / heterotroph-mediated H₂O₂ handling is where the action is, and the KG does not surface that for either condition in this study.

### Cross-omics concordance (axis × condition)

| axis | axenic | coculture |
|---|---|---|
| n_stress | RNA UP + Prot UP — concordant | **RNA DOWN + Prot UP — DIVERGENT** (homeostatic relief at RNA, stable proteome) |
| photo | RNA single-point UP + Prot UP-and-rising — concordant | RNA flat-to-down + Prot UP — discordant; coculture preserves photosynthesis at RNA |
| proteotoxic | RNA single-point UP + Prot mid — concordant low | RNA mid + Prot mid — concordant low |
| oxidative | RNA-Prot both flat — concordant null | RNA-Prot both flat — concordant null |
| cell_death | RNA strong UP + Prot weak-declining — discordant | RNA modest UP + Prot weak — partially concordant |

Cross-omics divergence is concentrated in the n_stress × coculture cell. That is exactly the point in the (axis, condition) matrix where the multi-omics framing of the paper is most informative.

### Panel sensitivity check

[`figures/panel_comparison.png`](figures/panel_comparison.png) — scatter of positive-panel `axis_score` vs cyanorak-panel `axis_score`, one point per (axis, condition, omics, TP) cell where both panels exist. Findings:

- **Both panels agree qualitatively** — points cluster around the y = x line; sign of axis_score flips in only a few cells with small magnitudes (proteotoxic and oxidative, where neither signal is strong)
- **Positive panel amplifies the signal** — points lie systematically above-and-right of y = x for engaged cells (n_stress and photo): the validated positive panel produces larger absolute axis_scores than the broader cyanorak panel for the same (cell × TP). This is because the positive panel concentrates strongly-responding genes; the cyanorak panel includes many flat or weakly-responding genes that dilute the mean.
- **Positive panel detects the n_stress coculture-RNA divergence more strongly** — positive panel reads -1.7 at coculture-RNA d60; cyanorak reads -0.4 there. Both negative; positive panel quantifies the homeostatic relief more clearly.

The two-panel comparison is itself a methodology finding: the focused, validated positive panel is the more sensitive instrument; the cyanorak panel is the conservative confirmation. Both are reported in step-6 evaluation; the headline figures use the positive panel.

## Surprises

**S1. The n_stress coculture-RNA cross-omics divergence (RNA DOWN, Prot UP) was missed in step 3's heatmap interpretation.** The data was always there — N-stress genes show consistent negative log2FC in coculture-RNA across day 18, 31, 60, 89. Re-reading the step-3 heatmap (and the step 3 control_de_long.csv directly) confirms: glnB at -2.6 d89; ntcA at -2.1 d89; glnA at -2.2 d60; etc. My step-3 narrative ("All 5 N-stress positives clearly UP across both omics in both conditions") was wrong. The score-based step-5 view exposed it because the negative axis_mean is the score's primary input. **Methodology lesson:** heatmap reads are easy to anchor on the strongest cells and miss the muted-but-systematic ones. Step 3 should have included a per-(omics × condition) summary of axis-mean log2FC, not just per-gene visualization. Logged as F5.

**S2. The largest single-cell axis_score is photo × axenic-Prot × death-phase day 31 = +3.15.** PSII disassembly accelerates abruptly between day 14 (NL, +1.17) and day 31 (death, +3.15). This cleanly reproduces the published axenic-collapse narrative and is one of the few cells that exceeds z = +2 in this analysis.

**S3. The coculture-RNA photo signal goes negative at day 89 (-0.76).** All other coculture-RNA photo cells are around 0. The dip means HLI proteins are *not* upregulated and PSII proteins are *not* downregulated at day 89 in coculture — i.e., coculture cells are still photosynthesizing normally even at 90 days of N-deprivation. This is a *positive* finding for the coculture-survives narrative.

**S4. The cell-death axis is the only axis where RNA and protein disagree on direction.** RNA shows engagement (axenic single-point +1.51, coculture trajectory +0.83-1.24); proteomics drops over time (axenic 0.87 → 0.48; coculture 1.10 → 0.54). Translation: the cell-death-signaling RNAs (spoT, lrtA, isiB) are engaged but the corresponding proteins are degrading or below MS detection at later TPs. Cell-death is the *only* axis where the protein-level read is *lower* than RNA, in both conditions. Plausibly this is because cell-death markers themselves degrade in dying cells.

**S5. Both panels agree, but positive panel is more sensitive.** The cyanorak comparison was set up to test "does gene-set choice change the conclusions?" — the answer is no qualitatively, yes quantitatively. The validated 3-5-gene positive panel gives larger absolute axis_scores because it concentrates strongly-responding genes; the broader cyanorak panel is more conservative because it includes flat genes. For step-6 reporting, positive panel is primary; cyanorak is the robustness check.

## Decisions

**2026-04-27 — Positive panel is the primary report; cyanorak is the robustness check.** Both have qualitative agreement; positive panel has higher sensitivity (S5). Step-6 evaluation against the framing uses positive-panel scores; cyanorak is a sensitivity control.

**2026-04-27 — The n_stress coculture-RNA divergence is reported as a key finding, with the explicit framing that it is "transcriptional relief while proteome remains engaged."** S1 is a step-3 reading error, not new biology, but the *interpretation* — that coculture cells reduce N-scavenging transcription because Alteromonas-supplied N relieves the need to make more, while the existing N-scavenging proteome stays stable — is a step-5-derived synthesis. Step 6 will frame this against the published narrative.

**2026-04-27 — Photo and n_stress are the two informative axes; oxidative and proteotoxic are not separately reported beyond noting their weakness.** Both engaged axes (photo, n_stress) cross meaningful thresholds (raw axis_mean ≥ 1.0 or axis_score ≥ 1.5) at multiple cells. Oxidative and proteotoxic do not. Cell-death sits in between — engaged at RNA, weak at protein.

**2026-04-27 — No formal p-values in step 5.** Step 4 deferred permutation tests to step 6. Step 5 reports continuous trajectories with the conventional z = +2 line as a visual heuristic only.

## Decide-gate checklist

- **Outputs produced** —
  - `scripts/01_score_all_axes.py` (run as `uv run analyses/.../5_analyze/scripts/01_score_all_axes.py`)
  - `data/all_axes_scores.csv` (108 rows, 5 axes × 2 panels × 4 experiments × variable TPs)
  - `data/panel_definitions.csv` (gene set + direction per (axis × panel))
  - `data/genome_de_proteomics_axenic.csv`, `..._proteomics_coculture.csv`, `..._rnaseq_axenic.csv`, `..._rnaseq_coculture.csv` — full DE per experiment
  - `data/01_score_all_axes.log`
  - `scripts/02_plot_trajectories.py` (run as `uv run analyses/.../5_analyze/scripts/02_plot_trajectories.py`)
  - `figures/trajectories_positive_panel.png` + `.pdf` — **headline trajectory figure**
  - `figures/trajectories_positive_panel_raw.png` + `.pdf` — raw axis-mean (un-normalized) view
  - `figures/trajectories_cyanorak_panel.png` + `.pdf` — sensitivity check (broader gene set)
  - `figures/panel_comparison.png` + `.pdf` — positive vs cyanorak agreement scatter
- **Results presented** — per-axis tables (n_stress, photo, cell-death detailed; proteotoxic, oxidative summarized as weak); cross-omics concordance matrix; panel-comparison summary all inline above; full data and figures linked.
- **QC gate** — coverage table verifies axis-gene presence per cell (HLI proteins missing in proteomics; small cyanorak gaps in proteomics; all positive-panel cell_death and n_stress at full coverage); the driving-example axenic-prot day 14 N-stress score (+1.71) reproduces exactly what step 4 produced (cross-step consistency); the panel-comparison scatter shows qualitative agreement (both panels point the same direction in every cell with strong signal).
- **Decisions made this step** — four logged above (positive panel primary; n_stress divergence reported; only photo + n_stress informative; no p-values yet).
- **Friction logged this step** — F5 added (step-3 heatmap reading error: "all UP" claim contradicted by the actual cocult-RNA negatives; methodology lesson: include per-(omics × condition) axis-mean summary alongside the heatmap in future step-3 work).
- **Advance rationale** — every cell of the locked (5 axes × 2 omics × 2 conditions × variable TPs) matrix is scored with two panel kinds; the headline trajectory figure makes the per-cell story legible; the central new finding (n_stress coculture cross-omics divergence) is quantified and consistent across both panel kinds; step 6 can evaluate the results against step-3 framing, harvest caveats, run the deferred permutation tests if needed, and finalize the paper.
