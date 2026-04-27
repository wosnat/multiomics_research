# Step 2 — KG entries

## Context

Step 1 locked the question: within-condition trajectories of five stress axes for *Prochlorococcus* MED4 from Weissberg 2025, scored independently in transcriptome and proteome, with exponential phase as the within-condition baseline. F1 (logged at step 1 close) noted that cross-condition contrast at matched TPs is not supported by the data; the specific cause was deferred to this step.

This step (a) selects the experiments that the question requires, (b) characterizes the TP / omics / replicate structure, and (c) closes out F1's cause.

## What I did

Wrote and ran `scripts/01_select_experiments.py`:

```bash
.venv/bin/python analyses/2026-04-27-1117-prochlorococcus_stress_axenic_vs_coculture/2_kg_selection/scripts/01_select_experiments.py
```

Calls `list_experiments(publication_doi=["10.1101/2025.11.24.690089"], verbose=True, limit=None)`, post-hoc filters on `organism_name == "Prochlorococcus MED4"`, attaches per-TP growth phase by zipping `timepoints` with the parallel `time_point_growth_phases` list, and emits three CSVs plus a log.

**Filter funnel:**

| Filter | Experiments |
|---|---|
| `publication_doi == 10.1101/2025.11.24.690089` | 10 |
| `organism_name == "Prochlorococcus MED4"` (post-hoc) | 5 |

The Python API's `organism=` parameter does partial-match on profiled organism *and* coculture partner, so passing it directly returned 6 experiments (one HOT1A3 experiment with MED4 as partner leaked in). Logged as F2 in `../gaps_and_friction.md`; the script uses post-hoc filtering instead.

## Results

### Selected experiments (5 / 10 in the publication)

| Experiment ID (suffix) | Omics | Condition | Time-course? | TPs | Distinct genes | `table_scope` |
|---|---|---|---|---|---|---|
| `..._coculture_alteromonas_hot1a3_med4_rnaseq` | RNA | coculture **vs** axenic (direct contrast at exponential) | no | single (exp day 11) | 1849 | `all_detected_genes` |
| `..._med4_rnaseq_coculture` | RNA | coculture (vs own exponential) | yes | day 18, 31, 60, 89, [60+89] | 1849 | `all_detected_genes` |
| `..._med4_rnaseq_axenic` | RNA | axenic (vs own exponential) | **no** | single (nutrient_limited, no per-TP) | 1849 | `all_detected_genes` |
| `..._med4_proteomics_coculture` | Prot | coculture (vs own exponential) | yes | day 18, 31, 60, 89, [60+89] | 1424 | `all_detected_genes` |
| `..._med4_proteomics_axenic` | Prot | axenic (vs own exponential) | yes | day 14, 31, 89 | 1424 | `all_detected_genes` |

Full table → [`data/experiments_pro_med4.csv`](data/experiments_pro_med4.csv); per-TP table with growth phases → [`data/timepoints_pro_med4.csv`](data/timepoints_pro_med4.csv); raw flat dump → [`data/experiments_pro_med4_full.csv`](data/experiments_pro_med4_full.csv); script log → [`data/01_select_experiments.log`](data/01_select_experiments.log).

### How each experiment maps to the locked question

| Question cell | Source experiment | Notes |
|---|---|---|
| axenic × RNA × trajectory | `..._med4_rnaseq_axenic` | **degenerate**: single nutrient_limited point, no per-TP breakdown |
| axenic × Prot × trajectory | `..._med4_proteomics_axenic` | 3 TPs spanning nutrient_limited (day 14) and death (day 31, 89) |
| coculture × RNA × trajectory | `..._med4_rnaseq_coculture` | 4 TPs, all nutrient_limited |
| coculture × Prot × trajectory | `..._med4_proteomics_coculture` | 4 TPs, all nutrient_limited |
| (cross-condition contrast at exponential) | `..._coculture_alteromonas_hot1a3_med4_rnaseq` | available in KG but out of scope per step 1; would test "is coculture different from axenic at exponential" — uninformative for stress |

The 4 trajectory experiments are structured as **"PRO99-lowN nutrient starvation" treatment vs "PRO99-lowN exponential growth" control**, both within the same condition. So each TP entry already encodes log2FC vs that condition's own exponential baseline. This is exactly the within-condition trajectory framing locked at step 1 — the data structure matches our framing without further reshaping.

The "`days 60+89`" rows present in the coculture experiments are pooled TPs (no `timepoint_hours`); the original analysis combined day 60 and day 89 reads into a single contrast, presumably for statistical power reasons. We have day 60 and day 89 separately, so we will **exclude the pooled `days 60+89` row** from trajectory scoring (no information gain over the separate days; treating it as an extra TP would mis-weight the late timepoints). Decision logged below.

### Per-TP growth phase, with cross-condition overlap

| Omics | Axenic (TP, phase) | Coculture (TP, phase) | Calendar-shared TPs | Phase-shared (TP + phase) |
|---|---|---|---|---|
| RNA-seq | — (no per-TP data) | day 18 (NL), day 31 (NL), day 60 (NL), day 89 (NL) | none | none |
| Proteomics | day 14 (NL), day 31 (**death**), day 89 (**death**) | day 18 (NL), day 31 (NL), day 60 (NL), day 89 (NL) | day 31, day 89 | **none** |

(NL = nutrient_limited.)

This is the F1 closeout. The two structural problems for cross-condition contrast are:

1. **Axenic RNA-seq is not a time-course** — 1 collapsed point, so per-TP RNA contrast across conditions is structurally impossible.
2. **Phase mismatch on calendar-shared proteomics TPs** — at day 31 and day 89, axenic is in *death* phase while coculture is still *nutrient_limited*. A DE contrast at those days would not isolate "Alteromonas mitigates stress" from "axenic and coculture are in fundamentally different physiological states at the same calendar day." That phase mismatch is itself biologically informative — it shows that axenic *Prochlorococcus* enters death substantially earlier than coculture under N-deprivation — but it precludes a clean matched-TP statistical contrast.

Full text in `../gaps_and_friction.md` F1.

## Surprises

**S1. The axenic RNA-seq experiment is not a time-course.** Going in, I assumed the 4 trajectory experiments would be symmetric (axenic × {RNA, Prot} and coculture × {RNA, Prot}, all multi-TP). They are not: axenic-RNA collapses all TPs into a single nutrient_limited contrast vs exponential. This degrades the locked plan to **3 viable per-TP trajectories** (axenic-Prot 3-point, coculture-Prot 4-point, coculture-RNA 4-point) plus 1 single-point summary (axenic-RNA). Affects step 5 figure design and step 6 evaluation: any cross-omics concordance check for axenic is forced to compare a single-point RNA summary against a 3-point protein trajectory.

**S2. Calendar-shared proteomics TPs are in different growth phases.** Day 31 axenic is already in *death* phase; day 31 coculture is in *nutrient_limited*. Same for day 89. This is biologically meaningful — the time-to-death difference between conditions is part of the published story — but it means "shared TP" is a misleading concept; the right alignment for cross-condition reading is by *physiological state*, not calendar day. Step 5 plots will be against `timepoint_hours` with phase annotations.

**S3. The pooled `days 60+89` TP exists in both coculture experiments.** Looks like the original analysis added a combined-late-timepoint contrast for statistical power. We have the individual days, so we drop the pooled row from trajectory scoring (decision below). Worth mentioning in the methods write-up so readers don't ask why our trajectories are 4-point and not 5-point.

## Decisions

**2026-04-27 — Drop the pooled `days 60+89` rows from trajectory scoring.** Both coculture experiments contain a fifth pooled "days 60+89" row alongside the individual day 60 and day 89 rows. Including it would double-count the late-timepoint signal (or worse, downweight it in some normalizations). The individual TPs carry the same information at finer resolution. Used only when explicitly noting it; not part of the per-TP trajectory.

**2026-04-27 — Treat axenic-RNA as a single-point summary, not a trajectory.** The KG records `is_time_course=False` for `..._med4_rnaseq_axenic`. We will report the axenic RNA stress score as one number (single nutrient_limited contrast vs exponential), not a trajectory. The coculture-RNA, coculture-Prot, and axenic-Prot trajectories remain multi-point. Cross-omics concordance check for axenic must compare 1 RNA point against the 3 protein TPs explicitly (we are not interpolating an axenic-RNA trajectory).

**2026-04-27 — TP alignment for cross-condition visual reading uses `timepoint_hours` and phase annotation, not calendar-day labels.** S2 makes calendar-day alignment misleading. Step 5 trajectory plots will use `timepoint_hours` on the x-axis and annotate growth phase per point.

## Decide-gate checklist

- **Outputs produced** —
  - `scripts/01_select_experiments.py` (run as: `.venv/bin/python analyses/2026-04-27-1117-prochlorococcus_stress_axenic_vs_coculture/2_kg_selection/scripts/01_select_experiments.py`)
  - `data/experiments_pro_med4.csv` (5 rows, compact)
  - `data/experiments_pro_med4_full.csv` (5 rows, all metadata, debug)
  - `data/timepoints_pro_med4.csv` (15 rows: 5 experiments × variable TPs, with per-TP growth phase)
  - `data/01_select_experiments.log` (filter funnel + per-TP listing + cross-condition overlap)
- **Results presented** — selected experiments table (5 rows), question-cell mapping table, per-TP growth-phase + cross-condition overlap table all inline above; full tables linked.
- **QC gate** — filter funnel matches expectation post-fix (`publication=10 → organism_name=5`); each TP's `tp_gene_count` equals the experiment's `distinct_gene_count` (1424 for proteomics, 1849 for RNA-seq) — no per-TP gene-count surprises; `table_scope == "all_detected_genes"` for all 5 experiments → fair cross-experiment comparison; per-TP growth phase attached correctly (axenic-Prot day 31 = death, day 89 = death; all coculture TPs = nutrient_limited).
- **Decisions made this step** — drop pooled `days 60+89` rows; treat axenic-RNA as single-point summary; TP alignment uses `timepoint_hours` + phase annotation in plots.
- **Friction logged this step** — F1 fully characterized; F2 added (API `organism=` filter ambiguity); F3 added (`experiments_to_dataframe` does not expand `time_point_growth_phases` per TP row — workaround via local helper `attach_per_tp_phase()`). All in `../gaps_and_friction.md`.
- **Advance rationale** — selection complete and validated; structure and limits of the data are now explicit; step 3 framing can build directly on the per-experiment DE-vs-exponential structure with a clear set of trajectory cells to score.
