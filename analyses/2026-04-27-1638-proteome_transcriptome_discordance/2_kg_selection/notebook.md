# Step 2 — KG entries

## Context

Step 1 locked the question: per-gene paired-FC discordance between RNA-seq and proteomics in *Prochlorococcus* MED4 and *Alteromonas macleodii* HOT1A3 under PRO99-lowN N-starvation, axenic and coculture, in Weissberg 2025 (DOI `10.1101/2025.11.24.690089`). Eleven classification axes were retained for downstream stratification of discordant genes. F1 was logged at step-1 close (KG metadata gap — MED4 axenic RNA-seq has no timepoint label; researcher confirmed pairing to day 14).

This step (a) selects the experiments the question requires, (b) builds the paired-observation table by joining RNA-seq and proteomics on matched timepoints, (c) closes the F1 data-pipeline issue by asserting the MED4 axenic RNA-seq → day 14 mapping in code, and (d) QCs gene-detection coverage per (organism × condition) pair.

## What I did

Wrote and ran three scripts:

```bash
.venv/bin/python analyses/2026-04-27-1638-proteome_transcriptome_discordance/2_kg_selection/scripts/01_select_experiments.py
.venv/bin/python analyses/2026-04-27-1638-proteome_transcriptome_discordance/2_kg_selection/scripts/02_build_paired_observations.py
.venv/bin/python analyses/2026-04-27-1638-proteome_transcriptome_discordance/2_kg_selection/scripts/qc_gene_coverage.py
```

`01_select_experiments.py` — calls `list_experiments(publication_doi=[DOI], verbose=True, limit=None)`, filters to `omics_type ∈ {RNASEQ, PROTEOMICS}` AND `treatment_type` containing "nitrogen", attaches per-TP growth phase by zipping `timepoints` with the parallel `time_point_growth_phases` list (per-tp phase is dropped by `experiments_to_dataframe`; local helper `attach_per_tp_phase` re-attaches).

`02_build_paired_observations.py` — joins RNA-seq and proteomics on (organism, condition, timepoint_hours), inner-merge so only matched-TP pairs survive. Applies the F1 fix in code (MED4 axenic RNA-seq → 336 h / day 14 / nutrient_limited / 1849 genes). Drops the pooled "days 60+89" rows (timepoint_hours = NaN) — they are statistical-power pooling of d60 and d89 and would double-count the late-timepoint signal.

`qc_gene_coverage.py` — for each of the 4 (organism × condition) pairings, queries distinct locus_tags per experiment via `differential_expression_by_gene(experiment_ids=[exp_id], limit=None)` (which accepts `experiment_ids` without `locus_tags` and returns one row per gene × timepoint; we collapse to the locus_tag set). Reports |both|, |RNA-only|, |Protein-only|, plus a cross-condition consistency check within each organism. Cross-checked against the envelope's `matching_genes` field per call.

**Filter funnel:**

| Filter | Experiments |
|---|---|
| `publication_doi == 10.1101/2025.11.24.690089` | 10 |
| `omics_type ∈ {RNASEQ, PROTEOMICS}` ∧ `treatment_type` contains "nitrogen" | **8** |

The 2 experiments dropped at the second filter are the single-timepoint coculture-vs-axenic RNA-seq contrasts (`treatment_type=["coculture"]`, no matching proteomics — out of scope per step 1).

## Results

### Selected experiments (8 / 10 in the publication)

| Organism | Omics | Condition | Time-course? | Timepoints (h) | Distinct genes |
|---|---|---|---|---|---|
| MED4 | RNA-seq | axenic | **No** (single contrast) | (336, asserted via F1) | 1849 |
| MED4 | RNA-seq | coculture | Yes | 432, 744, 1440, 2136 (+ pooled) | 1849 |
| MED4 | Proteomics | axenic | Yes | 336, 744, 2136 | 1424 |
| MED4 | Proteomics | coculture | Yes | 432, 744, 1440, 2136 (+ pooled) | 1424 |
| HOT1A3 | RNA-seq | axenic | Yes | 432, 744 (+ pooled) | 3950 |
| HOT1A3 | RNA-seq | coculture | Yes | 432, 744, 1440, 2136 (+ pooled) | 3950 |
| HOT1A3 | Proteomics | axenic | Yes | 432, 744 | 2226 |
| HOT1A3 | Proteomics | coculture | Yes | 432, 744, 1440, 2136 (+ pooled) | 2226 |

All 8 experiments use `table_scope == "all_detected_genes"` → fair cross-experiment comparison. The single contrast on MED4 axenic RNA-seq is asserted to map to 336 h (day 14) per F1; this is documented in the script and re-stated in `../gaps_and_friction.md` F1.

Full tables: [`data/experiments_weissberg2025.csv`](data/experiments_weissberg2025.csv) (compact), [`data/experiments_weissberg2025_full.csv`](data/experiments_weissberg2025_full.csv) (full metadata), [`data/timepoints_weissberg2025.csv`](data/timepoints_weissberg2025.csv) (per-TP with growth phase).

### Paired observations (11 rows, all `nutrient_limited` phase)

| Organism | Condition | TP label | Hours | RNA-seq genes | Proteomics genes |
|---|---|---|---|---|---|
| MED4 | axenic | day 14 | 336 | 1849 | 1424 |
| MED4 | coculture | day 18 | 432 | 1849 | 1424 |
| MED4 | coculture | day 31 | 744 | 1849 | 1424 |
| MED4 | coculture | day 60 | 1440 | 1849 | 1424 |
| MED4 | coculture | day 89 | 2136 | 1849 | 1424 |
| HOT1A3 | axenic | day 18 | 432 | 3950 | 2226 |
| HOT1A3 | axenic | day 31 | 744 | 3950 | 2226 |
| HOT1A3 | coculture | day 18 | 432 | 3950 | 2226 |
| HOT1A3 | coculture | day 31 | 744 | 3950 | 2226 |
| HOT1A3 | coculture | day 60 | 1440 | 3950 | 2226 |
| HOT1A3 | coculture | day 89 | 2136 | 3950 | 2226 |

Coverage breakdown by (organism, condition): MED4 axenic 1, MED4 coculture 4, HOT1A3 axenic 2, HOT1A3 coculture 4. Total 11. All paired observations are in `nutrient_limited` growth phase — the `death`-phase MED4 axenic proteome timepoints (d31, d89) are excluded because no matching RNA-seq exists (see F2 — measurement failure ≠ biological absence; we are not interpreting protein-only timepoints as discordance).

Full table: [`data/paired_observations.csv`](data/paired_observations.csv).

### Gene-coverage QC (per organism × condition pair)

| Organism | Condition | RNA-seq | Proteomics | Both (paired pool) | RNA only | Protein only | % protein in RNA |
|---|---|---|---|---|---|---|---|
| MED4 | axenic | 1849 | 1424 | **1424** | 425 | 0 | 100.0 % |
| MED4 | coculture | 1849 | 1424 | **1424** | 425 | 0 | 100.0 % |
| HOT1A3 | axenic | 3950 | 2226 | **2221** | 1729 | 5 | 99.8 % |
| HOT1A3 | coculture | 3950 | 2226 | **2221** | 1729 | 5 | 99.8 % |

The paired-discordance analysis pool is **1424 MED4 genes** and **2221 HOT1A3 genes** (proteomics-detected ∩ RNA-seq-detected). Cross-condition consistency: within each organism the RNA-seq detected set is identical between axenic and coculture, and likewise for proteomics — so we can pool axenic and coculture observations on the gene-set side without condition-specific detection bias.

The remaining 425 MED4 RNA-only and 1729 HOT1A3 RNA-only genes are *not* available for paired-FC discordance because protein FC is missing. Whether to leave them out entirely or treat them as a "protein-undetected" category for axis-3 (architecture / size) analysis is a step-3 framing decision.

Full table: [`data/qc_gene_coverage.csv`](data/qc_gene_coverage.csv); script log: [`data/qc_gene_coverage.log`](data/qc_gene_coverage.log).

## Surprises

**S1. The 5 HOT1A3 "protein only" locus tags are all `istB` paralogs.** All 5 are annotated as IS21-like element helper ATPase IstB (locus tags `ACZ81_17130`, `ACZ81_09590`, `ACZ81_07515`, `ACZ81_04425`, `ACZ81_01340`). This is the classic paralog-resolution mismatch between MS proteomics (which assigns shared peptides to multiple genome-redundant copies) and short-read RNA-seq (which discards multi-mappers). The interpretation is methodological, not biological — these 5 transposase paralogs are removed from the paired-discordance pool by the inner-join naturally; they're flagged here so the absence is not mis-read as biology. Logged as F3.

**S2. MED4 has no protein-only genes; HOT1A3 has 5.** Closely related to S1. MED4's smaller genome and lack of recently expanded paralogous gene families produces a clean 100% protein-in-RNA inclusion. HOT1A3 has the istB family, hence the 0.2% protein-only.

**S3. Cross-condition gene-set consistency is perfect.** Both organisms' RNA-seq and proteomics detected sets are *identical* between axenic and coculture (intersection = union for each (organism × omics)). This is what we want — it means the paired-discordance pool is condition-independent, and the cross-condition consistency axis (locked axis 6) compares the same genes' fold-changes, not different gene sets.

**S4. The pooled "days 60+89" rows exist in 4 of the 8 experiments.** All 4 (multi-TP) coculture-side experiments and the HOT1A3 axenic RNA-seq carry a fifth pooled "days 60+89" row alongside the individual day-60 and day-89 rows. The pooled row presumably exists for statistical-power reasons in the original analysis. We drop the pooled rows from paired observations (decision below) — the individual TPs carry the same information at finer resolution and pooling would double-count.

## Decisions

**2026-04-27 — F1 hardcoded fix: MED4 axenic RNA-seq → 336 h / day 14.** The KG records `is_time_course=False` for this experiment with no per-TP metadata. Researcher-confirmed pairing (step 1 dialogue): the single RNA-seq contrast corresponds to day 14, the only timepoint where RNA was extractable from axenic MED4 cells. The mapping is asserted in `02_build_paired_observations.py` constants (`F1_TIMEPOINT_HOURS = 336.0`, `F1_TIMEPOINT_LABEL = "day 14"`). This effectively closes F1's data-pipeline aspect for this analysis; the underlying KG metadata gap remains and the suggested KG fix is in `../gaps_and_friction.md` F1.

**2026-04-27 — Drop pooled "days 60+89" rows from paired observations.** The pooled row is a statistical-power summary of d60 and d89 reads. We have d60 and d89 individually and prefer the finer resolution; including the pooled row would double-count the late-timepoint signal.

**2026-04-27 — Exclude the 5 HOT1A3 `istB` paralogs from the paired pool.** They drop out naturally via the inner join; documented as a known paralog-resolution mismatch (S1, F3). No analytical action — just transparency about why the paired pool size is 2221 not 2226.

## Decide-gate checklist

- **Outputs produced** —
  - `scripts/01_select_experiments.py`
  - `scripts/02_build_paired_observations.py`
  - `scripts/qc_gene_coverage.py`
  - `data/experiments_weissberg2025.csv` (8 rows, compact)
  - `data/experiments_weissberg2025_full.csv` (8 rows, full metadata, debug)
  - `data/timepoints_weissberg2025.csv` (29 TP rows, with growth phase)
  - `data/paired_observations.csv` (11 paired rows)
  - `data/qc_gene_coverage.csv` (4 organism × condition pair rows)
  - `data/01_select_experiments.log`, `data/02_build_paired_observations.log`, `data/qc_gene_coverage.log`
- **Results presented** — selected-experiments table (8 rows), paired-observations table (11 rows), gene-coverage QC table (4 (org × cond) rows) all inline above; full tables linked.
- **QC gate** — filter funnel matches expectation (10 → 8 with the 2 dropped explicitly listed); all 8 experiments use `table_scope == "all_detected_genes"`; per-TP gene counts equal each experiment's `distinct_gene_count` (no per-TP gene-count variation); F1 fix verified by inspecting the row before/after the patch in the script log; gene-coverage cross-condition consistency check confirmed identical detected sets axenic↔coculture for both organisms; the 5 HOT1A3 protein-only locus tags identified as istB paralogs (paralog-resolution mismatch, methodological).
- **Decisions made this step** — F1 hardcoded; drop pooled "days 60+89" rows; exclude 5 istB paralogs (passive — natural inner-join exclusion).
- **Friction logged this step** — F1 closed (data-pipeline aspect); F3 added (HOT1A3 istB paralog ambiguity between MS and RNA-seq); F4 added (initial `qc_gene_coverage.py` used `run_cypher` when `differential_expression_by_gene(experiment_ids=[...])` covered the case — script fixed, identical results, generalizable lesson saved as feedback memory).
- **Advance rationale** — selection complete and validated; the paired-discordance pool is well-characterized (1424 MED4 + 2221 HOT1A3 genes); 11 paired (organism × condition × timepoint) observations identified; cross-condition gene-set consistency confirmed; ready for step 3 (analysis framing) — pick a driving example and frame hypothesis, controls, expected outcome, prioritize among the 11 axes.
