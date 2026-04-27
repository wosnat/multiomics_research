# Step 2 — KG entries

## Context

Step 1 locked the question: per-gene paired-FC discordance between RNA-seq and proteomics in *Prochlorococcus* MED4 and *Alteromonas macleodii* HOT1A3 under PRO99-lowN N-starvation, axenic and coculture, in Weissberg 2025 (DOI `10.1101/2025.11.24.690089`). Eleven classification axes were retained for downstream stratification of discordant genes. F1 was logged at step-1 close (KG metadata gap — MED4 axenic RNA-seq has no timepoint label; researcher confirmed pairing to day 14).

This step (a) selects the experiments the question requires, (b) builds the paired-observation table by joining RNA-seq and proteomics on matched timepoints, (c) closes the F1 data-pipeline issue by asserting the MED4 axenic RNA-seq → day 14 mapping in code, and (d) QCs gene-detection coverage per (organism × condition) pair.

## What I did

Wrote and ran four scripts:

```bash
.venv/bin/python analyses/2026-04-27-1638-proteome_transcriptome_discordance/2_kg_selection/scripts/01_select_experiments.py
.venv/bin/python analyses/2026-04-27-1638-proteome_transcriptome_discordance/2_kg_selection/scripts/02_build_paired_observations.py
.venv/bin/python analyses/2026-04-27-1638-proteome_transcriptome_discordance/2_kg_selection/scripts/qc_gene_coverage.py
.venv/bin/python analyses/2026-04-27-1638-proteome_transcriptome_discordance/2_kg_selection/scripts/03_paired_fc_overview.py
```

`01_select_experiments.py` — calls `list_experiments(publication_doi=[DOI], verbose=True, limit=None)`, filters to `omics_type ∈ {RNASEQ, PROTEOMICS}` AND `treatment_type` containing "nitrogen", attaches per-TP growth phase by zipping `timepoints` with the parallel `time_point_growth_phases` list (per-tp phase is dropped by `experiments_to_dataframe`; local helper `attach_per_tp_phase` re-attaches).

`02_build_paired_observations.py` — joins RNA-seq and proteomics on (organism, condition, timepoint_hours), inner-merge so only matched-TP pairs survive. Applies the F1 fix in code (MED4 axenic RNA-seq → 336 h / day 14 / nutrient_limited / 1849 genes). Drops the pooled "days 60+89" rows (timepoint_hours = NaN) — they are statistical-power pooling of d60 and d89 and would double-count the late-timepoint signal.

`qc_gene_coverage.py` — for each of the 4 (organism × condition) pairings, queries distinct locus_tags per experiment via `differential_expression_by_gene(experiment_ids=[exp_id], limit=None)` (which accepts `experiment_ids` without `locus_tags` and returns one row per gene × timepoint; we collapse to the locus_tag set). Reports |both|, |RNA-only|, |Protein-only|, plus a cross-condition consistency check within each organism. Cross-checked against the envelope's `matching_genes` field per call.

`03_paired_fc_overview.py` — for each of the 11 paired observations, fetches per-gene log2FC + significance for both omics, inner-joins on locus_tag (the paired-detection pool), and produces (a) a per-observation correlation summary (Pearson r, Spearman ρ across all paired genes), (b) a per-gene paired category by significance × direction (`concordant_up`, `concordant_down`, `rna_only_up`, `rna_only_down`, `prot_only_up`, `prot_only_down`, `opposite_rna_up`, `opposite_rna_down`, `both_ns`), collapsed for display into compact 4-bucket categories (`concordant_sig`, `asymmetric_sig`, `opposite_sig`, `both_ns`), (c) an 11-panel scatter grid of (log2FC_RNA vs log2FC_protein), and (d) a stacked bar of compact-category percentages per observation. F1 fix re-applied (override `timepoint_hours = 336.0` on rows from the MED4 axenic RNAseq experiment). gene_name is coalesced from both omics (`gene_name_rna`, fallback to `gene_name_prot`); locus_tag remains the primary identifier (Rule 2). Note: 30 % of MED4 paired-pool genes and 38 % of HOT1A3 paired-pool genes have null `gene_name` in the KG annotations — never filter or join by `gene_name`.

`qc_distributions_and_markers.py` — four step-2 QC checks: (1) per (experiment, timepoint) log2FC distribution stats (median, IQR, |max|, σ); (2) NaN log2FC counts in the paired pool, per observation; (3) per-TP significance crosscheck — row-level `significant_up` / `significant_down` counts vs the experiment-node `tp_significant_up` / `tp_significant_down` fields, with fallback to experiment-level `genes_by_status_*` for single-contrast experiments where per-TP fields are NaN by structural design; (4) canonical N-stress marker panel (ntcA, glnA, glnB, urtA, amt1) plus housekeeping references (rpoB, atpA), resolved via `resolve_gene` to locus_tags then filtered by locus_tag (not gene_name; many genes have null gene_name).

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

### Paired RNA / protein agreement overview

A purely descriptive look at how mRNA fold-changes track protein fold-changes per gene, before any discordance metric is defined. Per paired observation, every gene in the paired-detection pool gets a category by (RNA significance × protein significance × direction). Genes where neither omic significantly responds (`both_ns`) form the bulk of every observation. Genes where one omic significantly responds and the other does not (`asymmetric_sig`) form the second-largest bucket and **dominate over genes where both omics significantly respond in the same direction (`concordant_sig`) by an order of magnitude across the full data**.

#### Pearson / Spearman correlations of per-gene log2FC

| Organism | Condition | TP | n (paired) | Pearson r | Spearman ρ |
|---|---|---|---|---|---|
| HOT1A3 | axenic | day 18 | 2221 | −0.01 | −0.05 |
| HOT1A3 | axenic | day 31 | 2221 | −0.03 | −0.07 |
| HOT1A3 | coculture | day 18 | 2221 | −0.09 | −0.11 |
| HOT1A3 | coculture | day 31 | 2221 | +0.03 | −0.01 |
| HOT1A3 | coculture | day 60 | 2221 | −0.03 | −0.03 |
| HOT1A3 | coculture | day 89 | 2221 | +0.00 | −0.03 |
| MED4 | axenic | day 14 | 1424 | **+0.23** | +0.23 |
| MED4 | coculture | day 18 | 1424 | −0.06 | −0.08 |
| MED4 | coculture | day 31 | 1424 | +0.12 | +0.10 |
| MED4 | coculture | day 60 | 1424 | **+0.18** | +0.18 |
| MED4 | coculture | day 89 | 1424 | **+0.22** | +0.21 |

#### Compact-category breakdown (% of paired pool)

| Observation | concordant_sig | asymmetric_sig | opposite_sig | both_ns |
|---|---:|---:|---:|---:|
| HOT1A3 axenic day 18 | 0.0 % | 16.3 % | 0.0 % | 83.7 % |
| HOT1A3 axenic day 31 | 3.2 % | **46.1 %** | 3.1 % | 47.6 % |
| HOT1A3 coculture day 18 | 0.8 % | 21.7 % | 0.7 % | 76.8 % |
| HOT1A3 coculture day 31 | 3.1 % | 32.5 % | 1.9 % | 62.5 % |
| HOT1A3 coculture day 60 | 2.8 % | 35.5 % | 3.2 % | 58.6 % |
| HOT1A3 coculture day 89 | 3.2 % | 35.9 % | 3.7 % | 57.1 % |
| MED4 axenic day 14 | 4.5 % | 45.7 % | 1.0 % | 48.8 % |
| MED4 coculture day 18 | 0.1 % | 5.4 % | 0.0 % | **94.5 %** |
| MED4 coculture day 31 | 2.5 % | 28.5 % | 1.3 % | 67.7 % |
| MED4 coculture day 60 | 2.8 % | 24.4 % | 1.3 % | 71.6 % |
| MED4 coculture day 89 | 2.5 % | 23.7 % | 0.6 % | 73.2 % |

**Cross-observation pooled (all 20 446 paired-gene-observations):** concordant_sig = 2.3 %, asymmetric_sig = 29.3 %, opposite_sig = 1.6 %, both_ns = 66.8 %.

Figures:
- [`figures/paired_fc_scatter_grid.png`](figures/paired_fc_scatter_grid.png) — 11-panel scatter of log2FC mRNA vs log2FC protein per observation, color by compact category, with y = x reference and Pearson r in panel title.
- [`figures/paired_fc_category_bars.png`](figures/paired_fc_category_bars.png) — stacked bars of compact-category percentages per observation.

Long-format per-gene table: [`data/paired_fc_long.csv`](data/paired_fc_long.csv); per-observation summary: [`data/paired_fc_summary.csv`](data/paired_fc_summary.csv).

### QC #1 — log2FC distribution per (experiment, timepoint)

29 distribution rows (one per experiment × per-TP, including the NaN-hours pooled rows). All distributions are roughly symmetric around zero (medians in [−0.16, +0.07]). RNAseq experiments have wider IQR than proteomics experiments (RNAseq IQR ≈ 0.7–1.5; proteomics IQR ≈ 0.3–0.9), consistent with the known dynamic-range compression of bottom-up MS proteomics relative to RNAseq counts. Maximum |log2FC| reaches +9.6 / −9.6 in some RNAseq TPs (HLI family, known biology). No scaling pathologies; no experiment is on a different log-base.

Full table: [`data/qc_fc_distributions.csv`](data/qc_fc_distributions.csv); figure: [`figures/qc_fc_distribution_box.png`](figures/qc_fc_distribution_box.png).

### QC #2 — NaN log2FC counts in the paired pool

Zero NaN log2FC across every (organism × condition × timepoint) paired observation, on both omics sides. The 1424 / 2221 paired-pool sizes therefore equal the n_with_fc used for correlation; the NaN handling caveat does not apply to this dataset.

Full table: [`data/qc_paired_nan_summary.csv`](data/qc_paired_nan_summary.csv).

### QC #3 — per-TP significance crosscheck

Row-level `expression_status` counts (`significant_up`, `significant_down`) compared against the experiment-node `tp_significant_up` / `tp_significant_down` fields. For the MED4 axenic RNAseq experiment (`is_time_course = False`, per-TP node fields NaN by structural design), the comparison falls back to the experiment-level `genes_by_status_*` fields. **29 / 29 ✓** — all per-TP and the one fallback case agree exactly. Catches API-vs-node-metadata drift.

Full table: [`data/qc_significance_crosscheck.csv`](data/qc_significance_crosscheck.csv).

### QC #4 — canonical N-stress markers (MED4 only)

Resolved via `resolve_gene(name, organism="MED4")`, filtered downstream by locus_tag (not gene_name) per Rule 2.

| Marker | Locus tag | Product | Role |
|---|---|---|---|
| ntcA | PMM0246 | global nitrogen regulatory protein | concordance positive control |
| glnA | PMM0920 | glutamine synthetase, type I | concordance positive control |
| glnB | PMM1463 | nitrogen regulatory protein P-II | concordance positive control |
| urtA | PMM0970 | ABC-type urea transporter, substrate binding component | concordance positive control |
| amt1 | PMM0263 | ammonium transporter | concordance positive control |
| rpoB | PMM1485 | DNA-directed RNA polymerase, beta subunit | housekeeping reference |
| atpA | PMM1451 | ATP synthase alpha chain | housekeeping reference (imperfect) |

#### Marker log2FC across the 5 MED4 paired observations

| Marker | ax d14 RNA / Prot | cc d18 RNA / Prot | cc d31 RNA / Prot | cc d60 RNA / Prot | cc d89 RNA / Prot |
|---|---|---|---|---|---|
| **ntcA** | +3.4 / +1.9 (concordant_up) | −1.5 / +1.8 (prot_only_up) | −1.6 / +2.2 (prot_only_up) | −1.7 / +1.3 (prot_only_up) | **−2.1 / +1.2** (opposite_rna_down) |
| **glnA** | +3.4 / +2.1 (concordant_up) | −0.5 / +2.6 (prot_only_up) | −1.7 / +2.9 (prot_only_up) | **−2.2 / +2.6** (opposite_rna_down) | −1.3 / +2.2 (prot_only_up) |
| **glnB** | −0.5 / +1.3 (prot_only_up) | −1.0 / +1.3 (prot_only_up) | **−2.0 / +2.1** (opposite_rna_down) | **−1.7 / +3.0** (opposite_rna_down) | **−2.6 / +3.3** (opposite_rna_down) |
| **urtA** | +2.7 / +0.8 (rna_only_up) | −0.3 / +1.7 (prot_only_up) | −1.3 / +2.3 (prot_only_up) | **−2.0 / +2.5** (opposite_rna_down) | −0.6 / +2.5 (prot_only_up) |
| **amt1** | +3.9 / +0.3 (rna_only_up) | +0.3 / +0.3 (both_ns) | +0.6 / +1.0 (prot_only_up) | **−1.5 / +2.3** (opposite_rna_down) | +0.1 / +2.3 (prot_only_up) |
| rpoB | −0.3 / −0.4 (both_ns) | +0.3 / −0.3 (both_ns) | +0.4 / −0.6 (both_ns) | +0.1 / −0.4 (both_ns) | +0.4 / −0.3 (both_ns) |
| atpA | **−4.7 / −0.6 (rna_only_down)** | −0.1 / +0.0 (both_ns) | −1.0 / −0.7 (both_ns) | −1.2 / −0.5 (both_ns) | −0.4 / −0.3 (both_ns) |

Factual observations [KG]:
- The five proposed concordance positive controls (ntcA, glnA, glnB, urtA, amt1) are concordant in only 2 of 25 (5 markers × 5 paired observations) cells — `ntcA` and `glnA` at MED4 axenic d14. They do not behave as concordance positive controls at any MED4 coculture timepoint.
- At MED4 coculture days 31 / 60 / 89, all five canonical markers have protein log2FC > 0 (between +1.0 and +3.3) while their RNA log2FC is ≤ 0 (between −2.6 and +0.6). This produces `prot_only_up` (RNA non-significant + protein significantly up) or `opposite_rna_down` (RNA significantly down + protein significantly up) categories.
- rpoB log2FC is in [−0.6, +0.4] across all 5 observations, both omics; `both_ns` everywhere. rpoB behaves as a clean housekeeping reference for this dataset.
- atpA log2FC_RNA is −4.73 at MED4 axenic d14 (significant_down); protein log2FC = −0.58 (not significant). At MED4 coculture timepoints, atpA is `both_ns`. atpA is not flat at axenic d14, so it is not usable as a uniformly-flat reference; for the discordance analysis it is itself a `rna_only_down` data point at axenic d14.

[interpretation, to be tested at step 5] The MED4 coculture marker pattern is consistent with several non-exclusive mechanisms: (a) post-transcriptional persistence — translation rate or protein stability dominates over transcription at coculture timepoints; (b) mRNA degradation / sampling artefact at coculture; (c) a proteomics-side baseline drift (the exponential control was sampled before coculture-stable steady state). Step 5 testing — including stratification by gene class and stability checks across the 11 paired observations — is required before any one mechanism can be claimed.

Full table: [`data/qc_canonical_markers.csv`](data/qc_canonical_markers.csv); figure: [`figures/qc_canonical_markers.png`](figures/qc_canonical_markers.png).

## Surprises

**S1. The 5 HOT1A3 "protein only" locus tags are all `istB` paralogs.** All 5 are annotated as IS21-like element helper ATPase IstB (locus tags `ACZ81_17130`, `ACZ81_09590`, `ACZ81_07515`, `ACZ81_04425`, `ACZ81_01340`). This is the classic paralog-resolution mismatch between MS proteomics (which assigns shared peptides to multiple genome-redundant copies) and short-read RNA-seq (which discards multi-mappers). The interpretation is methodological, not biological — these 5 transposase paralogs are removed from the paired-discordance pool by the inner-join naturally; they're flagged here so the absence is not mis-read as biology. Logged as F3.

**S2. MED4 has no protein-only genes; HOT1A3 has 5.** Closely related to S1. MED4's smaller genome and lack of recently expanded paralogous gene families produces a clean 100% protein-in-RNA inclusion. HOT1A3 has the istB family, hence the 0.2% protein-only.

**S3. Cross-condition gene-set consistency is perfect.** Both organisms' RNA-seq and proteomics detected sets are *identical* between axenic and coculture (intersection = union for each (organism × omics)). This is what we want — it means the paired-discordance pool is condition-independent, and the cross-condition consistency axis (locked axis 6) compares the same genes' fold-changes, not different gene sets.

**S4. The pooled "days 60+89" rows exist in 4 of the 8 experiments.** All 4 (multi-TP) coculture-side experiments and the HOT1A3 axenic RNA-seq carry a fifth pooled "days 60+89" row alongside the individual day-60 and day-89 rows. The pooled row presumably exists for statistical-power reasons in the original analysis. We drop the pooled rows from paired observations (decision below) — the individual TPs carry the same information at finer resolution and pooling would double-count.

**S5. Per-gene paired log2FC correlations are very low across the entire dataset.** Pearson r between per-gene log2FC_RNA and log2FC_protein is in the range −0.09 to +0.23 across all 11 paired observations. HOT1A3 is essentially decorrelated at every timepoint (Pearson r ≈ 0); MED4 is weakly positive (best at axenic d14 and coculture d60 / d89, r ≈ 0.18–0.23). This is much lower than the 0.4–0.6 typically reported for steady-state paired transcriptome-proteome correlations in bacteria (Maier et al. 2009; Schwanhäusser et al. 2011 [literature, intrinsic — to be cited from KG-resolved publications at step 6]). Two interpretations, both compatible with the data: (a) under acute / chronic N-stress, the *response* (fold-change) layer is post-transcriptionally remodeled — different from the steady-state-abundance correlation that prior literature typically reports; (b) MS detection bias and transcriptome dynamic-range differences add noise that compresses the apparent correlation. Either way, this is the strongest possible motivation for the discordance analysis: there's a lot of decoupling to explain.

**S6. Concordant-significant genes are rare (2.3 % pooled), asymmetric-significant genes are common (29.3 % pooled).** Across all 20 446 paired-gene-observations, fewer than 1 in 40 genes show coordinated transcriptional + protein response in the same direction. About 1 in 3 are caught by exactly one omic. The dominant story in this data is *what one omic sees that the other misses*, not *how mRNA and protein agree*. The rare opposite-direction class (1.6 %) is the most biologically suspicious — opposite-sign FC at matched timepoints implies either active post-transcriptional reversal or a measurement artefact; both warrant scrutiny at step 5.

**S7. Coculture day 18 in MED4 is essentially silent (94.5 % both_ns).** At the earliest time after coculture establishment, almost no genes significantly respond in *either* omic — the early-coculture stress signal is below the per-omic significance threshold for nearly all 1424 paired-detected genes. This is consistent with d18 being a pre-stress baseline within the coculture lineage; meaningful response emerges only by d31 onward.

**S8. HOT1A3 axenic d31 is the most discordant single observation (46.1 % asymmetric).** Nearly half of HOT1A3's paired-detected genes show one-omic-only sig response at the second axenic timepoint. This is the most extreme single-observation discordance in the dataset and likely a strong driving example for step 4 method development.

**S9. The proposed concordance positive controls do not behave as concordance positive controls in MED4 coculture.** ntcA / glnA / glnB / urtA / amt1 — chosen for the panel because they are concordant up in MED4 axenic d14 (where canonical N-stress response is observed) — show `prot_only_up` or `opposite_rna_down` at every MED4 coculture timepoint (d18, d31, d60, d89). RNA log2FC is ≤ 0 while protein log2FC ≥ +1.0 in 18 of 20 (5 markers × 4 coculture TPs) cells. This affects step 3: these genes cannot be used as concordance positive controls for the coculture observations; they are instead **discordance examples** in coculture. Step 3 must either (a) split the control panel by condition (axenic = concordance pos ctrl, coculture = discordance reference) or (b) find different concordance positive controls for coculture, or (c) acknowledge that the dataset may not contain genes that are concordant under coculture and use synthetic / spike-in references if needed.

**S10. atpA log2FC_RNA = −4.73 at MED4 axenic d14.** atpA was proposed as a housekeeping reference; this single value (`rna_only_down` category, RNA significant down + protein not significant) means atpA is not uniformly flat. rpoB is the cleaner housekeeping candidate for this dataset (|log2FC| ≤ 0.6 in both omics across all 5 MED4 observations).

**S11. 30 % of MED4 and 38 % of HOT1A3 paired-pool genes have null `gene_name` in the KG annotations.** Locus_tag is the only safe identifier for filtering and joining (Rule 2). Confirmed: `gene_name` is consistent across omics where present (no rows have RNA-side gene_name set but proteomics-side null, or vice versa) — defensive `fillna` coalesce in `03_paired_fc_overview.py` recovered nothing, but the code is correct against future cases.

## Decisions

**2026-04-27 — F1 hardcoded fix: MED4 axenic RNA-seq → 336 h / day 14.** The KG records `is_time_course=False` for this experiment with no per-TP metadata. Researcher-confirmed pairing (step 1 dialogue): the single RNA-seq contrast corresponds to day 14, the only timepoint where RNA was extractable from axenic MED4 cells. The mapping is asserted in `02_build_paired_observations.py` constants (`F1_TIMEPOINT_HOURS = 336.0`, `F1_TIMEPOINT_LABEL = "day 14"`). This effectively closes F1's data-pipeline aspect for this analysis; the underlying KG metadata gap remains and the suggested KG fix is in `../gaps_and_friction.md` F1.

**2026-04-27 — Drop pooled "days 60+89" rows from paired observations.** The pooled row is a statistical-power summary of d60 and d89 reads. We have d60 and d89 individually and prefer the finer resolution; including the pooled row would double-count the late-timepoint signal.

**2026-04-27 — Exclude the 5 HOT1A3 `istB` paralogs from the paired pool.** They drop out naturally via the inner join; documented as a known paralog-resolution mismatch (S1, F3). No analytical action — just transparency about why the paired pool size is 2221 not 2226.

## Decide-gate checklist

- **Outputs produced** —
  - `scripts/01_select_experiments.py`
  - `scripts/02_build_paired_observations.py`
  - `scripts/qc_gene_coverage.py`
  - `scripts/03_paired_fc_overview.py`
  - `scripts/qc_distributions_and_markers.py`
  - `data/experiments_weissberg2025.csv` (8 rows, compact)
  - `data/experiments_weissberg2025_full.csv` (8 rows, full metadata, debug)
  - `data/timepoints_weissberg2025.csv` (29 TP rows, with growth phase)
  - `data/paired_observations.csv` (11 paired rows)
  - `data/qc_gene_coverage.csv` (4 organism × condition pair rows)
  - `data/paired_fc_long.csv` (20 446 paired-gene-observation rows)
  - `data/paired_fc_summary.csv` (11 per-observation summary rows)
  - `data/qc_fc_distributions.csv` (29 per-experiment per-TP distribution rows)
  - `data/qc_paired_nan_summary.csv` (11 paired-pool NaN summary rows)
  - `data/qc_significance_crosscheck.csv` (29 sig-crosscheck rows)
  - `data/qc_canonical_markers.csv` (35 marker × observation rows)
  - `figures/paired_fc_scatter_grid.png`, `figures/paired_fc_category_bars.png`, `figures/qc_fc_distribution_box.png`, `figures/qc_canonical_markers.png`
  - All `*.log` script logs in `data/`
- **Results presented** — selected-experiments table (8 rows), paired-observations table (11 rows), gene-coverage QC table (4 (org × cond) rows) all inline above; full tables linked.
- **QC gate** —
  - Filter funnel matches expectation (10 → 8 with the 2 dropped explicitly listed).
  - All 8 experiments use `table_scope == "all_detected_genes"`.
  - Per-TP gene counts equal each experiment's `distinct_gene_count` (no per-TP gene-count variation).
  - F1 fix verified by inspecting the row before/after the patch in the script log.
  - Gene-coverage cross-condition consistency check confirmed identical detected sets axenic↔coculture for both organisms.
  - 5 HOT1A3 protein-only locus tags identified as istB paralogs (paralog-resolution mismatch, methodological).
  - QC #1 distributions: medians ∈ [−0.16, +0.07]; RNAseq IQR > proteomics IQR (compression); no scaling pathologies.
  - QC #2 NaN log2FC: zero NaN across all 11 paired observations both omics.
  - QC #3 sig crosscheck: 29 / 29 ✓ row-level vs experiment-node sig counts match (with experiment-level fallback for the single-contrast MED4 axenic RNAseq).
  - QC #4 markers: ntcA / glnA / glnB / urtA / amt1 are concordant_up only at MED4 axenic d14 (2 of 25 cells); at MED4 coculture they are `prot_only_up` or `opposite_rna_down` (S9). rpoB is the clean housekeeping reference; atpA crashes at axenic d14 (S10).
- **Decisions made this step** — F1 hardcoded; drop pooled "days 60+89" rows; exclude 5 istB paralogs (passive — natural inner-join exclusion).
- **Friction logged this step** — F1 closed (data-pipeline aspect); F3 added (HOT1A3 istB paralog ambiguity between MS and RNA-seq); F4 added (initial `qc_gene_coverage.py` used `run_cypher` when the API covered the case — fixed, identical results, generalizable lesson saved as feedback memory); F5 added (editorializing data observations at step 2 before the formal analysis runs — re-wrote notebook + paper text in factual terms, saved as feedback memory).
- **Advance rationale** — selection complete and validated; the paired-discordance pool is well-characterized (1424 MED4 + 2221 HOT1A3 genes); 11 paired (organism × condition × timepoint) observations identified; cross-condition gene-set consistency confirmed; data overview (Pearson correlations + categorical breakdown) and QC (#1–#4) all done; ready for step 3 (analysis framing) — pick a driving example and frame hypothesis, controls, expected outcome, prioritize among the 11 axes. The QC #4 marker observations (S9, S10) directly inform step-3 control-panel selection: the proposed concordance positive controls do not behave concordantly under coculture, so step 3 must redesign the control structure for coculture observations.
