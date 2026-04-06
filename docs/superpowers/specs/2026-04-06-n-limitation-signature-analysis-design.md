# N-Limitation Signature Analysis: Weissberg 2025

**Date:** 2026-04-06
**Status:** Approach A spec (B and C in backlog)

## Question

Is Prochlorococcus nitrogen-limited in Weissberg 2025's axenic and coculture experiments, and how does the severity differ over time?

Both conditions use PRO99-lowN medium, so both are N-limited by design. The real question: can we quantify the **degree** of N-limitation molecularly, and show that coculture alleviates it over time (consistent with Alteromonas N-recycling)?

## Three approaches (execute sequentially)

| Approach | Core idea | Status |
|----------|-----------|--------|
| **A: Reference signature scoring** | Build N-limitation gene set from Tolonen 2006 + Read 2017, score Weissberg conditions against it | This spec |
| **B: Pathway-first, then validate** | Define pathways a priori, test enrichment per condition/timepoint, validate against references | Backlog |
| **C: Hybrid reference-anchored pathways** | Combine A + B: reference-validated markers as core, Weissberg-specific as extended, multi-omics integration | Backlog |

## Deliverable

Publication-ready analysis section showing:
- N-limitation gene signature derived from two independent reference studies
- Signature activation scores for axenic vs coculture across all timepoints
- Pathway-grouped interpretation
- Targeted Alteromonas N-recycling supporting evidence

All findings tagged: `[KG]` for data, `[interpretation]` for biological reasoning, `[gap]` for missing data.

## Artifact structure

```
analyses/2026-04-06-1432-n_limitation_signature/
├── exploration/           # Exploration logs (one per iteration)
│   └── 2026-04-06-signature-building.md
├── data/                  # Staged KG extracts (CSV)
│   ├── de_reference_tolonen_ndep.csv
│   ├── de_reference_read_ndep.csv
│   ├── core_signature_genes.csv
│   ├── extended_signature_genes.csv
│   ├── de_weissberg_med4_signature.csv
│   └── de_weissberg_hot1a3_nrecycling.csv
├── scripts/
│   ├── extract_reference_de.py      # Extract DE from reference experiments
│   ├── extract_weissberg_de.py      # Extract DE from Weissberg experiments
│   ├── build_signature.py           # Intersect references → core signature
│   ├── score_signature.py           # Compute activation metrics (uses sig_utils)
│   ├── plot_trajectories.py         # Trajectory figures
│   └── explore_alteromonas.py       # HOT1A3 N-recycling check
├── sig_utils/                 # Reusable scoring building blocks
│   ├── __init__.py
│   ├── metrics.py             # hit_rate(), mean_signed_log2fc(), rank_score(),
│   │                          #   mean_signed_normalized_rank()
│   ├── signature.py           # build_core_signature(), intersect_de_lists(),
│   │                          #   assign_direction(), assign_ranks()
│   └── io.py                  # load_de_csv(), save_scores_csv(),
│                              #   load_signature_csv()
├── results/
│   ├── signature_scores_core.csv
│   ├── signature_scores_extended.csv
│   ├── core_vs_extended_comparison.csv
│   ├── reference_baseline_scores.csv
│   ├── pathway_summary.csv
│   ├── trajectory_rnaseq.png
│   ├── trajectory_proteomics.png
│   └── alteromonas_nrecycling.csv
├── superpowers/           # Snapshot of design spec and plan at analysis start
│   └── 2026-04-06-n-limitation-signature-analysis-design.md
├── README.md
├── methods.md             # Publication-ready methods (living document)
├── gaps_and_friction.md   # Issues log
└── references.md          # DOIs and citations
```

The `sig_utils/` package contains reusable building blocks for signature-based scoring. Designed so that Approaches B and C can import the same metric functions with different gene sets (pathway gene sets, tiered gene sets) without duplicating code. Each function operates on a DataFrame of DE results + a signature definition; no KG or MCP dependency — pure computation on staged data.

The `superpowers/` subdirectory holds a copy of the design spec (and later the implementation plan) as they were when the analysis started. This preserves the original intent even if the specs in `docs/superpowers/specs/` evolve.

## Process: logging as you go

Per the research-methodology skill (artifacts guide), logging happens **during** the analysis, not reconstructed afterward:

- **Exploration logs** (`exploration/YYYY-MM-DD-{topic}.md`): one per research iteration. Each includes: Question, Approach, Findings (tagged `[KG]`/`[interpretation]`/`[gap]`), Assessment (established/preliminary/speculative), Gaps and friction, Next. Written in real time as queries are run and results come in.
- **methods.md**: living document, updated incrementally as findings become established. Publication-ready. Required sections: Research question, Data scope (DOIs, experiment IDs, inclusions/exclusions), Gene selection (filters, counts at each step), Statistical methods (tests, correction, thresholds), Results summary (effect sizes, p-values, file references), Limitations.
- **gaps_and_friction.md**: running log of KG data bugs, KG gaps, MCP friction, skill/methodology friction. Each entry recorded both in the exploration log and here. This feeds back to the multiomics_explorer repo.
- **references.md**: DOIs for all publications, tool citations, KG version, software versions.
- **README.md**: summary of key findings, file index, links between exploration logs.

These are not optional — they are hard deliverables alongside the data and figures.

## Data extraction: MCP vs Python API

- **MCP:** Use for initial scoping, gene lookups, and small queries (<50 results)
- **Python API:** Use for full DE extraction (hundreds of genes x multiple experiments x timepoints). Steps 1 and 3 will exceed MCP limits and must use `extract_*.py` scripts with the multiomics_explorer Python API.
- Staged data files in `data/` serve as the single source of truth for all downstream scripts.

## Data sources

All in the multiomics KG. All MED4 nitrogen treatment experiments.

### Reference studies (signature building)

| Study | Experiment ID | Platform | Design |
|-------|-------------|----------|--------|
| Tolonen 2006 | `10.1038/msb4100087_nitrogen_nitrogen_deprivation_med4_med4_microarray` | Microarray | Time course 0-48h, N-deprivation vs N-replete |
| Tolonen 2006 | `10.1038/msb4100087_growth_medium_growth_on_cyanate_as_med4_microarray` | Microarray | Steady state, cyanate as sole N source |
| Tolonen 2006 | `10.1038/msb4100087_growth_medium_growth_on_urea_as_med4_microarray` | Microarray | Steady state, urea as sole N source |
| Read 2017 | `10.1038/ismej.2017.88_nitrogen_stress_ndepleted_pro99_medium_med4_rnaseq` | RNA-seq | Time course 3-24h, N-depleted vs N-replete |

### Target experiments (Weissberg 2025)

| Experiment | Organism | Platform | Design |
|-----------|----------|----------|--------|
| `..._med4_rnaseq_axenic` | MED4 | RNA-seq | Single point, starvation vs exponential |
| `..._med4_rnaseq_coculture` | MED4 | RNA-seq | Time course d18, d31, d60, d89, d60+89 |
| `..._med4_proteomics_axenic` | MED4 | Proteomics | Time course d14, d31, d89 |
| `..._med4_proteomics_coculture` | MED4 | Proteomics | Time course d18, d31, d60, d89, d60+89 |

Note: RNA-seq axenic is a single timepoint (not a time course), so it gives one score to compare against coculture trajectory.

## Approach A: Reference signature scoring

### Step 1 — Extract DE genes from reference studies

For each reference experiment, query all significant DE genes across timepoints.

**Tolonen 2006 N-deprivation:**
- Use timepoints 6h, 12h, 24h, 48h (skip 0h and 3h — nearly zero DE genes)
- Record per gene: direction, log2FC, rank, rank_up, rank_down, timepoint of peak |log2FC|

**Read 2017 N-depleted:**
- Use timepoints 3h, 12h, 24h
- Same fields recorded
- Note: this dataset is filtered to top 50% of genes by expression level (table_scope = filtered_subset), so absence from this set does not mean non-responsive

**Tolonen 2006 alternative N sources (cyanate, urea):**
- Steady-state experiments, few DE genes (20-43 each)
- Used as supplementary validation, not for building the core signature

### Step 2 — Define core signature by intersection

**Core signature:** genes DE in BOTH Tolonen N-deprivation AND Read N-depleted, with concordant direction (both up or both down).

For each core signature gene, record:
- Direction (up or down)
- **Per-study best rank:** lowest rank within Tolonen and lowest rank within Read, separately. Preserves per-study signal strength for later cross-study comparison.
- **Cross-study best rank:** lowest rank across both studies — used as the primary priority ordering for Approach A scoring.
- **Directional rank (rank_up or rank_down):** use the direction-specific rank matching the gene's concordant direction, not the global rank. Rationale: a gene that is rank 50 overall but rank_up 3 is a top N-response gene that happens to have a modest absolute fold change relative to down-regulated genes. The directional rank better captures its importance within the N-limitation response.
- **Global rank:** record alongside directional rank for completeness, but directional rank is primary.
- log2FC at peak timepoint in each study (for within-study magnitude context only — not for cross-study comparison)
- Timepoint of peak response in each study

**Extended signature:** genes DE in only one study. Scored in parallel with the core signature using the same metrics (hit rate, mean signed log2FC, rank score). Comparing core vs extended results tests whether the core is representative or whether one-study-only genes add information.

Tag each extended gene with its source:
- `tolonen_only_read_absent` — not in Read's filtered subset (top 50% by expression). Absence may be due to low expression, not lack of response. These are the genes affected by limitation 3.
- `tolonen_only_read_ns` — present in Read but not significant. Genuine non-replication.
- `read_only` — significant in Read but not in Tolonen.

These tags feed into Approach C's Tier 2 analysis.

**Expected size:** Tolonen has ~400-500 significant genes across timepoints, Read has ~300-400 (from filtered subset). Core intersection likely 50-150 genes. If <30, flag this and rely more heavily on the extended signature results for interpretation.

### Step 3 — Query Weissberg DE data for signature genes

Pull DE results for all core + extended signature locus tags across all 4 MED4 Weissberg experiments.

Record per gene x experiment x timepoint: log2FC, padj, rank, rank_up, rank_down, expression_status.

Note: Weissberg coculture experiments include a `days 60+89` combined timepoint alongside individual d60 and d89. Use individual timepoints for trajectory plots. Include d60+89 in data tables for completeness but do not plot it (it would duplicate information and confuse the time axis).

### Step 4 — Compute signature activation scores

Three complementary metrics per condition x timepoint:

| Metric | What it measures | Computation |
|--------|-----------------|-------------|
| **Signature hit rate** | Breadth of N-response | Fraction of core signature genes that are significant in Weissberg **in the expected direction** (concordant with reference). Genes significant in the opposite direction are counted separately as "reversed" — biologically interesting, not noise. |
| **Mean signed log2FC** | Magnitude + direction (within-platform only) | See formula below. Do NOT compare log2FC magnitudes between RNA-seq and proteomics or between Weissberg and reference studies. |

**Mean signed log2FC formula:**

Each core signature gene has a known direction from the reference studies (up or down under N-deprivation). We define a sign factor per gene:
- `sign_i = +1` if the gene is UP in the references (e.g., glnA, ntcA)
- `sign_i = -1` if the gene is DOWN in the references (e.g., photosynthesis genes)

Then for each Weissberg condition x timepoint:

```
score = mean( sign_i * log2FC_weissberg_i )   for all detected signature genes i
```

**Interpretation:**
- score > 0 → signature genes are moving in the N-limitation direction (up genes are up, down genes are down)
- score ≈ 0 → no coherent N-limitation signal
- score < 0 → signature genes are moving *opposite* to N-limitation (recovery or reversal)

**Example:** If glnA (reference: up) has log2FC = +3.0 in Weissberg, contribution = +1 * 3.0 = +3.0. If a photosynthesis gene (reference: down) has log2FC = -2.0, contribution = -1 * -2.0 = +2.0. Both contribute positively because both are responding in the expected N-limitation direction.

**Open design decision: log2FC vs rank for this metric.** Mean signed log2FC is meaningful within a single Weissberg experiment (same platform, same pipeline), but rank may be more robust. Compute both and compare.

**Rank-based formula:**

Each signature gene has a directional rank in Weissberg (rank_up if the gene is significant up, rank_down if significant down). We normalize by total genes and apply a concordance sign:

```
concordance_i = +1  if Weissberg direction matches reference direction
concordance_i = -1  if Weissberg direction is opposite to reference direction
concordance_i =  0  if gene is not significant in Weissberg (no directional rank)

normalized_rank_i = 1 - (directional_rank_i / total_genes)
    # top-ranked gene → near 1, lowest-ranked → near 0

rank_score = mean( concordance_i * normalized_rank_i )   for all detected signature genes i
```

**Interpretation:**
- rank_score > 0 → signature genes are highly ranked in the expected N-limitation direction
- rank_score ≈ 0 → no coherent signal (genes not significant, or mixed directions)
- rank_score < 0 → signature genes are highly ranked in the opposite direction (reversal)

**Example:** glnA (reference: up) is significant_up in Weissberg at rank_up 3 out of 1849 genes → concordance = +1, normalized_rank = 1 - 3/1849 = 0.998, contribution = +0.998. A photosynthesis gene (reference: down) is significant_up in Weissberg (opposite!) at rank_up 50 → concordance = -1, normalized_rank = 1 - 50/1849 = 0.973, contribution = -0.973.

Compare both metrics in the exploration phase; use whichever produces clearer separation between axenic and coculture trajectories. Report the decision and reasoning in the exploration log.

| Metric | What it measures | Computation |
|--------|-----------------|-------------|
| **Rank correlation** | Pattern similarity (cross-platform safe) | Spearman rho between reference **directional rank** (rank_up for up genes, rank_down for down genes) and Weissberg directional rank, for genes significant in the expected direction; yields rho + p-value. Rank is the primary cross-platform metric. |

**Reference baselines:** Compute all three metrics (hit rate, mean signed log2FC, rank score) for the reference experiments themselves, per timepoint:
- Tolonen N-deprivation: 6h, 12h, 24h, 48h
- Read N-depleted: 3h, 12h, 24h

This provides an "acute N-deprivation trajectory" as a comparison baseline. Allows questions like: does Weissberg day-18 coculture resemble Tolonen 6h (early, mild) or 24h (severe)? Does axenic exceed all reference timepoints?

Note on circularity: the signature was built from genes significant at *any* timepoint across these studies, so per-timepoint hit rate is NOT 100% — many signature genes are not yet significant at early timepoints (e.g., 3h, 6h). Mean signed log2FC and rank score vary by timepoint regardless. This is valid, not circular.

**Permutation test for mean signed log2FC:** If core signature has >30 genes, shuffle gene labels 1000x, recompute mean signed log2FC each time, derive empirical p-value. If <30 genes, report without p-value (underpowered).

**All metrics computed in `scripts/score_signature.py`**, not in chat. Summary statistics cited from script output.

**Score both core and extended signatures:** Run all metrics for both gene sets in parallel. Present as side-by-side comparison. If core and extended tell the same story → core is representative. If they diverge → the extended adds information and should be discussed (especially the `tolonen_only_read_absent` subset, which may contain genuine N-response genes filtered out by Read's expression threshold).

### Step 5 — Visualize trajectories

- X-axis: timepoint (days for Weissberg, hours for references). Y-axis: activation score.
- Weissberg lines: axenic (dashed) vs coculture (solid).
- Reference baselines: Tolonen and Read per-timepoint scores plotted on a separate x-axis or as horizontal bands/annotations, showing what "acute N-deprivation at Xh" looks like on the same metric scale.
- Separate panels for RNA-seq and proteomics (different scales, different gene coverage).
- One panel per metric, or combined figure.
- Core and extended signature scores as separate line sets (or main + shaded region).

**Expected result:** Axenic shows strong activation at its measured timepoint(s). Coculture shows activation that plateaus or partially reverses at later timepoints. Reference baselines provide anchoring: "coculture day 60 resembles acute N-dep at ~6h" vs "axenic resembles 24h+".

### Step 6 — Pathway annotation (descriptive)

Group core signature genes by functional theme (locus tags are primary identifiers; gene names as labels):
- **N-assimilation:** PMM0920 (glnA), PMM1463 (glnB), PMM0246 (ntcA)
- **N-transport/scavenging:** PMM0970-0974 (urtABCDE), PMM0371 (cynB), PMM0373 (cynS), PMM0964-0965 (ureAB)
- **Photosynthesis:** (genes downregulated under N-stress — identify from signature)
- **Carbon metabolism:** (Tolonen found C/N integration differences — identify from signature)
- **Protein synthesis:** (ribosomal proteins, translation — identify from signature)
- **Stress response:** (general stress markers — identify from signature)

Pathway membership assigned from KG annotations (`gene_category`, `gene_ontology_terms`), not from intrinsic knowledge.

For each theme, report: locus tags, gene names, gene count, direction, which are active in Weissberg axenic vs coculture. This is a table — no statistical test. Sets up Approach B.

### Step 7 — Alteromonas supporting check

Define a targeted HOT1A3 gene set:
- Ammonium transporters
- Amino acid degradation / deaminase enzymes
- Organic nitrogen mineralization
- Peptidase / protease genes

Query DE status in Weissberg HOT1A3 coculture time course experiments (RNA-seq + proteomics).

Report: at which timepoints are these genes upregulated? Does timing align with when MED4's N-limitation score stabilizes?

This is descriptive — no formal test between organisms (different gene spaces).

## Backlog: Approach B additions

On top of Approach A:

- **A priori pathway definitions** from KG ontologies (KEGG, GO BP, CyanoRak roles)
- **Fisher's exact test** per pathway per timepoint: are pathway genes disproportionately DE vs genome background? Direction-aware (separate up/down tests).
- **FDR correction** (Benjamini-Hochberg) across all pathway x timepoint tests
- **Mann-Whitney U** to compare log2FC distributions of pathway members between axenic and coculture
- **Heatmap visualization:** pathways x timepoints, colored by activation level, split by condition

## Backlog: Approach C additions

On top of A + B:

- **Tier 1/Tier 2 partitioning** per pathway: reference-validated vs Weissberg-specific genes
- **Hypergeometric test** per pathway: is Tier 1 overlap with Weissberg greater than chance?
- **Pathway activation trajectory:** mean signed log2FC of Tier 1 genes over time, with Spearman trend test (monotonic worsening vs plateau)
- **Multi-omics concordance:** Spearman correlation between transcript and protein log2FC at matched timepoints for Tier 1 genes
- **Integrated figure:** heatmap with Tier 1/2 distinction + trajectory plots + Alteromonas overlay

## Known limitations

- **Timescale mismatch:** References measure hours (acute stress), Weissberg measures days-to-months (chronic). Some acute-response genes may have returned to baseline by day 18. Approach A will detect this as low hit rate — which is itself informative.
- **Platform differences:** Microarray (Tolonen) vs RNA-seq (Read, Weissberg) vs proteomics (Weissberg). Intersection requirement mitigates this. Per statistical rigor rules: compare direction and rank across platforms, not log2FC magnitude. Mean signed log2FC is valid only within a single experiment; rank correlation is the cross-platform metric.
- **Read 2017 gene coverage:** Filtered to top 50% by expression. Some N-responsive genes may be absent simply because they're lowly expressed. The extended signature (one-study-only genes) partially addresses this.
- **Axenic RNA-seq is a single timepoint:** Cannot track trajectory. Proteomics axenic has 3 timepoints (d14, d31, d89) — use this for temporal comparison.
- **Small sample sizes:** With 4-5 timepoints per condition, trajectory analysis is descriptive, not inferential. Spearman trend test is the most we can do statistically.
