# Step 5 — Analyze (full sweep + aggregate views)

## Context

Step 4 built and validated the dossier methodology module on 4 anchor genes. Step 5 sweeps it across the full 116-candidate set, produces one card per candidate, and computes aggregate views over the 116-card corpus to identify candidate functional bins (the stretch-(b) question from step 1) and characterize the floor / broad-stress subsets.

## What I did

```
uv run python analyses/2026-04-29-1025-axenic_up_hypotheticals_med4/5_analyze/scripts/01_sweep_dossiers.py
uv run python analyses/2026-04-29-1025-axenic_up_hypotheticals_med4/5_analyze/scripts/02_aggregate_views.py
```

**Sweep:** loaded `2_kg_selection/data/candidate_set.csv` (116 candidates), seeded the group probe cache from step 4 (12 entries from the anchors), iterated `build_dossier` over each candidate sorted by log2fc desc. Per-gene markdown card written to `data/cards/<locus_tag>.md`; structured cards aggregated to `data/dossiers.json`; `data/index.md` regenerated as a 1-line summary per candidate (sorted by log2fc desc, with link to each card).

Sweep timing: **31.8 s for 116 cards (avg 0.27 s/card)**. The cache grew from 12 → 431 unique ortholog group probes. Of the 431 group probes, 27 are singletons (group has only the candidate gene itself; no orthologs to probe).

**Aggregates:** `02_aggregate_views.py` reads `dossiers.json` and emits per-axis fill rate, cluster pivot (per-cluster candidate count + curated description), treatment pivot, breadth distribution, floor-genes table (F1 no-ontology + F3 no-cluster-no-homolog), broad-stress-responders table (breadth ≥ 3), top-FC signal summary (top 10 by log2fc with cluster-anchor takeaway), and ortholog-probe headline summary across the candidate set's groups.

## Results

### Per-axis fill rate (n=116; re-confirms step 2/3)

| Axis | n_with | n_total | fraction |
|---|---:|---:|---:|
| clusters | 113 | 116 | 0.974 |
| ortholog_groups | 114 | 116 | 0.983 |
| ontology_terms | 99 | 116 | 0.853 |
| derived_metrics | 9 | 116 | 0.078 |

Figure: `figures/agg_axis_fill.png`. The cross-study response-profile axis is universal (every candidate has at least the locked experiment row) and is therefore not in the table.

### Cluster pivot — top 15 clusters by candidate count

26 unique clusters touched by candidates across 4 clustering analyses. Saved to `data/agg_clusters_pivot.csv`; figure `figures/agg_clusters_pivot.png`.

| cluster_id | type | members | candidates | share | curated functional description |
|---|---|---:|---:|---:|---|
| `…:med4_expression_level:VEG` | condition_comparison | 1073 | **72** | 0.067 | (N/A — top-quartile RPKM across 10 conditions) |
| `…:med4_kmeans_nstarvation:4` | time_course | 86 | 20 | 0.233 | Enriched for amino acid synthesis genes and genes involved in carbon metabolism… |
| `…:med4_expression_level:MEG` | condition_comparison | 444 | 17 | 0.038 | Enriched for essential genes with homologs in DEG… |
| `…:med4_expression_level:HEG` | condition_comparison | 290 | 15 | 0.052 | Enriched for energy production / translation… |
| `…:med4_diel_clusters:18` | diel | 180 | 14 | 0.078 | Enriched for menaquinone and ubiquinone genes (6/9 genes, p=0.00045) |
| `…:med4_diel_clusters:7` | diel | 121 | **12** | 0.099 | **Enriched for nitrogen metabolism genes (4/8 genes, p=0.087)** |
| `…:med4_diel_clusters:6` | diel | 137 | 10 | 0.073 | Enriched for purine ribonucleotide biosynthesis genes (7/18 genes, p=0.0049) |
| `…:med4_diel_clusters:17` | diel | 136 | 9 | 0.066 | (N/A) |
| `…:med4_phage_transcription_groups:2` | time_course | 25 | 8 | **0.32** | Includes genes involved in RNA degradation and modification such as rne, rnhB… |
| `…:med4_kmeans_nstarvation:5` | time_course | 73 | 7 | 0.096 | Contains two rpoD-like sigma factors that are upregulated during nitrogen stress… |
| `…:med4_diel_clusters:9` | diel | 99 | 6 | 0.061 | Enriched for RNA synthesis, modification, and DNA transcription genes… |
| `…:med4_diel_clusters:5` | diel | 120 | 6 | 0.050 | Enriched for respiratory terminal oxidases (3/3 genes, p=0.019)… |
| `…:med4_diel_clusters:8` | diel | 87 | 5 | 0.057 | Enriched for chaperone genes (7/14 genes, p=0.00023) |
| `…:med4_expression_level:LEG` | condition_comparison | 82 | 5 | 0.061 | (N/A) |
| `…:med4_expression_level:NEG` | condition_comparison | 66 | 4 | 0.061 | (N/A) |

[KG] Plain-numbers reading: 72 / 116 candidates fall in the VEG (top-quartile RPKM) cluster — most candidates are highly expressed at baseline. The K-means N-starvation analysis (msb4100087) carries 100% curated descriptions, and clusters 4 (n=20) and 5 (n=7) are the two largest candidate-touching K-means N-stress clusters. The diel analysis carries N/A for some clusters but 4 of the top-15 by candidate count have curated descriptions including a nitrogen-metabolism-enriched cluster (12 candidates) and a menaquinone/ubiquinone cluster (14 candidates). The phage cluster 2 has the highest candidate-share-of-cluster (8 / 25 = 32 %).

### Treatment pivot — cross-treatment response counts (n=116)

| treatment | candidates_responded | candidates_with_treatment_data | fraction_responded (of n=116) |
|---|---:|---:|---:|
| nitrogen | 116 | 116 | 1.00 |
| coculture | 42 | 116 | 0.36 |
| carbon | 27 | 29 | 0.23 (of 116); 27 / 29 = 0.93 of those with carbon data |
| iron | 14 | 14 | 0.12 |
| light | 9 | 9 | 0.08 |
| viral | 9 | 9 | 0.08 |
| salt | 6 | 6 | 0.05 |
| phosphorus | 4 | 4 | 0.03 |

Figure: `figures/agg_treatments_pivot.png`. The 100 % nitrogen response is by definition (the candidate set IS sig_up in MED4 axenic N-stress). [interpretation] The carbon column is a coverage-limited finding: of the 29 candidates that have any carbon-experiment data in the KG, 27 respond — but the other 87 are `not_known` (no carbon data). Cannot generalize.

### Breadth distribution — number of treatments per candidate

| n_treatments_responded | n_candidates |
|---:|---:|
| 1 | 43 |
| 2 | 40 |
| 3 | 28 |
| 4 | 5 |

Figure: `figures/agg_breadth_distribution.png`. 43 / 116 (37 %) of candidates respond only to nitrogen in the cross-study profile; 73 / 116 (63 %) respond to nitrogen plus at least one other treatment.

### Floor genes (F1 + F3)

`data/agg_floor_genes.csv`. **17 F1** (no ontology terms; AQ=0 in all cases). **2 F3** strict (no cluster AND no homolog) — the two RefSeq-only candidates `TX50_RS09500` and `TX50_RS09520`. A third RefSeq-only candidate `TX50_RS09860` has 1 ortholog group but no clusters and no ontology, so it falls outside the strict F3 definition (qualifies for F1 only).

### Broad-stress responders (breadth ≥ 3)

`data/agg_broad_responders.csv`. **33 candidates** with breadth ≥ 3 (28 % of n=116). The 5 candidates with breadth = 4 (responding to 4 distinct treatments):

| locus_tag | log2fc | AQ | product | groups_responded |
|---|---:|---:|---|---|
| PMM0684 | 4.74 | 1 | conserved hypothetical protein (DUF1651) | carbon, coculture, nitrogen, viral |
| PMM0819 | 4.66 | 1 | conserved hypothetical protein | carbon, coculture, nitrogen, viral |
| PMM0348 | 4.18 | 1 | conserved hypothetical protein family PM-17 | carbon, light, nitrogen, salt |
| PMM0731 | 2.10 | 1 | uncharacterized conserved membrane protein | coculture, iron, nitrogen, salt |
| PMM0087 | 1.87 | 1 | conserved hypothetical protein | carbon, light, nitrogen, salt |

[KG] All 5 are AQ=1 (have description, no named product). Two pairs share the same response axis combination: PMM0684 + PMM0819 both respond to carbon + coculture + nitrogen + viral; PMM0348 + PMM0087 both respond to carbon + light + nitrogen + salt.

### Top 10 by log2fc — split by cluster-anchor availability

`data/agg_top_fc_signals.csv`.

| locus_tag | log2fc | rank_up | n_clusters | first informative cluster (curated description) |
|---|---:|---:|---:|---|
| PMM1828 | 6.72 | 3 | 1 | (none — only VEG with N/A description) |
| PMM1813 | 6.28 | 5 | 1 | (none — only VEG) |
| PMM0958 | 5.74 | 7 | 3 | [diel] cluster 18 menaquinone/ubiquinone; [time_course] K-means N-stress cluster 1 (urtA, cynA) |
| PMM0030 | 5.22 | 9 | 4 | [condition_comparison] MEG (essential genes with DEG homologs) |
| PMM1427 | 4.94 | 12 | 2 | [time_course] K-means N-stress cluster 2 (late sustained) |
| PMM0684 | 4.74 | 14 | 4 | [condition_comparison] HEG (energy / translation); also broad responder breadth=4 |
| PMM1966 | 4.71 | 15 | 1 | (none — only VEG) |
| PMM1898 | 4.68 | 16 | 1 | (none — only VEG; F1 no-ontology) |
| PMM0819 | 4.66 | 17 | 4 | [time_course] K-means cluster 3 (19 genes, late sustained); also broad responder breadth=4 |
| PMM1939 | 4.65 | 18 | 1 | (none — only VEG; F1 no-ontology) |

[KG] 5 of top-10 by log2fc have at least one cluster with a curated functional description; the other 5 are in VEG only (top-quartile RPKM) which carries no curated description. Two of the top-10 are F1 no-ontology candidates (PMM1898, PMM1939). Two of the top-10 are also breadth=4 broad responders (PMM0684, PMM0819).

### Ortholog probe headline (across 431 unique groups touched by candidate set)

From `data/agg_ortholog_probe_summary.csv` derived from `group_probe_cache.json`:
- 431 unique ortholog groups probed (including the 12 seeded from anchors).
- 27 singleton groups (no orthologs beyond the candidate itself; probe returned `None`).
- 404 non-singleton groups: median `has_expression` = 2 of N members; median `has_clusters` = 1; median `has_derived_metrics` = 0.

[interpretation] Most candidate ortholog groups have a small number of orthologs with cross-study expression data and at most one with cluster context; derived-metric data on orthologs is rare. The probe surface gives a per-group "any other organism's KG data?" signal rather than rich cross-organism dossier content per ortholog.

## Surprises

**S1 — VEG cluster contains 72 / 116 candidates (62 %) but carries no curated description.** VEG is the top RPKM-quartile classification across 10 growth conditions — most of the 116 N-stress-up hypotheticals are also "highly expressed at steady state". This finding is descriptive (these genes have been observed at high levels of expression in studies that aren't N-stress) but does not by itself suggest a functional role. The cluster's lack of curated description means VEG cannot serve as a "potential role" anchor for these 72 candidates.

**S2 — A diel-clustering N-metabolism enriched cluster (cluster 7, n=121) contains 12 candidates.** [KG] The cluster's curated description reads "Enriched for nitrogen metabolism genes (4/8 genes, p=0.087)". 12 of our 116 N-stress-up hypotheticals fall in this diel cluster. The p-value (0.087) is borderline by traditional thresholds and the statistic counts only 4 of 8 (50 %) annotated nitrogen-metabolism genes in the cluster. [interpretation] One reading: a sub-set of the candidate hypotheticals participate in a diel-rhythmic N-metabolism program that goes beyond the binary axenic-d14 contrast that defined the candidate set. Alternative readings include: (a) this cluster is large (121 members) and overlap with our 116-set may be statistical (12 / 121 = 10 % is not concentrated by candidate-share), (b) the diel rhythm and the N-stress response could be independent expression programs that happen to share members. Treat as a step-6 caveat-bearing finding rather than a conclusion.

**S3 — Phage cluster 2 (n=25, "RNA degradation and modification") contains 8 candidates = 32 % of the cluster.** Highest candidate-share-of-cluster value in the candidate set. [interpretation] 8 candidate hypotheticals are co-clustered with rne (RNase E) and rnhB (RNase H) under phage induction — they may be RNA-processing factors recruited under stress. Step-6 caveat: this comes from a single phage-induction time-course (Lindell 2007, nature06130); cross-condition replication is absent.

**S4 — 5 of the top-10 by log2fc (PMM1828, PMM1813, PMM1966, PMM1898, PMM1939) sit only in VEG with no other informative cluster anchor.** [interpretation] The very-highest-FC hypotheticals are also the most context-poor on the cluster axis. The dossier card for these relies more on the response-profile and ortholog-probe axes than on the cluster axis. Specifically PMM1898 and PMM1939 are also F1 no-ontology candidates — for them, the dossier surfaces effectively only DE evidence, response profile, and the eggNOG ortholog probes. This is the dossier's weakest content slice.

**S5 — Two pairs of breadth=4 broad responders share the SAME treatment-axis combination.** PMM0684 + PMM0819 both respond to carbon, coculture, nitrogen, viral; PMM0348 + PMM0087 both respond to carbon, light, nitrogen, salt. [interpretation] The matching axes within each pair could indicate (a) shared regulatory network or (b) coincidence given that the candidate set was selected for nitrogen response. Not a conclusion; flagged as a step-6 caveat slot.

## Decisions

**2026-04-29 — Three candidate functional bins surfaced for step-6 evaluation (the stretch-(b) question).** From the cluster-pivot results, three putative bins emerge: (1) the **diel-N-metabolism bin** — 12 candidates in `med4_diel_clusters:7` (curated description: nitrogen metabolism enriched, 4/8 genes); (2) the **phage-RNA-processing bin** — 8 candidates (32 % of cluster) in `med4_phage_transcription_groups:2` (RNA degradation / modification); (3) the **broad-stress-responders bin** — 5 candidates with breadth=4 in the cross-treatment response profile (PMM0684, PMM0819, PMM0348, PMM0731, PMM0087). These three bins are documented as candidate functional bins for step-6 evaluation; whether each represents a coherent functional grouping or coincidental overlap with a broader category is a step-6 question, not a step-5 conclusion.

**2026-04-29 — VEG-only top-FC subset flagged as the dossier's weakest slice for step 6.** 5 of top-10 by log2fc (PMM1828, PMM1813, PMM1966, PMM1898, PMM1939) sit only in VEG with no curated cluster anchor. PMM1898 and PMM1939 also have no ontology (F1). The dossier card content for these 5 reduces to identity + DE + response profile + ortholog probe envelopes. Step 6 evaluation should call out that the dossier's promise to support "potential role" reasoning is weakest precisely for the highest-FC hypotheticals.

**2026-04-29 — Floor genes confirm the F1 / F3 splits from step 2/3.** 17 F1 no-ontology candidates as enumerated in step 2 (matches exactly). F3 strict (no-cluster + no-homolog) = 2 candidates; the third RefSeq-only candidate `TX50_RS09860` has 1 ortholog group so falls outside strict F3 — sub-floor-but-not-floor. F3 entry stands as written; the analysis-level distinction is captured in the floor-genes CSV.

## Decide-gate checklist

- **Outputs produced** —
  - Scripts: `5_analyze/scripts/01_sweep_dossiers.py`, `02_aggregate_views.py`
  - Data: `data/dossiers.json` (116 entries), `data/cards/<locus_tag>.md` (116 cards), `data/index.md`, `data/group_probe_cache.json` (431 entries), `data/agg_axis_fill.csv`, `data/agg_clusters_pivot.csv` (26 clusters), `data/agg_treatments_pivot.csv`, `data/agg_breadth_distribution.csv`, `data/agg_floor_genes.csv` (17 + 2), `data/agg_broad_responders.csv` (33), `data/agg_top_fc_signals.csv`, `data/agg_ortholog_probe_summary.csv`, `data/01_sweep_dossiers.log`, `data/02_aggregate_views.log`
  - Figures: `figures/agg_axis_fill.png`, `figures/agg_clusters_pivot.png`, `figures/agg_treatments_pivot.png`, `figures/agg_breadth_distribution.png`
  - Reproducibility: scripts run via `uv run python analyses/<analysis_dir>/5_analyze/scripts/<name>.py` from repo root.
- **Results presented** — per-axis fill rates, cluster pivot top-15 (with curated descriptions), treatment pivot, breadth distribution, floor-genes table (full 17), 5 broad-responders, top-10-by-log2fc with cluster-anchor takeaway. All shown inline + in chat.
- **QC gate** — sweep envelope `total_matching` cross-checks confirmed against step-2 `gene_overview` envelope (ortholog/cluster/ontology/DM coverage matches expected counts); `not_found` empty for all 116 locus tags; truncation = false on every paginated call (limit=None); cache final size 431 entries grew monotonically through the sweep; toy-verified adjustment math (step 4) applies to every nitrogen row in dossiers.json without negative-count violations.
- **Decisions made this step** — three candidate functional bins identified for step-6 evaluation; VEG-only top-FC subset flagged as weakest dossier slice; F1 / F3 floor counts re-confirmed.
- **Advance rationale** — full 116-card dossier produced and aggregated; per-axis fill rates re-confirmed; three candidate functional bins documented for step-6 evaluation; floor and top-FC weak-slice cases identified; ready for step 6 (evaluate against framing, harvest caveats, consolidate KG-enhancement proposals from F1/F2/F3, finalize paper).
