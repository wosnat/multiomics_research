# Step 2 — KG entries: candidate gene set + per-facet feasibility

## Context

Step 1 locked the question: per-gene KG dossier of upregulated hypotheticals in the MED4 axenic RNA-seq experiment of Weissberg 2025; candidate set = `significant_up` × `annotation_quality ≤ 1`. Step 2 quantifies the candidate set, snapshots per-gene facet availability via `gene_overview`, and smoke-tests each per-gene KG endpoint on the top-3 candidates by log2fc to confirm content (not just counts) before step 3 framing and step 4 dossier construction.

## What I did

Plan was refined against the actual MCP API (read tool docs for `differential_expression_by_gene`, `gene_overview`, `gene_details`, `gene_ontology_terms`, `gene_clusters_by_gene`, `gene_derived_metrics`, `gene_homologs`, `gene_response_profile`). Key API insight: `gene_overview` already exposes `annotation_quality` plus per-gene facet-availability counts (`expression_edge_count`, `closest_ortholog_group_size`, `cluster_membership_count`, `derived_metric_count`, `annotation_types`). One call on the 405 sig_up locus tags both filters the candidate set and snapshots facet availability — no separate "is each facet populated?" probe needed at step 2.

### Interactive MCP probe (sub-step A — schema confirmation)

```
differential_expression_by_gene(experiment_ids=[locked], direction="up", significant_only=True, summary=True)
gene_overview(locus_tags=[5 sample sig_up locus tags], verbose=True)
```

Confirmed `matching_genes = 405`, `top_categories[0] = ('Unknown', 133)`, `annotation_quality` populated as int 0–3, `annotation_types` per-row, `derived_metric_count` per-row.

### Scripts run (sub-steps B–D)

```
uv run python 2_kg_selection/scripts/01_extract_candidate_set.py
uv run python 2_kg_selection/scripts/qc_filter_funnel.py
uv run python 2_kg_selection/scripts/qc_facet_smoke_test.py
```

`01_extract_candidate_set.py` — pulls all sig_up DE rows for the locked experiment via `differential_expression_by_gene(experiment_ids=[locked], direction="up", significant_only=True, limit=None, verbose=True)`; pulls `gene_overview` for the 405 distinct sig_up locus tags via `gene_overview(locus_tags=[…], limit=None, verbose=True)`; filters to `annotation_quality ≤ 1`; outputs `data/sig_up_de_rows.csv`, `data/sig_up_405_overview.csv`, `data/candidate_set.csv` (sorted by log2fc desc).

`qc_filter_funnel.py` — loads the artifacts; produces (a) annotation_quality bar chart across the 405 sig_up; (b) candidate-set log2fc histogram split by AQ=0 vs AQ=1; (c) facet-availability stratified by FC bucket (top quartile vs rest of candidate set), saved to `data/qc_facet_availability_by_fc.csv` and `figures/qc_facet_availability_by_fc.png`.

`qc_facet_smoke_test.py` — picks the top-3 candidates by log2fc; calls each per-gene MCP function ONCE with `locus_tags=[g1, g2, g3]` (batch); writes full responses to `data/qc_facet_smoke_test.json` and a 3-gene × 7-facet summary with `n_rows_for_gene` + `content_sample` per cell to `data/qc_facet_smoke_test_summary.csv`.

## Results

### Filter funnel

| Stage | Count | Source |
|---|---:|---|
| Detected genes (experiment) | 1849 | `list_experiments` envelope `gene_count` |
| Significant up (DESeq2, publication thresholds) | 405 | `differential_expression_by_gene` envelope `matching_genes` |
| Candidate set (annotation_quality ≤ 1) | **116** | `01_extract_candidate_set.py` filter |

### annotation_quality distribution among the 405 sig_up genes

| annotation_quality | n | included in candidate set |
|---:|---:|---|
| 0 (pure hypothetical) | 28 | yes |
| 1 (has description, no named product) | 88 | yes |
| 2 (named product, no full description) | 30 | no |
| 3 (well-annotated) | 259 | no |
| **Total** | **405** | **116 in / 289 out** |

Figure: `figures/qc_annotation_quality_dist.png`.

### Candidate-set log2fc summary

| stat | value |
|---|---:|
| min log2fc | 1.018 |
| median log2fc | 2.061 |
| mean log2fc | 2.451 |
| max log2fc | 6.716 |
| top-quartile cutoff (FC ≥ Q3) | 3.122 |
| top quartile size (FC ≥ 3.122) | 29 |

For context, the all-405-sig_up max |log2fc| is 9.583 (PMM1404, hli, AQ=2 — well-annotated, out of candidate set). The candidate set has the strong-FC tail truncated by the AQ filter at the very top, but a substantial top-quartile remains. Figure: `figures/qc_log2fc_by_annotation_quality.png`.

### gene_category distribution: full sig_up vs candidate vs non-candidate

Saved to `data/qc_gene_category_dist.csv`. Categories sorted by sig_up count.

| gene_category | sig_up (n=405) | candidate AQ ≤ 1 (n=116) | non-candidate AQ ≥ 2 (n=289) | candidate share of sig_up |
|---|---:|---:|---:|---:|
| Unknown | 133 | 103 | 30 | 0.77 |
| Stress response and adaptation | 53 | 6 | 47 | 0.11 |
| Coenzyme metabolism | 37 | 2 | 35 | 0.05 |
| Carbohydrate metabolism | 25 | 0 | 25 | 0.00 |
| Replication and repair | 21 | 0 | 21 | 0.00 |
| Post-translational modification | 19 | 1 | 18 | 0.05 |
| Amino acid metabolism | 18 | 0 | 18 | 0.00 |
| Cell wall and membrane | 16 | 2 | 14 | 0.13 |
| Translation | 14 | 1 | 13 | 0.07 |
| Transport | 14 | 0 | 14 | 0.00 |
| Nucleotide metabolism | 13 | 0 | 13 | 0.00 |
| Transcription | 9 | 1 | 8 | 0.11 |
| Central intermediary metabolism | 7 | 0 | 7 | 0.00 |
| Photosynthesis | 5 | 0 | 5 | 0.00 |
| Lipid metabolism | 4 | 0 | 4 | 0.00 |
| Regulatory functions | 4 | 0 | 4 | 0.00 |
| Cellular processes | 3 | 0 | 3 | 0.00 |
| Energy production | 3 | 0 | 3 | 0.00 |
| Cell cycle and division | 3 | 0 | 3 | 0.00 |
| Secondary metabolites | 2 | 0 | 2 | 0.00 |
| Signal transduction | 2 | 0 | 2 | 0.00 |

103 / 116 candidate genes (89 %) carry gene_category "Unknown"; the remaining 13 are spread across 6 non-Unknown categories: 6 Stress response and adaptation, 2 Coenzyme metabolism, 2 Cell wall and membrane, and 1 each in Post-translational modification, Translation, Transcription. Conversely, 11 of the 21 distinct sig_up gene_categories contribute zero candidates (Carbohydrate, Replication, Amino acid, Transport, Nucleotide, Photosynthesis, Lipid, Regulatory, Cellular, Energy, Cell cycle, Secondary, Signal — all metabolic / well-described biology, all well-annotated).

Within the "Unknown" category, candidate share is 77 % (103 / 133) — most "Unknown"-category sig_up genes are AQ ≤ 1, but 30 are still AQ ≥ 2 (named product despite the catch-all category).

### Top 10 candidates by log2fc (from `data/candidate_set.csv`)

| locus_tag | de_product | AQ | log2fc | annotation_types | OG size | clusters | DMs |
|---|---|---:|---:|---|---:|---:|---:|
| PMM1828 | Hypothetical protein | 1 | 6.72 | cog_category | 6 | 1 | 0 |
| PMM1813 | conserved hypothetical protein | 1 | 6.28 | cog_category, cyanorak_role, tigr_role | 4 | 1 | 0 |
| PMM0958 | conserved hypothetical protein (DUF1830) | 1 | 5.74 | pfam, cog_category, cyanorak_role, tigr_role | 11 | 3 | 1 |
| PMM0030 | conserved hypothetical protein | 0 | 5.22 | tigr_role | 4 | 4 | 0 |
| PMM1427 | conserved hypothetical protein | 1 | 4.94 | cog_category, cyanorak_role, tigr_role | 14 | 2 | 0 |
| PMM0684 | conserved hypothetical protein (DUF1651) | 1 | 4.74 | pfam, cog_category, cyanorak_role, tigr_role | 7 | 4 | 0 |
| PMM1966 | conserved hypothetical protein | 1 | 4.71 | cog_category, tigr_role | 13 | 1 | 0 |
| PMM1898 | Conserved hypothetical protein | 0 | 4.68 | (none) | 1 | 1 | 0 |
| PMM0819 | conserved hypothetical protein | 1 | 4.66 | pfam, cog_category | 4 | 4 | 0 |
| PMM1939 | Conserved hypothetical protein | 0 | 4.65 | (none) | 6 | 1 | 0 |

Two top-10 candidates carry an explicit DUF Pfam hit (PMM0958 = DUF1830, PMM0684 = DUF1651). Two top-10 candidates have **no ontology data at all** (PMM1898, PMM1939 — both AQ=0). PMM1898 has `closest_ortholog_group_size = 1` (a singleton — only this MED4 gene in its curated group).

### Per-facet availability across the candidate set (n=116; from `gene_overview`)

| Facet | Genes with data | Fraction |
|---|---:|---:|
| `annotation_types` (any ontology row, any source) | 99 | 0.85 |
| `expression_edge_count > 0` (any cross-study DE) | 116 | 1.00 |
| `significant_up_count > 0` (sig_up in another study) | 116 | 1.00 |
| `closest_ortholog_group_size > 0` | 114 | 0.98 |
| `cluster_membership_count > 0` | 113 | 0.97 |
| `derived_metric_count > 0` | 9 | 0.08 |

The "any ontology row" headline of 99 / 116 is descriptive and does not by itself imply content-bearing functional knowledge — many of those rows for hypothetical genes carry catch-all term names ("Function unknown", "Conserved hypothetical proteins", "Hypothetical proteins / Conserved"). The next sub-section reports the candidate set's ontology landscape at the term level, directly from the `gene_ontology_terms` envelope. Per-gene term-level classification (which gene is "informative" vs "uninformative") is deferred to step 4/5 dossier construction; at step 2 we report what terms appear and at what frequency, without committing to a per-gene informativeness label.

### Per-source ontology presence on the candidate set (descriptive only; n=116; from `data/qc_ontology_per_source_flags.csv`)

| Source | Genes with at least one term from this source | Fraction |
|---|---:|---:|
| tigr_role | 90 | 0.78 |
| cog_category | 88 | 0.76 |
| cyanorak_role | 64 | 0.55 |
| pfam | 35 | 0.30 |
| go_mf | 6 | 0.05 |
| kegg | 5 | 0.04 |
| go_cc | 5 | 0.04 |
| go_bp | 3 | 0.03 |
| brite | 1 | 0.01 |
| ec | 0 | 0.00 |

Per-source presence is a coverage description, not a content claim. The same gene can have a tigr_role row populated by "Hypothetical proteins / Conserved" (uninformative) and a separate informative term in another source.

### Candidate-set ontology landscape at the term level (n=116; from `gene_ontology_terms` summary envelope)

`gene_ontology_terms(locus_tags=[…116 candidates…], organism="MED4", summary=True, limit=None)` returns the envelope with `by_term` (term × count of genes) and `by_ontology` (source × n_terms × n_genes). Summary stats from the envelope:

- `total_matching` (gene × term rows): **309**
- `total_terms` (distinct terms across the candidate set): **69**
- `total_genes` (candidates with at least one term): **99**
- `terms_per_gene` min / median / max: **1 / 3 / 10**

#### Top 20 terms in the candidate set (full table in `data/qc_candidate_ontology_terms_by_term.csv`, figure `figures/qc_candidate_ontology_terms_top.png`)

| ontology_type | n genes | term_name |
|---|---:|---|
| tigr_role | 84 | Hypothetical proteins / Conserved |
| cog_category | 83 | Function unknown |
| cyanorak_role | 57 | Other > Conserved hypothetical proteins |
| cyanorak_role | 4 | Cellular processes > Adaptation/acclimation to atypical conditions and detoxification > Other |
| cyanorak_role | 4 | Other > Conserved hypothetical domains |
| go_cc | 4 | membrane |
| tigr_role | 4 | Not Found |
| pfam | 3 | Domain of unknown function (DUF3303) |
| pfam | 3 | Protein of unknown function (DUF1651) |
| pfam | 2 | Protein of unknown function (DUF2839) |
| cog_category | 2 | Coenzyme transport and metabolism |
| cyanorak_role | 2 | Cellular processes > Adaptation/acclimation to atypical conditions and detoxification > Phosphorus |
| cog_category | 1 | Cell wall/membrane/envelope biogenesis |
| cog_category | 1 | Translation, ribosomal structure and biogenesis |
| cog_category | 1 | Post-translational modification, protein turnover, chaperones |
| cyanorak_role | 1 | Cellular processes > Adaptation/acclimation to atypical conditions and detoxification > Iron |
| cyanorak_role | 1 | Transcription > Other |
| cyanorak_role | 1 | Cell envelope > Surface structures |
| pfam | 1 | Rho termination factor, N-terminal domain |
| pfam | 1 | Ycf66 protein N-terminus |

The three highest-frequency terms — TIGR "Hypothetical proteins / Conserved" (84 genes), COG "Function unknown" (83 genes), Cyanorak "Other > Conserved hypothetical proteins" (57 genes) — are catch-all "this is a hypothetical" terms. They cover almost the entire candidate set on those three sources. **However**, COG / Cyanorak / TIGR also carry a small tail of informative terms — e.g., COG "Coenzyme transport and metabolism" (2 genes), Cyanorak "Cellular processes > Adaptation/acclimation > Phosphorus" (2 genes), Cyanorak "Cell envelope > Surface structures" (1 gene), TIGR "Transcription / Other" (1 gene). So the source itself is not a reliable proxy for content informativeness; the term name is. The previous source-level "hypclass-only" / "informative" / "no ontology" 3-bucket classification is retracted (`gaps_and_friction.md` F2) — at step 2 we report the term-level landscape and let dossier construction at step 4/5 evaluate informativeness per gene.

Pfam terms in the candidate set are dominated by DUF / "Domain of unknown function" / "Protein of unknown function" assignments, but a handful are non-DUF (Rho termination factor, Ycf66, CopG-like RHH, etc., visible in the full term list). KEGG carries 5 distinct rows, 4 of which are "uncharacterized protein". GO_BP / GO_MF / GO_CC are sparse (3, 6, 5 genes total, with GO_CC concentrated on "membrane").

### Per-source ontology presence stratified by FC bucket (top quartile log2fc ≥ 3.12, n=29; rest n=87)

| Source | top-quartile (n=29) | rest (n=87) |
|---|---:|---:|
| tigr_role | 23 (79 %) | 67 (77 %) |
| cog_category | 22 (76 %) | 66 (76 %) |
| cyanorak_role | 11 (38 %) | 53 (61 %) |
| pfam | 5 (17 %) | 30 (34 %) |
| go_mf, go_cc, kegg, go_bp, brite, ec | small | small |

Pfam coverage drops from 34 % in the rest of the candidate set to 17 % in the top-FC quartile; cyanorak_role drops from 61 % to 38 %. The highest-FC candidates therefore have shallower coverage from the sources that carry actual term content (Pfam in particular). This is consistent with the "highest-FC are most under-annotated" observation but is now phrased as per-source coverage rather than as a per-gene informativeness bucket.

### Per-facet availability stratified by FC bucket (from `data/qc_facet_availability_by_fc.csv`)

Top-quartile = candidate set with `log2fc ≥ 3.122` (n=29); rest = remaining 87.

| facet | top-quartile (n=29) | rest (n=87) |
|---|---:|---:|
| has any ontology | 23 (79 %) | 76 (87 %) |
| has expression | 29 (100 %) | 87 (100 %) |
| has orthologs | 29 (100 %) | 85 (98 %) |
| has clusters | 29 (100 %) | 84 (97 %) |
| has derived metrics | 1 (3 %) | 8 (9 %) |

Top-quartile FC has slightly *lower* ontology coverage than the rest (79 % vs 87 %) — the highest-FC hypotheticals are also the most poorly annotated. Expression / orthologs / clusters are essentially universal in the candidate set; derived metrics are sparse everywhere. Figure: `figures/qc_facet_availability_by_fc.png`.

### Per-gene facet smoke test on top-3 (PMM1828, PMM1813, PMM0958)

Each per-gene MCP function called once with `locus_tags=[PMM1828, PMM1813, PMM0958]` (batch). One row per (facet × gene), `n_rows_for_gene` + `content_sample`. Full responses in `data/qc_facet_smoke_test.json`.

| facet | PMM1828 (log2fc 6.72) | PMM1813 (log2fc 6.28) | PMM0958 (log2fc 5.74) |
|---|---|---|---|
| `gene_overview` | 1 row — AQ=1, product='Hypothetical protein', annotation_types=[cog_category] | 1 row — AQ=1, conserved hypothetical, [cog_category, cyanorak_role, tigr_role] | 1 row — AQ=1, conserved hypothetical (DUF1830), [pfam, cog_category, cyanorak_role, tigr_role] |
| `gene_details` | 1 row — no sparse fields populated | 1 row — no sparse fields populated | 1 row — no sparse fields populated |
| `gene_ontology_terms` | 1 term: cog_category 'Function unknown' | 3 terms: cog 'Function unknown', cyanorak 'Other > Conserved hypothetical proteins', tigr 'Hypothetical proteins / Conserved' | 4 terms incl. **pfam 'Domain of unknown function (DUF1830)'** |
| `gene_clusters_by_gene` | 1 membership: **'Prochlorococcus MED4 cluster VEG (variably expressed)'** | 1 membership: **'Prochlorococcus MED4 cluster VEG (variably expressed)'** | 3 memberships: VEG cluster + 'Prochlorococcus cluster 18 (menaquinone and ubiquinone)' + 'Prochlorococcus MED4 cluster 1 (5 genes)' |
| `gene_derived_metrics` | 0 rows | 0 rows | 1 DM: 'MED4 vesicle DNA average read coverage (Biller 2014 Table S4, top-50 ORFs)' |
| `gene_homologs` | 4 groups (cyanorak curated CK_00039912 + 3 eggNOG levels) | 4 groups (cyanorak CK_00054066 + 3 eggNOG) | 4 groups (cyanorak CK_00001751 + 3 eggNOG) |
| `gene_response_profile` | groups_responded = [coculture, iron, nitrogen]; groups_not_responded = [] | groups_responded = [carbon, coculture, nitrogen]; groups_not_responded = [] | groups_responded = [coculture, iron, nitrogen]; groups_not_responded = [] |

Every facet returns content (not just counts) for at least one of the 3 top candidates, except `gene_details`, which adds nothing beyond `gene_overview` for hypothetical genes (the sparse fields it exposes — `ec_numbers`, `signal_peptide`, `transmembrane_regions`, `cazy_ids`, `transporter_classification`, etc. — are by definition not populated for "Hypothetical protein"-class genes).

## Surprises

**S1 — All 3 top-FC candidates share the same MED4 cluster ('VEG (variably expressed)').** PMM1828, PMM1813, and PMM0958 all belong to a cluster named "Prochlorococcus MED4 cluster VEG (variably expressed)". PMM0958 additionally sits in two other clusters (a Prochlorococcus-wide menaquinone/ubiquinone cluster and a MED4-specific 5-gene cluster). [interpretation] The shared VEG cluster is a candidate axis for step 5 grouping, but the cluster's actual member overlap with the wider candidate set is not yet measured — defer to step 4/5.

**S2 — The candidate-set ontology landscape is dominated by three catch-all "hypothetical-class" terms, but the source-level proxy that produced an earlier "61 hypclass-only" classification was over-reaching.** From `gene_ontology_terms` envelope: 309 gene × term rows, 69 distinct terms across 99 candidates with at least one term (17 candidates have none — see F1). The three highest-frequency terms — TIGR "Hypothetical proteins / Conserved" (84 genes), COG "Function unknown" (83 genes), Cyanorak "Other > Conserved hypothetical proteins" (57 genes) — are catch-all "this is a hypothetical" terms and cover most of the candidate set on those three sources. **However**, those same sources also carry a small tail of informative terms in this candidate set — e.g., COG "Coenzyme transport and metabolism" (2 genes), Cyanorak "Cellular processes > Adaptation/acclimation > Phosphorus / Iron / Other" (~7 genes total), Cyanorak "Cell envelope > Surface structures" (1 gene), TIGR "Transcription / Other" (1 gene). So source ≠ content informativeness; per-gene informativeness must be judged term-by-term, not source-by-source. The earlier 3-bucket source-level classification ("informative / hypclass-only / no ontology") is retracted (F2). At step 2 we publish the by-term landscape; at step 4/5 the dossier construction evaluates informativeness per gene.

**S2b — Two of the top-10 by log2fc have no ontology data at all** (PMM1898, PMM1939, both AQ=0). They sit among the candidate genes the KG can say nothing about on the ontology axis (the 17-gene set logged as F1). PMM1898 also has `closest_ortholog_group_size = 1` (singleton in its curated group).

**S3 — 9 / 116 candidates carry derived-metric data; 1 of those was the smoke-test PMM0958 with a vesicle-payload signal** ("MED4 vesicle DNA average read coverage (Biller 2014 Table S4, top-50 ORFs)"). The 9-gene set is small but the metric type (vesicle DNA payload) is biologically distinct from anything the per-gene tools surface — worth a deliberate look at step 5 even though sparseness limits the cross-gene comparison.

**S4 — `gene_details` adds no incremental information for the candidate set.** All three smoke-tested top-FC candidates returned "no sparse fields populated" (no ec_numbers, no signal_peptide, no transmembrane_regions, no cazy_ids, no transporter_classification). The dossier construction at step 4 can probably skip `gene_details` and rely on `gene_overview` + `gene_ontology_terms` for the gene-side data; gene_details is a deep-dive surface designed for well-annotated genes.

**S5 — DE edge timepoint metadata IS present** (`timepoint: "day 14"`, `timepoint_hours: 336.0` in every sig_up DE row), even though the experiment node carries `is_time_course=false` and `timepoints=null`. This refines the F1 entry from the discordance analysis (`gaps_and_friction.md` of `analyses/2026-04-27-1638-…`): the per-edge label exists; only the experiment-level timepoints rollup is missing. For this analysis (single contrast, no pairing across organisms), no workaround is needed.

## Decisions

**2026-04-29 — Keep AQ ≤ 1 cut; candidate set = 116 genes.** Set size 116 is workable for per-gene dossier construction (above the 10–100 target band but inside the order of magnitude). 28 AQ=0 + 88 AQ=1 split lets us slice the dossier "pure hypothetical" vs "has-description" later if the dossier reveals the two slices behave differently. No need to re-tighten or re-loosen the cut at this step.

**2026-04-29 — Drop `gene_details` from the planned dossier-construction surface.** All three smoke-tested top-FC hypotheticals returned no incremental sparse fields beyond what `gene_overview` already exposes. Step 4 dossier construction will rely on `gene_overview` + `gene_ontology_terms` + `gene_clusters_by_gene` + `gene_derived_metrics` + `gene_homologs` + `gene_response_profile` for per-gene data; if a specific subset of candidates turns out to need the deep-dive (e.g. transporter sub-groups), re-add `gene_details` then.

**2026-04-29 — Retract the source-level "informative / hypclass-only / no ontology" 3-bucket split; report ontology coverage at the term level instead.** Initial attempt classified per-gene by source: terms in `pfam`, `kegg`, `go_*`, `brite`, `ec` were treated as informative; terms in `cog_category`, `cyanorak_role`, `tigr_role` were treated as hypothetical-class-only. Empirical check via `gene_ontology_terms` envelope showed this proxy is over-reaching: those three "hypothetical-class" sources also carry a small tail of genuinely informative terms (e.g. COG "Coenzyme transport and metabolism", Cyanorak "Cellular processes > Adaptation/acclimation > Phosphorus / Iron", TIGR "Transcription / Other"). At step 2 we publish: (a) per-source presence flags as descriptive coverage, (b) the term-level landscape from the envelope (`by_term`, `by_ontology`). Per-gene informativeness assignment is deferred to step 4/5 dossier construction, where each candidate's actual term names are inspected. The candidate_set.csv keeps annotation_types as a pipe-delimited string; the auxiliary file `qc_ontology_per_source_flags.csv` carries per-source booleans for diagnostic queries.

## Decide-gate checklist

- **Outputs produced** —
  - Scripts: `2_kg_selection/scripts/01_extract_candidate_set.py`, `qc_filter_funnel.py`, `qc_facet_smoke_test.py`, `qc_categories_and_ontology.py`
  - Data: `data/sig_up_de_rows.csv` (405 × 23), `data/sig_up_405_overview.csv` (405 × 26), `data/candidate_set.csv` (116 × ~30, sorted by log2fc desc), `data/qc_facet_availability_by_fc.csv`, `data/qc_facet_smoke_test_summary.csv`, `data/qc_facet_smoke_test.json`, `data/qc_gene_category_dist.csv`, `data/qc_ontology_per_source_flags.csv`, `data/qc_candidate_ontology_terms.csv` (309 gene × term rows), `data/qc_candidate_ontology_terms_by_term.csv` (69 unique terms ranked by gene-frequency), `data/qc_candidate_ontology_by_source.csv`, `data/01_extract_candidate_set.log`, `data/qc_facet_smoke_test.log`
  - Figures: `figures/qc_annotation_quality_dist.png`, `figures/qc_log2fc_by_annotation_quality.png` (revised: full sig_up n=405 as background + candidate-set foreground split by AQ), `figures/qc_facet_availability_by_fc.png`, `figures/qc_gene_category_dist.png` (new: stacked bar of all 21 sig_up gene_categories, candidate vs non-candidate share), `figures/qc_candidate_ontology_terms_top.png` (new: top-20 candidate-set terms by gene-frequency, colored by ontology source)
  - Reproducibility: all four scripts run via `uv run python analyses/<analysis_dir>/2_kg_selection/scripts/<name>.py` from repo root; defaults pin the locked experiment id, AQ cut = 1, top-3 by log2fc, top-20 ontology terms.
- **Results presented** — filter funnel, AQ distribution, candidate-set log2fc summary, gene_category distribution, top-10 candidates table, per-facet availability table, FC-stratified availability table, smoke-test 3 × 7 matrix — all shown inline above and in chat.
- **QC gate** —
  - DE call envelope `matching_genes = 405` matches the `list_experiments` `genes_by_status.significant_up = 405` from step 1; `truncated = false` confirms full extraction at `limit=None`.
  - `gene_overview` envelope `total_matching = 405`, `not_found = []` confirms every sig_up locus tag resolves.
  - `annotation_quality` populated for all 405 (no nulls in candidate-set construction; the AQ-null filter is a no-op here).
  - Smoke-test confirmed every facet endpoint returns informative content for at least one of the three top-FC candidates; `gene_details` flagged as redundant for the candidate set (S4).
- **Decisions made this step** — keep AQ ≤ 1 cut (set size 116); drop `gene_details` from planned dossier surface.
- **Advance rationale** — candidate set quantified and frozen on disk; per-facet availability snapshotted across the set; smoke test confirms step 4 has working endpoints to build the dossier from. Ready for step 3 (analysis framing — pick the driving example, lock the dossier layout's hypothesis-bearing axes, define controls).
