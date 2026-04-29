# Step 3 — Analysis framing

## Context

Step 2 froze a 116-gene candidate set with per-gene-axis facet availability snapshotted at the gene_overview level. Step 3 here takes a "discovery first, framing second" stance: rather than pre-architect the dossier surface, run each per-gene MCP tool in verbose mode against three anchor genes (the 2 driving examples + 1 positive control) and aggregate views across the candidate set, then design the dossier card from what the KG actually returns.

Anchors:
- **Driving example A — PMM0958** (log2fc 5.74, AQ=1, "conserved hypothetical protein (DUF1830)") — rich content; exercises every dossier axis.
- **Driving example B — PMM1898** (log2fc 4.68, AQ=0, "Conserved hypothetical protein", singleton ortholog group) — F1-floor case; exercises empty-axis behavior.
- **Positive control — PMM0246 / ntcA** (log2fc 3.35, AQ=3, "global nitrogen regulatory protein") — well-annotated canonical N-regulator; validates the surface against established biology.

## What I did

```
uv run python 3_analysis_framing/scripts/01_discover_dossier_axes.py
```

Six per-axis MCP discovery calls on the 3 anchors (verbose where applicable) plus four aggregate-view calls across the 116-gene candidate set:

```
gene_overview            (anchors, verbose=True)
gene_ontology_terms      (anchors, organism="MED4", limit=None)
gene_clusters_by_gene    (anchors and candidate_set, verbose=True, limit=None)
gene_derived_metrics     (anchors and candidate_set, verbose=True, limit=None)
gene_homologs            (anchors and candidate_set, verbose=True, limit=None)
gene_response_profile    (anchors and candidate_set, group_by="treatment_type")
```

Outputs in `data/`:
- `discover_anchors.json` — full verbose responses per anchor × axis
- `discover_anchors_summary.csv` — anchor × axis row counts
- `discover_clusters_aggregate.csv` — per-clustering-analysis × candidate-count + curated-description coverage
- `discover_homologs_aggregate.csv` — per-source × specificity-rank × candidate count
- `discover_dms_aggregate.csv` — per-metric_type × candidate count
- `discover_response_aggregate.csv` — per-treatment_type × candidate response classification
- `01_discover_dossier_axes.log`

## Results

### Per-axis anchor findings

| Axis | PMM0958 (driver A) | PMM1898 (driver B) | PMM0246 (ntcA) |
|---|---|---|---|
| Identity (gene_overview) | AQ=1, product DUF1830, 4 annotation_types, 31 expression edges, 17 sig_up cross-study | AQ=0, product "Conserved hypothetical protein", no annotation_types, 7 expression edges, singleton ortholog group | AQ=3, gene_name ntcA, product "global nitrogen regulatory protein", 9 annotation_types, 27 expression edges, 14 sig_up cross-study |
| Ontology (gene_ontology_terms) | 4 terms (cog "Function unknown", cyanorak "Conserved hypothetical proteins", **pfam DUF1830**, tigr "Hypothetical proteins / Conserved") | **0 terms** | 16 terms incl. cog "Transcription", cyanorak "Cellular processes > Adaptation/acclimation > Nitrogen", "Central intermediary metabolism > Nitrogen metabolism", "Regulatory functions > Other"; go_mf "DNA-binding transcription factor activity" |
| Cluster (gene_clusters_by_gene verbose) | 3 memberships: VEG (N/A description), **diel cluster 18 "Enriched for menaquinone and ubiquinone genes"**, **K-means N-stress cluster 1 "Contains nitrogen transport genes such as urtA and cynA. Enriched for transport and binding category"** (early transient, fastest-responding 5-gene cluster) | 1 membership: VEG (N/A description) | 3 memberships: VEG, diel cluster 18, **K-means N-stress cluster 2 "Contains the rapidly-responding subset of hli genes (including hli10) and ntcA. The most highly upregulated hli genes have the strongest NtcA binding sites"** |
| Ortholog (gene_homologs verbose) | 4 groups: cyanorak curated (CK_00001751, 11 members, cross-genus Pro+Syn) + 3 eggNOG levels; consensus_product all "conserved hypothetical protein (DUF1830)"; cog "Function unknown" | 4 groups: curated (CK_00047606, **singleton — 1 member, 1 organism**) + 3 eggNOG levels (Prochloraceae+; 6 members, 5 organisms each, single_genus); no cog_categories on curated; cyanorak_role "Conserved hypothetical proteins" everywhere | 4 groups: curated (CK_00000468, ntcA, 13 members, cross-genus); Prochloraceae (10 members, "CRP family of transcriptional regulators"); Cyanobacteria (19 members, cross-genus to Thermosynechococcus); Bacteria (eggnog COG0664, 129 members across 27 organisms incl. Alteromonas, Marinobacter, etc.); functional_description consistent across groups: "N regulatory functions, transcription, N metabolism" |
| Derived metrics (gene_derived_metrics verbose) | 1 DM: "MED4 vesicle DNA average read coverage (Biller 2014)", value 2210, rank 28, percentile 40, bucket mid, vesicle compartment | **0 DMs** | **0 DMs** |
| Cross-study response profile | groups_responded=[coculture, iron, nitrogen]; nitrogen 8 exps tested 4 up, up_max_log2fc 6.02, up_best_rank 1 | groups_responded=[coculture, nitrogen]; nitrogen 2 of 8 tested, 1 up, up_max_log2fc 4.68 (from this experiment) | groups_responded=[coculture, nitrogen]; nitrogen 8 of 8 tested, 5 up + 1 down across 25 TPs, up_max_log2fc 4.45 |

### Aggregate view across the candidate set (n=116)

**Cluster aggregate** (113/116 have memberships; 238 gene-cluster rows; 26 distinct clusters across 4 clustering analyses):

| analysis_id | cluster_type | candidates | distinct clusters | rows with curated functional_description |
|---|---|---:|---:|---:|
| `med4_expression_level` (1471-2180-14-11) | condition_comparison | 113 | 5 | 32 / 113 |
| `med4_diel_clusters` (journal.pone.0005135) | diel | 80 | 14 | 58 / 80 |
| `med4_kmeans_nstarvation` (msb4100087) | time_course | 37 | 6 | 37 / 37 |
| `med4_phage_transcription_groups` (nature06130) | time_course | 8 | 1 | 8 / 8 |

The K-means N-starvation analysis has 100% curated-description coverage and is the most informative cluster axis for "potential role" interpretation; the expression-level (VEG/HEG/...) analysis is the noisiest (most "N/A" descriptions).

**Ortholog aggregate** (114/116 have at least one group; 437 gene-group rows; 4 specificity ranks all populated):

| source | taxonomic_level | specificity_rank | gene-group rows | distinct candidates |
|---|---|---:|---:|---:|
| cyanorak | curated | 0 | 114 | 114 |
| eggnog | Prochloraceae | 1 | 102 | 102 |
| eggnog | Cyanobacteria | 2 | 110 | 110 |
| eggnog | Bacteria | 3 | 111 | 111 |

**Derived-metric aggregate** (9/116 candidates; 29 rows; all numeric; 7 distinct metric types from 2 publications):

| metric_type | publication | compartment | rows | distinct candidates |
|---|---|---|---:|---:|
| `vesicle_dna_avg_read_coverage` | 10.1126/science.1243457 (Biller 2014) | vesicle | 5 | 5 |
| `damping_ratio` / `diel_amplitude_*_log2` / `peak_time_*_h` / `protein_transcript_lag_h` (6 metrics) | 10.1371/journal.pone.0043432 (Waldbauer 2012) | whole_cell | 4 each (24 total) | 4 |

The 9 candidates with DMs split into 5 with vesicle-payload signal and 4 with diel-rhythmicity signal (no overlap between the 2 sets in the candidate-set DM rows).

**Response-profile aggregate** (116/116 have response data; nitrogen treatment by definition responds for 100%):

| treatment_type | responded | not responded | tested but not responded | not known (no data) |
|---|---:|---:|---:|---:|
| nitrogen | 116 | 0 | 0 | 0 |
| coculture | 42 | 74 | 0 | 0 |
| carbon | 27 | 2 | 0 | 87 |
| iron | 14 | 0 | 102 | 0 |
| light | 9 | 0 | 0 | 107 |
| viral | 9 | 0 | 107 | 0 |
| salt | 6 | 0 | 110 | 0 |
| phosphorus | 4 | 0 | 0 | 112 |

Breadth distribution (n_treatments responded per candidate): 1 → 43 candidates, 2 → 40, 3 → 28, **4 → 5 candidates** (PMM0087, PMM0348, PMM0684 (DUF1651), PMM0731 "uncharacterized conserved membrane protein", PMM0819). The 5 broad-stress responders are a step-5 stretch lead.

### Dossier surface (locked at step 3)

Seven axes per gene (six per-gene + one per-ortholog-group probe), each with its own field list and explicit empty-handling. Card layout (ordering, prominence, collapsing) is deferred to step 4 method design.

| Axis | Fields surfaced | When empty |
|---|---|---|
| **Identity & DE evidence (this experiment)** | locus_tag, gene_name, product, gene_category, annotation_quality, log2fc, padj, rank_up, all_identifiers | Never empty (every candidate has these) |
| **Cluster** | per cluster: cluster_id, cluster_name, cluster_type, cluster_method, member_count, cluster_functional_description, cluster_expression_dynamics, cluster_temporal_pattern, treatment, analysis_name | "no cluster memberships" — 3 candidates (TX50_RS09500/09520/09860) |
| **Ortholog group** (per group, sorted by specificity_rank ascending) | from `gene_homologs(verbose=True)`: group_id, source, taxonomic_level, member_count, organism_count, genera, has_cross_genus_members, consensus_gene_name, consensus_product, group description / functional_description, consensus cyanorak_roles + cog_categories | "no homolog groups" — 2 candidates (TX50_RS09500/09520) |
| **Ortholog group signal probe** (per group, all 4 specificity ranks; cached by group_id) | 2 calls per group: (a) `genes_by_homolog_group(group_ids=[g], limit=None)` → list of member locus_tags; (b) `gene_overview(locus_tags=[members], summary=True)` → envelope only: `total_matching`, `by_organism`, `by_category`, `by_annotation_type`, `has_expression`, `has_significant_expression`, `has_clusters`, `has_derived_metrics`, `has_orthologs`. Surfaces "of the N ortholog members of this group, K have expression / L have clusters / P have DMs / by_organism breakdown". Per-ortholog full mini-cards are deferred to on-demand drill-down at step 5/6 only when the probe surfaces a worth-pursuing signal. | Group has no members beyond the candidate (singleton, e.g. PMM1898's curated CK_00047606) — surface "singleton group, no orthologs to probe" |
| **Cross-study response profile (2-call pattern)** | (a) `gene_response_profile(locus_tags=[g], experiment_ids=[LOCKED_EXP])` for the locked-experiment-only profile; (b) `gene_response_profile(locus_tags=[g])` for the full cross-study profile spanning all MED4 experiments. Per treatment: experiments_total/tested/up/down, timepoints_total/tested/up/down, up_best_rank, up_max_log2fc, down_*. Plus groups_responded / tested_not_responded / not_known | Never empty (the locked experiment is always 1 row) |
| **Ontology** (one row per term, leaves only — see caveats below) | term_id, term_name, ontology_type, level, BRITE tree/tree_code when present | "no ontology terms" — 17 candidates (F1) |
| **Derived metrics** (one row per DM) | derived_metric_id, name, value, value_kind, rank_by_metric, metric_percentile, metric_bucket, metric_type, compartment, treatment_type, publication_doi | "no derived metrics" — 107 candidates |

### Things to be aware of (caveats noted at framing)

1. **`gene_ontology_terms` returns leaves only by default; rollup deferred.** The dossier's ontology axis uses leaf annotations only (`mode="leaf"`, the default). For this candidate set this is sufficient — cyanorak/TIGR `term_name` already encodes ancestors via `>` path syntax (e.g., ntcA's `'Cellular processes > Adaptation/acclimation to atypical conditions and detoxification > Nitrogen'` carries all ancestor levels in the leaf string); COG categories ARE level-0 terms (no rollup possible / meaningful); Pfam has no traditional hierarchy; GO BP/MF/CC and BRITE are too sparse in this candidate set (3+6+5+1 = 15 rows total) for rollup to surface a meaningful pattern. **If a step-5 stretch-(b) question requires ancestor-level aggregation** (e.g. "are candidates concentrated in any broad GO BP ancestor?"), add `gene_ontology_terms(mode="rollup", level=N)` or `cluster_enrichment(ontology=..., level=N)` then — driven by a specific question, not pre-emptively.

2. **`gene_response_profile` aggregates the locked driving experiment as one of N nitrogen experiments.** Without separation, the cross-study `nitrogen` row's counts double-count this analysis's contribution. Mitigated by the 2-call pattern (locked-only + full) and the adjust+flag rules locked below. The reader sees both perspectives and can compare.

3. **Ortholog group signal probe runs at all 4 specificity ranks (curated cyanorak + 3 eggNOG levels).** Cost is modest after caching by group_id (~300–600 unique groups across the candidate set, 2 calls each). The per-ortholog full mini-card design (one card per ortholog member, with cluster + response + ontology of the member) is **deferred to on-demand drill-down at step 5/6** — only when the signal probe surfaces a finding worth pursuing. The probe's `has_expression: 0/N` result on a Bacteria-level group is itself a meaningful answer (KG carries no data on distant orthologs) — honest emptiness, not a surface gap.

### Response profile — adjust + flag rules (locked at step 3)

The full-profile call (b) includes the locked experiment as one of N nitrogen experiments, so a naïve "nitrogen: 4 of 8 experiments up" reading double-counts this analysis's contribution. The dossier card adjusts the cross-study panel as follows:

- **Cleanly subtract** the locked experiment's contribution from the `nitrogen` row's count fields: `experiments_total`, `experiments_tested`, `experiments_up`, `experiments_down`, `timepoints_total`, `timepoints_tested`, `timepoints_up`, `timepoints_down`. Emit a derived sub-row "other N-stress experiments" with these subtracted counts.
- **Flag** (do not subtract) extremes / ranks: `up_best_rank`, `up_max_log2fc`, `down_best_rank`, `down_max_log2fc`. Append a note "may include or be from this analysis's locked experiment (locked log2fc = X)" so the reader can compare with the locked-experiment-only stats.
- **No adjustment** on other treatment_types (coculture, iron, etc.) — the locked experiment is not in those groups.

This split lets the reader see: (i) what this analysis's experiment shows, (ii) what the rest of the KG agrees on across other N-stress experiments and other treatments, (iii) where extremes might be driven by this analysis vs by a separate study.

## Surprises

**S1 — Cluster axis is the strongest "potential role" anchor when curated descriptions are present.** The K-means N-starvation analysis carries 100% curated descriptions; the curated text often *literally* says what the cluster is about (e.g., for PMM0958's K-means cluster 1: "Contains nitrogen transport genes such as urtA and cynA"; for PMM0246's K-means cluster 2: "the most highly upregulated hli genes have the strongest NtcA binding sites"). When such text is populated, it carries the role suggestion directly — no synthesis needed.

**S2 — 5 candidates respond to 4 treatments** (PMM0087, PMM0348, PMM0684 with DUF1651, PMM0731 uncharacterized membrane protein, PMM0819) → broad-stress hypothetical responders. Step-5 stretch lead, not the primary deliverable.

**S3 — `gene_response_profile` aggregates the locked experiment as one of N nitrogen experiments.** Without separating the locked experiment, "nitrogen: 4 of 8 up" double-counts. Locked + adjust-and-flag rules added to the dossier surface.

## Decisions

**2026-04-29 — Drop tier designation; describe the dossier as a flat list of 6 axes with field lists + empty-handling.** The earlier "Tier 1 / 2 / 3" split was layout machinery being designed before the data arrived. Card layout (ordering, prominence, collapsing) is a step-4 method-design decision driven by what produces a readable per-gene card, not a step-3 framing decision.

**2026-04-29 — Driving examples = PMM0958 (rich, A) + PMM1898 (F1 floor, B); positive control = PMM0246 / ntcA.** glnB dropped from the control set as out-of-scope. PMM1898 absorbs the dossier-floor demonstration role that was previously a separate negative control.

**2026-04-29 — Drop `gene_details` from the planned dossier surface, confirmed.** Step 2 already showed gene_details adds nothing for hypothetical-class genes; step 3 anchors confirmed (no sparse fields for any anchor including the rich PMM0246 — even ntcA's gene_details is mostly already-surfaced content).

**2026-04-29 — Cluster axis surfaced as the dossier's strongest "potential role" anchor when curated descriptions are populated.** No representative-member synthesis (retracted earlier); the KG-native `cluster_functional_description` + `cluster_expression_dynamics` + `cluster_temporal_pattern` fields carry the cluster's curated semantics directly.

**2026-04-29 — Response-profile axis uses the 2-call pattern (locked experiment alone + full cross-study), with adjust-and-flag rules on the cross-study panel.** Cleanly-subtractable counts get derived "other-than-this-analysis" sub-rows; extremes / ranks are flagged in place rather than subtracted.

**2026-04-29 — F3 added** for the 2 RefSeq-style locus tags (TX50_RS09500, TX50_RS09520) with no cluster, no homolog, no ontology. Most-stripped floor case in the candidate set.

## Decide-gate checklist

- **Outputs produced** —
  - Script: `3_analysis_framing/scripts/01_discover_dossier_axes.py`
  - Data: `data/discover_anchors.json` (verbose responses per anchor × axis), `data/discover_anchors_summary.csv` (anchor × axis row counts), `data/discover_clusters_aggregate.csv` (per-analysis_id × candidate count + curated-description coverage), `data/discover_homologs_aggregate.csv` (per-source × specificity rank), `data/discover_dms_aggregate.csv` (per-metric_type × candidate count), `data/discover_response_aggregate.csv` (per-treatment classification), `data/01_discover_dossier_axes.log`
  - No new figures at step 3 — this step's outputs are the locked dossier surface in this notebook + the per-axis data files.
  - Reproducibility: `uv run python analyses/<analysis_dir>/3_analysis_framing/scripts/01_discover_dossier_axes.py` from repo root.
- **Results presented** — per-axis anchor findings table, 4 aggregate-view tables, locked dossier surface (6 axes × fields × empty-handling), response-profile adjust-and-flag rules, no-data candidate identification.
- **QC gate** — anchor-level discovery confirmed: PMM0958 returns rich content on all 6 axes; PMM1898 returns ontology empty (F1-floor case confirmed); PMM0246 returns content-rich N-regulator card consistent with established biology in cyanorak roles + cluster co-membership with hli/ntcA. Aggregate views: 113/116 cluster, 114/116 homolog, 99/116 ontology, 9/116 DM, 116/116 response coverage cross-checked against gene_overview's facet-availability flags from step 2.
- **Decisions made this step** — six locked decisions above (tier designation dropped; driving examples; control set; gene_details dropped; cluster axis prioritized; response-profile 2-call pattern; F3 added).
- **Advance rationale** — dossier surface frozen on disk; per-axis fields explicit; empty-handling rules per axis; response-profile self-experiment overlap handled with explicit pattern; ready for step 4 (build the methods module that materializes one card per candidate gene from the locked surface, drive on PMM0958 + PMM1898, validate against PMM0246).
