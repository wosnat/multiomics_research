# Step 4 — Methods (dossier construction module)

## Context

Step 3 locked the 7-axis dossier surface and the 4 anchor genes. Step 4 builds the methodology module that materializes one dossier card per candidate gene from the locked surface. The module is driven on PMM0958 (rich-content driver A) and PMM1898 (F1-floor driver B), validated on PMM0246 (ntcA control), and stress-tested on TX50_RS09500 (F3 floor — RefSeq-only, no clusters / homologs / ontology).

## What I did

```
analyses/2026-04-29-1025-axenic_up_hypotheticals_med4/4_methods/
  dossier.py                        ← methodology utility module
  scripts/
    01_build_dossier_anchors.py     ← driver: build cards for the 4 anchors
    qc_toy_verify_adjust.py          ← hand-verified toy cases for the response-profile adjust math
  data/
    dossiers.json                    ← aggregated structured cards (one entry per gene)
    cards/<locus_tag>.md             ← per-gene markdown cards
    index.md                         ← 1-line summary per gene with link to its card
    group_probe_cache.json           ← cached signal-probe envelopes (12 groups; reused at step 5)
```

**File layout decision (step 4):** per-gene MD files in `cards/` + one aggregated `dossiers.json` for cross-gene work + one `index.md` for navigation. This scales cleanly to step 5's 116-card sweep — readers can link to a specific card without grepping a single 10K-line file, while the aggregated JSON supports cross-gene aggregation and CSV / pivot views at step 5/6.

Run commands (from repo root):

```
uv run python analyses/2026-04-29-1025-axenic_up_hypotheticals_med4/4_methods/scripts/01_build_dossier_anchors.py
uv run python analyses/2026-04-29-1025-axenic_up_hypotheticals_med4/4_methods/scripts/qc_toy_verify_adjust.py
```

### Module surface (`dossier.py`)

Public API (imported by step 5 sweep script):

| Object | Signature | What it does |
|---|---|---|
| `GroupProbeCache` | `GroupProbeCache(cache_path: Path \| None)` with `get_or_compute(group_id, candidate_locus_tag) -> dict \| None` and `save()` | Caches the per-group signal probe envelope by `group_id`. Returns `None` for singleton groups (no orthologs to probe). Persistable JSON. |
| `build_dossier(locus_tag, locked_experiment_id, organism, group_cache) -> dict` | Returns a nested dict with one section per axis | Issues one MCP call per axis (gene_overview verbose, differential_expression_by_gene scoped to locked exp, gene_clusters_by_gene verbose, gene_homologs verbose, gene_ontology_terms, gene_derived_metrics verbose, gene_response_profile × 2). Computes the response-profile adjust+flag inline. |
| `render_card_markdown(card) -> str` | Walks the structured dict and emits a markdown card | Handles every empty-axis branch (no clusters, no homologs, no ontology, no DMs, singleton ortholog group). |

### Per-call cost per dossier card

| Source | Calls (per gene) |
|---|---:|
| gene_overview (verbose, candidate only) | 1 |
| differential_expression_by_gene (locked experiment only) | 1 |
| gene_clusters_by_gene (verbose) | 1 |
| gene_homologs (verbose) | 1 |
| Per-ortholog-group signal probe (cached): genes_by_homolog_group + gene_overview(summary=True) | 0–8 (cache hits reduce on repeat calls) |
| gene_ontology_terms (leaves only, organism-scoped) | 1 |
| gene_derived_metrics (verbose) | 1 |
| gene_response_profile (locked-only + full) | 2 |
| **Total per card (no cache hits)** | **8 + 2 × n_groups** |

For the 4 anchors with cold cache: 4 × 8 + 12 unique groups × 2 = **56 calls total**. Driver completed in under 5 seconds.

## Results

### 1. Toy-data verification of response-profile adjust+flag math

`qc_toy_verify_adjust.py` runs 5 hand-calculated cases against `_adjust_nitrogen_row`:

| Case | Setup | Expected `other_experiments` | Result |
|---|---|---|---|
| T1 (PMM0958-like) | full: 8 N exps, 4 up, 14/25 TPs up; locked: 1/1 up, 1/1 TP up | `{exps_total: 7, tested: 7, up: 3, down: 0, TPs_total: 24, tested: 24, up: 13, down: 0}` | OK |
| T2 (PMM1898-like) | full: 2 of 8 tested, 1 up + 1 down, 1/6 TPs up + 1 down; locked: 1/1 up, 1/1 TP up | `{exps_total: 7, tested: 1, up: 0, down: 1, TPs_total: 24, tested: 5, up: 0, down: 1}` | OK |
| T3 (only-locked-tested edge) | full: 1 of 8 tested, 1 up; locked: 1/1 up | `{all zeros for 'other'}` | OK |
| T4 (full missing) | `full=None` | function returns `None` | OK |
| T5 (locked missing) | `locked=None` | returns full unchanged with `extremes_may_include_locked=False` | OK |

The math is sound for both populated and edge cases.

### 2. Anchor dossier card observations

#### PMM0958 (driver A — rich content)

The card produces a content-rich, interpretable surface: identity (DUF1830, AQ=1, 4 KG identifiers including TX50/UniProt/RefSeq), DE evidence (log2fc=5.74, rank_up=7), 3 cluster memberships (the K-means N-stress cluster 1 carries the curated description "Contains nitrogen transport genes such as urtA and cynA. Enriched for transport and binding category" — early-transient dynamics, fastest-responding 5-gene cluster), 4 ortholog groups (cyanorak curated cross-genus 11-member group + 3 eggNOG levels with 6-member single-genus groups; signal probes show 2–3 of 5 / 10 ortholog members carry expression data, and 1 carries derived-metric data — likely the candidate's own vesicle-payload signal generalizing to a sibling), full response profile (responds in coculture / iron / nitrogen; the cross-study nitrogen panel adjusts to "other 7 N-stress experiments: 7 tested, 3 up across 13/24 TPs" with `up_best_rank=1` flagged as may-include-locked), 4 ontology terms (3 catch-all + 1 informative DUF1830), 1 vesicle-payload DM (rank 28, percentile 40, bucket mid).

[interpretation] The cluster axis carries the strongest "potential role" anchor here: the K-means N-stress cluster co-locates PMM0958 with urtA + cynA + 2 other genes as the fastest-upregulated N-transport-binding subset of the MED4 N-stress response. This is the form the dossier was designed to deliver.

#### PMM1898 (driver B — F1 floor, singleton ortholog group)

Card surfaces the floor case correctly:
- Ontology section: "No ontology terms in any source. (Member of the F1 17-gene no-ontology subset.)" — explicit empty-row.
- Curated cyanorak group: "singleton group (only the candidate gene itself) — no orthologs to probe" — explicit no-probe-possible.
- eggNOG groups (Prochloraceae / Cyanobacteria / Bacteria) all lump PMM1898 into the same 6-member single-genus Prochlorococcus group; signal probes show 2/5 with expression, 2/5 with clusters, 0/5 with DMs across those 5 ortholog members.
- VEG cluster (only cluster membership) carries "N/A" for all curated description fields → rendered as em-dash.

[interpretation] Even at the F1 floor, the dossier surfaces non-trivial cross-organism context via the eggNOG-level ortholog probes. PMM1898's same-family Prochlorococcus orthologs do have some expression data + cluster memberships, even though MED4-PMM1898 itself has no ontology and a curated-cyanorak singleton group. This is more informative than the "nothing on any axis" card a naive design would have produced for an F1-floor candidate.

#### PMM0246 / ntcA (positive control)

Card validates the surface against established N-regulator biology:
- Identity: AQ=3, gene_name=ntcA, 9 annotation_types, function_description "Required for full expression of proteins subject to ammonium repression. Transcriptional activator of genes subject to nitrogen control".
- Cluster: K-means N-stress cluster 2 carries "Contains the rapidly-responding subset of hli genes (including hli10) and ntcA. The most highly upregulated hli genes have the strongest NtcA binding sites." — a literal mechanism-bearing description.
- Ortholog groups: 4 groups; curated cyanorak has 13 members across Pro+Syn (12 in probe); eggNOG Bacteria-level COG0664 has 129 members across 27 organisms incl. Alteromonas, Marinobacter, Pseudomonas, Shewanella, etc., with by_category breakdown 52 Transcription / 34 Stress response / 28 Signal transduction / 8 Transport / 4 Regulatory functions / 2 Energy production. Probe shows 52/128 of distantly-related orthologs have expression data, 17/128 are in clusters, 2/128 carry DMs.
- Ontology: 16 terms across 9 sources, including BRITE "CRP/FNR family", KEGG "ntcA; CRP/FNR family transcriptional regulator, global nitrogen regulator", Pfam "Cyclic nucleotide-binding domain" + "Crp-like helix-turn-helix domain", GO_BP "regulation of nitrogen utilization", GO_MF "DNA-binding transcription factor activity".
- Response profile: nitrogen 8/8 tested, 5 up + 1 down across 13/25 TPs; adjusted to "other 7 N-stress experiments: 7 tested, 4 up + 1 down across 12/24 TPs"; flagged extremes `up_best_rank=4, up_max_log2fc=4.451` may include locked.

The positive control passes — the dossier surfaces well-established N-regulator biology consistently across cluster, ortholog, ontology, and response axes.

#### TX50_RS09500 (F3 floor — RefSeq-only, no cluster / no homolog / no ontology)

Card surfaces 4 explicit empty rows + the response profile (which still contains the locked-experiment row + the cross-study row). Full output:
- Cluster section: "No cluster memberships in any clustering analysis."
- Ortholog section: "No homolog groups in any source. (See `gaps_and_friction.md` F3 if RefSeq-only locus tag.)"
- Ontology section: "No ontology terms in any source. (Member of the F1 17-gene no-ontology subset.)"
- DM section: "No derived metrics."
- Response profile: nitrogen — 1/8 tested in cross-study (only the locked experiment + 1 other), with 6/6 TPs up in the other experiment. After adjustment: "other 7 N-stress experiments: 1 tested, 1 up across 5/24 TPs."

The card runs to completion without error and reports the floor honestly. The response profile axis still surfaces a non-trivial signal (1 other N-experiment showed this gene up across 5/5 timepoints) — even at the F3 floor the dossier captures what is actually there.

### 3. Sample card excerpts (in chat — full cards in `data/anchor_dossiers.md`)

PMM0958's K-means N-stress cluster row (the strongest "potential role" anchor in the candidate set):
```
### `cluster:msb4100087:med4_kmeans_nstarvation:1`
- name: Prochlorococcus MED4 cluster 1 (5 genes)
- type: time_course / method: K-means (K=9)
- member_count: 5
- cluster_functional_description: Contains nitrogen transport genes such as urtA and cynA. Enriched for transport and binding category.
- cluster_expression_dynamics: early transient
- cluster_temporal_pattern: Most rapidly and highly upregulated genes in MED4 during nitrogen starvation, with expression first appearing within 6 hours and peaking early in the time course.
- treatment: N-starvation time course (0, 3, 6, 12, 24, 48h)
- analysis: MED4 K-means N-starvation clusters
```

PMM0246/ntcA's eggNOG Bacteria signal probe (largest probe in the anchor set):
```
- **signal probe** (envelope from `gene_overview(locus_tags=[128 ortholog members], summary=True)`):
    - by_organism: Alteromonas (MarRef v6)=2, Alteromonas macleodii EZ55=3, … Synechococcus elongatus PCC 7942=11, Synechococcus elongatus UTEX 2973=11, …, Thermosynechococcus vestitus BP-1=4
    - by_category: Energy production=2, Regulatory functions=4, Signal transduction=28, Stress response and adaptation=34, Transcription=52, Transport=8
    - by_annotation_type: brite=68, cog_category=128, cyanorak_role=44, ec=3, go_bp=73, go_cc=65, go_mf=117, kegg=68, pfam=127, tigr_role=53
    - has_expression: 52 / 128
    - has_significant_expression: 27 / 128
    - has_clusters: 17 / 128
    - has_derived_metrics: 2 / 128
```

## Surprises

**S1 — Group probe cache hit count = 0 across the 4 anchors.** Each ortholog group_id was unique to its anchor; no cache hits at this scale. Step 5 sweep across all 116 candidates is expected to yield more cache hits, particularly at deeper eggNOG levels (Bacteria-level groups can lump many distantly-related cyanobacterial genes into a single COG-style group).

**S2 — Even F1 / F3 floor candidates produce informative cards via the response-profile axis.** Both PMM1898 (F1 floor, singleton curated group, no ontology) and TX50_RS09500 (F3 floor, no cluster + no homolog + no ontology) get a non-trivial response profile entry beyond the locked experiment. The dossier surface degrades gracefully — the card never collapses to "no data anywhere" because the response-profile axis is universal (every candidate has at least the locked-experiment row).

**S3 — Specificity_rank ordering bug.** Initial sort key `g.get("specificity_rank") or 99` placed rank-0 (curated cyanorak) groups LAST instead of FIRST because Python treats `0` as falsy. Fixed to `g.get("specificity_rank") if g.get("specificity_rank") is not None else 99`. Caught by visual inspection of PMM0958's first card output. Logged here for awareness; not severe enough to be an F-entry friction.

## Decisions

**2026-04-29 — Dossier card structure: nested dict + markdown renderer.** Structured dict per gene saved as JSON for cross-gene aggregation at step 5; markdown rendered separately via `render_card_markdown(card)`. Keeps data and presentation decoupled — step 5 can also produce CSV / per-cluster pivot views from the same JSON.

**2026-04-29 — File layout: per-gene MDs + aggregated JSON + index.md.** Each candidate's card lands in `data/cards/<locus_tag>.md`. The structured cards aggregate into a single `data/dossiers.json` (one entry per candidate). A `data/index.md` provides a 1-line summary per gene (locus_tag + gene_name + log2fc + AQ + product + per-axis count + role) with a link to the card. Scales cleanly to step 5's 116-card sweep without producing one 10K-line MD; the aggregated JSON supports cross-gene aggregation and pivot views at step 5/6.

**2026-04-29 — Group probe cache persisted to disk between runs.** `GroupProbeCache(cache_path=...)` reads at construction, writes at `save()`. Step 5's full-sweep run will benefit from anchor probes computed at step 4 (12 cached entries seeded). Cache is plain JSON with `null` entries for singleton groups, so it's auditable.

**2026-04-29 — Response-profile adjustment math verified by toy cases.** 5 hand-calculated cases pass (`qc_toy_verify_adjust.py`); the math is sound for the populated case (PMM0958-like), the partial-coverage case (PMM1898-like), the only-locked-tested edge, and the missing-input edges (full or locked None).

**2026-04-29 — Empty-axis branches handled per axis.** No clusters → "No cluster memberships in any clustering analysis." No homologs → "No homolog groups in any source. (See `gaps_and_friction.md` F3 if RefSeq-only locus tag.)" No ontology → "No ontology terms in any source. (Member of the F1 17-gene no-ontology subset.)" No DMs → "No derived metrics." Singleton curated cyanorak group → "singleton group (only the candidate gene itself) — no orthologs to probe." All four anchor cards exercise at least one of these branches.

## Decide-gate checklist

- **Outputs produced** —
  - Module: `4_methods/dossier.py` (single-file; ~370 lines; importable from step-5 driver)
  - Driver script: `4_methods/scripts/01_build_dossier_anchors.py` (builds cards for the 4 anchors)
  - QC script: `4_methods/scripts/qc_toy_verify_adjust.py` (5 toy cases for `_adjust_nitrogen_row`)
  - Data: `data/dossiers.json` (aggregated), `data/cards/PMM0958.md`, `data/cards/PMM1898.md`, `data/cards/PMM0246.md`, `data/cards/TX50_RS09500.md` (per-gene cards), `data/index.md` (1-line summary per gene with link), `data/group_probe_cache.json` (cached probes)
  - No figures — step 4 is methodology, not visualization. Step 5 produces aggregate visualizations.
- **Results presented** — anchor card excerpts (PMM0958 K-means cluster, PMM0246 ortholog probe, PMM1898 floor + F3 floor handling) shown inline above and in chat; full cards in `data/anchor_dossiers.md`. Toy-verification 5/5 PASS table shown.
- **QC gate** — toy-data math verification passed (5/5 cases); 4 anchors successfully produced cards spanning rich-content + F1 floor + positive control + F3 floor; ortholog sort bug found via visual inspection and fixed; group probe cache structure validated (12 entries seeded for step 5 reuse).
- **Decisions made this step** — four locked decisions above (data + presentation decoupling; persisted probe cache; toy-verified adjust math; per-axis empty branches).
- **Advance rationale** — methodology module materializes one dossier card per candidate gene from the locked 7-axis surface, validates against rich + floor + positive-control anchors, with toy-verified adjustment math and persisted cache. Ready for step 5 (sweep across all 116 candidates; build the per-cluster / per-treatment aggregate views; produce the README index of the dossier gallery).
