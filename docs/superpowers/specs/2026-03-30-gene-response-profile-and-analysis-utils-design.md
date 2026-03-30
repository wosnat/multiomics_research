# Gene response profile MCP tool and analysis utilities

## Problem

During the nitrogen stress analysis of MED4, two recurring friction points emerged:

1. **No cross-experiment gene-level summary in MCP.** The existing `differential_expression_by_gene` tool returns one row per gene x experiment x timepoint. Answering "which stresses does this gene respond to?" requires querying all experiments, collecting all rows, and manually aggregating. For 17 genes across 30 experiments, this exceeds the LLM context window and forces ad hoc jq/python workarounds.

2. **No reusable analysis utilities.** Common operations (build a response matrix, compare gene sets, filter by specificity) were done as one-off scripts during the analysis. These patterns recur across research questions.

## Scope

This spec covers two components across three repos:

| Component | Repo | Purpose |
|-----------|------|---------|
| `gene_response_profile` MCP tool | multiomics_explorer | Cross-experiment gene-level summary for interactive exploration |
| Directional rank on expression edges | biocypher_kg | Precomputed rank among up/down genes per experiment x timepoint |
| Analysis utilities module | multiomics_explorer | Reusable functions that consume the API for common analysis patterns |

Out of scope (separate brainstorm):
- Skill/methodology updates (trigger, state tracking, truncation-to-API switching)
- Time course visualization tools
- GO/functional enrichment

## Component 1: `gene_response_profile` MCP tool

### What it does

Takes a gene list and returns one row per gene summarizing its response across all experiments. Answers: "what does this gene respond to, how strongly, and how consistently?"

### Parameters

| Name | Type | Default | Required | Description |
|------|------|---------|----------|-------------|
| locus_tags | list[string] | — | yes | Gene locus tags |
| organism | string \| null | null | no | Organism name (inferred from genes if omitted) |
| treatment_types | list[string] \| null | null | no | Filter to specific treatment types |
| experiment_ids | list[string] \| null | null | no | Restrict to specific experiments |
| group_by | string | "treatment_type" | no | How to group the response summary. `"treatment_type"` (default): aggregates across experiments by treatment, response_summary keys are treatment type strings. `"experiment"`: one entry per experiment with no cross-experiment aggregation, response_summary keys are experiment IDs. Top-level `treatments_responded` / `treatments_not_responded` / `treatments_not_known` lists are likewise keyed by experiment ID when `group_by="experiment"`. |

### Response format

#### Envelope

```
{
  "organism_name": "Prochlorococcus marinus subsp. pastoris str. CCMP1986",
  "genes_queried": 17,
  "genes_with_response": 15,
  "not_found": [],
  "no_expression": ["PMM1234"],
  "returned": 15,
  "results": [ ... ]
}
```

#### Per-gene result

```json
{
  "locus_tag": "PMM0370",
  "gene_name": "cynA",
  "product": "cyanate transporter",
  "gene_category": "Inorganic ion transport",
  "treatments_responded": ["nitrogen_stress", "coculture"],
  "treatments_not_responded": ["light_stress", "iron_stress", "carbon_stress"],
  "treatments_not_known": ["salt_stress"],
  "response_summary": {
    "nitrogen_stress": {
      "experiments_tested": 4,
      "experiments_up": 3,
      "experiments_down": 0,
      "timepoints_tested": 14,
      "timepoints_up": 8,
      "timepoints_down": 0,
      "up_best_rank": 3,
      "up_median_rank": 8,
      "up_max_log2fc": 5.7,
      "down_best_rank": null,
      "down_median_rank": null,
      "down_max_log2fc": null,
    },
    "coculture": {
      "experiments_tested": 2,
      "experiments_up": 0,
      "experiments_down": 2,
      "timepoints_tested": 6,
      "timepoints_up": 0,
      "timepoints_down": 5,
      "up_best_rank": null,
      "up_median_rank": null,
      "up_max_log2fc": null,
      "down_best_rank": 12,
      "down_median_rank": 15,
      "down_max_log2fc": -13.0,
    }
  }
}
```

### Field definitions

#### Top-level gene fields

| Field | Type | Description |
|-------|------|-------------|
| locus_tag | string | Gene identifier |
| gene_name | string \| null | Gene name (null if unannotated) |
| product | string \| null | Gene product description |
| gene_category | string \| null | Functional category |
| treatments_responded | list[string] | Treatment types where gene is significant in at least one timepoint |
| treatments_not_responded | list[string] | Treatment types where a not_significant expression edge exists (gene was measured, did not change) |
| treatments_not_known | list[string] | Treatment types where no expression edge exists for this gene (gene was not profiled or not reported) |
| response_summary | object | Per-treatment detail (keys = treatment type strings) |

#### Per-treatment fields

| Field | Type | Description |
|-------|------|-------------|
| experiments_tested | int | Experiments for this treatment where the gene was profiled |
| experiments_up | int | Experiments where gene was significant_up in at least one timepoint |
| experiments_down | int | Experiments where gene was significant_down in at least one timepoint |
| timepoints_tested | int | Total timepoints across all experiments for this treatment |
| timepoints_up | int | Timepoints where gene was significant_up |
| timepoints_down | int | Timepoints where gene was significant_down |
| up_best_rank | int \| null | Best (lowest) directional rank across all significant_up timepoints. Null if experiments_up = 0 |
| up_median_rank | float \| null | Median directional rank across all significant_up timepoints. Null if experiments_up = 0 |
| up_max_log2fc | float \| null | Maximum log2FC across all significant_up timepoints (positive). Null if experiments_up = 0 |
| down_best_rank | int \| null | Best (lowest) directional rank across all significant_down timepoints. Null if experiments_down = 0 |
| down_median_rank | float \| null | Median directional rank across all significant_down timepoints. Null if experiments_down = 0 |
| down_max_log2fc | float \| null | Most negative log2FC across all significant_down timepoints (negative). Null if experiments_down = 0 |

### Semantics

- **"responded"** = significant (significant_up or significant_down) in at least one timepoint in at least one experiment for that treatment.
- **"not responded" vs "not known"** is determined by edge existence. If a `not_significant` expression edge exists for the gene in that treatment, it's "not responded" (the gene was measured and did not change). If no edge exists at all, it's "not known" (the gene was not profiled or not reported).
- **Rank and log2FC stats** are computed across all timepoints where the gene was significant in the given direction, across all experiments for that treatment. They are per-timepoint stats aggregated across timepoints, not per-experiment.
- **Directional ranks** (`up_best_rank`, `down_best_rank`, etc.) use the precomputed `rank_up` / `rank_down` properties on expression edges (see Component 2). These rank a gene among all genes going the same direction within an experiment x timepoint.
- **Cross-platform caveat:** rank is comparable across platforms (rank 5 out of 1800 genes means the same thing on microarray or RNA-seq). log2FC magnitudes are not directly comparable across platforms but are included for biological interpretation (e.g., log2FC 5 = 32-fold).

### What the tool does NOT do

- Average log2FC across experiments (different platforms, different baselines)
- Combine p-values across experiments (assumes independence, same test)
- Compute a single effect size across studies
- Classify genes as "N-specific" vs "general stress" (research-specific logic belongs in analysis scripts)
- Report time course dynamics (use `differential_expression_by_gene` to drill into temporal patterns per experiment)

## Component 2: Directional rank on expression edges (biocypher_kg)

### What it adds

Two new precomputed properties on each `AFFECTS_EXPRESSION` edge in the KG:

| Property | Type | Description |
|----------|------|-------------|
| rank_up | int \| null | Rank by \|log2FC\| among significant_up genes within the same experiment x timepoint. Null if the gene is not significant_up. 1 = strongest upregulated gene. |
| rank_down | int \| null | Rank by \|log2FC\| among significant_down genes within the same experiment x timepoint. Null if the gene is not significant_down. 1 = strongest downregulated gene. |

The existing `rank` property (rank among all genes regardless of direction) is retained.

### Computation

During KG build in biocypher_kg, after expression edges are created:

1. For each experiment x timepoint, collect all edges with `expression_status = 'significant_up'`
2. Sort by `|log2fc|` descending
3. Assign `rank_up` = 1, 2, 3, ...
4. Repeat for `expression_status = 'significant_down'` → `rank_down`
5. Non-significant edges get `rank_up = null, rank_down = null`

### Impact on existing tools

- `differential_expression_by_gene` should include `rank_up` and `rank_down` in its per-result fields (alongside existing `rank`)
- `gene_response_profile` uses these for its `up_best_rank`, `up_median_rank`, `down_best_rank`, `down_median_rank` fields

## Component 3: Analysis utilities (multiomics_research)

### Module location

`multiomics_explorer/multiomics_explorer/utils/` — importable Python module alongside the existing API. No CLI wrappers.

All functions return DataFrames. No side effects (no file saving). The caller decides whether and where to persist output.

### Utility 1: `response_matrix`

Builds a gene x treatment matrix showing response status.

```python
from multiomics_explorer.utils import response_matrix

df = response_matrix(
    genes=["PMM0370", "PMM0920", "PMM0468"],
    organism="MED4",
    # optional:
    # experiment_ids=["tolonen_...", "read_..."],  # scope to specific experiments
    # conn=GraphConnection for reuse
)
```

**Returns:** DataFrame with one row per gene, one column per treatment type.

```
             nitrogen_stress  coculture  light_stress  iron_stress  ...
locus_tag
PMM0370      up               down       not_responded not_responded
PMM0920      up               down       not_responded not_responded
PMM0468      down             not_known   down          not_responded
```

Cell values:
- `"up"` — significant_up in at least one timepoint (experiments_up > 0, experiments_down = 0)
- `"down"` — significant_down (experiments_down > 0, experiments_up = 0)
- `"mixed"` — significant in both directions across experiments
- `"not_responded"` — a not_significant edge exists (gene was measured, did not change)
- `"not_known"` — no expression edge exists (gene was not profiled or not reported)

Additional columns: `gene_name`, `product`, `gene_category`.

Optional `group_map: dict[str, str]` parameter maps experiment_id → custom group label. When provided, calls `gene_response_profile` with `group_by="experiment"`, then re-aggregates by the group labels. Matrix columns become the custom group names instead of treatment types.

```python
df = response_matrix(
    genes=["PMM0370", "PMM0920"],
    organism="MED4",
    group_map={
        "tolonen_ndep_...": "axenic_early",
        "weissberg_axenic_...": "axenic_late",
        "weissberg_cocult_...": "coculture",
    }
)
# columns: axenic_early, axenic_late, coculture (instead of treatment types)
```

**Implementation:** Calls `gene_response_profile` API, pivots the response_summary into a matrix.

### Utility 2: `gene_set_compare`

Compares two gene sets by their response profiles.

```python
from multiomics_explorer.utils import gene_set_compare

result = gene_set_compare(
    set_a=["PMM0370", "PMM0920", "PMM0687"],  # e.g., 3h responders
    set_b=["PMM0468", "PMM0552", "PMM0370"],  # e.g., 48h responders
    organism="MED4",
    set_a_name="3h_responders",  # optional labels
    set_b_name="48h_responders",
    # experiment_ids=["tolonen_..."],  # optional: scope to specific experiments
    # group_map={"exp_id": "group_label"},  # optional: custom grouping
)
```

**Returns:** dict with:

```python
{
  "overlap": DataFrame,        # genes in both sets, with response profiles
  "only_a": DataFrame,         # genes only in set A
  "only_b": DataFrame,         # genes only in set B
  "shared_groups": [...],      # groups where both sets respond
  "divergent_groups": [...],   # groups where sets respond differently
  "summary_per_group": DataFrame  # per-group breakdown (see below)
}
```

The `overlap`, `only_a`, `only_b` DataFrames include the response matrix columns (one per group) plus gene metadata. When `group_map` is provided, groups are the custom labels; otherwise they are treatment types.

The `summary_per_group` DataFrame has one row per group, showing how many genes from each set respond. **Note: exact columns are tentative — will refine based on real usage.**

```
                    3h_responders  48h_responders  overlap  shared
group
nitrogen_stress     3              2               1        True
coculture           2              2               1        True
light_stress        0              2               0        False
iron_stress         0              1               0        False
carbon_stress       1              2               1        True
```

Columns:
- `3h_responders` / `48h_responders` (named from `set_a_name` / `set_b_name`): count of genes from that set responding to this group
- `overlap`: count of genes in both sets that respond to this group
- `shared`: True if both sets have at least one gene responding

**Implementation:** Calls `response_matrix` (with `group_map` if provided) for the union of both sets, then partitions.

### Dependency chain

```
gene_response_profile (API)
    └── response_matrix (utility)
            └── gene_set_compare (utility)
```

All three utilities live in multiomics_explorer alongside the API. They call `gene_response_profile` directly as a Python function, not via MCP.

## Implementation order

1. **biocypher_kg:** Add `rank_up` / `rank_down` to expression edges, rebuild KG
2. **multiomics_explorer:** Add `gene_response_profile` API function + MCP tool, update `differential_expression_by_gene` to include directional ranks
3. **multiomics_explorer:** Implement utilities module (`response_matrix`, `gene_set_compare`)

Steps 1 and 2 are sequential (2 depends on 1). Step 3 depends on 2.

## Testing strategy

### MCP tool (multiomics_explorer)

- Unit tests: mock Neo4j, verify response envelope structure, field types, null handling
- Integration tests: real KG queries for known genes (cynA, psaJ) with verified expected values
- Edge cases: gene with no expression data, gene in `significant_only` experiments only, gene with mixed directions

### Utilities (multiomics_explorer)

- Unit tests: mock `gene_response_profile` return values, verify matrix shape, set comparison
- Integration tests: real API calls for the N-stress gene set from the nitrogen analysis (known expected classifications)

### Directional rank (biocypher_kg)

- Verify rank_up + rank_down partition: every significant edge should have exactly one directional rank set
- Verify rank 1 exists for each direction per experiment x timepoint
- Verify consistency: rank_up values are contiguous 1..N where N = count of significant_up edges
