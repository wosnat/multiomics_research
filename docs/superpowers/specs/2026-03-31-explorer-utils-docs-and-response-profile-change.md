# Explorer: utility reference docs + gene_response_profile spec change

Standalone spec for changes to `multiomics_explorer`. Can be loaded
and implemented independently of the research skill changes.

## Context

The N-stress markers analysis in `multiomics_research` uncovered two
issues in the explorer:

1. The analysis utility functions (`response_matrix`,
   `gene_set_compare` in `multiomics_explorer.analysis`) have no
   reference docs. Claude doesn't know they exist or how to use them.
   MCP tool functions have per-tool guides in
   `skills/multiomics-kg-guide/references/tools/` and are loaded as
   MCP resources — the utilities need equivalent docs.

2. `gene_response_profile` conflates "not measured" with "measured
   but not significant" in its `groups_not_known` field. For
   experiments with `significant_only` table scope on full-genome
   platforms, absence of an expression edge means the gene was
   measured and did not respond — not that it's unknown.

---

## 1. Utility function reference docs

### What

Reference docs for Python-only analysis functions that don't have
a corresponding MCP tool.

### Location

`skills/multiomics-kg-guide/references/utils/`

Mirrors `references/tools/`. One file per function, same format
as tool guides (what it does, parameters, return format, examples,
common mistakes).

### MCP resource loading

Register as MCP resources at `docs://utils/{name}`, parallel to
`docs://tools/{tool_name}`.

### Functions to document

Both are in `multiomics_explorer/analysis/expression.py`, exported
via `multiomics_explorer.analysis`:

#### `response_matrix`

```python
from multiomics_explorer.analysis import response_matrix

df = response_matrix(
    genes=["PMM0370", "PMM0920"],
    organism="MED4",
)
```

Parameters:
| Name | Type | Default | Description |
|------|------|---------|-------------|
| genes | list[str] | required | Locus tags to query |
| organism | str \| None | None | Organism filter (fuzzy match) |
| experiment_ids | list[str] \| None | None | Experiment filter (ignored when group_map set) |
| group_map | dict[str, str] \| None | None | experiment_id -> group label for custom grouping |
| conn | GraphConnection \| None | None | Reuse existing Neo4j connection |

Returns: DataFrame with index=locus_tag, one column per treatment
type (or custom group label), plus metadata columns (gene_name,
product, gene_category).

Cell values: `"up"`, `"down"`, `"mixed"`, `"not_responded"`,
`"not_known"`.

Key behavior:
- Without group_map: calls `gene_response_profile(group_by="treatment_type")`
- With group_map: calls `gene_response_profile(group_by="experiment")`,
  re-aggregates by label
- Uses `groups_not_responded` and `groups_not_known` from API result
  to fill non-response cells

#### `gene_set_compare`

```python
from multiomics_explorer.analysis import gene_set_compare

result = gene_set_compare(
    set_a=["PMM0370", "PMM0920"],
    set_b=["PMM0468", "PMM0552"],
    organism="MED4",
    set_a_name="early_responders",
    set_b_name="late_responders",
)
```

Parameters:
| Name | Type | Default | Description |
|------|------|---------|-------------|
| set_a | list[str] | required | First gene set (locus tags) |
| set_b | list[str] | required | Second gene set (locus tags) |
| organism | str \| None | None | Organism filter |
| set_a_name | str | "set_a" | Label for set A in summary |
| set_b_name | str | "set_b" | Label for set B in summary |
| experiment_ids | list[str] \| None | None | Experiment filter |
| group_map | dict[str, str] \| None | None | Custom grouping |
| conn | GraphConnection \| None | None | Reuse connection |

Returns: dict with keys:
- `overlap` — DataFrame of genes in both sets
- `only_a` — DataFrame of genes only in set_a
- `only_b` — DataFrame of genes only in set_b
- `shared_groups` — list[str] of groups where both sets respond
- `divergent_groups` — list[str] of groups where only one set responds
- `summary_per_group` — DataFrame indexed by group with columns:
  {set_a_name}, {set_b_name}, overlap, shared

### Doc format

Each reference file should follow the same structure as the tool
guides:

```markdown
# {function_name}

## What it does
{Brief description}

## Parameters
{Table}

## Response format
{Description of return type and fields}

## Few-shot examples
{2-3 examples showing common usage patterns}

## Chaining patterns
{How this function relates to other API functions}

## Common mistakes
{Mistake/correction pairs}
```

---

## 2. `gene_response_profile` change

### Problem

`groups_not_known` currently lists all treatment groups where no
expression edge exists for a gene. But for experiments using
`significant_only`, `significant_any_timepoint`, or
`filtered_subset` table scopes, absence means the gene was measured
and did not reach significance — especially on microarray platforms
that cover the full genome (~1700 genes for MED4).

The tool reports edge absence correctly but doesn't interpret what
absence means given the experiment's `table_scope`.

### Concrete example

In the N-stress analysis, `gene_response_profile` for ureA
(PMM0965) reports `groups_not_known: ["carbon_stress",
"iron_stress", "light_stress", "phosphorus_stress", "viral",
"salt_stress"]`. But ALL non-N MED4 experiments are microarray
with `significant_only` or `filtered_subset` scope. ureA was
measured in all of them — it just wasn't significant. The correct
interpretation is "tested, not responded," not "unknown."

The workaround required manually checking `list_experiments` for
table_scope per treatment type and reinterpreting the results. This
consumed an entire research iteration.

### Proposed changes

**New field: `groups_tested_not_responded`**

A treatment group belongs here when:
- Experiments exist for that organism + treatment group
- Those experiments use `significant_only`,
  `significant_any_timepoint`, or `filtered_subset` table_scope
- The gene has no expression edge in those experiments
- The platform plausibly covers the gene (implementation decision —
  see "not defined" below)

This field sits between `groups_not_responded` (has a
`not_significant` edge) and `groups_not_known` (truly no data).

**Narrowed `groups_not_known`**

Groups that move to `groups_tested_not_responded` are removed from
`groups_not_known`. After the change, `groups_not_known` contains
only groups where absence is truly uninformative.

**Updated per-result fields:**

| Field | Before | After |
|-------|--------|-------|
| `groups_responded` | list[string] | unchanged |
| `groups_not_responded` | list[string] — has `not_significant` edge | unchanged |
| `groups_tested_not_responded` | (does not exist) | **new** — no edge, but table_scope implies measured |
| `groups_not_known` | list[string] — no edge at all | **narrowed** — only truly uninformative absence |

**Optional: `table_scopes` in response_summary**

Add a `table_scopes` field to each group's response_summary entry:

```json
"carbon_stress": {
  "experiments_total": 4,
  "experiments_tested": 0,
  "table_scopes": ["filtered_subset"],
  ...
}
```

This lets users see at a glance whether absence is informative
without a separate `list_experiments` call.

### Impact on `response_matrix`

When `groups_tested_not_responded` is available, `response_matrix`
should treat it the same as `groups_not_responded` — cell value
`"not_responded"`. Currently, these genes show as `"not_known"`.

### What this spec does NOT define

- How to determine platform coverage (microarray gene count vs
  RNA-seq detection threshold) — implementation decision
- Edge cases for `filtered_subset` with very narrow criteria
  (e.g., "top 34 genes") — may need a coverage threshold
- Whether `groups_tested_not_responded` should appear in the
  `response_summary` detail or only in the top-level list
- Performance implications of cross-referencing table_scope

### Testing

- Use the N-stress markers gene set (21 genes) as a test case:
  ureA, glnB, glsF, PMM1462 should appear in
  `groups_tested_not_responded` for carbon, iron, light, P, viral,
  salt — not in `groups_not_known`
- cynA should remain in `groups_responded` for nitrogen, coculture,
  carbon, iron (it has edges there)
- Genes with `all_detected_genes` experiments that have
  `not_significant` edges should remain in `groups_not_responded`
  (no change)

---

## 3. DataFrame conversion utilities (ADDENDUM)

> Added after initial spec. Does not block items 1-2 above —
> implement whenever convenient.

### Problem

Converting API results to DataFrames requires boilerplate for
functions with nested return fields (lists, dicts). The pattern is
identical every time but easy to get wrong (broken CSV export,
silent object-in-cell issues). The API guide in `multiomics_research`
currently documents this as copy-paste boilerplate — it should be
a utility function instead.

### Proposed functions

Add to `multiomics_explorer/analysis/expression.py` (or a new
`multiomics_explorer/analysis/frames.py` if expression.py is
getting large):

#### `profile_to_dataframe`

Converts a `gene_response_profile` result dict into two DataFrames.

```python
from multiomics_explorer.analysis import profile_to_dataframe

result = gene_response_profile(locus_tags=["PMM0370", "PMM0920"])
genes_df, summary_df = profile_to_dataframe(result)
```

**Returns:** tuple of two DataFrames:

1. `genes_df` — one row per gene, list columns joined as
   semicolon-separated strings, `response_summary` dropped.
   Columns: `locus_tag`, `gene_name`, `product`, `gene_category`,
   `groups_responded`, `groups_not_responded`, `groups_not_known`
   (all as strings, not lists).

2. `summary_df` — one row per gene x group, flattened from
   `response_summary`. Columns: `locus_tag`, `gene_name`, `group`,
   plus all stats fields (`experiments_total`, `experiments_tested`,
   `experiments_up`, `experiments_down`, `timepoints_total`,
   `timepoints_tested`, `timepoints_up`, `timepoints_down`,
   `up_best_rank`, `up_median_rank`, `up_max_log2fc`,
   `down_best_rank`, `down_median_rank`, `down_max_log2fc`).

Both DataFrames are CSV-safe (no nested objects in cells).

#### `results_to_dataframe` (optional)

Generic helper for the simple case. Thin wrapper that warns if
nested fields are detected.

```python
from multiomics_explorer.analysis import results_to_dataframe

result = differential_expression_by_gene(organism="MED4", limit=None)
df = results_to_dataframe(result)
# Same as pd.DataFrame(result["results"]) but warns if any column
# contains list/dict values
```

This is lower priority — `pd.DataFrame(result["results"])` works
for most functions. The value is the warning on nested fields.

### Reference docs

Add `references/utils/profile_to_dataframe.md` following the same
format as the other utility reference docs (Section 1 above).

### Impact on research repo

The Python API guide in `multiomics_research` should reference
`profile_to_dataframe` instead of embedding the flattening
boilerplate. The guide keeps the pre-script checklist, common
mistakes, and import paths — but the conversion pattern becomes
a function call, not copy-paste code.

---

## Implementation order

1. **Utility reference docs** — write `references/utils/response_matrix.md`
   and `references/utils/gene_set_compare.md`, register as MCP
   resources
2. **`gene_response_profile` change** — add `groups_tested_not_responded`,
   narrow `groups_not_known`, optionally add `table_scopes` to
   response_summary
3. **Update `response_matrix`** — handle `groups_tested_not_responded`
   as `"not_responded"`
4. **Update tool guide** — update
   `references/tools/gene_response_profile.md` with new fields and
   updated common mistakes section
5. **DataFrame conversion utilities** (addendum) — implement
   `profile_to_dataframe`, optionally `results_to_dataframe`, add
   reference doc `references/utils/profile_to_dataframe.md`

Steps 1 and 5 are independent of each other. Steps 2-4 are
sequential. Step 5 can be done at any time.
