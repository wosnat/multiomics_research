# Explorer change: DataFrame conversion utilities

**Date:** 2026-03-31
**Explorer commits:** `8f9a29f`..`641a3b1` (on main)

## What changed

Three new functions in `multiomics_explorer.analysis`:

```python
from multiomics_explorer.analysis import (
    to_dataframe,                    # any API result → flat CSV-safe DataFrame
    profile_summary_to_dataframe,    # gene_response_profile → gene × group detail
    experiments_to_dataframe,        # list_experiments → experiment × timepoint
)
```

### `to_dataframe(result)`

Universal converter. Pass any API function's return dict:

```python
from multiomics_explorer import genes_by_function
from multiomics_explorer.analysis import to_dataframe

result = genes_by_function("nitrogen")
df = to_dataframe(result)
df.to_csv("output.csv", index=False)  # always safe
```

What it does automatically (no hardcoded column names):
- List columns → joined with `" | "`
- Dict columns (e.g. `genes_by_status`) → inlined as `genes_by_status_significant_up`, etc.
- Nested structures (e.g. `response_summary`, `timepoints`) → dropped with `UserWarning` naming the dedicated function

### `profile_summary_to_dataframe(result)`

For `gene_response_profile` results only. Returns one row per
gene × treatment group with all stats columns.

```python
result = gene_response_profile(locus_tags=["PMM0370", "PMM0920"])
genes_df = to_dataframe(result)              # gene-level flat table
summary_df = profile_summary_to_dataframe(result)  # gene × group detail
```

### `experiments_to_dataframe(result)`

For `list_experiments` results only. Expands time-course experiments
to one row per timepoint; non-time-course get single row with NaN
timepoint fields.

```python
result = list_experiments(organism="MED4")
tp_df = experiments_to_dataframe(result)
```

## Impact on python-api-guide

Section 4 ("Handling Nested Fields") should be replaced. The 60-line
boilerplate pattern is now three lines:

**Before:**
```python
# 20+ lines of manual list joining, response_summary extraction...
```

**After:**
```python
df = to_dataframe(result)
# For gene_response_profile detail:
summary_df = profile_summary_to_dataframe(result)
# For list_experiments timepoints:
tp_df = experiments_to_dataframe(result)
```

The "Note: `profile_to_dataframe` will handle this conversion
automatically once available" stub can be removed — it's done.

The common mistake "Treating `not_known` cells as not measured" is
now partially addressed by `groups_tested_not_responded` in the
API (shipped earlier today).

## MCP resource

Reference doc available at `docs://analysis/to_dataframe`.

## Action needed

1. Update `skills/research/references/python-api-guide.md` section 4
   to use the new utilities instead of manual boilerplate
2. Run `uv sync --extra analysis` to pick up the new explorer version
