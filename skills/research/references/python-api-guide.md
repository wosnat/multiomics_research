# Python API Guide: Scripting with multiomics_explorer

Use this guide when writing extraction or analysis scripts against the multiomics knowledge graph. The Python API is the right tool for bulk extraction, computed aggregations, and reproducible scripting — use MCP interactively when exploring, use the API when building.

---

## 1. Installation and Import

The `multiomics-explorer` package is declared as a dependency in this project's `pyproject.toml` under the `[analysis]` extras group.

**Install (run once from `multiomics_research` root):**
```bash
uv sync --extra analysis
# or
pip install -e ".[analysis]"
```

**Import MCP tool functions (22 functions):**
```python
from multiomics_explorer import (
    resolve_gene,
    gene_overview,
    gene_details,
    genes_by_function,
    gene_homologs,
    list_organisms,
    list_publications,
    list_experiments,
    list_filter_values,
    search_ontology,
    genes_by_ontology,
    gene_ontology_terms,
    search_homolog_groups,
    genes_by_homolog_group,
    differential_expression_by_gene,
    differential_expression_by_ortholog,
    gene_response_profile,
    list_gene_clusters,
    gene_clusters_by_gene,
    genes_in_cluster,
    kg_schema,
    run_cypher,
)
```

**Import analysis utilities:**
```python
from multiomics_explorer.analysis import (
    response_matrix,
    gene_set_compare,
    to_dataframe,                    # any API result → flat CSV-safe DataFrame
    profile_summary_to_dataframe,    # gene_response_profile → gene × group detail
    experiments_to_dataframe,        # list_experiments → experiment × timepoint
)
```

**Run scripts from the `multiomics_research` project root** — the package is installed in this project's environment. Use `uv run script.py` or `.venv/bin/python script.py`. Do not use `uv run --directory /path/to/multiomics_explorer` — that runs in the wrong environment.

**Neo4j must be running** (same connection as the MCP server, configured via `multiomics_explorer/.env`).

---

## 2. Core Principle: MCP Tools = Python API Functions

The MCP tools and the Python API functions are the same code. Same name, same parameters, same return dict structure.

| Context | Use |
|---------|-----|
| Interactive exploration, ad-hoc questions | MCP tools (via Claude conversation) |
| Scripted extraction, bulk queries | Python API functions |
| Computed aggregations (matrices, comparisons) | Python API + analysis utilities |
| Reproducible pipelines | Python API in scripts |

**Tool documentation** is available as MCP resources at `docs://tools/{tool_name}` — e.g., `docs://tools/gene_response_profile`. Read these to understand available parameters and return fields.

---

## 3. Return Structure

Every function returns a dict with:
- **Envelope fields** (metadata): `total_matching`, `truncated`, `limit`, and function-specific summary fields
- **`results`**: list of dicts, one per row

**Standard extraction pattern:**
```python
import pandas as pd
from multiomics_explorer import list_experiments

result = list_experiments()
df = pd.DataFrame(result["results"])
```

**Retrieve all rows — use `limit=None`:**
```python
result = gene_response_profile(organism="MED4", limit=None)
df = pd.DataFrame(result["results"])
```
No pagination needed. `limit=None` returns all matching rows in one call.

**Always check these envelope fields before proceeding:**
```python
print(result["total_matching"])   # how many rows exist in total
print(result["truncated"])        # True if you got fewer rows than total
```

If `truncated` is `True` and you needed all rows, re-run with `limit=None`.

---

## 4. Handling Nested Fields

Some result fields are lists, dicts, or nested objects. `pd.DataFrame()` stores these as Python objects in cells — silently breaking `.to_csv()`, filtering, and groupby operations.

Use the DataFrame conversion utilities in `multiomics_explorer.analysis` instead of manual flattening:

### `to_dataframe(result)` — universal converter

Pass any API function's return dict. Automatically handles:
- **List columns** → joined with `" | "`
- **Dict columns** (e.g. `genes_by_status`) → inlined as `genes_by_status_significant_up`, etc.
- **Nested structures** (e.g. `response_summary`, `timepoints`) → dropped with a `UserWarning` naming the dedicated function

```python
from multiomics_explorer import genes_by_function
from multiomics_explorer.analysis import to_dataframe

result = genes_by_function("nitrogen")
df = to_dataframe(result)
df.to_csv("output.csv", index=False)  # always safe
```

Works for all API functions. For functions with deeply nested fields, `to_dataframe` gives you the flat gene-level table and warns you to use a dedicated converter for the nested parts.

### `profile_summary_to_dataframe(result)` — gene × treatment detail

For `gene_response_profile` results. Returns one row per gene × treatment group with all stats columns.

```python
from multiomics_explorer import gene_response_profile
from multiomics_explorer.analysis import to_dataframe, profile_summary_to_dataframe

result = gene_response_profile(locus_tags=["PMM0370", "PMM0920"])
genes_df = to_dataframe(result)                   # gene-level flat table
summary_df = profile_summary_to_dataframe(result)  # gene × group detail
```

### `experiments_to_dataframe(result)` — experiment × timepoint

For `list_experiments` results. Expands time-course experiments to one row per timepoint; non-time-course experiments get a single row with NaN timepoint fields.

```python
from multiomics_explorer import list_experiments
from multiomics_explorer.analysis import to_dataframe, experiments_to_dataframe

result = list_experiments(organism="MED4")
exp_df = to_dataframe(result)               # one row per experiment, flat
tp_df = experiments_to_dataframe(result)     # one row per timepoint
```

Reference: `docs://analysis/to_dataframe`

---

## 5. Pre-Script Checklist

Before writing any extraction script, run these checks interactively:

**1. Verify import works:**
```python
from multiomics_explorer import gene_response_profile
```

**2. Test return schema with a small call:**
```python
result = gene_response_profile(locus_tags=["PMM0370"], limit=1)
print(result.keys())                   # envelope fields
print(result["results"][0].keys())     # per-row fields
print(result["results"][0])            # inspect actual values and types
```

**3. Check column names before filtering:**
```python
df = pd.DataFrame(result["results"])
print(df.dtypes)
print(df.columns.tolist())
```
Never guess column names. Always verify against actual output.

**4. Check `total_matching` to understand dataset size:**
```python
print(result["total_matching"])
```
Use this to decide whether `limit=None` is safe or whether you need to filter upstream.

---

## 6. Common Mistakes

| Mistake | Fix |
|---------|-----|
| `uv run --directory /path/to/multiomics_explorer script.py` | Run from `multiomics_research` root — package is installed in this project's environment |
| Guessing column names (`product`, `is_significant`) | Test with `result["results"][0].keys()` first |
| `pd.DataFrame(result["results"])` then `.to_csv()` on nested data | Use `to_dataframe(result)` — handles all flattening automatically (see section 4) |
| Hardcoding the full organism name | Fuzzy matching works: `"MED4"` resolves correctly |
| Writing raw Neo4j or `requests` calls | Use the API functions — they handle connection, auth, and serialization |
| Setting an arbitrarily high `limit=10000` | Use `limit=None` to get all results without guessing |
| Treating `not_known` cells as "not measured" | `not_known` may mean "measured but not significant" — check `groups_tested_not_responded` for genes confirmed tested with no response, and `list_experiments` for the treatment type's `table_scope` |
