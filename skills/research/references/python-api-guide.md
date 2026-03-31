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

**Import MCP tool functions (19 functions):**
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
    kg_schema,
    run_cypher,
)
```

**Import analysis utilities:**
```python
from multiomics_explorer.analysis import response_matrix, gene_set_compare
```

**Run scripts from the `multiomics_research` project root** — the package is installed in this project's virtual environment (`.venv`). Use `.venv/bin/python script.py`, not `uv run --directory /path/to/multiomics_explorer`.

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

This is the most common source of bugs. Some result fields are lists, dicts, or nested objects. `pd.DataFrame()` will store these as Python objects in cells — this silently breaks `.to_csv()`, filtering, and groupby operations.

**Fields that require special handling:**

| Function | Field | Type | How to handle |
|----------|-------|------|---------------|
| `gene_response_profile` | `groups_responded` | list[string] | Join to string or explode |
| `gene_response_profile` | `groups_not_responded` | list[string] | Join to string or explode |
| `gene_response_profile` | `groups_not_known` | list[string] | Join to string or explode |
| `gene_response_profile` | `response_summary` | nested dict (keys=group names, values=stats dicts) | Extract separately — do not put in main DataFrame |
| `gene_overview` | `annotation_types` | list[string] | Join to string |
| `gene_overview` | `closest_ortholog_genera` | list[string] or None | Join to string, handle None |
| `list_experiments` | `timepoints` | list[dict] or None | Extract separately |
| `list_experiments` | `genes_by_status` | object | Extract separately or flatten |
| `differential_expression_by_gene` | `rows_by_status` | object | Envelope field — not in `results` |
| `differential_expression_by_gene` | `experiments` | list[nested objects] | Envelope field — not in `results` |
| `differential_expression_by_gene` | `top_categories` | list[object] | Envelope field — not in `results` |

### Worked Example: `gene_response_profile`

This function has both list columns and a deeply nested `response_summary` field. Here is the full extraction pattern.

```python
import pandas as pd
from multiomics_explorer import gene_response_profile

result = gene_response_profile(locus_tags=["PMM0370", "PMM0920"])

# 1. Flat gene-level table (join list columns, drop nested response_summary)
genes_df = pd.DataFrame(result["results"])

list_cols = ["groups_responded", "groups_not_responded", "groups_not_known"]
for col in list_cols:
    if col in genes_df.columns:
        genes_df[col] = genes_df[col].apply(
            lambda x: "; ".join(x) if isinstance(x, list) else x
        )

genes_df = genes_df.drop(columns=["response_summary"])
# genes_df is now safe to .to_csv()

# 2. Per-gene x per-treatment summary (extracted from nested response_summary)
summary_rows = []
for gene in result["results"]:
    for group, stats in gene["response_summary"].items():
        summary_rows.append({
            "locus_tag": gene["locus_tag"],
            "gene_name": gene["gene_name"],
            "group": group,
            **stats,
        })

summary_df = pd.DataFrame(summary_rows)
# Each row: one gene x one treatment group, with stats columns from the stats dict
```

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
| `uv run --directory /path/to/multiomics_explorer script.py` | Run from `multiomics_research` root — package is installed in this project's environment. Use `.venv/bin/python script.py`. |
| Guessing column names (`product`, `is_significant`) | Test with `result["results"][0].keys()` first |
| `pd.DataFrame(result["results"])` then `.to_csv()` on nested data | Flatten list/dict columns first (see section 4) |
| Hardcoding the full organism name | Fuzzy matching works: `"MED4"` resolves correctly |
| Writing raw Neo4j or `requests` calls | Use the API functions — they handle connection, auth, and serialization |
| Setting an arbitrarily high `limit=10000` | Use `limit=None` to get all results without guessing |
| Treating `not_known` cells as "not measured" | `not_known` may mean "measured but not significant" — check `list_experiments` for the treatment type's `table_scope` to understand what was covered |
