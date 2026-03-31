# Skill hardening, Python API guide, and gene_response_profile spec change

## Problem

The N-stress markers analysis (second iteration) revealed that the
research skill redesign from 2026-03-31 did not prevent the same
process failures that motivated it:

1. **Claude jumped straight into analysis without planning or
   artifacts.** No methods.md, no exploration log, no directory
   structure. The skill describes a checklist, but Claude treats
   checklist items as suggestions. The "orientation outputs" section
   is buried after the orientation queries, so Claude starts querying
   before noticing it should have created files first.

2. **Claude couldn't use the Python API.** Didn't know how to import
   `multiomics_explorer`, fought with column names, guessed at return
   structure, fell back to ad-hoc scripts. The tool guides have a
   3-line `## Package import equivalent` section but no guidance on
   scripting patterns, return structure handling, or common pitfalls
   (list fields in DataFrames, pagination).

3. **Claude didn't check prior work.** The user deliberately re-ran
   the same analysis to test whether skill improvements stuck. Claude
   started from scratch instead of checking `analyses/` for existing
   work on the same topic.

4. **`gene_response_profile` conflates "not measured" with "measured
   but not significant."** The `groups_not_known` field reports
   treatment groups where no expression edge exists, but for
   `significant_only` experiments, absence means the gene was measured
   and not significant — evidence of non-response, not unknown.

### Relationship to prior specs

- **2026-03-31 iterative redesign spec:** Introduced the 2-stage
  workflow (Orientation + Research loop). This spec adds Phase 0
  (pre-flight) as a hard gate before orientation, strengthens gate
  enforcement, and adds the prior-work check. The Stage 1/2 structure
  is unchanged.
- **2026-03-30 gene-response-profile spec:** Defined the tool's
  parameters and response format. This spec proposes an additive
  change: a new `groups_tested_not_responded` field and narrowed
  `groups_not_known` semantics.

## Scope

| Deliverable | Repo | Type |
|-------------|------|------|
| Phase 0 pre-flight + hard gates in research checklist | multiomics_research | Skill change |
| Prior-work check in Phase 0 | multiomics_research | Skill change |
| Process retrospective as mandatory exit artifact | multiomics_research | Skill change |
| Python API guide (`references/python-api-guide.md`) | multiomics_research | New reference |
| `gene_response_profile` spec change | multiomics_explorer | API spec (input/output only) |
| Utility function reference docs | multiomics_explorer | New reference files for `analysis` module functions |

Out of scope:
- Implementation details for `gene_response_profile` change
- Changes to other MCP tools
- Eval framework changes
- Skill-creator-based eval of the skill (deferred to after next iteration)

---

## Design

### 1. Research skill restructure: state machine with hard gates

Replace the current 2-stage workflow with a 3-phase state machine.
Each phase has a hard gate — a set of file-existence and content
checks that must pass before proceeding.

```
Phase 0: Pre-flight
  |
  GATE 0: methods.md exists, gaps_and_friction.md exists, exploration/ exists
  |
Phase 1: Orientation
  |
  GATE 1: exploration log has ## Findings, methods.md has ## Data scope,
          gaps_and_friction.md updated
  |
Phase 2: Research loop (iterate)
  |
  PER-ITERATION EXIT: exploration log written, scripts saved,
                       gaps_and_friction.md updated
  |
  EXIT GATE: README.md written, all scripts in scripts/,
             gaps_and_friction.md has ## Process retrospective
```

#### Phase 0: Pre-flight (new)

Entirely mechanical — no biology, no MCP calls. This phase exists
to make artifact-first discipline structural, not aspirational.

Steps:

1. **Check for prior work.** Search `analyses/` for existing analyses
   on the same topic. If found:
   - Read the prior README.md and gaps_and_friction.md
   - Note in the new methods.md: "Builds on / diverges from
     `analyses/{prior}/`" with rationale
   - Reference confirmed findings rather than re-deriving them
   - If the user explicitly asks to start fresh, note that too
2. **Create analysis directory structure:**
   ```
   analyses/{analysis_name}/
   ├── exploration/
   ├── data/
   ├── scripts/
   ├── results/
   ├── methods.md          (stub: ## Research question only)
   ├── gaps_and_friction.md (stub: headers only)
   ```
3. **Verify gate:** methods.md exists with ## Research question,
   gaps_and_friction.md exists, exploration/ directory exists.

**Gate language in the skill:**
> "Do not make any MCP tool calls or KG queries until Phase 0 is
> complete. The first tool call in any analysis must be preceded by
> file creation."

#### Phase 1: Orientation (unchanged from iterative redesign)

Scope KG, identify genes, resolve locus tags. Write the first
exploration log. Update methods.md with ## Data scope and ## Gene
identification.

**Gate 1 checks:**
- exploration/YYYY-MM-DD-orientation*.md exists and has ## Findings
- methods.md has ## Data scope section
- gaps_and_friction.md updated (even if "none yet")

#### Phase 2: Research loop (minor additions)

The Question -> Explore -> Log -> Assess loop from the iterative
redesign spec is unchanged. Two additions:

**Per-iteration exit checks:**
- Exploration log written with all required sections
- Any scripts created during this iteration saved to scripts/
  (not left as inline chat code)
- gaps_and_friction.md updated (even if "no new issues")

**Exit gate (end of analysis):**
- README.md written with file index and key findings
- All scripts in scripts/ (no orphaned code in chat)
- gaps_and_friction.md has ## Process retrospective section with:
  - What worked
  - What didn't work
  - Proposed changes (to skill, MCP, KG)

#### Gate enforcement language

The skill should use imperative, non-negotiable language for gates:

> **GATE 0 — Pre-flight complete.**
> Before any MCP call: methods.md, gaps_and_friction.md, and
> exploration/ must exist. If they do not, create them now.
> Do not proceed until these files are on disk.

This is stronger than the current checklist style ("- [ ] Create
methods.md stub"). The gates are not items to check off — they are
conditions that block progress.

---

### 2. Python API guide

New file: `skills/research/references/python-api-guide.md`

Referenced from the research checklist's Explore step and from the
artifacts guide's data extraction section.

#### Contents

**2.1 Installation and import**

- Package: `multiomics-explorer`
- Installed via this repo's `pyproject.toml` `[analysis]` extras:
  `pip install -e ".[analysis]"` or `uv sync --extra analysis`
- MCP tool functions: `from multiomics_explorer import <function_name>`
- Analysis utilities: `from multiomics_explorer.analysis import response_matrix, gene_set_compare`
- Run scripts from the `multiomics_research` project root
- Requires Neo4j running (same connection as MCP server)

**2.2 Core principle: MCP tools = Python API functions**

Every MCP tool has an identical Python function: same name, same
parameters, same return dict structure. The tool guides loaded as
MCP resources (`docs://tools/{tool_name}`) document both interfaces.

- Use MCP for interactive exploration and reasoning
- Use Python API for scripting, bulk extraction, and computation
- When a research iteration needs data reshaped or aggregated,
  switch to a Python script

**2.3 Return structure**

Every function returns a dict with:
- Envelope fields (metadata): `total_matching`, `returned`,
  `truncated`, `offset`, etc.
- `results`: list of dicts, one per result row

Standard pattern for tabular data:
```python
result = differential_expression_by_gene(
    organism="MED4", locus_tags=["PMM0370"], limit=None
)
df = pd.DataFrame(result["results"])
```

Use `limit=None` to retrieve all matching rows (no pagination
needed).

**2.4 Handling nested fields**

Some result fields are lists or dicts. `pd.DataFrame()` will store
these as Python objects in cells, which breaks CSV export and
filtering.

Fields that need special handling:

| Function | Field | Type | How to handle |
|----------|-------|------|---------------|
| `gene_response_profile` | `groups_responded` | list | Join to string or explode |
| `gene_response_profile` | `groups_not_responded` | list | Join to string or explode |
| `gene_response_profile` | `groups_not_known` | list | Join to string or explode |
| `gene_response_profile` | `response_summary` | nested dict | Extract separately, don't put in main DataFrame |
| `gene_overview` | `annotation_types` | list | Join to string |
| `list_experiments` | `timepoints` | list of dicts | Extract separately |

Example — extracting response_summary into a flat table:
```python
from multiomics_explorer import gene_response_profile

result = gene_response_profile(locus_tags=["PMM0370", "PMM0920"])

# Flat gene-level table (drop nested fields)
genes_df = pd.DataFrame(result["results"])
list_cols = ["groups_responded", "groups_not_responded", "groups_not_known"]
for col in list_cols:
    if col in genes_df.columns:
        genes_df[col] = genes_df[col].apply(lambda x: "; ".join(x) if isinstance(x, list) else x)
genes_df = genes_df.drop(columns=["response_summary"])

# Per-gene x per-treatment summary (from nested response_summary)
summary_rows = []
for gene in result["results"]:
    for group, stats in gene["response_summary"].items():
        summary_rows.append({"locus_tag": gene["locus_tag"],
                             "gene_name": gene["gene_name"],
                             "group": group, **stats})
summary_df = pd.DataFrame(summary_rows)
```

**2.5 Pre-script checklist**

Before writing any extraction script:
1. Verify import: `from multiomics_explorer import <func>`
2. Test return schema with a minimal call:
   ```python
   result = func(..., limit=1)
   print(result.keys())
   print(result["results"][0].keys())
   ```
3. Check column names before filtering — don't assume field names
4. Check `total_matching` to understand dataset size

**2.6 Common mistakes**

| Mistake | Fix |
|---------|-----|
| `uv run --directory /path/to/multiomics_explorer script.py` | Run from `multiomics_research` root. The package is installed in this project's environment. |
| Guessing column names (`product`, `is_significant`) | Test with `result["results"][0].keys()` first |
| `pd.DataFrame(result["results"])` on nested data then `.to_csv()` | Flatten list/dict columns first (see 2.4) |
| Hardcoding full organism name | Fuzzy matching works: `"MED4"` resolves correctly |
| Writing raw Neo4j/requests calls | Use the API functions — they handle connection, serialization, error handling |
| Setting arbitrary high `limit=10000` | Use `limit=None` to get all results |

---

### 3. Explorer changes (utility docs + gene_response_profile)

Detailed in a standalone spec that can be loaded in the explorer
repo independently:

[`2026-03-31-explorer-utils-docs-and-response-profile-change.md`](2026-03-31-explorer-utils-docs-and-response-profile-change.md)

Summary:

- **Utility reference docs:** Write `references/utils/response_matrix.md`
  and `references/utils/gene_set_compare.md` in the explorer's skill
  directory, register as MCP resources at `docs://utils/{name}`.
  Functions already exist at `multiomics_explorer.analysis`.
- **`gene_response_profile` change:** Add `groups_tested_not_responded`
  field for treatment groups where absence means "measured, not
  significant" (based on table_scope). Narrow `groups_not_known` to
  truly uninformative absence only. Optionally add `table_scopes` to
  per-group response_summary.
- **`response_matrix` update:** Handle `groups_tested_not_responded`
  as `"not_responded"`.

The Python API guide (Section 2) should note that
`groups_tested_not_responded` exists and is the correct field for
specificity analyses. Until implemented, the workaround is to
manually check table_scope via `list_experiments`.

---

## Validation

### Skill changes

The next analysis iteration serves as the test. Success criteria:

- Phase 0 completes before any MCP call (methods.md and
  gaps_and_friction.md exist before first query)
- Prior work is checked and referenced if relevant
- All scripts saved to scripts/ (no orphaned code in chat)
- gaps_and_friction.md updated per iteration
- Process retrospective written at the end
- The skill works without project memory loaded (fresh session)

### Python API guide

- Claude uses the correct import path without trial and error
- Claude tests return schema before writing extraction scripts
- No "column not found" or "import not found" errors
- Nested fields handled correctly in DataFrames

### `gene_response_profile` change

- Validated when implemented in multiomics_explorer
- The N-stress markers analysis specificity table should be
  reproducible without manual table_scope correction

## Implementation order

1. **Python API guide** — new reference file in this repo, no code changes
2. **Skill restructure** — update research-checklist.md with Phase 0
   and hard gates, update SKILL.md to reference the new structure
3. **Utility function reference docs** — new files in
   multiomics_explorer `references/utils/`, register as MCP resources
4. **`gene_response_profile` spec change** — implement in
   multiomics_explorer

Steps 1 and 2 are independent. Step 3 depends on the utility
functions existing (from the prior spec). Step 4 is independent
of 1-3 but lower priority (workaround exists).
