# Skill Hardening and Python API Guide Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Harden the research skill with Phase 0 pre-flight hard gates and add a Python API guide so Claude follows the methodology and can script with the API without trial-and-error.

**Architecture:** Three files changed (research-checklist.md, SKILL.md, artifacts-guide.md) and one new file (python-api-guide.md). The checklist gets a new Phase 0 section with hard gate language. SKILL.md gets updated to reference Phase 0 and the API guide. The artifacts guide gets a pointer to the API guide for its data extraction section.

**Tech Stack:** Markdown skill files (no code changes)

**Spec:** [`docs/superpowers/specs/2026-03-31-skill-hardening-and-api-guide.md`](../specs/2026-03-31-skill-hardening-and-api-guide.md)

---

## File Structure

| File | Action | Responsibility |
|------|--------|----------------|
| `skills/research/references/python-api-guide.md` | Create | Import paths, return structure, nested field handling, pre-script checklist, common mistakes |
| `skills/research/references/research-checklist.md` | Modify | Add Phase 0 pre-flight section with hard gates, add per-iteration exit checks, add exit gate with process retrospective |
| `skills/research/SKILL.md` | Modify | Add Phase 0 description to the workflow section, reference python-api-guide.md, update Stage 1/2 gate language |
| `skills/research/references/artifacts-guide.md` | Modify | Add pointer to python-api-guide.md in data extraction section, add process retrospective to gaps_and_friction.md format |

---

### Task 1: Create Python API guide

**Files:**
- Create: `skills/research/references/python-api-guide.md`

- [ ] **Step 1: Write the API guide**

Create `skills/research/references/python-api-guide.md` with the following content:

```markdown
# Python API guide

How to use the `multiomics_explorer` Python package in analysis
scripts. For MCP tool reference (parameters, return format,
examples), read the MCP resource `docs://tools/{tool_name}`.
For analysis utility reference, read `docs://utils/{name}`.

---

## Installation and import

The `multiomics-explorer` package is installed as an optional
dependency of this project:

```bash
# Install (if not already)
uv sync --extra analysis
# or: pip install -e ".[analysis]"
```

Import MCP tool functions directly:
```python
from multiomics_explorer import differential_expression_by_gene
from multiomics_explorer import gene_response_profile
from multiomics_explorer import list_experiments
# ... any of the 19 MCP tool functions
```

Import analysis utilities from the `analysis` module:
```python
from multiomics_explorer.analysis import response_matrix
from multiomics_explorer.analysis import gene_set_compare
```

Run scripts from the `multiomics_research` project root.
Requires Neo4j running (same connection as the MCP server).

---

## Core principle: MCP tools = Python API functions

Every MCP tool has an identical Python function: same name, same
parameters, same return dict structure. The tool guides loaded as
MCP resources (`docs://tools/{tool_name}`) document both.

| Use MCP for | Use Python API for |
|---|---|
| Interactive exploration | Bulk data extraction |
| Browsing and quick lookups | Computation and reshaping |
| Reasoning about results | Producing CSV/tables/figures |

When a research iteration needs data reshaped or aggregated,
switch to a Python script. Do not parse MCP output as text.

---

## Return structure

Every function returns a dict with:
- **Envelope fields** (metadata): `total_matching`, `returned`,
  `truncated`, `offset`, etc.
- **`results`**: list of dicts, one per result row

### Standard extraction pattern

```python
import pandas as pd
from multiomics_explorer import differential_expression_by_gene

result = differential_expression_by_gene(
    organism="MED4",
    locus_tags=["PMM0370", "PMM0920"],
    limit=None,  # get ALL rows, no pagination needed
)
df = pd.DataFrame(result["results"])
df.to_csv("data/output.csv", index=False)
```

Use `limit=None` to retrieve all matching rows.

---

## Handling nested fields

Some result fields are lists or dicts. `pd.DataFrame()` stores
these as Python objects in cells, which breaks CSV export and
filtering. Flatten them before saving.

### Fields that need special handling

| Function | Field | Type | How to handle |
|----------|-------|------|---------------|
| `gene_response_profile` | `groups_responded` | list | Join: `"; ".join(x)` |
| `gene_response_profile` | `groups_not_responded` | list | Join: `"; ".join(x)` |
| `gene_response_profile` | `groups_not_known` | list | Join: `"; ".join(x)` |
| `gene_response_profile` | `response_summary` | nested dict | Extract separately (see below) |
| `gene_overview` | `annotation_types` | list | Join: `"; ".join(x)` |
| `list_experiments` | `timepoints` | list of dicts | Extract separately |

### Example: extracting gene_response_profile into flat tables

```python
import pandas as pd
from multiomics_explorer import gene_response_profile

result = gene_response_profile(locus_tags=["PMM0370", "PMM0920"])

# 1. Flat gene-level table (drop nested fields)
genes_df = pd.DataFrame(result["results"])
list_cols = ["groups_responded", "groups_not_responded", "groups_not_known"]
for col in list_cols:
    if col in genes_df.columns:
        genes_df[col] = genes_df[col].apply(
            lambda x: "; ".join(x) if isinstance(x, list) else x
        )
genes_df = genes_df.drop(columns=["response_summary"])

# 2. Per-gene x per-treatment summary (from nested response_summary)
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
```

---

## Pre-script checklist

Before writing any extraction script:

1. **Verify import:**
   ```python
   from multiomics_explorer import <function_name>
   ```

2. **Test return schema** with a minimal call:
   ```python
   result = func(..., limit=1)
   print(result.keys())
   print(result["results"][0].keys())
   ```

3. **Check column names** before filtering — don't assume field
   names like `product` or `is_significant` exist.

4. **Check `total_matching`** to understand dataset size before
   setting `limit=None`.

---

## Common mistakes

| Mistake | Fix |
|---------|-----|
| `uv run --directory /path/to/multiomics_explorer script.py` | Run from `multiomics_research` root — the package is installed in this project's environment |
| Guessing column names (`product`, `is_significant`) | Test with `result["results"][0].keys()` first |
| `pd.DataFrame(result["results"])` then `.to_csv()` on nested data | Flatten list/dict columns first (see above) |
| Hardcoding full organism name | Fuzzy matching works: `"MED4"` resolves correctly |
| Writing raw Neo4j or requests calls | Use the API functions — they handle connection and serialization |
| Setting arbitrary high `limit=10000` | Use `limit=None` to get all results |
| Using `response_matrix()` output without checking `not_known` cells | `not_known` may mean "measured but not significant" depending on table_scope — check `list_experiments` for the treatment type |
```

- [ ] **Step 2: Verify the file reads correctly**

Read the file back and check for any markdown formatting issues.

- [ ] **Step 3: Commit**

```bash
git add skills/research/references/python-api-guide.md
git commit -m "Add Python API guide for multiomics_explorer scripting"
```

---

### Task 2: Add Phase 0 pre-flight to research checklist

**Files:**
- Modify: `skills/research/references/research-checklist.md`

The current checklist has: "Before starting", "Orientation", "Research loop", "Synthesis", "Review". We need to insert Phase 0 between "Before starting" and "Orientation", add per-iteration exit checks to the research loop, and add the exit gate with process retrospective.

- [ ] **Step 1: Add Phase 0 section after "Before starting"**

Insert the following after the "Before starting" section (after line 20, before the `## Orientation` section):

```markdown
---

## Phase 0: Pre-flight (REQUIRED before any MCP call)

> **GATE 0 — Pre-flight complete.**
> Before any MCP call: the analysis directory, methods.md,
> gaps_and_friction.md, and exploration/ must exist on disk.
> Do not proceed until these files are created.

### Prior work check

- [ ] Search `analyses/` for existing analyses on the same topic
- [ ] If prior work exists:
  - Read the prior README.md and gaps_and_friction.md
  - Note in methods.md: "Builds on / diverges from
    `analyses/{prior}/`" with rationale
  - Reference confirmed findings rather than re-deriving them
- [ ] If user explicitly asks to start fresh, note that in
  methods.md

### Create analysis directory

- [ ] Create the directory structure:
  ```
  analyses/{analysis_name}/
  ├── exploration/
  ├── data/
  ├── scripts/
  ├── results/
  ├── methods.md          (stub: ## Research question only)
  ├── gaps_and_friction.md (stub: category headers only)
  ```
- [ ] Verify: methods.md exists with ## Research question
- [ ] Verify: gaps_and_friction.md exists with category headers
- [ ] Verify: exploration/ directory exists
```

- [ ] **Step 2: Add per-iteration exit checks to the Research loop section**

After the existing "### Assess" checklist items in the Research loop section (after the line `- [ ] Decide: next question, or conclude?`), add:

```markdown

### Iteration exit checks

- [ ] Exploration log written with all required sections
- [ ] Any scripts created this iteration saved to `scripts/`
  (not left as inline chat code)
- [ ] `gaps_and_friction.md` updated (even if "no new issues
  this iteration")
```

- [ ] **Step 3: Add exit gate to the Synthesis section**

After the existing Synthesis checklist items, add:

```markdown

### Exit gate

- [ ] All scripts in `scripts/` (no orphaned code in chat)
- [ ] `gaps_and_friction.md` has `## Process retrospective`
  section with:
  - What worked
  - What didn't work
  - Proposed changes (to skill, MCP, KG)
```

- [ ] **Step 4: Add gate language to the Orientation section**

Before the existing Orientation checklist items (before `### Scope`), add:

```markdown

> **GATE 1 — Orientation complete.**
> Before entering the research loop: an exploration log with
> ## Findings must exist, methods.md must have ## Data scope,
> and gaps_and_friction.md must be updated.
```

- [ ] **Step 5: Verify the checklist reads as a coherent document**

Read the full file back and check that the flow is:
Before starting → Phase 0 (pre-flight) → Orientation → Research loop → Synthesis → Review

- [ ] **Step 6: Commit**

```bash
git add skills/research/references/research-checklist.md
git commit -m "Add Phase 0 pre-flight hard gates and iteration exit checks to research checklist"
```

---

### Task 3: Update SKILL.md with Phase 0 and API guide reference

**Files:**
- Modify: `skills/research/SKILL.md`

- [ ] **Step 1: Rename "The 2-stage workflow" to "The 3-phase workflow" and add Phase 0**

Replace the line `## The 2-stage workflow` (line 74) with:

```markdown
## The 3-phase workflow
```

Insert the following Phase 0 section before `### Stage 1: Orientation` (before line 76):

```markdown
### Phase 0: Pre-flight (before any MCP call)

No biology. No MCP calls. Create the analysis directory and
artifact stubs first. Check for prior analyses on the same topic.

See [research checklist — Phase 0](references/research-checklist.md)
for the full checklist.

> **GATE 0:** Do not make any MCP tool calls or KG queries until
> methods.md, gaps_and_friction.md, and exploration/ exist on disk.

```

- [ ] **Step 2: Add gate language to Stage 1**

After the existing Gate section in Stage 1 (after line 106, "If key data is missing, report the gap before proceeding."), add:

```markdown

> **GATE 1:** Do not enter the research loop until an exploration
> log with ## Findings exists, methods.md has ## Data scope, and
> gaps_and_friction.md is updated.
```

- [ ] **Step 3: Add reference to python-api-guide.md in the header**

Update the opening reference block (lines 9-11) from:

```markdown
See [research checklist](references/research-checklist.md) for the
step-by-step protocol, [artifacts guide](references/artifacts-guide.md)
for output structure, and [anti-hallucination](references/anti-hallucination.md)
for concrete failure modes to avoid.
```

to:

```markdown
See [research checklist](references/research-checklist.md) for the
step-by-step protocol, [artifacts guide](references/artifacts-guide.md)
for output structure, [Python API guide](references/python-api-guide.md)
for scripting with the multiomics_explorer package, and
[anti-hallucination](references/anti-hallucination.md) for concrete
failure modes to avoid.
```

- [ ] **Step 4: Add per-iteration and exit gate language to Stage 2**

After the Assess section's last bullet ("Decide: next question, or conclude?", line 172), add:

```markdown

> **Per-iteration exit:** Before moving to the next iteration,
> verify: exploration log written, scripts saved to `scripts/`,
> gaps_and_friction.md updated.

```

After the Synthesis section (after line 180, "Populate `references.md`"), add:

```markdown

> **Exit gate:** Before concluding the analysis: README.md written,
> all scripts in `scripts/`, gaps_and_friction.md has
> `## Process retrospective` (what worked, what didn't, proposed
> changes).

```

- [ ] **Step 5: Add API guide reference to the Explore section**

In the Stage 2 Explore section (around line 126), after the bullet about writing scripts, add:

```markdown
- See [Python API guide](references/python-api-guide.md) for
  import paths, return structure, and common scripting patterns.
  **Before writing any script:** test the return schema with a
  minimal call (see guide §Pre-script checklist).
```

- [ ] **Step 6: Verify the file reads coherently**

Read the full SKILL.md back and check the flow:
Phase 0 → Stage 1 (Orientation) → Stage 2 (Research loop) → Synthesis

- [ ] **Step 7: Commit**

```bash
git add skills/research/SKILL.md
git commit -m "Add Phase 0 pre-flight, gate language, and API guide reference to research skill"
```

---

### Task 4: Update artifacts guide with process retrospective and API guide pointer

**Files:**
- Modify: `skills/research/references/artifacts-guide.md`

- [ ] **Step 1: Add process retrospective to gaps_and_friction.md format**

In the `## gaps_and_friction.md` section (around line 109), after the format example block that ends with `## Skill/methodology friction`, add a new section to the format:

```markdown

```markdown
## Process retrospective

### What worked
1. ...

### What didn't work
1. ...

### Proposed changes

**To the skill:**
- ...

**To the MCP/KG:**
- ...
```

Written at the end of the analysis (exit gate). Captures lessons
for skill and tool improvement.
```

- [ ] **Step 2: Add API guide pointer to data extraction section**

In the `## data/` section (around line 140), after the "**How to extract:**" code example, add:

```markdown

See [Python API guide](python-api-guide.md) for import paths,
return structure handling (especially nested fields), and common
mistakes to avoid.
```

- [ ] **Step 3: Verify the file reads coherently**

Read the full artifacts-guide.md back.

- [ ] **Step 4: Commit**

```bash
git add skills/research/references/artifacts-guide.md
git commit -m "Add process retrospective format and API guide pointer to artifacts guide"
```

---

### Task 5: Final review

- [ ] **Step 1: Read all four files end-to-end**

Read these files and verify cross-references are correct:
- `skills/research/SKILL.md` — references python-api-guide.md, research-checklist.md
- `skills/research/references/research-checklist.md` — Phase 0/1/2 flow
- `skills/research/references/python-api-guide.md` — standalone, references MCP resources
- `skills/research/references/artifacts-guide.md` — references python-api-guide.md

- [ ] **Step 2: Check that gate language is consistent**

Verify that GATE 0, GATE 1, per-iteration exit, and exit gate
language matches between SKILL.md and research-checklist.md.

- [ ] **Step 3: Fix any inconsistencies found**

If cross-references are broken or gate language diverges, fix
inline.

- [ ] **Step 4: Final commit if any fixes were made**

```bash
git add -u skills/research/
git commit -m "Fix cross-references in research skill files"
```
