# Research Notebook Discipline — Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Update the research-methodology skill to enforce interactive notebook discipline with do→show→explore→decide checkpoints on every analytical step.

**Architecture:** One new reference file (`references/research-notebook.md`) containing the notebook format, step cycle, QC checkpoint types, and code lifecycle rules. Targeted edits to 5 existing files. No code changes — all documentation/skill content.

**Spec:** `docs/superpowers/specs/2026-04-07-research-notebook-discipline-design.md`

---

## File Map

| File | Action | Responsibility |
|------|--------|---------------|
| `skills/research-methodology/references/research-notebook.md` | Create | Notebook format, step cycle, QC checkpoints, manifest rules, code lifecycle |
| `skills/research-methodology/SKILL.md` | Modify | Add load-timing advisory + Rule 8 |
| `skills/research-methodology/references/artifacts.md` | Modify | Add required docs to directory structure, add templates, replace exploration log with notebook pointer |
| `skills/research-methodology/references/statistical-rigor.md` | Modify | Add worked-example requirement for formulas |
| `skills/research-methodology/references/python-api-guide.md` | Modify | Add verbose join-early rule + single-timepoint check |
| `CLAUDE.md` | Modify | Add skill load-timing instruction + subagent review override |

---

### Task 1: Create `references/research-notebook.md`

**Files:**
- Create: `skills/research-methodology/references/research-notebook.md`

This is the main deliverable — the full notebook discipline reference.

- [ ] **Step 1: Write the research-notebook.md file**

```markdown
# Research notebook

Every analysis is driven by an interactive research notebook — a
chronological log where each analytical step is recorded, inspected,
explored with the researcher, and approved before the next step
proceeds. Implementation can be fast or delegated; quality control
and exploration are always interactive.

## The step cycle: do → show → explore → decide

Every step that produces data or analytical output follows this cycle:

1. **Do** — run the script, produce artifacts (CSV, figures,
   summaries). Can be delegated to a subagent.
2. **Show** — present QC diagnostics: what went in, what came out,
   what's missing, sanity checks.
3. **Explore** — interactive walkthrough in chat: pick specific
   genes/conditions, verify the logic on concrete examples, ask and
   answer questions. Capture the key data points and conclusions in
   the notebook.
4. **Decide** — researcher says "continue", "explain X", or "redo
   with Y". The decision is logged. Only then does the next step
   begin.

**No step that produces data or analytical output may proceed
without completing the full cycle.**

Mechanical tasks (formatting, plotting from existing data, file
reorganization) can skip the explore phase but still need
show + decide.

## What a "step" is

A step is:
- A Python script with explicit command-line invocation
- Output artifacts (CSV, figures, summary files)
- A notebook entry capturing all four phases of the cycle

## Notebook format

One notebook per analysis: `exploration/YYYY-MM-DD-notebook.md`

Chronological and append-only. Reruns add new entries, not
overwrites. This is a messy lab notebook, not a publication.

### Entry template

~~~markdown
---

## YYYY-MM-DD HH:MM — Step N: {description}

### Command
```bash
uv run scripts/script_name.py --arg value --output data/output.csv
```

### Inputs
- [input_file.csv](../data/input_file.csv) — N rows, N genes, N timepoints

### Outputs
- [output_file.csv](../data/output_file.csv) — N rows, description
- [figure.png](../results/figure.png) — description

### QC
- Summary statistics, counts, sanity checks
- What's present, what's missing, what's surprising

### Exploration
- Walked through gene X (PMM0xxx): values, logic, conclusion
- Checked gene Y: values — expected / unexpected because [reason]
- Asked: [question] → [answer with data]

### Decision
What was decided and why. Proceed / redo / adjust.
~~~

### What the notebook captures

Everything:
- **Spec walkthrough** — the initial section-by-section review of
  the analysis design (questions, decisions, rationale)
- **Each analytical step** — command, inputs, outputs, QC,
  exploration, decision
- **Reruns** — new entry with "Why" section explaining what changed
  and the new output
- **Chat-based follow-up** — interpretation discussions, "what
  about gene X?" explorations, with actual data points
- **Questions and decisions** — including dead ends and rejected
  alternatives
- **Relative links** to all artifacts produced

### Rerun entries

When a step is rerun after a change, add a new entry (don't
overwrite the original):

~~~markdown
---

## YYYY-MM-DD HH:MM — Step N rerun: {what changed}

### Why
What the researcher flagged and why the step needed to change.

### Changes
- What was modified (file, line, logic change)
- How many genes/rows/results were affected

### Command
```bash
uv run scripts/script_name.py --arg new_value --output data/output.csv
```

### Outputs
- [output_file.csv](../data/output_file.csv) — N rows (updated)

### Decision
Proceed with updated output / need another iteration.
~~~

## QC checkpoint types

What to show depends on the step type.

### Data extraction
- Row count, gene count, timepoint/condition count
- Sample rows (first 5 or curated edge cases)
- What's missing: expected genes absent, unexpected metadata
  (e.g., `timepoint=single`)
- One-sentence summary: "This is X genes across Y conditions
  from Z experiments"

### Gene selection / filtering
- Counts at each filter step: started with → filter 1 →
  filter 2 → final
- Sample excluded genes with reason
- Sample included genes — known markers present?
- Flag surprises: "Gene X was excluded because Y — expected?"

### Computation / metric
- Worked example: 2-3 genes through the formula with actual
  numbers, step by step
- Summary statistics of the output (distribution, range, NaNs)
- Sanity check against known biology: "glnA should score
  high — does it?"

### Scoring / comparison
- Full results table in markdown (not prose summary — the
  actual numbers)
- Best/worst scores, surprises, anything unexpected
- Cross-condition comparison with expectation check

## Manifest updates

Every step that produces output files updates the relevant
manifest (`data/DATA_MANIFEST.md` or `results/RESULTS_MANIFEST.md`)
immediately — not retroactively at the end.

## Code lifecycle: analysis-first, productize later

Research code has two phases. Follow phase 1 during analysis;
flag phase 2 candidates in the notebook.

### Phase 1 — Analysis code (methodology-first)

Code is written to solve the immediate research problem. It lives
in the analysis directory. The goal is correct methodology, not
good software engineering.

**Separate reusable logic from scripts.** When an analysis
introduces a new methodology (signature scoring, enrichment,
concordance metrics), put the reusable logic in a utility package
within the analysis directory (`analyses/*/sig_utils/`,
`analyses/*/enrich_utils/`). Scripts call the utilities for
specific data; utilities contain the methodology. This separation
is what makes toy-testing possible and productization
straightforward later.

**Specs describe methodology, not implementation.** For novel
utilities (scoring functions, metrics, gene set operations), the
spec should contain formulas with worked examples, expected I/O,
and pseudocode — not Python implementation. Regular
extraction/plotting scripts are straightforward enough to go
directly in the plan.

**Toy-data verification before real data.** When building a
reusable utility, verify with hand-calculated toy examples first.
This is a notebook step: create small synthetic input, compute
expected output by hand, run the utility, compare, log the
verification. Applies to anything in a shared `*_utils/` package.
One-off scripts don't need it.

**Refine through the notebook QC cycle.** The interactive
do→show→explore→decide loop is how the methodology gets validated.
Formula corrections, edge cases, direction logic — all discovered
through the researcher walking through concrete examples.

### Phase 2 — Productization (software-first)

After the analysis, if a utility proves reusable (used across
multiple analyses or conditions), flag it for productization — a
separate brainstorm with proper API design, tests, and
documentation. It moves from the analysis directory to a shared
package (e.g., `multiomics_explorer/analysis/`).

Don't productize speculatively. Wait for proven reuse — the
analysis notebook is the evidence.
```

- [ ] **Step 2: Verify the file reads correctly**

Run: read the file, confirm no formatting issues or broken markdown.

- [ ] **Step 3: Commit**

```bash
git add skills/research-methodology/references/research-notebook.md
git commit -m "feat: add research notebook discipline reference"
```

---

### Task 2: Update SKILL.md

**Files:**
- Modify: `skills/research-methodology/SKILL.md`

Add load-timing advisory at the top (before Rule 1) and Rule 8 at the bottom.

- [ ] **Step 1: Add load-timing advisory**

After the frontmatter and before Rule 1, add:

```markdown
> **Load this skill BEFORE brainstorming.** The notebook structure,
> artifact requirements, and methodology rules shape the entire
> analysis design. Loading after the spec is written means
> retrofitting — which is what happened in the N-limitation analysis.
```

- [ ] **Step 2: Add Rule 8**

After Rule 7 (Don't hallucinate), add:

```markdown
## Rule 8: Research notebook, not pipeline

Every analysis is driven by an interactive notebook where each
step is recorded, explored with the researcher, and approved
before the next step. Implementation can be fast or delegated;
quality control is always interactive.

See [Research notebook](references/research-notebook.md) for the
notebook format, QC checkpoint requirements, step cycle, and code
lifecycle rules.
```

- [ ] **Step 3: Add research-notebook.md to the References list**

Add to the References section at the bottom:

```markdown
- [Research notebook](references/research-notebook.md) — notebook format, step cycle, QC checkpoints, code lifecycle
```

- [ ] **Step 4: Commit**

```bash
git add skills/research-methodology/SKILL.md
git commit -m "feat: add Rule 8 (research notebook) and load-timing advisory to SKILL.md"
```

---

### Task 3: Update `references/artifacts.md`

**Files:**
- Modify: `skills/research-methodology/references/artifacts.md`

Three changes: update directory structure, add document templates, replace exploration log section.

- [ ] **Step 1: Update directory structure**

Replace the directory structure block at lines 23-33 of `artifacts.md` (the fenced code block under `## Directory structure`) with:

```
analyses/{analysis_name}/
├── exploration/       # Research notebook (one per analysis)
├── data/              # Staged data from KG (CSV/TSV)
│   └── DATA_MANIFEST.md  # What each data file contains
├── scripts/           # Python scripts (extract, analyze, explore)
├── results/           # Outputs: tables, figures, statistics
│   └── RESULTS_MANIFEST.md  # What each result file contains
├── README.md          # Summary, key findings, file index
├── methods.md         # Publication-ready methods (living document)
├── decisions.md       # Design decisions with rationale (why, not what)
├── caveats.md         # Interpretation caveats for readers of results
├── gaps_and_friction.md  # KG/MCP/methodology issues log
└── references.md      # Data sources and citations
```

- [ ] **Step 2: Replace exploration log section**

Replace lines 39-77 of `artifacts.md` — the entire `## Exploration logs` section including the `~~~markdown` template block and the "Chain between iterations" paragraph. Replace with:

```markdown
## Research notebook

One notebook per analysis: `exploration/YYYY-MM-DD-notebook.md`.
A chronological, append-only lab notebook logging every analytical
step as it happens. Not a retrospective summary.

See [Research notebook](research-notebook.md) for the full format,
step cycle, QC checkpoint types, and code lifecycle rules.

The README should link to the notebook.
```

- [ ] **Step 3: Add document templates after the `## results/` section**

Insert after line 157 of `artifacts.md` (end of `## results/` section, before `## References and citations`):

```markdown
## Required analysis documents

### `data/DATA_MANIFEST.md`

Updated with every extraction step — not retroactively.

~~~markdown
# Data Manifest

All files produced by extraction scripts from the multiomics KG.
Run scripts from `multiomics_research` root with `uv run`.

## {Section name}

| File | Rows | Genes | Timepoints | Produced by | Description |
|------|------|-------|------------|-------------|-------------|
| `filename.csv` | N | N | list | `script.py` | One-line description. Note any filtering or scope limitations. |

## Common columns (all DE files)

List shared columns with brief descriptions. Note any non-obvious
semantics (e.g., `not_significant` vs `not_known`, directional
rank nulls).
~~~

### `results/RESULTS_MANIFEST.md`

Updated with every scoring, plotting, or analysis step.

~~~markdown
# Results Manifest

All files produced by scoring, plotting, and analysis scripts.

## {Section name}

| File | Rows | Produced by | Description |
|------|------|-------------|-------------|
| `filename.csv` | N | `script.py` | One-line description. List key columns. |

## Figures

| File | Produced by | Description |
|------|-------------|-------------|
| `figure.png` | `script.py` | What it shows, what axes/panels represent. |
~~~

### `decisions.md`

Design decisions with rationale. Written during brainstorming and
updated during analysis when new decisions are made.

~~~markdown
# Decision Log

Design decisions with rationale — WHY the analysis was done this
way, not what it does (that's in methods.md).

## {Topic area}

### {Decision title}

**Decision:** What was chosen.

**Rationale:** Why. What alternatives were considered and rejected.

**Status:** Implemented / Open question / Superseded by [X].
~~~

### `caveats.md`

Interpretation caveats for readers of results. Distinct from
gaps_and_friction.md (which is about process/tooling). This is
the "fine print" for any figure or claim.

~~~markdown
# Caveats for Interpretation

Things a reader of these results needs to know before drawing
conclusions.

## {Caveat title}

- Bullet points explaining the limitation
- How it affects interpretation
- What it means for specific claims or comparisons
~~~
```

- [ ] **Step 4: Commit**

```bash
git add skills/research-methodology/references/artifacts.md
git commit -m "feat: add required docs, templates, and notebook pointer to artifacts guide"
```

---

### Task 4: Update `references/statistical-rigor.md`

**Files:**
- Modify: `skills/research-methodology/references/statistical-rigor.md`

- [ ] **Step 1: Add worked-example requirement**

After the `## What you must compute in scripts` section, add:

```markdown
## Worked examples in specs

Every formula in a spec must include at least one worked example
with concrete numbers showing input → computation → output.

Shorthand like "reference direction as expected sign" is not
sufficient — a researcher must be able to verify the formula
produces the expected result by reading the example. If there are
edge cases (e.g., what happens when a gene is not significant?),
include an example for those too.
```

- [ ] **Step 2: Commit**

```bash
git add skills/research-methodology/references/statistical-rigor.md
git commit -m "feat: require worked examples for every formula in specs"
```

---

### Task 5: Update `references/python-api-guide.md`

**Files:**
- Modify: `skills/research-methodology/references/python-api-guide.md`

- [ ] **Step 1: Add verbose join-early rule and single-timepoint check**

Add to the `## 6. Common Mistakes` table:

```markdown
| Extracting with `verbose=True` but not joining `product`/`gene_category` into gene-level outputs | Join metadata immediately — don't defer annotation to a separate step. Unnamed genes in signature CSVs are opaque. |
| Not checking for `timepoint=single` or null after extraction | After extraction, explicitly check for experiments with single/null timepoints. Confirm they're handled correctly in downstream scripts (plotting, groupby). |
```

- [ ] **Step 2: Commit**

```bash
git add skills/research-methodology/references/python-api-guide.md
git commit -m "feat: add verbose-join-early and single-timepoint-check rules"
```

---

### Task 6: Update CLAUDE.md

**Files:**
- Modify: `CLAUDE.md`

- [ ] **Step 1: Update the research methodology section**

Replace the existing `## Research methodology` section with:

```markdown
## Research methodology

**Load the `research-methodology` skill BEFORE brainstorming.** It
contains the rules for KG usage, gene identity, artifact structure,
notebook discipline, and anti-hallucination. These rules shape the
analysis design — loading after the spec means retrofitting.

For specific analysis types, invoke the corresponding recipe skill
(e.g., `enrichment`, `response-matrix`, `conservation`).

Use the superpowers workflow for research: brainstorm the question,
write a plan, execute with checkpoints. The methodology skill
provides the domain rules; superpowers provides the process
discipline.

### Process overrides

- **Don't skip subagent reviews** for tasks that produce data
  outputs. At minimum, run spec compliance review. The tasks that
  seem simple are where silent bugs hide.
```

- [ ] **Step 2: Commit**

```bash
git add CLAUDE.md
git commit -m "feat: add skill load-timing and subagent review override to CLAUDE.md"
```

---

### Task 7: Final verification

- [ ] **Step 1: Verify all cross-references**

Check that these links resolve:
- `SKILL.md` → `references/research-notebook.md`
- `artifacts.md` → `research-notebook.md` (relative within references/)
- `SKILL.md` References section lists research-notebook.md

- [ ] **Step 2: Read through the complete updated SKILL.md**

Verify rules 1-8 flow logically, no duplication, load-timing advisory is visible at the top.

- [ ] **Step 3: Commit any fixes**

```bash
git add -A skills/research-methodology/ CLAUDE.md
git commit -m "fix: resolve cross-reference issues"
```

(Skip this step if no fixes needed.)
