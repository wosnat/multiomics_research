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
