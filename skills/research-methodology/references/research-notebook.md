# Research notebook

## Step protocol and enforcement

For per-step commit timing, hard gates, and the chat-capture
pattern, see [Step protocol](step-protocol.md). This document
owns the notebook **format and content**; step-protocol owns
**when things happen and what gates enforce them**.

---

Every analysis is driven by an interactive research notebook — a
chronological log where each analytical step is recorded, inspected,
explored with the researcher, and approved before the next step
proceeds. Implementation can be fast or delegated; quality control
and exploration are always interactive.

## The step cycle: do → show → explore → decide

Every step that produces data or analytical output follows this
cycle:

1. **Do** — write/update script, run, update manifests, commit.
2. **Show** — present QC diagnostics, begin notebook entry with
   summary tables, figure links, counts.
3. **Explore** — interactive walkthrough, append chat-capture
   section to notebook entry.
4. **Decide** — researcher says continue/redo/adjust. Decision
   logged, notebook entry committed.

Per-phase obligations, commit timing, and hard gates are in
[Step protocol](step-protocol.md). What follows here is the
notebook **format**.

Mechanical tasks (formatting, plotting from existing data, file
reorganization) can skip the explore phase but still need
show + decide.

## What a "step" is

A step is:
- A Python script with explicit command-line invocation
- Output artifacts (CSV, figures, summary files)
- A notebook entry capturing all four phases of the cycle

## Interactive discovery steps

Some steps are naturally exploratory — browsing experiments,
classifying conditions, discussing what to include. These don't
fit neatly into "run a script, check the output." The pattern for
interactive steps:

1. Explore interactively via MCP queries and chat discussion
2. Produce a **frozen output file** (CSV, table) that downstream
   scripts consume — this is the reproducible artifact
3. Write a **notebook entry** documenting the reasoning,
   classifications, and decisions — detailed enough that a
   different researcher could follow the logic
4. Optionally, write a script that reproduces the frozen output
   from KG queries

The frozen output + notebook entry are the minimum. A script is
preferred but not required for discovery steps where the process
is inherently iterative and conversational.

Rule 5 (scripts over chat reasoning) still applies to
computations. Interactive steps are for *discovery and
classification*, not for computing statistics in chat.

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

### Chat exploration
**Q: [researcher's question, as asked]**
Data: [what was looked up / computed to answer it]
Finding: [what the data showed, with concrete numbers]
Impact: [how this affects interpretation or next steps]

**Q: [next question]**
...

### Decision
What was decided and why. Proceed / redo / adjust.
~~~

Logs must be sufficient to verify the step without rerunning —
see [Artifacts guide — Log verbosity](artifacts.md) for the
standard.

### What the notebook captures

Everything:
- **Spec walkthrough** — the initial section-by-section review of
  the analysis design (questions, decisions, rationale)
- **Each analytical step** — command, inputs, outputs, QC,
  exploration, decision
- **Reruns** — see [Rerun and revision workflow](#rerun-and-revision-workflow)
- **Chat-based follow-up** — interpretation discussions, "what
  about gene X?" explorations, with actual data points
- **Questions and decisions** — including dead ends and rejected
  alternatives
- **Relative links** to all artifacts produced

## Rerun and revision workflow

When a step is rerun (due to a formula correction, normalization
change, filter adjustment, or any methodology change mid-analysis):

1. **Update the utility code** — change the shared methodology
   (`*_utils/`), not just the script
2. **Run tests** — verify the change is correct against toy data
3. **Rerun all downstream steps** that consumed the changed
   output — every step after the change gets a fresh run
4. **New notebook entry per rerun** — each rerun gets its own
   entry explaining what changed and what the impact was
5. **Commit each rerun's outputs** — see
   [Step protocol](step-protocol.md) for commit timing during
   redos

The rerun entry template:

~~~markdown
---

## YYYY-MM-DD HH:MM — Step N rerun: {what changed}

### Why
What was wrong with the previous output and how it was discovered.

### Changes
- What code was modified (file, function, logic change)
- How tests were affected (new tests, changed expected values)
- How many downstream steps needed rerunning

### Outputs
- Updated files with new row counts / statistics

### Impact
- What changed in the results (direction, magnitude, interpretation)
- What didn't change (confirming robustness where applicable)

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
