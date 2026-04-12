# Artifacts guide

Research questions produce files, not chat messages. Chat is for
reasoning, planning, and interpretation. Data, statistics, figures,
and exploration logs go to disk.

## Contents

1. [When to create artifacts](#when-to-create-artifacts)
2. [Directory structure](#directory-structure)
3. [Superpowers artifact capture](#superpowers-artifact-capture)
4. [Research notebook](#research-notebook)
5. [gaps_and_friction.md](#gaps_and_frictionmd)
6. [methods.md](#methodsmd)
7. [Script naming conventions](#script-naming-conventions)
8. [data/](#data)
9. [results/](#results)
10. [Git tracking convention](#git-tracking-convention)
11. [Log verbosity](#log-verbosity)
12. [Required analysis documents](#required-analysis-documents) — DATA_MANIFEST, RESULTS_MANIFEST, decisions, caveats
13. [References and citations](#references-and-citations)

## When to create artifacts

Create an analysis directory when:
- The question requires statistical computation
- Multiple data extractions feed into an analysis
- Results need to be shared or revisited
- The analysis takes more than 2-3 tool calls
- Figures or publication-ready tables are needed

Simple lookups and quick answers don't need artifacts:
- "What is gene PMM0120?" — MCP lookup, answer in chat
- "How many genes respond to nitrogen?" — MCP summary, answer in chat
- "List the catalase genes" — MCP query, table in chat

## Directory structure

```
analyses/YYYY-MM-DD-HHMM-{name}/
├── exploration/           # Research notebook
├── data/                  # Staged KG data (CSV/TSV)
│   └── DATA_MANIFEST.md
├── scripts/               # Pipeline scripts (numbered)
├── {name}_utils/          # Reusable methodology package (when applicable)
│   ├── __init__.py
│   └── tests/
├── logs/                  # Per-step diagnostic logs
├── results/               # Tables, figures, statistics
│   └── RESULTS_MANIFEST.md
├── superpowers/           # Spec, plan, brainstorm-log (copied here)
│   ├── spec.md
│   ├── plan.md
│   └── brainstorm-log.md
├── README.md
├── methods.md
├── decisions.md
├── caveats.md
├── gaps_and_friction.md
└── references.md
```

**Naming:** `YYYY-MM-DD-HHMM-{name}` — timestamp prefix for
chronological ordering. `{name}` is descriptive and lowercase
with underscores (e.g., `2026-04-08-1038-n_limitation_signature_v2`).

**New directories:**
- **`{name}_utils/`** — reusable methodology package when the
  analysis introduces novel methods (signature scoring, enrichment,
  etc.). Contains pure methodology code with tests. Only created
  when the analysis needs it — not every analysis gets one. See
  [Research notebook — Code lifecycle](research-notebook.md) for
  when this applies.
- **`logs/`** — one log per pipeline step, capturing diagnostics
  sufficient to verify the step without rerunning. See "Log
  verbosity" below.
- **`superpowers/`** — local copies of the spec, plan, and
  brainstorm-log. See "Superpowers artifact capture" below.

## Superpowers artifact capture

After brainstorming and planning, copy the spec, plan, and
brainstorm-log into the analysis's `superpowers/` directory. The
analysis should be self-contained — a reader should find the full
decision history (why this design, what alternatives were
considered, what the plan was) alongside the data and code, not
scattered across `docs/superpowers/`.

The canonical copies in `docs/superpowers/specs/` and
`docs/superpowers/plans/` remain the repo-level index. The copies
in the analysis are the local record.

Files to capture:
- `superpowers/spec.md` — the design spec
- `superpowers/plan.md` — the implementation plan
- `superpowers/brainstorm-log.md` — the Q&A from brainstorming
  (questions asked, options considered, decisions made, rationale)

## Research notebook

One notebook per analysis: `exploration/YYYY-MM-DD-notebook.md`.
A chronological, append-only lab notebook logging every analytical
step as it happens. Not a retrospective summary.

See [Research notebook](research-notebook.md) for the full format,
step cycle, QC checkpoint types, and code lifecycle rules.

The README should link to the notebook.

## gaps_and_friction.md

A running log of KG, MCP, and methodology issues. Updated during
every research iteration, not written retroactively.

Each entry is recorded in two places:
- The exploration log entry (for context)
- `gaps_and_friction.md` (aggregated backlog)

**Format:**

```markdown
# Gaps and friction points

## KG data bugs
1. **{Short description}** — {Details. Action: ...}

## KG gaps
1. **{Short description}** — {What's missing and why it matters.}

## MCP friction
1. **{Short description}** — {What was hard and what would help.}

## Skill/methodology friction
1. **{Short description}** — {What the skill didn't guide well.}

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

## methods.md

A living document updated incrementally as findings become
established. Publication-ready.

**Required sections:**
- **Research question** — operationalized, not the user's original phrasing
- **Data scope** — organisms, publications (DOIs), experiments (IDs, conditions, timepoints), inclusions/exclusions, KG version
- **Gene selection** — how identified, filters, counts at each step, paralog handling
- **Statistical methods** — tests, correction, background set, software versions, thresholds
- **Results summary** — key findings with effect sizes and p-values, references to output files
- **Limitations** — KG gaps, statistical caveats, what the analysis can and cannot conclude

## Script naming conventions

- `scripts/extract_*.py` — data extraction from KG (reusable)
- `scripts/analyze_*.py` — computation that produces results (reusable)
- `scripts/explore_*.py` — ad hoc iteration scripts (kept for reproducibility)

Scripts must be self-contained, documented, versioned, deterministic,
and have clear I/O (no hardcoded absolute paths).

## data/

KG extracts staged as files:
- Format: CSV (preferred) or TSV
- One file per extraction query
- Include metadata columns (organism, experiment_id, locus_tag)
- Name descriptively: `de_genes_med4_nitrogen.csv`, `catalase_orthologs.csv`

## results/

- **Tables:** CSV with clear headers and units
- **Figures:** PNG (300 DPI min) and PDF/SVG for publication; include axis labels, legends, titles
- **Statistics:** summary stats as CSV or in methods.md; full test results (statistic, p-value, effect size, CI)

## Git tracking convention

For per-step commit timing, manifest-commit coupling, and hard
gates, see [Step protocol](step-protocol.md).

### Per-analysis .gitignore

Each analysis gets its own `.gitignore` created during
scaffolding (before step 1). No blanket repo-level rules for
`analyses/*/data/` or `analyses/*/results/`. Tracking decisions
are explicit and logged in the notebook.

Default `.gitignore` template for a new analysis:
```
# Large intermediate data reproducible from KG
# (list specific files here, not blanket patterns)
__pycache__/
```

Everything else tracked by default. If a file should be ignored,
list it explicitly with a comment explaining why.

### What to track vs gitignore

| Track in git | Gitignore (list explicitly) |
|---|---|
| Signatures, scores, applied subsets (small CSVs) | Raw DE extracts (large CSVs, reproducible from KG) |
| Plots (PNG) | `__pycache__/` directories |
| Logs | |
| Manifests | |
| Notebook entries | |

Rule of thumb: if it's small and captures analytical decisions,
track it. If it's large and reproducible by rerunning a script,
gitignore it explicitly.

### Reproduction instructions

When large files are gitignored, the README must include
instructions for regenerating them (which scripts to run, in
what order). A reader who clones the repo should be able to
reproduce the full analysis.

## Log verbosity

Logs capture everything a researcher needs to verify the step
without rerunning it. This includes:
- Summary statistics (row counts, gene counts, filter funnel)
- Diagnostic traces (marker gene values at each stage, QC checks)
- Edge case results (tie-breaking outcomes, genes at
  classification boundaries)

Scripts with `--explore` flags should write their diagnostic
output to the log file as well as stdout. The log is the
persistent record; stdout is for the interactive session.

A log that says "N genes built" is insufficient. A log that shows
the filter funnel (how many started → how many passed each
filter → how many in the final set), marker gene traces, and
classification edge cases is sufficient.

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

Design-time decisions with rationale — from brainstorming and
planning. Execution-time decisions (continue/redo/adjust at the
end of each step) go in the notebook's `### Decision` section.

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

## References and citations

### Publication tracking
- Record publication DOIs from `list_publications` output
- Every experiment traces to a publication — maintain this link
- Cite the original study when presenting its expression data

### Methods referencing
- KG version/build date (from `kg_schema` or connection metadata)
- Tool versions and parameters used
- Statistical software versions
- Python and package versions

### Reference library
For analyses that produce a manuscript or report, maintain
`references.bib` or `references.md` with: original data
publications, methods references, KG/tool citations.
