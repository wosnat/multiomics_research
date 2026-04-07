# Artifacts guide

Research questions produce files, not chat messages. Chat is for
reasoning, planning, and interpretation. Data, statistics, figures,
and exploration logs go to disk.

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

**Naming:** `{analysis_name}` should be descriptive and lowercase
with underscores: `catalase_expression`,
`nitrogen_stress_enrichment`, `photosystem_conservation`.

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
