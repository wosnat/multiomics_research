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
├── exploration/       # Exploration logs (one per research iteration)
├── data/              # Staged data from KG (CSV/TSV)
├── scripts/           # Python scripts (extract, analyze, explore)
├── results/           # Outputs: tables, figures, statistics
├── README.md          # Summary, key findings, file index
├── methods.md         # Publication-ready methods (living document)
├── gaps_and_friction.md  # KG/MCP/methodology issues log
└── references.md      # Data sources and citations
```

**Naming:** `{analysis_name}` should be descriptive and lowercase
with underscores: `catalase_expression`,
`nitrogen_stress_enrichment`, `photosystem_conservation`.

## Exploration logs

One file per research iteration. Written **during** the session,
not reconstructed afterward.

**Naming:** `YYYY-MM-DD-{topic}.md`

**Format:**

~~~markdown
# YYYY-MM-DD: {Topic}

## Question
{Testable question or comparison. Mark: hypothesis / exploratory / follow-up}

## Approach
{What queries/scripts will be used. Brief — not a methods section.}

## Findings
{Results, tables, observations. Each finding tagged:}
- [KG] cynA (PMM0370) is significant only under nitrogen_stress across all 30 MED4 experiments
- [interpretation] This suggests cynA is a reliable N-specific marker
- [gap] No physiological data (Fv/Fm) to confirm cell viability at 48h

## Assessment
{What did we learn? Classify findings:}
- **Established:** {consistent across experiments, statistically supported}
- **Preliminary:** {single experiment, or no stats}
- **Speculative:** {interpretation beyond data}

## Gaps and friction
{Any issues encountered — also appended to gaps_and_friction.md}

## Next
{What question does this lead to? Or: ready to conclude.}
~~~

**Chain between iterations:** The `## Next` section links iterations.
The README should list all exploration logs in order.

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
