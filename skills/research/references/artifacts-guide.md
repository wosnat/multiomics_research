# Artifacts guide

Every non-trivial analysis produces a directory under `analyses/`
with a standard structure. This ensures reproducibility and makes
it easy for collaborators (or future Claude sessions) to understand
what was done.

---

## Directory structure

```
analyses/{analysis_name}/
├── exploration/       # Exploration logs (one per research loop iteration)
├── data/              # Staged data from KG (CSV/TSV)
├── scripts/           # Python scripts (extract, analyze, explore)
├── results/           # Outputs: tables, figures, statistics
├── README.md          # Summary, key findings, file index
├── methods.md         # Publication-ready methods (living document)
├── gaps_and_friction.md  # KG/MCP/methodology issues log
└── references.md      # Data sources and citations
```

### Naming convention

`{analysis_name}` should be descriptive and lowercase with
underscores: `catalase_expression`, `nitrogen_stress_enrichment`,
`photosystem_conservation`.

---

## exploration/

Exploration logs — one file per iteration of the research loop.
Written **during** the session, not reconstructed afterward.

### Naming

`YYYY-MM-DD-{topic}.md` — one file per research question/iteration.
The topic should be descriptive: `orientation-and-scope`,
`cross-stress-specificity`, `timepoint-classification`.

### Format

Each exploration log follows this structure:

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
{Any KG/MCP/methodology issues encountered — also appended to gaps_and_friction.md}

## Next
{What question does this lead to? Or: ready to conclude.}
~~~

### Source tags

Every finding must be tagged with its source:
- `[KG]` — data directly from a KG query or script output
- `[interpretation]` — biological reasoning using intrinsic
  knowledge (literature context, mechanistic inference)
- `[gap]` — something the KG can't answer (missing data, missing
  annotations, missing organisms)

### Chain between iterations

The `## Next` section links iterations. File A's Next says
"classify all timepoints by stress specificity" → file B is
`2026-03-30-timepoint-classification.md`. The README should
list all exploration logs in order.

---

## gaps_and_friction.md

A running log of KG, MCP, and methodology issues discovered during
the analysis. This is a first-class output — updated during every
iteration of the research loop, not written retroactively.

Each entry is recorded in two places:
- The exploration log entry (for context — what were you doing
  when you hit this?)
- `gaps_and_friction.md` (aggregated backlog for tool/KG
  development)

### Format

Organize by category:

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
```

Add a process retrospective at the end of the analysis (exit gate):

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

---

## data/

KG extracts staged as files. These are the inputs to analysis
scripts.

- Format: CSV (preferred) or TSV
- One file per extraction query
- Include metadata columns (organism, experiment_id, locus_tag)
- Name files descriptively: `de_genes_med4_nitrogen.csv`,
  `catalase_orthologs.csv`

**How to extract:**

```python
#!/usr/bin/env python3
"""Extract differential expression data for analysis."""
import pandas as pd
from multiomics_explorer import differential_expression_by_gene

data = differential_expression_by_gene(
    organism="Prochlorococcus marinus MED4",
    experiment_ids=["exp_001", "exp_002"],
    significant_only=True,
)
df = pd.DataFrame(data["results"])
df.to_csv("data/de_genes_med4.csv", index=False)
print(f"Extracted {len(df)} rows (total_matching: {data['total_matching']})")
```

See [Python API guide](python-api-guide.md) for import paths,
return structure handling (especially nested fields), and common
mistakes to avoid.

Scripts in `data/` or `scripts/` — either works, but extraction
scripts should be clearly labeled.

---

## scripts/

Analysis scripts that read from `data/` and write to `results/`.

Requirements:
- **Self-contained** — run with `python scripts/analysis.py` from
  the analysis directory
- **Documented** — docstring explaining what the script does
- **Versioned** — log library versions in output or methods.md
- **Deterministic** — set random seeds where applicable
- **Clear I/O** — explicit input/output file paths, no hardcoded
  absolute paths

### Naming convention

- `scripts/extract_*.py` — data extraction from KG (reusable)
- `scripts/analyze_*.py` — computation that produces results
  (reusable)
- `scripts/explore_*.py` — ad hoc iteration scripts (may be
  throwaway, kept for reproducibility)

### Examples

- `scripts/extract_de_med4.py` — KG extraction
- `scripts/analyze_enrichment.py` — pathway enrichment analysis
- `scripts/analyze_volcano.py` — visualization
- `scripts/explore_cross_stress.py` — one-off cross-experiment
  comparison from a research loop iteration

Explore scripts are referenced from the exploration log
(`## Approach: ran scripts/explore_cross_stress.py`). If a pattern
recurs, promote the script to `extract_*` or `analyze_*`.

---

## results/

Outputs from analysis scripts.

### Tables
- CSV format for data tables
- Include headers with clear column names
- Include units where applicable (log2FC, padj, etc.)
- Name: `enrichment_results.csv`, `top_genes.csv`

### Figures
- PNG (for quick viewing) and PDF/SVG (for publication)
- Resolution: 300 DPI minimum for PNG
- Include axis labels, legends, titles
- Name: `volcano_nitrogen.png`, `heatmap_catalase.pdf`
- Consistent style across figures in the same analysis

### Statistics
- Summary statistics as CSV or in methods.md
- Full test results (test statistic, p-value, effect size,
  confidence interval)

---

## methods.md

A living document that grows alongside exploration. Updated
incrementally as findings become established — not written
retroactively at the end.

Publication-ready methods document. **Required sections:**

### Research question
What was asked, in precise terms. Not the user's original phrasing
but the operationalized question.

### Data scope
- Organisms and strains (with KG identifiers)
- Publications (DOIs)
- Experiments (IDs, conditions, time points)
- What was included and what was excluded, with reasons
- KG build date / version

### Gene selection
- How genes of interest were identified (search terms, ontology
  terms, homolog groups)
- Filters applied (significance, fold-change, quality)
- Gene counts at each filtering step
- Paralog handling

### Statistical methods
- Tests used (with parameters)
- Multiple testing correction method
- Background set definition (for enrichment)
- Software versions
- Significance thresholds

### Results summary
- Key findings with effect sizes and p-values
- Reference to specific output files in `results/`
- What was significant, what wasn't

### Limitations
- KG gaps (missing organisms, conditions, annotations)
- Statistical caveats (missing padj, small samples, borderline
  significance)
- What the analysis can and cannot conclude
- Known issues (duplicate entries, paralog conflation risks)

---

## references.md

Track all data sources and citations.

```markdown
# References

## Data sources
- KG build: [date], multiomics_biocypher_kg [version/commit]
- Neo4j: [connection info]

## Publications (data)
- [DOI] — [short description of what data was used from this study]
- [DOI] — ...

## Methods references
- DESeq2: Love et al. (2014) Genome Biology 15:550
- [enrichment tool]: [citation]
- [clustering method]: [citation]

## Software
- Python [version]
- multiomics_explorer [version]
- pandas [version]
- scipy [version]
- matplotlib [version]
```

---

## README.md

Brief summary for someone encountering the analysis directory.

```markdown
# {Analysis title}

## Question
[One-sentence research question]

## Key findings
- [Finding 1]
- [Finding 2]

## Files
| File | Description |
|---|---|
| `data/de_genes.csv` | Differential expression data from KG |
| `scripts/enrichment.py` | Pathway enrichment analysis |
| `results/volcano.png` | Volcano plot of DE genes |
| `methods.md` | Full methods documentation |

## How to reproduce
```bash
cd analyses/{analysis_name}
python scripts/extract_data.py   # stage data from KG
python scripts/enrichment.py     # run analysis
python scripts/volcano_plot.py   # generate figures
```

## Date
[Analysis date]
```

---

## When NOT to create artifacts

Simple lookups and quick answers don't need the full artifact
structure:
- "What is gene PMM0120?" — MCP lookup, answer in chat
- "How many genes respond to nitrogen?" — MCP summary, answer in
  chat
- "List the catalase genes" — MCP query, table in chat

Create artifacts when:
- The question requires statistical computation
- Multiple data extractions feed into an analysis
- Results need to be shared or revisited
- The analysis takes more than 2-3 tool calls
- Figures or publication-ready tables are needed
