---
name: research
description: Academic research methodology for multi-omics analysis. Activate when answering biological questions, analyzing expression data, or producing research artifacts.
argument-hint: "[research question or 'review' to audit an existing analysis]"
---

# Research methodology

See [research checklist](references/research-checklist.md) for the
step-by-step protocol, [artifacts guide](references/artifacts-guide.md)
for output structure, and [anti-hallucination](references/anti-hallucination.md)
for concrete failure modes to avoid.

## KG is the sole data source

This is the top-level rule. Everything else follows from it.

Every claim must trace to a KG query. Never rely on intrinsic
knowledge for data — gene names, expression values, experiment
details, organism properties, ortholog assignments. Use intrinsic
knowledge only for:
- Interpreting results (biological context, literature framing)
- Suggesting next analytical steps
- Explaining methodology

**When the KG is insufficient**, say so explicitly and flag it as
a gap. Do not fill gaps with assumptions, web searches, or general
knowledge. Common gaps:
- Missing organisms or strains not yet in the KG
- Annotations not loaded (e.g., some GO terms, regulatory elements)
- Expression data not available for a condition/study
- No adjusted p-values in the original dataset
- Missing metadata (timepoints, replicates, normalization method)
- No physiological data (growth curves, Fv/Fm, cell counts)

This rule is made operational through **source tagging** in the
research loop: every finding is tagged `[KG]`, `[interpretation]`,
or `[gap]` (see Stage 2).

## Core principles

### 1. Reproducibility over convenience

Every analysis must be reproducible without access to the Claude
conversation. This means:
- **Scripts over chat reasoning** — computations go in `.py` files,
  not in chat responses
- **Data staged to files** — KG extracts saved as CSV/TSV before
  analysis, using the Python package (not MCP tool output)
- **Methods documented** — every decision recorded in `methods.md`
  as it happens, not retroactively
- **No manual steps** — if it can't be scripted, document it as a
  limitation

### 2. Locus tags, not gene names

Gene names are ambiguous. Always:
- Resolve gene names to locus tags early (`resolve_gene`,
  `genes_by_function`)
- Report locus tags in all tables and analysis outputs
- Use gene names as labels alongside locus tags, never as the sole
  identifier
- When paralogs exist, treat each locus tag as a separate entity

### 3. Artifacts, not answers

Research questions produce files, not chat messages. Chat is for
reasoning, planning, and interpretation. Data, statistics, figures,
and **exploration logs** go to disk. See the
[artifacts guide](references/artifacts-guide.md).

---

## The 2-stage workflow

### Stage 1: Orientation (linear, run once)

Establish scope and identify the initial gene set before entering
the research loop.

**Scope:**
- `list_organisms` — what species/strains are available
- `list_publications` — what studies cover the question
- `list_experiments` — which experiments, conditions, time points

**Gene identification:**
- Discovery: `resolve_gene`, `genes_by_function`, `genes_by_ontology`,
  `genes_by_homolog_group`
- Details: `gene_overview`, `gene_ontology_terms`, `gene_homologs`
- Annotation: `search_ontology`, `search_homolog_groups`

**Outputs:**
- `methods.md` stub with: research question, data scope, gene
  identification approach
- First exploration log entry
  (`exploration/YYYY-MM-DD-orientation.md`)
- `gaps_and_friction.md` initialized (even if empty)

**Gate:** Can the KG answer this question?
- Are the needed organisms present?
- Are there experiments for the conditions of interest?
- Are all genes resolved to locus tags? Are paralogs identified?
- Is the gene set complete (check `total_matching` vs `returned`)?
- If key data is missing, report the gap before proceeding.

### Stage 2: Research loop (iterative, run N times)

Each iteration follows: **Question → Explore → Log → Assess**.
The loop exits when the research question is answered or the
remaining gaps can't be addressed with available data.

#### Question

- State what you're testing as a testable claim or comparison
  (not open-ended)
- Mark the question type:
  - `hypothesis` — have a prediction to test
  - `exploratory` — no prediction, looking for patterns
  - `follow-up` — triggered by a previous iteration's findings

#### Explore

- Use MCP tools for browsing and quick lookups
- Use `gene_response_profile` for cross-experiment gene-level
  summaries (when available; fall back to manual aggregation via
  API extraction otherwise)
- Write scripts for any computation or reshaping — save to
  `scripts/explore_*.py` immediately rather than doing one-offs
  in chat
- **Cross-experiment comparison caveat:** when comparing across
  platforms (microarray vs RNA-seq), compare direction (up/down)
  and rank, not fold-change magnitude. Include a platform column
  in cross-experiment tables. See
  [anti-hallucination §2.5](references/anti-hallucination.md).

#### Log

Write to `exploration/YYYY-MM-DD-{topic}.md` **during** the
iteration, not after. One file per iteration.

Tag every finding with its source:
- `[KG]` — data from KG queries
- `[interpretation]` — biological reasoning using intrinsic
  knowledge
- `[gap]` — things the KG can't answer

See exploration log format in the
[artifacts guide](references/artifacts-guide.md).

#### Assess

Classify each finding:
- **Established** — consistent across experiments, statistically
  supported
- **Preliminary** — single experiment, or no statistical support
- **Speculative** — interpretation beyond data

Self-check: "Did I use any knowledge that didn't come from the KG?
If yes, is it tagged `[interpretation]`?"

- Update the working hypothesis
- Log any gaps/friction encountered → append to
  `gaps_and_friction.md`
- When findings are classified as `established`, update the
  relevant section of `methods.md`
- Decide: next question, or conclude?

### Synthesis (after the loop)

- Finalize `methods.md` — ensure coherence and publication-readiness
- Write or update `README.md` — summary, key findings, file index,
  exploration log links
- Produce publication artifacts (figures, tables) in `results/`
- Populate `references.md`

---

## MCP vs Python package

| Use MCP when | Use Python package when |
|---|---|
| Browsing and orienting | Extracting full datasets |
| Checking gene counts and summaries | Running enrichment or statistics |
| Navigating between tools | Producing tables and figures |
| Quick lookups (< 50 results) | Any computation on gene lists |
| Working in chat mode | Result set > MCP limit |
| Cross-experiment overview (`gene_response_profile`) | Reshaping or aggregating data |

**When to switch:** If you need to reshape, aggregate, or compute
on MCP results, write a script immediately. Save it to
`scripts/explore_*.py` and reference it from the exploration log.
Do not attempt to work around MCP limits in chat with one-off
parsing.

**Script naming convention:**
- `scripts/extract_*.py` — data extraction from KG (reusable)
- `scripts/analyze_*.py` — computation that produces results
  (reusable)
- `scripts/explore_*.py` — ad hoc iteration scripts (may be
  throwaway, kept for reproducibility)

Explore scripts can be promoted into the analysis utilities module
if the pattern recurs.

```python
from multiomics_explorer import differential_expression_by_gene
data = differential_expression_by_gene(
    experiment_ids=["..."], significant_only=True)
# limit=None by default in api/ -> all rows
pd.DataFrame(data["results"]).to_csv("data/de_genes.csv", index=False)
```

---

## Statistical rigor

### What the KG provides
- log2FC and padj from the original DESeq2 analysis
- Per-study statistics — do not combine p-values across studies
- Precomputed summary counts (gene counts, significant counts)

### What you must do in scripts
- Enrichment analysis (Fisher's exact, hypergeometric) with proper
  background set from the KG
- Multiple testing correction when comparing across conditions
- Effect size reporting (not just significance)
- Volcano plots, heatmaps with proper clustering

### What to flag
- Studies without adjusted p-values (report as fold-change only)
- Borderline significance (padj near 0.05)
- Duplicate KG entries for the same gene/condition (different
  contrasts — investigate, don't cherry-pick)
- Small sample sizes or missing replicates
- Conflation of paralogs or orthologs in aggregated results

---

## References and citations

### Publication tracking
- Record publication DOIs from `list_publications` output
- Every experiment traces to a publication — maintain this link
  throughout the analysis
- Cite the original study when presenting its expression data

### Methods referencing
- KG version/build date (from `kg_schema` or connection metadata)
- Tool versions and parameters used
- Statistical software versions (scipy, statsmodels, etc.)
- Python and package versions

### Reference library
For analyses that produce a manuscript or report:
- Maintain a `references.bib` or `references.md` in the analysis
  directory
- Include: original data publications, methods references (DESeq2,
  enrichment tools), KG/tool citations
- Link each claim in the analysis to its data source (KG query)
  and original publication

---

## Review protocol

Every analysis should be reviewed. Use
`/multiomics-research:research review` or apply this checklist:

1. **Data provenance** — Can every number be traced to a KG query
   or script output?
2. **Source tagging** — Are findings tagged `[KG]`,
   `[interpretation]`, or `[gap]` in exploration logs?
3. **Gene identity** — Are locus tags used throughout? Are paralogs
   separated?
4. **Completeness** — Were full datasets extracted for statistical
   analyses, or were truncated MCP results used?
5. **Confidence classification** — Are findings marked as
   established, preliminary, or speculative?
6. **Statistics** — Were p-values taken from the KG (original study)
   or computed? If computed, is the method documented?
7. **Gaps declared** — Are KG limitations and missing data explicitly
   stated in both exploration logs and `gaps_and_friction.md`?
8. **Reproducibility** — Can someone run `python scripts/*.py` and
   reproduce the results?
9. **Exploration trail** — Do the exploration logs form a followable
   chain (each entry's `## Next` points to the next)?
10. **Artifacts** — Do all claimed files exist and contain what
    they should?
