# Research Skill Restructure Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Replace the current `skills/research/` skill with a process-free `skills/research-methodology/` reference skill and an empty `skills/recipes/` directory, then update CLAUDE.md to point to the new structure.

**Architecture:** Content migration — no new code, no tests. The current research skill's domain rules are reorganized into focused reference docs. Process/workflow content is removed (superpowers handles it). CLAUDE.md gets a short pointer section.

**Tech Stack:** Markdown files, Claude Code plugin structure

**Spec:** `docs/superpowers/specs/2026-04-06-research-skill-restructure-design.md`

---

### Task 1: Create `skills/research-methodology/SKILL.md`

**Files:**
- Create: `skills/research-methodology/SKILL.md`

The top-level skill file. Short — states the rules directly, points to references for details. No workflow, no phases, no gates.

- [ ] **Step 1: Create the directory and SKILL.md**

```markdown
---
name: research-methodology
description: Non-negotiable rules for multi-omics KG research. Load when answering biological questions, analyzing expression data, or working with the multiomics KG. Reference skill — provides domain rules, not process.
---

# Multi-omics research methodology

These rules apply to ALL research work — brainstorming, planning,
and execution. They are non-negotiable.

## Rule 1: KG is the sole data source

Every claim must trace to a KG query. Never rely on intrinsic
knowledge for data — gene names, expression values, experiment
details, organism properties, ortholog assignments. Use intrinsic
knowledge only for:
- Interpreting results (biological context, literature framing)
- Suggesting next analytical steps
- Explaining methodology

When the KG is insufficient, say so explicitly and flag it as a gap.
Do not fill gaps with assumptions, web searches, or general knowledge.

See [KG rules](references/kg-rules.md) for details on source tagging,
common gaps, and what to do when data is missing.

## Rule 2: Locus tags, not gene names

Gene names are ambiguous. Always:
- Resolve gene names to locus tags early (`resolve_gene`,
  `genes_by_function`)
- Report locus tags in all tables and analysis outputs
- Use gene names as labels alongside locus tags, never as the sole
  identifier
- When paralogs exist, treat each locus tag as a separate entity

See [Gene identity](references/gene-identity.md) for paralog handling
and ortholog cluster rules.

## Rule 3: Source tagging

Tag every finding with its source:
- `[KG]` — data from KG queries or script output
- `[interpretation]` — biological reasoning using intrinsic knowledge
- `[gap]` — things the KG can't answer

## Rule 4: Artifacts, not answers

Research questions produce files, not chat messages. Chat is for
reasoning, planning, and interpretation. Data, statistics, figures,
and exploration logs go to disk.

See [Artifacts guide](references/artifacts.md) for directory structure,
exploration log format, and file naming conventions.

## Rule 5: Scripts over chat reasoning

Computations go in `.py` files, not in chat responses. Data staged
to files before analysis. No manual steps — if it can't be scripted,
document it as a limitation.

See [Python API guide](references/python-api-guide.md) for imports,
return structure, and scripting patterns.

## Rule 6: Statistical rigor

See [Statistical rigor](references/statistical-rigor.md) for what
the KG provides, what you must compute in scripts, and what to flag.

## Rule 7: Don't hallucinate

See [Anti-hallucination](references/anti-hallucination.md) for
concrete failure modes and prevention patterns.

## References

- [KG rules](references/kg-rules.md) — source tagging, common gaps, KG insufficiency
- [Gene identity](references/gene-identity.md) — locus tags, paralogs, ortholog clusters
- [Artifacts guide](references/artifacts.md) — directory structure, exploration logs, methods.md
- [Anti-hallucination](references/anti-hallucination.md) — failure modes and prevention
- [Python API guide](references/python-api-guide.md) — scripting with multiomics_explorer
- [Statistical rigor](references/statistical-rigor.md) — KG statistics, what to compute, what to flag
```

- [ ] **Step 2: Verify file exists**

Run: `cat skills/research-methodology/SKILL.md | head -5`
Expected: the frontmatter header with `name: research-methodology`

- [ ] **Step 3: Commit**

```bash
git add skills/research-methodology/SKILL.md
git commit -m "feat: create research-methodology skill with top-level rules"
```

---

### Task 2: Create `references/kg-rules.md`

**Files:**
- Create: `skills/research-methodology/references/kg-rules.md`

Extracted from current SKILL.md §KG-is-the-sole-data-source and §Source-tagging. Also pulls domain rules from research-checklist.md (scope check items, what to verify about organisms/experiments/genes).

- [ ] **Step 1: Create the file**

```markdown
# KG rules

## KG is the sole data source

This is the top-level rule. Everything else follows from it.

Every claim must trace to a KG query. Never rely on intrinsic
knowledge for data — gene names, expression values, experiment
details, organism properties, ortholog assignments. Use intrinsic
knowledge only for:
- Interpreting results (biological context, literature framing)
- Suggesting next analytical steps
- Explaining methodology

## When the KG is insufficient

Say so explicitly and flag it as a gap. Do not fill gaps with
assumptions, web searches, or general knowledge. Common gaps:
- Missing organisms or strains not yet in the KG
- Annotations not loaded (e.g., some GO terms, regulatory elements)
- Expression data not available for a condition/study
- No adjusted p-values in the original dataset
- Missing metadata (timepoints, replicates, normalization method)
- No physiological data (growth curves, Fv/Fm, cell counts)

## Source tagging

Tag every finding with its source:
- `[KG]` — data directly from a KG query or script output
- `[interpretation]` — biological reasoning using intrinsic
  knowledge (literature context, mechanistic inference)
- `[gap]` — something the KG can't answer (missing data, missing
  annotations, missing organisms)

### Separating KG from literature

Clearly separate KG-derived findings from literature context:
- Use a distinct format: "**From the KG:** ... **From the
  literature (not verified in this KG):** ..."
- Never cite a specific paper from intrinsic knowledge as if
  it were a KG reference
- If the user needs literature support, say so — this system is
  a KG explorer, not a literature review tool

## Scoping: what to verify

Before diving into analysis, confirm the KG has what you need:
- Are the needed organisms present? (`list_organisms`)
- Are there experiments for the conditions of interest?
  (`list_experiments`)
- Are relevant publications available? (`list_publications`)
- Is the gene set complete? (check `total_matching` vs `returned`)
- Document any gaps: "The KG does not contain [X] data for [Y]"

## MCP vs Python API

| Use MCP when | Use Python API when |
|---|---|
| Browsing and orienting | Extracting full datasets |
| Checking gene counts and summaries | Running enrichment or statistics |
| Navigating between tools | Producing tables and figures |
| Quick lookups (< 50 results) | Any computation on gene lists |
| Working in chat mode | Result set > MCP limit |
| Cross-experiment overview (`gene_response_profile`) | Reshaping or aggregating data |
| | Building response matrices (`response_matrix()`) |
| | Comparing gene sets (`gene_set_compare()`) |

### When to switch from MCP to utilities

If you need to reshape, aggregate, or compute on MCP results, use
the analysis utilities first. Write scripts only for project-specific
logic not covered by the utilities. Do not attempt to work around
MCP limits in chat with one-off parsing.
```

- [ ] **Step 2: Commit**

```bash
git add skills/research-methodology/references/kg-rules.md
git commit -m "feat: add KG rules reference (source tagging, gaps, MCP vs API)"
```

---

### Task 3: Create `references/gene-identity.md`

**Files:**
- Create: `skills/research-methodology/references/gene-identity.md`

Extracted from current SKILL.md §Locus-tags and anti-hallucination.md §Category-1. Combines the rule with the failure modes.

- [ ] **Step 1: Create the file**

```markdown
# Gene identity rules

## Locus tags, not gene names

Gene names are ambiguous. Always:
- Resolve gene names to locus tags early (`resolve_gene`,
  `genes_by_function`)
- Report locus tags in all tables and analysis outputs
- Use gene names as labels alongside locus tags, never as the sole
  identifier

### Gene discovery tools

- By name/keyword: `resolve_gene`, `genes_by_function`
- By pathway/ontology: `search_ontology` → `genes_by_ontology`
- By homology: `search_homolog_groups` → `genes_by_homolog_group`
- By prior knowledge: start with `resolve_gene` for known names

## Paralog handling

When paralogs exist (multiple locus tags share a gene name):
- Treat each locus tag as a separate entity throughout
- Never aggregate expression data across locus tags that share
  a name without explicitly noting they are paralogs
- In tables: always include locus_tag column; add a suffix like
  "katA-1", "katA-2" when presenting results
- Use `gene_overview` on all locus tags to confirm they are
  distinct genes

**Example failure:** "katA is downregulated by starvation (log2FC
-7.8) then partially recovers at day 89 (log2FC +1.8)." Wrong —
the downregulation is ACZ81_02025, the "recovery" is ACZ81_11985.
Different genes, opposite behaviors, same name.

## Ortholog cluster rules

Genes in the same ortholog cluster are NOT interchangeable:
- Ortholog clusters are for cross-organism comparison of
  *conserved function*, not for within-organism gene equivalence
- When reporting per-gene expression, use locus tags, not cluster
  membership
- If a gene has no expression data, say so — don't substitute
  a cluster-mate's data

**Example failure:** MIT1002_02513 (katE) and MIT1002_03530 (katB)
are both in cluster 4644E. The analysis reported katE expression
values as "katB" because the cluster matched.

## Gene function claims

Gene function claims must come from `gene_overview`,
`gene_ontology_terms`, or `genes_by_function` output — not from
training data. If intrinsic knowledge disagrees with the KG
annotation, report the KG annotation and note the discrepancy.

Prefix intrinsic knowledge with: "Based on general knowledge
(not from this KG)..."
```

- [ ] **Step 2: Commit**

```bash
git add skills/research-methodology/references/gene-identity.md
git commit -m "feat: add gene identity reference (locus tags, paralogs, orthologs)"
```

---

### Task 4: Create `references/artifacts.md`

**Files:**
- Create: `skills/research-methodology/references/artifacts.md`

Reorganized from current `references/artifacts-guide.md`. Keeps artifact structure and format specs. Removes process language (when to create them, gates). Also absorbs references/citations section from current SKILL.md.

- [ ] **Step 1: Create the file**

```markdown
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
```

- [ ] **Step 2: Commit**

```bash
git add skills/research-methodology/references/artifacts.md
git commit -m "feat: add artifacts reference (directory structure, logs, methods.md)"
```

---

### Task 5: Create `references/statistical-rigor.md`

**Files:**
- Create: `skills/research-methodology/references/statistical-rigor.md`

Extracted from current SKILL.md §Statistical-rigor. Also absorbs relevant anti-hallucination patterns (cross-study p-value comparison, fabricating summary statistics).

- [ ] **Step 1: Create the file**

```markdown
# Statistical rigor

## What the KG provides

- log2FC and padj from the original DESeq2 analysis
- Per-study statistics — do not combine p-values across studies
- Precomputed summary counts (gene counts, significant counts)

## What you must compute in scripts

- Enrichment analysis (Fisher's exact, hypergeometric) with proper
  background set from the KG
- Multiple testing correction when comparing across conditions
- Effect size reporting (not just significance)
- Volcano plots, heatmaps with proper clustering

## What to flag

- Studies without adjusted p-values (report as fold-change only)
- Borderline significance (padj near 0.05)
- Duplicate KG entries for the same gene/condition (different
  contrasts — investigate, don't cherry-pick)
- Small sample sizes or missing replicates
- Conflation of paralogs or orthologs in aggregated results

## Strength of language vs strength of evidence

- padj < 0.001: "strongly significant" / "highly significant"
- padj < 0.01: "significant"
- padj < 0.05: "nominally significant" or "borderline significant"
- padj = 0.05: "at the conventional threshold" — flag explicitly
- No padj: "fold-change of X (no statistical test available)"
- Never write "consistent with [mechanism]" for data without
  statistical support — use "suggestive of" or "potentially
  related to"

## Cross-experiment comparison caveats

- Compare direction (up/down) and rank across platforms, not
  magnitude
- Compare log2FC magnitudes only within the same platform and study
- When presenting cross-experiment tables, include a platform
  column and note the caveat
- Use `gene_response_profile` rank fields (comparable across
  platforms) instead of raw fold changes for cross-study comparisons
- Expression data from different platforms (microarray, RNA-seq)
  or different statistical tests (Goldenspike, Rockhopper, DESeq2)
  are NOT directly comparable

## Cross-study p-value comparison

P-values from different studies have different designs, sample
sizes, and statistical power. They are not comparable:
- Never compare p-values across studies
- Compare effect sizes (log2FC) across studies if needed, but
  note the caveat
- For cross-study claims, use direction and rough magnitude,
  not precise p-value ranking

## Summary statistics

Summary statistics must be computed in scripts over complete data,
not eyeballed from MCP output:
- If you report a summary statistic, cite the script that computed it
- For quick estimates in chat, explicitly say "rough estimate from
  the top N results" and never present as a finding
```

- [ ] **Step 2: Commit**

```bash
git add skills/research-methodology/references/statistical-rigor.md
git commit -m "feat: add statistical rigor reference"
```

---

### Task 6: Move `anti-hallucination.md` and `python-api-guide.md`

**Files:**
- Move: `skills/research/references/anti-hallucination.md` → `skills/research-methodology/references/anti-hallucination.md`
- Move: `skills/research/references/python-api-guide.md` → `skills/research-methodology/references/python-api-guide.md`

These files are stable and move unchanged.

- [ ] **Step 1: Copy files to new location**

```bash
cp skills/research/references/anti-hallucination.md skills/research-methodology/references/anti-hallucination.md
cp skills/research/references/python-api-guide.md skills/research-methodology/references/python-api-guide.md
```

- [ ] **Step 2: Verify files match**

```bash
diff skills/research/references/anti-hallucination.md skills/research-methodology/references/anti-hallucination.md
diff skills/research/references/python-api-guide.md skills/research-methodology/references/python-api-guide.md
```

Expected: no output (files are identical)

- [ ] **Step 3: Commit**

```bash
git add skills/research-methodology/references/anti-hallucination.md skills/research-methodology/references/python-api-guide.md
git commit -m "feat: move anti-hallucination and python-api-guide to research-methodology"
```

---

### Task 7: Create empty recipes directory

**Files:**
- Create: `skills/recipes/.gitkeep`

Recipes don't exist yet. Create the directory with a `.gitkeep` so git tracks it.

- [ ] **Step 1: Create directory**

```bash
mkdir -p skills/recipes
touch skills/recipes/.gitkeep
```

- [ ] **Step 2: Commit**

```bash
git add skills/recipes/.gitkeep
git commit -m "feat: add empty recipes directory for future analysis protocols"
```

---

### Task 8: Update CLAUDE.md

**Files:**
- Modify: `CLAUDE.md`

Update the architecture section and key files table. Add the research methodology section. Update eval references.

- [ ] **Step 1: Update the Architecture section**

Replace the current architecture block:

```
multiomics_research/            # This repo (consumer)
  skills/research/              # Research methodology skill
  hooks/                        # Usage logging hooks
  evals/                        # Evaluation cases and runners
  scripts/                      # Analysis and eval scripts
  benchmarks/                   # Eval results (gitignored)
```

With:

```
multiomics_research/            # This repo (consumer)
  skills/research-methodology/  # Domain rules for KG research (reference skill)
  skills/recipes/               # On-demand analysis protocols (one skill per method)
  hooks/                        # Usage logging hooks
  evals/                        # Evaluation cases and runners
  scripts/                      # Analysis and eval scripts
  benchmarks/                   # Eval results (gitignored)
```

- [ ] **Step 2: Add Research methodology section**

Add after the Plugin Structure section, before Evaluation Framework:

```markdown
## Research methodology

When answering biological questions, analyzing expression data, or working
with the multiomics KG, invoke the `research-methodology` skill. It contains
the rules for KG usage, gene identity, artifact structure, and
anti-hallucination. These rules apply to all research work — brainstorming,
planning, and execution.

For specific analysis types, invoke the corresponding recipe skill (e.g.,
`enrichment`, `response-matrix`, `conservation`).

Use the superpowers workflow for research: brainstorm the question, write a
plan, execute with checkpoints. The methodology skill provides the domain
rules; superpowers provides the process discipline.
```

- [ ] **Step 3: Update Key Files table**

Replace:

```markdown
| `skills/research/SKILL.md` | Research methodology skill |
```

With:

```markdown
| `skills/research-methodology/SKILL.md` | Research methodology rules (reference skill) |
| `skills/recipes/` | On-demand analysis protocols |
```

- [ ] **Step 4: Update eval trigger path**

Replace:

```bash
python scripts/run_trigger_eval.py evals/trigger_evals.json --skill-path skills/research
```

With:

```bash
python scripts/run_trigger_eval.py evals/trigger_evals.json --skill-path skills/research-methodology
```

- [ ] **Step 5: Commit**

```bash
git add CLAUDE.md
git commit -m "docs: update CLAUDE.md for research-methodology skill restructure"
```

---

### Task 9: Remove old `skills/research/` directory

**Files:**
- Delete: `skills/research/` (entire directory)

Only do this after all new files are confirmed in place.

- [ ] **Step 1: Verify new skill is complete**

```bash
ls skills/research-methodology/SKILL.md
ls skills/research-methodology/references/kg-rules.md
ls skills/research-methodology/references/gene-identity.md
ls skills/research-methodology/references/artifacts.md
ls skills/research-methodology/references/statistical-rigor.md
ls skills/research-methodology/references/anti-hallucination.md
ls skills/research-methodology/references/python-api-guide.md
```

Expected: all files exist

- [ ] **Step 2: Remove old directory**

```bash
git rm -r skills/research/
```

- [ ] **Step 3: Commit**

```bash
git commit -m "refactor: remove old skills/research/ (replaced by research-methodology)"
```

---

### Task 10: Update plugin.json description

**Files:**
- Modify: `.claude-plugin/plugin.json`

Update the description to reflect the new structure.

- [ ] **Step 1: Update description**

Change the `description` field from:

```
"Research methodology and evaluation tools for the Prochlorococcus/Alteromonas multi-omics knowledge graph. Provides a research skill that guides Claude through systematic multi-omics analysis: orientation, gene work, expression analysis, with anti-hallucination safeguards and reproducible artifact generation."
```

To:

```
"Research methodology and evaluation tools for the Prochlorococcus/Alteromonas multi-omics knowledge graph. Provides domain rules for KG-based research (source tagging, gene identity, anti-hallucination, artifact structure) and on-demand analysis recipes. Designed to work with the superpowers framework for process discipline."
```

- [ ] **Step 2: Commit**

```bash
git add .claude-plugin/plugin.json
git commit -m "docs: update plugin.json description for restructured skills"
```
