# Research Skill Iterative Redesign — Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Restructure the research skill from linear phases to orientation + iterative research loop with structured exploration logging.

**Architecture:** Four markdown files are modified: SKILL.md (full rewrite), research-checklist.md (rewrite), artifacts-guide.md (extend), anti-hallucination.md (extend). No code, no tests — this is a skill/documentation change. Each task modifies one file and commits.

**Tech Stack:** Markdown (Claude Code skill format)

**Spec:** `docs/superpowers/specs/2026-03-31-research-skill-iterative-redesign.md`

---

### Task 1: Rewrite SKILL.md

**Files:**
- Modify: `skills/research/SKILL.md` (full rewrite, preserve frontmatter)

The current SKILL.md has: frontmatter, core principles (4), 3-phase workflow, MCP vs Python, statistical rigor, references/citations, review protocol. The new structure replaces the principles and workflow sections while keeping the technical reference sections largely intact.

- [ ] **Step 1: Rewrite SKILL.md with the new structure**

Replace the entire file content (keeping the frontmatter lines 1-5 unchanged) with:

```markdown
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
```

- [ ] **Step 2: Read the new file and verify structure**

Read `skills/research/SKILL.md` and verify:
- Frontmatter is intact (name, description, argument-hint)
- Sections appear in spec order: KG rule → Core principles → 2-stage workflow → MCP vs Python → Statistical rigor → References → Review protocol
- No leftover content from the old file
- All cross-references to reference docs are correct paths

- [ ] **Step 3: Commit**

```bash
git add skills/research/SKILL.md
git commit -m "Rewrite research skill: linear phases → iterative research loop

Restructures from 5 linear phases to orientation + Question/Explore/Log/Assess
loop. Elevates KG-as-sole-source to top-level rule. Adds source tagging,
confidence classification, and exploration logging to the core workflow."
```

---

### Task 2: Rewrite research-checklist.md

**Files:**
- Modify: `skills/research/references/research-checklist.md` (full rewrite)

The current checklist has: Before starting, Phase 1-5 with checkboxes, Review. The new checklist has: Before starting (unchanged), Orientation checklist (merged Phase 1+2), Research loop checklist (per iteration), Synthesis checklist, Review (updated).

- [ ] **Step 1: Rewrite research-checklist.md**

Replace the entire file content with:

```markdown
# Research checklist

Step-by-step protocol for answering a biological research question
using the multiomics KG.

---

## Before starting

- [ ] Restate the research question precisely
- [ ] Identify what kind of question it is:
  - **Lookup** — gene identity, function, annotations → MCP only
  - **Survey** — one gene across conditions, or one condition across
    genes → MCP + maybe extraction
  - **Comparative** — two conditions, two organisms, pathway
    comparison → extraction + statistics required
  - **Discovery** — enrichment, clustering, unbiased screen →
    full extraction + computational pipeline required
- [ ] For comparative/discovery: create the analysis directory
  structure (see artifacts guide)

---

## Orientation (run once)

### Scope

- [ ] `list_organisms` — confirm target organisms are in the KG
- [ ] `list_publications` — find relevant studies
- [ ] `list_experiments(summary=true)` — understand available
  conditions, time points, omics types
- [ ] `list_experiments` with filters — get specific experiment IDs

### Gene identification

- [ ] Identify genes of interest using appropriate discovery tools:
  - By name/keyword: `resolve_gene`, `genes_by_function`
  - By pathway/ontology: `search_ontology` → `genes_by_ontology`
  - By homology: `search_homolog_groups` → `genes_by_homolog_group`
  - By prior knowledge: start with `resolve_gene` for known names
- [ ] Resolve all gene names to locus tags
- [ ] Check for paralogs — multiple locus tags per gene name?
  - If yes: `gene_overview` on all locus tags, treat separately
- [ ] `gene_overview` — confirm annotation types, expression data
  availability, ortholog coverage

### Orientation outputs

- [ ] Create `methods.md` stub (research question, data scope, gene
  set)
- [ ] Write first exploration log
  (`exploration/YYYY-MM-DD-orientation.md`)
- [ ] Initialize `gaps_and_friction.md`

### Gate: scope check

- [ ] Are the needed organisms present?
- [ ] Are there experiments for the conditions of interest?
- [ ] Are all genes resolved to locus tags? Paralogs identified?
- [ ] Is the gene set complete? (check `total_matching` vs
  `returned`)
- [ ] Document any gaps: "The KG does not contain [X] data for [Y]"

---

## Research loop (per iteration)

Run this checklist for each iteration of the research loop.

### Question

- [ ] State a testable claim or comparison (not open-ended)
- [ ] Mark type: `hypothesis`, `exploratory`, or `follow-up`

### Explore

- [ ] Use MCP for browsing, Python scripts for computation
- [ ] For cross-experiment comparison: use `gene_response_profile`
  or extract via API and aggregate in a script
- [ ] Save any ad hoc computation to `scripts/explore_*.py`
- [ ] If comparing across platforms: note caveat, compare direction
  and rank not magnitude

### Log

- [ ] Write exploration log entry
  (`exploration/YYYY-MM-DD-{topic}.md`) **during** the iteration
- [ ] Tag every finding: `[KG]`, `[interpretation]`, or `[gap]`

### Assess

- [ ] Classify findings: `established`, `preliminary`, or
  `speculative`
- [ ] Self-check: did I use knowledge not from the KG? Is it tagged?
- [ ] Append any gaps/friction to `gaps_and_friction.md`
- [ ] Update `methods.md` if established findings changed scope
- [ ] Decide: next question, or conclude?

---

## Synthesis (after the loop)

- [ ] Finalize `methods.md` — coherent, publication-ready
- [ ] Write or update `README.md` — summary, findings, file index,
  exploration log links
- [ ] Produce publication artifacts in `results/`
- [ ] Populate `references.md` or `references.bib`

---

## Review (self or peer)

- [ ] Every number traces to a named KG query or script output
- [ ] Findings are tagged with source (`[KG]`, `[interpretation]`,
  `[gap]`) in exploration logs
- [ ] No gene name used without its locus tag
- [ ] Paralogs and orthologs handled correctly
- [ ] Truncation metadata respected (not counting from limited rows)
- [ ] Findings classified (established / preliminary / speculative)
- [ ] Statistical methods appropriate and documented
- [ ] KG gaps and limitations stated in exploration logs and
  `gaps_and_friction.md`
- [ ] All artifact files exist and are consistent with the narrative
- [ ] Exploration logs form a followable chain (`## Next` links)
- [ ] Claims match the data — no over-interpretation
```

- [ ] **Step 2: Read the new file and verify**

Read `skills/research/references/research-checklist.md` and verify:
- Before starting section preserved
- Orientation section merges old Phase 1 + Phase 2
- Research loop section has per-iteration checklist matching spec (Question, Explore, Log, Assess)
- Synthesis section covers post-loop finalization
- Review section includes new items (source tagging, confidence, exploration trail)

- [ ] **Step 3: Commit**

```bash
git add skills/research/references/research-checklist.md
git commit -m "Rewrite research checklist: orientation + per-iteration loop checklist

Replaces linear Phase 1-5 checklist with orientation (run once) and
research loop (per iteration) structure. Adds source tagging, confidence
classification, and exploration log checks."
```

---

### Task 3: Extend artifacts-guide.md

**Files:**
- Modify: `skills/research/references/artifacts-guide.md`

Three additions: updated directory structure (with `exploration/` and `gaps_and_friction.md`), exploration log format section, scripts naming convention. The rest of the file is unchanged.

- [ ] **Step 1: Update the directory structure**

In `artifacts-guide.md`, find the directory structure block (lines 13-19):

```
analyses/{analysis_name}/
├── data/              # Staged data from KG (CSV/TSV)
├── scripts/           # Python analysis scripts
├── results/           # Outputs: tables, figures, statistics
├── README.md          # Summary, key findings, file index
├── methods.md         # Publication-ready methods document
└── references.md      # Data sources and citations
```

Replace with:

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

- [ ] **Step 2: Add exploration/ section after the directory structure**

After the `### Naming convention` subsection (line 26), insert a new section:

```markdown
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
```

- [ ] **Step 3: Add gaps_and_friction.md section**

After the new `exploration/` section, insert:

```markdown
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
```

- [ ] **Step 4: Update the scripts/ section**

Find the `## scripts/` section. After the existing "Example patterns" list (lines 77-80):

```
- `scripts/extract_data.py` — KG extraction
- `scripts/enrichment.py` — pathway enrichment analysis
- `scripts/volcano_plot.py` — visualization
- `scripts/clustering.py` — gene clustering
```

Replace with:

```markdown
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
```

- [ ] **Step 5: Update methods.md section to note it's a living document**

Find the `## methods.md` section header (line 109). After the header line, before "Publication-ready methods document", insert:

```markdown
A living document that grows alongside exploration. Updated
incrementally as findings become established — not written
retroactively at the end.
```

- [ ] **Step 6: Read the modified file and verify**

Read `skills/research/references/artifacts-guide.md` and verify:
- Directory structure includes `exploration/`, `gaps_and_friction.md`
- Exploration log format section present with source tags and chain description
- gaps_and_friction.md section present
- Scripts naming convention updated
- methods.md described as living document
- Rest of file unchanged (data/, results/, README.md, references.md, "When NOT to create artifacts")

- [ ] **Step 7: Commit**

```bash
git add skills/research/references/artifacts-guide.md
git commit -m "Extend artifacts guide: exploration logs, gaps log, script naming

Adds exploration/ directory format with source tagging, gaps_and_friction.md
as standard artifact, and extract/analyze/explore script naming convention."
```

---

### Task 4: Extend anti-hallucination.md

**Files:**
- Modify: `skills/research/references/anti-hallucination.md`

Add pattern 2.5 (cross-experiment aggregation without caveats) after the existing pattern 2.4.

- [ ] **Step 1: Add pattern 2.5**

In `anti-hallucination.md`, find the end of section 2.4 (after line 131, the last line of the "Prevention" list for 2.4). Before the `---` separator and `## Category 3`, insert:

```markdown

### 2.5 Cross-experiment aggregation without caveats

**What happens:** Expression data from different platforms
(microarray, RNA-seq) or different statistical tests (Goldenspike,
Rockhopper, DESeq2) are compared as if they are directly
comparable. "Gene X has log2FC 3.2 in study A and 1.1 in study B,
so the response is stronger in A."

**Real example (nitrogen stress analysis):** Tolonen 2006
(microarray, Goldenspike) and Read 2017 (RNA-seq, Rockhopper)
fold changes for the same genes were compared. Different platforms
have different dynamic ranges — a log2FC of 3 on microarray is
not the same as log2FC 3 on RNA-seq.

**Prevention:**
- Compare direction (up/down) and rank across platforms, not
  magnitude
- Compare log2FC magnitudes only within the same platform and study
- When presenting cross-experiment tables, include a platform
  column and note the caveat
- Use `gene_response_profile` rank fields (comparable across
  platforms) instead of raw fold changes for cross-study
  comparisons
```

- [ ] **Step 2: Read the modified file and verify**

Read `skills/research/references/anti-hallucination.md` and verify:
- Pattern 2.5 appears between 2.4 and Category 3
- The real example references the nitrogen stress analysis
- Prevention bullets include the rank vs magnitude guidance
- Rest of file unchanged

- [ ] **Step 3: Commit**

```bash
git add skills/research/references/anti-hallucination.md
git commit -m "Add anti-hallucination pattern: cross-experiment aggregation caveats

New pattern 2.5 warns against comparing fold-change magnitudes across
platforms. Based on nitrogen stress analysis experience with microarray
vs RNA-seq comparisons."
```

---

### Task 5: Final verification

- [ ] **Step 1: Verify all four files are internally consistent**

Read all four files and check cross-references:
- `SKILL.md` references `references/research-checklist.md` — path correct?
- `SKILL.md` references `references/artifacts-guide.md` — path correct?
- `SKILL.md` references `references/anti-hallucination.md` (§2.5) — section exists?
- `research-checklist.md` mentions exploration log format — consistent with artifacts guide?
- Source tags (`[KG]`, `[interpretation]`, `[gap]`) described consistently across all files?
- Confidence levels (`established`, `preliminary`, `speculative`) described consistently?

- [ ] **Step 2: Verify the skill loads correctly**

```bash
# Check skill structure is intact
ls -la skills/research/SKILL.md skills/research/references/
# Check frontmatter is valid
head -5 skills/research/SKILL.md
```

Expected output:
```
---
name: research
description: Academic research methodology for multi-omics analysis. Activate when answering biological questions, analyzing expression data, or producing research artifacts.
argument-hint: "[research question or 'review' to audit an existing analysis]"
---
```
