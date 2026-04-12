# Step Protocol Enforcement — Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Consolidate per-step process obligations into a single step-protocol reference, add a chat-capture gate, and fix git commit/ignore discipline in the research-methodology skill.

**Architecture:** New `references/step-protocol.md` becomes the single source of truth for per-step timing and gates. Existing references (research-notebook.md, artifacts.md) have enforcement content removed to avoid duplication. SKILL.md Rule 8 updated to point to both step-protocol and research-notebook.

**Spec:** [2026-04-12-step-protocol-enforcement-design.md](../specs/2026-04-12-step-protocol-enforcement-design.md)

**Note:** This is a skill text update — no code, no tests. All files are markdown in the research-methodology skill directory. "Verify" steps mean reading the file and confirming content.

**Skill root:** `.claude/skills/research-methodology/`

---

### Task 1: Create `references/step-protocol.md` [Spec §1, §2, §3]

**Files:**
- Create: `.claude/skills/research-methodology/references/step-protocol.md`

- [ ] **Step 1: Create step-protocol.md**

```markdown
# Step protocol

The step cycle (do → show → explore → decide) is defined in
[research-notebook.md](research-notebook.md). This document owns
**when things happen and what gates enforce them**. Follow it for
every step that produces data or analytical output.

## Commit structure

Two commits per step:
- **Commit 1** (end of "do"): script + outputs + log + manifest
  updates
- **Commit 2** (end of "explore"): notebook entry with QC tables
  and chat-capture section

Show and explore produce ONE notebook entry together — QC section
from show, then chat-capture section from explore — committed
once at the end of explore.

If the researcher requests a redo, the exploration reasoning is
preserved separately from the artifacts that get replaced.

## Before starting a step

- Confirm previous step's artifacts are committed (script +
  outputs + log + manifests — Commit 1)
- Confirm previous step's notebook entry is committed (including
  chat-capture — Commit 2)

## "do" phase

1. Write or update the script
2. Run script, capture outputs and log
3. Update `DATA_MANIFEST.md` or `RESULTS_MANIFEST.md` with new
   files (see [artifacts.md](artifacts.md) for manifest format)
4. Git commit: script + outputs + log + manifest updates

## "show" phase (begins the notebook entry)

Present QC diagnostics to the researcher. Write the QC section
of the notebook entry, including:

- Summary tables of outputs — the same tables shown in chat,
  written to the notebook as markdown tables, not prose
  paraphrases
- Links to figures produced
  (`![description](../results/fig.png)`)
- Row counts, gene counts, key statistics

See [research-notebook.md — QC checkpoint types](research-notebook.md)
for what to show per step type.

## "explore" phase (continues the same notebook entry)

Interactive walkthrough with the researcher. Append the
chat-capture section to the notebook entry AS exploration
happens, not after.

### Chat-capture format

```markdown
## Chat exploration

**Q: [researcher's question, as asked]**
Data: [what was looked up / computed to answer it]
Finding: [what the data showed, with concrete numbers]
Impact: [how this affects interpretation or next steps]

**Q: [next question]**
...
```

**What gets captured:** Questions that surfaced data points,
changed interpretation, or influenced decisions. Any exchange
that produced a number, a comparison, or a "therefore."

**What does NOT get captured:** Clarifications ("what does column
X mean"), formatting/typo fixes, tool-call mechanics.

Why not verbatim chat? Chat is noisy with tool calls, retries,
and formatting. The value is the Q → data → finding → impact
chain, not the raw transcript.

**GATE 3** applies here — see below.

## "decide" phase (closes the notebook entry)

Researcher says continue, redo, or adjust. Decision logged in
notebook entry. Then git commit (Commit 2): the completed
notebook entry.

## Redo path

When the researcher says "redo with X":

1. **do:** Update script, rerun, update manifests → NEW commit
   (never amend previous)
2. **show:** New QC presentation — same requirements as first
   pass (summary tables, figure links, counts)
3. **explore:** New walkthrough with chat-capture → APPEND new
   notebook entry (never revise previous entries)
4. **decide:** Researcher reviews. The notebook entry explicitly
   lists downstream steps that consumed the changed output.
   Researcher decides whether to cascade redo or accept as-is.
   Decision logged in notebook.

The failed attempt's notebook entry stays as a record of what
was tried and why it was rejected. Append-only.

## Hard gates

Each gate explains what goes wrong without it (citing real
failures from past analyses), then states the stop condition.

### GATE 1: Step boundary

B1's notebook was written retroactively — exploration reasoning
was lost and couldn't be verified against the actual data state
at the time.

**Do not start step N+1 until step N's notebook entry (including
chat-capture section) is committed.**

### GATE 2: Manifest currency

B1's manifests were updated in bulk at the end. By that point,
file descriptions were reconstructed from memory rather than
written when the data was fresh.

**Do not commit a script's outputs without updating the relevant
manifest in the same commit.**

### GATE 3: Chat-capture

B1's most valuable findings came from researcher questions during
explore ("is this dominated by catch-all categories?", "what
about Steglich's low power?"). These drove analytical decisions
but were lost from the chat context.

**Do not enter "decide" phase until the chat-capture section
exists in the notebook entry.**

## Git discipline

### Per-analysis .gitignore

Each analysis gets its own `.gitignore` created during
scaffolding (before step 1), with explicit decisions logged in
the notebook. Default template:

```
# Large intermediate data reproducible from KG
# (list specific files here, not blanket patterns)
__pycache__/
```

Everything else tracked by default. If a file should be ignored,
list it explicitly with a comment explaining why.

### Scaffolding commit

When creating an analysis directory, the initial scaffolding
(directory structure, empty manifests, `.gitignore`, notebook
stub) gets its own commit before step 1 begins. This establishes
tracking decisions upfront.

### Commit discipline

- Two commits per step (do + explore)
- Redo: new commits, never amend
- Each commit is self-contained: includes the artifacts it
  produces AND the manifest/notebook updates for those artifacts
```

- [ ] **Step 2: Verify step-protocol.md**

Read the file. Confirm:
- All three gates present with "why" and stop condition
- Chat-capture format includes Q/Data/Finding/Impact
- Per-phase obligations match spec (do, show, explore, decide)
- Redo path includes downstream dependency flagging
- Git discipline section covers per-analysis .gitignore and scaffolding commit

---

### Task 2: Update `references/research-notebook.md` [Spec §2, §4]

Add chat-capture entry format. Remove step-cycle enforcement content that moved to step-protocol.md. Keep notebook format, entry templates, QC checkpoint types, code lifecycle.

**Files:**
- Modify: `.claude/skills/research-methodology/references/research-notebook.md`

- [ ] **Step 1: Replace the notebook-commit gate section (lines 1-19)**

Replace the opening "Notebook-commit gate" section with a pointer to step-protocol.md:

Old (lines 1-19):
```markdown
# Research notebook

## Notebook-commit gate

In the N-limitation v2 analysis, Steps 1-3 had rich notebook
entries because they were written as part of the step. Steps 4-6
lost their entries to chat history — the analysis gained momentum
and the notebook fell behind. The richest findings (proteomics
coverage bias, RNA/protein discordance test, rank normalization
fix) were discovered in those later steps and exist only in the
conversation. A future reader has methods.md but not the
investigative trail.

The fix: if a step produces data or analytical output, its
notebook entry must be committed to git before the next step
begins. The commit includes the notebook entry and that step's
generated artifacts (see [Artifacts guide — Git tracking
convention](artifacts.md)). This keeps the notebook as a live
record rather than a retroactive summary.
```

New:
```markdown
# Research notebook

## Step protocol and enforcement

For per-step commit timing, hard gates, and the chat-capture
pattern, see [Step protocol](step-protocol.md). This document
owns the notebook **format and content**; step-protocol owns
**when things happen and what gates enforce them**.
```

- [ ] **Step 2: Replace the step cycle section (lines 29-50)**

The step cycle definition stays but becomes a brief summary pointing to step-protocol for enforcement:

Old (lines 29-50):
```markdown
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
```

New:
```markdown
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
```

- [ ] **Step 3: Add chat-capture section to the entry template (after line 119)**

Insert a new section between `### Exploration` and `### Decision` in the entry template:

After the existing Exploration section:
```markdown
### Exploration
- Walked through gene X (PMM0xxx): values, logic, conclusion
- Checked gene Y: values — expected / unexpected because [reason]
- Asked: [question] → [answer with data]
```

Insert:
```markdown

### Chat exploration
**Q: [researcher's question, as asked]**
Data: [what was looked up / computed to answer it]
Finding: [what the data showed, with concrete numbers]
Impact: [how this affects interpretation or next steps]

**Q: [next question]**
...
```

- [ ] **Step 4: Remove the "Manifest updates" section (lines 213-217)**

This content moved to step-protocol.md. Remove:

```markdown
## Manifest updates

Every step that produces output files updates the relevant
manifest (`data/DATA_MANIFEST.md` or `results/RESULTS_MANIFEST.md`)
immediately — not retroactively at the end.
```

- [ ] **Step 5: Verify research-notebook.md**

Read the full file. Confirm:
- Opens with pointer to step-protocol.md for enforcement
- Step cycle is summarized, not duplicated
- Chat-capture format present in entry template
- Manifest updates section removed
- QC checkpoint types, code lifecycle, interactive discovery, rerun workflow all preserved
- No orphaned references to removed content

- [ ] **Step 6: Commit**

```bash
git add .claude/skills/research-methodology/references/step-protocol.md
git add .claude/skills/research-methodology/references/research-notebook.md
git commit -m "feat(skill): add step-protocol reference, update research-notebook

New step-protocol.md consolidates per-step commit timing, hard
gates (step boundary, manifest currency, chat-capture), redo
path, and git discipline into one reference.

research-notebook.md updated to point to step-protocol for
enforcement, adds chat-capture format to entry template, removes
manifest-timing section (moved to step-protocol)."
```

---

### Task 3: Update `references/artifacts.md` [Spec §3, §4]

Remove per-step commit rules and manifest-timing rules (moved to step-protocol.md). Replace blanket gitignore guidance with per-analysis .gitignore. Keep directory structure, file naming, manifest format.

**Files:**
- Modify: `.claude/skills/research-methodology/references/artifacts.md`

- [ ] **Step 1: Replace "Git tracking convention" section (lines 189-217)**

Old (lines 189-217):
```markdown
## Git tracking convention

Commit generated artifacts with the step that produces them, not
retroactively at the end. Each step's commit includes its outputs,
its log, and updated manifests.

**What to track vs gitignore:**

| Track in git | Gitignore |
|---|---|
| Signatures, scores, applied subsets (small CSVs) | Raw DE extracts (large CSVs, reproducible from KG) |
| Plots (PNG) | `__pycache__/` directories |
| Logs | |
| Manifests | |
| Notebook entries | |

Rule of thumb: if it's small and captures analytical decisions
(a signature gene list, a score table), track it. If it's large
and reproducible by rerunning a script against the KG, gitignore
it.

**Reproduction instructions:** When large files are gitignored,
the README must include instructions for regenerating them (which
scripts to run, in what order). A reader who clones the repo
should be able to reproduce the full analysis.

**Manifest timing:** The manifest (`DATA_MANIFEST.md` or
`RESULTS_MANIFEST.md`) is updated in the same commit that adds
the new data or result file. Not retroactively.
```

New:
```markdown
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
```

- [ ] **Step 2: Verify artifacts.md**

Read the full file. Confirm:
- Git tracking section points to step-protocol for commit timing
- Per-analysis .gitignore replaces blanket guidance
- Track-vs-gitignore table preserved
- Reproduction instructions preserved
- Manifest timing removed (now in step-protocol)
- "Commit generated artifacts with the step" wording removed (now in step-protocol)
- Directory structure, file naming, manifest format all preserved

- [ ] **Step 3: Commit**

```bash
git add .claude/skills/research-methodology/references/artifacts.md
git commit -m "refactor(skill): move commit timing to step-protocol, add per-analysis gitignore

artifacts.md now points to step-protocol.md for commit timing
and manifest-commit coupling. Replaces blanket gitignore guidance
with per-analysis .gitignore created at scaffolding."
```

---

### Task 4: Update SKILL.md [Spec §4]

Update Rule 8 pointer and add step-protocol to the references list.

**Files:**
- Modify: `.claude/skills/research-methodology/SKILL.md`

- [ ] **Step 1: Update Rule 8 (lines 95-104)**

Old:
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

New:
```markdown
## Rule 8: Research notebook, not pipeline

Every analysis is driven by an interactive notebook where each
step is recorded, explored with the researcher, and approved
before the next step. Implementation can be fast or delegated;
quality control is always interactive.

See [Step protocol](references/step-protocol.md) for per-step
commit timing, hard gates, the chat-capture pattern, and redo
path. See [Research notebook](references/research-notebook.md)
for notebook format, entry templates, QC checkpoint types, and
code lifecycle rules.
```

- [ ] **Step 2: Add step-protocol to references list (after line 114)**

Old:
```markdown
## References — read on demand, not all at once

- [KG rules](references/kg-rules.md) — read when scoping a new analysis or uncertain about data sourcing vs literature
- [Gene identity](references/gene-identity.md) — read when working with gene names, paralogs, or orthologs
- [Artifacts guide](references/artifacts.md) — read when setting up an analysis directory or unsure about file structure
- [Anti-hallucination](references/anti-hallucination.md) — read before presenting findings or when self-checking results
- [Python API guide](references/python-api-guide.md) — read before writing extraction or analysis scripts
- [Statistical rigor](references/statistical-rigor.md) — read when computing enrichment, comparing across studies, or reporting p-values
- [Research notebook](references/research-notebook.md) — read when starting or resuming an interactive analysis
```

New:
```markdown
## References — read on demand, not all at once

- [Step protocol](references/step-protocol.md) — read at the start of every analysis execution; owns commit timing, hard gates, chat-capture, redo path
- [KG rules](references/kg-rules.md) — read when scoping a new analysis or uncertain about data sourcing vs literature
- [Gene identity](references/gene-identity.md) — read when working with gene names, paralogs, or orthologs
- [Artifacts guide](references/artifacts.md) — read when setting up an analysis directory or unsure about file structure
- [Anti-hallucination](references/anti-hallucination.md) — read before presenting findings or when self-checking results
- [Python API guide](references/python-api-guide.md) — read before writing extraction or analysis scripts
- [Statistical rigor](references/statistical-rigor.md) — read when computing enrichment, comparing across studies, or reporting p-values
- [Research notebook](references/research-notebook.md) — read when starting or resuming an interactive analysis
```

- [ ] **Step 3: Verify SKILL.md**

Read the full file. Confirm:
- Rule 8 points to both step-protocol and research-notebook with clear ownership split
- Step-protocol is first in the references list with "read at the start of every analysis execution"
- Rules 1-7 unchanged
- No orphaned references

- [ ] **Step 4: Commit**

```bash
git add .claude/skills/research-methodology/SKILL.md
git commit -m "feat(skill): add step-protocol to Rule 8 and references list

Rule 8 now points to step-protocol.md for enforcement and
research-notebook.md for format. Step-protocol listed first
in references with guidance to read at every analysis start."
```

---

### Task 5: Update repo-level .gitignore [Spec §3]

Remove the blanket `analyses/*/data/` and `analyses/*/results/` rules and per-analysis exceptions. Each analysis will manage its own .gitignore going forward.

**Files:**
- Modify: `.gitignore`

- [ ] **Step 1: Replace analyses section (lines 12-22)**

Old:
```
# Analysis outputs (per-analysis, large CSVs not tracked by default)
analyses/*/data/
analyses/*/results/
# v2 analysis: track data/results (small, artifact-focused)
!analyses/2026-04-08-1038-n_limitation_signature_v2/data/
!analyses/2026-04-08-1038-n_limitation_signature_v2/results/
# v2: skip large DE extracts (reproducible from scripts 02 + 04)
analyses/2026-04-08-1038-n_limitation_signature_v2/data/de_*.csv
# B1 pathway enrichment: track data and results (all small files)
!analyses/2026-04-09-1713-pathway_enrichment_b1/data/
!analyses/2026-04-09-1713-pathway_enrichment_b1/results/
```

New:
```
# Analysis outputs: each analysis manages its own .gitignore
# (see step-protocol.md — per-analysis .gitignore created at scaffolding)
```

- [ ] **Step 2: Add per-analysis .gitignore to existing analyses**

Create `.gitignore` in each existing analysis directory so their tracking behavior is preserved.

For v2 (`analyses/2026-04-08-1038-n_limitation_signature_v2/.gitignore`):
```
# Large DE extracts (reproducible from scripts 02 + 04)
data/de_*.csv
__pycache__/
```

For B1 (`analyses/2026-04-09-1713-pathway_enrichment_b1/.gitignore`):
```
# All data and results are small files — track everything
__pycache__/
```

- [ ] **Step 3: Verify git status**

Run `git status` and confirm:
- No previously tracked files became untracked
- No previously ignored files suddenly appear as untracked
- Both per-analysis .gitignore files show as new

- [ ] **Step 4: Commit**

```bash
git add .gitignore
git add analyses/2026-04-08-1038-n_limitation_signature_v2/.gitignore
git add analyses/2026-04-09-1713-pathway_enrichment_b1/.gitignore
git commit -m "refactor: replace blanket analysis gitignore with per-analysis rules

Each analysis now manages its own .gitignore (created at
scaffolding per step-protocol.md). Removes repo-level blanket
rules and per-analysis exceptions."
```

---

### Task 6: Final verification

- [ ] **Step 1: Read all modified files end-to-end**

Read in order:
1. `.claude/skills/research-methodology/SKILL.md`
2. `.claude/skills/research-methodology/references/step-protocol.md`
3. `.claude/skills/research-methodology/references/research-notebook.md`
4. `.claude/skills/research-methodology/references/artifacts.md`
5. `.gitignore`

Confirm:
- No duplicated enforcement content between step-protocol and other files
- All cross-references resolve (step-protocol → research-notebook, artifacts; research-notebook → step-protocol; artifacts → step-protocol)
- Gate numbering consistent (1, 2, 3)
- Chat-capture format identical in step-protocol and research-notebook entry template

- [ ] **Step 2: Verify git log**

Run `git log --oneline -5` and confirm 4 clean commits from tasks 2-5.
