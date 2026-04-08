# Research Methodology Skill Improvements — Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Incorporate lessons from the N-limitation v2 analysis into the research-methodology skill, CLAUDE.md, and reference docs.

**Architecture:** Four file edits — artifacts guide (directory tree, git tracking, log verbosity, superpowers capture), research notebook (commit gate, interactive steps, rerun workflow), SKILL.md (Rule 5 amendment), CLAUDE.md (plan gates, numbering convention). All changes are documentation; no code or tests.

**Spec:** `docs/superpowers/specs/2026-04-08-research-methodology-v2-improvements-design.md`

---

## File Structure

All modifications to existing files:

```
skills/research-methodology/
├── SKILL.md                          # Modify: Rule 5 amendment (line ~65)
└── references/
    ├── artifacts.md                  # Modify: directory tree, naming, 4 new sections
    └── research-notebook.md          # Modify: commit gate, interactive steps, rerun workflow
CLAUDE.md                             # Modify: 3 additions to Research methodology section
```

---

### Task 1: Update artifacts guide — directory tree and naming

**Files:**
- Modify: `skills/research-methodology/references/artifacts.md:36-55`

- [ ] **Step 1: Replace the directory tree and naming paragraph**

Replace the current directory tree (lines 37-50) and naming paragraph (lines 52-54) with:

```markdown
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
```

- [ ] **Step 2: Verify the edit**

Read the file and confirm the tree, naming, and new directory explanations are in place. Confirm surrounding sections (Research notebook, gaps_and_friction.md) are intact.

---

### Task 2: Add superpowers artifact capture section

**Files:**
- Modify: `skills/research-methodology/references/artifacts.md` — insert after the directory tree section, before "Research notebook"

- [ ] **Step 1: Insert the superpowers artifact capture section**

Insert after the new directory explanations, before the `## Research notebook` heading:

```markdown
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
```

- [ ] **Step 2: Verify the edit**

Read the file and confirm the new section appears between the directory tree and the Research notebook section.

---

### Task 3: Add git tracking convention and log verbosity sections

**Files:**
- Modify: `skills/research-methodology/references/artifacts.md` — insert before `## Required analysis documents`

- [ ] **Step 1: Insert the git tracking convention section**

Insert before the `## Required analysis documents` heading:

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
```

- [ ] **Step 2: Verify the edit**

Read the file and confirm both new sections appear before the Required analysis documents section.

---

### Task 4: Update research notebook — commit gate and interactive steps

**Files:**
- Modify: `skills/research-methodology/references/research-notebook.md:6-16` (step cycle section) and after line 16 (insert new subsection)

- [ ] **Step 1: Add the notebook-commit gate**

Insert at line 6, before "Every analysis is driven by...", as the opening of the file body:

```markdown
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

---

```

- [ ] **Step 2: Add the interactive discovery steps subsection**

Insert after the "What a 'step' is" section (after line ~30):

```markdown
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
```

- [ ] **Step 3: Add log cross-reference**

Insert one line after the step entry template (after the `~~~` closing the template):

```markdown
Logs must be sufficient to verify the step without rerunning —
see [Artifacts guide — Log verbosity](artifacts.md) for the
standard.
```

- [ ] **Step 4: Verify the edits**

Read the full file and confirm the gate, interactive steps section, and log cross-reference are in place and the existing content is intact.

---

### Task 5: Expand rerun/revision workflow in research notebook

**Files:**
- Modify: `skills/research-methodology/references/research-notebook.md` — replace the "Rerun entries" section (~lines 96-120)

- [ ] **Step 1: Replace the rerun section**

Replace the current "Rerun entries" heading and content with:

```markdown
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
5. **Commit each rerun's outputs** — the commit-per-step rule
   still applies; don't batch rerun commits

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
```

- [ ] **Step 2: Verify the edit**

Read the file and confirm the rerun section is replaced cleanly, with the original rerun template gone.

---

### Task 6: Amend Rule 5 in SKILL.md

**Files:**
- Modify: `skills/research-methodology/SKILL.md:64-69` (after Rule 5 text)

- [ ] **Step 1: Add the interactive steps exception**

Insert after the existing Rule 5 content ("See [Python API guide]..." line), before Rule 6:

```markdown

**Exception: interactive discovery steps.** Steps that are
inherently exploratory (browsing available data, classifying
experiments, scoping what the KG contains) may be done
interactively rather than scripted. These steps must still produce
a frozen output file (CSV) and a notebook entry documenting the
reasoning. Computations — statistics, scores, enrichment — always
go in scripts. See [Research notebook — Interactive discovery
steps](references/research-notebook.md) for the pattern.
```

- [ ] **Step 2: Verify the edit**

Read SKILL.md and confirm Rule 5 now has the exception paragraph and Rule 6 follows cleanly.

---

### Task 7: Update CLAUDE.md

**Files:**
- Modify: `CLAUDE.md:72-76` (Process overrides section) and `CLAUDE.md:57-71` (Research methodology section)

- [ ] **Step 1: Add notebook-commit gates to Process overrides**

Add after the existing "Don't skip subagent reviews" bullet (line 75):

```markdown

- **Notebook-commit gates in plans.** When writing implementation
  plans for research analyses, every plan must include a
  notebook-commit gate between data-producing steps. The gate is:
  "commit notebook entry for Step N before beginning Step N+1."
  This ensures the executing agent treats the notebook as a
  blocking dependency, not a nice-to-have.
```

- [ ] **Step 2: Add spec steps vs plan tasks**

Add after the "Use the superpowers workflow..." paragraph (line 70), before "### Process overrides":

```markdown

### Naming conventions

- **Spec steps ≠ plan tasks.** Spec steps describe the analytical
  pipeline (Step 1, Step 2, ...). Plan tasks describe
  implementation work and may split one spec step into multiple
  tasks. Specs must not reference files by task number. Plans must
  cross-reference the spec step in each task description (e.g.,
  "Task 5: Build signature utilities [Step 3]"). This prevents
  confusion when the plan introduces scaffolding or verification
  tasks that don't correspond to a spec step.

- **Interactive discovery steps.** Rule 5 (scripts over chat
  reasoning) has a carve-out for interactive discovery/scoping
  steps — see the research-methodology skill for the
  frozen-output + notebook-entry pattern.
```

- [ ] **Step 3: Verify the edits**

Read CLAUDE.md and confirm both additions are in the right sections and the file reads coherently.

---

### Task 8: Commit

- [ ] **Step 1: Commit all changes**

```bash
git add \
  skills/research-methodology/SKILL.md \
  skills/research-methodology/references/artifacts.md \
  skills/research-methodology/references/research-notebook.md \
  CLAUDE.md \
  docs/superpowers/specs/2026-04-08-research-methodology-v2-improvements-design.md \
  docs/superpowers/plans/2026-04-08-research-methodology-v2-improvements.md

git commit -m "docs: improve research methodology skill from N-lim v2 retrospective

Add notebook-commit gate, interactive discovery step pattern,
git tracking convention, log verbosity standard, superpowers
artifact capture, and expanded rerun workflow to the research
methodology skill. Update CLAUDE.md with plan gates and
spec-vs-plan numbering convention."
```
