# Research Methodology Skill Improvements (post N-lim v2)

**Date:** 2026-04-08
**Status:** Active spec
**Source:** Retrospective from `analyses/2026-04-08-1038-n_limitation_signature_v2/gaps_and_friction.md`

## Motivation

The N-limitation signature v2 analysis applied the research-methodology skill from the start and produced the best-quality output so far. The retrospective identified concrete gaps where the skill didn't guide well enough, leading to drift under momentum. This spec describes targeted improvements to close those gaps.

## Scope

Four files edited, no new files:

| File | What changes |
|------|-------------|
| `skills/research-methodology/references/artifacts.md` | Directory tree, superpowers capture, git tracking, log verbosity, manifest timing |
| `skills/research-methodology/references/research-notebook.md` | Notebook-commit gate, interactive discovery steps, log cross-reference, expanded rerun workflow |
| `skills/research-methodology/SKILL.md` | Rule 5 amendment for interactive steps |
| `CLAUDE.md` | Plan notebook gates, spec-vs-plan numbering, Rule 5 cross-reference |

## Out of scope

- MCP/KG improvements (matching_genes in list_experiments, resolve_gene parameter naming)
- Plotting cosmetics
- Superpowers skill changes (external, not controlled)
- New rules or references — all changes are additions to existing sections

---

## Change 1: Artifacts guide — directory tree and conventions

**File:** `skills/research-methodology/references/artifacts.md`

### 1a. Updated directory tree

Replace the current directory tree with:

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

Replace the existing naming paragraph ("Naming: `{analysis_name}` should be descriptive...") with:

> **Naming:** `YYYY-MM-DD-HHMM-{name}` — timestamp prefix for chronological ordering. `{name}` is descriptive and lowercase with underscores (e.g., `2026-04-08-1038-n_limitation_signature_v2`).

New entries explained:
- **`YYYY-MM-DD-HHMM-{name}`** — as above.
- **`{name}_utils/`** — reusable methodology package when the analysis introduces novel methods (signature scoring, enrichment, etc.). Contains pure methodology code with tests. Only created when the analysis needs it — not every analysis gets one. See research-notebook code lifecycle for when this applies.
- **`logs/`** — one log per pipeline step, capturing diagnostics sufficient to verify the step without rerunning.
- **`superpowers/`** — local copies of the spec, plan, and brainstorm-log. See "Superpowers artifact capture" below.

### 1b. Superpowers artifact capture

New section after the directory tree:

> **Superpowers artifact capture**
>
> After brainstorming and planning, copy the spec, plan, and brainstorm-log into the analysis's `superpowers/` directory. The analysis should be self-contained — a reader should find the full decision history (why this design, what alternatives were considered, what the plan was) alongside the data and code, not scattered across `docs/superpowers/`.
>
> The canonical copies in `docs/superpowers/specs/` and `docs/superpowers/plans/` remain the repo-level index. The copies in the analysis are the local record.
>
> Files to capture:
> - `superpowers/spec.md` — the design spec
> - `superpowers/plan.md` — the implementation plan
> - `superpowers/brainstorm-log.md` — the Q&A from brainstorming (questions, options, decisions, rationale)

### 1c. Git tracking convention for analysis outputs

New section:

> **Git tracking convention**
>
> Commit generated artifacts with the step that produces them, not retroactively at the end. Each step's commit includes its outputs, its log, and updated manifests.
>
> **What to track vs gitignore:**
>
> | Track in git | Gitignore |
> |---|---|
> | Signatures, scores, applied subsets (small CSVs) | Raw DE extracts (large CSVs, reproducible from KG) |
> | Plots (PNG) | `__pycache__/` directories |
> | Logs | |
> | Manifests | |
> | Notebook entries | |
>
> Rule of thumb: if it's small and captures analytical decisions (signature gene list, score table), track it. If it's large and reproducible by rerunning a script against the KG, gitignore it.
>
> **Reproduction instructions:** When large files are gitignored, the README must include instructions for regenerating them (which scripts to run, in what order). A reader who clones the repo should be able to reproduce the full analysis.
>
> **Manifest timing:** The manifest (`DATA_MANIFEST.md` or `RESULTS_MANIFEST.md`) is updated in the same commit that adds the new data or result file. Not retroactively.

### 1d. Log verbosity standard

New section (or addition to existing logs guidance):

> **Log verbosity**
>
> Logs capture everything a researcher needs to verify the step without rerunning it. This includes:
> - Summary statistics (row counts, gene counts, filter funnel)
> - Diagnostic traces (marker gene values at each stage, QC checks)
> - Edge case results (tie-breaking outcomes, genes at classification boundaries)
>
> Scripts with `--explore` flags should write their diagnostic output to the log file as well as stdout. The log is the persistent record; stdout is for the interactive session.
>
> A log that says "N genes built" is insufficient. A log that shows the filter funnel (how many started → how many passed each filter → how many in the final set), marker gene traces, and classification edge cases is sufficient.

---

## Change 2: Research notebook — gates and patterns

**File:** `skills/research-methodology/references/research-notebook.md`

### 2a. Notebook-commit gate

Add at the top of "The step cycle" section, before the do→show→explore→decide description:

> **Notebook entry committed before next step begins.**
>
> In the N-limitation v2 analysis, Steps 1-3 had rich notebook entries because they were written as part of the step. Steps 4-6 lost their entries to chat history — the analysis gained momentum and the notebook fell behind. The richest findings (proteomics coverage bias, RNA/protein discordance test, rank normalization fix) were discovered in those later steps and exist only in the conversation. A future reader has methods.md but not the investigative trail.
>
> The fix: if a step produces data or analytical output, its notebook entry must be committed to git before the next step begins. The commit includes the notebook entry and that step's generated artifacts (see artifacts guide — git tracking convention). This keeps the notebook as a live record rather than a retroactive summary.

### 2b. Interactive discovery steps

New subsection after "What a 'step' is":

> **Interactive discovery steps**
>
> Some steps are naturally exploratory — browsing experiments, classifying conditions, discussing what to include. These don't fit neatly into "run a script, check the output." The pattern for interactive steps:
>
> 1. Explore interactively via MCP queries and chat discussion
> 2. Produce a **frozen output file** (CSV, table) that downstream scripts consume — this is the reproducible artifact
> 3. Write a **notebook entry** documenting the reasoning, classifications, and decisions — this must be detailed enough that the reasoning could be reproduced by a different researcher
> 4. Optionally, write a script that reproduces the frozen output from KG queries
>
> The frozen output + notebook entry are the minimum. A script is preferred but not required for discovery steps where the process is inherently iterative and conversational.
>
> Rule 5 (scripts over chat reasoning) still applies to computations. Interactive steps are for *discovery and classification*, not for computing statistics in chat.

### 2c. Log cross-reference

Add one line near the step entry template:

> Logs must be sufficient to verify the step without rerunning — see artifacts guide for the verbosity standard.

### 2d. Expanded rerun/revision workflow

Replace or expand the current "Rerun entries" section:

> **Rerun and revision workflow**
>
> When a methodology change is needed mid-analysis (formula correction, normalization change, filter adjustment):
>
> 1. **Update the utility code** — change the shared methodology (`*_utils/`), not just the script
> 2. **Run tests** — verify the change is correct against toy data
> 3. **Rerun all downstream steps** that consumed the changed output — every step after the change gets a fresh run
> 4. **New notebook entry per rerun** — each rerun gets its own entry with a "Why" section explaining what changed and what the impact was
> 5. **Commit each rerun's outputs** — the commit-per-step rule still applies; don't batch rerun commits
>
> The rerun entry template:
>
> ~~~markdown
> ---
>
> ## YYYY-MM-DD HH:MM — Step N rerun: {what changed}
>
> ### Why
> What was wrong with the previous output and how it was discovered.
>
> ### Changes
> - What code was modified (file, function, logic change)
> - How tests were affected (new tests, changed expected values)
> - How many downstream steps needed rerunning
>
> ### Outputs
> - Updated files with new row counts / statistics
>
> ### Impact
> - What changed in the results (direction, magnitude, interpretation)
> - What didn't change (confirming robustness where applicable)
>
> ### Decision
> Proceed with updated output / need another iteration.
> ~~~

---

## Change 3: SKILL.md — Rule 5 amendment

**File:** `skills/research-methodology/SKILL.md`

Add one paragraph to Rule 5:

> **Exception: interactive discovery steps.** Steps that are inherently exploratory (browsing available data, classifying experiments, scoping what the KG contains) may be done interactively rather than scripted. These steps must still produce a frozen output file (CSV) and a notebook entry documenting the reasoning. Computations — statistics, scores, enrichment — always go in scripts. See [Research notebook — Interactive discovery steps](references/research-notebook.md) for the pattern.

---

## Change 4: CLAUDE.md

**File:** `CLAUDE.md`

### 4a. Notebook gate in plans

Add to the "Process overrides" section:

> - **Notebook-commit gates in plans.** When writing implementation plans for research analyses, every plan must include a notebook-commit gate between data-producing steps. The gate is: "commit notebook entry for Step N before beginning Step N+1." This ensures the executing agent treats the notebook as a blocking dependency, not a nice-to-have.

### 4b. Spec steps vs plan tasks

Add to the "Research methodology" section:

> - **Spec steps ≠ plan tasks.** Spec steps describe the analytical pipeline (Step 1, Step 2, ...). Plan tasks describe implementation work and may split one spec step into multiple tasks. Specs must not reference files by task number. Plans must cross-reference the spec step in each task description (e.g., "Task 5: Build signature utilities [Step 3]"). This prevents confusion when the plan introduces scaffolding or verification tasks that don't correspond to a spec step.

### 4c. Interactive discovery steps cross-reference

Add to the "Research methodology" section:

> - **Interactive discovery steps.** Rule 5 (scripts over chat reasoning) has a carve-out for interactive discovery/scoping steps — see the research-methodology skill for the frozen-output + notebook-entry pattern.

---

## Evaluation

These are incremental improvements to an existing reference skill, not a new skill. Formal test cases (skill-creator style) are not warranted — the changes are documentation, not behavioral logic.

The next research analysis using this skill is the real evaluation. Specifically, watch for:
- Does the generated plan include notebook-commit gates between steps?
- Does the notebook stay current through the later steps?
- Does the first step use the interactive discovery pattern (frozen CSV + notebook entry) without trying to force a script?
- Are artifacts committed per step, not batched at the end?
- Do logs contain enough detail to verify steps without rerunning?

If any of these fail, the skill language needs to be stronger or restructured — revisit with another retrospective.
