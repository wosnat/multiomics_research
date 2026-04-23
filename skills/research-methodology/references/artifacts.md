# Artifacts guide

Research questions produce files, not chat messages. Chat is for reasoning and dialogue; data, statistics, figures, and step narratives go to disk.

## Contents

1. [When to create an analysis directory](#when-to-create-an-analysis-directory)
2. [Directory structure per analysis](#directory-structure-per-analysis)
3. [Scaffold creation](#scaffold-creation)
4. [Per-step folders](#per-step-folders)
5. [paper.md and gaps_and_friction.md](#papermd-and-gaps_and_frictionmd)
6. [Data and figure conventions](#data-and-figure-conventions)
7. [Git tracking](#git-tracking)
8. [Log verbosity](#log-verbosity)
9. [References](#references)

## When to create an analysis directory

Create when:
- The question requires statistical computation
- Multiple data extractions feed into an analysis
- Results need to be shared or revisited
- The work takes more than 2–3 tool calls
- Figures or publication-ready tables are needed

Simple lookups stay in chat:
- "What is gene PMM0120?" — MCP lookup, answer in chat
- "How many genes respond to nitrogen?" — MCP summary, answer in chat
- "List the catalase genes" — MCP query, table in chat

## Directory structure per analysis

```
analyses/YYYY-MM-DD-<slug>/
  paper.md                  # academic-style writeup, grows from step 1
  gaps_and_friction.md      # transitional — methodology/KG/tooling friction log
  .gitignore                # per-analysis, created at scaffolding
  1_question/
    notebook.md
  2_kg_selection/
    notebook.md, scripts/, data/, figures/
  3_analysis_framing/
    notebook.md, scripts/, data/, figures/
  4_methods/
    notebook.md
    <module_name>.py        # ad-hoc methods module
    scripts/, data/, figures/
  5_analyze/
    notebook.md, scripts/, data/, figures/
  6_evaluate/
    notebook.md, scripts/, data/, figures/
```

**Slug format:** `YYYY-MM-DD-<descriptor>` where descriptor is snake_case, 2–4 words reflecting the research question (e.g., `2026-04-23-n_limitation_signature`). Same-day collisions: append a letter suffix (`-a`, `-b`).

**Why per-step folders:** keeps each step's artifacts self-contained and reviewable; avoids an analysis-wide `scripts/` or `results/` that becomes a dumping ground and loses step context.

## Scaffold creation

**Who:** Claude proposes the slug; the user approves or counter-proposes.
**When:** at the start of step 1, before any folder or file is written.

1. **Propose a slug** based on the user's initial prompt. Confirm with the user before writing anything — no scaffold exists until the slug is confirmed.
2. **Create the minimal scaffold:**
   ```
   analyses/YYYY-MM-DD-<slug>/
     paper.md              # skeleton with empty section headers: Question, Background, Methods, Results, Discussion, References
     gaps_and_friction.md  # header only; entries added as friction occurs
     .gitignore            # per-analysis (template below)
     1_question/
       notebook.md         # empty; populated by the step 1 dialogue
   ```
3. **Begin step 1 dialogue.** `1_question/notebook.md` grows through the clarifying exchange.
4. **At step 1 decide:** `paper.md`'s Question section is populated; the first git commit includes the scaffold + step 1 artifacts.

Step folders for 2–6 are created **progressively** — each folder (and its `scripts/` / `data/` / `figures/` subfolders as needed) is created when that step begins, not pre-created at scaffold time. This avoids empty folders cluttering the analysis during work in progress.

## Per-step folders

Each step folder except `1_question/` contains:
- **`notebook.md`** — narrative + decide-gate checklist (see [research-notebook.md](research-notebook.md))
- **`scripts/`** — Python scripts that do the step's work
- **`data/`** — script outputs (CSV/TSV); inputs from prior steps referenced by relative path
- **`figures/`** — PNG/PDF/SVG outputs

Step 1's folder contains only `notebook.md` — step 1 is a conversation, not a computation.

### QC artifacts

QC artifacts live in the same `scripts/` / `data/` / `figures/` subfolders with a **`qc_` filename prefix** (e.g., `scripts/qc_check_tp_coverage.py`, `data/qc_tp_coverage.csv`, `figures/qc_tp_coverage_histogram.png`). No separate `qc/` subfolder. The prefix sorts QC files together in listings without creating a binary main/QC split at the folder level.

### Methods module (step 4 only)

`4_methods/<module_name>.py` holds the ad-hoc method module (scoring function, enrichment utility, etc.). Subsequent steps import from it by relative path. See [research-notebook.md — Code lifecycle](research-notebook.md) for the methodology-first → productization-later pattern.

### Script requirements

Scripts must be:
- Self-contained
- Documented at the top (purpose, inputs, outputs)
- Deterministic (seeded if random)
- Have clear I/O — no hardcoded absolute paths; arguments via CLI

### Script naming

- **Pipeline scripts (numbered within the step):** `01_select_experiments.py`, `02_compute_scores.py`
- **QC scripts:** `qc_<what_is_checked>.py`
- **Ad-hoc exploration scripts:** `explore_<what>.py` — kept in `scripts/` for reproducibility even if one-off

## paper.md and gaps_and_friction.md

Both live at the analysis root. Content rules are in [research-notebook.md](research-notebook.md):
- `paper.md` — academic-style writeup, skeleton sections from day 1, populated at each step's decide phase
- `gaps_and_friction.md` — transitional friction log, distinct from decisions (which live in step notebooks)

## Data and figure conventions

### data/
- Format: CSV (preferred) or TSV
- One file per extraction query or script output
- Include metadata columns (organism, experiment_id, locus_tag) when relevant
- Name descriptively: `de_genes_med4_nitrogen.csv`, not `output.csv`

### figures/
- PNG at 300 DPI minimum for notebook inclusion
- PDF or SVG for publication-ready figures
- Axis labels, legends, and titles mandatory

## Git tracking

For per-step commit timing, hard gates, and redo commits, see [step-protocol.md](step-protocol.md).

### Per-analysis .gitignore

Created during scaffolding (at start of step 1). Default template:
```
# Large intermediate data reproducible from KG
# (list specific files here, not blanket patterns)
__pycache__/
```

No blanket repo-level rules for `analyses/*/data/` or `analyses/*/figures/`. Tracking decisions are explicit, per-file, and logged in the step's `notebook.md`.

### What to track vs gitignore

| Track in git | Gitignore (list explicitly) |
|---|---|
| Signatures, scores, applied subsets (small CSVs) | Raw DE extracts (large CSVs, reproducible from KG) |
| Figures (PNG/PDF/SVG) | `__pycache__/` |
| `notebook.md`, `paper.md`, `gaps_and_friction.md` | |
| Scripts (all) | |

Rule of thumb: if it's small and captures analytical decisions, track it. If it's large and reproducible by rerunning a script, gitignore it explicitly and document the regeneration command in the step's `notebook.md` or in `paper.md` Methods.

## Log verbosity

Each step's scripts should emit enough diagnostic output to verify the step without rerunning. Two capture paths:

- **Short logs:** paste into the step's `notebook.md` (Results or Surprises sections)
- **Long logs:** per-script log file at `<step>/data/<script_name>.log`, referenced from `notebook.md`

**Minimum content:** summary statistics (row counts, gene counts, filter funnel), diagnostic traces (marker gene values at key stages), edge-case results (tie-breaking outcomes, classification boundary cases).

A log that says "N genes built" is insufficient. A log that shows the filter funnel (how many started → passed each filter → final count), marker gene traces, and classification edge cases is sufficient.

Scripts with `--explore` or `--verbose` flags should write their diagnostic output to the log file as well as stdout. The log is the persistent record; stdout is for the interactive session.

## References

References (publications cited in the analysis) live in **`paper.md`'s References section**, not a separate file. Every reference must be resolved through `list_publications` and cited by DOI or KG experiment ID — never from intrinsic knowledge (see [anti-hallucination.md — Category 5.2](anti-hallucination.md)).

KG version, tool versions, statistical software, and Python package versions go in `paper.md`'s Methods section.
