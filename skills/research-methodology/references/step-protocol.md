# Step protocol

Each analysis step advances through the rhythm **do → show → explore → decide** (phase content defined in [research-notebook.md](research-notebook.md)). This document owns **when things happen and what gates enforce them**.

The **just-in-time formalization** principle applies throughout: terms, metrics, stability checks, decisions, and caveats enter the analysis only when the data demands them. Do not pre-specify framework inventories. See [research-notebook.md](research-notebook.md) for the full principle and its application to step 3 framing.

## Commit structure

**One commit per step**, at the end of the decide phase. The commit includes everything the step produced:

- `notebook.md` (narrative + decide-gate checklist)
- `scripts/`, `data/`, `figures/` (main + `qc_*` files)
- updates to `paper.md` (the step's synthesis section)
- updates to `gaps_and_friction.md` (if the step encountered friction)

No mid-step commits, no separate "do commit + decide commit" pattern. The decide phase is the atomic step boundary.

Step 1 is a special case: the commit includes both the scaffold (analysis folder, `paper.md` skeleton, `gaps_and_friction.md` header, `1_question/notebook.md`) and step 1's own outputs. See [artifacts.md](artifacts.md) for scaffold creation.

## Before starting a step

- Previous step's commit exists
- Previous step's `notebook.md` has the decide-gate checklist populated
- Previous step's `paper.md` section is populated

If any of these is missing, close the previous step first.

## "do" phase

Do the step's work — scope depends on the step:

- **Step 1:** clarifying dialogue with the researcher (see [research-notebook.md — Using brainstorming for step 1](research-notebook.md))
- **Steps 2, 4, 5:** write and run scripts; produce data and figures
- **Steps 3, 6:** select, validate, or evaluate; produce scripts + data + figures for the QC side; write prose for the framing or evaluation side

Outputs land wherever the step naturally produces them: a conversation lands in `notebook.md`; scripts land in `scripts/`; their outputs land in `data/` and `figures/`. QC artifacts use the `qc_` filename prefix (see [artifacts.md](artifacts.md)).

No commit yet. The step's outputs are uncommitted working-tree state until decide.

## "show" phase

Populate `notebook.md` with what was produced. Recommended sections (see [research-notebook.md](research-notebook.md) for full content):

- **Context** — what this step is for; what the prior step decided
- **What I did** — scripts run with their command-line invocation for non-trivial cases; KG queries issued
- **Results** — summary tables shown inline (as markdown tables, not prose paraphrases); links to full tables in `data/` and figures in `figures/`; cited publications resolved via `list_publications` (never from memory — see [anti-hallucination.md — Category 5](anti-hallucination.md))

Summary tables in **Results** are the same tables presented to the researcher in chat — copied as markdown into the notebook, not paraphrased.

## "explore" phase

Investigate anomalies, surprises, or gaps:

- Ask follow-up clarifying questions (step 1)
- Add `qc_*.py` checks; run sensitivity analyses; cross-validate against controls (steps 2–6)

Capture anomalies worth flagging as **Surprises** in `notebook.md`. If a researcher question during this phase produces a data point or changes interpretation, both the prose and the data live in the notebook — the narrative IS the exploration record. No separate chat-capture section.

## "decide" phase

1. **Finalize `notebook.md`:**
   - Ensure Context / What I did / Results / Surprises are populated as applicable
   - Add **Decisions** section if any forks were taken (prose + date; see [research-notebook.md](research-notebook.md))
   - Write the **decide-gate checklist** at the end of notebook.md:
     - **Outputs produced** — filenames in `scripts/`, `data/`, `figures/`, with command lines for non-trivial scripts
     - **Results presented** — summary tables shown inline; links to full tables and figures generated this step
     - **QC gate** — what was checked → result (one line per check)
     - **Decisions made this step** — prose + date, if any; omit the section if none
     - **Advance rationale** — one line, why this step is ready to close

2. **Update `paper.md`:** write the synthesis paragraph (or figure inclusion, or methods sub-section) for the section this step populates — see [research-notebook.md — paper.md growth](research-notebook.md) for the section-to-step mapping.

3. **Append to `gaps_and_friction.md`** if friction was encountered this step (KG issues, MCP schema mismatches, methodology gaps, anti-hallucination corrections).

4. **Present state to researcher:** show the `notebook.md` content, the `paper.md` diff, and any `gaps_and_friction.md` additions. Wait for explicit approval or redirect.

5. **On approval, commit.** One commit, containing all of the step's changes.

6. Begin next step (create its folder as needed — see [artifacts.md](artifacts.md) for progressive folder creation).

## Redo path

When the researcher says "redo step N with X":

1. **do:** update script or framing; rerun; regenerate outputs. New artifacts overwrite old in the step folder.
2. **show / explore:** new tables, figures, Results; update Surprises if changed.
3. **decide:** new decide-gate checklist, new `paper.md` synthesis, new `gaps_and_friction.md` entry if the redo surfaced friction. **New commit (never amend the previous).**

The previous commit remains in git history as the record of what was tried. The working-tree `notebook.md` is overwritten because it now describes what actually happened in the successful attempt — it is not a log of prior attempts.

If the redo invalidates downstream steps, the redo's `notebook.md` must list the downstream steps that consumed its outputs. The researcher decides whether to cascade the redo.

`gaps_and_friction.md` is append-only: redo friction entries accumulate.

## Hard gates

### GATE 1: Step boundary

B1 and B2 partially wrote notebooks retroactively — exploration reasoning was lost and couldn't be verified against the actual data state at the time.

**Do not start step N+1 until step N is committed, including `notebook.md`, `paper.md` updates, and `gaps_and_friction.md` updates if applicable.**

### GATE 2: Researcher approval

B2's scope drift (decisions D5–D8 added mid-execution) slipped past because there was no atomic gate between "I finished some work" and "I'm advancing." The decide phase presents state to the researcher; the researcher approves, requests a redo, or redirects.

**Do not commit the step without explicit researcher approval of the decide-gate state.**

### GATE 3: Results presented, not paraphrased

Summary tables shown in chat must also appear as markdown tables in `notebook.md`. Prose paraphrases of numbers lose precision and are unreviewable.

**Do not close the step if the Results section in `notebook.md` is prose where a table belongs.**

## Git discipline

### Per-analysis .gitignore

Created during scaffolding (at start of step 1). Default:
```
# Large intermediate data reproducible from KG
# (list specific files here, not blanket patterns)
__pycache__/
```
Everything else tracked by default. Explicit entries with a comment explaining why.

### Scaffolding

Scaffolding and step 1 land in the same commit. Claude creates the scaffold during step 1's do phase, before the dialogue begins (see [artifacts.md — Scaffold creation](artifacts.md)). No separate scaffolding commit.

### Redo commits

Redo produces new commits, not amendments. The failed attempt's commit stays in history. Within the working tree, the step's `notebook.md` is overwritten to reflect the successful attempt.
