# Step protocol

The step cycle (do → show → explore → decide) is defined in
[research-notebook.md](research-notebook.md). This document owns
**when things happen and what gates enforce them**. Follow it for
every step that produces data or analytical output.

## Commit structure

Two commits per step:
- **Commit 1** (end of "do"): script + outputs + log + manifest
  updates
- **Commit 2** (end of "decide"): notebook entry with QC tables,
  chat-capture section, and decision

Show, explore, and decide produce ONE notebook entry together —
QC section from show, chat-capture section from explore, decision
from decide — committed once at the end of decide.

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
the notebook entry, closing it. Then git commit (Commit 2): the
completed notebook entry (QC + chat-capture + decision).

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
