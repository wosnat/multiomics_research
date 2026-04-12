# Step Protocol Enforcement — Design Spec

## Goal

Harden the research-methodology skill's process enforcement by
consolidating per-step obligations into a single step protocol,
adding a chat-capture gate for interactive exploration, and fixing
git commit/ignore discipline.

**Scope:** Process enforcement only. No new analytical content
(ontology selection, enrichment patterns). KG/MCP tool changes out
of scope.

**Motivation:** B1 pathway enrichment analysis revealed recurring
violations despite existing rules:
- Notebook entries written retroactively (steps 1-2 explored, then
  bulk-written)
- Manifests updated after the fact, not with producing commit
- Valuable chat exploration lost — researcher's questions drove key
  findings but only fragments reached the notebook
- Results blanket-gitignored, requiring manual exceptions
- Per-step commits not enforced

These aren't new rules — they exist across research-notebook.md,
artifacts.md, and CLAUDE.md. The problem is that they're scattered
across 3 files, organized by artifact type rather than by when they
happen. The executing agent reconciles them imperfectly and drops
obligations when optimizing for throughput.

## Design

### 1. New reference: `references/step-protocol.md`

Single source of truth for what must happen at each phase of the
step cycle. Organized by **when**, not by **what**.

#### Step cycle

```
do      → write/update script → run → update manifests → COMMIT 1
show    → present QC to researcher (captured in notebook entry)
explore → interactive walkthrough (captured in notebook with chat-capture) → COMMIT 2
decide  → continue / redo / adjust
```

Two commits per step:
- **Commit 1** (end of "do"): script + outputs + log + manifest
  updates
- **Commit 2** (end of "explore"): notebook entry with chat-capture
  section

Show and explore produce ONE notebook entry together (QC section
from show, then chat-capture section from explore), committed once
at the end of explore.

Rationale for two commits: if the researcher requests a redo, the
exploration reasoning is preserved separately from the artifacts
that get replaced.

#### Per-phase obligations

**Before starting a step:**
- Confirm previous step's artifacts are committed (script +
  outputs + log + manifests — Commit 1)
- Confirm previous step's notebook entry is committed (including
  chat-capture — Commit 2)

**"do" phase:**
1. Write or update the script
2. Run script, capture outputs and log
3. Update DATA_MANIFEST.md or RESULTS_MANIFEST.md with new files
4. Git commit: script + outputs + log + manifest updates

**"show" phase** (begins the notebook entry for this step):
- Present QC diagnostics to researcher
- Write QC section of notebook entry, including:
  - Summary tables of outputs (same tables shown in chat, written
    to notebook — not prose paraphrases)
  - Links to figures produced (`![description](../results/fig.png)`)
  - Row counts, gene counts, key statistics

**"explore" phase** (continues the same notebook entry):
- Interactive walkthrough with researcher
- Append chat-capture section to notebook entry AS exploration
  happens, not after (see Section 2)
- GATE 3: chat-capture section must exist before moving to decide

**"decide" phase** (closes the notebook entry):
- Researcher says continue, redo, or adjust
- Decision logged in notebook entry
- Git commit (Commit 2): completed notebook entry

#### Redo path

When the researcher says "redo with X":

1. **do:** Update script, rerun, update manifests → NEW commit
   (never amend previous)
2. **show:** New QC presentation — same requirements as first pass
   (summary tables, figure links, counts)
3. **explore:** New walkthrough with chat-capture → APPEND new
   notebook entry (never revise previous entries)
4. **decide:** Researcher reviews. Protocol explicitly lists
   downstream steps that consumed changed output. Researcher
   decides whether to cascade redo or accept as-is. Decision
   logged in notebook.

The failed attempt's notebook entry stays as a record of what was
tried and why it was rejected. Append-only.

#### Hard gates

Each gate includes a "why" citing the real failure it prevents,
plus a stop condition.

**GATE 1: Step boundary**
> B1's notebook was written retroactively — exploration reasoning
> was lost and couldn't be verified against the actual data state
> at the time.
>
> Do not start step N+1 until step N's notebook entry (including
> chat-capture section) is committed.

**GATE 2: Manifest currency**
> B1's manifests were updated in bulk at the end. By that point,
> file descriptions were reconstructed from memory rather than
> written when the data was fresh.
>
> Do not commit a script's outputs without updating the relevant
> manifest in the same commit.

**GATE 3: Chat-capture**
> B1's most valuable findings came from researcher questions during
> explore ("is this dominated by catch-all categories?", "what
> about Steglich's low power?"). These drove analytical decisions
> but were lost from the chat context.
>
> Do not enter "decide" phase until the chat-capture section exists
> in the notebook entry.

### 2. Chat-capture pattern

New section in notebook entry format. Written during the "explore"
phase, not retroactively.

**Format:**

```markdown
## Chat exploration

**Q: [researcher's question, as asked]**
Data: [what was looked up / computed to answer it]
Finding: [what the data showed, with concrete numbers]
Impact: [how this affects interpretation or next steps]

**Q: [next question]**
...
```

**What gets captured:**
- Questions that surfaced data points, changed interpretation, or
  influenced decisions
- Any exchange that produced a number, a comparison, or a
  "therefore"

**What does NOT get captured:**
- Clarifications ("what does column X mean")
- Formatting/typo fixes
- Tool-call mechanics

**Why not verbatim chat?** Chat is noisy with tool calls, retries,
formatting. The value is the Q -> data -> finding -> impact chain.

### 3. Git discipline

**Per-analysis .gitignore:**

No blanket `analyses/*/results/` rule at the repo level. Each
analysis gets its own `.gitignore` created during scaffolding
(before step 1), with explicit decisions logged in the notebook.

Default template:
```
# Large intermediate data reproducible from KG
# (list specific files here, not blanket patterns)
__pycache__/
```

Everything else tracked by default. If a file should be ignored,
it's listed explicitly with a comment explaining why.

**Commit discipline:**
- Two commits per step (do + explore), as defined above
- Redo: new commits, never amend
- Each commit is self-contained: includes the artifacts it
  produces AND the manifest/notebook updates for those artifacts

**Scaffolding commit:**
When creating an analysis directory, the initial scaffolding
(directory structure, empty manifests, .gitignore, notebook stub)
gets its own commit before step 1 begins. This establishes the
tracking decisions upfront.

### 4. Changes to existing skill files

**SKILL.md:**
- Rule 8 pointer updated: "See step-protocol.md for the per-step
  protocol and enforcement gates. See research-notebook.md for
  notebook entry format and templates."
- Add step-protocol.md to the references list with guidance:
  "read at the start of every analysis execution"

**references/research-notebook.md:**
- Add chat-capture entry format (the Q/Data/Finding/Impact
  template)
- Add redo append-only entry template
- REMOVE step-cycle enforcement content (gates, commit timing) —
  that moves to step-protocol.md
- KEEP: notebook format, entry templates, QC checkpoint types,
  code lifecycle, interactive discovery pattern

**references/artifacts.md:**
- REMOVE per-step commit rules and manifest-timing rules (move to
  step-protocol.md)
- REPLACE blanket gitignore guidance with "per-analysis .gitignore
  created at scaffolding"
- KEEP: directory structure, file naming conventions, what-to-track
  guidance, manifest format

**Unchanged files:**
- Rules 1-7 in SKILL.md
- gene-identity.md
- kg-rules.md
- anti-hallucination.md
- statistical-rigor.md
- python-api-guide.md

## What success looks like

In the next analysis:
1. Every step has two commits — one for work product, one for
   notebook
2. Manifests are current at every commit, not retroactively patched
3. Chat exploration is captured in notebook with
   Q/Data/Finding/Impact structure
4. No blanket gitignore exceptions needed
5. Redo attempts are append-only in the notebook with downstream
   dependency flagged
