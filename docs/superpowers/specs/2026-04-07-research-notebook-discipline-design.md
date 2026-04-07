# Research Notebook Discipline — Design Spec

**Date:** 2026-04-07
**Motivation:** The N-limitation signature analysis (2026-04-06) produced good results but the researcher couldn't keep up with the implementation pace. Claude optimized for throughput when the researcher needed comprehension. The result was a black box — correct outputs without understanding of the intermediate steps.

**Core insight:** The brainstorming process (one section at a time, researcher approves before proceeding) is the right model for analysis execution too. The skill should enforce this.

## Problem Statement

The current research-methodology skill says "artifacts, not answers" — produce files, not chat. This is correct but incomplete. It optimizes for *output* when research requires *understanding*. Specifically:

1. Implementation runs too fast. Once a plan exists, Claude dispatches subagents and produces 12 tasks of output without researcher checkpoints.
2. The exploration log is a retrospective summary, not a chronological record of what happened.
3. Analytical decisions (what's in a gene set, why a gene was excluded, whether a score is reasonable) happen inside scripts or subagent contexts, invisible to the researcher.
4. QC and interactive exploration (walking through specific genes, verifying logic) aren't part of the defined workflow.

## Design

### New Rule 8: Research notebook, not pipeline

Added to SKILL.md as a top-level rule with a new reference file `references/research-notebook.md`.

**Rule text:**

> Every analysis is driven by an interactive research notebook — a chronological log where each analytical step is recorded, inspected, explored with the researcher, and approved before the next step proceeds. Implementation can be fast or delegated; quality control and exploration are always interactive.

### The step cycle: do → show → explore → decide

Every step that produces data or analytical output follows this cycle:

1. **Do** — run the script, produce artifacts (CSV, figures, summaries). Can be delegated to a subagent.
2. **Show** — present QC diagnostics: what went in, what came out, what's missing, sanity checks.
3. **Explore** — interactive walkthrough in chat: pick specific genes/conditions, verify the logic on concrete examples, ask and answer questions.
4. **Decide** — researcher says "continue", "explain X", or "redo with Y". The decision is logged. Only then does the next step begin.

**No step that produces data or analytical output may proceed without completing the full cycle.**

Mechanical tasks (formatting, plotting from existing data, file reorganization) can skip the explore phase but still need show + decide.

### What a "step" is

A step is:
- A Python script with explicit command-line invocation
- Output artifacts (CSV, figures, summary files)
- A log entry in the notebook capturing all four phases

### Notebook format

One notebook per analysis: `exploration/YYYY-MM-DD-notebook.md`

Chronological and append-only. Reruns add new entries, not overwrites.

Each entry:

```markdown
---

## YYYY-MM-DD HH:MM — Step N: {description}

### Command
```bash
uv run scripts/build_signature.py --tolonen data/de_reference_tolonen_ndep.csv --read data/de_reference_read_ndep.csv --output data/core_signature_genes.csv
​```

### Inputs
- [de_reference_tolonen_ndep.csv](../data/de_reference_tolonen_ndep.csv) — 1,247 genes, 6 timepoints
- [de_reference_read_ndep.csv](../data/de_reference_read_ndep.csv) — 892 genes, 4 timepoints

### Outputs
- [core_signature_genes.csv](../data/core_signature_genes.csv) — 198 genes (83 up, 115 down)
- [extended_signature_genes.csv](../data/extended_signature_genes.csv) — 367 genes

### QC
- Gene counts at each filter step: 1,247 Tolonen → 892 Read → 210 in both → 198 concordant direction → 12 discordant excluded
- Known markers present: glnA ✓, amt1 ✓, nirA ✓
- Read table_scope=filtered_subset — 72 genes in Tolonen but absent from Read

### Exploration
- Walked through glnA (PMM0370): Tolonen rank_up=5 at 24h, Read rank_up=12 at 12h. Both up. Core ✓.
- Checked psbA (PMM1547): rank_down=3 (Tolonen), rank_down=8 (Read). Photosystem shutdown under N-stress — expected.
- Asked: why is PMM0842 (hypothetical) in core? → Significant up in both, rank_up=45/30. No annotation. Keep — signature is data-driven.
- Looked at excluded genes: PMM0501 up in Tolonen, down in Read. Genuinely discordant — correct to exclude.

### Decision
Core signature validated. The 72 tolonen-only/read-absent genes go to extended, flagged as caveat. Proceed to scoring.
```

### What the notebook captures

Everything:
- **Spec walkthrough** — the initial section-by-section review of the analysis design (questions, decisions, rationale)
- **Each analytical step** — command, inputs, outputs, QC, exploration, decision
- **Reruns** — new entry with "Why" section explaining what changed and the new output
- **Chat-based follow-up** — interpretation discussions, "what about gene X?" explorations, with actual data points
- **Questions and decisions** — including dead ends and rejected alternatives
- **Relative links** to all artifacts

### Rerun entries

When a step is rerun after a change:

```markdown
---

## YYYY-MM-DD HH:MM — Step 2 rerun: updated direction logic

### Why
Researcher flagged that majority-direction tie-breaking was arbitrary. Changed to: exclude ties.

### Changes
- Modified `build_signature.py` line 45: ties now excluded instead of broken by idxmax
- 3 genes affected (PMM0xxx, PMM0yyy, PMM0zzz)

### Command
```bash
uv run scripts/build_signature.py --tolonen data/de_reference_tolonen_ndep.csv --read data/de_reference_read_ndep.csv --output data/core_signature_genes.csv --exclude-ties
​```

### Outputs
- [core_signature_genes.csv](../data/core_signature_genes.csv) — 195 genes (updated)

### Decision
Proceed with 195.
```

### QC checkpoint types

What to show depends on the step type:

**Data extraction:**
- Row count, gene count, timepoint/condition count
- Sample rows (first 5 or curated edge cases)
- What's missing: expected genes absent, unexpected metadata (e.g., `timepoint=single`)
- One-sentence summary: "This is X genes across Y conditions from Z experiments"

**Gene selection / filtering:**
- Counts at each filter step: started with → filter 1 → filter 2 → final
- Sample excluded genes with reason
- Sample included genes — known markers present?
- Flag surprises

**Computation / metric:**
- Worked example: 2-3 genes through the formula with actual numbers
- Summary statistics (distribution, range, NaNs)
- Sanity check against known biology

**Scoring / comparison:**
- Full results table in markdown (not prose summary)
- Best/worst scores, surprises
- Cross-condition comparison with expectation check

### Manifest updates

Every step that produces output files updates the relevant manifest (`data/DATA_MANIFEST.md` or `results/RESULTS_MANIFEST.md`) immediately — not retroactively at the end.

## Changes to Existing Files

### SKILL.md
- Add load-timing advisory at top: "Load this skill BEFORE brainstorming, not after."
- Add Rule 8: Research notebook, not pipeline (pointing to `references/research-notebook.md`)

### references/artifacts.md
- Update directory structure to include required docs: `data/DATA_MANIFEST.md`, `results/RESULTS_MANIFEST.md`, `decisions.md`, `caveats.md`
- Add concise templates for each
- Replace exploration log template with pointer to notebook format in `references/research-notebook.md`
- Add rule: manifests updated with every step

### references/statistical-rigor.md
- Add: "Every formula in a spec must include a worked example with concrete numbers: input → computation → output."

### references/python-api-guide.md
- Add: "When extracting DE data with `verbose=True`, join `product` and `gene_category` into gene-level outputs immediately."
- Add: "After extraction, check for experiments with `timepoint=single` or null — confirm handled correctly downstream."

### CLAUDE.md
- Update research methodology section: load skill before brainstorming.

## New Files

### references/research-notebook.md
The full notebook discipline reference: step cycle, notebook format, checkpoint types, manifest rules. Content as designed above.

## Out of Scope

- MCP/KG changes (table_scope in DE rows, background_gene_count, score_gene_set utility) — tracked in gaps_and_friction.md for the explorer repo
- Superpowers process changes (no code in plans, mandatory toy-data verification, subagent review enforcement) — separate from this skill
- Recipe-level checkpoint definitions — can be added per-recipe later, building on the general notebook discipline
