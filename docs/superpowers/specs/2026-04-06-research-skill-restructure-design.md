# Research skill restructure: superpowers-driven methodology

**Date:** 2026-04-06
**Status:** Draft

## Problem

The current research skill (`skills/research/`) encodes two things:

1. **Domain rules** — KG-as-sole-source, source tagging, locus tags, anti-hallucination, artifact structure, statistical rigor
2. **Process workflow** — 3-phase structure (Pre-flight → Orientation → Research loop → Synthesis), gates, checklists

Claude doesn't follow the process reliably. It treats the checklist as suggestions, skips gates, and goes off into ad hoc scripts. Meanwhile, the superpowers framework (brainstorm → spec → plan → execute) works well — Claude follows it because the staged approach has hard gates, interactive checkpoints, and concrete artifacts at each stage.

## Design

**Principle:** Superpowers owns the process. The research skill owns the domain knowledge. Don't reinvent process machinery — leverage what works.

### Two layers

| Layer | What | When loaded | How |
|-------|------|-------------|-----|
| **Research methodology** | Domain rules that apply to all research work | Always (triggered by CLAUDE.md) | Reference skill — no process, just rules |
| **Recipes** | Concrete protocols for specific analysis types | On demand | Individual skills, one per method |

### Layer 1: `skills/research-methodology/`

A reference skill containing non-negotiable rules for working with the multiomics KG. No process steps, no checklists — just the rules themselves.

```
skills/research-methodology/
  SKILL.md                        # Top-level rules + pointers to references
  references/
    kg-rules.md                   # KG-as-sole-source, source tagging, what to do when KG is insufficient
    gene-identity.md              # Locus tags over gene names, paralog handling, ortholog cluster rules
    artifacts.md                  # Directory structure, exploration logs, methods.md, gaps_and_friction.md
    anti-hallucination.md         # Concrete failure modes and prevention (exists today, mostly unchanged)
    python-api-guide.md           # Scripting guide: imports, return structure, pre-script checklist (exists today, unchanged)
    statistical-rigor.md          # What KG provides, what you compute in scripts, what to flag
```

**SKILL.md content:** Short. States the top-level rules directly (KG is sole source, locus tags not gene names, source tagging, artifacts not answers). Points to reference docs for details. No workflow, no phases, no gates.

**Trigger:** CLAUDE.md instructs Claude to load this skill whenever working with the multiomics KG or answering biological questions. This means the rules are in context during brainstorming, planning, and execution — not just during a `/research` invocation.

### Layer 2: `skills/recipes/`

Individual skills for specific analysis methods. Each recipe is self-contained: what the method does, when to use it, step-by-step protocol, expected inputs/outputs.

```
skills/recipes/
  enrichment/
    SKILL.md                      # Pathway/GO enrichment analysis protocol
  response-matrix/
    SKILL.md                      # Gene × treatment comparison using response_matrix()
  conservation/
    SKILL.md                      # Cross-organism ortholog conservation analysis
  ...                             # More recipes added over time
```

**Recipes don't exist yet.** They'll be developed as analysis patterns stabilize. Each recipe references the methodology skill for general rules and adds method-specific instructions.

**Trigger:** Invoked by name when needed — either by the user (`/enrichment`) or referenced in a plan step ("invoke the enrichment recipe for this step").

### CLAUDE.md changes

Replace the current research skill references with:

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

### How it flows

1. User asks a biological question
2. CLAUDE.md triggers loading the research-methodology skill → domain rules in context
3. User runs `/brainstorm` → superpowers brainstorming, but with methodology rules loaded. Scoping naturally includes KG-specific steps (check organisms, experiments, gene resolution) because the rules say to
4. User runs `/write-plan` → plan steps reference methodology rules as constraints. Plan identifies which recipe skills are needed for specific steps
5. User runs `/execute-plan` → each step follows the plan. Methodology rules are in context. Recipe skills loaded as needed per step

### What gets removed

The current `skills/research/` directory is replaced entirely:

- **Process workflow** (Phase 0, Stage 1, Stage 2, Synthesis) → removed. Superpowers handles this
- **Gates and checklists** → removed from the skill. The plan's checkpoints serve this role, and they're enforceable because the user approves them
- **Research checklist** (`references/research-checklist.md`) → removed. Its content either moves to methodology references (the domain rules parts) or becomes unnecessary (the process parts)
- **Domain rules** → reorganized into `research-methodology/references/`
- **Anti-hallucination** → moves to `research-methodology/references/` (mostly unchanged)
- **Python API guide** → moves to `research-methodology/references/` (unchanged)
- **Artifacts guide** → reorganized into `research-methodology/references/artifacts.md`

### What's new

- **Recipes directory** (`skills/recipes/`) — empty initially, populated as analysis patterns emerge
- **CLAUDE.md pointer** — short section triggering methodology skill loading
- **Separation of concerns** — methodology never tells you what step to do next; superpowers never tells you what a locus tag is

## Migration

### Content mapping from current skill

| Current location | Destination | Changes |
|---|---|---|
| `SKILL.md` §KG is the sole data source | `research-methodology/SKILL.md` + `references/kg-rules.md` | Remove workflow references, keep rules |
| `SKILL.md` §Core principles | `research-methodology/SKILL.md` | Condensed as top-level rules |
| `SKILL.md` §3-phase workflow | Removed | Superpowers handles process |
| `SKILL.md` §MCP vs Python package | `references/python-api-guide.md` | Merge with existing guide |
| `SKILL.md` §Statistical rigor | `references/statistical-rigor.md` | Standalone reference |
| `SKILL.md` §References and citations | `references/artifacts.md` | Merge with artifacts guide |
| `SKILL.md` §Review protocol | Removed | Becomes a plan checkpoint or a recipe |
| `references/research-checklist.md` | Removed | Domain parts → methodology refs; process parts → superpowers |
| `references/artifacts-guide.md` | `references/artifacts.md` | Reorganized, process language removed |
| `references/anti-hallucination.md` | `references/anti-hallucination.md` | Mostly unchanged |
| `references/python-api-guide.md` | `references/python-api-guide.md` | Unchanged |

### What stays the same

- Anti-hallucination patterns — proven, stable
- Python API guide — proven, stable
- Artifact directory structure (`analyses/{name}/exploration/`, `data/`, `scripts/`, `results/`) — this is a rule, not a process
- Source tagging (`[KG]`, `[interpretation]`, `[gap]`) — this is a rule, not a process
- Exploration log format — this is an artifact spec, not a process

## Success criteria

1. Research work follows the superpowers flow (brainstorm → plan → execute) without a separate process encoded in the skill
2. Domain rules (locus tags, source tagging, KG-as-sole-source) are consistently followed because they're in context during all phases
3. Recipes can be added incrementally without touching the methodology skill
4. CLAUDE.md is concise — pointer to methodology, not pages of rules
5. No duplicate process machinery between research skill and superpowers

## Open questions

1. **Review protocol** — currently a checklist in the research skill. Should it become a recipe (`/research-review`), a plan checkpoint template, or something else?
2. **Recipe granularity** — how specific should recipes be? "Enrichment analysis" or "GO enrichment with Fisher's exact test on MED4 nitrogen stress genes"? Probably the former, with the plan adding specificity.
3. **Methodology skill trigger** — CLAUDE.md says "invoke when working with the KG." This relies on Claude reading CLAUDE.md and invoking the skill itself (same mechanism superpowers uses via its `using-superpowers` skill). No hook needed — CLAUDE.md is always in context.
