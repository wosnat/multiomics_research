---
name: research-methodology
description: Use when answering biological questions, analyzing expression data, planning or brainstorming a research analysis, reviewing results, or working with the multiomics KG in any capacity. Non-negotiable domain rules and research process for multi-omics KG work. CRITICAL — load BEFORE invoking brainstorming for step 1 of an analysis; loading after the step-1 dialogue means retrofitting.
---

# Multi-omics research methodology

These rules apply to ALL research work — dialogue, analysis,
execution. They are non-negotiable.

> **Load this skill BEFORE invoking `superpowers:brainstorming` for step 1
> of an analysis.** Step 1's capture location and terminal behavior are
> overridden by this skill (see Rule 8 and
> [research-notebook.md — Using brainstorming for step 1](references/research-notebook.md)).
> Loading after step 1 is committed means retrofitting.

## Rule 1: KG is the sole data source

Every claim must trace to a KG query. Never rely on intrinsic
knowledge for data — gene names, expression values, experiment
details, organism properties, ortholog assignments. Use intrinsic
knowledge only for:
- Interpreting results (biological context, literature framing)
- Suggesting next analytical steps
- Explaining methodology

When the KG is insufficient, say so explicitly and flag it as a gap.
Do not fill gaps with assumptions, web searches, or general knowledge.

See [KG rules](references/kg-rules.md) for common gaps, scoping
checks, and when to use MCP vs Python API.

## Rule 2: Locus tags, not gene names

Gene names are ambiguous. Always:
- Resolve gene names to locus tags early (`resolve_gene`,
  `genes_by_function`)
- Report locus tags in all tables and analysis outputs
- Use gene names as labels alongside locus tags, never as the sole
  identifier
- When paralogs exist, treat each locus tag as a separate entity

See [Gene identity](references/gene-identity.md) for paralog handling
and ortholog cluster rules.

## Rule 3: Source tagging

Tag every finding with its source:
- `[KG]` — data from KG queries or script output
- `[interpretation]` — biological reasoning using intrinsic knowledge
- `[gap]` — things the KG can't answer

## Rule 4: Artifacts, not answers

Research questions produce files, not chat messages. Chat is for
reasoning, planning, and interpretation. Data, statistics, figures,
and exploration logs go to disk.

See [Artifacts guide](references/artifacts.md) for directory structure,
exploration log format, and file naming conventions.

## Rule 5: Scripts over chat reasoning

Chat-computed statistics are unreproducible — they can't be rerun,
verified, or shared. Computations go in `.py` files, not in chat
responses. Data staged to files before analysis. No manual steps —
if it can't be scripted, document it as a limitation.

See [Python API guide](references/python-api-guide.md) for imports,
return structure, and scripting patterns.

**Exception: interactive discovery steps.** Steps that are
inherently exploratory (browsing available data, classifying
experiments, scoping what the KG contains) may be done
interactively rather than scripted. These steps must still produce
a frozen output file (CSV) and a notebook entry documenting the
reasoning. Computations — statistics, scores, enrichment — always
go in scripts. See [Research notebook — Interactive discovery
steps](references/research-notebook.md) for the pattern.

## Rule 6: Statistical rigor

Expression data without proper statistical treatment produces
misleading claims — wrong background sets, uncorrected multiple
testing, and magnitude comparisons across platforms are common
failure modes.

See [Statistical rigor](references/statistical-rigor.md) for what
the KG provides, what you must compute in scripts, and what to flag.

## Rule 7: Don't hallucinate

See [Anti-hallucination](references/anti-hallucination.md) for
concrete failure modes and prevention patterns.

## Rule 8: The 6-step research flow

Every research analysis is structured as 6 steps, each advancing
through the rhythm **do → show → explore → decide** with an atomic
researcher-approval gate at decide close.

1. **Research question** — user prompt + clarifying questions →
   locked question
2. **KG entries** — relevant publications, experiments, organisms,
   data types identified from the KG
3. **Analysis framing** — (a) selection: publications / experiments /
   organisms / data types; (b) framing: hypothesis, target, positive
   and negative controls, expected outcome — all in KG terms
4. **Methods** — pick one item from step 3 as driving example;
   produce an ad-hoc Python module
5. **Analyze** — run the method; produce scored outputs, figures,
   tables
6. **Evaluate** — assess against framing; harvest caveats; finalize
   paper

Steps 1–3 form the research proposal (locked at end of step 3
decide). Steps 4–6 execute against it. Step 1 is a conversation
(typically uses `superpowers:brainstorming` with overrides); steps
2–6 produce scripts, data, figures, and QC alongside `notebook.md`.

### Just-in-time formalization

Plans, framings, predictions, metrics, decisions, and caveats enter
the analysis **only when the data demands them**. Look at the data
before drafting the plan; start with the simplest plan; expand only
when a specific finding forces it. Step-1 sub-narratives belong to
step 3, after step 2 has shown what's in the KG.

See [Research notebook — Just-in-time formalization](references/research-notebook.md#just-in-time-formalization) for per-step concrete rules.

### What replaces the old spec/plan/decisions machinery

- **No `spec.md` / `plan.md`** — the 6-step flow replaces the
  monolithic plan; each step's intent is captured in its own
  `notebook.md`
- **No global `decisions.md`** — decisions live in the step's
  `notebook.md` where they were forced
- **No `hypotheses.md`** — the hypothesis lives in step 3's
  `notebook.md` (the framing)
- **No cross-file labels** (H1/P3/D8/T4) — labels are
  document-scoped; labels used must be paired with a short readable
  name in the same paragraph
- **Per-step `notebook.md` + single `paper.md`** at the analysis
  root replace the single exploration notebook + methods.md +
  caveats.md
- **`gaps_and_friction.md`** (transitional) captures methodology /
  KG / tooling friction, distinct from decisions

See [Step protocol](references/step-protocol.md) for commit timing,
decide-gate checklist, hard gates, and redo path. See
[Research notebook](references/research-notebook.md) for notebook
format, `paper.md` growth, `gaps_and_friction.md`, and the
brainstorming override for step 1. See
[Artifacts guide](references/artifacts.md) for directory structure
and scaffold creation.

## Rule 9: Plain language; describe before interpreting

Write in plain English. Numbers and direction first, interpretation
second. Don't reach for fancy vocabulary when ordinary words work.

**Banned in steps 1–5** (these signal commitment before the analysis
has earned it): "striking", "massive", "central finding",
"biologically loaded", "biologically explosive", "reframes", "rich",
"hand-wavy" as praise. Reserve interpretive vocabulary for step 6.

**Describe before interpreting.** Write the numbers and direction
first ("5 markers show RNA log2FC < 0 and protein log2FC > 0 at
coculture days 31–89"). Tag interpretation `[interpretation]` and
list plausible alternatives.

**Verify before generalizing.** Words like "all", "every",
"systematically", "primarily", "no genes" require a query against
the data file, not a heatmap glance. If you write "all 5 are UP",
check all 5 in the data first.

**Proposals must cite specific data or friction.** When suggesting
a fix, framework, or methodology change, point to the data or
incident that motivates it. One occurrence is a note; process
change needs the same friction in two analyses. Don't enumerate
speculative proposals.

See [Anti-hallucination — 2.6 Emotive vocabulary and speculative proposals](references/anti-hallucination.md#26-emotive-vocabulary-and-speculative-proposals).

## References — read on demand, not all at once

- [Step protocol](references/step-protocol.md) — read at the start of every analysis execution; owns commit timing, decide-gate checklist, hard gates, redo path
- [Research notebook](references/research-notebook.md) — read when starting or resuming an analysis; owns notebook format, paper.md growth, gaps_and_friction.md, brainstorming override for step 1
- [Artifacts guide](references/artifacts.md) — read at scaffold time or when unsure about directory structure, per-step folders, QC naming
- [Anti-hallucination](references/anti-hallucination.md) — read before presenting findings; covers tool-schema claims, publication attribution, field semantics, and the other anti-hallucination patterns
- [KG rules](references/kg-rules.md) — read when scoping a new analysis or uncertain about data sourcing vs literature
- [Gene identity](references/gene-identity.md) — read when working with gene names, paralogs, or orthologs
- [Python API guide](references/python-api-guide.md) — read before writing extraction or analysis scripts
- [Statistical rigor](references/statistical-rigor.md) — read when computing enrichment, comparing across studies, or reporting p-values
