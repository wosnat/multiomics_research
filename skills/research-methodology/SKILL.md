---
name: research-methodology
description: Non-negotiable rules for multi-omics KG research. Load when answering biological questions, analyzing expression data, planning or brainstorming a research analysis, reviewing results, or working with the multiomics KG in any capacity. Must be loaded BEFORE brainstorming — the rules shape analysis design. Reference skill — provides domain rules, not process.
---

# Multi-omics research methodology

These rules apply to ALL research work — brainstorming, planning,
and execution. They are non-negotiable.

> **Load this skill BEFORE brainstorming.** The notebook structure,
> artifact requirements, and methodology rules shape the entire
> analysis design. Loading after the spec is written means
> retrofitting — which is what happened in the N-limitation analysis.

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

## Rule 8: Research notebook, not pipeline

Every analysis is driven by an interactive notebook where each
step is recorded, explored with the researcher, and approved
before the next step. Implementation can be fast or delegated;
quality control is always interactive.

See [Research notebook](references/research-notebook.md) for the
notebook format, QC checkpoint requirements, step cycle, and code
lifecycle rules.

## References — read on demand, not all at once

- [KG rules](references/kg-rules.md) — read when scoping a new analysis or uncertain about data sourcing vs literature
- [Gene identity](references/gene-identity.md) — read when working with gene names, paralogs, or orthologs
- [Artifacts guide](references/artifacts.md) — read when setting up an analysis directory or unsure about file structure
- [Anti-hallucination](references/anti-hallucination.md) — read before presenting findings or when self-checking results
- [Python API guide](references/python-api-guide.md) — read before writing extraction or analysis scripts
- [Statistical rigor](references/statistical-rigor.md) — read when computing enrichment, comparing across studies, or reporting p-values
- [Research notebook](references/research-notebook.md) — read when starting or resuming an interactive analysis
