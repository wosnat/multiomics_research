---
name: research-methodology
description: Non-negotiable rules for multi-omics KG research. Load when answering biological questions, analyzing expression data, or working with the multiomics KG. Reference skill — provides domain rules, not process.
---

# Multi-omics research methodology

These rules apply to ALL research work — brainstorming, planning,
and execution. They are non-negotiable.

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

See [KG rules](references/kg-rules.md) for details on source tagging,
common gaps, and what to do when data is missing.

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

Computations go in `.py` files, not in chat responses. Data staged
to files before analysis. No manual steps — if it can't be scripted,
document it as a limitation.

See [Python API guide](references/python-api-guide.md) for imports,
return structure, and scripting patterns.

## Rule 6: Statistical rigor

See [Statistical rigor](references/statistical-rigor.md) for what
the KG provides, what you must compute in scripts, and what to flag.

## Rule 7: Don't hallucinate

See [Anti-hallucination](references/anti-hallucination.md) for
concrete failure modes and prevention patterns.

## References

- [KG rules](references/kg-rules.md) — source tagging, common gaps, KG insufficiency
- [Gene identity](references/gene-identity.md) — locus tags, paralogs, ortholog clusters
- [Artifacts guide](references/artifacts.md) — directory structure, exploration logs, methods.md
- [Anti-hallucination](references/anti-hallucination.md) — failure modes and prevention
- [Python API guide](references/python-api-guide.md) — scripting with multiomics_explorer
- [Statistical rigor](references/statistical-rigor.md) — KG statistics, what to compute, what to flag
