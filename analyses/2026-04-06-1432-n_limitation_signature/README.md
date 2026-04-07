# N-Limitation Signature Analysis

**Question:** Is Prochlorococcus nitrogen-limited in Weissberg 2025's axenic
and coculture experiments, and how does the severity differ over time?

**Approach:** Reference signature scoring (Approach A of 3).

## Key findings

- **Core signature: 198 genes** (83 up, 115 down) from the concordant
  intersection of Tolonen 2006 (microarray) and Read 2017 (RNA-seq). Extended
  signature adds 169 lower-confidence genes for a total of 367.

- **RNA-seq axenic: strongly N-limited.** Day 14 axenic MED4 scores
  hit_rate=0.818, rank_score=0.743, permutation p=0.000. The transcriptome
  mirrors the reference N-deprivation response. [KG]

- **RNA-seq coculture: no transcriptomic N-limitation signal.** All
  coculture timepoints (days 18–89) score near zero and are not significant
  (p=0.062–0.672). The presence of Alteromonas eliminates the MED4
  transcriptional stress response. [interpretation]

- **Proteomics: both conditions show N-limitation at protein level.**
  Axenic MED4 is significant at day 31 and day 89 (rank scores 0.23–0.42).
  Coculture MED4 is significant at all timepoints (rank scores 0.41–0.76),
  with consistently higher scores than axenic. [interpretation]

- **RNA/protein discordance in coculture is the main biological finding.**
  The coculture transcriptome shows no N-limitation signature while the
  proteome shows persistent, strong N-limitation. This suggests Alteromonas
  suppresses transcriptional N-stress signaling in MED4 while the underlying
  low-N environment continues to shape the proteome. [interpretation]

- **Alteromonas N-recycling: 65 of 171 candidate genes significantly DE.**
  N-recycling activity increases from day 18 through day 60, then plateaus.
  Both up- and down-regulated genes contribute. [KG]

## What to do next

**Priority 1: Walkthrough and verification.** The analysis was built fast. Before extending it, walk through each step: verify the signature makes biological sense, check the scoring outputs against known genes (e.g., does glnA score high in axenic? does it score low in coculture?), and confirm the RNA/protein discordance isn't an artifact of different gene coverage between platforms.

**Priority 2: Investigate the RNA/protein discordance.** This is the novel finding and wasn't in the original plan. Key questions:
- Which specific genes drive the proteomics signal in coculture? Are they the same genes that drive the RNA-seq signal in axenic, or a different subset?
- Is the discordance uniform across pathways, or concentrated in specific functional categories (e.g., translation machinery has slow protein turnover)?
- Can we quantify per-gene RNA/protein concordance at matched timepoints?

**Priority 3: Approach B (pathway enrichment).** The backlog in the spec lists the statistical tests (Fisher's exact, FDR, Mann-Whitney) but needs a full brainstorm given Approach A's results. Key decisions needed: which ontologies for pathways, what background set for enrichment, how to handle the RNA/protein split.

**Priority 4: Approach C (hybrid).** Builds on A + B. The Tier 1/Tier 2 partitioning and multi-omics concordance are well-defined in the spec but should be re-scoped to focus on the discordance finding.

**Spec:** `docs/superpowers/specs/2026-04-06-n-limitation-signature-analysis-design.md`
**Plan (Approach A):** `docs/superpowers/plans/2026-04-06-n-limitation-signature-approach-a.md`
**Gaps and retrospective:** `gaps_and_friction.md`

## File index

- `methods.md` — publication-ready methods
- `gaps_and_friction.md` — issues log
- `references.md` — citations and versions
- `data/` — staged KG extracts
- `scripts/` — extraction and analysis scripts
- `sig_utils/` — reusable scoring utilities
- `results/` — output tables and figures
- `exploration/` — exploration logs

## Exploration logs

- [2026-04-06: Signature building](exploration/2026-04-06-signature-building.md) — Tolonen/Read intersection, core and extended signature construction
