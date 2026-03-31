# Research skill redesign: iterative workflow with structured exploration logging

## Problem

The research skill assumes a linear 5-phase workflow (Orientation → Gene ID → Expression → Analysis → Documentation). The nitrogen stress MED4 analysis revealed that real research is iterative: orient once, then cycle through hypotheses, cross-experiment comparisons, and hypothesis revision. The linear structure caused several friction points:

1. **No exploration state management** — findings lived in chat context and had to be manually checkpointed to exploration logs halfway through the session.
2. **No cross-experiment comparison workflow** — the skill guides per-gene, per-experiment queries but the actual work needed gene-set-level cross-stress classification.
3. **No confidence/maturity marking** — no way to distinguish established findings from preliminary observations or speculation.
4. **MCP → API switching is vague** — the skill says "if truncated, switch to Python" but doesn't guide the common case where results need reshaping/aggregation.
5. **Platform comparability caveats get lost** — mentioned in anti-hallucination but not surfaced at the point of cross-experiment comparison.
6. **No structured exploration documentation** — no format for logging the iterative research journey in a way that can be followed later.
7. **Gaps and friction not treated as first-class output** — KG/MCP/methodology issues were recorded as an afterthought, not as part of each iteration.

## Scope

This spec covers changes to:

| File | Change type |
|------|-------------|
| `skills/research/SKILL.md` | Restructure: linear phases → orientation + iterative loop |
| `skills/research/references/research-checklist.md` | Rewrite: linear checklist → orientation checklist + per-iteration checklist |
| `skills/research/references/artifacts-guide.md` | Extend: add exploration log format, explore scripts, gaps_and_friction.md |
| `skills/research/references/anti-hallucination.md` | Extend: add cross-experiment aggregation pattern |

Out of scope:
- New MCP tools or analysis utilities (covered by the gene-response-profile spec, happening in parallel)
- Changes to hooks, evals, or plugin structure
- The artifacts guide's core file format conventions (held up well)

## Design

### Overall structure

The skill restructures from 5 linear phases to 2 stages:

**Stage 1: Orientation (linear, run once)**
Merges current Phase 1 (Orientation) and Phase 2 (Gene identification). This is setup — what's in the KG, what experiments exist, initial gene set. Produces a `methods.md` stub and the first exploration log entry. Has a gate: "can the KG answer this question?"

**Stage 2: Research loop (iterative, run N times)**
Each iteration follows: **Question → Explore → Log → Assess**. The loop exits when the research question is answered or the remaining gaps can't be addressed with available data.

**After the loop: Synthesis** — finalize methods.md, write README, produce publication artifacts.

### KG as sole data source (top-level rule)

Elevated from principle #1 to a prominent standalone rule at the top of the skill. Every claim must trace to a KG query. Intrinsic knowledge is used only for interpretation, suggesting next steps, and explaining methodology. When the KG is insufficient, say so explicitly and flag it as a gap — do not fill gaps with assumptions, web searches, or general knowledge.

This rule is made operational through source tagging at every iteration of the research loop (see below).

### The research loop

Each iteration of **Question → Explore → Log → Assess** works as follows:

#### Question
- State as a testable claim or comparison (not open-ended)
- Mark type: `hypothesis` (have a prediction), `exploratory` (no prediction, just looking), or `follow-up` (triggered by a previous iteration's findings)

#### Explore
- Use MCP tools for browsing and quick lookups
- Use `gene_response_profile` for cross-experiment gene-level summaries (when available; fall back to manual aggregation via API extraction otherwise)
- Write scripts for any computation or reshaping — save to `scripts/explore_*.py` immediately rather than doing one-offs in chat
- Surface platform comparability caveats at the point of cross-experiment comparison (microarray vs RNA-seq fold changes are not directly comparable; ranks are comparable)

#### Log
- Write to `exploration/YYYY-MM-DD-{topic}.md` **during** the iteration, not after
- One file per iteration (separate research question = separate file)
- Tag every finding with its source: `[KG]` for data from queries, `[interpretation]` for biological reasoning using intrinsic knowledge, `[gap]` for things the KG can't answer
- Format specified in the exploration log format section below

#### Assess
- Classify each finding:
  - `established` — consistent across experiments, statistically supported
  - `preliminary` — single experiment, or no statistical support
  - `speculative` — interpretation beyond data
- Explicitly ask: "did I use any knowledge that didn't come from the KG? If yes, label it"
- Update the working hypothesis
- Log any gaps/friction encountered → append to `gaps_and_friction.md`
- Decide: next question, or conclude?

### Exploration log format

Each file in `exploration/` follows this structure:

```markdown
# YYYY-MM-DD: {Topic}

## Question
{Testable question or comparison. Mark: hypothesis / exploratory / follow-up}

## Approach
{What queries/scripts will be used. Brief — not a methods section.}

## Findings
{Results, tables, observations. Each finding tagged:}
- [KG] cynA (PMM0370) is significant only under nitrogen_stress across all 30 MED4 experiments
- [interpretation] This suggests cynA is a reliable N-specific marker
- [gap] No physiological data (Fv/Fm) to confirm cell viability at 48h

## Assessment
{What did we learn? Classify findings:}
- **Established:** {consistent across experiments, statistically supported}
- **Preliminary:** {single experiment, or no stats}
- **Speculative:** {interpretation beyond data}

## Gaps and friction
{Any KG/MCP/methodology issues encountered — also appended to gaps_and_friction.md}

## Next
{What question does this lead to? Or: ready to conclude.}
```

The `## Next` section creates a chain between iterations. File A's Next says "classify all timepoints by stress specificity" → file B is `2026-03-30-timepoint-classification.md`.

### Methods.md as a living document

Methods.md grows alongside exploration but only absorbs settled conclusions:

- **Orientation** → creates methods.md stub with: research question, data scope, gene identification approach
- **Each iteration** → when the Assess step classifies findings as `established`, the relevant methods.md section is updated (data scope expands, gene set changes, caveats added)
- **Synthesis** → methods.md gets a final pass for coherence and publication-readiness

Exploration logs preserve the journey. Methods.md preserves the destination.

### Scripts convention

Scripts go into `scripts/` with a naming convention:

- `scripts/extract_*.py` — data extraction from KG (reusable)
- `scripts/analyze_*.py` — computation that produces results (reusable)
- `scripts/explore_*.py` — ad hoc iteration scripts (may be throwaway, kept for reproducibility)

Explore scripts are referenced from the exploration log (`## Approach: ran scripts/explore_cross_stress.py`). They can be promoted into the analysis utilities module if the pattern recurs.

### Gaps and friction as first-class output

`gaps_and_friction.md` is a standard artifact for every analysis, not an afterthought. Every iteration's Assess step asks whether any KG/MCP/methodology issues were encountered and appends them immediately. Each entry is recorded in both places:
- The exploration log entry (for context — what were you doing when you hit this?)
- `gaps_and_friction.md` (for the aggregated backlog that feeds into tool/KG development)

## Changes to reference docs

### research-checklist.md

Replace the current linear checklist with two sections:

1. **Orientation checklist** (run once) — largely unchanged from current Phase 1 + Phase 2. Gates: organisms present? experiments available? gene set resolved to locus tags?
2. **Research loop checklist** (per iteration):
   - [ ] Question stated as testable claim, type marked
   - [ ] Queries/scripts used (not chat reasoning for computation)
   - [ ] Findings tagged with source ([KG], [interpretation], [gap])
   - [ ] Exploration log entry written before moving on
   - [ ] Findings classified (established / preliminary / speculative)
   - [ ] Gaps appended to gaps_and_friction.md
   - [ ] Methods.md updated if established findings changed scope

### anti-hallucination.md

Add one pattern under Category 2 (Narrative fabrication):

**2.5 Cross-experiment aggregation without caveats**

*What happens:* Expression data from different platforms (microarray, RNA-seq) or different statistical tests (Goldenspike, Rockhopper, DESeq2) are compared as if they are directly comparable. "Gene X has log2FC 3.2 in study A and 1.1 in study B, so the response is stronger in A."

*Real example (nitrogen stress analysis):* Tolonen 2006 (microarray, Goldenspike) and Read 2017 (RNA-seq, Rockhopper) fold changes for the same genes were compared. Different platforms have different dynamic ranges — a log2FC of 3 on microarray is not the same as log2FC 3 on RNA-seq.

*Prevention:*
- Compare direction (up/down) and rank across platforms, not magnitude
- Compare log2FC magnitudes only within the same platform and study
- When presenting cross-experiment tables, include a platform column and note the caveat
- Use `gene_response_profile` rank fields (comparable across platforms) instead of raw fold changes for cross-study comparisons

### artifacts-guide.md

Add to the standard directory structure:

```
analyses/{analysis_name}/
├── exploration/       # Exploration logs (one per iteration)
├── data/              # Staged data from KG (CSV/TSV)
├── scripts/           # Python scripts (extract, analyze, explore)
├── results/           # Outputs: tables, figures, statistics
├── README.md          # Summary, key findings, file index
├── methods.md         # Publication-ready methods (living document)
├── gaps_and_friction.md  # KG/MCP/methodology issues log
└── references.md      # Data sources and citations
```

Add the exploration log format as a subsection (or pointer to the skill's description of it).

Add the scripts naming convention: `extract_*`, `analyze_*`, `explore_*`.

## Changes to SKILL.md structure

The restructured SKILL.md has these sections:

1. **KG is the sole data source** (top-level rule, standalone — not one of several principles)
2. **Core principles** (reproducibility, locus tags, artifacts not answers — reframed to include exploration logs as artifacts)
3. **Stage 1: Orientation** (linear — current phases 1+2 merged, with gate)
4. **Stage 2: Research loop** (Question → Explore → Log → Assess, with source tagging and confidence classification)
5. **Synthesis** (finalize methods.md, README, publication artifacts)
6. **MCP vs Python package** (kept, with added explore script convention and clearer switching guidance)
7. **Statistical rigor** (unchanged)
8. **References and citations** (unchanged)
9. **Review protocol** (updated to cover exploration logs, not just final artifacts)

## Dependencies

- The `gene_response_profile` MCP tool and analysis utilities (separate spec, parallel development) will improve the Explore step but are not required. The skill works with current tools and notes where the new tools will help.
- No changes to hooks, evals, or plugin structure.

## Validation

The next analysis session (coculture comparison / paper rewrite continuation) serves as the test. Success criteria:
- Exploration logs written during the session, not reconstructed
- Every finding tagged with source
- Gaps recorded as they're encountered
- Methods.md updated incrementally
- The session produces a clear chain of exploration log entries that someone could follow later
