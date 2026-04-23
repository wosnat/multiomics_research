# Analysis flow redesign

**Date:** 2026-04-23
**Status:** design
**Motivation:** The pathway_enrichment_b2 analysis was stopped mid-execution because the workflow around it — not the data — was unusable. Steps 1–4 produced solid artifacts; steps 4 decide and 5 were never written, and the researcher got lost in a 104KB plan with H1–H8 / P1–P4 / D1–D8 / T1–T4 labels accumulated before the data arrived. This spec replaces the superpowers-based plan/execute workflow with a leaner academic research methodology.

## Headline principle: just-in-time formalization

Terms, predictions, metrics, stability checks, decisions, and caveats enter the analysis **only when the data demands them**. Nothing is enumerated in advance "just in case."

Concrete rules that follow from this principle:

1. **Step 3 framing has a floor, not a template.** Minimum required content: the hypothesis in prose, what success means operationally, positive and negative controls selected from the KG. Nothing else is required.
2. **Preregistration is optional and minimal.** If confirmation bias is a real risk, preregister 1–3 named predictions — not a 4×4 matrix.
3. **Stability checks are added when a specific result triggers them**, not planned up front.
4. **Decisions are written when forced, not anticipated.** No "we may decide X" placeholders.
5. **Methods (step 4) stay minimal** — ad-hoc Python module with exactly what this analysis needs.
6. **Caveats are harvested at step 6**, from what actually happened; not pre-cataloged.

## Changes at a glance

**Dropped from default workflow:**
- `superpowers:writing-plans` — no more 100KB 14-task plans
- `superpowers:executing-plans` — replaced by per-step decide-gate
- `superpowers:subagent-driven-development`
- Separate `spec.md` / `plan.md` / `decisions.md` / `hypotheses.md` documents — content merges into per-step `notebook.md` + single `paper.md`
- Global symbol taxonomies (H1–H8, P1–P4, T1–T4, M1–M4, D1–D8, C1–C9)

**Kept as optional, invoked on demand:**
- `superpowers:brainstorming` — for exploring the research question (step 1)
- `superpowers:verification-before-completion` — before claiming a step is done
- `superpowers:systematic-debugging` — when code breaks
- `superpowers:requesting-code-review` — for substantial method modules

**Replaced by:** the 6-step academic research flow below.

## The 6-step flow

1. **Research question** — formulate the question in prose
2. **KG entries** — identify relevant publications, experiments, organisms, data types in the KG
3. **Analysis framing** — (a) selection: publications, experiments, organisms, data types (DE, cluster, ontologies); (b) framing: hypothesis, target, positive and negative controls, expected outcome — all in KG terms
4. **Methods** — pick one item from step 3 as a driving example; select or generate an analysis method (statistical test, score, ...); produce an ad-hoc Python module
5. **Analyze** — run the method; produce scored outputs, figures, tables
6. **Evaluate** — assess results against framing; write up caveats; finalize paper

Steps 1–3 are collectively the **research proposal**. Locked at the end of step 3 decide. Steps 4–6 execute against it.

## Intra-step rhythm: do → show → explore → decide

Every step — including step 1 — advances through these four phases:

- **do** — run the step's work; produce outputs in `scripts/`, `data/`, `figures/`
- **show** — populate `notebook.md` with tables, figures, observations; display to researcher
- **explore** — investigate anomalies or surprises; add `qc_*.py` checks as needed; add `qc_*.csv` / `qc_*.png` as needed
- **decide** — finalize notebook.md; minimal checklist (below); **researcher approval gate**; git commit; advance

The decide phase is the only formal gate between steps.

### Decide-gate checklist (radically minimal)

Each step's decide phase must produce, in its `notebook.md`:

- **Outputs produced** — filenames in `scripts/`, `data/`, `figures/`
- **QC gate** — what was checked → result (one line per check)
- **Decisions made this step** — prose + date, if any; omit the section if none
- **Advance rationale** — one line, why this step is ready to close

Claude then pauses, displays the state, waits for researcher approval, commits the step, and begins the next.

The checklist must stay this minimal. It is not a template to extend with optional fields. Inflation of this list would reintroduce the same premature-formalization failure the redesign exists to prevent.

## Directory structure per analysis

```
analyses/YYYY-MM-DD-<slug>/
  paper.md                        # single academic-style writeup, grows from step 1
  1_question/
    notebook.md
    scripts/
    data/
    figures/
  2_kg_selection/
    notebook.md, scripts/, data/, figures/
  3_analysis_framing/
    ...
  4_methods/
    notebook.md
    <module_name>.py              # ad-hoc method module
    scripts/, data/, figures/
  5_analyze/
    ...
  6_evaluate/
    ...
```

**QC artifacts** live in the same `scripts/` / `data/` / `figures/` folders as main outputs, using a `qc_` filename prefix (e.g., `scripts/qc_check_tp_coverage.py`, `figures/qc_tp_coverage_histogram.png`). No separate `qc/` subfolder.

## notebook.md narrative (per step)

Freeform prose. Recommended sections (not a required template — include what applies):

- **Context** — what this step is for; what the prior step decided
- **What I did** — work performed; scripts run; KG queries issued
- **What I saw** — tables and figure links inline; cited publications from the KG (by DOI / experiment ID — resolved via `list_publications`, never from memory)
- **Surprises** — anomalies, data oddities, unexpected distributions
- **Decisions** — in prose with dates, if any were forced this step
- **Advance rationale** — one line at the end

Figures produced in this step must be linked from `notebook.md`. Publications referenced must be cited with their KG-resolved identifiers (DOI or experiment ID), not free-text names.

## paper.md growth pattern

Single `paper.md` at the analysis root. Skeleton sections exist from day 1 and fill in as steps complete:

| Section | Populated from |
|---|---|
| Question | step 1 |
| Background | step 2 (KG entries and prior work) |
| Methods | steps 3 (framing) and 4 (implementation) |
| Results | step 5 |
| Discussion | step 6 |

When the analysis ends, the paper ends. Avoids the B2 pattern where the write-up was deferred to a final step that never happened.

## Labels

Labels are permitted within a single document when each label is paired with a short readable name in the same paragraph on first mention. No cross-file labels.

Permitted: "the coculture-stronger prediction (2) was not supported; (2) rests on the assumption that axenic cells shut down before fully engaging response"

Not permitted: "P2 failed" (no name; requires a glossary) or "see D7 in decisions.md" (cross-file label reference).

If in practice paired-label-plus-name still obscures content, the convention tightens to prose-only.

## Skill changes

### `skills/research-methodology/SKILL.md` — rewrite in place

**Keep (domain rules):**
- Rule 1: KG is the sole data source
- Rule 2: Locus tags, not gene names
- Rule 3: Source tagging (`[KG]` / `[interpretation]` / `[gap]`)
- Rule 4: Artifacts, not answers
- Rule 5: Scripts over chat reasoning (with interactive-discovery carve-out)
- Rule 6: Statistical rigor
- Rule 7: Don't hallucinate (strengthen: mandatory `ToolSearch` before claims about MCP tool capabilities; mandatory `list_publications` before naming authors)

**Rewrite (workflow section):**
- Rule 8 currently describes "Research notebook, not pipeline" with spec/plan/decisions/hypotheses documents and recipes. Rewrite to describe the 6-step flow, the do-show-explore-decide rhythm, the decide-gate checklist, and the just-in-time formalization principle.
- References to `spec.md`, `plan.md`, `decisions.md`, `hypotheses.md` → replaced by per-step `notebook.md` + `paper.md`.
- Reference to recipes as prerequisites → recipes as optional references (see below).

**Reference docs affected** (under `skills/research-methodology/references/`):
- `step-protocol.md` — rewrite to describe do-show-explore-decide + decide-gate checklist
- `research-notebook.md` — rewrite for the new per-step notebook.md structure + paper.md growth
- `artifacts.md` — update directory layout to match new structure
- Others (`kg-rules.md`, `gene-identity.md`, `anti-hallucination.md`, `python-api-guide.md`, `statistical-rigor.md`) — light edits only, if any

### Recipe skills — kept as optional references

`skills/recipes/enrichment/`, `response-matrix/`, `conservation/` remain in the repo. They are **not** invoked automatically and are **not** prerequisites for analyses. If step 4 methods for a new analysis happens to match an existing recipe pattern, the researcher may point Claude at the recipe as a starting-point reference.

Any text in the `research-methodology` skill stating that recipes are loaded after the methodology skill must be removed or rephrased as optional.

### `CLAUDE.md` — narrow edit

Replace the "Research methodology" section (currently describes loading the research-methodology skill and invoking recipe skills via superpowers workflow) and the "Process overrides" section (currently lists notebook-commit gates, subagent reviews, interactive discovery steps) with a concise description of:

- The 6-step flow
- The do-show-explore-decide intra-step rhythm
- The decide-gate + researcher approval
- The just-in-time formalization principle

Leave unchanged: Project Overview, Architecture, Dependencies, Plugin Structure, Naming conventions, Evaluation Framework, Logs, Key Files.

## Out of scope

- **Retrofitting B2** to the new flow. B2's data is salvageable — someone could write the missing step 4 decide + step 5 — but that's a follow-up effort, not part of this redesign.
- **Evals and hooks** (`evals/`, `hooks/`, `scripts/analyze_usage.py`, `scripts/run_*.py`). Unaffected by this change.
- **MCP tool development** in the sibling `multiomics_explorer` repo.

## Open questions for the first analysis under the new flow

These are intentionally left for the first real analysis to answer, rather than pre-decided:

- Does "label + name in same paragraph" read well in practice, or does it need to tighten to prose-only?
- Does the rewritten `research-methodology` skill stay short, or does it drift long when the workflow section is added back in a new form?
- Does any recipe skill get used often enough that "optional reference" should be promoted back to "default load"?

## Success criteria

The redesign works if, one analysis from now:

- The researcher can read the analysis folder end-to-end without needing a glossary of symbolic labels
- Each step's `notebook.md` is written and committed before the next step begins (the B2 failure mode where steps 3 and 4 had no notebook entries is prevented)
- The final step produces a completed `paper.md` (the B2 failure mode where step 5 never started is prevented)
- The workflow is small enough to hold in working memory without a 100KB plan document
