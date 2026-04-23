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

1. **Research question** — the user states the question; Claude asks clarifying questions; they converge on a formulated question locked in `notebook.md`. This step is a conversation, not a computation. `superpowers:brainstorming` is the natural tool.
2. **KG entries** — identify relevant publications, experiments, organisms, data types in the KG
3. **Analysis framing** — (a) selection: publications, experiments, organisms, data types (DE, cluster, ontologies); (b) framing: hypothesis, target, positive and negative controls, expected outcome — all in KG terms
4. **Methods** — pick one item from step 3 as a driving example; select or generate an analysis method (statistical test, score, ...); produce an ad-hoc Python module
5. **Analyze** — run the method; produce scored outputs, figures, tables
6. **Evaluate** — assess results against framing; write up caveats; finalize paper

Steps 1–3 are collectively the **research proposal**. Locked at the end of step 3 decide. Steps 4–6 execute against it.

Step 1 is the only conversation-only step — it converges on a formulated question through dialogue and typically produces only `notebook.md`. Every other step (2–6) involves computation: step 2 queries the KG and filters entries; step 3 validates selection and control choices with QC scripts; step 4 builds the method module; step 5 runs the analysis; step 6 performs sensitivity / evaluation analyses against the framing. Steps 2–6 all produce scripts, data, figures, and QC alongside `notebook.md`.

### Using `superpowers:brainstorming` for step 1

The brainstorming skill's dialogue pattern (clarifying questions one at a time, proposing approaches, converging) fits step 1 well. Two overrides apply:

- **Capture location.** The skill's default is to write a design doc to `docs/superpowers/specs/YYYY-MM-DD-<topic>-design.md`. Override this — the output lands in `analyses/<slug>/1_question/notebook.md` instead. The formulated research question, the clarifying dialogue (in summary form), any rejected alternatives, and the converged scope all live there.
- **Terminal action.** The skill's terminal state is to invoke `superpowers:writing-plans`. **Skip this.** Step 1 decide advances to step 2 (KG entries), not to implementation-plan writing. There is no monolithic plan for the analysis — the 6-step flow replaces it.

If the skill's preamble nudges toward writing a spec doc outside the analysis folder or invoking writing-plans at the end, treat those as defaults being overridden by the research methodology skill (which has higher priority in this repo).

## Intra-step rhythm: do → show → explore → decide

Every step — including step 1 — advances through these four phases:

- **do** — do the step's work; outputs land wherever the step naturally produces them (a conversation lands in `notebook.md`; scripts land in `scripts/` with outputs in `data/` and `figures/`)
- **show** — populate `notebook.md` with what was produced: the question (step 1), the KG entries (step 2), the framing and selection QC (step 3), the method module and test outputs (step 4), the analysis tables and figures (step 5), the evaluation results (step 6)
- **explore** — investigate anomalies, surprises, or gaps: ask follow-up clarifying questions (step 1); add `qc_*.py` checks, run sensitivity analyses, cross-validate against controls, produce exploratory figures (steps 2–6)
- **decide** — finalize notebook.md; update the relevant section of paper.md with this step's synthesis; minimal checklist (below); **researcher approval gate**; git commit; advance

The decide phase is the only formal gate between steps.

### Decide-gate checklist (radically minimal)

Each step's decide phase must produce, in its `notebook.md`:

- **Outputs produced** — filenames in `scripts/`, `data/`, `figures/`, with command lines for non-trivial scripts (for reproducibility)
- **Results presented** — summary tables shown inline in `notebook.md`; links to full tables and figures generated this step
- **QC gate** — what was checked → result (one line per check)
- **Decisions made this step** — prose + date, if any; omit the section if none
- **Advance rationale** — one line, why this step is ready to close

Claude then pauses, displays the state, waits for researcher approval, commits the step, and begins the next.

The checklist must stay this minimal. It is not a template to extend with optional fields. Inflation of this list would reintroduce the same premature-formalization failure the redesign exists to prevent.

## Directory structure per analysis

```
analyses/YYYY-MM-DD-<slug>/
  paper.md                        # single academic-style writeup, grows from step 1
  gaps_and_friction.md            # transitional — methodology/KG/tooling friction log
  1_question/
    notebook.md
  2_kg_selection/
    notebook.md, scripts/, data/, figures/
  3_analysis_framing/
    notebook.md, scripts/, data/, figures/
  4_methods/
    notebook.md
    <module_name>.py              # ad-hoc method module
    scripts/, data/, figures/
  5_analyze/
    notebook.md, scripts/, data/, figures/
  6_evaluate/
    notebook.md, scripts/, data/, figures/
```

**QC artifacts** live in the same `scripts/` / `data/` / `figures/` folders as main outputs, using a `qc_` filename prefix (e.g., `scripts/qc_check_tp_coverage.py`, `figures/qc_tp_coverage_histogram.png`). No separate `qc/` subfolder.

### Scaffold creation

Claude creates the scaffold at the **start of step 1**, before the brainstorming dialogue begins. The researcher does not pre-create folders.

1. **Propose a slug.**
   - **Who:** Claude proposes; the user approves or counter-proposes.
   - **When:** at the start of step 1, before any folder or file is written. No scaffold is created until the slug is confirmed.
   - **Source:** the user's initial prompt / topic.
   - **Format:** `YYYY-MM-DD-<descriptor>` where descriptor is snake_case, 2–4 words reflecting the research question (e.g., `2026-04-23-n_limitation_signature`). Same-day collision: append a letter suffix (`-a`, `-b`).
2. **Create the minimal scaffold:**
```
analyses/YYYY-MM-DD-<slug>/
  paper.md              # skeleton with empty section headers: Question, Background, Methods, Results, Discussion, References
  gaps_and_friction.md  # header only; entries added as friction occurs
  1_question/
    notebook.md         # empty; populated by the step 1 dialogue
```
3. **Begin step 1 dialogue.** `1_question/notebook.md` grows through the clarifying exchange.
4. **At step 1 decide:** `paper.md`'s Question section is populated from the locked research question; the first git commit includes the scaffold plus step 1 artifacts.

Step folders for 2–6 are created **progressively** — each folder (and its `scripts/` / `data/` / `figures/` subfolders as needed) is created when that step begins, not pre-created at scaffold time. This avoids empty folders cluttering the analysis during work in progress.

## notebook.md narrative (per step)

Freeform prose. Recommended sections (not a required template — include what applies):

- **Context** — what this step is for; what the prior step decided
- **What I did** — work performed; scripts run with their command-line invocation for non-trivial cases; KG queries issued
- **Results** — summary tables shown inline; links to full tables in `data/` and figures in `figures/` produced this step; cited publications from the KG (by DOI or experiment ID — resolved via `list_publications`, never from memory)
- **Surprises** — anomalies, data oddities, unexpected distributions
- **Decisions** — in prose with dates, if any were forced this step
- **Advance rationale** — one line at the end

Figures produced in this step must be linked from `notebook.md`. Publications referenced must be cited with their KG-resolved identifiers (DOI or experiment ID), not free-text names.

## paper.md growth pattern

Single `paper.md` at the analysis root. Skeleton sections exist from day 1 and fill in during each step's **decide** phase — after the step's notebook is finalized but before the commit. Each decide produces a synthesis paragraph (or figure inclusion, or methods sub-section) for the relevant part of paper.md:

| Section | Populated from |
|---|---|
| Question | step 1 |
| Background | step 2 (KG entries and prior work) |
| Methods | steps 3 (framing) and 4 (implementation) |
| Results | step 5 |
| Discussion | step 6 |
| References | accumulates across all steps that cite publications |

References are populated as publications are cited in any step's `notebook.md` or in `paper.md` sections. Every reference must be resolved through `list_publications` and cited by DOI or KG experiment ID — never drafted from intrinsic knowledge (per Rule 7, anti-hallucination). Citation format inside paper prose can be short (author-year or numeric); the References section at the end carries the resolved DOI / experiment ID for each.

When the analysis ends, the paper ends. Avoids the B2 pattern where the write-up was deferred to a final step that never happened.

## gaps_and_friction.md (transitional)

A top-level file at `analyses/<slug>/gaps_and_friction.md` captures friction encountered during the analysis, distinct from the decisions that land in per-step `notebook.md`:

- **Decision** = a fork the analysis had to take, based on data (logged in the relevant step's `notebook.md`)
- **Friction** = a problem that slowed us down, surprised us, or revealed a gap in methodology / KG / tooling (logged in `gaps_and_friction.md`)

These can co-occur — a methodology gap can force both a decision and a friction entry. Log in both places when that happens.

**What goes in `gaps_and_friction.md`:**

- KG data issues and bugs encountered (e.g., B2's background-collapse bug)
- MCP tool schema or capability mismatches
- Methodology gaps discovered during execution (nuances the framing didn't anticipate)
- Anti-hallucination corrections (claims from memory caught by verification)
- Process friction (things that slowed the work)

Each entry is prose with a date, a short name, what happened, and — if relevant — the workaround and the downstream impact on methodology/KG/tooling.

**Why this is transitional.** The 6-step methodology itself is under development. Every analysis we run teaches us something about what's missing or awkward. `gaps_and_friction.md` is the learning record that feeds back into methodology and KG/tooling improvements. Mandatory while the methodology is being stabilized — likely the first 3–5 analyses. Once the pattern settles, we revisit: keep it, make it optional, or fold into `notebook.md`. Decision criterion: when two consecutive analyses produce a near-empty `gaps_and_friction.md` (only incidental friction, no methodology gaps), propose retiring the requirement.

The decide-gate checklist gains an implicit fifth item: if friction was encountered this step, it's appended to `gaps_and_friction.md` before the commit.

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
