# Research notebook

## Step protocol and enforcement

For per-step commit timing, hard gates, and the decide-gate checklist, see [Step protocol](step-protocol.md). This document owns the notebook **format and content**; step-protocol owns **when things happen and what gates enforce them**.

---

Every analysis is driven by **per-step interactive notebooks** — one `notebook.md` per step folder — plus a single `paper.md` at the analysis root that grows as a polished write-up. Implementation can be fast or delegated; quality control, exploration, and decision-making are always interactive.

## Just-in-time formalization

Plans, framings, predictions, metrics, terms, decisions, and caveats enter the analysis **only when the data demands them**. Nothing is enumerated in advance "just in case."

Look at the data before drafting the plan. Pull what's in the KG — counts, fields, coverage — then propose the minimal plan that fits. Start simple; expand only when a specific finding forces it.

Concrete rules:
- **Step 1 grounds the dialogue in MCP queries, not assumptions.** Before locking scope or proposing sub-questions, run `list_publications`, `list_experiments`, `list_organisms` filtered to the prompt's context. Surface counts and structural surprises (e.g. "axenic RNA-seq is single-contrast, not time-course") in the dialogue. Capture the queries and key counts in step-1 `notebook.md` under a **"KG context"** sub-section. Step 1 grounds; step 2 enumerates.
- **Step 1 locks the question, not sub-questions.** Sub-narratives that form during dialogue (e.g., "protein persists while mRNA is gone") belong to step 3. Defer; flag as step-3 candidate.
- **Step 3 framing has a floor, not a template.** Minimum: hypothesis in prose, what success means operationally, positive and negative controls from the KG. Nothing else required.
- **Preregistration is optional and minimal.** If confirmation bias is a real risk, preregister 1–3 named predictions — not a 4×4 matrix of ordering and thresholds.
- **Stability checks are added when a specific result triggers them**, not planned up front.
- **Decisions are written when forced by data**, not anticipated. No "we may decide X" placeholders.
- **Methods (step 4) stay minimal** — ad-hoc Python module with exactly what this analysis needs.
- **Caveats are harvested at step 6** from what actually happened; not pre-cataloged.

This principle governs every step. If you find yourself listing or architecting things the analysis might need before the data has arrived, stop.

## The 6-step flow

1. **Research question** — conversation: user prompt + Claude clarifying questions → locked question
2. **KG entries** — identify relevant publications, experiments, organisms, data types
3. **Analysis framing** — (a) selection: publications/experiments/organisms/data types; (b) framing: hypothesis, target, positive and negative controls, expected outcome — all in KG terms
4. **Methods** — pick one item from step 3 as a driving example; select or generate an analysis method; produce an ad-hoc Python module
5. **Analyze** — run the method; produce scored outputs, figures, tables
6. **Evaluate** — assess results against framing; harvest caveats; finalize paper

Steps 1–3 are the **research proposal**. Locked at end of step 3. Steps 4–6 execute against it.

Step 1 produces `notebook.md` (including a "KG context" sub-section that captures grounding queries, counts, and any structural surprises encountered during the dialogue). Step 1 is interactive, not scripted, but it is not assumption-driven — MCP queries ground the conversation. Steps 2–6 add `scripts/`, `data/`, `figures/`, and QC alongside `notebook.md`.

## The intra-step rhythm: do → show → explore → decide

Every step advances through these four phases (see [step-protocol.md](step-protocol.md) for commit timing and gates):

- **do** — do the step's work; outputs land wherever the step naturally produces them
- **show** — populate `notebook.md` with what was produced
- **explore** — investigate anomalies, surprises, or gaps; add `qc_*.py` checks, sensitivity analyses, or follow-up clarifying questions as needed
- **decide** — finalize notebook, update paper.md, present state to researcher, commit on approval

## notebook.md format

One `notebook.md` per step folder. Freeform prose, not a rigid template — include what applies to the step.

### Recommended sections

- **Context** — what this step is for; what the prior step decided
- **What I did** — work performed; scripts run with their command-line invocation for non-trivial cases; KG queries issued
- **Results** — summary tables shown inline (as markdown tables, not prose paraphrases); links to full tables in `data/` and figures in `figures/` produced this step; cited publications from the KG by DOI or experiment ID — resolved via `list_publications`, never from memory (see [anti-hallucination.md — Category 5](anti-hallucination.md))
- **Surprises** — anomalies, data oddities, unexpected distributions worth flagging
- **Decisions** — in prose with dates, if any forks were taken this step; omit if none
- **Advance rationale** — one line at the end

### Decide-gate checklist (end of notebook.md)

At step close, the notebook ends with this checklist (see [step-protocol.md](step-protocol.md) for the approval gate):

- **Outputs produced** — filenames in `scripts/`, `data/`, `figures/`, with command lines for non-trivial scripts (for reproducibility)
- **Results presented** — summary tables shown inline in `notebook.md`; links to full tables and figures generated this step
- **QC gate** — what was checked → result (one line per check)
- **Decisions made this step** — prose + date, if any; omit the section if none
- **Advance rationale** — one line, why this step is ready to close

The checklist must stay this minimal. It is not a template to extend with optional fields. Inflation of this list reintroduces the premature-formalization failure that the per-step notebook is designed to prevent.

### Labels

Labels are permitted within a single document when each label is paired with a short readable name in the same paragraph on first mention. **No cross-file labels.**

- OK: "the coculture-stronger prediction (2) was not supported; (2) rests on the assumption that axenic cells shut down before fully engaging response"
- Not OK: "P2 failed" (no name) or "see D7 in decisions.md" (cross-file label)

If paired-label-plus-name still obscures content, drop to prose-only — labels are a convenience, not a requirement.

### Overwrite vs append

A step's `notebook.md` represents what-we-now-believe-happened in that step. On redo, the notebook is **overwritten**, not appended — the narrative reflects the successful attempt, not a log of past attempts. The prior attempt lives in git history (see [step-protocol.md — Redo path](step-protocol.md)).

`gaps_and_friction.md` is **append-only** (see below).

## paper.md growth pattern

Single `paper.md` at the analysis root. Skeleton sections exist from day 1 (seeded by Claude during scaffolding at start of step 1) and fill in during each step's **decide** phase — after the step's notebook is finalized but before the commit.

| paper.md section | Populated from |
|---|---|
| Question | step 1 |
| Background | step 2 (KG entries and prior work) |
| Methods | steps 3 (framing) and 4 (implementation) |
| Results | step 5 |
| Discussion | step 6 |
| References | accumulates across all steps that cite publications |

References are populated as publications are cited. Every reference must be resolved through `list_publications` and cited by DOI or KG experiment ID — never drafted from intrinsic knowledge (see [anti-hallucination.md — Category 5.2](anti-hallucination.md)). Citation format inside prose can be short (author-year or numeric); the References section at the end carries the resolved DOI or experiment ID for each.

When the analysis ends, the paper ends. The B2 failure mode was deferring write-up to a final step that never happened — `paper.md`'s incremental growth prevents this.

## gaps_and_friction.md (transitional)

A top-level file at `analyses/<slug>/gaps_and_friction.md` captures friction encountered during the analysis, distinct from decisions:

- **Decision** = a fork the analysis had to take, based on data → logged in the relevant step's `notebook.md`
- **Friction** = a problem that slowed us down, surprised us, or revealed a gap in methodology / KG / tooling → logged in `gaps_and_friction.md`

These can co-occur — a methodology gap can force both a decision and a friction entry. Log in both places when that happens.

**What goes in `gaps_and_friction.md`:**
- KG data issues and bugs encountered
- MCP tool schema or capability mismatches (see [anti-hallucination.md — Category 5.1](anti-hallucination.md))
- Methodology gaps discovered during execution (nuances the framing didn't anticipate)
- Anti-hallucination corrections (claims from memory caught by verification)
- Process friction (things that slowed the work)

Each entry is prose with a date, a short name, what happened, and — if relevant — the workaround and the downstream impact on methodology/KG/tooling.

**Why this is transitional.** The 6-step methodology itself is under development. Every analysis teaches us something about what's missing or awkward. `gaps_and_friction.md` is the learning record that feeds back into methodology and KG/tooling improvements. Mandatory while the methodology is being stabilized — likely the first 3–5 analyses. Once the pattern settles, revisit: keep, make optional, or fold into `notebook.md`. Retirement criterion: when two consecutive analyses produce a near-empty `gaps_and_friction.md` (only incidental friction, no methodology gaps), propose retiring the requirement.

Append-only. Redo friction entries accumulate.

## Using `superpowers:brainstorming` for step 1

The brainstorming skill's dialogue pattern (clarifying questions one at a time, proposing approaches, converging) fits step 1 well. **Three overrides apply:**

1. **Capture location.** The skill's default writes a design doc to `docs/superpowers/specs/YYYY-MM-DD-<topic>-design.md`. Override: output lands in `analyses/<slug>/1_question/notebook.md`. The formulated research question, clarifying dialogue in summary form, rejected alternatives, and converged scope all live there.
2. **KG grounding alongside the dialogue.** Brainstorming is pure conversation by default. Step 1 here adds lightweight MCP grounding (`list_publications`, `list_experiments`, `list_organisms` filtered to the prompt's context) before scope or sub-question decisions. Capture queries and counts in a "KG context" sub-section in `notebook.md`. See the Just-in-time concrete rules above for the rationale.
3. **Terminal action.** The skill's terminal state is invoking `superpowers:writing-plans`. **Skip this.** Step 1 decide advances to step 2 (KG entries), not to implementation-plan writing. There is no monolithic plan for the analysis — the 6-step flow replaces it.

If the skill's preamble nudges toward writing a spec doc outside the analysis folder or invoking writing-plans at the end, treat those as defaults being overridden by this methodology skill (which has higher priority for research work in this repo).

## QC checkpoint types

What to show depends on the step type.

### KG selection (step 2)
- Row counts per filter: experiments at each stage (started with → filter by organism → filter by assay → final)
- Sample rows of selected entries (experiment ID, publication, TPs, omics)
- Per-TP gene counts (`tp_gene_count`), **not** cumulative `gene_count` (see [anti-hallucination.md — Category 5.3](anti-hallucination.md))
- Publication attributions resolved via `list_publications`

### Analysis framing (step 3)
- Controls selected from the KG with validation: what distinguishes positive from negative; coverage of TPs/conditions; distributional QC
- Hypothesis statement in prose
- Expected outcome phrased in KG-operational terms (what table / metric / direction will change)

### Computation / metric (step 4)
- Worked example: 2–3 genes or clusters through the formula with actual numbers, step by step
- Summary statistics of the output (distribution, range, NaNs)
- Sanity check against known biology ("glnA should score high — does it?")

### Scoring / comparison (step 5)
- Full results table in markdown (the actual numbers, not prose summary)
- Best/worst scores, surprises, anything unexpected
- Cross-condition comparison with expectation check against the step 3 framing

### Evaluation (step 6)
- Preregistered predictions held / not held
- Sensitivity / LOO stability where applicable
- Harvested caveats

## Code lifecycle: analysis-first, productize later

Research code has two phases. Follow phase 1 during analysis; flag phase 2 candidates.

### Phase 1 — Analysis code (methodology-first)

Code lives in the analysis directory. Goal: correct methodology, not good software engineering.

**Separate reusable logic from scripts.** When step 4 introduces a new method, put the reusable logic in a utility package within the step folder (`4_methods/<module_name>.py`) and let analysis scripts in later steps import it. Scripts call utilities for specific data; utilities contain the methodology. This separation is what makes toy-testing possible and productization straightforward later.

**Methods modules describe methodology, not implementation scaffolding.** For novel utilities (scoring functions, metrics, gene-set operations), the module should be minimal — the formula with a worked example, expected I/O, and the minimum code to compute it. Regular extraction/plotting scripts are straightforward enough not to need a separate utility.

**Toy-data verification before real data.** When building a reusable utility, verify with hand-calculated toy examples first. Create small synthetic input, compute expected output by hand, run the utility, compare, log the verification in `4_methods/notebook.md`. Applies to anything in a shared utility; one-off scripts don't need it.

**Refine through the notebook QC cycle.** The do → show → explore → decide loop is how methodology gets validated. Formula corrections, edge cases, direction logic — all discovered through the researcher walking through concrete examples.

### Phase 2 — Productization (software-first)

After the analysis, if a utility proves reusable (used across multiple analyses), flag it for productization — a separate brainstorm with API design, tests, and documentation. It moves from the analysis directory to a shared package (e.g., `multiomics_explorer/analysis/`).

Don't productize speculatively. Wait for proven reuse — the analysis notebooks across multiple analyses are the evidence.
