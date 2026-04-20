# Research-methodology skill — v3 improvements proposed from B2 brainstorm

**Date:** 2026-04-18
**Source:** Brainstorming session for [pathway-enrichment-b2 design](2026-04-18-pathway-enrichment-b2-design.md)
**Status:** Living document — accretes findings during B2 execution; consolidates into skill v3 when B2 closes.

---

## 1. Purpose

Capture skill-improvement proposals that surfaced during B2 brainstorming, so future brainstorms don't need the same reminders. Two categories:

1. **Content gaps** — missing rules or guidance on specific analytical patterns (enrichment scoring, signature derivation, etc.).
2. **Meta-process gaps** — structural/framing issues in the skill that let the assistant drift from interactive-research pace into pipeline thinking.

This doc is the input to the next skill revision. It lives here rather than in B2's `gaps_and_friction.md` because:
- It aggregates pre-execution findings (before B2's analysis dir exists).
- It's a skill-change proposal, not an analysis-local friction log.
- During B2 execution, `gaps_and_friction.md` will accrete new friction that gets merged back here at the end.
- Prior precedent: `2026-04-06-research-skill-restructure-design.md`, `2026-04-08-research-methodology-v2-improvements-design.md`, `2026-04-07-research-notebook-discipline-design.md`.

## 2. Content gaps — concrete rules to add

Each item is a place where the researcher had to explicitly prompt the assistant during brainstorming for a methodological detail that should have come from the skill.

### Statistical rigor (owner: `statistical-rigor.md`)

**C1. Cap signed_score magnitude in alignment scoring.** Uncapped `−log10(padj)` lets one extreme pathway (padj → 0) dominate a mean across ~10–20 signature pathways. Propose: any spec defining a score of the form `mean(f(signed_score))` must specify a cap (`SCORE_CAP = 10` / padj = 1e-10 is reasonable). Document the rationale: beyond the cap, magnitude differences don't carry useful information for alignment, and padj underflow to 0 produces NaN/inf otherwise.

**C2. Split multi-experiment enrichment by `table_scope`.** Experiments with restricted `table_scope` (`significant_only`, `top_n`, `significant_any_timepoint`, `filtered_subset`) cannot use `background="table_scope"` — the DE table equals or contains the foreground, collapsing the Fisher 2×2. They must use `background="organism"` as a separate pathway_enrichment call; results concatenated with a `background_used` column. The skill mentions "only significant genes as background" as a pitfall but doesn't prescribe the two-call split. Add worked example covering mixed-table_scope experiment sets.

**C3. Bidirectional-pathway tie-breaker for signature derivation.** When a pathway is enriched both up (≥N clusters) and down (≥N clusters) in the reference set, current skill has no guidance. Propose: exclude from signature, log to a separate "ambiguous/dropped" file for transparency. A pathway without clean direction signal across references can't anchor alignment scoring.

**C4. Temporal filter for time-course signature derivation.** B1 observed that pathways can flip direction between early timepoints (<3h) and settled timepoints. Current skill has no guidance on this. Propose: when deriving a signature from a time-course reference, exclude clusters below a biologically-justified early threshold (default: `timepoint_hours < 3` for MED4-scale organisms). Early clusters stay in analysis — just not in signature derivation.

**C5. Per-experiment dominance check for signature robustness.** A signature pathway supported by ≥N clusters all from one experiment is fragile. The skill currently doesn't prescribe computing per-experiment contribution counts. Propose: signature CSVs must include a `per_experiment_breakdown` column; show/explore phase must surface single-experiment dominance as a flag.

**C6. Cross-background calibration.** When scores come from mixed backgrounds (C2), NC floor and PC range must be computed per `(ontology, background_used)` group, not pooled. Skill currently doesn't address this. Propose: when a spec uses `background="organism"` for some clusters and `background="table_scope"` for others, require per-group calibration and flag any T cluster whose matched-background calibration set is empty.

**C7. Dropped-but-notable pathway log.** Not every excluded pathway is equally uninteresting. Bidirectional pathways with high enrichment, or below-threshold pathways with extreme magnitude in 1–2 clusters, represent interesting dropped biology. Propose: add a "transparency of exclusion" rule — any pathway with `max |signed_score| ≥ 3` (padj ≤ 1e-3) in any reference cluster gets logged to a dropped-notable file even if it doesn't make the signature.

**C8. Edge-case checklist for formula-based derivations.** For every score / filter / derivation step defined in a spec, prompt for: magnitude outliers (cap?), direction ties (tie-breaker?), temporal edges (time-course flip?), support concentration (single-experiment dominance?). Make this a section in `statistical-rigor.md` so spec-writers check each box.

### Methodology and scoping (owner: `kg-rules.md` or new `scoping.md`)

**C9. Organism-matched reference signature.** For "is organism A in condition X?" questions, derive the reference signature from organism A's own experiments, not pooled cross-organism. Non-matched organism experiments go in as context only. Current skill lacks this principle; during B2 brainstorm the researcher had to pull me back from a cross-organism pooling approach.

**C10. Experimental-design classification vocabulary.** Target / Reference / Positive Control / Negative Control / Context is a generic pattern that applies across analyses. Propose: add standard vocabulary to the skill (maybe in `kg-rules.md`), with a CSV schema for `experiments_classified.csv` (columns: experiment_id, class, rationale).

**C11. "Open-mind reset" for follow-up analyses.** When doing a redo/follow-up (B2 is a redo of B1), the skill should prompt: "what's changed in the KG, the tool surface, and the methodology since the predecessor?" The researcher had to say "open mind — not trying to replicate B1" to break me out of B1's framing.

### Artifacts and presentation (owner: `artifacts.md`)

**C12. CSV preferred over JSON for LLM-readable caches.** Skill mentions CSV preference but doesn't explicitly frame as "tabular, greppable, token-efficient for Claude Code inspection." Propose: add a subsection "Artifact format convention for Claude Code" that explicitly prescribes CSV via `to_dataframe`, pickle only for stateful API objects, envelope metadata in notebook.

**C13. Figure scope: show context, highlight what matters.** Don't filter a figure to just the pathways in the score; show all enriched pathways and highlight signature members. Principle: researcher needs to see what moved *outside* the score, not just what was scored. Add to `artifacts.md`.

**C14. Time-course visualization convention.** For time-course experiments, lineplot trajectories (panel per experiment) complement heatmaps. Skill has no figure-type-to-experiment-type guidance.

**C19. Figures are iterative; plans/specs should not over-specify pre-run.** Figure polish (fontsize, row density, panel grid, highlight styling, label rotation) can only be judged against real data — spec'ing 200 lines of matplotlib upfront produces false confidence and rot on mismatch. Pattern surfaced in B2 planning: the figures step was the only one where "write full inline template" felt wrong. Proposed rule (for `artifacts.md` or new `figures.md`): separate **structural choices** (axis semantics, row selection, column ordering, value caps — spec-locked) from **cosmetic choices** (sizes, highlights, rotations — marked `# TWEAK:` in the script, iterated in show/explore/decide). Plans should instruct executors to treat first-run output as a draft, not a deliverable. `step-protocol.md` should explicitly allow "figures iteration" as a sub-loop within Step-N show/explore/decide — multiple commits acceptable as `fig 1 iter 1`, `fig 1 iter 2`, each with a notebook entry.

### Python API patterns (owner: `python-api-guide.md`)

**C15. Pickling stateful API result objects.** `EnrichmentResult`-style dataclasses can pickle individually but fail in dict aggregation. Propose: add a section on "persisting stateful API return objects," prescribing (1) round-trip verification before full run, (2) split-to-subdir fallback, (3) constituent-state reconstruction, (4) inline re-run as last resort.

**C16. Prefer `search_text` fuzzy catch for experiment discovery.** The skill's KG-rules examples use explicit filter values (`treatment_type=["nitrogen_stress"]`) which requires pre-discovering the tag vocabulary. `search_text="nitrogen"` is often cleaner — Lucene full-text match across name/treatment/control/context fields. Add as idiomatic discovery pattern.

### Process — step-protocol enforcement (owner: `step-protocol.md`)

**C17. Plan must materialize each step-protocol phase as a discrete task.** The step-protocol describes the cycle but doesn't instruct plan-writers to create separate tasks for do/show/explore/decide. Without this, a plan-agent collapses the cycle into one "run Step 2 and report" task. Propose: add a "Plan-authoring rules" section that requires plans to list show/explore/decide as explicit tasks, and require spec-compliance review to check for collapsed phases.

**C18. Notebook entry is a per-step deliverable, explicit in every spec.** Skill references notebook writing in `research-notebook.md` but specs can reference it generically ("follows the cycle") rather than explicitly list "write notebook entry" per step. Propose: add a spec-template rule that every step's description explicitly mentions the notebook entry as a deliverable of show/explore/decide.

---

## 3. Meta-process gaps — framing issues that cause drift

These aren't content; they're structural weaknesses in how the skill is framed or referenced, causing the assistant to drift from interactive-research pace into pipeline thinking. Each was observed multiple times in this brainstorm.

**M-A. Drift toward "pipeline" framing in first-pass specs.** First-pass spec structure is sequential data-producing phases; interactive gates get added only after prompting. Root cause: `research-notebook.md` describes the cycle but doesn't frame a step as "interactive walkthrough supported by a script" — which is the actual intent. Proposed change: rewrite the "what a step is" subsection to lead with "a step is an interactive walkthrough; a script is what makes it reproducible."

**M-B. Parallel-terminology invention.** I invented "QC phase" when show/explore/decide already covered that function. Loaded skill didn't prevent this. Proposed change: add an explicit warning in `step-protocol.md` — "do not introduce phase vocabulary outside do/show/explore/decide. If you find yourself defining 'QC phase' / 'review step' / 'validation pass,' you are reinventing show+explore."

**M-C. Under-emphasis of notebook writing per step in specs.** Implicit reference to the skill isn't enough; specs tend to drop explicit notebook obligations. Proposed change: `artifacts.md` / spec-template guidance should require explicit "Notebook entry: [what goes in]" bullet per step.

**M-D. Plan-level pace enforcement is not automatic after spec.** Even with the skill loaded, the spec → plan → execution handoff loses the interactive pace unless the plan explicitly enforces it. Proposed change: add a "Plan-authoring rules" section (see C17); brainstorming skill should invoke this when handing off to writing-plans.

**M-E. Anchoring on predecessor analysis.** When a brainstorm is framed as "redo of X," I tunnel on X's choices. Proposed change: add a "follow-up analysis reset" prompt in `kg-rules.md` — before adopting any of the predecessor's choices, explicitly check what's new in the KG, tools, and skill.

**M-F. Script-first default rather than interactive-first.** I treat "do" as the main step and show/explore/decide as adjuncts. Skill's intent is the inverse. Proposed change: `step-protocol.md` should emphasize that the *interactive walkthrough is the primary deliverable* — the researcher learning from the data — and the script is the reproducibility backstop.

**M-G. Batch-present drift in brainstorming.** I push toward presenting whole designs in chunks even though the brainstorming skill says "one question at a time." Proposed change: the brainstorming skill should add a red-flag self-check: "am I about to present more than one decision that could be separately iterated?"

**M-H. Tool-doc re-read reminder missing.** I added a BRITE tree iteration because I didn't re-read the tool docs carefully; prior knowledge was stale. Proposed change: `kg-rules.md` should include — "when using a tool in a spec or script, re-read `docs://tools/<name>` even if familiar. Prior knowledge rots."

**M-I. Pre-registration not surfaced as structured.** I initially treated pre-registration as free text rather than a structured object with dimensions (condition × background × expected score). Proposed change: `statistical-rigor.md` should prescribe pre-registration as a named artifact with required fields.

**M-J. Edge-case blind spots for formulas.** For every formula-based step, I needed a prompt for the edge cases. Proposed change: formalize as "edge-case checklist" in `statistical-rigor.md` — see C8.

**M-K. Analysis classification vocab not inherited.** T/R/PC/NC/CTX is reusable; had to specify manually. Proposed change: see C10.

---

## 4. B2 plan-review conclusions (added 2026-04-20)

Before execution began, the B2 plan went through a multi-pass review that surfaced 20+ issues in the inline Python / bash templates. The issues cluster into a small number of patterns that aren't spec-specific — they're hygiene failures common to any research-extraction codebase. Capturing them here so (a) the plan-authoring skill can get a hygiene checklist, (b) future plan reviews don't need to rediscover the patterns, (c) the research-methodology skill v3 can promote the reliable ones to hard rules.

### 4.1 Four recurring failure modes in research-code templates

These are the "crooked code" patterns the plan review caught repeatedly. Propose adding a dedicated `research-code-hygiene.md` (or a subsection of `python-api-guide.md`) that enumerates them with corrections.

**F1. Heuristic string parsing on structured IDs.** Example caught: `organisms = {row["experiment_id"].split("_")[-1] for row in CLASSIFICATIONS}` — using the suffix of a DOI-derived ID to infer organism. Failure mode: IDs are opaque keys; any token-position heuristic is silently wrong when the ID schema evolves. **Rule:** carry every field you need as an explicit column. If you find yourself parsing an ID string to recover a field, fix the upstream data structure instead.

**F2. Silent-last-wins on multi-row frames.** Example caught: `class_map = dict(zip(classified["experiment_id"], classified["class"]))` applied to a DataFrame with one row per (experiment × timepoint). `dict(zip(...))` keeps only the last value per key; if class varies (even by typo) across timepoints, the lookup is silently wrong. Same pattern: `exp_ids = df["experiment_id"].tolist()` when `df` has repeated IDs → passes duplicates to APIs. **Rule:** dedupe to the natural primary key (`df.drop_duplicates("experiment_id")`) BEFORE casting a DataFrame into a scalar lookup, a list, or a set intended to describe experiments.

**F3. Hardcoded condition-specific strings as equality checks.** Example caught: `classified["organism_name"] == "Prochlorococcus MED4"` — breaks silently if the KG display name drifts, or if an upstream call returns `"MED4"` (substring form). **Rule:** match by case-insensitive substring (`.str.contains(MED4, case=False)`), by validated enum drawn from `list_organisms()`, or by a parameter loaded from config — never by a hardcoded full-name equality unless validated non-empty immediately after.

**F4. Empty-intermediate blindness.** Every research step has a DataFrame that could be empty: filter removed everything, enrichment returned zero significant rows, calibration has no NCs in this `(ontology, background)` group. Without a guard, the pipeline runs to completion with an empty output and the failure surfaces three tasks later at a downstream `.iloc[0]`. **Rule:** every `groupby` / `concat` / `read_csv` / post-filter frame gets either (a) an explicit empty-guard (`if df.empty: log.error(...); return 1`), or (b) a `log.info(f"rows={len(df)}")` sanity trace. Prefer (a) for anything that would break downstream; (b) is acceptable only for inspection.

### 4.2 Empirical probe before elaborate contingency

Spec §8 risk 7 enumerated four escalation paths for `EnrichmentResult` pickle failure (split-to-subdir, `dill`/`cloudpickle`, constituent-state reconstruction, inline re-run). During review, a 10-minute probe against the live KG confirmed standard pickle round-trips cleanly at B2 scale. **The four escalation paths remain valid as precautionary fallbacks, but they never would have been needed.**

**Rule for specs and plans:** when a risk enumerates ≥3 escalation paths, the plan must include an empirical probe before committing to contingency planning. Probes are typically <15 minutes to write and <5 minutes to run, and they either rule out the risk entirely or surface a concrete failure mode that sharpens the fallback choice. Speculative contingency planning without a probe is a distractor.

Propose: add this rule to `step-protocol.md` as a sub-section "Before contingency planning, probe." Also see C15 (pickling-specific) — this is the generalization.

### 4.3 Spec → plan parameterization drift

Example caught: spec §5.0 hardcoded `exploration/2026-04-18-notebook.md` (the brainstorm date). The plan was executed on 2026-04-19, copied the spec verbatim into the analysis dir via `cp`, and the copied spec then referenced a file that wasn't there. Same pattern: organism names, ontology keys, experiment IDs, anything date-stamped.

**Rule:** specs reference parameters by role (`<target_organism>`, `<notebook_filename>`), not by concrete values that bind at spec-write time. The plan (or a preamble in the spec) maps parameters → concrete values. If the spec must use concrete values (e.g., "MED4" because the entire analysis is MED4-specific), call them out in a "Concrete parameters" subsection that's easy to find and update.

Propose: add a spec-template rule to `artifacts.md` — specs declare their parameters explicitly; plans resolve them. When a spec is copied into an analysis dir, the parameter resolution travels with it.

### 4.4 Figures are iterative within a step

Already captured in C19 — but reinforced by the review. The figures script was the only Step where "write the full inline code upfront" felt wrong. The resolved approach (skeleton + `# TWEAK:` markers + explicit iteration protocol in show/explore/decide) is the right pattern for any visualization step. Do not over-specify figure cosmetics in plans; do specify structural choices (axis semantics, row selection, color cap).

### 4.5 multiomics_explorer API surface gaps

Confirmed during review; flagged for the explorer team.

**A1. `list_experiments` lacks `experiment_ids` filter.** All other "fetch by experiment" tools take one (`pathway_enrichment`, `ontology_landscape`, `gene_overview`, `gene_response_profile`, `differential_expression_by_gene`). For the Task-1 pattern "classify N experiments, fetch metadata for those N," the workaround is pull unfiltered + local filter. Cheap today (~76 experiments), wasteful as KG grows. Propose: add `experiment_ids: list[str] | None = None` parameter mirroring the filter shape of the sibling tools.

**A2. `pathway_enrichment` "Package import equivalent" docs block is wrong.** The MCP tool doc says "returns dict with keys..." but the Python `pathway_enrichment` actually returns an `EnrichmentResult` object (wrapping `.results` as a DataFrame, plus accessor methods). This was the single largest source of confusion during review — I initially thought `experiment_id` wasn't a column on `result.results`. Propose: fix the doc block to say "returns `EnrichmentResult`; call `.to_envelope(...)` for the MCP dict shape" and cross-reference the Per-result fields table as the authoritative column list.

**A3. Cluster naming convention (`{experiment_id}|{timepoint}|{direction}`, NaN → `"NA"`) is not prominent in the docs.** Every drill-down pattern using `result.explain(cluster, term_id)` depends on this convention. Currently it's a bullet in Common mistakes. Propose: move to the top of the Response format section or create a dedicated "Cluster naming" subsection.

**A4. `ontology_landscape` default `limit=10` silent truncation.** Easy to assume the full landscape is returned. Propose: either change the default to `None`, or surface the truncation with a `truncated=True` signal that's checked at the top of the doc (it IS returned in the envelope but not flagged prominently).

### 4.6 Plan-authoring hygiene checklist (new)

Beyond the existing "spec-completeness" self-review in plans, a separate "execution-correctness" checklist would catch most of what this review found. Propose adding to `writing-plans` skill (or `step-protocol.md`):

- [ ] Every inline Python script template has an empty-guard on every DataFrame read / filter / groupby result.
- [ ] No string-parsing heuristics on structured IDs (experiment_id, cluster, term_id).
- [ ] No hardcoded condition-specific strings in equality checks without an adjacent non-empty validation.
- [ ] Every `dict(zip(df.col, df.col))` or `df.col.tolist()` targeted at a per-experiment lookup has a `drop_duplicates` first, or is explicitly safe because the frame is already at experiment granularity.
- [ ] Every `except Exception: continue` has a justification comment or is upgraded to fail-fast.
- [ ] Bash blocks: every task's first two lines are `cd <repo root>` + `export ANALYSIS_DIR="$(cat .analysis_dir)"`. No inherited-shell assumptions.
- [ ] `git add` uses explicit file lists or bash arrays, never raw shell globs.
- [ ] Long matplotlib / visualization code in plans uses skeleton + `# TWEAK:` pattern, not full polished code.
- [ ] Cross-script shared utilities live in a named module importable via `sys.path` (no `importlib.util.spec_from_file_location` hacks for 04-style filenames).
- [ ] For every risk with ≥3 escalation paths in the spec, the plan has an empirical probe task.

### 4.7 Review-process meta-observation

This review was ~80% pattern-matching for hygiene issues that aren't domain-specific and ~20% analysis-specific findings (pickle empirical validation, API gap, spec-plan drift). The hygiene 80% is mechanizable: a pre-execution review sub-skill (or the code-reviewer subagent primed with §4.6 checklist) could catch most of it automatically. Worth experimenting with for future multi-step research plans before they go to execution.

---

## 5. Consolidation and next steps

1. During B2 execution, `analyses/<date>-<slug>/gaps_and_friction.md` captures live friction.
2. At B2 completion, merge the "Skill/methodology friction" and "Proposed changes to the skill" sections from `gaps_and_friction.md` into §2 and §3 here.
3. Items in §4 (plan-review conclusions) are already generalizable — they can feed skill v3 directly without waiting for B2 execution to finish. But validate during B2 that the hygiene checklist (§4.6) actually prevents the patterns it targets, and add any new patterns surfaced during execution.
4. When this doc has enough mass (end of B2 + maybe one more analysis), it becomes the input to a skill-v3 brainstorm with the same research-methodology skill restructuring pattern as `2026-04-06-research-skill-restructure-design.md`.

---

## 6. Cross-references

- [B2 design spec](2026-04-18-pathway-enrichment-b2-design.md) — the analysis this brainstorm produced.
- [B2 execution plan](../plans/2026-04-18-pathway-enrichment-b2-plan.md) — the plan the §4 review hardened.
- [B1 gaps_and_friction](../../analyses/2026-04-09-1713-pathway_enrichment_b1/gaps_and_friction.md) — prior analysis friction log; some items overlap and strengthen the case for v3 changes.
- [Prior skill-improvement specs](.) — `2026-04-06-research-skill-restructure-design.md`, `2026-04-07-research-notebook-discipline-design.md`, `2026-04-08-research-methodology-v2-improvements-design.md` — establish the pattern for skill-change proposals.
