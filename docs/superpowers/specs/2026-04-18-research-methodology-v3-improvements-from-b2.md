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

## 4. Consolidation and next steps

1. During B2 execution, `analyses/<date>-<slug>/gaps_and_friction.md` captures live friction.
2. At B2 completion, merge the "Skill/methodology friction" and "Proposed changes to the skill" sections from `gaps_and_friction.md` into §2 and §3 here.
3. When this doc has enough mass (end of B2 + maybe one more analysis), it becomes the input to a skill-v3 brainstorm with the same research-methodology skill restructuring pattern as `2026-04-06-research-skill-restructure-design.md`.

---

## 5. Cross-references

- [B2 design spec](2026-04-18-pathway-enrichment-b2-design.md) — the analysis this brainstorm produced.
- [B1 gaps_and_friction](../../analyses/2026-04-09-1713-pathway_enrichment_b1/gaps_and_friction.md) — prior analysis friction log; some items overlap and strengthen the case for v3 changes.
- [Prior skill-improvement specs](.) — `2026-04-06-research-skill-restructure-design.md`, `2026-04-07-research-notebook-discipline-design.md`, `2026-04-08-research-methodology-v2-improvements-design.md` — establish the pattern for skill-change proposals.
