# Brainstorming Log: N-Limitation Signature v2

**Date:** 2026-04-07 to 2026-04-08
**Participants:** Researcher + Claude

## Context

The original analysis (v1) at `analyses/2026-04-06-1432-n_limitation_signature/` completed Approach A but:
- Methodology skill loaded too late (retrofitted)
- Results are a black box — researcher can't follow or defend the method
- Process shortcuts (skipped reviews, untested plan code)

Goals for v2:
1. Clean method suitable for packaging
2. Deep understanding of the signature scoring approach
3. Correct methodology from the start

## Decisions

### Q1: Redo or walk-through?
**Options:** (A) Walk through existing scripts/outputs, fixing as we go. (B) Rebuild from scratch, reusing concepts.
**Decision:** B — rebuild from scratch. Cleaner for packaging, avoids anchoring to potentially wrong decisions.

### Q2: Scope — which parts?
**Options:** All four parts (extraction, signature, scoring, Alteromonas) or just 1-3.
**Decision:** Steps 1-3 only (reference extraction, signature building, scoring). Alteromonas N-recycling is a separate analysis.

### Q3: Study selection approach?
**Options:** (A) Keep same Tolonen+Read pair. (B) Start from KG discovery, let selection emerge. (C) Other.
**Decision:** B — data-driven from KG. Part of understanding the method from the ground up.

### Q4: Scoring metrics?
**Options:** (A) One metric, understood thoroughly. (B) Carry forward all three. (C) Decide after seeing data.
**Decision:** A — rank score as sole primary metric. One understood metric > three black boxes.

### Q5: How interactive?
**Options:** (A) Deep — trace 3-5 genes at each step. (B) Moderate — QC + deep-dive on surprises. (C) Light.
**Decision:** A — deep. This is how we break the black box.

### Q6: Where does reusable logic live?
**Options:** (A) Inside analysis dir, productize later. (B) Shared location from the start.
**Decision:** A — analysis-first, per the skill's phase 1 rule.

### Q7: Signature application in utils?
**Decision:** Yes — `apply_signature()` is reusable methodology, not script glue. Returns the inspectable intermediate (which genes matched, concordance, values).

### Q8: Per-step scripts with logging?
**Decision:** Yes. Scripts write diagnostic logs to `logs/`. Explore phase split: script diagnostics (reproducible) + chat reasoning (notebook).

### Q9: Where do logs sit?
**Decision:** `logs/` inside the analysis directory, alongside data and results.

### Q10: Companion Jupyter notebooks?
**Decision:** One optional notebook after step 6 (not per-step). Decided at that point whether it's worth creating.

### Q11: Table scope handling?
**Decision:** Capture `table_scope` from DE response envelope, add as column on every row in CSVs. Critical for interpreting absent genes.

### Q12: Toy tests as scripts?
**Decision:** Yes — saved as `sig_utils/tests/test_*.py` to seed the test suite for productization.

### Q13: Shared extraction utility?
**Decision:** `sig_utils/extraction.py` wraps DE extraction logic. All download scripts use it.

### Q14: Permutation testing scope?
**Decision:** Part of scoring, not limited to specific experiments. `score_with_significance()` wraps rank_score + permutation into one call.

### Q15: Superpowers products in analysis dir?
**Decision:** Yes — spec, plan, and brainstorm log copied into `superpowers/` so decision history is self-contained.

### Q16: Positive/negative controls?
**Decision:** Every discovered experiment gets classified. Positive controls (N-stress) must score high, negative controls (non-N or early timepoints) must score low. Gate: controls must separate before scoring targets.

### Q17: Additional experiments as controls?
**Decision:** Yes — select additional experiments from the KG as positive and negative controls beyond the reference studies.

### Q18: Experiment scoping table columns?
**Decision:** Added treatment_type and background_factors to the scoping output.
