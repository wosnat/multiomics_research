# Gaps and friction — Prochlorococcus stress, axenic vs co-culture (Weissberg 2025)

Append-only log of methodology / KG / tooling friction encountered during this analysis.
Decisions live in each step's `notebook.md`; this file captures friction (gaps, schema mismatches, anti-hallucination corrections, process slowdowns).

---

## F1 — Cross-condition statistical contrast at matched TPs not supported by the data

**Date:** 2026-04-27 (logged in step 1)

**What happened.** During step-1 brainstorming Q3 (comparison structure), the researcher selected option (a) within-condition trajectory and explicitly noted: *"would like b [across-condition contrast at matched TPs] but we don't have the data for it."* The exact reason — whether it is mismatched sampling schedules between axenic and coculture, axenic perishing before coculture reaches late TPs, replicate structure, batch / normalization barriers, or another study-specific issue — was not pinned down in the dialogue and is to be characterized in step 2 when the actual KG entries (TPs, sample counts, replicate structure) are inspected.

**Workaround.** Scope of step 1 narrowed to within-condition trajectories only (option (a) at Q3). The question still answers "is each condition stressed?" but cannot directly say "is axenic more stressed than coculture in axis X at TP n?" Cross-condition reading is by visual inspection of curve shapes, not by statistical contrast. This is a real and acknowledged ceiling on what the analysis can claim.

**Downstream impact.**
- *Methodology:* the analysis cannot statistically support the central comparative claim of the abstract (that the heterotroph mitigates *specific* stress axes). It can support the within-condition claim that one condition accumulates stress and the other does not.
- *Step 2 (KG entries):* must surface and verify the actual TP / sample structure for both conditions and confirm what specifically blocks (b). If (b) turns out to be partially feasible (e.g., feasible at one or two TPs), that may motivate a redo of step 1's comparison-structure decision.
- *Possible follow-up analysis:* a separate, narrower analysis at a single matched TP (if any exists with comparable replicate structure) could later attempt the cross-condition contrast. Out of scope for this analysis.
