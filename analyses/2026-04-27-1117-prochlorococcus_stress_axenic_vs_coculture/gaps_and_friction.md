# Gaps and friction — Prochlorococcus stress, axenic vs co-culture (Weissberg 2025)

Append-only log of methodology / KG / tooling friction encountered during this analysis.
Decisions live in each step's `notebook.md`; this file captures friction (gaps, schema mismatches, anti-hallucination corrections, process slowdowns).

---

## F1 — Cross-condition statistical contrast at matched TPs not supported by the data

**Date:** 2026-04-27 (logged in step 1; cause characterized in step 2 on 2026-04-27)

**What happened.** During step-1 brainstorming Q3 (comparison structure), the researcher selected option (a) within-condition trajectory and explicitly noted: *"would like b [across-condition contrast at matched TPs] but we don't have the data for it."* The specific cause was deferred to step 2.

**Step 2 characterization (2026-04-27).** The KG selection in step 2 (see `2_kg_selection/notebook.md`) makes the structural cause explicit:

| Omics | Axenic experiment | Coculture experiment | Calendar TPs shared | Phys.-shared (same TP + same phase) |
|---|---|---|---|---|
| RNA-seq | `..._med4_rnaseq_axenic` — **not time-course**, single nutrient_limited contrast | `..._med4_rnaseq_coculture` — 4 TPs (day 18 / 31 / 60 / 89, all nutrient_limited) | none (axenic has no per-TP RNA) | none |
| Proteomics | `..._med4_proteomics_axenic` — 3 TPs: day 14 (nutrient_limited), day 31 (death), day 89 (death) | `..._med4_proteomics_coculture` — 4 TPs: day 18, 31, 60, 89, all nutrient_limited | day 31, day 89 | **none** — calendar-shared days are in different physiological states |

Two structural issues:
1. **No per-TP axenic RNA-seq.** The axenic RNA experiment is recorded as a single `is_time_course=False` contrast against exponential, with no per-TP breakdown — so a per-TP RNA contrast across conditions is not possible.
2. **Phase mismatch on calendar-shared TPs.** Day 31 and day 89 exist in both axenic and coculture proteomics, but axenic is already in *death phase* by day 31, while coculture is still *nutrient_limited* at the same calendar day (and at all sampled TPs). A direct DE contrast at those days would conflate "stress mitigation by Alteromonas" with "axenic and coculture being in fundamentally different physiological states at the same calendar day."

The phase mismatch is biologically informative — it is itself a sign that axenic *Prochlorococcus* enters death much earlier than coculture under N-deprivation — but it precludes a clean matched-TP statistical contrast for the cross-condition claim.

**Workaround.** Scope locked at step 1 to within-condition trajectories only. The cross-condition reading remains visual / descriptive: shape and timing of trajectories side by side, with the implicit acknowledgment that "calendar-aligned" axenic and coculture TPs are not physiologically comparable.

**Downstream impact.**
- *Step 3 framing:* must build the within-condition trajectory framing on top of these per-experiment DE-vs-exponential structures (each TP is already log2FC vs the same condition's exponential phase). The within-condition baseline is therefore implicitly defined by the experiments themselves.
- *Step 5 visualization:* trajectories should plot against `timepoint_hours`, not just calendar day labels, and should annotate physiological phase per point so viewers can tell which TPs are nutrient_limited vs death.
- *Statistical claims this analysis cannot make:* "coculture has less stress than axenic at TP X." It *can* make: "axenic accumulates X-axis stress earlier (e.g., entering death by day 31) while coculture stays in nutrient_limited at all sampled TPs."

---

## F2 — `list_experiments(organism=...)` filter is ambiguous (matches profiled organism OR coculture partner)

**Date:** 2026-04-27 (encountered in step 2)

**What happened.** Initial selection script called `list_experiments(publication_doi=[DOI], organism="Prochlorococcus MED4", limit=None)` expecting 5 *Prochlorococcus* MED4 experiments. The API returned 6 — the extra one was an HOT1A3 experiment (HOT1A3 profiled, *Prochlorococcus* MED4 as the coculture partner). The API's docstring states the filter is a *partial match on profiled organism AND coculture partner*, so this is documented behavior, but it is silently wrong for the common case "give me experiments where MED4 is the organism being profiled."

**Workaround.** Drop the `organism` parameter from the API call; post-hoc filter on `organism_name == "Prochlorococcus MED4"` in pandas. Logged in `2_kg_selection/scripts/01_select_experiments.py` with a comment block explaining why.

**Downstream impact.**
- *KG / API gap:* the same risk applies to any analysis that scopes by organism via `list_experiments` and is unaware of the coculture-partner match behavior. A safer API design would separate `profiled_organism=` from `coculture_partner=` (the latter already exists as a parameter; the former does not).
- *Methodology:* every analysis should verify experiment counts against expectation immediately after the filter call. The filter funnel pattern (logged here as `publication filter -> N1 / organism filter -> N2`) catches this kind of leak.
- *Suggested KG-side fix:* add a `profiled_organism=` parameter that matches only `organism_name`, and keep `organism=` as the current behavior for backward compatibility — or rename `organism=` to `any_organism=` to make the OR semantics explicit. To be raised against `multiomics_explorer`.

---

## F3 — `experiments_to_dataframe` does not expand `time_point_growth_phases` per TP row

**Date:** 2026-04-27 (encountered in step 2)

**What happened.** `experiments_to_dataframe(result)` correctly expands the `timepoints` list to one row per TP, populating per-TP columns (`timepoint`, `timepoint_order`, `timepoint_hours`, `tp_gene_count`, `tp_significant_up/down/not_significant`). But the parallel array `time_point_growth_phases` (one phase per timepoint, ordered by `timepoint_order`) is **not zipped against the TP rows.** Instead, the whole list is flattened to a `" | "`-joined string and copied identically onto every TP row of that experiment.

Concretely, for `..._med4_proteomics_axenic`:

| timepoint | timepoint_order | `time_point_growth_phases` (function output) |
|---|---|---|
| day 14 | 1 | `"nutrient_limited \| death \| death"` |
| day 31 | 2 | `"nutrient_limited \| death \| death"` |
| day 89 | 3 | `"nutrient_limited \| death \| death"` |

To recover the per-TP phase (day 14 → nutrient_limited; day 31 → death; day 89 → death), the consumer must zip `result["results"][i]["timepoints"]` with `result["results"][i]["time_point_growth_phases"]` themselves, by `timepoint_order` index.

**Workaround in this analysis.** A local helper `attach_per_tp_phase()` in `2_kg_selection/scripts/01_select_experiments.py` builds a `(experiment_id, timepoint_order) -> phase` lookup from the raw envelope and adds a `tp_growth_phase` column to the DataFrame returned by `experiments_to_dataframe`.

**Downstream impact.**
- *Methodology:* every analysis that uses `experiments_to_dataframe` to ask "what was the growth phase at this TP?" silently gets the *experiment-level joined string* instead of the per-TP value. If a downstream script does `df[df["time_point_growth_phases"] == "death"]` or similar, it filters on the join string, not on the per-TP phase — silently wrong.
- *Suggested KG-side fix (deferred to later):* in `experiments_to_dataframe`, when expanding `timepoints` per row, also emit a `tp_growth_phase` column populated from `time_point_growth_phases[timepoint_order - 1]`. The experiment-level string column can stay (or be renamed `experiment_growth_phases_joined`) for backward compatibility. To be raised against `multiomics_explorer` later — user has indicated this is a gap to handle later, not now.
