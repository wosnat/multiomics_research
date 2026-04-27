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

---

## F4 — Annotation gaps for stress axes in MED4

**Date:** 2026-04-27 (encountered in step 3)

**What happened.** While defining stress-axis controls in step 3, several canonical ontology terms have **0 MED4 genes** in this KG, even though the underlying biology is well-established:

| Axis | Term | Ontology | MED4 genes | Comment |
|---|---|---|---|---|
| N-stress | `go:0006995` cellular response to nitrogen starvation | GO BP | **0** | The hand-curated `cyanorak.role:D.1.3` (16 genes) and `E.4` (28 genes) cover this well, so the GO BP gap is workaround-able. |
| Cell death | `go:0008219` cell death | GO BP | **0** | No MED4 genes annotated to "cell death." |
| Cell death | `go:0012501` programmed cell death | GO BP | **0** | Same. |
| Cell death | `cyanorak.role:D.3` Cell growth and death | cyanorak | **0** | Even the hand-curated cyanorak ontology has no MED4 genes in "Cell growth and death." Cell death is structurally absent from MED4 annotations. |
| Photo-stress | photosystem-II-related GO BP terms | GO BP | psbA (PMM0223) is **not** in `go:0009769`, `go:0010206`, `go:0019684` despite being the canonical D1 protein. | psbA is properly annotated in `cyanorak.role:J.8` (Photosystem II). GO BP coverage of psbA is the gap. |
| Oxidative | catalase / `katG` | (n/a) | **0** | Not a KG gap — MED4 *biologically lacks* catalase (Black Queen Hypothesis). The gap is the absence of the gene from the genome, not from the KG. Logged here so future axes don't expect it to be present. |
| HLI proteins | `cyanorak.role:J.8` (PSII) | cyanorak | not in J.8 | The 23 high-light-inducible (HLI) proteins are not in the cyanorak PSII role; they live in `gene_category="Stress response and adaptation"` only. Cross-reference is via gene_category, not ontology. |

**Workaround in this analysis.** Step 3 control panel uses cyanorak (hand-curated, 73% MED4 coverage) where it has gene sets, falls back to GO BP where cyanorak is empty, and uses canonical literature markers (psbA, lrtA) for genes that aren't in either ontology. The cell-death axis is reframed as **"late-stationary / starvation response"** because formal cell-death annotations are unavailable; positive controls are spoT (stringent response), lrtA (ribosome hibernation), isiB (flavodoxin / starvation marker).

**Downstream impact.**
- *Methodology:* the cell-death axis cannot use a single canonical gene-set definition. It must rely on a small hand-picked set of starvation-response markers, and the axis label should be qualified ("late-stationary / starvation response") rather than promising programmed-cell-death-detection capability.
- *Step 4 (methods):* gene-set definitions per axis should mix evidence sources transparently — cyanorak primary, GO BP supplementary, literature markers explicitly tagged. The mixture must be visible in any methods write-up.
- *KG / annotation suggestion:* a future cyanorak update could populate `D.3` (Cell growth and death) with starvation/stationary-phase markers — even minimal coverage (spoT, lrtA, isiB and similar) would be enough for axis-level analyses to ground in a single ontology. Out of scope for this analysis.

---

## F5 — Step-3 heatmap reading error: "all UP" claim contradicted by data

**Date:** 2026-04-27 (caught in step 5)

**What happened.** Step 3's notebook described the step-3 control-validation heatmap with: *"All 5 N-stress positives are clearly UP across both omics in both conditions. Magnitudes biggest in axenic late TPs."* That description is wrong for **coculture-RNA**. Re-reading the underlying `3_analysis_framing/data/control_de_long.csv` directly: in the coculture-RNA experiment, glnB at day 89 has log2FC = -2.58 (significant_down); ntcA at day 89 = -2.08 (significant_down); glnA at day 60 = -2.17 (significant_down). The 5 N-stress positives are systematically **DOWN** in coculture-RNA across all 4 TPs (day 18, 31, 60, 89). The step-3 narrative anchored on the strongest cells (axenic-Prot death-phase) and missed the muted-but-systematic cocult-RNA negative pattern.

This was caught in step 5 only because the score-based view computes a quantitative `axis_mean_signed_lfc` per cell, which immediately surfaced the negative values for cocult-RNA n_stress (-1.18 at day 31, -1.82 at day 60, etc.).

**Workaround / correction.** Step 5's notebook now reports the cocult-RNA n_stress DOWN-regulation as a **central finding** (transcriptional relief in coculture: cells reduce N-scavenging transcription while the proteome stays engaged, consistent with Alteromonas-supplied N relieving the need to make more N-scavenging proteins). The step-3 claim "all UP" remains in step 3's notebook as an artifact of how the analysis evolved — fixing it retroactively would erase the methodology learning. Instead, this entry plus the step 5 notebook serve as the corrective record.

**Methodology lesson for future analyses.** Heatmap reads are easy to anchor on the strongest cells and miss the muted-but-systematic ones. **Step-3 control validation should include a per-(omics × condition) summary table of axis-mean log2FC alongside the per-gene heatmap.** A two-line table per axis showing per-condition-per-omics direction would have surfaced the cocult-RNA negative pattern at step 3 close instead of step 5. Recommend adopting this for future analyses.

**Downstream impact.**
- *Step 5:* turned what could have been a rote "apply scoring across all cells" into a major-finding-discovery step. Net positive.
- *Step 3 retrofit:* not retrofitted. The narrative-overwrite anti-pattern (overwriting step-3 retroactively) would erase the methodology lesson. Step 3's narrative stands; this entry + step 5's notebook hold the corrected reading.
- *General:* the validation-summary numerical table proposal will be adopted for future step-3 work.
