# Gaps and friction points

## KG gaps

1. **KG conflates "not measured" with "measured but not significant"** — Non-N experiments use table_scope `significant_only` or `filtered_subset`, so only significant genes have KG edges. Genes absent from these tables were measured (all are microarray with ~1700 genes ≈ full genome) but not significant — evidence of non-response, not "untested." However, the KG has no way to distinguish this from "not on the platform." The `gene_response_profile` tool reports these as `groups_not_known`, which is technically correct (no edge) but biologically misleading. **Corrected interpretation:** most N-metabolism genes that appear "untested" in other stresses were actually measured and did not respond — supporting N-specificity.

2. **No `all_detected_genes` experiments for non-N stresses** — All 10 `all_detected_genes` experiments are nitrogen_stress (7) or coculture (3). Every carbon, iron, light, P, viral, and salt experiment uses `significant_only`, `filtered_subset`, or `significant_any_timepoint`. This means there is no direct way to verify whether N-metabolism genes are expressed-but-not-significant vs truly absent in those conditions without going back to the original data.

## KG data bugs

(none yet)

## MCP friction

1. **`gene_response_profile` conflates "not measured" with "not significant"** — The `groups_not_known` field reports treatment types where no expression edge exists. But for `significant_only` experiments (all non-N MED4 experiments), absence means "measured and not significant" — evidence of non-response, not unknown. The tool should ideally cross-reference experiment table_scope to distinguish `not_known` (truly no data) from `not_responded` (data exists, not significant). Workaround: manually check table_scope of experiments in each treatment type via `list_experiments`.

## Skill/methodology friction

1. **Started analysis without writing artifacts** — Ran orientation and gene identification in MCP without creating methods.md or exploration log first. The skill says to write these during the iteration, not after. Caught by user.

---

## Process retrospective

### What worked

1. **MCP → Python API escalation path.** Using `gene_response_profile` for a quick cross-treatment overview, then switching to the Python API for full extraction, worked exactly as the skill prescribes. The MCP overview surfaced the right 21 genes fast; the API extraction gave the full 586-row dataset needed for proper analysis.

2. **The iterative loop caught a real error.** The first iteration produced a specificity table with `?` (not tested) everywhere. Rather than accepting this and concluding "we can't determine specificity," the loop structure pushed toward a follow-up question ("WHY are they not tested?"). This uncovered the table_scope issue — the most important finding of the analysis. A single-pass approach would have missed it.

3. **`gene_response_profile` is the right entry point.** One call for 21 genes × all treatments gave a complete summary. Compared to the previous analysis (nitrogen_stress_med4/) which did this with jq one-offs and oversized MCP results, this was dramatically faster and more reliable.

4. **Source tagging kept claims honest.** Marking findings `[KG]` vs `[interpretation]` vs `[gap]` forced the distinction between "the data shows X" and "I think X means Y." The gap tag caught the table_scope issue before it became a false conclusion.

### What didn't work

1. **Artifacts came too late.** Ran all of orientation and half of the first research iteration before writing any files. The user had to interrupt and ask "are you following the skill guidelines?" The skill says "write during, not after" — I treated it as "write when convenient." This is the single biggest process failure. Fix: create methods.md stub and first exploration log BEFORE the first MCP call, even if they're mostly empty.

2. **Script debugging wasted time.** Three failed script runs due to: wrong organism name (`Prochlorococcus marinus MED4` vs `Prochlorococcus MED4`), wrong column name (`product` — not in API output), wrong significance check (`is_significant` vs `expression_status`). Fix: check API output schema before writing scripts. A 2-line test (`data = func(...); print(data["results"][0].keys())`) would have prevented all three errors.

3. **Didn't check what the previous analysis had already found.** Read the existing nitrogen_stress_med4/ analysis files but didn't leverage what it already established (e.g., the 5-phase model, the 108 Tolonen genes). Started from scratch instead of building on confirmed findings. The skill doesn't address how to handle prior analyses on the same topic.

4. **venv confusion.** Tried to run with `uv run --directory` instead of the project venv. Small but avoidable — should check the project's Python environment once at the start.

### Proposed changes

**To the skill:**
- Add a "pre-flight check" step: before the first MCP call, create the analysis directory and write stub files (methods.md with the question, empty gaps_and_friction.md, empty exploration/ dir). This makes the artifact-first discipline structural, not aspirational.
- Add guidance on reusing prior analyses: "Check for existing analyses on the same topic. Reference confirmed findings rather than re-deriving them. Note where you diverge."
- Add to the "Explore" section: "When using the Python API, verify the output schema with a test call before writing a full extraction script."

**To the MCP/KG:**
- `gene_response_profile` should add a `groups_tested_not_responded` field that cross-references experiment table_scope. For `significant_only` experiments, genes absent from the table ARE tested (assuming the platform covers them). This is the highest-value improvement — it would have saved an entire iteration of this analysis.
- Consider adding a `table_scope` field to the `gene_response_profile` response_summary per treatment group, so users can see at a glance whether absence is informative.

**To the Python API extraction workflow:**
- The multiomics_explorer API already returns `treatment_type`, `omics_type`, `publication_doi` in DE results — no need to merge from experiment metadata. The extraction script had unnecessary code for this. Document the API output schema somewhere accessible.
