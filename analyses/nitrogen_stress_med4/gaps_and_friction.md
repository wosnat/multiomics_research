# Gaps and friction points

Issues discovered during this analysis that feed back into KG/MCP/skill development.

## KG data bugs

1. **Axenic RNA-seq "days 60+89" is misattributed** — The experiment `10.1101/2025.11.24.690089_growth_state_pro99lown_nutrient_starvation_med4_rnaseq_axenic` has 1849 expression edges at timepoint "days 60+89". This is actually late coculture vs exponential growth axenic — a cross-condition contrast loaded under the wrong experiment. **Action: remove from KG in biocypher_kg repo.** The only valid axenic RNA-seq timepoint is day 14.

## KG gaps

1. **No physiological data** — Fv/Fm, chlorophyll fluorescence, growth curves, cell counts from publications are not in the KG. These are critical for interpreting whether cells are alive or dead at late timepoints. (Tolonen 2006 reports Fv/Fm declining to near-zero by 48h; this context changes interpretation of the 48h transcriptome.)

2. **No diel cycle data** — Expression patterns over light/dark cycles would help distinguish N-stress from light-dependent effects.

3. **No clustering/co-expression** — Gene co-expression clusters and dynamic pattern classifications from the time courses. Would enable "which cluster does this gene belong to?" queries.

4. **No enrichment results** — GO/functional enrichment per gene set is not queryable. Must be computed in scripts.

5. **Papers themselves not in KG** — Full text, figures, supplementary data from publications not accessible. Would help with interpretation.

## Missing middle: MCP → API gap

The workflow has a gap between "browse 5 results" (MCP) and "write a full script" (Python API). Repeated need for cross-experiment gene-level aggregations that neither handles well:
- "For these N genes, which stresses do they respond to?" — needs one row per gene, not per edge
- "Show me the time course for a gene set across conditions" — needs a tidy summary table

**Workaround used:** jq/python one-offs to parse oversized MCP results saved to temp files. This is fragile and not reproducible.

**Proposed fixes:**
1. **MCP: gene-level summary mode** — Add a `summary_by="gene"` option to `differential_expression_by_gene` that returns one row per gene with: treatments responded to, directions, max/min FC, experiment count. **Change in multiomics_explorer.**
2. **This repo: reusable analysis scripts** — Python scripts using the multiomics_explorer API for common cross-experiment queries:
   - `scripts/cross_stress_profile.py` — gene list → classification table (N-specific vs nutrient vs general)
   - `scripts/timecourse_extract.py` — experiment IDs + gene list → tidy CSV for plotting
   - These save output to `data/` and are the bridge between MCP exploration and publication figures.

## MCP friction

1. **Large results exceed token limits** — Querying 17 N-metabolism genes across all experiments returns 84KB+ JSON. Had to save to file and parse with jq. The MCP tool itself works fine, but the LLM context can't hold the full result.

2. **No cross-experiment gene-level comparison tool** — Asking "does gene X respond to any non-nitrogen stress?" requires querying all experiments then filtering. A dedicated "gene specificity" or "gene response profile" tool would be useful.

3. **No batch gene query with summary** — Querying a gene set and getting a per-gene summary (which treatments, directions, max FC) in one call would reduce round-trips.

## Skill/methodology friction

1. **Capturing analysis state** — The interactive exploration generates insights but they live in chat context. Need a systematic way to checkpoint findings as we go, not just at the end.

2. **Distinguishing "explored" from "concluded"** — Some observations are preliminary (e.g., "48h might be dying") vs established (e.g., "cynA is N-specific across all experiments"). Need a way to mark confidence levels.

3. **Platform comparability caveat** — The cross-experiment comparison (microarray vs RNA-seq, different stats tests) needs prominent caveats. The skill's anti-hallucination guide covers this but it's easy to forget in the flow of analysis.
