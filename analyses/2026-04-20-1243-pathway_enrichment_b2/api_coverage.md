# API coverage — pathway enrichment B2

Running log of MCP tools and Python API primitives exercised during the analysis, with what worked and where friction surfaced. Seeded during Step 2; completed at Task 12 (Step 5 do) with the final close-out.

For the comprehensive methodology-level friction analysis, see `gaps_and_friction.md`. This file is the lightweight coverage-and-observations log for API evaluation (goal #3 in spec §2).

---

## MCP tools

| Tool | Step(s) used | Purpose | Worked | Friction / gaps |
|---|---|---|---|---|
| `list_organisms()` | 1a MCP orientation | Surface available organisms for class-selection. | ✅ | None. |
| `list_experiments(summary=True)` | 1a MCP orientation | Top-level breakdown of experiments by organism / treatment / table_scope. | ✅ | None. |
| `list_experiments(search_text="nitrogen", verbose=True)` | 1a MCP orientation | Cross-organism N-related experiment discovery via Lucene relevance. | ✅ | Useful; surfaced cross-organism CTX candidates. |
| `list_experiments(organism="MED4", verbose=True)` | 1a | Full MED4 landscape for T/R/PC/NC classification. | ✅ | `gene_count` field is cumulative across TPs, which is misleading for detection-power interpretation — see `gaps_and_friction.md` "gene_count misreported as cumulative instead of per-timepoint" (skill-lesson + API-docs refinement candidate). Use `tp_gene_count` from `experiments_to_dataframe` instead. |
| `list_publications(organism="MED4", verbose=True)` | 1a | Attribution of experiment_id → first_author / title for the classification table. | ✅ | Discovery needed — the `first_author` join isn't automatic on `list_experiments` output. Skill-friction: `list_experiments` provides `publication_doi` but not authors; author attribution requires a separate `list_publications` call. Documented as anti-hallucination rule: never attribute authors from memory. |
| `list_filter_values("treatment_type")` | 1a (didn't end up using) | Sanity-check N-tag vocabulary. | n/a | Skipped — Lucene search on "nitrogen" covered the space adequately for B2. |
| `ontology_landscape(organism=..., experiment_ids=[...], verbose=True)` | 1b | Per-(ontology, level) coverage / relevance ranking to pick Pick 1 (cyanorak_role L1) and Pick 2 (kegg L2). | ✅ | `example_terms` column (nested list) dropped by `to_dataframe`. `by_ontology` envelope summary only available in raw response — captured in notebook entry rather than CSV. |
| `search_ontology("nitrogen", ontology=X)` | 1b | Cross-ontology nitrogen-term density check across 7 ontologies. | ✅ | BRITE excluded — requires `tree` filter per BRITE subtree; worth flagging as UX improvement candidate. |
| `pathway_enrichment(organism, experiment_ids, ontology, level, background, direction, significant_only)` | 2 do | Main enrichment run — Fisher ORA per cluster with per-(organism, ontology, bg) calls. | ✅ (after upstream fix) | **KG-data bug: `de_enrichment_inputs` background collapsed to foreground** (fixed upstream 2026-04-20; see `gaps_and_friction.md` KG data bugs section). Initial run produced 0 significant results across all MED4 × table_scope clusters; root-caused in `multiomics_explorer.de_enrichment_inputs`; fixed in sibling repo; re-ran successfully. **Tool-schema hallucination** — initial claim that `term_ids` was not supported was wrong; schema inspection showed `term_ids` is fully supported for DAG ontologies. Anti-hallucination rule: always `ToolSearch select:<tool_name>` before asserting limitations. |
| `EnrichmentResult.explain(cluster, term_id)` | 2 explore | Gene-level drill-down on specific cluster × pathway combinations (`glnA` / `atpA-I` / ribosome gene presence). | ✅ | Output is markdown-formatted and returns a namespace object — `._repr_markdown_()` for rendering. Useful for chat-capture. |
| `EnrichmentResult.overlap_genes(cluster, term_id)` | 2 explore | List of `GeneRef` objects for foreground ∩ pathway — structured access to locus_tag, gene_name, log2fc, padj, rank per gene. | ✅ | `GeneRef` is unhashable — can't put directly into a `set()`. Extract `.locus_tag` for set-based overlap comparisons. |
| `EnrichmentResult.term2gene` | 2 explore (redundancy cross-check) | MED4 gene membership of any pathway term — used for cross-pathway Jaccard. | ✅ | Is a **long DataFrame** (locus_tag × term_id rows), not the dict-like structure initially expected. Column names: `locus_tag, gene_name, product, gene_category, term_id, term_name, level`. |
| `EnrichmentResult.inputs.gene_sets` / `inputs.background` | — | Per-cluster gene_set + background dicts. Not yet used. | — | Available but not exercised in B2 Step 2. |
| `EnrichmentResult.to_compare_cluster_frame()` | — | Pivoted DataFrame. Not used in Step 2. | — | Listed as available in `.explain` preview; may be useful in Step 4/5. |
| `EnrichmentResult.generate_summary()` | — | Pre-built summary. Not used. | — | Available; not exercised. |
| `EnrichmentResult.to_envelope()` | — | Envelope metadata. Not used. | — | Available; not exercised. |

---

## Python API / utility observations

| Primitive | Observation |
|---|---|
| `experiments_to_dataframe(...)` | Produces one row per (experiment × timepoint) — correct for per-TP gene_count inspection. Top-level `gene_count` on the raw response object is cumulative (sum across TPs), which is wrong for per-TP interpretation. Flagged in gaps_and_friction.md Skill/methodology friction. |
| `omics_type` column in `enrichment_all.csv` | Populated for non-Weissberg experiments but **NaN for Weissberg T experiments** — source-of-truth is `experiments_classified.csv` (or upstream `list_experiments`). Discovered during heatmap v2 construction when T panel dividers displayed "nan" instead of Prot/RNA labels. Worked around by merging from classified; worth surfacing as an upstream API consistency gap. |
| Pickle round-trip of `EnrichmentResult` dict | Verified stage-1 (single object) + stage-2 (multi-key dict). `.explain()` works on objects retrieved from the loaded dict. Precautionary escalation paths in spec §8 Risk 7 (split pickle files / dill / state-only pickle / inline re-run) were not needed at this scale (8-key dict, ~0.5 MB per entry). |
| Fisher ORA `signed_score = sign × −log10(padj)` | Behaves as log-scale evidence strength. Meaningful below \|s\|≈10; statistical saturation above — differences reflect integer-gene-count detection precision, not biology. Informed the display-cap choice in the QC heatmap. |

---

## MCP doc / example improvements noted during B2 (candidate API-docs revisions)

1. **`pathway_enrichment` `term_ids` discoverability for DAG ontologies.** Example at `docs://examples/pathway_enrichment.py` shows `level`-only usage. For DAG ontologies (go_*), a `level`-only call silently drops biologically-meaningful terms at other depths. Propose adding an example showing hand-curated `term_ids` panel paired with `search_ontology` for DAG cases. See gaps_and_friction.md.

2. **`ontology_landscape` flag DAG-vs-flat ontologies explicitly.** Current landscape returns per-(ontology, level) stats without distinguishing flat (cyanorak_role, tigr_role, kegg) from DAG (go_*) — researcher has to know this from outside the tool. Proposing an `ontology_kind` field: `flat` / `hierarchical-tree` / `dag`.

3. **`enrichment_all.csv` schema consistency for `omics_type`.** Populated inconsistently (null for Weissberg T). Either fix upstream to populate consistently, or document that downstream code should source from `experiments_classified.csv`.

4. **`list_experiments` top-level `gene_count` naming.** Rename to `total_row_count` or add sibling `genes_per_timepoint` / `n_timepoints` so the cumulative-sum-across-TPs semantics are self-documenting. Current name invites misinterpretation.

Defer the "propose doc-level patches" action to Task 12 Step 5 when the full API usage is known.

---

## Completion (Task 12 Step 5)

_To be filled at Task 12 — final close-out after Step 4/5 exercises additional primitives (`.generate_summary`, possibly `.to_compare_cluster_frame`, score-computation helpers)._
