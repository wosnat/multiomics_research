# Brainstorming Log: Pathway Enrichment B1

**Date:** 2026-04-09
**Participants:** Researcher + Claude

## Context

Approach A (v2) complete — 189-gene N-limitation signature, rank-score metric, clean control separation. Key findings: axenic RNA-seq scores 0.58, coculture near zero; genuine RNA/protein discordance. Now starting Approach B: pathway-level enrichment as a complementary lens.

## Decisions

### Q1: What biological questions should Approach B answer?
**Options:** (a) Which pathways drive the signature? (b) Pathway differences between axenic/coculture? (c) RNA/protein discordance at pathway level? (d) Something else?
**Decision:** All of the above, but spanning several analyses — too much for one cycle.

### Q2: What should the first B analysis (B1) tackle?
**Options:** (a) B1: pathway enrichment, B2: pathway-level discordance, B3: pathway trajectories. (b) Different decomposition.
**Decision:** B1 = annotation landscape + enrichment with one ontology. B2/B3 follow with additional ontologies and cross-platform comparison. Start with B1.

### Q3: How is B connected to A?
**Discussion:** B was designed as an independent, complementary lens — not built on A's signature. Shared DE data and experiment/control framework, but different analytical approach. A is gene-centric + reference-anchored; B is pathway-centric + ontology-anchored. They converge in interpretation and eventually in Approach C.
**Decision:** Parallel analysis, not extension.

### Q4: Enrichment background set — what is it?
**Discussion:** Fisher's exact needs a 2×2 table. Background = universe of genes that could have been detected. Data varies by experiment and table_scope.
**Decision:** Background = all genes returned for that experiment × timepoint. General rule, no special cases per experiment.

### Q5: Background per timepoint, experiment, or paper?
**Discussion:** Checked actual data — most experiments have stable gene counts across timepoints (Tolonen 1,697, Weissberg RNA-seq 1,849, Weissberg proteomics 1,424). Read varies slightly (840-853). Steglich/cyanate/urea are single-timepoint.
**Decision:** Per experiment × timepoint. For stable experiments, equivalent to per-experiment. For Read, respects per-timepoint filtering.

### Q6: How to handle different table_scopes?
**Discussion:** Four table_scope types in the KG data:
- `all_detected_genes`: clean background
- `filtered_subset`: valid background with caveat
- `significant_any_timepoint`: valid per-timepoint test (tested vs all responsive)
- `significant_only`: no background, descriptive only
**Decision:** All testable except `significant_only`. Flag test_type in output column: `vs_genome`, `vs_filtered_genome`, `vs_all_responsive`, `descriptive_only`.

### Q7: Where do pathway definitions come from?
**Discussion:** KG has 9 ontology types with varying coverage and hierarchy depth: GO BP (3,052 terms, 15 levels deep), KEGG (4,742 terms, 4 levels), CyanoRak (173 terms, 3 levels), COG (26, flat), TIGR (114), EC, Pfam.
**Decision:** Answering "which ontology" should be a large part of the analysis. First two steps are annotation landscape survey and ranking.

### Q8: Hierarchy level selection?
**Discussion:** GO has 15 levels — testing every level creates massive redundancy. CyanoRak has 3 clean levels. KEGG has 4 explicit levels. The right level balances power (enough genes per term) and specificity (interpretable biology).
**Decision:** Hierarchy level selection is a discovery step, not a preset. Determined during step 2 interactive exploration.

### Q9: MCP tool limitations for hierarchy?
**Discussion:** `gene_ontology_terms` returns leaf terms only. `genes_by_ontology` does hierarchy expansion but term-by-term. No tool for hierarchy traversal or bulk "genes at level X."
**Decision:** Use `run_cypher` for hierarchy queries. Document query patterns as input for future MCP ontology hierarchy tool.

### Q10: Relationship of enrich_utils to sig_utils?
**Options:** (a) Sibling package, independent. (b) Merge into shared `analysis_utils`.
**Decision:** (a) Sibling. Share DE data format but not logic. Productization decision comes later.

### Q11: Toy-data verification for enrich_utils?
**Decision:** Yes — same pattern as v2. Hand-calculated synthetic data with assertions. Test scripts in `enrich_utils/tests/` seed the test suite. Verification is a notebook step (pipeline step 4) before running on real data.

### Q12: MCP tool limitation in pathway definition step?
**Discussion:** Step 3 needs to get genes at a chosen hierarchy level. No bulk tool for this — need `genes_by_ontology` per term or `run_cypher`.
**Decision:** Script handles it. Document the query pattern for MCP tool requirements output.

### Q13: Scope of enrich_utils — which pipeline steps?
**Decision:** Steps 1-3 are reusable (annotation extraction, landscape survey, pathway definition). These are questions any enrichment analysis needs to answer. Step 4+ (Fisher's exact, FDR) is also reusable. All go in enrich_utils.
