# Gaps and friction points

## KG data bugs

## KG gaps

1. **Pfam hierarchy missing** — `Pfam_in_pfam_clan` relationship returned 0 edges for MED4 genes. Pfam domains exist (1,492 genes annotated) but no clan grouping in the KG. Low priority — Pfam is not a candidate ontology for enrichment.

2. **KEGG genes disconnected from pathway hierarchy** — 300 of 1,065 KEGG-annotated MED4 genes (28%) have KO annotations but don't connect to any pathway via `Kegg_term_is_a_kegg_term`. Genome coverage flat at 41% across all hierarchy levels. May be KOs without pathway assignments, or missing edges in the KG.

## MCP friction

1. **No ontology landscape tool** — Characterizing an ontology's suitability for enrichment required a 4-step pipeline: extract annotations → extract hierarchy → roll up per level → compute stats. This should be one MCP call. Proposed tool:

   `ontology_landscape(ontology, organism)` → per level:
   - `genome_coverage`: genes reachable at this level / total organism genes
   - `n_terms`: terms with ≥1 gene at this level
   - `median_genes`, `max_genes`: term-size distribution
   - `example_terms`: top N terms by gene count with names

   This was the critical query pattern for selecting the right ontology and hierarchy level. Without it, the initial analysis incorrectly selected CyanoRak level 2 (18% genome coverage) over level 1 (80%).

2. **No bulk "genes at hierarchy level" tool** — Building pathway definitions at a chosen level required either 69 separate `genes_by_ontology` calls (one per term) or a custom roll-up pipeline. Proposed tool:

   `genes_at_ontology_level(ontology, level, organism)` → all genes rolled up to the specified hierarchy level, grouped by term.

3. **Neo4j memory limits on bulk extraction** — `gene_ontology_terms` for all 1,976 MED4 genes hit the 1.4 GiB transaction memory limit for GO MF. Fixed by batching in 500-gene chunks in `extract_annotations`. The MCP tool or Python API should handle batching internally.

4. **`gene_ontology_terms` returns leaf terms only** — No way to get annotations at a specific hierarchy level via MCP. Must extract leaves then roll up in Python. The proposed `genes_at_ontology_level` tool would solve this.

## Skill/methodology friction

## Process retrospective

### What worked

### What didn't work

### Proposed changes

**To the skill:**

**To the MCP/KG:**
