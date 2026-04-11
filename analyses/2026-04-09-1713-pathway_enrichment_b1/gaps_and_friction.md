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

- **genome_coverage metric** — Adding genome_coverage (genes reachable at this level / total genome genes) to the ontology characterization was the single most important methodological contribution of this analysis. Without it, CyanoRak level 2 appeared to be the best choice based on term-size statistics alone. The metric revealed that level 2 covers only 18% of the genome and that level 1 is the correct selection. This metric should be a required output of any ontology landscape characterization.

- **Signed enrichment score** — Collapsing the up/down enrichment into a single signed -log10(padj) score enabled a single heatmap to show direction, magnitude, and experiment × timepoint structure simultaneously. The visualization immediately revealed the RNA/protein discordance and control separation without requiring side-by-side comparison of separate figures.

- **Marker gene tracing** — Tracing known N-limitation markers (glnA, cynA, rbcL, atpD) through each pipeline step (annotation extraction → pathway definition → enrichment results) provided a running sanity check. When glnA and cynA appeared in E.4, rbcL in J.2, and atpD in J.1, those were early signals that the pipeline was biologically coherent before the enrichment tests ran.

- **enrich_utils package** — Structuring the analysis utilities as a local importable package (`enrich_utils/`) rather than inline scripts made each subsequent script cleaner and the logic reusable. The `run_enrichment_all_timepoints`, `apply_bh_correction`, and `signed_enrichment_score` functions encapsulated the core logic and were easy to test in isolation.

### What didn't work

- **Initial CyanoRak level 2 selection** — The first pass of `02_survey_landscape.py` ranked CyanoRak level 2 by term-size statistics (median 12 genes, attractive range) without computing genome_coverage. This led to `03_define_pathways.py` being called with `--level 2`, which produced 0 marker gene matches (level 2 only covers niche subcategories). The bug was caught before enrichment was run, but it cost one iteration of scripts 2 and 3. Root cause: the survey script did not compute genome_coverage in its initial version.

- **NaN timepoint bug in enrichment runner** — The initial `run_enrichment_all_timepoints` used `groupby("timepoint")` which silently drops groups where `timepoint` is NaN. This affected Tolonen cyanate, Tolonen urea, and Steglich (single-timepoint experiments stored with NaN timepoint). The fix was to replace `groupby` with explicit NaN-aware grouping. The bug was caught immediately on first-pass QC (those experiments showed 0 tests run), but it required a full rerun of script 04.

### Proposed changes

**To the skill:**

- Add genome_coverage to the ontology selection checklist in research-methodology. Any enrichment analysis using a hierarchical ontology must compute and compare genome_coverage per level before selecting a level. This is not optional — term-size statistics alone are insufficient and misleading.
- Add a note that single-timepoint experiments may require NaN-aware grouping when DE tables use NaN for timepoint rather than a sentinel value.

**To the MCP/KG:**

- `ontology_landscape` tool (already filed as MCP friction item 1) — would have eliminated scripts 01 and 02 entirely for the landscape characterization phase.
- `genes_at_ontology_level` tool (MCP friction item 2) — would eliminate the 69-call loop in script 03 validation.
- Batch support in `gene_ontology_terms` (MCP friction item 3) — would eliminate the chunking workaround in script 01 for large gene sets.
