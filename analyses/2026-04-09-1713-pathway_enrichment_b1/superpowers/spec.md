# Pathway Enrichment Analysis B1: Annotation Landscape and First Enrichment

**Date:** 2026-04-09
**Status:** Active spec
**Predecessor:** `analyses/2026-04-08-1038-n_limitation_signature_v2/` — completed Approach A (reference signature scoring). This is the first Approach B analysis: pathway-level enrichment.

## Question

Which functional pathways are enriched among differentially expressed genes in each Weissberg 2025 condition/timepoint, and how does the enrichment pattern differ between axenic and coculture, RNA-seq and proteomics?

Sub-questions:
- What ontology annotation landscape exists for MED4 in the KG?
- Which ontology gives the best pathway definitions for enrichment testing (coverage, granularity, hierarchy level)?
- Which pathways are enriched per condition/timepoint?
- Do pathway-level results explain the RNA/protein discordance found in Approach A?
- Are there condition-specific pathway patterns the global signature score missed?

## Relationship to Approach A

Approach B is a **parallel analysis of the same data through a different lens**, not an extension of A.

- **A asks:** do these specific reference-validated genes respond? (gene-centric, reference-anchored)
- **B asks:** do coherent functional groups respond, regardless of whether they were in the references? (pathway-centric, ontology-anchored)

**Shared inputs:**
- Same `de_*.csv` DE extracts from v2 (no re-extraction)
- Same experiment/control classification from v2's `experiment_scoping.csv`

**Independent outputs:** enrichment results can find biology that A's signature missed (genes outside Tolonen + Read), and can decompose A's global score into pathway-level contributions.

They converge in interpretation, and eventually in Approach C (hybrid).

## Goals

1. **Research notebook and analysis artifacts** — interactive exploration notebook guiding and recording every step (do→show→explore→decide), plus full artifact set per the research-methodology skill: methods.md, decisions.md, caveats.md, gaps_and_friction.md, manifests, references.md
2. **Annotation landscape characterization** — systematic survey of all KG ontologies for MED4, scoped to the gene universes in the v2 experiments
3. **Ranked ontology recommendations** — scored by coverage, granularity, biological relevance, with rationale. Reusable for future analyses.
4. **Pathway enrichment results** — using the top-ranked ontology, across all v2 experiments
5. **Clean `enrich_utils` package** — reusable enrichment utilities with toy-tested functions, sibling to `sig_utils`, candidate for later productization
6. **MCP tool requirements** — document what hierarchy queries this analysis needed via `run_cypher`, as input for a future dedicated ontology hierarchy tool

## Out of scope

- Enrichment with second/third ontologies (future B analyses)
- Approach C hybrid integration
- Alteromonas pathway analysis
- GSEA / rank-based enrichment (different method, different analysis)
- Productization of `enrich_utils` into `multiomics_explorer`

## Data sources

Reuses all DE data from v2. No new extraction needed.

| File | Experiment | Genes | table_scope |
|------|-----------|------:|-------------|
| `de_ref_tolonen_ndep.csv` | Tolonen N-deprivation | 1,697 | all_detected_genes |
| `de_ref_read_ndep.csv` | Read N-depleted | 840-853 | filtered_subset |
| `de_ctrl_tolonen_cyanate.csv` | Tolonen cyanate | 1,697 | all_detected_genes |
| `de_ctrl_tolonen_urea.csv` | Tolonen urea | 1,697 | all_detected_genes |
| `de_ctrl_aharonovich_coculture.csv` | Aharonovich coculture | 1,714 | all_detected_genes |
| `de_ctrl_steglich_high_white_light.csv` | Steglich high light | 198 | filtered_subset |
| `de_weissberg_rnaseq_axenic.csv` | Weissberg RNA-seq axenic | 1,849 | all_detected_genes |
| `de_weissberg_rnaseq_coculture.csv` | Weissberg RNA-seq coculture | 1,849 | all_detected_genes |
| `de_weissberg_proteomics_axenic.csv` | Weissberg proteomics axenic | 1,424 | all_detected_genes |
| `de_weissberg_proteomics_coculture.csv` | Weissberg proteomics coculture | 1,424 | all_detected_genes |

Gene universe varies per experiment. Read varies slightly per timepoint (840-853). All others are stable across timepoints.

## Enrichment test design

### The 2×2 table

For each pathway × experiment × timepoint × direction:

|  | In pathway | Not in pathway |
|--|-----------|---------------|
| DE in this direction | a | b |
| Not DE in this direction | c | d |

- **Background** (a+b+c+d) = all genes returned for that experiment × timepoint
- **Foreground** (a+b) = genes significant in the tested direction (up or down), or significant in either direction (combined test)
- **Pathway membership** (a+c) = genes in the pathway that are present in the background

**One-sided Fisher's exact test** for over-representation. Report odds ratio alongside p-value.

### Direction-aware testing

Three tests per pathway × experiment × timepoint:
- **Up:** DE genes = significant_up only
- **Down:** DE genes = significant_down only
- **Combined:** DE genes = significant in either direction

The up/down split is essential for biological interpretation — a pathway with 10 up and 10 down genes may not be "enriched" in the combined test but is clearly responding in both directions.

### Background set rule by table_scope

| table_scope | Background | Test type | Interpretation |
|------------|-----------|-----------|----------------|
| `all_detected_genes` | All genes in experiment × timepoint | `vs_genome` | Pathway enriched vs full measured genome |
| `filtered_subset` | All genes in experiment × timepoint | `vs_filtered_genome` | Pathway enriched vs filtered genome (caveat: genes excluded before testing) |
| `significant_any_timepoint` | All genes in experiment × timepoint | `vs_all_responsive` | Pathway enriched at this timepoint vs all genes responsive at any timepoint |
| `significant_only` | No valid background | `descriptive_only` | Overlap reported but no statistical test |

The `test_type` is recorded as a column in the output, not a branching code path — the Fisher's exact computation is identical; only the interpretation differs.

### Multiple testing correction

Benjamini-Hochberg FDR correction applied within each experiment × timepoint × direction group. Not across experiments — each experiment is an independent test.

### Pathway coverage annotation

For each pathway × experiment: what fraction of the pathway's genes are in this experiment's gene universe? Report as `pathway_coverage` column. Flag tests where coverage < 50% as underpowered.

### Worked example (Fisher's exact, direction-aware)

Suppose the "Nitrogen metabolism" pathway has 15 genes annotated in the KG. Of those, 14 are in the Weissberg RNA-seq axenic experiment (1,849 genes total). At day 14, 450 genes are significant_up.

Of the 14 pathway genes in the background, 8 are significant_up.

```
                  In pathway    Not in pathway
Significant up        8              442         = 450
Not sig up            6             1393         = 1399
                     14             1835         = 1849

odds_ratio, pvalue = scipy.stats.fisher_exact([[8, 442], [6, 1393]], alternative='greater')
odds_ratio = (8 × 1393) / (442 × 6) = 11144 / 2652 = 4.20
Fold enrichment (observed/expected): 8 / (450 × 14 / 1849) = 8 / 3.41 = 2.35×

test_type: vs_genome
pathway_coverage: 14/15 = 0.93
```

Note: odds ratio (4.20) and fold enrichment (2.35×) measure different things. Odds ratio compares odds of being in the pathway given DE vs not DE. Fold enrichment compares observed count to expected under independence. Both are reported; fold enrichment is more intuitive for interpretation.

Report: "Nitrogen metabolism" enriched among upregulated genes at day 14 (Fisher's exact p = X, BH q = Y, odds ratio = 4.20, fold enrichment = 2.35×, 8/14 pathway genes significant_up vs 3.4 expected, pathway coverage 93%).

### Worked example (filtered_subset background)

Read N-dep at 12h: 840 genes in the universe, 200 significant_up. "Nitrogen metabolism" has 15 KG genes, but only 9 are in Read's filtered_subset.

```
                  In pathway    Not in pathway
Significant up        6              194         = 200
Not sig up            3              637         = 640
                      9              831         = 840

odds_ratio, pvalue = scipy.stats.fisher_exact([[6, 194], [3, 637]], alternative='greater')
odds_ratio = (6 × 637) / (194 × 3) = 3822 / 582 = 6.57
Fold enrichment: 6 / (200 × 9 / 840) = 6 / 2.14 = 2.80×

test_type: vs_filtered_genome
pathway_coverage: 9/15 = 0.60
```

The caveat: 6 pathway genes were excluded by Read's expression filter. They might be DE but were never tested. This inflates fold enrichment (smaller denominator) but the p-value is computed honestly against the measured universe.

### Worked example (significant_any_timepoint)

If an experiment has table_scope `significant_any_timepoint` with 300 genes total, and at timepoint T, 80 are significant_up. "Nitrogen metabolism" has 7 genes in this set.

```
                  In pathway    Not in pathway
Significant up        5               75         = 80
Not sig up            2              218         = 220
                      7              293         = 300

odds_ratio, pvalue = scipy.stats.fisher_exact([[5, 75], [2, 218]], alternative='greater')
odds_ratio = (5 × 218) / (75 × 2) = 1090 / 150 = 7.27
Fold enrichment: 5 / (80 × 7 / 300) = 5 / 1.87 = 2.67×

test_type: vs_all_responsive
pathway_coverage: 7/15 = 0.47 (flagged: < 50%)
```

Interpretation: the pathway is enriched *at this timepoint* relative to all genes that respond at *any* timepoint. Not the same as enrichment vs the genome.

## Ontology landscape

### Available ontologies in the KG

From the KG schema, the following ontology types have gene annotations and hierarchy relationships:

| Ontology | Total terms | Gene→term edge | Hierarchy edge(s) | Levels | Notes |
|----------|----------:|----------------|-------------------|--------|-------|
| `go_bp` | 3,052 | `Gene_involved_in_biological_process` | `is_a`, `part_of` | 15 deep | Also has `regulates` edges (positive/negative) |
| `go_mf` | ? | `Gene_enables_molecular_function` | `is_a`, `part_of` | ? | |
| `go_cc` | ? | `Gene_located_in_cellular_component` | `is_a`, `part_of` | ? | |
| `kegg` | 4,742 | `Gene_has_kegg_ko` | `Kegg_term_is_a_kegg_term` | 4 (`category`, `subcategory`, `pathway`, `ko`) | Genes connect at `ko` level only; `KeggTerm.level` property |
| `cyanorak_role` | 173 | `Gene_has_cyanorak_role` | `Cyanorak_role_is_a_cyanorak_role` | 3 (letter, X.Y, X.Y.Z) | Cyanobacteria-curated |
| `tigr_role` | 114 | `Gene_has_tigr_role` | — (flat) | 1 | No hierarchy in KG |
| `cog_category` | 26 | `Gene_in_cog_category` | — (flat) | 1 | Single-letter codes |
| `ec` | ? | `Gene_catalyzes_ec_number` | `Ec_number_is_a_ec_number` | 4 (X, X.Y, X.Y.Z, X.Y.Z.W) | |
| `pfam` | ? | `Gene_has_pfam` | `Pfam_in_pfam_clan` | 2 (domain, clan) | |

MCP tools: `gene_ontology_terms` returns leaf annotations (gene→term direction). `genes_by_ontology` finds genes for a term with hierarchy expansion (term→gene direction, one call per term). For bulk extraction: use Python API imports or `run_cypher` with the edges above.

Steps 1-2 will fill in the "?" entries and assess MED4-specific coverage.

### MCP tool capabilities and gaps

**Available:**
- `gene_ontology_terms(locus_tags, ontology)` — leaf annotations per gene (gene→terms)
- `genes_by_ontology(term_ids, ontology, organism)` — genes for a term with hierarchy expansion (term→genes), but one call per term
- `search_ontology(search_text, ontology)` — text search for term IDs
- `run_cypher(query)` — raw Cypher, escape hatch for hierarchy queries

**Direct graph edges per ontology:** The KG has specific relationships for each ontology type, defined in `ONTOLOGY_CONFIG` in `multiomics_explorer/kg/queries_lib.py`. Scripts can use `run_cypher` with these edges directly instead of going through the `gene_ontology_terms` tool:

| Ontology | Gene→term edge | Hierarchy edge(s) | Notes |
|----------|---------------|-------------------|-------|
| `go_bp` | `Gene_involved_in_biological_process` | `is_a`, `part_of` | Two hierarchy edge types |
| `go_mf` | `Gene_enables_molecular_function` | `is_a`, `part_of` | Two hierarchy edge types |
| `go_cc` | `Gene_located_in_cellular_component` | `is_a`, `part_of` | Two hierarchy edge types |
| `kegg` | `Gene_has_kegg_ko` | `Kegg_term_is_a_kegg_term` | Genes connect at `ko` level only; `KeggTerm.level` property for hierarchy |
| `cyanorak_role` | `Gene_has_cyanorak_role` | `Cyanorak_role_is_a_cyanorak_role` | 3 levels via `code` (A, A.1, A.1.1) |
| `tigr_role` | `Gene_has_tigr_role` | — (flat) | No hierarchy edges |
| `cog_category` | `Gene_in_cog_category` | — (flat) | No hierarchy edges |
| `ec` | `Gene_catalyzes_ec_number` | `Ec_number_is_a_ec_number` | 4-level numeric hierarchy |
| `pfam` | `Gene_has_pfam` | `Pfam_in_pfam_clan` | 2 levels: domain → clan |

**Gap: no hierarchy traversal tool.** To survey the hierarchy structure (depth, parent-child relationships, terms per level), we must use `run_cypher` with the edges above. This is a known limitation. The analysis will document what queries were needed, as input for a future dedicated tool.

**Gap: no bulk gene→term extraction at scale.** `gene_ontology_terms` paginates (default limit 5). For extracting all annotations for ~1,800 genes, use the Python API import:
```python
from multiomics_explorer import gene_ontology_terms
result = gene_ontology_terms(locus_tags=[...], ontology="cyanorak_role", limit=None)
```
Or use `run_cypher` with the direct gene→term edges for maximum control.

**Gap: no "genes at hierarchy level X" tool.** To get all genes at CyanoRak level 2 (X.Y), we need either:
- `genes_by_ontology` called once per level-2 term (up to ~70 calls for CyanoRak)
- `run_cypher` with hierarchy expansion: `MATCH (g:Gene)-[:Gene_has_cyanorak_role]->(leaf)-[:Cyanorak_role_is_a_cyanorak_role*0..]->(parent) WHERE parent.code =~ '^[A-Z]\\.[0-9]+$'`
Scripts will use the Python API + `run_cypher` for bulk extraction.

**Cypher verification rule:** Every `run_cypher` query must be tested interactively before embedding in a script — run it with a small limit, inspect the output, verify it returns the expected shape and values. Cypher is easy to get subtly wrong (wrong traversal direction, missing deduplication, unintended Cartesian products). Treat each query as a mini do→show→verify step.

**Explorer layer rules:** When writing Cypher queries or using the Python API, follow the architecture conventions in `multiomics_explorer/.claude/skills/layer-rules`. Key points:
- `queries_lib.py` has `ONTOLOGY_CONFIG` with the exact edge names and hierarchy relationships per ontology — use these, don't guess
- Python API functions (e.g., `gene_ontology_terms()`, `genes_by_ontology()`) return `dict` with `total_matching`, `results`, etc. — use `limit=None` for full extraction
- Cypher conventions: `$param` placeholders (never f-string user input), `AS snake_case` aliases, APOC available (`apoc.coll.frequencies()` etc.)
- `collect()` silently drops NULLs from scalar values — use `collect({key: val})` maps instead
- `WITH-as-rename` drops other variables — carry UNWIND loop vars explicitly

**Complexity escape hatch:** If the hierarchy queries become too complex or error-prone for `run_cypher`, take a 2-step approach:
1. Document the requirement as an MCP tool spec for `multiomics_explorer` (what the tool should accept and return)
2. Pause this analysis step, wait for the tool to be delivered
3. Resume using the new/updated MCP tool

This is preferable to writing fragile Cypher that produces wrong results silently. The MCP tool requirements doc (goal 6) captures these needs regardless.

## Pipeline: two phases

### Phase 0: Build, test, and validate `enrich_utils` (before full analysis)

Build the utility functions, verify with toy data, then test with real marker genes/terms to discover MCP gaps early. This phase ends with either "tools work, proceed to phase 1" or "file MCP enhancement, wait for delivery."

#### Step 0a: Define and toy-test enrichment functions

**Goal:** Implement and verify the core enrichment computation (Fisher's exact, FDR, pathway coverage) with hand-calculated synthetic data.

**enrich_utils functions:**
- `run_enrichment(de_df, pathway_defs, gene_universe, table_scope)` → DataFrame with one row per pathway × direction, with p_value, padj, odds_ratio, counts, test_type, pathway_coverage
- `run_enrichment_all_timepoints(de_df, pathway_defs)` → runs per timepoint, concatenates

**Test cases:**

1. **Strong enrichment (up):** 10-gene pathway, 8 DE up out of 100 total, 30 DE up total. Hand-calculated expected odds ratio, p-value.
2. **No enrichment:** 10-gene pathway, 1 DE up out of 100 total, 30 DE up total. p-value should be non-significant.
3. **Down enrichment:** same as case 1 but testing down direction.
4. **Combined test:** pathway with 5 up + 5 down out of 10 members. Combined test vs separate up/down.
5. **Empty pathway in universe:** pathway has 10 KG genes but 0 in this experiment's universe. Coverage = 0, test skipped.
6. **Underpowered pathway:** pathway has 10 KG genes but 3 in universe. Coverage = 0.3, flagged.
7. **All pathway genes are DE:** edge case — all 10 in pathway are DE up. Maximum enrichment.
8. **FDR correction:** 20 pathways tested, 2 nominally significant. Verify BH correction produces expected q-values.
9. **test_type assignment:** verify that table_scope maps correctly to test_type in output.

**Outputs:**
- `enrich_utils/tests/test_enrichment.py` — test script with synthetic data and assertions

**Explore:** Walk through each test case. Verify edge cases.

**Decision:** Enrichment functions correct? Proceed to survey/hierarchy functions.

#### Step 0b: Define and toy-test survey and hierarchy functions (pure DataFrame)

**Goal:** Implement hierarchy roll-up (gene→non-leaf term), survey profiling, and pathway definition functions. All pure DataFrame operations — no KG, no MCP. Tested entirely with synthetic data.

**enrich_utils functions:**
- `roll_up_to_level(annotations_df, hierarchy_df, level)` → DataFrame of gene × rolled-up term at the chosen level
- `survey_ontology(annotations_df, hierarchy_df, gene_universe)` → OntologyProfile with coverage, term stats, per-level stats
- `rank_ontologies(profiles)` → ranked list with scores
- `build_pathway_definitions(annotations_df, hierarchy_df, level, min_genes)` → DataFrame of pathway_id, pathway_name, locus_tags, gene_count
- `scope_pathways_to_universe(pathway_defs, gene_universe)` → per-universe gene sets and coverage fractions

**Test cases (toy):**
- Synthetic 3-level hierarchy (5 roots, 15 mid-level, 40 leaves) with 100 genes. Genes annotated to **any level** (some to leaves, some to mid-level, some to roots — reflecting real KG where annotation depth varies). Genes may have **multiple annotations** (e.g., gene in both "Nitrogen metabolism" and "Transport"). Terms may share genes.
- Roll-up verification:
  - Gene annotated to a leaf → appears in leaf's parent and grandparent
  - Gene annotated directly to a mid-level term → appears in that term and its parent, but NOT in child terms
  - Gene with annotations to two different branches → appears in both roll-ups
  - Deduplication (same term): gene annotated to both a leaf and its parent → counted once at the parent level
  - Deduplication (DAG convergence): gene annotated to two child terms that share a parent → counted once in the parent. The hierarchy is a DAG (especially GO), not a tree — multiple paths converge
- Survey on synthetic data: verify coverage, term-size distribution, per-level stats (accounting for multi-annotation genes)
- Pathway definitions: verify min-gene filter, coverage calculation with multi-annotation genes

**Outputs:**
- `enrich_utils/tests/test_survey.py` — test script for survey and hierarchy functions

**Decision:** Roll-up logic correct? Survey stats match expectations?

#### Step 0c: Implement KG extraction and validate (marker genes/terms)

**Goal:** Implement the extraction functions that produce DataFrames from the KG, and validate the full pipeline (extraction → roll-up → survey → pathway defs) with real marker genes. This is where MCP friction will surface.

**enrich_utils functions:**
- `extract_annotations(locus_tags, ontology, api)` → DataFrame of gene × term annotations
- `extract_hierarchy(ontology, api, max_depth=4)` → DataFrame of parent-child edges with level assignments

**Method:**
- Select marker genes that cover the edge cases from 0b's toy tests:
  - Genes with multiple annotations in the same ontology (e.g., multiple GO BP terms)
  - Genes with annotations at different hierarchy levels (leaf vs mid-level)
  - Genes with no annotation in a given ontology
  - Genes annotated to sibling terms sharing a parent (DAG convergence)
  - Starting set from v2: glnA/PMM0920, cynA/PMM0370, rbcL/PMM0550, atpD/PMM1452, PMM0030, PMM0346 — extend during exploration to cover missing edge cases
- Implement `extract_annotations` using Python API `gene_ontology_terms(locus_tags, ontology, limit=None)`
- Implement `extract_hierarchy` using `run_cypher` with `ONTOLOGY_CONFIG` edges, capped at `max_depth`
- For CyanoRak and KEGG (shallow hierarchies): extract, roll up with 0b functions, verify against `genes_by_ontology` spot-checks
- For GO BP: attempt the same. This is where the DAG complexity will likely cause friction.
- Run `survey_ontology` and `build_pathway_definitions` on the marker gene set end-to-end

**Verify against MCP:** For each rolled-up term, call `genes_by_ontology(term_ids=[term], ontology=X, organism="MED4")` and compare gene sets. Discrepancies reveal bugs in extraction, roll-up, or gaps in the hierarchy data.

**Outputs:**
- `enrich_utils/extraction.py` — extraction functions calling existing MCP/Python API/`run_cypher`, or calling new MCP tools once delivered
- `enrich_utils/tests/test_extraction.py` — integration tests using the marker genes: extraction, roll-up, spot-check against `genes_by_ontology`
- `logs/00_validate_utils.log` — marker gene annotations per ontology, roll-up verification, discrepancies
- `mcp_tool_requirements.md` — concrete list of what queries were painful or impossible, proposed tool specs with examples

**Decision gate:**
- **Extraction + roll-up work for shallow ontologies (CyanoRak, KEGG)?** → Proceed to phase 1.
- **Works but GO is problematic?** → Proceed to phase 1 excluding GO. File MCP enhancement for GO hierarchy support. GO becomes a future B analysis.
- **Fundamental friction (even shallow ontologies fail)?** → **Stop. File MCP enhancement with concrete spec. Wait for delivery before phase 1.**

---

### Phase 1: Run the analysis (steps 1-6)

Requires: `enrich_utils` validated in phase 0, MCP tools sufficient (or enhanced).

Each step follows the do→show→explore→decide cycle. No step proceeds without the researcher's decision. Notebook entry committed before the next step begins.

### Step 1: Extract MED4 annotations from KG

**Goal:** Pull all gene-to-term mappings and hierarchies for every ontology, scoped to the gene universes from the v2 experiments.

**Method:**
- Collect gene universes: union of locus tags per experiment from v2 `de_*.csv` files. Record per-experiment universe sizes.
- Extract annotations for the **full MED4 genome**, not just experiment gene universes. Use `extract_annotations` with all MED4 locus tags (from `genes_by_function("*", organism="MED4", limit=None)` or equivalent). Filter to experiment-specific universes downstream. This gives the true genome-level denominator for coverage calculations (e.g., "78% of MED4 genes have CyanoRak annotations; RNA-seq covers 95% of those").
- For hierarchical ontologies: use `extract_hierarchy` (max_depth=4) to get hierarchy edges
- Record which genes have no annotations in any ontology (the "dark genome")

**Outputs:**
- `data/annotations_*.csv` — one per ontology, gene × leaf term mappings
- `data/hierarchy_*.csv` — one per hierarchical ontology, parent-child edges (depth ≤ 4)
- `data/gene_universes.csv` — per-experiment gene universe with counts
- `logs/01_extract_annotations.log` — per-ontology: gene count, term count, genes with no terms, hierarchy depth

**Explore:** Coverage per gene universe? RNA-seq vs proteomics? Unannotated genes? Trace the marker genes from step 0c through the full extraction — do their annotations match what was validated in 0c?

**Decision:** Any ontologies to drop?

### Step 2: Characterize annotation landscape

**Goal:** For each ontology, assess suitability for enrichment testing across all v2 experiment gene universes.

**Method — `enrich_utils.survey`:**
- Per ontology, per gene universe:
  - Coverage: fraction of genes with ≥1 annotation
  - Term count at each hierarchy level (using `roll_up_to_level` from phase 0)
  - Term-size distribution per level: min, Q1, median, Q3, max genes per term
  - Terms with <3 genes (underpowered), terms with >200 genes (too broad)
- Hierarchy analysis:
  - For flat ontologies (COG, TIGR): single level, stats complete
  - For hierarchical ontologies: terms per level, gene-set size distribution per level, redundancy
  - For GO (if included): top 4 levels only
- Identify "sweet spot" level per ontology
- Annotation bias check: are unannotated genes disproportionately DE?

**Outputs:**
- `data/ontology_profiles.csv` — per ontology × level: coverage, term count, term-size stats, sweet-spot flag
- `results/ontology_ranking.csv` — ranked ontologies with scores
- `results/annotation_landscape.png` — visual comparison across ontologies
- `logs/02_characterize_landscape.log`

**Explore:** Which ontologies give good power at which levels? CyanoRak level 2 vs KEGG pathway vs GO depth 4? What's the per-term coverage distribution across experiments — are some terms well-covered in RNA-seq but missing in proteomics?

**Decision:** Rank ontologies. Select top choice. Justify.

### Step 3: Define pathway gene sets

**Goal:** Build pathway definitions from the selected ontology at the chosen hierarchy level.

**Method:**
- Use `roll_up_to_level` + `build_pathway_definitions` on the selected ontology
- Use `scope_pathways_to_universe` per experiment
- Filter by minimum gene set size (decided during explore)

**Outputs:**
- `data/pathway_definitions.csv` — pathway ID, name, gene count, gene list
- `data/pathway_coverage_per_experiment.csv` — pathway × experiment: genes in universe, coverage fraction
- `logs/03_define_pathways.log` — filtering funnel, size distribution, N-related pathways check

**Explore:** How many pathways? Median gene set size? N-related pathways present? Coverage in proteomics vs RNA-seq?

**Decision:** Pathway definitions reasonable? Proceed to enrichment?

### Step 4: Run enrichment tests

**Goal:** Fisher's exact test per pathway × experiment × timepoint × direction, across all v2 experiments.

**Method — `enrich_utils.enrichment`:**
- For each experiment × timepoint (single-timepoint experiments like Tolonen cyanate/urea, Steglich, and Weissberg RNA-seq axenic are treated as one timepoint — no special-casing needed):
  - Load DE data, determine background (all genes at this experiment × timepoint)
  - Determine `test_type` from `table_scope`
  - For each pathway:
    - Scope pathway gene set to this background
    - Compute pathway coverage
    - If coverage = 0: skip (no genes testable)
    - Run Fisher's exact (one-sided) for up, down, combined
    - Record: odds_ratio, p_value, a/b/c/d counts, pathway_coverage
  - BH correction within this experiment × timepoint × direction group
  - Flag underpowered tests (pathway_coverage < 0.50)

**enrich_utils functions:**
- `run_enrichment(de_df, pathway_defs, gene_universe, table_scope)` → DataFrame with one row per pathway × direction, with p_value, padj, odds_ratio, counts, test_type, pathway_coverage
- `run_enrichment_all_timepoints(de_df, pathway_defs)` → runs per timepoint, concatenates

**Outputs:**
- `results/enrichment_all.csv` — full results: pathway, experiment, timepoint, direction, p_value, padj, odds_ratio, fold_enrichment, observed (a), expected, a, b, c, d, test_type, pathway_coverage, n_pathway_genes_in_universe, n_pathway_genes_total
- `results/enrichment_significant.csv` — filtered to padj < 0.05
- `logs/04_run_enrichment.log` — total tests, tests per experiment, significant counts, N-pathway traces

**Explore:** Which pathways are enriched in RNA-seq axenic? In coculture? Do they differ? Does proteomics show the same or different pathways? Trace nitrogen-related and photosynthesis pathways specifically. How many tests total, how many survive FDR? Compare enrichment in positive controls (reference studies) vs negative controls — do the right pathways light up?

**Decision:** Results make biological sense? Surprises? Proceed to visualization?

### Step 5: Visualize

**Goal:** Heatmaps and summary figures.

**Outputs:**
- `results/heatmap_enrichment_up.png` — pathways × experiments/timepoints, colored by -log10(padj), for upregulated genes
- `results/heatmap_enrichment_down.png` — same for downregulated
- `results/pathway_comparison_axenic_vs_coculture.png` — side-by-side for RNA-seq and proteomics
- `results/pathway_comparison_rnaseq_vs_proteomics.png` — same conditions, different platforms
- `results/control_enrichment.png` — reference and negative control enrichment patterns

### Step 6: Interpret and document

**Goal:** Update analysis documents, assess enrich_utils for packaging, produce MCP tool requirements.

**Outputs:**
- Updated: `methods.md`, `decisions.md`, `caveats.md`, `gaps_and_friction.md`
- Updated: `data/DATA_MANIFEST.md`, `results/RESULTS_MANIFEST.md`
- `mcp_tool_requirements.md` — what `run_cypher` queries were needed, proposed tool API for ontology hierarchy traversal
- Ranked ontology list as backlog for next B analysis
- Package-readiness assessment for `enrich_utils`

## Utility package: enrich_utils

### Boundary

Sibling to `sig_utils` — independent package in the analysis directory. Shares the DE data format but not logic. No imports from `sig_utils`.

Four modules:
- **Extraction** (`extraction.py`) — KG data extraction: annotations and hierarchies. Calls Python API / `run_cypher`. Only module that touches the KG.
- **Hierarchy** (`hierarchy.py`) — hierarchy roll-up (`roll_up_to_level`). Pure DataFrame operations — takes annotations + hierarchy edges, produces gene × rolled-up term.
- **Survey** (`survey.py`) — annotation landscape characterization, ontology profiling, ranking. Pure DataFrame operations.
- **Enrichment** (`enrichment.py`) — Fisher's exact, FDR correction, pathway coverage. Pure DataFrame-in, DataFrame-out.

### Modules

```
enrich_utils/
├── __init__.py
├── extraction.py      # KG data extraction: annotations and hierarchies (only KG-touching module)
├── hierarchy.py       # Roll-up: gene × leaf term → gene × level term (pure DataFrame)
├── survey.py          # Annotation landscape: coverage, term stats, ranking (pure DataFrame)
├── enrichment.py      # Fisher's exact, FDR, pathway coverage (pure DataFrame)
├── io.py              # Load/save helpers
└── tests/
    ├── test_enrichment.py   # Toy data tests for Fisher's exact + FDR
    ├── test_survey.py       # Toy data tests for coverage, roll-up, and profiling
    └── test_extraction.py   # Integration tests against KG with marker genes
```

### Flow

```
Gene universes + KG
    → extract_annotations()           # gene × leaf term per ontology (step 1)
    → extract_hierarchy()             # parent-child edges, depth ≤ 4 (step 1)

Annotations + hierarchy + gene universe
    → survey_ontology()               # per-ontology profile (step 2)
    → rank_ontologies()               # ranked list with scores

Selected ontology + hierarchy
    → build_pathway_definitions()     # pathway gene sets at chosen level (step 3)
    → scope_pathways_to_universe()    # per-experiment coverage

DE data + pathway definitions + gene universe
    → run_enrichment()                # Fisher's exact + FDR per pathway × timepoint × direction
    → enrichment results DataFrame
```

### Toy-tested

Both survey and enrichment functions are verified with hand-calculated synthetic data before touching real KG data. Test scripts in `enrich_utils/tests/` serve as the seed for a real test suite during productization.

## Scripts

```
scripts/
├── 01_extract_annotations.py    # Pull annotations + hierarchies from KG
├── 02_survey_landscape.py       # Characterize ontology suitability
├── 03_define_pathways.py        # Build pathway gene sets at chosen level
├── 04_run_enrichment.py         # Fisher's exact across all experiments
└── 05_plot_results.py           # Heatmaps and comparison figures
```

Each script:
- Uses `enrich_utils` for all reusable logic
- Writes outputs to `data/` or `results/`
- Writes diagnostic log to `logs/`
- Has a `--explore` flag for marker pathway traces and QC diagnostics

Phase 0 (build/test/validate) uses the test suite in `enrich_utils/tests/` and interactive MCP exploration — no numbered pipeline scripts. Phase 1 scripts are numbered 01-05 mapping to steps 1-5. Step 6 (interpret) produces documentation, not script output.

## Marker genes and pathways

Traced through every step from phase 0c onward. Marker genes are carried from v2; marker pathways are assigned specific term IDs after step 2 selects the ontology.

### Marker genes

From v2, selected to cover edge cases (see step 0c):

| Gene | Locus tag | Direction (v2) | Why trace it |
|------|-----------|:-:|--------------|
| glnA | PMM0920 | up | Canonical N-limitation marker, likely multi-annotated (N-metabolism + glutamate family) |
| cynA | PMM0370 | up | N-scavenging transporter, should appear in transport pathways |
| rbcL | PMM0550 | down | Photosynthesis/CO2 fixation, tests down-enrichment |
| atpD | PMM1452 | down | Energy metabolism, tests cross-pathway annotation |
| PMM0030 | PMM0030 | up | Unnamed, rank 1 — may have no annotation in some ontologies |
| PMM0346 | PMM0346 | down | Unnamed, proteomics edge case — may have sparse annotation |

Additional genes added during step 0c to cover: multiple terms in same ontology, annotations at non-leaf level, no annotation, DAG convergence (sibling terms sharing parent).

### Marker pathways

| Expected pathway | Axenic RNA-seq | Coculture RNA-seq | Negative controls | Why |
|-----------------|:-:|:-:|:-:|-----|
| Nitrogen metabolism | Enriched up | Low/absent | Low (except cyanate) | Core N-limitation response |
| Photosynthesis | Enriched down | Low/absent | Low | N-stress suppresses photosynthesis |
| Protein synthesis / ribosomal | Enriched down | Low/absent | Low | Growth shutdown under N-stress |
| Transport | Enriched up | Low/absent | Low | Nutrient scavenging |

Specific term IDs depend on the selected ontology — assigned in step 3. Expected behaviors based on v2 signature results (axenic score 0.58, coculture near zero).

**Proteomics expectations:** May differ from RNA-seq — the RNA/protein discordance from v2 should manifest as different pathway enrichment patterns between platforms. This is a key question for B1.

## Directory structure

```
analyses/YYYY-MM-DD-HHMM-pathway_enrichment_b1/
├── exploration/
│   └── YYYY-MM-DD-notebook.md
├── data/
│   ├── DATA_MANIFEST.md
│   ├── gene_universes.csv
│   ├── annotations_*.csv          # One per ontology (leaf terms)
│   ├── hierarchy_*.csv            # One per hierarchical ontology (depth ≤ 4)
│   ├── ontology_profiles.csv
│   ├── pathway_definitions.csv
│   └── pathway_coverage_per_experiment.csv
├── scripts/
│   ├── 01_extract_annotations.py
│   ├── 02_survey_landscape.py
│   ├── 03_define_pathways.py
│   ├── 04_run_enrichment.py
│   └── 05_plot_results.py
├── enrich_utils/
│   ├── __init__.py
│   ├── extraction.py
│   ├── hierarchy.py
│   ├── survey.py
│   ├── enrichment.py
│   ├── io.py
│   └── tests/
│       ├── test_enrichment.py
│       ├── test_survey.py
│       └── test_extraction.py
├── logs/
│   ├── 01_extract_annotations.log
│   └── ...
├── results/
│   ├── RESULTS_MANIFEST.md
│   ├── ontology_ranking.csv
│   ├── annotation_landscape.png
│   ├── enrichment_all.csv
│   ├── enrichment_significant.csv
│   ├── heatmap_enrichment_up.png
│   ├── heatmap_enrichment_down.png
│   ├── pathway_comparison_*.png
│   └── control_enrichment.png
├── superpowers/
│   ├── spec.md
│   ├── plan.md
│   └── brainstorm-log.md
├── README.md
├── methods.md
├── decisions.md
├── caveats.md
├── gaps_and_friction.md
├── mcp_tool_requirements.md
└── references.md
```

## Known limitations

- **Ontology coverage bias:** Different ontologies annotate different subsets of genes. Unannotated genes are invisible to enrichment. The annotation bias check in step 2 quantifies this.
- **Hierarchy level choice is subjective:** The "sweet spot" balances power and specificity. Different levels may give different results. The ranked ontology output lets future analyses try different levels.
- **GO hierarchy complexity:** GO's 15-level DAG with `is_a`, `part_of`, and `regulates` edges creates massive redundancy. May require algorithmic pruning (elim algorithm) or limiting to a single edge type. Decided during step 2.
- **MCP tool pagination:** Bulk annotation extraction must use the Python API, not MCP. MCP is used for interactive exploration only.
- **Single ontology per analysis:** This analysis tests one ontology. Results may differ with another ontology — that's the point of the ranked backlog.
- **Read 2017 background caveat:** filtered_subset excludes low-expression genes, inflating fold enrichment for pathways with many lowly-expressed members.
- **Steglich unreliable:** Only 198 genes — most pathways will have zero or near-zero genes in the universe. Include for completeness but do not draw conclusions.
- **No correction across experiments:** FDR is within experiment × timepoint × direction. Cross-experiment comparison is descriptive, not corrected.

## Notebook discipline

One notebook: `exploration/YYYY-MM-DD-notebook.md`. Append-only, chronological. Each step gets an entry with: command, inputs, outputs, QC, exploration, decision. Notebook entry committed before the next step begins.

Steps 1-2 are interactive discovery steps — MCP queries + discussion produce frozen output files + notebook entries. Steps 3-5 are script-based with the full do→show→explore→decide cycle.

## Predecessor analysis

The v2 analysis at `analyses/2026-04-08-1038-n_limitation_signature_v2/` remains untouched. This analysis reads its `de_*.csv` files and `experiment_scoping.csv` as inputs. The notebook may reference v2 findings for comparison (e.g., "A's global score was 0.58 for axenic — B shows this is driven by N-metabolism and transport pathways").
