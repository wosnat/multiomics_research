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

Fisher's exact (one-sided, over-representation):
Expected: 450 × 14 / 1849 = 3.41
Observed: 8
Fold enrichment: 8 / 3.41 = 2.35×
p-value: scipy.stats.fisher_exact([[8, 442], [6, 1393]], alternative='greater')[1]

test_type: vs_genome
pathway_coverage: 14/15 = 0.93
```

Report: "Nitrogen metabolism" enriched 2.35-fold among upregulated genes at day 14 (Fisher's exact p = X, BH q = Y, 8/14 pathway genes significant_up vs 3.4 expected, pathway coverage 93%).

### Worked example (filtered_subset background)

Read N-dep at 12h: 840 genes in the universe, 200 significant_up. "Nitrogen metabolism" has 15 KG genes, but only 9 are in Read's filtered_subset.

```
                  In pathway    Not in pathway
Significant up        6              194         = 200
Not sig up            3              637         = 640
                      9              831         = 840

Expected: 200 × 9 / 840 = 2.14
Observed: 6
Fold enrichment: 6 / 2.14 = 2.80×

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

test_type: vs_all_responsive
pathway_coverage: 7/15 = 0.47 (flagged: < 50%)
```

Interpretation: the pathway is enriched *at this timepoint* relative to all genes that respond at *any* timepoint. Not the same as enrichment vs the genome.

## Ontology landscape

### Available ontologies in the KG

From the KG schema, the following ontology types have gene annotations and hierarchy relationships:

| Ontology | Total terms | Hierarchy | Levels | Gene→term | Term→genes |
|----------|----------:|-----------|--------|-----------|-----------|
| `go_bp` | 3,052 | `is_a`, `part_of`, `regulates` | 15 deep | `gene_ontology_terms` (leaf only) | `genes_by_ontology` (with expansion) |
| `go_mf` | ? | `is_a`, `part_of` | ? | same | same |
| `go_cc` | ? | `is_a`, `part_of` | ? | same | same |
| `kegg` | 4,742 | `is_a` | 4 explicit (`category`, `subcategory`, `pathway`, `ko`) | `Gene_has_kegg_ko` | `genes_by_ontology` |
| `cyanorak_role` | 173 | `is_a` | 3 (letter, X.Y, X.Y.Z) | `Gene_has_cyanorak_role` | `genes_by_ontology` |
| `tigr_role` | 114 | ? | ? | `Gene_has_tigr_role` | `genes_by_ontology` |
| `cog_category` | 26 | flat | 1 | `Gene_in_cog_category` | `genes_by_ontology` |
| `ec` | ? | `is_a` | 4 (X, X.Y, X.Y.Z, X.Y.Z.W) | `Gene_catalyzes_ec_number` | `genes_by_ontology` |
| `pfam` | ? | `Pfam_in_pfam_clan` | 2 (domain, clan) | `Gene_has_pfam` | `genes_by_ontology` |

Steps 1-2 will fill in the "?" entries and assess MED4-specific coverage.

### MCP tool capabilities and gaps

**Available:**
- `gene_ontology_terms(locus_tags, ontology)` — leaf annotations per gene (gene→terms)
- `genes_by_ontology(term_ids, ontology, organism)` — genes for a term with hierarchy expansion (term→genes), but one call per term
- `search_ontology(search_text, ontology)` — text search for term IDs
- `run_cypher(query)` — raw Cypher, escape hatch for hierarchy queries

**Gap: no hierarchy traversal tool.** To survey the hierarchy structure (depth, parent-child relationships, terms per level), we must use `run_cypher`. This is a known limitation. The analysis will document what queries were needed, as input for a future dedicated tool.

**Gap: no bulk gene→term extraction.** `gene_ontology_terms` paginates (default limit 5). For extracting all annotations for ~1,800 genes, use the Python API import:
```python
from multiomics_explorer import gene_ontology_terms
result = gene_ontology_terms(locus_tags=[...], ontology="cyanorak_role", limit=None)
```

**Gap: no "genes at hierarchy level X" tool.** To get all genes at CyanoRak level 2 (X.Y), we need either:
- `genes_by_ontology` called once per level-2 term (up to ~70 calls for CyanoRak)
- `run_cypher` to batch the query
Scripts will use the Python API + `run_cypher` for bulk extraction.

## Pipeline: 7 steps

Each step follows the do→show→explore→decide cycle. No step proceeds without the researcher's decision. Notebook entry committed before the next step begins.

### Step 1: Extract MED4 annotations from KG (interactive discovery)

**Goal:** Pull all gene-to-term mappings for every ontology, scoped to the gene universes from the v2 experiments.

**Method:**
- Collect gene universes: union of locus tags per experiment from v2 `de_*.csv` files. Record per-experiment universe sizes.
- For each ontology type: use Python API `gene_ontology_terms(locus_tags=universe, ontology=X, limit=None)` to get leaf annotations
- For hierarchical ontologies: use `run_cypher` to extract the full hierarchy (parent-child edges, level assignments)
- Record which genes have no annotations in any ontology (the "dark genome" for enrichment purposes)

**enrich_utils function:**
- `extract_annotations(locus_tags, ontology, api)` → DataFrame of gene × term leaf annotations
- `extract_hierarchy(ontology, api)` → DataFrame of parent-child edges with level assignments

**Outputs:**
- `data/annotations_*.csv` — one per ontology, gene × leaf term mappings
- `data/hierarchy_*.csv` — one per hierarchical ontology, parent-child edges
- `data/gene_universes.csv` — per-experiment gene universe with counts
- `logs/01_extract_annotations.log` — per-ontology: gene count, term count, genes with no terms

**Explore:** Which ontologies have the best coverage for each gene universe? How does coverage differ between RNA-seq (1,849 genes) and proteomics (1,424 genes)? How many genes are unannotated in all ontologies? Are unannotated genes enriched among DE genes (annotation bias)?

**Decision:** Any ontologies to drop as obviously unsuitable (empty or near-empty for MED4)?

### Step 2: Characterize annotation landscape

**Goal:** For each ontology, assess suitability for enrichment testing across all v2 experiment gene universes.

**Method — `enrich_utils.survey`:**
- Per ontology, per gene universe:
  - Coverage: fraction of genes with ≥1 annotation
  - Term count (at each hierarchy level for hierarchical ontologies)
  - Term-size distribution: min, Q1, median, Q3, max genes per term
  - Terms with <3 genes (underpowered), terms with >200 genes (too broad)
- Hierarchy analysis (for GO, KEGG, CyanoRak, TIGR, EC):
  - Depth/levels, terms per level
  - Gene-set size distribution at each level
  - Redundancy: fraction of genes annotated to both a term and its parent
- Identify "sweet spot" level per ontology: terms large enough for power (≥5 genes in the universe), specific enough for interpretation
- Annotation bias check: are unannotated genes disproportionately DE in any experiment?

**enrich_utils functions:**
- `survey_ontology(annotations_df, hierarchy_df, gene_universe)` → OntologyProfile with coverage, term stats, per-level stats
- `rank_ontologies(profiles)` → ranked list with scores and rationale

**Outputs:**
- `data/ontology_profiles.csv` — per ontology × level: coverage, term count, term-size stats, sweet-spot flag
- `results/ontology_ranking.csv` — ranked ontologies with scores
- `results/annotation_landscape.png` — visual comparison across ontologies
- `logs/02_characterize_landscape.log`

**Explore:** Walk through profiles. Which ontologies give good power at which levels? CyanoRak level 2 vs KEGG pathway vs GO at depth 5 — what do the term sizes look like for our gene universes? Is annotation biased?

**Decision:** Rank ontologies. Select top choice for enrichment. Justify.

### Step 3: Define pathway gene sets

**Goal:** Extract pathway definitions from the selected ontology at the chosen hierarchy level, for each experiment's gene universe.

**Method:**
- From step 1 annotations + hierarchy, filter to selected ontology + level
- For hierarchical ontologies at a rolled-up level: a gene belongs to a parent term if it is annotated to any descendant of that term
  - For CyanoRak/KEGG: hierarchy is in `hierarchy_*.csv`, roll-up is a join
  - For GO: use `genes_by_ontology(term_ids=[parent_term], organism="MED4")` per parent term, or `run_cypher` for batch extraction
- Per pathway definition: pathway ID, pathway name, full gene set (all MED4 genes), per-experiment gene set (intersection with that experiment's universe)
- Filter: minimum gene set size in the experiment universe (threshold decided during explore — likely ≥5)

**MCP tool limitation:** Getting genes at a non-leaf level requires either:
- Calling `genes_by_ontology` once per term at the chosen level (works for CyanoRak ~70 terms, expensive for GO ~400 terms)
- Using `run_cypher` for batch queries
The script will use whichever is appropriate for the selected ontology. Document the query pattern for the MCP tool requirements output.

**enrich_utils function:**
- `build_pathway_definitions(annotations_df, hierarchy_df, level, min_genes)` → DataFrame of pathway_id, pathway_name, locus_tags (set), gene_count
- `scope_pathways_to_universe(pathway_defs, gene_universe)` → same DataFrame with per-universe gene sets and coverage fractions

**Outputs:**
- `data/pathway_definitions.csv` — pathway ID, name, gene count (full genome), gene list
- `data/pathway_coverage_per_experiment.csv` — pathway × experiment: genes in universe, coverage fraction
- `logs/03_define_pathways.log` — filtering funnel, size distribution, N-related pathways check

**Explore:** How many pathways survive the size filter? What's the median gene set size? Do known N-related pathways appear (nitrogen metabolism, photosynthesis, transport)? Any pathways with poor coverage in proteomics but good coverage in RNA-seq?

**Decision:** Pathway definitions look reasonable? Minimum size threshold correct? Proceed to enrichment?

### Step 4: Toy-data verification of enrichment functions

**Goal:** Verify Fisher's exact enrichment, FDR correction, and pathway coverage functions with hand-calculated synthetic data before real data.

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

Each test: synthetic DE DataFrame + synthetic pathway definitions → expected Fisher's exact result (computed by hand) → assert match.

**enrich_utils tests:**
- `enrich_utils/tests/test_enrichment.py` — all cases above
- `enrich_utils/tests/test_survey.py` — synthetic annotation data, verify coverage and term-size stats

**Outputs:**
- `enrich_utils/tests/test_enrichment.py` — test script with synthetic data and assertions
- `enrich_utils/tests/test_survey.py` — test script for survey functions
- `logs/04_verify_enrichment.log` — expected vs actual for each test case

**Explore:** Walk through each test case. Verify edge cases. Does the enrichment function have the properties we want?

**Decision:** Tests pass? Functions correct? Proceed to real data?

### Step 5: Run enrichment tests

**Goal:** Fisher's exact test per pathway × experiment × timepoint × direction, across all v2 experiments.

**Method — `enrich_utils.enrichment`:**
- For each experiment × timepoint:
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
- `results/enrichment_all.csv` — full results: pathway, experiment, timepoint, direction, p_value, padj, odds_ratio, fold_enrichment, a, b, c, d, test_type, pathway_coverage, n_pathway_genes_in_universe, n_pathway_genes_total
- `results/enrichment_significant.csv` — filtered to padj < 0.05
- `logs/05_run_enrichment.log` — total tests, tests per experiment, significant counts, N-pathway traces

**Explore:** Which pathways are enriched in RNA-seq axenic? In coculture? Do they differ? Does proteomics show the same or different pathways? Trace nitrogen-related and photosynthesis pathways specifically. How many tests total, how many survive FDR? Compare enrichment in positive controls (reference studies) vs negative controls — do the right pathways light up?

**Decision:** Results make biological sense? Surprises? Proceed to visualization?

### Step 6: Visualize

**Goal:** Heatmaps and summary figures.

**Outputs:**
- `results/heatmap_enrichment_up.png` — pathways × experiments/timepoints, colored by -log10(padj), for upregulated genes
- `results/heatmap_enrichment_down.png` — same for downregulated
- `results/pathway_comparison_axenic_vs_coculture.png` — side-by-side for RNA-seq and proteomics
- `results/pathway_comparison_rnaseq_vs_proteomics.png` — same conditions, different platforms
- `results/control_enrichment.png` — reference and negative control enrichment patterns

### Step 7: Interpret and document

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

Two layers:
- **Survey** (`survey.py`) — annotation landscape characterization, ontology profiling, ranking. Operates on annotation DataFrames. May use KG Python API for extraction.
- **Enrichment** (`enrichment.py`) — Fisher's exact, FDR correction, pathway coverage. Pure DataFrame-in, DataFrame-out. Never calls the KG.
- **Hierarchy** (`hierarchy.py`) — hierarchy extraction and level roll-up. Uses KG Python API / `run_cypher` for extraction, pure DataFrame operations for roll-up.
- **I/O** (`io.py`) — load/save helpers.

### Modules

```
enrich_utils/
├── __init__.py
├── survey.py          # Annotation landscape: coverage, term stats, ranking
├── enrichment.py      # Fisher's exact, FDR, pathway coverage (pure computation)
├── hierarchy.py       # Hierarchy extraction and level roll-up
├── io.py              # Load/save helpers
└── tests/
    ├── test_enrichment.py   # Toy data tests for Fisher's exact + FDR
    └── test_survey.py       # Toy data tests for coverage and profiling
```

### Flow

```
Gene universes + KG
    → extract_annotations()           # gene × leaf term per ontology
    → extract_hierarchy()             # parent-child edges per ontology

Annotations + hierarchy + gene universe
    → survey_ontology()               # per-ontology profile
    → rank_ontologies()               # ranked list with scores

Selected ontology + hierarchy
    → build_pathway_definitions()     # pathway gene sets at chosen level
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

Note: pipeline step 4 (toy-data verification) runs the test suite in `enrich_utils/tests/`, not a numbered pipeline script. The test scripts are the artifact. Pipeline scripts are numbered 01-05 mapping to steps 1-3 and 5-6.

## Marker pathways

Starting set for tracing through every step (confirmed after step 2 selects the ontology):

| Expected pathway | Expected behavior | Why |
|-----------------|-------------------|-----|
| Nitrogen metabolism | Enriched up in axenic RNA-seq | Core N-limitation response |
| Photosynthesis | Enriched down in axenic RNA-seq | N-stress suppresses photosynthesis |
| Protein synthesis / ribosomal | Enriched down in N-stress | Growth shutdown |
| Transport | Enriched up in N-stress | Nutrient scavenging |

Specific term IDs depend on the selected ontology — assigned in step 3.

## Directory structure

```
analyses/YYYY-MM-DD-HHMM-pathway_enrichment_b1/
├── exploration/
│   └── YYYY-MM-DD-notebook.md
├── data/
│   ├── DATA_MANIFEST.md
│   ├── gene_universes.csv
│   ├── annotations_*.csv          # One per ontology
│   ├── hierarchy_*.csv            # One per hierarchical ontology
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
│   ├── survey.py
│   ├── enrichment.py
│   ├── hierarchy.py
│   ├── io.py
│   └── tests/
│       ├── test_enrichment.py
│       └── test_survey.py
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
