# Add Cluster Analysis Info to List Tools

## Goal

Surface clustering analysis summary data in `list_publications`, `list_organisms`, and `list_experiments` so researchers see cluster coverage without needing to call `list_clustering_analyses` separately. Mirrors how expression data (experiment counts, treatment types, omics types) is already surfaced in these tools.

## What changes for the researcher

Before: a researcher calling `list_publications()` sees experiment counts and treatment types but has no idea whether clustering analyses exist for a study. They must separately call `list_clustering_analyses` to discover this.

After: each publication, organism, and experiment result includes `clustering_analysis_count` and `cluster_types` (always), plus `cluster_count` and `clustered_gene_count` (verbose). Top-level summaries include a `by_cluster_type` breakdown.

## New fields

### Per-result fields

| Field | Type | Mode | Description |
|---|---|---|---|
| `clustering_analysis_count` | int | always | Number of linked ClusteringAnalysis nodes |
| `cluster_types` | list[str] | always | Distinct cluster types (e.g. `["response_pattern", "diel_expression_pattern"]`) |
| `cluster_count` | int | verbose | Total GeneCluster nodes across linked analyses |
| `clustered_gene_count` | int | verbose | Total genes assigned to clusters across linked analyses |

Added to all three tools: `list_publications`, `list_organisms`, `list_experiments`.

### Top-level summary breakdown

| Field | Type | Description |
|---|---|---|
| `by_cluster_type` | list[{cluster_type, count}] | Count of matching results (publications/organisms/experiments) that have at least one clustering analysis of that type. Same semantics as `by_treatment_type`. |

Added to all three tools. Same pattern as existing `by_treatment_type`, `by_omics_type`.

## Approach: materialize on nodes

Cluster summary fields are materialized as properties on `OrganismTaxon`, `Publication`, and `Experiment` nodes during KG build, matching the existing pattern for `experiment_count`, `treatment_types`, `omics_types`, etc.

### Graph edges used

| Node | Edge to ClusteringAnalysis | Direction |
|---|---|---|
| `OrganismTaxon` | `ClusteringanalysisBelongsToOrganism` | CA → Organism |
| `Publication` | `PublicationHasClusteringAnalysis` | Pub → CA |
| `Experiment` | `ExperimentHasClusteringAnalysis` | Exp → CA |

All are direct one-hop edges (no chaining through intermediate nodes).

### Computation

For each node, aggregate over its linked ClusteringAnalysis nodes:

- `clustering_analysis_count` = count of linked CAs
- `cluster_types` = distinct `CA.cluster_type` values
- `cluster_count` = sum of `CA.cluster_count` across linked CAs
- `clustered_gene_count` = sum of `CA.total_gene_count` across linked CAs

Nodes with no linked CAs get `clustering_analysis_count=0`, `cluster_types=[]`, `cluster_count=0`, `clustered_gene_count=0`.

## Implementation layers

### Layer 0: Validate Cypher against live KG

Run the proposed enrichment Cypher (read-only, no SET) against the live Neo4j to confirm edge directions, property names, and expected counts before writing any code. Verify that the aggregated values match `list_clustering_analyses` summary output.

### Layer 1: KG enrichment (biocypher_kg)

Add a post-import enrichment step (or extend the existing one that computes `experiment_count`, `treatment_types`, etc.) to materialize the four properties on `OrganismTaxon`, `Publication`, and `Experiment` nodes.

Example Cypher for OrganismTaxon:
```cypher
MATCH (o:OrganismTaxon)
OPTIONAL MATCH (ca:ClusteringAnalysis)-[:ClusteringanalysisBelongsToOrganism]->(o)
WITH o,
     count(ca) AS ca_count,
     collect(DISTINCT ca.cluster_type) AS ctypes,
     sum(coalesce(ca.cluster_count, 0)) AS total_clusters,
     sum(coalesce(ca.total_gene_count, 0)) AS total_genes
SET o.clustering_analysis_count = ca_count,
    o.cluster_types = ctypes,
    o.cluster_count = total_clusters,
    o.clustered_gene_count = total_genes
```

Similar for Publication (via `PublicationHasClusteringAnalysis`) and Experiment (via `ExperimentHasClusteringAnalysis`).

### Layer 2: Cypher queries (multiomics_explorer `queries_lib.py`)

**Detail queries** for all three tools: add the four new properties to the RETURN clause. `cluster_count` and `clustered_gene_count` are always returned from Cypher; verbose gating happens in Python.

**Summary queries** for all three tools: collect `cluster_types` arrays, flatten, and compute `by_cluster_type` via `apoc.coll.frequencies()`. Same pattern as existing `by_treatment_type`.

Files:
- `queries_lib.py` lines ~586-699 (list_publications detail + summary)
- `queries_lib.py` lines ~712-752 (list_organisms)
- `queries_lib.py` lines ~827-1007 (list_experiments detail + summary)

### Layer 3: API + MCP (multiomics_explorer `functions.py`, `tools.py`)

**`functions.py`:**
- All three API functions: include `clustering_analysis_count` and `cluster_types` in results always
- Gate `cluster_count` and `clustered_gene_count` behind `verbose=True`
- Add `by_cluster_type` to summary output

**`tools.py`:**
- Update Pydantic response models for all three tools to include the new fields
- Update tool docstrings to mention clustering info availability

Files:
- `functions.py` lines ~521-548 (list_organisms), ~551-636 (list_publications), ~639-817 (list_experiments)
- `tools.py` lines ~135-172 (list_organisms), ~1024-1098 (list_publications), ~1187-1325 (list_experiments)

### Layer 4: About YAML + regenerate (multiomics_explorer)

**Input YAML updates** (`inputs/tools/`):

`list_publications.yaml`:
- Add `cluster_count` and `clustered_gene_count` to `verbose_fields`
- Update example response to include `clustering_analysis_count`, `cluster_types`
- Add chaining: `"list_publications → list_clustering_analyses(publication_doi=[...])"`

`list_organisms.yaml`:
- Add `cluster_count` and `clustered_gene_count` to `verbose_fields`
- Update example response to include `clustering_analysis_count`, `cluster_types`
- Add chaining: `"list_organisms → list_clustering_analyses(organism=...)"`

`list_experiments.yaml`:
- Add `cluster_count` and `clustered_gene_count` to `verbose_fields`
- Update example response to include `clustering_analysis_count`, `cluster_types`
- Add chaining: `"list_experiments → list_clustering_analyses(experiment_ids=[...])"`

**Regenerate:**
```bash
cd /home/osnat/github/multiomics_explorer
uv run python scripts/build_about_content.py list_publications list_organisms list_experiments
```

### Layer 5: Documentation (multiomics_research)

Update `docs/superpowers/specs/2026-03-31-gene-cluster-tools-what-changed.md` to note the new fields on list tools.

### Layer 6: KG spec doc (biocypher_kg)

Create a spec/task doc in biocypher_kg describing the new materialized properties and the enrichment Cypher. This is the handoff for the KG rebuild — implementation of Layer 1 waits for the KG update.

## Current KG data (as of 2026-04-05)

15 clustering analyses across 7 organisms:
- By organism: NATL2A (4), MED4 (4), BP-1 (2), M. ruber (2), MIT1002 (1), MIT9313 (1), MIT9301 (1)
- By cluster type: response_pattern (7), periodicity_classification (2), diel_expression_pattern (2), expression_classification (1), diel_cycling (1), expression_pattern (1), expression_level (1)

After materialization, e.g.:
- `OrganismTaxon "Prochlorococcus MED4"` would get `clustering_analysis_count=4`, `cluster_types=["response_pattern", "diel_expression_pattern", ...]`
- Publications without cluster data get `clustering_analysis_count=0`, `cluster_types=[]`

## Out of scope

- No new MCP tools — this adds fields to existing tools only
- No changes to `list_clustering_analyses`, `gene_clusters_by_gene`, or `genes_in_cluster`
- No filtering by cluster_type on the list tools (can be added later if needed)
- No changes to `gene_overview` or `gene_response_profile`
