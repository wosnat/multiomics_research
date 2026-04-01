# Gene Cluster Tools: What Changed in multiomics_explorer

## Summary

Three new MCP tools and corresponding Python API functions for querying GeneCluster nodes in the knowledge graph. These address the "No clustering/co-expression" gap from the nitrogen stress MED4 analysis.

## New MCP Tools

| Tool | Purpose | Entry point |
|---|---|---|
| `list_gene_clusters` | Browse, search (Lucene), and filter clusters by organism, cluster_type, treatment_type, omics_type, publication | Text search or filter-driven |
| `gene_clusters_by_gene` | Batch gene lookup: "what clusters are these genes in?" Single organism enforced. | Gene locus tags |
| `genes_in_cluster` | Drill into cluster membership: "what genes are in this cluster?" | Cluster IDs from tools above |

### Typical workflow

```
list_gene_clusters(search_text="photosynthesis")
    → cluster IDs

genes_in_cluster(cluster_ids=["cluster:msb4100087:med4:down_photosynthesis_psi"])
    → member genes (locus_tag, gene_name, product, gene_category)

gene_clusters_by_gene(locus_tags=["PMM0370", "PMM0920"])
    → which clusters these genes belong to
```

## New Python API Functions

```python
from multiomics_explorer import (
    list_gene_clusters,        # NEW
    gene_clusters_by_gene,     # NEW
    genes_in_cluster,          # NEW
)
```

Same name, same parameters, same return dict as the MCP tools.

### Examples

```python
# Search for clusters
result = list_gene_clusters(search_text="N transport", verbose=True)
# Returns: total_entries, total_matching, by_organism, by_cluster_type, ...
# Per result: cluster_id, name, organism_name, functional_description, ...

# Gene-centric lookup
result = gene_clusters_by_gene(locus_tags=["PMM0370", "PMM0920", "PMM0958"])
# Returns: total_clusters, genes_with_clusters, genes_without_clusters,
#          not_found, not_matched, results (one row per gene × cluster)

# Drill into cluster members
result = genes_in_cluster(
    cluster_ids=["cluster:msb4100087:med4:up_n_transport"],
    verbose=True, limit=None,
)
# Returns: by_cluster, top_categories, genes_per_cluster_max/median,
#          results (one row per gene with gene-level AND cluster-level descriptions)
```

## MCP Resource Documentation

Tool docs are available as MCP resources at `docs://tools/{tool_name}`:

```
docs://tools/list_gene_clusters
docs://tools/gene_clusters_by_gene
docs://tools/genes_in_cluster
```

Analysis utility docs at `docs://analysis/{name}`:

```
docs://analysis/response_matrix
docs://analysis/gene_set_compare
docs://analysis/to_dataframe
```

**Access via MCP:** These are resource templates (parameterized URIs). They don't appear in `ListMcpResourcesTool` but can be read directly with `ReadMcpResourceTool`:

```
ReadMcpResourceTool(server="multiomics-kg", uri="docs://tools/list_gene_clusters")
```

**Known limitation:** Resource templates don't auto-list. A future fix will register them as static resources so they appear in discovery. For now, the full list is documented here and in `CLAUDE.md`.

## Current KG Data (as of 2026-03-31)

16 GeneCluster nodes from Tolonen 2006:
- 9 MED4 clusters (5–124 member genes) — N-stress K-means
- 7 MIT9313 clusters (6–128 member genes) — N-stress K-means
- All `stress_response` type, `nitrogen_stress` treatment, `MICROARRAY` platform
- Single publication: DOI 10.1038/msb4100087

More clusters coming in KG Phase 1 (Zinser 2009 diel, ~18 clusters) and Phase 2 (Alonso-Saez, Bagby, Steglich).

### Known test genes for scripts

| Gene | Cluster | Members |
|---|---|---|
| PMM0370 (cynA) | cluster:msb4100087:med4:up_n_transport | 5 |
| PMM0920 (glnA) | cluster:msb4100087:med4:up_n_transport | 5 |
| PMM0297 | cluster:msb4100087:med4:down_translation | 124 |
| PMT0992 | cluster:msb4100087:mit9313:up_n_transport | 6 |

## Key Design Points

- **Single organism enforced** on `gene_clusters_by_gene` and `genes_in_cluster` — don't mix PMM and PMT locus tags in one call
- **`membership_score` and `p_value` are null** for Tolonen data (K-means). Will be populated for Mfuzz/soft clustering papers (Zinser 2009)
- **Verbose mode** on `genes_in_cluster` returns both gene-level fields (`function_description`, `gene_summary`) AND cluster-level fields (`functional_description`, `behavioral_description`) — clearly labeled in the docs
- **Publication filter** available on all 3 tools — useful once multiple papers contribute clusters

## Updates Needed in multiomics_research

### python-api-guide.md

Add to the import block (section 1):

```python
from multiomics_explorer import (
    # ... existing 19 functions ...
    list_gene_clusters,          # NEW — 22 total
    gene_clusters_by_gene,       # NEW
    genes_in_cluster,            # NEW
)
```

Update "19 functions" → "22 functions" in the comment.

### research skill

The `/research` skill's tool awareness should include the cluster tools. Relevant for analysis questions like:
- "Do these genes co-cluster?" → `gene_clusters_by_gene`
- "What genes cluster with X?" → `gene_clusters_by_gene` → `genes_in_cluster`
- "What clusters exist for this treatment?" → `list_gene_clusters`
