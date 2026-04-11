# Data Manifest

All files produced by extraction scripts from the multiomics KG.
Run scripts from `multiomics_research` root with `uv run`.

---

## Gene universes

| File | Script | Description |
|------|--------|-------------|
| `gene_universes.csv` | `01_extract_annotations.py` | Per-experiment gene universe (genes in each DE table). 10 experiments from v2 `de_*.csv` files. Columns: `experiment`, `gene_id`. |

---

## Annotations

One file per ontology. Extracted from the multiomics KG via `gene_ontology_terms` and `genes_by_function`. All files have columns: `gene_id`, `term_id`, `term_name`, `level` (where applicable).

| File | Ontology | Genes annotated | Terms |
|------|----------|---:|---:|
| `annotations_cyanorak_role.csv` | CyanoRak functional role | 1,504 | 141 |
| `annotations_cog_category.csv` | COG category | 1,751 | 19 |
| `annotations_ec.csv` | EC numbers | 828 | 700 |
| `annotations_go_bp.csv` | GO Biological Process | 1,122 | 669 |
| `annotations_go_cc.csv` | GO Cellular Component | 805 | 121 |
| `annotations_go_mf.csv` | GO Molecular Function | 1,200 | 886 |
| `annotations_kegg.csv` | KEGG pathway / KO | 1,065 | 1,017 |
| `annotations_pfam.csv` | Pfam domain | 1,492 | 1,527 |
| `annotations_tigr_role.csv` | TIGR role | 1,765 | 106 |

---

## Hierarchy files

Term-to-parent relationships extracted from the KG. Used for rolling up leaf annotations to the chosen hierarchy level. Columns: `child_id`, `child_name`, `parent_id`, `parent_name`, `level`.

| File | Ontology | Note |
|------|----------|------|
| `hierarchy_cyanorak_role.csv` | CyanoRak | 3 levels (0=root, 1=subcategory, 2=leaf) |
| `hierarchy_ec.csv` | EC numbers | 4 levels |
| `hierarchy_go_bp.csv` | GO Biological Process | 12 levels |
| `hierarchy_go_cc.csv` | GO Cellular Component | 7 levels |
| `hierarchy_go_mf.csv` | GO Molecular Function | 10 levels |
| `hierarchy_kegg.csv` | KEGG | 4 levels (category→subcategory→pathway→ko) |
| `hierarchy_pfam.csv` | Pfam | 0 edges (Pfam clan relationships not present in KG for MED4) |

---

## Ontology characterization

| File | Script | Description |
|------|--------|-------------|
| `ontology_profiles.csv` | `02_survey_landscape.py` | Per-ontology × level statistics. Columns: `ontology`, `level`, `genome_coverage`, `n_terms`, `median_genes`, `max_genes`. 43 rows. |

---

## Pathway definitions

| File | Script | Description |
|------|--------|-------------|
| `pathway_definitions.csv` | `03_define_pathways.py` | 69 CyanoRak level 1 pathways with ≥5 genes. Columns: `term_id`, `term_name`, `n_genes`, `gene_ids`. |
| `pathway_coverage_per_experiment.csv` | `03_define_pathways.py` | Per-pathway coverage in each experiment's gene universe. Columns: `term_id`, `experiment`, `n_genes_in_universe`, `coverage_fraction`. |
