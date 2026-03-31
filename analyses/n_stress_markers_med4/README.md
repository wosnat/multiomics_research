# How do we know when Prochlorococcus is nitrogen stressed?

## Question

What are the reliable transcriptional markers of nitrogen stress in Prochlorococcus MED4, and how specific are they to N-stress vs other stresses?

## Key findings

- **5 genes are N-specific** across all tested stresses (7 treatment types): ureA (PMM0965), ureB (PMM0964), glnB (PMM1463), glsF (PMM1512), PMM1462 (hypothetical adjacent to glnB). These never respond to carbon, iron, light, phosphorus, viral, or salt stress.
- **ureA is the single best N-stress marker**: UP in 5/8 N experiments (13/25 timepoints), never responds to any other stress, peaks at 6h post-deprivation, stays elevated for 48h+.
- **PMM1462 is the fastest responder**: log2FC +4.7 at 3h, +5.6 at 6h — stronger than any named gene. Unknown function (hypothetical, adjacent to glnB).
- **6 more genes are N+coculture specific** (ntcA, cynB, cynD, pipX, PMM0687, urtD) — coculture response is biologically coherent (N-related).
- **A critical KG gap was identified**: `gene_response_profile` reports N-metabolism genes as `not_known` in other stresses, but they were actually measured and not significant (evidence of non-response). The tool conflates "not measured" with "measured but not significant" for experiments using `significant_only` table scope.
- **Coculture proteins persist for months** even when mRNA declines — post-transcriptional maintenance of N-scavenging machinery in the survival state.

## Files

| File | Description |
|------|-------------|
| [methods.md](methods.md) | Data sources, gene identification, approach |
| [gaps_and_friction.md](gaps_and_friction.md) | KG/MCP/skill issues found during analysis |
| [data/n_marker_de_all.csv](data/n_marker_de_all.csv) | Full DE data for 21 genes × all MED4 experiments (586 rows) |
| [data/n_marker_specificity.csv](data/n_marker_specificity.csv) | Per-gene × per-treatment specificity summary |
| [data/n_marker_timecourse.csv](data/n_marker_timecourse.csv) | Temporal profiles for Tier 1-2 genes (250 rows) |
| [scripts/extract_n_marker_profiles.py](scripts/extract_n_marker_profiles.py) | Extract DE data and build specificity table |
| [scripts/extract_n_marker_timecourse.py](scripts/extract_n_marker_timecourse.py) | Extract temporal profiles for best markers |

## Exploration log

1. [2026-03-31: Orientation — N-marker specificity](exploration/2026-03-31-orientation-n-marker-specificity.md) — KG scope, gene identification, initial specificity table, critical coverage gap identified
2. [2026-03-31: Corrected specificity](exploration/2026-03-31-corrected-specificity.md) — Table_scope reinterpretation, corrected tier classification, temporal dynamics, coculture protein persistence

## Data sources

8 nitrogen experiments + 22 other-stress experiments for MED4, from 11 publications across microarray, RNA-seq, and proteomics. See [methods.md](methods.md) for details.

## How to reproduce

```bash
cd analyses/n_stress_markers_med4
python scripts/extract_n_marker_profiles.py    # extracts DE data, builds specificity table
python scripts/extract_n_marker_timecourse.py  # extracts temporal profiles for best markers
```

Requires: multiomics_explorer Python package, Neo4j running with multiomics KG loaded.

## Date

2026-03-31
