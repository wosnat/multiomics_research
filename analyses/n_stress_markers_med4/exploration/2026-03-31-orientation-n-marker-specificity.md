# 2026-03-31: Orientation — N-stress marker specificity

## Question

What are reliable transcriptional markers of nitrogen stress in Prochlorococcus MED4, and how specific are they to N-stress vs other stresses? (exploratory)

## Approach

1. `list_organisms`, `list_publications(nitrogen, MED4)`, `list_experiments(MED4, summary)` — scope available data
2. `genes_by_function` — identify candidate N-metabolism genes (21 genes)
3. `gene_response_profile` — cross-treatment response for all 21 genes
4. Python API extraction — full DE data for 21 genes across all MED4 experiments
5. Build specificity table: gene × treatment_type showing response status

## Findings

### Scope

- [KG] MED4 has 46 experiments across 8 treatment types, 11 publications, 3 omics platforms
- [KG] 8 nitrogen_stress experiments from 3 publications (Tolonen 2006, Read 2017, Weissberg 2025)
- [KG] 3 publications with N-stress data: microarray time course (0-48h), RNA-seq time course (3-24h), long-term multi-omics (14-89 days)

### Gene identification

- [KG] 21 candidate N-metabolism genes identified by functional search: ntcA, glnA, glnB, pipX, glsF, amt1, cynABDS, urtABCDE, ureABC, and 3 hypotheticals (PMM0958, PMM1462, PMM0687)
- [KG] All 21 genes resolved to unique locus tags. No paralogs identified.

### Critical gap: gene coverage across treatment types

- [KG] DE data extracted for 21 genes: 586 rows total across all MED4 experiments
- [KG] Gene coverage per treatment type:
  - nitrogen_stress: 21/21 genes, 8 experiments, 497 rows
  - coculture: 21/21 genes, 2 experiments, 42 rows
  - carbon_stress: 9/21 genes, 3 experiments, 20 rows
  - iron_stress: 2/21 genes, 3 experiments, 8 rows
  - viral: 2/21 genes, 2 experiments, 11 rows
  - light_stress: 2/21 genes, 2 experiments, 3 rows
  - phosphorus_stress: 1/21 genes, 1 experiments, 5 rows
- [gap] 12 of 21 N-marker genes have NO data in iron, light, viral, phosphorus, or salt stress experiments. These genes cannot be classified as "N-specific" — they are "N-responsive with unknown specificity."
- [interpretation] This gap likely arises because non-N experiments used platforms or table scopes that excluded these genes. Microarray experiments would include them (all genes on chip), but RNA-seq experiments with `significant_only` or `filtered_subset` scopes would exclude non-significant genes from the results table.

### Specificity table (from gene_response_profile MCP tool)

| Gene | Name | N-stress | Coculture | Carbon | Iron | Light | P | Viral | Salt | Specificity |
|------|------|----------|-----------|--------|------|-------|---|-------|------|-------------|
| PMM0246 | ntcA | MIX | MIX | ? | ? | ? | ? | ? | ? | N-only(untested:5) |
| PMM1463 | glnB | MIX | - | ? | ? | ? | ? | ? | ? | N-only(untested:5) |
| PMM0393 | pipX | UP | DN | ? | ? | ? | ? | ? | ? | N-only(untested:5) |
| PMM0920 | glnA | MIX | MIX | MIX | ? | DN | ? | ? | ? | broad |
| PMM1512 | glsF | UP | - | ? | ? | ? | ? | ? | ? | N-only(untested:5) |
| PMM0263 | amt1 | MIX | - | DN | ? | ? | ? | UP | ? | broad |
| PMM0370 | cynA | UP(15) | UP | DN | MIX | ? | ? | ? | ? | broad |
| PMM0371 | cynB | UP(11) | MIX | ? | ? | ? | ? | ? | ? | N-only(untested:5) |
| PMM0372 | cynD | UP(11) | UP | ? | ? | ? | ? | ? | ? | N-only(untested:5) |
| PMM0373 | cynS | UP(10) | MIX | DN | ? | ? | ? | ? | ? | broad |
| PMM0970 | urtA | MIX | UP | MIX | ? | DN | DN | UP | ? | broad |
| PMM0971 | urtB | MIX | - | MIX | ? | ? | ? | ? | ? | narrow |
| PMM0972 | urtC | UP | - | MIX | ? | ? | ? | ? | ? | narrow |
| PMM0973 | urtD | UP | UP | ? | ? | ? | ? | ? | ? | N-only(untested:5) |
| PMM0974 | urtE | UP | DN | UP | ? | ? | ? | ? | ? | broad |
| PMM0963 | ureC | UP | DN | DN | ? | ? | ? | ? | ? | broad |
| PMM0964 | ureB | MIX | - | ? | ? | ? | ? | ? | ? | N-only(untested:5) |
| PMM0965 | ureA | UP(13) | - | ? | ? | ? | ? | ? | ? | N-only(untested:5) |
| PMM0958 | (DUF1830) | UP(14) | MIX | ? | MIX | ? | ? | ? | ? | broad |
| PMM1462 | (hyp) | MIX | - | ? | ? | ? | ? | ? | ? | N-only(untested:5) |
| PMM0687 | (hyp) | MIX | MIX | ? | ? | ? | ? | ? | ? | N-only(untested:5) |

Key: UP(n)=sig up in n timepoints, MIX=up in some timepoints down in others, DN=sig down, -=tested not sig, ?=not tested

### Best N-stress marker candidates (from what data we have)

- [KG] **ureA (PMM0965)**: UP in 13/25 N-stress timepoints, never responds to coculture or carbon (the only other stresses where it was tested). Strongest case for N-specificity, but untested in 5 treatment types.
- [KG] **cynB (PMM0371)**: UP in 11/25 N-stress timepoints, never down. But mixed in coculture (UP in proteomics, DOWN in one RNA-seq timepoint). Untested in 5 types.
- [KG] **cynD (PMM0372)**: UP in 11/25 N-stress timepoints, UP in coculture proteomics. Untested in 5 types.
- [KG] **PMM0958 (DUF1830)**: UP in 14/25 N-stress timepoints, best rank 1. But also responds to iron_stress. Not purely N-specific.
- [KG] **cynA (PMM0370)**: UP in 15/25 N-stress timepoints, best rank 1, median rank 4. The most consistently responsive gene. But also responds to carbon, iron, and coculture — clearly not N-specific.

### MIX signals — what do they mean?

- [KG] ntcA, glnB, glnA, urtA, PMM1462 all show "MIX" under N-stress — UP in some timepoints, DOWN in others. This reflects the time course: these genes are initially upregulated (early N-stress response) then downregulated (late/dying).
- [interpretation] A gene that goes UP early and DOWN late is not a reliable steady-state marker. It's informative as a **temporal** marker: UP = early stress, DOWN = late/death.

## Assessment

- **Established:** No gene in this set can be confirmed as N-specific across all stresses. The data coverage is too sparse.
- **Preliminary:** ureA (PMM0965) is the strongest N-specificity candidate — UP in 5/8 N experiments (13 timepoints), no response in carbon or coculture. But untested in 5 treatment types.
- **Preliminary:** cynA (PMM0370) is the most *sensitive* N-stress indicator (responds in 6/8 experiments, best rank 1) but is NOT specific — also responds to carbon, iron, coculture.
- **Speculative:** The previous analysis's classification of genes as "N-specific" was based on absence of evidence, not evidence of absence. Many genes were simply not tested in non-N experiments.

## Gaps and friction

1. **Gene coverage gap across treatment types** — Most N-metabolism genes have no data in iron, light, viral, P, salt experiments. This makes specificity classification impossible. Root cause: different experiments have different table scopes (all_detected vs significant_only vs filtered_subset). See gaps_and_friction.md.

## Update: table_scope reinterpretation

After investigating the non-N experiments:
- [KG] ALL non-N MED4 experiments use `significant_only`, `significant_any_timepoint`, or `filtered_subset` table scopes. None use `all_detected_genes`.
- [KG] ALL non-N experiments are microarray-based (~1696 genes ≈ full genome) or RNA-seq with significance filtering.
- [interpretation] For these table scopes, gene **absence** = **measured but not significant** (the microarray covers essentially the full genome). This is evidence of **non-response**, not "untested."
- [gap] The KG and `gene_response_profile` tool cannot distinguish "not on platform" from "measured but not significant" — both result in no expression edge. The `groups_not_known` label is technically correct but biologically misleading for experiments with `significant_only` scope.

### Corrected specificity interpretation

The `?` (not_tested) entries in the specificity table should be reinterpreted as `-` (tested, not significant) for treatment types covered by microarray experiments. The table scopes are:
- carbon_stress: `filtered_subset` (166-285 genes per condition; "genes significant in any of 4 gas conditions")
- iron_stress: `significant_any_timepoint` (112 genes)
- light_stress: `filtered_subset` (198 genes; "genes significant in any of 6 light conditions") + `significant_only` (49-51 genes)
- phosphorus_stress: `filtered_subset` (34 genes; "genes with q<0.05 at 48h")
- viral: `significant_only` (29-32 genes)
- salt_stress: `significant_only` (27 genes)

**Revised conclusion:** Most N-metabolism genes ARE N-specific — they were measured in other stresses and did not respond. The previous analysis's "N-specific" classification was likely correct, but for the wrong reason (it assumed absence = not tested, when actually absence = not significant).

The few genes that DO appear in non-N experiments (cynA in iron/carbon, urtA in light/P/viral, glnA in light/carbon) are the ones where cross-reactivity is real.

## Next

- Produce a corrected specificity table that accounts for table_scope
- Rank genes by marker quality: consistency across N experiments × specificity against other stresses
- Consider temporal dynamics: which genes are early markers (hours) vs late markers (days)?
