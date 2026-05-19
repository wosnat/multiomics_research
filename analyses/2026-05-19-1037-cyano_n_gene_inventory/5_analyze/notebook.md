# Step 5 — Analyze

## Context

Step 4 built and toy-verified the `heatmap_clustering` module (Jaccard + UPGMA + clustered heatmap). Step 5 runs it on the real 19 × 61 inventory matrix and produces the primary figure.

## What I did

```bash
uv run python 5_analyze/scripts/01_cluster_and_plot.py
```

Loads `2_kg_selection/data/05_inventory_matrix.csv` (19 × 61), runs `cluster_inventory` (Jaccard distance on binarized matrix, UPGMA on both axes), writes reordered matrix + linkage CSVs, renders the heatmap.

Layout iterations: initial version had row labels cut off (long strain names); fixed by (i) shortening row labels to strain-only (genus is in the annotation strip), (ii) moving the colorbar to the bottom of the figure to free the right margin, (iii) widening `figsize` to 26 × 10.

## Results

### Strain clustering (row leaf order, top → bottom)

| # | Strain (label in figure) | Genus | Clade | Notes |
|---|---|---|---|---|
| 1 | WH8102 | Synechococcus | 5.1 / III | marine Syn cluster |
| 2 | CC9311 | Synechococcus | 5.1 / I | marine Syn cluster |
| 3 | BL107 | Synechococcus | 5.1 / IV | marine Syn cluster |
| 4 | PCC 7002 | Synechococcus | — | non-Cyanorak; rich N machinery |
| 5 | BP-1 (Thermosyn.) | Thermosynechococcus | — | thermophile |
| 6 | WH7803 | Synechococcus | 5.1 / V | marine Syn, but separates (lacks urease) |
| 7 | PCC 7942 (elong.) | Synechococcus elongatus | — | freshwater |
| 8 | UTEX 2973 (elong.) | Synechococcus elongatus | — | freshwater (paired with PCC 7942) |
| 9 | SS120 | Prochlorococcus | LLII | LL outlier (lacks urea pathway) |
| 10 | MIT9303 | Prochlorococcus | LLIV | paired with MIT9313 |
| 11 | MIT9313 | Prochlorococcus | LLIV | paired with MIT9303 |
| 12 | NATL2A | Prochlorococcus | LLI | LLI block |
| 13 | MIT0801 | Prochlorococcus | LLI | LLI block |
| 14 | NATL1A | Prochlorococcus | LLI | LLI block |
| 15 | MIT9301 | Prochlorococcus | HLII | HLII block |
| 16 | AS9601 | Prochlorococcus | HLII | HLII block |
| 17 | MIT9312 | Prochlorococcus | HLII | HLII block |
| 18 | MED4 | Prochlorococcus | HLI | HLI pair |
| 19 | RSP50 | Prochlorococcus | HLI | HLI pair |

**Hypothesis (taxa-clustering) was supported.** The strain dendrogram cleanly separates Prochlorococcus from non-Prochlorococcus at the top split. Within each cluster, taxonomy is largely recovered: LLIV strains paired (MIT9303 + MIT9313), three LLI strains together (NATL2A + MIT0801 + NATL1A), three HLII strains together (MIT9301 + AS9601 + MIT9312), HLI pair (MED4 + RSP50). Two strains break perfect taxonomy: SS120 (LLII) isolates from the other LL Pro because of urea-pathway loss; WH7803 separates from the other marine Syn because of urease loss. Both are interpretable as gene-content-driven (not noise).

### Column clustering — broad pattern

Left side of heatmap: anchor singletons + lineage-specific genes (very sparse rows). Right side: universally-conserved core (NtcA, GlnA, GlnB, GlsF, CarA, MoeB, MetC — the 7 universal groups from QC Q1) form the dense right block. Middle columns are intermediate-prevalence genes (urea pathway, cyanate operon, nitrate/nitrite pathway, Mo cofactor).

Notable column groupings:
- urea pathway (urtABCDE + ureABCDEFG) clusters tightly — strains either have the whole operon or none of it (consistent with operon-level horizontal transfer / loss).
- nitrate/nitrite pathway (nrtP, narB, narM, nirA, focA) clusters together — present in marine Syn + PCC 7002 + Thermosyn + LLI Pro NATL2A, absent from HL Pro.
- cyanate operon (cynABDSH) clusters together but splits into Pro-lineage vs Syn-lineage subclusters (the orthology bridge correctly preserves lineage-specific paralogs).

### Inventory totals

- 1159 strain × group cells total.
- 539 present (46.5%).
- 18 cells with paralogs (copy_count ≥ 2).

### Primary figure (`01_clustered_heatmap.png`)

200 DPI PNG + PDF. Final v2 layout (after researcher feedback): row dendrogram on left | strain labels | row annotation strip (genus abbreviated PRO/SYN/PARA/PICO/THERM/SYN-el, clade) | heatmap (19 × 61, discrete colormap capped at 3) on the right. Column dendrogram on top with a pathway-category color strip beneath it (10 categories, legend on the far right). Horizontal copy-count colorbar at the bottom. Pigment and n_genes annotations were dropped per researcher feedback.

### Pathway-collapsed views (Figures 2, 3, 4)

Built by `02_pathway_views.py`. The 61 ortholog groups are collapsed into 10 pathway categories (`4_methods/n_pathway_categories.py`):

| Pathway | n_groups in inventory |
|---|---|
| urease | 10 |
| Mo cofactor | 9 |
| cyanate | 8 |
| nitrate / nitrite | 7 |
| urea uptake | 7 |
| Met biosynthesis | 5 |
| regulation | 5 |
| assimilation core | 4 |
| ammonium uptake | 3 |
| other | 3 |

**Figure 2 — Strain × Pathway (19 × 10).** `figures/02_strain_x_pathway.png`. Cell = `# present groups / total groups in pathway` (with `k/N` text overlay). Strain order inherited from Figure 1 clustering for direct cross-figure reading. Reveals at a glance:
- SS120's row is mostly zeros (only regulation, assimilation core, ammonium uptake, Met biosynthesis, other).
- PCC 7002 has near-saturation across most pathways.
- *S. elongatus* PCC 7942 and UTEX 2973 are visually identical.

**Figure 3 — Clade-group × Pathway summary (8 × 10).** `figures/03_cladegroup_x_pathway.png`. Strains aggregated into 8 clade groups; cell = mean fraction across strains in the group. Reveals clade-level distinctions:
- **All Pro clades except LLII (= SS120) retain ~70% of urea uptake + ~70% of urease.**
- **Cyanate is HLI Pro-specific** (38% in HLI; 0% in HLII, LLII, LLIV; 12% in LLI from MIT0801 alone) — consistent with the MED4-lineage cyanate operon (cynABDS at PMM0370–0373).
- **Nitrate/nitrite is mostly Syn + Thermosyn**: marine Syn 68%, non-marine Syn 52%, Thermosyn 43%; absent in HL Pro; partial in LL Pro (29%, via NATL2A nirA + MIT9303/MIT9313 nirA).
- **Mo cofactor biosynthesis is a Syn/Thermosyn signature**: marine Syn 92%, non-marine Syn 59%, Thermosyn 56%, vs Pro 11% everywhere.
- **Met biosynthesis** is 100% in marine Syn, 80% in all Pro, 47% in non-marine Syn (S. elongatus lacks the metB-like paralog), 20% in Thermosyn (Thermosynechococcus uses a different Met biosynthesis route — flagged in step-2 surprises).

**Figure 4 — Pathway composition stacked bars per strain.** `figures/04_strain_pathway_composition.png`. Horizontal stacked bars, strains ordered by total N-group count descending. Pathway colors match Figure 1's top strip. Visualises N-machinery "investment" per strain: CC9311 (43 groups, top) and WH8102 (41) are the richest; SS120 (12, bottom) is sparest. Reveals that **strains with higher total N-group counts diversify across more pathways** — not just more copies of the same pathway.

## Surprises

- **The data-driven clustering recovered taxonomy almost exactly.** Both deviations (SS120 in Pro; WH7803 in marine Syn) are explained by single-pathway losses, not noise. Strong validation that the inventory captures real biology.
- **PCC 7002 clusters with the freshwater/non-marine group, not with marine Synechococcus** — despite being a *Synechococcus* by genus. Its N-machinery profile (rich nitrate + nitrite + urea + cyanate + Mo cofactor + ntrB/ntrX regulators + amt1 × 3) is more like the *S. elongatus* / Thermosynechococcus pattern than the streamlined marine *Synechococcus*. Suggests N-uptake is plastic across the *Synechococcus* genus and partly reflects habitat / lifestyle.
- **MIT9303 ntcA × 2 (paralog) is visible in the figure** as a darker cell in MIT9303's row at the ntcA column. The Q6 paralog finding from step 2 is now confirmed visually in the data-driven layout.

## Decisions

**2026-05-19 — Primary figure ships with horizontal colorbar at bottom.** Iterated through 3 layout versions to fit the strain labels; final layout puts colorbar below the heatmap (gs row 3) with extended `hspace=0.55` so column labels don't clip into the colorbar.

**2026-05-19 — Strain labels are strain-only (no "Pro_"/"Syn_" prefix).** The genus is already in the annotation strip; prefixes were redundant and consumed right-margin space.

## Advance rationale

Primary figure produced; structural hypothesis confirmed (taxa-clustering); deviations are explainable by specific gene-content losses. Step 6 can evaluate against the framing and write the Discussion. The secondary figure (A1 + B1 + C2 view) is still planned but deferred until after this primary figure is reviewed.

---

## Decide-gate checklist

- **Outputs produced.**
  - Scripts: `5_analyze/scripts/01_cluster_and_plot.py`, `5_analyze/scripts/02_pathway_views.py`
  - Data: `01_clustered_matrix.csv` (19 × 61), `01_row_linkage.csv`, `01_col_linkage.csv`, `01_row_order.txt`, `01_col_order.txt`, `02_strain_pathway_count.csv` (19 × 10), `02_strain_pathway_pct.csv` (19 × 10), `03_cladegroup_pathway_pct.csv` (8 × 10). Plus 2 log files.
  - Figures: `01_clustered_heatmap.png/.pdf` (Figure 1, primary clustered heatmap), `02_strain_x_pathway.png/.pdf` (Figure 2, strain × pathway collapsed view), `03_cladegroup_x_pathway.png/.pdf` (Figure 3, clade-group × pathway summary), `04_strain_pathway_composition.png/.pdf` (Figure 4, per-strain pathway composition stacked bars).
- **Results presented.** Strain leaf order table, column-clustering broad pattern, inventory totals, and the figure description all shown inline above.
- **QC gate.**
  - Clustering pipeline ran on real data without warning → PASS.
  - Strain dendrogram shows expected Pro/Syn split and within-genus clade-level grouping → PASS (and validates the inventory).
  - Two clustering deviations (SS120 within Pro; WH7803 within marine Syn) are explained by single-pathway gene-content losses, not noise → PASS.
  - Inventory total cells = 19 × 61 = 1159; 539 present (46.5%); 18 with paralogs → matches step 2 numbers.
- **Decisions made this step.** Two: primary figure layout (horizontal bottom colorbar); strain label form (strain-only).
- **Advance rationale.** Figure produced and biologically validated. Step 6 (evaluate + Discussion) is the natural next step.
