# Step 2 — Is the up-regulation specific to N stress? (pooled, 3 families)

## Question
Step-1 flowers showed these genes go UP under nitrogen stress, but they also move
under other conditions. **Is up-regulation enriched under N relative to other
stresses — and is that true for all three families, or only some?**

## Method (locked framing)
- **Self-contained, reproducible:** `scripts/pooled_n_specificity.py` downloads the
  data from the KG via the `multiomics_explorer` Python API
  (`differential_expression_by_ortholog` → `to_dataframe`), caches it to
  `data/pooled_de_downloaded.csv`, then computes the test and figure. Rerun:
  `uv run python analyses/.../scripts/pooled_n_specificity.py`.
- **Three families** (pooled across all genomes with DE data):
  - `urt*` = urea transporter, ortholog groups urtA–E (CK_00000076/1365/1366/1367/8074)
  - `amt*` = ammonium transporter, amt1 + amt2 (CK_00000244, CK_00008701)
  - `cynS` = cyanate hydratase (CK_00001552)
- **Test:** per family, Fisher's exact on `(experiment is N vs non-N) × (family UP vs not-up)`.
  - Unit = experiment (each experiment = one genome). Timepoints **and** subunits
    collapsed: family is "up" if **any** member is significantly up in **any** timepoint;
    "tested" if any member is detected.
  - **DE / not-DE call only — no fold-change magnitude** → robust to platform
    dynamic-range differences (RNA-seq / microarray / proteomics). [feedback: rank or DE/not-DE, not FC]
  - One-sided (`greater`) = up enriched under N. BH-FDR across the 3 families.
- Scope: **pooled only** (the MED4 per-gene panel was dropped — underpowered).

## Results [KG]  (`data/n_specificity_pooled_results.csv`, figure `figures/fig8_n_specificity_forest.*`)
Downloaded: 421 rows, 7 genomes, 41 experiments.

| family | up/tested (N) | up/tested (other) | OR | q (BH) |
|---|---|---|---|---|
| **urt\*** (urea) | **10/11** | **7/27** | **28.6** | **0.0010 \*** |
| amt\* (ammonium) | 6/11 | 6/17 | 2.2 | 0.27 |
| cynS (cyanate) | 4/8 | 1/9 | 8.0 | 0.17 |

## Interpretation
**Plain English:** Pooled across genomes, **urea-transporter (urt\*) up-regulation is
strongly and significantly nitrogen-specific** — up in 10 of 11 N experiments vs only
7 of 27 non-N (odds ratio ≈ 29, q = 0.001). **Cyanate hydratase (cynS)** points the same
way (4/8 vs 1/9, OR 8) but does **not** reach significance (q = 0.17, underpowered — only
8 N experiments test it). **Ammonium transporter (amt\*) is NOT N-specific** (6/11 vs 6/17,
OR 2.2, q = 0.27): it goes up about as readily under other conditions.

`[interpretation]` The three families behave differently, which sharpens the marker story:
- **urt\*** — the robust N-stress reporter in this dataset.
- **cynS** — probably N-specific (large OR) but **under-measured**; needs more N experiments
  across strains to confirm.
- **amt\*** — a **general-response** gene, not a clean N marker (ammonium is the *preferred*
  N source, so amt regulation tracks many physiological states, not N starvation alone).

`[interpretation]` Net: "are these up-regulations specific for N stress?" → **urea transport, yes;
cyanate, likely; ammonium transport, no.** A urt-anchored panel is the defensible N-stress detector.

## Caveats `[KG]` / `[interpretation]`
- **Family-level "any-subunit-up"** — urt\* has 5 subunits, so "any member up" is more permissive
  than a single gene (raises the up-rate in both N and non-N; the OR is still high, but the
  absolute rates should be read with this in mind).
- **Pseudoreplication** — timepoints collapsed to experiment; experiments within one study are
  not fully independent (no study-level random effect).
- **Per-publication thresholds** — "significant" uses each study's own DE call (this is the point —
  platform-robust — but the unit of "up" varies by study).
- **tested-absent / table_scope** — for `significant_only` / `top_n` source tables, a gene is only
  recorded when significant, so "tested-but-not-up" is undercounted there; can inflate up-rates.
- **"non-N" is a heterogeneous union** (carbon, P, light, viral, coculture, salt, iron, darkness,
  plastic) — the test asks "N vs everything else", not "N vs any one specific stress".
- **Direction = up only** (down and n.s. both count as "not-up").
- **cynS only** — the Synechococcus-only paralog group cynH (CK_00003051) is excluded; it has no DE data anyway.

## Decision
Slide headline: **"Up-regulation is significantly N-specific for urea transport (urt\*: 10/11 vs
7/27, OR≈29, q=0.001); cyanate (cynS) trends the same; ammonium transport (amt\*) is not N-specific."**
Use `fig8` alongside the flowers (fig6/fig7).
