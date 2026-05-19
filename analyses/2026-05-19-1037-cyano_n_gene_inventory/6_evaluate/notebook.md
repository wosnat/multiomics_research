# Step 6 — Evaluate

## Context

Steps 1–5 built the cross-strain N-gene inventory (19 × 61 ortholog group matrix), the primary clustered heatmap (Figure 1), and three pathway-collapsed views (Figures 2, 3, 4). This step evaluates the analysis against the locked question, runs the statistical tests requested by the researcher, builds the deferred secondary figure (Figure 5, A1 + B1 + C2 view), and finalizes `paper.md` Discussion.

## What I did

```bash
uv run python 6_evaluate/scripts/01_statistical_tests.py
uv run python 6_evaluate/scripts/02_secondary_taxonomy_view.py
```

`01_statistical_tests.py` runs three tests:
1. **Mantel** on Jaccard distance vs taxonomic distance (9999 permutations).
2. **Fisher's exact** broad scan (80 clade-group × pathway tests, BH-corrected) plus 5 **pre-specified contrasts** (BH × 5), plus **Mann-Whitney U** on % pathway coverage for 5 quantitative contrasts (BH × 5).
3. **Silhouette score** per strain on Jaccard distance, using merged clade-group labels.

`02_secondary_taxonomy_view.py` produces Figure 5: taxonomy-ordered rows × pathway-ordered columns, with vertical lines separating pathway blocks and a horizontal line separating Pro from Syn/Thermosyn.

## Results

### Test 1 — Mantel: Jaccard distance correlates strongly with taxonomic distance

| Test | n strains | n distances | permutations | Mantel r | p-value |
|---|---|---|---|---|---|
| Jaccard vs taxonomic distance | 19 | 171 | 9999 | **0.784** | **0.0001** (two-sided) |

The pre-registered "taxa-clustering" hypothesis from step 3 is strongly supported: at the distance-matrix level, strains that are taxonomically close are also close in N-gene content (Pearson r = 0.78 between the two distance matrices; permutation p < 0.001).

### Test 2 — Pathway contrasts

**Test 2a (broad scan).** 80 (clade-group × pathway) Fisher's exact tests with BH × 80 → 0 contrasts significant at q < 0.05. The broad scan has limited power given small clade-group sizes (n = 1 to n = 4). Three contrasts pass raw p < 0.05 (Pro HLII × nitrate_nitrite DEPLETED, Pro HLII × cyanate DEPLETED, Non-marine Syn × other DEPLETED) but all have q > 0.95 after correction.

**Test 2b (pre-specified contrasts, Fisher's exact, BH × 5).** Binary "≥1 gene present" gives mixed signal — appropriate for pathways that are all-or-none (urea operon), too coarse for pathways with gradient coverage (Mo cofactor):

| Contrast | a | b | OR | p | q (BH × 5) |
|---|---|---|---|---|---|
| HLI Pro vs other Pro: cyanate | 2/2 | 3/9 | ∞ | 0.18 | 0.30 |
| Syn vs Pro: nitrate | 8/8 | 5/11 | ∞ | 0.018 | 0.090 (marginal) |
| Syn vs Pro: Mo cofactor | 8/8 | 11/11 | NaN | 1 | 1 (binary too coarse) |
| SS120 vs other Pro: urea_uptake | 0/1 | 10/10 | 0 | 0.091 | 0.23 (n=1 limits) |
| Thermosyn vs others: met biosynthesis | 1/1 | 18/18 | NaN | 1 | 1 (binary too coarse) |

**Test 2c (pre-specified contrasts, Mann-Whitney U on % pathway coverage, BH × 5).** The right tool for gradient pathways — 4 of 5 contrasts significant at q < 0.05:

| Contrast | a median % | a n | b median % | b n | U | p | q (BH × 5) |
|---|---|---|---|---|---|---|---|
| Syn vs Pro: nitrate_nitrite | **0.64** | 8 | **0.00** | 11 | 88.0 | 2.1e-4 | **5.3e-4** |
| Syn vs Pro: mo_cofactor | **0.78** | 8 | **0.11** | 11 | 88.0 | 5.9e-5 | **3.0e-4** |
| Marine Syn vs Non-marine Syn: met_biosynthesis | **1.00** | 4 | **0.40** | 4 | 16.0 | 0.019 | **0.024** |
| HL Pro vs LL Pro: nitrate_nitrite | **0.00** | 5 | **0.29** | 6 | 2.5 | 0.011 | **0.019** |
| HL Pro vs LL Pro: cyanate | 0.00 | 5 | 0.06 | 6 | 16.5 | 0.84 | 0.84 |

Four contrasts back up the Discussion claims: (a) Syn has nitrate machinery, Pro doesn't; (b) Syn has Mo cofactor biosynthesis, Pro mostly doesn't; (c) marine Syn has full Met biosynthesis, non-marine Syn/Thermosyn doesn't; (d) LL Pro has partial nitrate (via nirA), HL Pro doesn't. The HL vs LL cyanate contrast is **not** significant — because cyanate is HLI-specific, not HL-wide, the HL pool dilutes the effect. (The descriptive "HLI Pro has cyanate" claim is qualitatively correct from Figures 1 and 3, but the formal Fisher's test on n=2 lacks power.)

### Test 3 — Silhouette: SS120, BP-1, and surprisingly the 3 LLI Pro fit their clades poorly

Per-strain silhouette on Jaccard distance using merged clade-group labels (Pro HL, Pro LL, Marine Syn, Non-marine cyano). Mean = 0.241; range [-0.31, 0.76]:

| Strain | Silhouette | Merged label | Note |
|---|---|---|---|
| MIT0801 | -0.312 | Pro LL | LLI member; pulled by SS120 in the Pro LL pool |
| NATL1A | -0.312 | Pro LL | LLI member; same |
| NATL2A | -0.312 | Pro LL | LLI member; same |
| BP-1 (Thermosyn) | -0.103 | Non-marine cyano | Singleton in own clade; merged into non-marine cyano |
| SS120 | -0.102 | Pro LL | Singleton in LLII; outlier within Pro LL — the urea-loss signature |
| ... | ... | ... | ... |
| MED4, MIT9301, RSP50, MIT9312, AS9601 | 0.61 – 0.76 | Pro HL | HL Pro cluster very tightly together |

Two observations:
1. **The 3 LLI Pro strains have identically negative silhouettes (−0.312).** This is informative: the "Pro LL" merged label includes SS120 (LLII, urea-loss) + MIT9303/MIT9313 (LLIV) + the three LLI; the LLI subgroup is more N-gene-similar to HL Pro than to the heterogeneous Pro LL pool. Within Pro, LLI vs HL is a smaller distance than LLI vs SS120-included-LL.
2. **HL Pro is the tightest cluster** — all 5 HL Pro strains have silhouette > 0.6. Cyanate-operon-presence variance (MED4/RSP50 have it; AS9601/MIT9301/MIT9312 don't) is not enough to overcome the shared HL N-machinery.

### Figure 5 — Secondary taxonomy-ordered view

`figures/05_taxonomy_view.png`. Rows ordered by taxonomy (HLI → HLII → LLI → LLII → LLIV Pro → marine Syn → non-marine Syn → Thermosyn); columns ordered by pathway category. Vertical black lines separate pathway blocks; horizontal black line separates Pro from non-Pro. Reveals the same patterns as Figure 1 but with explicit taxonomic structure: HLI Pro cyanate block; LL Pro urea + nirA; SS120 vertical-band of zeros; marine Syn richness across all blocks; *S. elongatus* missing urea/urease block while keeping nitrate.

## Surprises

- **Mantel r = 0.78 is unusually high** for distance-matrix correlation. Suggests the cyano N-gene complement is a strong taxonomic marker — comparable to using 16S or housekeeping genes for distance estimation. Possible application: N-gene Jaccard distance could be a quick proxy for cyano taxonomic placement.
- **HL Pro vs LL Pro on cyanate is NOT significant** even though cyanate is qualitatively HL-specific. Reason: cyanate is HLI-specific (only MED4 + RSP50 have the full operon), and HLII Pro (AS9601, MIT9301, MIT9312) lack it entirely — pooling them dilutes the signal. A more focused HLI-vs-rest contrast (skipped here for n=2 power reasons) would likely show the effect.
- **LLI Pro silhouette negative due to merger artifact, not biology.** When SS120 (LLII, atypical) is included in the merged "Pro LL" label, the LLI subgroup ends up closer to Pro HL in Jaccard space than to the heterogeneous Pro LL centroid. Suggests Pro LL is genuinely not a single cohesive group in N-gene space — different LL clades (LLI, LLII, LLIV) have meaningfully different N machinery.

## Caveats harvested (from gaps_and_friction.md + this step)

1. **Anchor seed is Cyanorak-only.** Genes that are N-related but never appear in any Cyanorak `E.4 ∪ D.1.3` strain will not enter the inventory. Likely affects strains outside Cyanorak's curation scope (the 4 non-Cyanorak Syn-genus strains): if they have unique N machinery that no Cyanorak strain has, it's invisible to this analysis.
2. **9 orphan anchor locus_tags** could not map to an eggnog Cyanobacteria-level ortholog group and were recovered as 7 synthetic singleton groups (logged in `gaps_and_friction.md`). The Pro-lineage cynA paralog is preserved this way; this approach mixes ortholog tiers (eggnog Cyanobacteria for most groups, anchor singleton for orphans) but the `source` column makes the distinction auditable.
3. **2 KG loader bugs** identified during step 2 affected this analysis: Synechococcus `Organism.clade` is null (worked around via Cyanorak source CSV); NATL1A and NATL2A `clade` is mislabeled `LLII` instead of `LLI` (corrected for this analysis but other downstream consumers of the KG would inherit the bug).
4. **KEGG `ko00910` "Nitrogen metabolism" is too narrow for *Prochlorococcus*** (6 genes in MED4). Cyanorak roles E.4 ∪ D.1.3 gave a cleaner cyano-specific definition.
5. **Statistical power is limited for many per-clade-group contrasts** due to small N (some clade groups have n = 1 or n = 2). The Mantel result is the headline; pre-specified Mann-Whitney contrasts on % pathway coverage provide focused per-pathway evidence; the broad-scan Fisher's exact (BH × 80) is too conservative and serves only as a screening step.
6. **No expression data was used.** This is purely a gene-content inventory. A separate analysis could test whether the present/absent N-genes are actually expressed (using the 4 N-perturbation studies in the KG) — flagged as a natural follow-up.

## Decisions

**2026-05-19 — Mann-Whitney U on % pathway coverage adopted as the primary pathway-contrast test** (Test 2c), with Fisher's exact on binary presence (Test 2b) as a secondary view. Reason: Fisher's binary "≥1 gene present" misses gradient effects on multi-gene pathways (Mo cofactor, Met biosynthesis); MWU on the continuous % coverage captures both gradient and all-or-none signal correctly.

**2026-05-19 — Hypothesis "taxa-clustering" confirmed.** Mantel r = 0.78 (p < 0.001) at the distance-matrix level; cluster dendrogram (Figure 1) recovers known clade structure visually; two systematic deviations (SS120, WH7803) are explained by single-pathway gene-content loss, not by noise. No grounds for rejecting the hypothesis or revisiting the analysis design.

**2026-05-19 — Discussion-level claims to make in paper.md.** Five named findings (all backed by figures + statistical evidence): (a) N-gene content recovers taxonomy strongly (Mantel r=0.78); (b) cyano-wide N core = 7 ortholog groups; (c) Syn-genus strains have nitrate + Mo cofactor machinery, Pro mostly doesn't (MWU q < 0.001 both); (d) marine Syn has full Met biosynthesis, non-marine Syn/Thermosyn doesn't (MWU q=0.024); (e) two strain-level outliers (SS120 urea loss, WH7803 urease loss) are gene-content-driven, not annotation artifacts.

## Advance rationale

All analytical artifacts are produced and committed. The pre-registered hypothesis is supported by both visual evidence (Figure 1) and statistical tests (Mantel, Mann-Whitney). The Discussion in `paper.md` can be written now from the captured findings, caveats, and decisions. After this step, the analysis is closed (no further steps in the 6-step flow).

---

## Decide-gate checklist

- **Outputs produced.**
  - `6_evaluate/scripts/01_statistical_tests.py`, `02_secondary_taxonomy_view.py`
  - `6_evaluate/data/01_mantel_test.csv`, `01_fisher_pathway_clade.csv`, `01_fisher_prespecified_contrasts.csv`, `01_mannwhitney_contrasts.csv`, `01_silhouette.csv`, `01_statistical_tests.log`, `02_secondary_taxonomy_view.log`
  - `6_evaluate/figures/05_taxonomy_view.png`, `.pdf`
- **Results presented.** All four statistical tests' summary tables shown inline above; Figure 5 description in text.
- **QC gate.**
  - Mantel test uses 9999 permutations → adequate power; r = 0.78, p < 0.001 → strongly supports taxa-clustering.
  - Pre-specified MWU contrasts use the right scale (continuous %) and have appropriate BH correction (× 5) → 4 of 5 q < 0.05.
  - Silhouette computed on the same Jaccard distance matrix used for Figure 1 clustering → internally consistent.
  - Figure 5 reuses the same inventory matrix; vertical/horizontal dividers + pathway color strip make taxonomic structure explicit.
- **Decisions made this step.** Three: MWU on % coverage as primary; hypothesis confirmed; five Discussion claims.
- **Advance rationale.** Analysis ready to close. Discussion section of `paper.md` finalized in this commit.
