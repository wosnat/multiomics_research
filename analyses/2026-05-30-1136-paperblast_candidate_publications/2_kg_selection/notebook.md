# Step 2 — KG entries: the prioritized seed-gene list

## Context

Step 1 locked the workflow: build a prioritized gene list (this step), the
researcher runs PaperBLAST offline (step 3), Claude analyzes the downloads into
a candidate-publication shortlist (later). This step produces the gene list.

Two project goals were separated here (the researcher chose to pursue them one at
a time): **(1) more papers** — rank toward genes most likely to return PaperBLAST
hits; **(2) info on low-info genes** — rank toward responding hypotheticals. This
step's deliverable is ranked for **goal 1**; the same pipeline produces a goal-2
ranking (`--mode jackpot`) when wanted. Pools are kept **separate** (no cross-pool
ortholog dedup), and **4 pools** are required: `pro`, `syn`, `alt`, `other_hetero`.

## What I did

A 5-stage pipeline, all typed `multiomics_explorer` tools (no raw Cypher), run
from the repo root with `env -u VIRTUAL_ENV uv run python <script>`:

1. **`top_genes.py`** — per-experiment top responders. A gene is a top responder
   of an experiment if its best (lowest) KG-supplied rank among significant rows
   is `<= TOP_N (=10)` in **any** timepoint, **either** direction (`rank_up` or
   `rank_down`, 1 = strongest |log2FC| within experiment×timepoint). **Conditional
   timepoint rule:** time-course experiments (≥2 timepoints) require the gene to
   be top-10 in `>= 2` timepoints; single-contrast experiments require 1. This
   drops transient single-timepoint spikes without zeroing out single-contrast
   experiments. Verified on the Tolonen MED4 N-deprivation time-course (50→21
   genes under the rule) and a single-contrast cyanate-growth experiment (stayed
   at 20). `to_dataframe()` + pandas does the timepoint collapse.
2. **`02_all_experiments.py`** — pool top responders across **all 197 experiments**
   (no cross-experiment dedup yet). Per-organism / per-genus distribution.
3. **`03_dedup_and_enrich.py`** — dedup to one row per locus_tag (freq = #
   experiments where top responder); enrich with `gene_overview` (annotation_state,
   closest_ortholog_group_size/genera) and `gene_aa_sequence` (protein_id).
4. **`04_ortholog_dedup.py`** — assign 4 pools by genus; attach the **finest-
   available ortholog group** as dedup key (cyanorak curated rank-0 when present,
   else finest eggNOG family→phylum→COG, else locus_tag); dedup **within pool**.
   Representative of a collapsed group = **best BLAST odds** (largest
   ortholog_group_size, then most genera, then freq, then rank). The group's
   response score aggregates across all members (freq=max, best_rank=min).
5. **`05_rank.py --mode papers --top 8 --cat-cap 2 --org-cap 3`** — score and
   select. `score = response × conservation` (info_gain dropped for goal 1).
   `response = 0.6·rank_norm + 0.4·freq_norm`; `conservation = 0.5·log(og_size) +
   0.5·n_genera` (PaperBLAST hit-likelihood proxy). Select top 8 per pool with
   soft caps: ≤2 per gene_category (multi-pathway breadth) and ≤3 per organism
   (so no single strain fills a pool).

## Results

### Funnel
| stage | count |
|---|---|
| experiments processed | 155 contributed / 42 errored (package bug) / 197 total |
| pooled top-responder rows (before dedup) | 2290 |
| distinct genes after locus_tag dedup | 1653 |
| seeds after within-pool ortholog dedup | **1430** (pro 575 / syn 308 / alt 282 / other_hetero 265) |
| **round-one selection** | **32** (8 per pool) |

The 42 errored experiments fail inside the package's diagnostic-query (an upstream
bug, not our call); they are mostly metabolomics / vesicle / exoproteomics, which
carry little or no gene-level DE anyway. Logged in `gaps_and_friction.md`.

### Distribution before dedup (by genus)
| genus | pooled rows | distinct genes | organisms | experiments |
|---|---|---|---|---|
| Prochlorococcus | 1070 | 726 | 9 | 68 |
| Alteromonas | 461 | 320 | 4 | 35 |
| Synechococcus | 252 | 218 | 5 | 17 |
| Marinobacter | 227 | 148 | 1 | 18 |
| Parasynechococcus | 100 | 90 | 1 | 6 |
| Ruegeria | 80 | 68 | 1 | 3 |
| Pseudomonas | 40 | 38 | 1 | 2 |
| Picosynechococcus | 30 | 22 | 1 | 3 |
| Shewanella | 30 | 23 | 1 | 3 |

(`Parasynechococcus` / `Picosynechococcus` are GTDB reclassifications of marine
*Synechococcus*; folded into the `syn` pool. Full table: `data/02_distribution_by_organism.csv`.)

### Round-one selected seeds (32 = 8 per pool) — papers mode

**pro** (5 organisms)
| # | locus | organism | gene | product | freq | og_size |
|---|---|---|---|---|---|---|
| 1 | PMM1028 | MED4 | – | uncharacterized conserved secreted | 9 | 28 |
| 2 | PMM1404 | MED4 | hli | high light inducible protein | 10 | 40 |
| 3 | PMM1135 | MED4 | hli | high light inducible protein | 10 | 31 |
| 4 | PMT9312_0550 | MIT9312 | rbcL | RuBisCO large subunit | 4 | 20 |
| 5 | PMT9312_0549 | MIT9312 | csoS1A | carboxysome shell protein | 3 | 20 |
| 6 | PMT1742 | MIT9313 | rplN | 50S ribosomal protein L14 | 2 | 20 |
| 7 | Pro0418 | marinus | mqoA | malate:quinone oxidoreductase | 2 | 20 |
| 8 | NATL1_02111 | NATL1A | tsaB | tRNA threonylcarbamoyladenosine | 1 | 20 |

**syn** (4 organisms): resP, outer-membrane porin, nusG, cyabrB2, rpoB, nucleoside 2-deoxyribosyltransferase, argB, psbO
**alt** (3 organisms): IS110 transposase, acnB, adk, flaA, cheW, fecA, slyA, lldP
**other_hetero** (4 organisms): gtsA (Marinobacter), clpB (Marinobacter), prlC (Marinobacter), SPO0362 (Ruegeria), idh (Pseudomonas), sucC (Ruegeria), tolB (Ruegeria), cytochrome-c biogenesis (Shewanella)

Full ranking: `data/05_ranked_seeds_papers.csv` (1430 seeds). Round-one set:
`data/05_selected_papers.csv` (32). PaperBLAST query = `rep_protein_id` (RefSeq).

## Surprises

- **Sanity check passes:** the Tolonen-only run surfaced the N-stress canon —
  cynA (rank 1), glnA, urtA, ntcA induced; ATP-synthase / RuBisCO / ribosomal
  proteins repressed. The pipeline finds real signal.
- **MED4 dominates the raw pool** (524/2290 rows, 34 experiments) — the per-organism
  cap is what keeps `pro` from being all-MED4 and `other_hetero` from being
  all-Marinobacter.
- **A family-level second ortholog dedup is a near-no-op** (1430→1418); the only
  meaningful broader collapse is COG/Bacteria level, and its effect is *cross-pool*
  (1388→1108 global) — e.g. amt1/amtB, proA, purB, psbA recur across pools. We
  deliberately did **not** apply it: keeping pools separate maximizes paper
  discovery (a Pro and an Alteromonas homolog can yield different papers).
- **HLIP paralog family** is the biggest single within-pool collapse (7 hli →
  1 via cyanorak curated); the category cap further limits HLIPs to 2 in pro's top-8.

## Decisions

**2026-05-30 — per-experiment cut: top-10, both directions, ≥2 timepoints (conditional).**
Chosen over top-5 and ≥1-tp after comparing set sizes on Tolonen. Conditional on
time-course vs single-contrast so single-contrast experiments aren't zeroed.

**2026-05-30 — 4 pools (pro/syn/alt/other_hetero), kept separate.** No cross-pool
ortholog dedup. Marine *Synechococcus* GTDB-renames folded into `syn`.

**2026-05-30 — ortholog dedup key = finest available, uniform across pools.**
cyanorak curated rank-0 when present (covers Pro/Syn), else finest eggNOG
(family→phylum→COG), else locus_tag. Heterotrophs have no cyanorak, so eggNOG
family is their finest. Representative = best BLAST odds.

**2026-05-30 — ranked for goal 1 (more papers): score = response × conservation.**
info_gain dropped (it was pulling toward hypotheticals, which serves goal 2).
Goal-2 ranking available via `--mode jackpot`.

**2026-05-30 — round one = 8 per pool (32), caps: ≤2/category, ≤3/organism.**

## Advance rationale

The 32-seed round-one list is produced, diversified across organisms and pathways,
and ranked toward PaperBLAST hit-likelihood. The pipeline is reproducible (one
script per stage) and parameterized (top-n, caps, mode) so it re-runs for goal 2
or larger batches. Ready for step 3: the researcher runs these 32 protein IDs
through PaperBLAST offline and saves the result pages.

---

## Decide-gate checklist

- **Outputs produced.**
  - `scripts/top_genes.py` (per-experiment top responders)
  - `scripts/02_all_experiments.py` → `data/02_pooled_top_responders.csv`, `data/02_distribution_by_{organism,genus}.csv`, `data/02_all_experiments.log`
  - `scripts/03_dedup_and_enrich.py` → `data/03_genes_dedup.csv`
  - `scripts/04_ortholog_dedup.py` → `data/04_genes_ortholog.csv`, `data/04_pool_representatives.csv`, `data/04_ortholog_dedup.log`
  - `scripts/05_rank.py` → `data/05_ranked_seeds_papers.csv` (1430), `data/05_selected_papers.csv` (32)
  - `scripts/qc_knobs.py`, `scripts/qc_second_dedup.py` (exploratory checks)
  - Run: `env -u VIRTUAL_ENV uv run python 02_all_experiments.py --top-n 10 --min-tps 2`; then `03`, `04`; then `05_rank.py --mode papers --top 8 --cat-cap 2 --org-cap 3`.
- **Results presented.** Funnel, before-dedup genus distribution, and the 32-seed selection shown inline above (same tables shown in chat).
- **QC gate.**
  - Per-experiment cut compared (top-5/10 × ≥1/2 tp) on Tolonen → top-10/≥2tp chosen.
  - Method sanity-checked on Tolonen N-deprivation → N-stress canon recovered.
  - Ortholog dedup verified collapsing HLIP/operon paralogs correctly.
  - Second-dedup levels probed (`qc_second_dedup.py`) → family no-op, COG cross-pool only → not applied.
  - Category + organism caps verified to hold in the 32-seed output.
- **Decisions made this step.** Five (logged above): per-experiment cut; 4 separate pools; finest-available ortholog dedup key + best-BLAST-odds rep; goal-1 ranking; round-one size + caps.
- **Advance rationale.** 32-seed list produced, diversified, reproducible; ready for the manual PaperBLAST step.
