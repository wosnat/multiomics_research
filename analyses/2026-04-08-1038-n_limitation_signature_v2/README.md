# N-Limitation Signature Analysis v2

**Question:** Can we quantify nitrogen limitation in Prochlorococcus MED4 molecularly, using a gene signature from independent reference studies?

**Approach:** Methodological rebuild of v1 with proper notebook discipline, single primary metric (rank score normalized by n_significant), positive/negative controls, three signature tiers (top/core/extended), and package-ready sig_utils.

## Key findings

- **Core signature: 189 genes** (74 up, 115 down) from concordant intersection of Tolonen 2006 (microarray, 6h-48h) and Read 2017 (RNA-seq, 12h-24h). Excludes early transient timepoints from both studies.

- **RNA-seq axenic: strongly N-limited.** Day 14 rank_score=0.583, hit_rate=0.85, p=0.000. 158 of 185 matched genes concordant.

- **RNA-seq coculture: no transcriptomic signal.** All timepoints score near zero (p>0.01 at most timepoints). Alteromonas eliminates MED4 transcriptional N-stress response.

- **Proteomics coculture: persistent protein-level signal.** Significant at all timepoints (p=0.000), peak at day 31 (score=0.217, 53 concordant, 1 reversed).

- **RNA/protein discordance is genuine.** Tested with 147 genes detectable on both platforms: RNA-seq scores 0.548 vs proteomics 0.066 for the same gene set. Not a coverage artifact.

- **Proteomics missing genes are biased.** 42 of 189 core genes undetected — predominantly small/membrane proteins (hli, photosystem subunits, ribosomal) with better-than-average signature ranks.

## Changes from v1

| Aspect | v1 | v2 |
|--------|----|----|
| Core genes | 198 | 189 |
| Metrics | 3 (hit rate, log2FC, rank score) | 1 (rank score) + hit rate as context |
| Rank normalization | by total_genes (~0.95 for all) | by n_significant (meaningful spread) |
| Read 3h | included | excluded (transient up-bias) |
| Tie-breaking | pandas idxmax (arbitrary) | best directional rank |
| Extended classification | asymmetric (B-only unclassified) | symmetric (_ns/_absent both ways) |
| Discordant genes | counted but lost | saved to file |
| Controls | none | 4 negative controls + reference self-scoring |
| Scoring tiers | core + extended | top (38) + core (189) + extended (531) |
| Toy tests | 0 | 31 |

## File index

- `methods.md` — publication-ready methods
- `decisions.md` — design decisions with rationale
- `caveats.md` — interpretation caveats for readers
- `gaps_and_friction.md` — KG/methodology issues
- `references.md` — citations and versions
- `data/` — staged KG extracts (see `data/DATA_MANIFEST.md`)
- `scripts/` — per-step pipeline scripts (01-06)
- `sig_utils/` — reusable scoring utilities (extraction, signature, scoring, I/O, tests)
- `results/` — score tables and figures (see `results/RESULTS_MANIFEST.md`)
- `exploration/` — research notebook
- `superpowers/` — spec, plan, brainstorm log

## How to reproduce

Requires Neo4j running + multiomics_explorer installed. Run from `multiomics_research` root.

```bash
# Step 1: Experiment scoping (was done interactively — CSV tracked in git)

# Step 2: Extract reference and control DE data (large CSVs, not tracked)
uv run analyses/2026-04-08-1038-n_limitation_signature_v2/scripts/02_extract_reference_de.py

# Step 3: Build signature (produces core/extended/discordant CSVs)
uv run analyses/2026-04-08-1038-n_limitation_signature_v2/scripts/03_build_signature.py

# Step 4: Extract target DE data (large CSVs, not tracked)
uv run analyses/2026-04-08-1038-n_limitation_signature_v2/scripts/04_extract_target_de.py

# Step 5: Score all experiments (produces scores_all.csv + applied subsets)
uv run analyses/2026-04-08-1038-n_limitation_signature_v2/scripts/05_score_experiments.py

# Step 6: Generate plots
uv run analyses/2026-04-08-1038-n_limitation_signature_v2/scripts/06_plot_results.py

# Run tests
uv run python -m pytest analyses/2026-04-08-1038-n_limitation_signature_v2/sig_utils/tests/ -v
```

Steps 2 and 4 produce large DE extract CSVs (`de_*.csv`) that are not tracked in git. All other outputs (signatures, scores, plots, logs) are tracked. Re-run steps 2+4 first to regenerate the extracts, then steps 3-6 will work.

## Spec and plan

- **Spec:** `superpowers/spec.md`
- **Plan:** `superpowers/plan.md`
- **Brainstorm:** `superpowers/brainstorm-log.md`
- **Predecessor:** `analyses/2026-04-06-1432-n_limitation_signature/`
