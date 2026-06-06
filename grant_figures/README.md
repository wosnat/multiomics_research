# Grant-proposal figures

Vector figures describing the project (knowledge graph + MCP server + AI research
agent) for a funding proposal. Regenerate with:

```bash
uv run grant_figures/make_figures.py
```

Each figure is written as **PDF** (drop into LaTeX/Word), **SVG** (edit in
Inkscape/Illustrator — text stays as text), and **PNG** (quick preview). Text is
embedded/editable in both PDF and SVG.

## The figures

Figures are written for a **biology / bioinformatics audience**: they lead with
the science, keep the software framing light, and show headline numbers rather
than full cardinalities.

| File | Title | Use it for |
|---|---|---|
| `fig1_architecture` | The system in 3 layers — integrated knowledge graph (with its **data model** shown) → ask questions (plain language or Python) → AI research assistant, researcher in the loop | The "what we built" / aims overview figure |
| `fig2_workflow` | From question to finished result in six steps, with a researcher-approved checkpoint at each step | The "how it produces rigorous, reproducible science" figure |
| `fig2a_workflow_example` | The same six steps **worked through one real study**, with a small real QC/artifact thumbnail on each relevant step (3 select experiments · 4 method QC · 5 results · 6 evaluation scorecard) | The "here's the process producing real artifacts" companion to Fig 2 |
| `fig3_exemplar` | Science win, single hero plot: MED4 dies axenically under N-starvation but survives in coculture; KG re-analysis reproduces the rescue and surfaces RNA/protein discordance | The compact "why it matters" figure |
| `fig3a_exemplar_panels` | The same science win as a **four-panel data figure** (A trajectory · B RNA-vs-protein scatter · C controls · D signature composition) — all real embedded plots | The detailed "why it matters / feasibility" figure |

**Embedded real plots (not mock-ups).** Every data chart is rendered live from
the committed N-limitation analysis
(`analyses/2026-04-08-1038-n_limitation_signature_v2/`), so the figures stay in
sync if the analysis is re-run:

- **Fig 3** & **Fig 3a panel A** — the **90-day signature-score trajectory**:
  axenic RNA-seq (single point inside the reference N-limited band), coculture
  RNA-seq (collapses to ≈0 — the rescue), coculture proteomics (peaks ≈0.22 at
  day 31), axenic proteomics (context). `draw_trajectory()`, from `scores_all.csv`.
- **Fig 3a panel B** — **RNA vs protein log₂FC** scatter for the 147 signature
  genes at coculture day 31; the highlighted quadrant (protein ↑, RNA flat/down)
  is the discordance. `draw_discordance()`, from the `applied_de_*` tables.
- **Fig 3a panel C** — **specificity check**: reference N-dep studies score high,
  negative controls score low (the high-light cross-responsive control is flagged,
  not hidden). `draw_controls()`, from `scores_all.csv` roles.
- **Fig 3a panel D** — **signature composition**: 189 genes (74↑ / 115↓) by
  functional category. `draw_composition()`, from `core_signature.csv`.
- **Fig 2a** — small per-step thumbnails: `mini=True` variants of the renderers
  above (steps 4–5), plus `draw_experiment_matrix()` (step 3 selection) and
  `draw_scorecard()` (step 6 — prediction-vs-outcome assessment + caveats), so
  the process is shown producing real artifacts step by step.

To retarget the charts, point `SCORES_CSV` / `SIG_CSV` / `RNA_COC_CSV` /
`PRO_COC_CSV` at different result files.

The tool layer (39 typed tools in 9 software families, exposed over MCP + a
Python API) is collapsed in Fig 1 into the **five question types** a biologist
actually asks: expression & DE, function & pathways, orthologs / cross-organism,
metabolism & chemistry, and sequence & genome context.

**KG data model in Fig 1.** The Gene-centered hub inside the Knowledge Graph
layer (Gene → Publication, Experiment, Organism, Ortholog group, Ontology,
Metabolite) is a native, on-brand recreation adapted from
`multiomics_biocypher_kg/docs/kg_concept_figure.svg` — drawn directly in
matplotlib (`draw_kg_model()`), so it matches the figure's style and stays
vector. No SVG-rasterizer dependency required.

## Provenance of the numbers (so you can defend them)

**Fig 1 — KG cardinalities** are the live release counts from
`kg_release_info` (KG build `2026-06-03`):

- 120,416 genes · 232,758 expression edges · 197 experiments · 43 publications
  · 45 organisms
- 46,438 ortholog groups · 3,230 metabolites — from `docs://guide/concepts`
  (current-snapshot table)
- 39 MCP tools in 9 families · 14 ontologies · 4 evidence layers — from
  `docs://guide/start_here` and `docs://guide/concepts`

> Note: `kg_release_info` reports the full multi-organism counts (120,416 genes,
> 45 organisms); the concepts doc quotes rounded "backbone" figures (~99,871
> genes, 37 organisms) for the curated core. The figure uses the live
> `kg_release_info` numbers. Pick whichever framing your proposal needs, but be
> consistent.

**Fig 3 — science numbers** come from committed analyses:

- Axenic-vs-coculture survival, 90-day time course, PRO99-lowN —
  `analyses/2026-04-27-1117-prochlorococcus_stress_axenic_vs_coculture/paper.md`
  (Weissberg et al. 2025, DOI 10.1101/2025.11.24.690089)
- 189-gene N-limitation signature (74 up / 115 down); axenic RNA-seq day-14
  rank_score 0.583, hit-rate 0.85, 158/185 concordant; coculture RNA-seq signal
  absent; coculture proteomics significant at every timepoint (peak day 31);
  genuine RNA/protein discordance —
  `analyses/2026-04-08-1038-n_limitation_signature_v2/README.md`

If you re-run after a KG rebuild, refresh the Fig 1 stats from `kg_release_info`
and the Fig 3 numbers from the analysis READMEs above.
