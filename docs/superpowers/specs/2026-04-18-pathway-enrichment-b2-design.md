# Pathway enrichment B2 — N-limitation score for Weissberg 2025 MED4 conditions

**Date:** 2026-04-18
**Predecessor:** `analyses/2026-04-09-1713-pathway_enrichment_b1/` (B1 — pathway enrichment survey)
**Tracks:** biology + skill evaluation + MCP/Python API evaluation

---

## 1. Research question

**Biological:** How N-limited is MED4 in each of the four Weissberg 2025 conditions — axenic RNA-seq, axenic proteomics, coculture RNA-seq, coculture proteomics — measured against a reference N-limitation signature derived from MED4 experiments in the KG?

**Meta:** Does the research-methodology skill guide this analysis end-to-end without retrofitting? Do the new MCP tools (`ontology_landscape`, `pathway_enrichment`, `search_ontology`) and Python primitives (`fisher_ora`, `EnrichmentResult` accessors) cover the analytical needs of a realistic multi-ontology, multi-experiment enrichment workflow?

## 2. Goals

1. **Biology.** Produce a per-condition N-limitation score for each of the four Weissberg conditions, per chosen ontology, with stability checks and narrative interpretation.
2. **Skill evaluation.** Run the analysis strictly following the research-methodology step protocol (do → show → explore → decide, two commits per step, three hard gates). Record every place the skill guided effectively or got in the way in `gaps_and_friction.md` — as it happens, not retroactively.
3. **API evaluation.** Exercise `list_experiments` (with `search_text`), `list_filter_values`, `ontology_landscape`, `search_ontology`, `pathway_enrichment` (MCP wrapper), `EnrichmentResult.explain()` / `.overlap_genes()`, and `fisher_ora` (direct Python primitive). Record API friction in `api_coverage.md`.

## 3. Scope

- **Organism of interest:** MED4 (Prochlorococcus MED4). Non-MED4 N-limitation experiments appear in figures as context but do **not** contribute to the reference signature (see §7.1).
- **Ontologies:** 1–3 chosen in Step 1b via `ontology_landscape` + `search_ontology("nitrogen")`, including BRITE subtrees as dedicated candidates.
- **Experiments:** classified into five classes (§4). The target set is the four Weissberg 2025 conditions.
- **Signal surface:** pathway enrichment (Fisher's exact + BH) on DE gene sets per `(experiment, timepoint, direction)` cluster. Per-timepoint granularity preserved throughout.
- **Out of scope:** ortholog-level analysis, custom signature gene sets, GSEA/rank-based enrichment, cluster-membership enrichment (these are different analyses for different questions).
- **Reuse of B1:** none. The B1 `enrich_utils/` package is explicitly not reused — this analysis exercises the new MCP / Python API.

## 4. Experiment classification

Every selected experiment is assigned one of five classes in Step 1a:

| Class | Meaning | Role in scoring |
|---|---|---|
| **T** (Target) | The four Weissberg 2025 MED4 conditions | Scored; the answer |
| **R** (Reference) | Unambiguous MED4 N-limitation experiments (e.g., Tolonen N-deplete) | Derive the signature; also scored (expected: high) |
| **PC** (Positive control) | N-stress-related MED4 experiments that aren't strict N-deprivation (alternative N sources, chronic low-N) | Scored; expected direction: high |
| **NC** (Negative control) | MED4 experiments unrelated to N (light, phage, salt, dark, glucose) | Scored; expected: ≈ 0; defines noise floor |
| **CTX** (Context) | Non-MED4 N-limitation experiments from other organisms | Scored against MED4 signature; emergent evidence of signature conservation (high) or MED4-specificity (low). Figure context only. |

Non-MED4 N-limitation experiments (CTX class) cannot contribute to the MED4 signature (see §7.1 for why). Soft cap on CTX selection: 3–5 experiments across 1–3 organisms, chosen for representativeness of conditions (e.g., one per cyanobacterial strain with an N experiment). Cap enforced for Fig 1 legibility, not for biological reasons.

## 5. Phase structure

Six steps. Each step follows the research-methodology skill's do → show → explore → decide cycle; each produces two commits (Commit 1 after do; Commit 2 after decide). Three hard gates (step boundary, manifest currency, chat-capture) enforced per step-protocol.

### Step 1a — Experiment discovery + classification

**do:**
- MCP orientation (interactive — notebook-captured, chat-capture style):
  - `list_organisms()`
  - `list_experiments(summary=True)` — overall by-organism, by-treatment, by-table_scope breakdowns.
  - `list_experiments(search_text="nitrogen", verbose=True)` — cross-organism fuzzy N catch with Lucene relevance scores.
  - `list_experiments(organism="MED4", verbose=True)` — full MED4 landscape.
  - Optional: `list_filter_values("treatment_type")` — sanity-check N-tag vocabulary; catch tags like `nutrient_starvation` that wouldn't match `"nitrogen"` text.
- Researcher selects experiment IDs per class (T, R, PC, NC, CTX). No pre-filter on `table_scope`; restricted-scope experiments accepted with explicit caveat when they are the only candidate for a class.
- `scripts/01_select_experiments.py` reproduces the selected `list_experiments` calls via the Python API (`from multiomics_explorer import list_experiments` — run from repo root, not `--directory` to `multiomics_explorer`), applies the classification labels (hard-coded from 1a's decisions), and writes `data/experiments_classified.csv` via `experiments_to_dataframe` or `to_dataframe`.

**show:** counts by class × organism × omics_type; `table_scope` and `gene_count` distributions; sample rows per class.
**explore:** spot-check close-call classifications; document any restricted-`table_scope` inclusions in the notebook entry with rationale.
**decide:** classification locked → `experiments_classified.csv` committed.

### Step 1b — Ontology & level selection + key-pathway panel

**do:**
- `scripts/02_ontology_landscape.py` — programmatic wrapper for reproducibility:
  - For each organism with selected experiments: `ontology_landscape(organism=<org>, experiment_ids=<selected for that org>, verbose=True)` via Python API — for MED4 the experiment set is `T∪R∪PC∪NC`; for non-MED4 organisms the experiment set is `CTX` for that organism. BRITE rows in the output carry `tree` / `tree_code` sparse fields, so per-tree candidates surface without extra calls.
  - `search_ontology("nitrogen")` across all ontologies.
  - Both calls flattened via `to_dataframe` and written as CSV: `data/landscape_<org>.csv`, `data/nitrogen_ontology_search.csv`. Envelope summary fields (totals, `truncated`, `by_ontology` breakdown) captured in the notebook entry, not the CSV.
- Researcher selects 1–3 ontologies with explicit justification (coverage pick, relevance pick, BRITE subtree). Hard cap: 3. Writes `ontology_selection.md` (ranked table + rationale) manually — this is the interactive decision, not a computation.
- Identify key-pathway term_ids in each selected ontology for the canonical N-response categories (N-metabolism, photosynthesis, amino-acid transport, ribosome). These are the biological anchor for every downstream QC. Write `exploration/key_pathways.csv` (ontology, term_id, term_name, expected_direction, canonical_gene_marker).

**show:** ranked tables per organism; organism-coverage gaps per selected ontology; top populated terms per pick; key-pathway term_ids with their member counts.
**explore:** sanity-check key-pathway term_ids via `genes_by_ontology(term_id=..., organism="MED4", limit=None)` — does the N-metabolism term contain `glnA`, `amt`, `cynA`? Does photosynthesis contain `rbcL`, `psbA`, `psbD`? If no, the ontology or level is wrong — redo Step 1b.
**decide:** ontology set + key-pathway panel locked → both files committed.

### Step 2 — Enrichment run

**do:**
- `scripts/03_run_enrichment.py`: for each (organism, ontology) pair from Step 1b, split the organism's selected experiments by `table_scope` and run **two** `pathway_enrichment` calls:
  - **Call A — detected-genes group.** Experiments with `table_scope == "all_detected_genes"` → `background="table_scope"` (the pathway_enrichment default). `table_scope` is a principled per-cluster denominator: the set of genes that were measured and could have been significant.
  - **Call B — restricted-scope group.** Experiments with `table_scope ∈ {significant_only, significant_any_timepoint, top_n, filtered_subset}` → `background="organism"`. For these experiments the DE table equals (or is a strict subset of) the significant gene set, so `table_scope` background would collapse the Fisher 2×2 and break enrichment. Organism background is the only principled alternative — with the caveat that it inflates `N` for experiments that only measured a proper subset of the genome.
  - Skip either call if its experiment group is empty.
- Each call: `pathway_enrichment(organism, experiment_ids=<group>, ontology, level, background=<as above>, direction="both", significant_only=True)` via the Python API.
- Collect each `EnrichmentResult.results` DataFrame; concat with added `organism`, `ontology`, `background_used` (`"table_scope"` or `"organism"`) columns into `data/enrichment_all.csv`. Also pickle the `EnrichmentResult` objects (`data/enrichment_results.pkl`) so downstream steps can use `.explain()` and `.overlap_genes()` for drill-down.
- Per-timepoint granularity preserved — cluster key = `experiment_id | timepoint | direction`. NaN timepoints appear as `"NA"` (`pathway_enrichment` tool handles this correctly).

**show:** per-(org, ontology) cluster/test/significance counts, split by `background_used`; key-pathway signed_score across all clusters as a diagnostic heatmap (`exploration/qc/step2_key_pathway_heatmap.png` — not a publication figure).
**explore:** do R clusters show expected key-pathway enrichment in the expected direction? `result.explain(R_cluster, key_pathway_term_id)` for canonical R clusters (e.g., Tolonen N-deplete, N-metabolism term) — do canonical marker genes (`glnA`, `amt`, `cynA`) appear in the overlap? NC clusters ≈ zero on key pathways?
**decide:** proceed / redo → `enrichment_all.csv` + pickle + QC figure committed.

### Step 3 — Reference signature derivation (MED4-only)

**do:**
- `scripts/04_derive_signature.py`: restrict `enrichment_all.csv` to MED4 R clusters (all ontologies). For each `(ontology, term_id)`, count `n_clusters_supporting = count(padj < 0.05)` per direction. Signature rule (M1 = b):
  - **Core rule:** pathway is in the signature if it's enriched (padj < 0.05) in **≥3 MED4 R clusters, same direction**.
  - **Fallback rule:** if fewer than 5 pathways per ontology survive the core rule, relax to **≥2 MED4 R clusters, same direction** for that ontology and record the fallback in `decisions.md`.
- Write `data/reference_signature.csv` with columns: `ontology, term_id, term_name, direction, n_clusters_supporting, contributing_clusters, rule_applied` (core or fallback).

**show:** signature size per ontology; distribution of `n_clusters_supporting` (histogram); which key pathways made it into the signature; any fallback application logged.
**explore:** signature composition — sensible? Are all key pathways in the signature? If a key pathway (e.g., N-metabolism) is NOT in the signature, investigate: is the MED4 R set too small? Does the reference experiment not behave as expected? Capture in chat-capture.
**decide:** proceed / adjust threshold → `reference_signature.csv` committed.

### Step 4 — Scoring (T + R + PC + NC + CTX)

**do:**
- **Pre-registration before computing:** update `decisions.md` with expected outcomes per T condition (e.g., "axenic-RNA: high, in PC range; coculture-RNA: unknown — B1 showed transcriptional suppression, may score ≈ NC; coculture-protein: may retain signal per B1"). This section is finalized and staged *before* the scoring script runs, then included in Commit 1 alongside the script and outputs — so the pre-registration commit is a single atomic commit, not separate.
- `scripts/05_compute_scores.py`: for each (cluster, ontology) pair, compute Layer A score using formula M2 (β):

  ```
  score_A(cluster, ontology) = mean(sign_ref_p * signed_score(cluster, p)
                                    for p in signature(ontology)
                                    if p present in cluster results)
  ```

  where `sign_ref_p` = `+1` if `direction(p) == "up"` in signature, `-1` if `"down"`; `signed_score(cluster, p)` is `pathway_enrichment`'s per-result `signed_score` column (= `sign × -log10(padj)`, sign from direction). Clusters not tested against pathway p contribute nothing (excluded from the mean).

- Stability checks (M3):
  - **LOO on signature pathways:** per T cluster, recompute score leaving out each pathway in turn. Flag if any single-pathway exclusion flips the score sign or drops it >50%.
  - **LOO on R experiments:** re-derive signature excluding each R experiment in turn (re-filter `enrichment_all.csv`, not re-enrich), re-score T. Flag if removing any single R experiment flips a T condition's classification.
  - **Cross-ontology agreement:** for each T cluster, compare scores across all selected ontologies. Flag disagreement in direction or major magnitude discrepancy.
- **NC noise floor (per ontology, per direction of T score):** for each ontology `o`, let `nc_scores_o` = the set of Layer A scores over all NC clusters × ontology `o`. Compute `nc_mean_o = mean(nc_scores_o)` and `nc_std_o = std(nc_scores_o, ddof=1)`. Thresholds:
  - `score ≥ nc_mean_o + 2·nc_std_o` → **detectable signature** in ontology `o`.
  - `|score − nc_mean_o| < 2·nc_std_o` → **no signal** (indistinguishable from NC).
  - `score ≥ pc_mean_o − 2·pc_std_o` (where `pc_*_o` computed analogously over PC clusters) → **PC-like strength**.
- Stability-check flags are evaluated per T cluster:
  - **LOO pathway flip:** score changes sign (positive → negative or vice versa) OR `|new_score| < 0.5·|original_score|` when any single signature pathway is removed.
  - **LOO R experiment flip:** T cluster's threshold classification (detectable / no-signal / PC-like) changes when any single R experiment is excluded from signature derivation.
  - **Cross-ontology disagreement:** T cluster's threshold classification differs across ontologies (e.g., "PC-like" in CyanoRak but "no signal" in KEGG).
- Write `results/scores_all.csv` (cluster × ontology score) and `results/score_summary.csv` (T condition × ontology: peak score, peak timepoint, noise-floor comparison, classification).

**show:** score summary table (T / R / PC / NC per ontology); NC calibration check; LOO results per T cluster; cross-ontology agreement.
**explore:** per-T contribution decomposition — which signature pathways contributed most (top-5 by absolute magnitude)? Which disagreed (opposite sign from reference)? For top contributors, `result.explain(T_cluster, contributing_pathway)` — are the driving genes canonical markers? Leave-one-out results — any fragile scores?
**decide:** scores final / redo → score CSVs and LOO outputs committed.

### Step 5 — Figures + write-up

**do:**
- `scripts/06_make_figures.py`:
  - **Fig 1** — unified signed-enrichment heatmap. Rows = visible pathways (signature + key-pathway panel + any pathway significantly enriched in ≥1 T cluster; capped to ~60–80 rows by aggregate magnitude; multi-panel by ontology if overflow). Columns = clusters ordered by class (T → R → PC → NC → CTX), then by experiment × timepoint. Annotations: ontology (row grouping), class (column color), organism (column color — MED4 bolded/highlighted), direction. Signature rows visually highlighted (bold labels + left annotation column + subtle background shade).
  - **Fig 2** — score-trajectory lineplots, one panel per experiment. x = timepoint (ordinal or hours); y = Layer A score; one line per ontology (color). Single-timepoint experiments = single point. Panel title = experiment_id + class + organism. MED4 T panels visually framed. NC panels included as baseline. Horizontal reference lines at `nc_mean ± 2σ` per ontology.
- Write `methods.md`, `caveats.md`, `gaps_and_friction.md`, `api_coverage.md`, `README.md`, update `DATA_MANIFEST.md` and `RESULTS_MANIFEST.md` with the final file set.

**show:** final figures inline; summary table in README with per-condition Layer A scores.
**explore:** legibility (MED4 distinguishable? trajectories readable?); story clarity (does Fig 1 + Fig 2 + summary table answer the research question?); anything missing or misleading.
**decide:** publish → Step 5 commit; close the analysis with a final commit of all docs.

## 6. Scaffolding commit (before Step 1a)

Per step-protocol: before Step 1 begins, create:
- `analyses/2026-04-18-HHMM-pathway_enrichment_b2/` directory with subdirectories per artifacts-guide §Directory structure.
- Empty `data/DATA_MANIFEST.md`, `results/RESULTS_MANIFEST.md` with header.
- `exploration/2026-04-18-notebook.md` stub with spec-walkthrough section.
- Per-analysis `.gitignore` with explicit-file entries only (no blanket `data/*` or `results/*`):
  ```
  # No large intermediates expected — enrichment results are small
  __pycache__/
  *.pyc
  ```
- Copy `spec.md` (this document) and `brainstorm-log.md` into `superpowers/`.

## 7. Methodology details

### 7.1 Why MED4-only references

The signature describes MED4 biology — pathways that respond to N-limitation in MED4 specifically. Weissberg is MED4; scoring Weissberg against a MED4-derived signature is an apples-to-apples test.

Non-MED4 N-limitation experiments carry two risks as signature contributors: (i) organism-specific pathway differences (e.g., CyanoRak category populations differ across strains), (ii) different experimental designs (growth conditions, statistical test, table_scope). Rather than dilute the signature, non-MED4 experiments appear as context in Fig 1 and get scored against the MED4 signature — their scores become emergent evidence of signature conservation (high) or MED4-specificity (low).

### 7.2 Worked example — signature derivation (M1 = b)

Suppose MED4 R set = Tolonen N-deplete with 2 replicates × 3 timepoints × 2 directions = 12 R clusters, enriched against CyanoRak level 1 (110 terms).

Consider pathway `cyanorak.role:E.4` (N-metabolism):
- padj < 0.05 in 8/12 R clusters, all with `direction="up"` (same direction).
- `n_clusters_supporting = 8 ≥ 3` ✓ — signature member, direction = "up".

Consider pathway `cyanorak.role:D.2` (catch-all adaptation):
- padj < 0.05 in 5/12 R clusters, split 3 up / 2 down.
- No consistent direction → not in signature.

Consider pathway `cyanorak.role:M.3` (some niche category):
- padj < 0.05 in 2/12 R clusters, both "down".
- `n_clusters_supporting = 2 < 3` → not in signature under core rule.
- If final signature for CyanoRak has <5 pathways → fallback to `≥2` → `M.3` admitted with `rule_applied="fallback"`.

### 7.3 Worked example — scoring (M2 = β)

Suppose signature for CyanoRak has 12 pathways: 8 "up", 4 "down".

For Weissberg coculture-RNA timepoint t3 cluster:
- Pathway `E.4` (ref direction "up"): Weissberg `signed_score = +4.2` → contribution = `(+1) * (+4.2) = +4.2` (agrees with reference).
- Pathway `J.2` (photosynthesis, ref direction "down"): Weissberg `signed_score = -3.1` → contribution = `(-1) * (-3.1) = +3.1` (agrees).
- Pathway `Q.1` (AA transport, ref direction "up"): Weissberg `signed_score = -0.8` → contribution = `(+1) * (-0.8) = -0.8` (disagrees).
- Suppose 9 other pathways average contribution `+1.5`:

```
score_A = mean([+4.2, +3.1, -0.8, +1.5, +1.5, ...]) for 12 pathways
        ≈ +1.8
```

Compare to:
- NC noise floor (CyanoRak): `0.1 ± 0.4` → T score `+1.8` is well above noise.
- PC mean: `+2.5`. T score `+1.8` is in "intermediate" range.
- R mean: `+3.3`. Not at full R magnitude.

**Interpretation:** coculture-RNA t3 shows a detectable but attenuated N-limitation signature. Report magnitude, flag attenuation relative to R.

### 7.4 Source tagging

Per skill Rule 3:
- `[KG]` — all numbers in notebook, methods.md, caveats.md, figures.
- `[interpretation]` — biological reasoning using intrinsic knowledge (e.g., "photosynthesis downregulation is canonical N-response behavior").
- `[gap]` — anything the KG can't answer (e.g., "no MED4 protein-level N-depletion time course outside Weissberg in this KG").

Tagging applies to notebook entries and analysis documents. Chat reasoning during exploration does not need tags (per B1's skill friction observation).

### 7.5 Locus tags, not gene names

Per skill Rule 2:
- Signature CSV columns: `term_id, term_name, ...` — term_ids are the primary key.
- When drilling into pathway composition (`result.explain`, `result.overlap_genes`), locus_tags are primary; gene_names (e.g., `glnA`) appear as labels.
- Figure row labels may use `term_name` for readability but `term_id` must be in the CSV.

## 8. Risks and contingencies

1. **Insufficient MED4 R experiments.** If Phase 1a finds only 1 MED4 N-limitation experiment, the ≥3-cluster signature threshold is barely meaningful (single experiment × a few timepoints). Contingency: relax to fallback `≥2`, document the weakness in `caveats.md`. If fewer than 3 MED4 R clusters exist, the analysis cannot proceed as designed — flag to researcher, reconsider scope.

2. **Ontology-organism mismatch.** CyanoRak is cyanobacteria-specific. If a non-MED4 context experiment uses a non-cyanobacterial organism (e.g., Alteromonas), that organism is silently absent from the CyanoRak enrichment output. Contingency: document per-ontology organism coverage in `ontology_selection.md`; if non-MED4 context becomes ontology-incompatible, show it only in compatible ontologies.

3. **Mixed-background clusters across experiments.** Restricted-`table_scope` experiments use `organism` background (Step 2 Call B); `all_detected_genes` experiments use `table_scope` background (Step 2 Call A). `fold_enrichment` and `bg_ratio` magnitudes are not directly comparable across the two background types — the organism background has a larger `N`, so `fold_enrichment` runs larger for the same pathway. Contingency: `signed_score` (= `sign × −log10(padj)`) is still comparable because BH is applied per cluster. Report `background_used` in `scores_all.csv`; flag in `caveats.md` if any T/R clusters use the organism background so interpretation accounts for it. Weissberg T tables with `significant_only` or similar are specifically handled via Call B, not silently dropped.

4. **Signature dominated by catch-all pathways.** Per B1 caveat C3, broad categories (CyanoRak D.2, R.2) routinely show high enrichment due to size. Contingency: flag catch-all signature members in `caveats.md`; consider computing Layer A with and without them as a stability check alongside LOO.

5. **No cross-ontology agreement.** If the 2–3 ontologies give discordant Layer A scores for a T condition, the claim must be qualified per ontology. Contingency: no single scalar answer; report the per-ontology vector with interpretation.

6. **LOO on R experiments collapses the signature.** If removing a single R experiment drops signature size below 5 per ontology, the signature is really just that one experiment. Contingency: document, treat the analysis as descriptive rather than reference-anchored.

## 9. Artifact plan

Analysis directory: `analyses/2026-04-18-HHMM-pathway_enrichment_b2/` with standard structure per artifacts-guide §Directory structure.

Files expected at completion:
- `data/experiments_classified.csv`, `data/enrichment_all.csv`, `data/enrichment_results.pkl`, `data/reference_signature.csv`, `data/DATA_MANIFEST.md`.
- `results/scores_all.csv`, `results/score_summary.csv`, `results/loo_signature.csv`, `results/loo_r_experiments.csv`, `results/fig1_heatmap.png`, `results/fig1_heatmap.pdf`, `results/fig2_trajectories.png`, `results/fig2_trajectories.pdf`, `results/RESULTS_MANIFEST.md`.
- `exploration/2026-04-18-notebook.md`, `exploration/key_pathways.csv`, `exploration/qc/` (diagnostic figures per step).
- `scripts/01_select_experiments.py`, `02_ontology_landscape.py`, `03_run_enrichment.py`, `04_derive_signature.py`, `05_compute_scores.py`, `06_make_figures.py`, plus any `explore_*.py` for ad-hoc investigation.
- Root: `README.md`, `methods.md`, `decisions.md`, `caveats.md`, `gaps_and_friction.md`, `api_coverage.md`, `references.md` (bibliographic — original data publications, methods references, tool/KG citations per artifacts-guide).
- `superpowers/spec.md`, `superpowers/plan.md`, `superpowers/brainstorm-log.md`.

No `{name}_utils/` package — this analysis exercises the new shared API rather than building reusable utilities locally.

**Artifact format convention:**
- **CSV** (via `to_dataframe`) for all API outputs — landscape rows, enrichment results, signature, scores, experiment metadata. Tabular, greppable, token-efficient for Claude Code inspection, one-liner to load via `pd.read_csv`.
- **Pickle** only for `EnrichmentResult` objects (`data/enrichment_results.pkl`) — downstream steps need the live `.explain()` and `.overlap_genes()` accessors, which CSV cannot preserve.
- **Notebook markdown** for envelope metadata (totals, `truncated`, per-ontology breakdowns) — captured at the `show` phase of each step, not kept as a separate file.
- **PNG + PDF/SVG** for figures per artifacts-guide.

## 10. Meta-deliverables (goals 2 + 3)

- **`gaps_and_friction.md`** — live document, updated at every step's decide phase. Sections: KG data bugs, KG gaps, MCP friction, Skill/methodology friction, Process retrospective. Format per artifacts-guide.
- **`api_coverage.md`** — table of MCP tools and Python primitives used, for what, where each worked or failed. Compared against B1's friction list to show resolved / new / persistent friction.

Format for `api_coverage.md`:

| Tool / Function | Step used | Purpose | Worked well | Friction / gaps |
|---|---|---|---|---|
| `list_experiments(search_text="nitrogen")` | 1a | Cross-organism N discovery | ... | ... |
| `ontology_landscape` | 1b | Ontology/level ranking | ... | ... |
| ... | | | | |

## 11. References

- B1 analysis: `analyses/2026-04-09-1713-pathway_enrichment_b1/` — prior pathway enrichment survey; reference for biological markers, ontology-selection methodology lessons, catch-all caveats.
- Research methodology skill: `.claude/skills/research-methodology/SKILL.md` and references — rules for KG usage, artifacts, notebook, step protocol.
- Pathway enrichment methodology: `docs://analysis/enrichment` (MCP resource) — Fisher ORA, BH per cluster, signed_score, clusterProfiler mapping.
- Pathway enrichment tool: `docs://tools/pathway_enrichment` — parameters, response format, gotchas.
- Ontology landscape tool: `docs://tools/ontology_landscape` — coverage × size_factor ranking, genome_coverage caveat.
- Weissberg 2025 bioRxiv (DOI `10.1101/2025.11.24.690089`) — target publication.
