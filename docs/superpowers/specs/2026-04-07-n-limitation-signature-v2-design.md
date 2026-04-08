# N-Limitation Signature Analysis v2: Methodological Rebuild

**Date:** 2026-04-07
**Status:** Active spec
**Predecessor:** `analyses/2026-04-06-1432-n_limitation_signature/` — completed Approach A, but methodology applied post-hoc, results are a black box to the researcher. This rebuild applies correct methodology from the start and aims to produce a clean, understood, package-ready method.

## Question

Can we quantify the degree of nitrogen limitation in Prochlorococcus MED4 molecularly, using a gene signature derived from independent N-stress reference studies in the KG?

Sub-questions:
- Which reference studies does the KG hold for MED4 N-limitation?
- What genes form a reproducible N-limitation signature across references?
- Does a rank-based score separate known positive and negative controls?
- What does the signature show for Weissberg 2025 axenic vs coculture conditions (RNA-seq + proteomics)?

## Goals

1. **Clean method for packaging:** A `sig_utils` utility package with well-defined interfaces, toy-tested functions, and clear separation from analysis-specific scripts. Candidate for later productization into `multiomics_explorer/analysis/`.
2. **Researcher understanding:** Every step walked through interactively. No black boxes. The researcher can explain the method and defend the results.
3. **Correct methodology:** Research-methodology skill loaded before brainstorming (done). Proper notebook discipline, source tagging, artifacts structure from step 1.

## Out of scope

- Alteromonas N-recycling check (separate analysis)
- Approach B (pathway enrichment) and Approach C (hybrid)
- Multiple scoring metrics — rank score is the sole primary metric; others added later if needed
- Productization — flag candidates in the notebook, don't move them

## Validation strategy: positive and negative controls

Every experiment discovered in step 1 gets classified:

| Role | What it tests | Expected score |
|------|--------------|----------------|
| **Positive control** | Reference N-deprivation studies (defined the signature) + any other MED4 N-stress experiments in the KG | High rank score |
| **Negative control** | Non-nitrogen stresses (light, phosphorus), alternative N sources (Tolonen cyanate, urea), early timepoints before response onset (Tolonen 0h, 3h) | Low/zero rank score |
| **Target** | Weissberg 2025 axenic + coculture (RNA-seq + proteomics) | Unknown — this is what we're measuring |

**Gate:** Positive/negative controls must separate cleanly before we score targets. If they don't, the metric or signature needs work.

## Pipeline: 7 steps

Each step follows the do→show→explore→decide cycle. No step proceeds without the researcher's decision.

### Step 1: KG discovery — what MED4 N-limitation data exists?

**Goal:** Map all available MED4 experiments relevant to nitrogen metabolism.

**Method:**
- `list_experiments` for Prochlorococcus MED4
- `list_publications` to link experiments to studies
- Classify each experiment: positive control / negative control / target / irrelevant

**Outputs:**
- `data/experiment_scoping.csv` — experiment ID, publication, conditions, timepoints, platform, treatment type, background factors, classification
- `logs/01_discover_experiments.log` — full query results, classification reasoning

**Explore:** Walk through the scoping table. Discuss classifications. Identify any experiments that could serve as additional controls.

**Decision:** Which experiments are references, which are controls, which are targets.

### Step 2: Reference DE extraction

**Goal:** Extract full DE data for all reference and control experiments.

**Method:**
- Build a reusable extraction utility in `sig_utils/extraction.py` that wraps the DE extraction logic. All scripts (steps 2, 4, and any future extraction) use this same function. The utility:
  - Calls `differential_expression_by_gene()` with `verbose=True`, `limit=None`
  - Extracts ALL genes (not just significant) — needed for background sets and for classifying genes as "present but not significant" vs "absent"
  - Captures `table_scope` from the response envelope and adds as a column to every row
  - Returns a DataFrame in a consistent schema
- One CSV per experiment

**Outputs:**
- `data/de_*.csv` — one per experiment, named descriptively
- `logs/02_extract_reference_de.log` — row counts, gene counts, timepoints, marker gene values

**Explore:** Check marker genes (glnA, cynA, rbcL, atpD, PMM0030) in each reference dataset. Are they significant? Correct direction? What rank? Flag any metadata issues (e.g., `timepoint=single`).

**Decision:** Are the reference datasets complete and usable? Any experiments to add or drop?

### Step 3: Signature building

**Goal:** Build an N-limitation gene signature from reference studies.

**Sub-steps:**
1. **Design** — define the intersection logic, concordance criteria, direction assignment, and per-gene summarization. Worked examples with concrete numbers for each operation.
2. **Implement** — build functions in `sig_utils/signature.py`
3. **Toy-test** — synthetic DE data with hand-calculated expected outputs
4. **Apply** — run on real reference DE data
5. **QC** — filter funnel, marker gene traces, edge case tracers

**sig_utils functions (signature construction):**
- `summarize_per_gene(de_df)` → one row per gene with direction, best rank
- `intersect_references(study_a, study_b)` → concordant intersection
- `classify_non_core(study_a, study_b, core)` → extended signature categories

**Outputs:**
- `data/core_signature.csv` — core signature genes with locus_tag, gene_name, direction, per-study ranks, product, gene_category
- `data/extended_signature.csv` — non-core genes with classification
- `logs/03_build_signature.log` — filter funnel, traced genes, edge cases

**Explore:** Trace 5 marker genes + edge cases through the intersection logic. Verify filter funnel counts. Check that known biology is represented.

**Decision:** Is the core signature biologically sensible? Proceed to scoring or adjust criteria?

### Step 4: Target DE extraction

**Goal:** Extract Weissberg 2025 MED4 data for scoring.

**Method:** Same as step 2, applied to Weissberg experiments.

**Outputs:**
- `data/de_weissberg_*.csv` — one per experiment/platform/condition
- `logs/04_extract_target_de.log` — gene counts per platform, timepoint coverage, signature gene overlap

**Explore:** Check marker genes in Weissberg data. What does glnA look like in axenic vs coculture? How many signature genes are detectable per platform?

**Decision:** Any coverage issues that affect interpretation? Proceed?

### Step 5: Metric design and verification

**Goal:** Define the rank score metric formally, implement, and verify with synthetic data.

**Sub-steps:**
1. **Define on paper** — formula, properties, edge cases, worked example with 5-6 genes showing input → computation → output
2. **Implement** — `sig_utils/scoring.py` with `rank_score()` and `permutation_test()`
3. **Toy-test** — synthetic data covering: all concordant, all discordant, mixed, genes absent from dataset, gene not significant
4. **Verify** — hand-calculated expected values for each toy case

**sig_utils functions (scoring):**
- `apply_signature(signature_df, de_df)` → applied subset DataFrame with concordance info, per-gene scores. The script writes this to CSV — the inspectable intermediate showing exactly which genes matched and how.
- `rank_score(applied_df)` → single score per condition/timepoint
- `permutation_test(signature_df, de_df, n_perms)` → observed score, p-value, null distribution
- `score_with_significance(signature_df, de_df, n_perms)` → wraps apply + rank_score + permutation_test into one call, returns score and p-value per condition/timepoint. This is the main entry point for scoring.

**Outputs:**
- `sig_utils/scoring.py` — tested utility module
- `logs/05_metric_verification.log` — toy data, expected vs actual, edge case results

**Explore:** Walk through each toy example. Verify edge cases. Discuss: does this metric have the properties we want?

**Decision:** Is the metric sound? Proceed to real data?

### Step 6: Score all experiments

**Goal:** Score references (positive controls), control experiments (negative controls), and Weissberg (targets).

**Method:**
- `apply_signature()` for each experiment → applied subset CSVs
- `score_with_significance()` per condition/timepoint — computes rank score + permutation p-value in one call
- Check control separation

**Outputs:**
- `data/applied_*.csv` — one per experiment, showing signature genes in that experiment's context
- `results/scores_all.csv` — all experiments scored with rank score, p-value, and role (positive/negative/target)
- Note: reference studies scoring their own signature will have inflated scores and low p-values by construction — this is expected and flagged in caveats, not a reason to skip the test
- `logs/06_score_experiments.log` — full score table, control separation summary, gene-level contributions for marker genes

**Explore:** Full results table. Do positive controls score high? Negative controls low? Trace marker genes through scoring for specific conditions. Identify surprises.

**Decision:** Does the metric separate controls? Are Weissberg results interpretable? Proceed to visualization?

### Step 7: Interpret and document

**Goal:** Visualize results, update analysis documents, flag package candidates.

**Outputs:**
- `results/trajectory_*.png` — score over time, axenic vs coculture, with control reference bands
- `results/control_separation.png` — positive vs negative controls
- Updated: `methods.md`, `decisions.md`, `caveats.md`, `gaps_and_friction.md`, manifests
- Notebook: package-readiness assessment for sig_utils

**Explore:** Discuss the biological interpretation. RNA-seq vs proteomics comparison. What's novel, what confirms expectations, what's unexplained.

**Optional: companion Jupyter notebook.** A single notebook (`exploration/notebooks/signature_explorer.ipynb`) that loads the produced data files and lets the researcher interactively explore: signature gene table, marker gene traces, score results, key figures. Created here if the researcher wants it — decided after seeing the results.

**Decision:** Analysis complete, or further investigation needed?

## Utility package: sig_utils

### Boundary

sig_utils has two layers:
- **Extraction** (`extraction.py`) — wraps the KG Python API to produce DataFrames in a consistent schema. This is the only module that calls the KG.
- **Methodology** (`signature.py`, `scoring.py`) — pure DataFrame-in, DataFrame-out. Never calls the KG, never knows about specific experiments or organisms. A script could swap in completely different data and these modules would work the same.
- **I/O** (`io.py`) — load/save helpers for standard formats.

### Modules

```
sig_utils/
├── __init__.py
├── extraction.py     # KG DE extraction wrapper (consistent schema, table_scope)
├── signature.py      # Signature construction: summarize, intersect, classify
├── scoring.py        # Signature application, rank score, permutation test
├── io.py             # Load/save helpers for standard formats
└── tests/
    ├── test_signature.py   # Toy data tests for signature construction
    └── test_scoring.py     # Toy data tests for scoring and permutation
```

### Flow

```
DE data (DataFrame)
    → summarize_per_gene()
    → intersect_references()
    → Signature (DataFrame)

Signature + any DE data
    → score_with_significance()          # main entry point
        → apply_signature()              # → Applied subset (DataFrame / CSV)
        → rank_score()                   # → score per condition/timepoint
        → permutation_test()             # → p-value per condition/timepoint
    → score + p-value per condition/timepoint
```

### Toy-tested

Both signature construction and scoring functions are verified with hand-calculated synthetic data before touching real KG data. Verification is a notebook step, not an afterthought.

Toy tests are saved as scripts (`sig_utils/tests/test_signature.py`, `sig_utils/tests/test_scoring.py`) so they serve as the seed for a real test suite during productization. Each test script contains synthetic data, hand-calculated expected values, and assertions.

## Scripts

```
scripts/
├── 01_discover_experiments.py
├── 02_extract_reference_de.py
├── 03_build_signature.py
├── 04_extract_target_de.py
├── 05_score_experiments.py
└── 06_plot_results.py
```

Each script:
- Uses sig_utils for all reusable logic (extraction, signature building, scoring)
- Writes outputs to `data/` or `results/`
- Writes diagnostic log to `logs/`
- Has a `--explore` flag that prints marker gene traces and QC diagnostics to stdout

## Marker genes

Starting set, confirmed via `resolve_gene` in step 1:

| Gene | Locus tag | Expected direction | Why |
|------|-----------|-------------------|-----|
| glnA | PMM0920 | up | Canonical N-limitation marker — glutamine synthetase |
| cynA | PMM0370 | up | Cyanate transporter — N scavenging |
| rbcL | PMM0550 | down | RuBisCO — carbon fixation reduced under N-stress |
| atpD | PMM1452 | down | ATP synthase — energy metabolism reduced |
| PMM0030 | PMM0030 | up | Unnamed, rank 1 in v1 — strong responder, unknown function |

Additional edge-case tracers added during exploration as needed and logged in the notebook.

## Directory structure

```
analyses/YYYY-MM-DD-HHMM-n_limitation_signature_v2/
├── exploration/
│   └── YYYY-MM-DD-notebook.md     # Research notebook (append-only)
├── data/
│   ├── DATA_MANIFEST.md
│   ├── experiment_scoping.csv
│   ├── core_signature.csv
│   ├── extended_signature.csv
│   ├── de_*.csv                   # One per experiment
│   └── applied_*.csv              # Signature applied to each experiment
├── scripts/
│   ├── 01_discover_experiments.py
│   ├── 02_extract_reference_de.py
│   ├── 03_build_signature.py
│   ├── 04_extract_target_de.py
│   ├── 05_score_experiments.py
│   └── 06_plot_results.py
├── sig_utils/
│   ├── __init__.py
│   ├── extraction.py
│   ├── signature.py
│   ├── scoring.py
│   ├── io.py
│   └── tests/
│       ├── test_signature.py
│       └── test_scoring.py
├── logs/
│   ├── 01_discover_experiments.log
│   ├── 02_extract_reference_de.log
│   └── ...
├── results/
│   ├── RESULTS_MANIFEST.md
│   ├── scores_all.csv             # Includes rank scores and permutation p-values
│   ├── trajectory_*.png
│   └── control_separation.png
├── superpowers/
│   ├── spec.md                    # Copy of the design spec
│   ├── plan.md                    # Copy of the implementation plan (when created)
│   └── brainstorm-log.md          # Q&A from the brainstorming session
├── README.md
├── methods.md
├── decisions.md
├── caveats.md
├── gaps_and_friction.md
└── references.md
```

## Notebook discipline

One notebook: `exploration/YYYY-MM-DD-notebook.md`. Append-only, chronological.

Each step gets an entry with: command, inputs, outputs, QC, exploration (marker gene traces + edge cases + questions), decision.

Reruns get new entries with "why" section — never overwrite.

The notebook is the primary record of the researcher's understanding. Methods.md is the publication-ready distillation.

## Superpowers products

The analysis directory includes a `superpowers/` folder with copies of the design spec, implementation plan, and brainstorming log. These live with the analysis so the full decision history is self-contained — not scattered across `docs/superpowers/` which serves the repo-level index.

- `superpowers/spec.md` — the design spec (this document)
- `superpowers/plan.md` — the implementation plan (created next)
- `superpowers/brainstorm-log.md` — the Q&A from the brainstorming session: questions asked, options considered, decisions made, and rationale. Captures the reasoning that shaped the spec.

## Predecessor analysis

The original analysis at `analyses/2026-04-06-1432-n_limitation_signature/` remains untouched. The notebook may reference it for comparison ("v1 found 198 core genes, we found N — difference is because...") but this analysis stands independently.

Key lessons from v1 incorporated:
- Methodology skill loaded before brainstorming (not retrofitted)
- Toy-data verification before real data (mandatory, not ad-hoc)
- Single primary metric (not three unresolved)
- Verbose DE data joined into outputs immediately (not deferred)
- Proper notebook discipline from step 1
- Subagent reviews not skipped for data-producing tasks
