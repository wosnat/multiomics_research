# Pathway enrichment B2 — Implementation plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.
>
> **This is a research analysis, not a software project.** TDD patterns are adapted — tests here mean "script produces the expected artifact with the expected shape and the biological anchors look right." Several tasks are inherently interactive (show / explore / decide per the research-methodology skill's step-protocol). The plan materializes each do / show / explore / decide phase as a discrete task — per spec §5.0, the plan cannot collapse them.

**Goal:** Produce a per-condition N-limitation score for each of the four Weissberg 2025 MED4 conditions (axenic-RNA, axenic-protein, coculture-RNA, coculture-protein), anchored in a MED4 reference signature derived from KG N-limitation experiments, with stability checks + two publication figures + interpretive write-up. Simultaneously evaluate the research-methodology skill and the new MCP/Python enrichment API.

**Architecture:** Six-step analysis pipeline following research-methodology's do → show → explore → decide cycle. Each step = two commits (script + notebook). Five interleaved scripts, two pathway_enrichment MCP calls per (organism, ontology) split by `background_used`, MED4-only reference signature with temporal + bidirectional filters, sign-weighted alignment scoring with magnitude cap, three stability checks, unified multi-ontology heatmap + per-experiment trajectory lineplots.

**Tech Stack:** Python (`multiomics-explorer` Python API, pandas, matplotlib, seaborn). Neo4j-backed KG via MCP server. `uv run` from `multiomics_research` project root.

---

## Before starting — required context loading

A fresh Claude session executing this plan **must** load the following before touching any code:

- [ ] **Load research-methodology skill** via Skill tool. This plan depends on the step-protocol, research-notebook, artifacts, python-api-guide, statistical-rigor, kg-rules, gene-identity, and anti-hallucination references being active. Without it, the interactive pace and documentation rules are lost.
- [ ] **Read the full spec:** [`docs/superpowers/specs/2026-04-18-pathway-enrichment-b2-design.md`](../specs/2026-04-18-pathway-enrichment-b2-design.md). Every task below refers back to numbered sections of the spec.
- [ ] **Read the B1 gaps:** [`analyses/2026-04-09-1713-pathway_enrichment_b1/gaps_and_friction.md`](../../../analyses/2026-04-09-1713-pathway_enrichment_b1/gaps_and_friction.md). B1's lessons (temporal flips, catch-all caveats, MCP friction) inform this analysis.
- [ ] **Read the enrichment example:** `docs://examples/pathway_enrichment.py` MCP resource. This is the template for Step 2 and drill-down scripts.
- [ ] **Verify MCP server is running:** `multiomics-kg` server should be connected. Check by calling `list_organisms()` — should return Prochlorococcus MED4 and other organisms.
- [ ] **Verify Neo4j is running:** at `bolt://localhost:7687` (configured via `multiomics_explorer/.env`).
- [ ] **Verify Python environment:** `uv sync --extra analysis` from the repo root. Imports `from multiomics_explorer import list_experiments, ontology_landscape, pathway_enrichment, ...` should work.

If any of these fail, stop and report to the researcher. Do not proceed.

**Pacing reminder:** The spec §5.0 mandates interactive show + explore + decide for every step. Do not batch-run multiple scripts and present results in bulk. After each script runs, the researcher walks through outputs before the next step begins. Two commits per step (Commit 1 after do; Commit 2 after decide).

---

## Task 0 — Scaffold the analysis directory

**Files:**
- Create: `analyses/2026-04-19-HHMM-pathway_enrichment_b2/` directory tree (substitute actual HHMM from `date "+%H%M"`)
- Create: `analyses/.../exploration/notebook.md`
- Create: `analyses/.../data/DATA_MANIFEST.md`
- Create: `analyses/.../results/RESULTS_MANIFEST.md`
- Create: `analyses/.../.gitignore`
- Create: `analyses/.../superpowers/spec.md` (copied from canonical)
- Create: `analyses/.../superpowers/plan.md` (copied from canonical — this file)
- Create: `analyses/.../superpowers/brainstorm-log.md` (stub — to be filled if brainstorm transcript is captured separately)

- [ ] **Step 1: Determine HHMM timestamp, create directory tree, persist path**

```bash
cd /home/osnat/github/multiomics_research
export ANALYSIS_DIR="analyses/$(date +%Y-%m-%d)-$(date +%H%M)-pathway_enrichment_b2"
mkdir -p "$ANALYSIS_DIR"/{exploration/qc,data,scripts,logs,results,superpowers}
# Persist the path so subsequent tasks (possibly in fresh shells) can recover it
# without guessing. See "ANALYSIS_DIR preamble" below every later bash block.
echo "$ANALYSIS_DIR" > .analysis_dir
echo "Created $ANALYSIS_DIR (path saved to .analysis_dir)"
```

**ANALYSIS_DIR preamble for all later tasks.** Each task's bash block should start by re-exporting the path — the env var from Task 0 won't survive a new shell:

```bash
cd /home/osnat/github/multiomics_research
export ANALYSIS_DIR="$(cat .analysis_dir)"
```

If `.analysis_dir` is missing, Task 0 didn't run — go back and scaffold. **Do not** fall back to `date +%Y-%m-%d-HHMM` — the timestamp won't match Task 0's, and a second directory with a different name will appear.

- [ ] **Step 2: Create `.gitignore` with explicit entries only**

Write `$ANALYSIS_DIR/.gitignore`:

```
# No large intermediates expected — enrichment results are small
__pycache__/
*.pyc
```

No blanket `data/*` or `results/*` ignores. Tracking decisions are explicit.

- [ ] **Step 3: Create empty manifests with headers**

Write `$ANALYSIS_DIR/data/DATA_MANIFEST.md`:

```markdown
# Data Manifest

All files produced by extraction scripts from the multiomics KG.
Run scripts from `multiomics_research` root with `uv run`.

## Experiment landscape

| File | Rows | Genes | Timepoints | Produced by | Description |
|------|------|-------|------------|-------------|-------------|

## Enrichment outputs

| File | Rows | Genes | Timepoints | Produced by | Description |
|------|------|-------|------------|-------------|-------------|

## Signature

| File | Rows | Genes | Timepoints | Produced by | Description |
|------|------|-------|------------|-------------|-------------|
```

Write `$ANALYSIS_DIR/results/RESULTS_MANIFEST.md`:

```markdown
# Results Manifest

All files produced by scoring, plotting, and analysis scripts.

## Scores

| File | Rows | Produced by | Description |
|------|------|-------------|-------------|

## Figures

| File | Produced by | Description |
|------|-------------|-------------|
```

- [ ] **Step 4: Create notebook stub with spec-walkthrough section**

Write `$ANALYSIS_DIR/exploration/notebook.md`:

```markdown
# Pathway enrichment B2 — research notebook

Analysis directory: `analyses/2026-04-19-HHMM-pathway_enrichment_b2/`
Spec: `superpowers/spec.md`
Plan: `superpowers/plan.md`

---

## 2026-04-19 HH:MM — Spec walkthrough

### Context
Redo of B1 with new MCP/Python enrichment API. Goals: (1) N-limitation score
per Weissberg 2025 MED4 condition, (2) evaluate research-methodology skill,
(3) evaluate new pathway_enrichment / ontology_landscape tools.

### Spec sections reviewed
- §1 Research question, §2 Goals, §3 Scope, §4 Classification (T/R/PC/NC/CTX)
- §5 Phase structure (6 steps with do/show/explore/decide per spec §5.0)
- §7 Methodology details (M1 signature, M2 scoring, M3 calibration)
- §8 Risks and contingencies
- §9 Artifact plan

### Open items to resolve in Step 1a
- Specific experiments to place in each class (T/R/PC/NC/CTX)
- Whether restricted-`table_scope` experiments will be accepted
- Non-MED4 CTX selection (3–5 experiments, 1–3 organisms)

### Decision
Proceed to Step 1a.

---
```

- [ ] **Step 5: Copy spec and plan into analysis superpowers/ directory**

```bash
cp docs/superpowers/specs/2026-04-18-pathway-enrichment-b2-design.md "$ANALYSIS_DIR/superpowers/spec.md"
cp docs/superpowers/plans/2026-04-18-pathway-enrichment-b2-plan.md "$ANALYSIS_DIR/superpowers/plan.md"
touch "$ANALYSIS_DIR/superpowers/brainstorm-log.md"
```

Write a one-line stub to `brainstorm-log.md`:

```markdown
# Brainstorm log — pathway enrichment B2

Brainstorming was conducted interactively on 2026-04-18. Canonical spec in `spec.md`.
Meta findings captured in `docs/superpowers/specs/2026-04-18-research-methodology-v3-improvements-from-b2.md`.
```

- [ ] **Step 6: Commit scaffolding**

```bash
git add "$ANALYSIS_DIR"
git commit -m "scaffold: pathway enrichment B2 analysis directory"
```

Expected: one commit with the full directory tree, empty manifests, notebook stub, copied spec/plan.

---

## Task 1 — Step 1a do: experiment discovery + classification script

**Spec reference:** §5 Step 1a do

**Files:**
- Create: `$ANALYSIS_DIR/scripts/01_select_experiments.py`
- Create: `$ANALYSIS_DIR/data/experiments_classified.csv`
- Create: `$ANALYSIS_DIR/logs/step1a.log`

**Pattern:** This step has an interactive MCP-orientation sub-phase BEFORE the script. The researcher explores via chat; the script locks the selected experiment IDs. Do not run the script until the researcher has classified experiments.

- [ ] **Step 1: Run MCP orientation queries and capture observations in notebook**

In chat (not the script), call:

```
list_organisms()
list_experiments(summary=True)
list_experiments(search_text="nitrogen", verbose=True)
list_experiments(organism="MED4", verbose=True)
list_filter_values("treatment_type")  # if needed to catch tags not matching "nitrogen" fuzzy
```

Capture observations in the notebook entry's "Show" section (will be written in Task 2):
- Total experiments, by_organism breakdown.
- MED4 experiment count by treatment_type / omics_type / table_scope.
- Non-MED4 organisms with N-related experiments.
- `table_scope` and `gene_count` distributions (informational, not filtering).

- [ ] **Step 2: Interactively classify with the researcher**

Researcher and assistant agree on T/R/PC/NC/CTX assignments. Document each assignment in the chat with reasoning. Close-calls and restricted-`table_scope` inclusions get explicit caveats.

- [ ] **Step 3: Write the selection script**

Write `$ANALYSIS_DIR/scripts/01_select_experiments.py`:

```python
"""Select and classify experiments for pathway_enrichment_b2 analysis.

Hard-codes the T/R/PC/NC/CTX classifications decided interactively in Step 1a.
Re-running this script reproduces experiments_classified.csv from the current KG.
"""
from __future__ import annotations

import logging
import sys
from pathlib import Path

import pandas as pd
from multiomics_explorer import list_experiments
from multiomics_explorer.analysis import experiments_to_dataframe, to_dataframe

ANALYSIS_DIR = Path(__file__).resolve().parent.parent

# FILL IN FROM STEP 1a INTERACTIVE CLASSIFICATION — one entry per selected experiment.
# class values: "T" (target), "R" (reference), "PC", "NC", "CTX".
CLASSIFICATIONS: list[dict[str, str]] = [
    # Example structure — replace with actual selections.
    # {"experiment_id": "10.1101/2025.11.24.690089_growth_state_pro99lown_nutrient_starvation_med4_rnaseq_axenic", "class": "T", "rationale": "Weissberg axenic RNA-seq"},
    # {"experiment_id": "10.1038/ismej.2017.88_nitrogen_stress_ndepleted_pro99_medium_med4_rnaseq", "class": "R", "rationale": "Tolonen N-depleted, unambiguous N-limitation"},
    # ...
]


def main() -> int:
    logging.basicConfig(
        filename=str(ANALYSIS_DIR / "logs" / "step1a.log"),
        level=logging.INFO,
        format="%(asctime)s [%(levelname)s] %(message)s",
    )
    log = logging.getLogger("step1a")
    log.info("Starting 01_select_experiments.py")

    if not CLASSIFICATIONS:
        log.error("CLASSIFICATIONS is empty — fill in from Step 1a.")
        sys.exit(1)

    classified_df = pd.DataFrame(CLASSIFICATIONS)

    # API gap workaround: list_experiments has no experiment_ids filter, so we
    # pull the full landscape and filter locally. See api_coverage.md.
    result = list_experiments(verbose=True, limit=None)
    combined = experiments_to_dataframe(result)
    log.info(f"Pulled {len(combined)} experiment × timepoint rows from list_experiments(limit=None)")

    # Filter to our selected experiment_ids, join classification.
    selected = combined[combined["experiment_id"].isin(classified_df["experiment_id"])].copy()
    selected = selected.merge(classified_df, on="experiment_id", how="left")

    # Any experiment in CLASSIFICATIONS that didn't match the KG's returned IDs is a
    # typo or a stale ID — fail loudly rather than silently writing an incomplete CSV.
    missing = sorted(set(classified_df["experiment_id"]) - set(selected["experiment_id"]))
    if missing:
        log.error(f"Missing experiments after metadata join ({len(missing)}): {missing}")
        print(
            f"ERROR: {len(missing)} classified experiment_ids not found in list_experiments: "
            f"{missing}",
            file=sys.stderr,
        )
        sys.exit(2)

    # Write output.
    out_path = ANALYSIS_DIR / "data" / "experiments_classified.csv"
    selected.to_csv(out_path, index=False)
    log.info(f"Wrote {len(selected)} rows (from {len(classified_df)} classifications) to {out_path}")

    # Summary to stdout (for chat inspection).
    print(selected.groupby(["class", "organism_name", "omics_type"]).size().to_string())
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
```

- [ ] **Step 4: Fill in CLASSIFICATIONS with the Step 1a-decided experiments**

Replace the example CLASSIFICATIONS list with the actual ones. Each entry: `experiment_id`, `class` (T/R/PC/NC/CTX), `rationale` (1-sentence).

- [ ] **Step 5: Run the script**

```bash
cd /home/osnat/github/multiomics_research
uv run "$ANALYSIS_DIR/scripts/01_select_experiments.py"
```

Expected: `data/experiments_classified.csv` created with N rows (= len(CLASSIFICATIONS)), stdout shows class × organism × omics breakdown.

- [ ] **Step 6: Update DATA_MANIFEST.md**

Add row to the "Experiment landscape" table:

```markdown
| `experiments_classified.csv` | N | - | - | `01_select_experiments.py` | T/R/PC/NC/CTX-classified experiments with list_experiments metadata |
```

- [ ] **Step 7: Commit 1 (do-phase)**

```bash
git add "$ANALYSIS_DIR/scripts/01_select_experiments.py" \
        "$ANALYSIS_DIR/data/experiments_classified.csv" \
        "$ANALYSIS_DIR/data/DATA_MANIFEST.md" \
        "$ANALYSIS_DIR/logs/step1a.log"
git commit -m "step 1a do: experiment discovery + classification"
```

Gate check: script + outputs + log + manifest all in one commit (Gate 2).

---

## Task 2 — Step 1a show / explore / decide

**Spec reference:** §5 Step 1a show/explore/decide

- [ ] **Step 1: Present show tables in chat**

In chat, show:
- Total classified: N experiments across classes X/Y/Z.
- By class × organism × omics table (from the script's stdout).
- `table_scope` distribution.
- `gene_count` distribution per omics type (inform eligibility).
- Any restricted-`table_scope` experiments included — with caveat.

- [ ] **Step 2: Append show section to notebook entry**

Open `$ANALYSIS_DIR/exploration/notebook.md` and append the following. **Note:** the outer code fence below uses four backticks so that the inner `` ```bash `` fence renders cleanly when copied into the notebook — keep that distinction when pasting.

````markdown

## YYYY-MM-DD HH:MM — Step 1a: experiment discovery + classification

### Command
```bash
uv run scripts/01_select_experiments.py
```

### Outputs
- [experiments_classified.csv](../data/experiments_classified.csv) — N rows

### QC
[Class × organism × omics breakdown table in markdown]
[table_scope distribution]
[gene_count distribution per omics]

### Exploration (agent-driven QC)
- ...
````

- [ ] **Step 3: Interactive walkthrough with researcher**

Researcher asks questions about close-call classifications, missing experiments, restricted-scope inclusions. Capture per spec §7.4 format:

```markdown
### Chat exploration

**Q: [researcher's question]**
Data: [what was looked up]
Finding: [concrete numbers]
Impact: [on downstream scope]
```

- [ ] **Step 4: Researcher decides**

Decision: proceed / redo classification / adjust scope. Append to notebook entry:

```markdown
### Decision
[Continue / redo with adjustments]
```

- [ ] **Step 5: Commit 2 (decide-phase)**

```bash
git add "$ANALYSIS_DIR/exploration/notebook.md"
git commit -m "step 1a decide: classification approved"
```

Gate check: Gate 1 (step boundary — cannot start Step 1b until this commit lands), Gate 3 (chat-capture section exists before decide).

---

## Task 3 — Step 1b do: ontology landscape + selection

**Spec reference:** §5 Step 1b do

**Files:**
- Create: `$ANALYSIS_DIR/scripts/02_ontology_landscape.py`
- Create: `$ANALYSIS_DIR/data/landscape_<org>.csv` (one per organism)
- Create: `$ANALYSIS_DIR/data/nitrogen_ontology_search.csv`
- Create: `$ANALYSIS_DIR/ontology_selection.md`
- Create: `$ANALYSIS_DIR/exploration/key_pathways.csv`
- Create: `$ANALYSIS_DIR/logs/step1b.log`

- [ ] **Step 1: Write the landscape script**

Write `$ANALYSIS_DIR/scripts/02_ontology_landscape.py`:

```python
"""Run ontology_landscape per organism + search_ontology('nitrogen').

Produces CSVs for each organism's landscape and the nitrogen-ontology search.
Selection of 1-3 ontologies is done interactively by the researcher — this
script only extracts the data.
"""
from __future__ import annotations

import logging
from pathlib import Path

import pandas as pd
from multiomics_explorer import ontology_landscape, search_ontology
from multiomics_explorer.analysis import to_dataframe

ANALYSIS_DIR = Path(__file__).resolve().parent.parent


def main() -> int:
    logging.basicConfig(
        filename=str(ANALYSIS_DIR / "logs" / "step1b.log"),
        level=logging.INFO,
        format="%(asctime)s [%(levelname)s] %(message)s",
    )
    log = logging.getLogger("step1b")

    # Load classifications.
    classified = pd.read_csv(ANALYSIS_DIR / "data" / "experiments_classified.csv")
    log.info(f"Loaded {len(classified)} classified experiments")

    # Per-organism landscape.
    for org, group in classified.groupby("organism_name"):
        exp_ids = group["experiment_id"].tolist()
        log.info(f"Organism {org}: {len(exp_ids)} experiments for landscape")
        try:
            result = ontology_landscape(organism=org, experiment_ids=exp_ids, verbose=True)
        except Exception as e:
            log.error(f"ontology_landscape failed for {org}: {e}")
            continue
        df = to_dataframe(result)
        slug = org.replace(" ", "_").replace("/", "_")
        out = ANALYSIS_DIR / "data" / f"landscape_{slug}.csv"
        df.to_csv(out, index=False)
        log.info(f"Wrote {len(df)} landscape rows to {out}")
        print(f"\n=== {org} landscape (top 10 by relevance_rank) ===")
        print(df.sort_values("relevance_rank").head(10).to_string())

    # search_ontology('nitrogen') across all ontologies.
    try:
        nitrogen_search = search_ontology("nitrogen")
    except Exception as e:
        log.error(f"search_ontology failed: {e}")
        return 1
    nitrogen_df = to_dataframe(nitrogen_search)
    out = ANALYSIS_DIR / "data" / "nitrogen_ontology_search.csv"
    nitrogen_df.to_csv(out, index=False)
    log.info(f"Wrote {len(nitrogen_df)} nitrogen-term rows to {out}")
    print(f"\n=== nitrogen search (top 20 by relevance) ===")
    print(nitrogen_df.head(20).to_string())

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
```

- [ ] **Step 2: Run the script**

```bash
uv run "$ANALYSIS_DIR/scripts/02_ontology_landscape.py"
```

Expected: `data/landscape_<org>.csv` per organism, `data/nitrogen_ontology_search.csv`, stdout shows top-ranked rows per organism.

- [ ] **Step 3: Interactively select 1–3 ontologies with researcher**

Per spec §5 Step 1b + §3: coverage pick, relevance pick (from nitrogen density), optional BRITE subtree. Hard cap: 3.

- [ ] **Step 4: Write `ontology_selection.md`**

Write `$ANALYSIS_DIR/ontology_selection.md`:

```markdown
# Ontology + level selection for pathway enrichment B2

## Ranked landscape (MED4)
[top-10 rows from landscape_MED4.csv with columns: ontology_type, level, relevance_rank, genome_coverage, median_genes_per_term, n_terms_with_genes, tree if BRITE]

## Ranked landscape (other organisms)
[similar tables for any non-MED4 organisms in CTX]

## Nitrogen term density per ontology/level
[counts of nitrogen-related terms per ontology, at the levels being considered]

## Selected ontologies (1–3)

### Pick 1 — <ontology>, level <N>
- Rationale: [coverage-based / relevance-based / BRITE subtree]
- Genome coverage (MED4): X%
- Median genes/term: Y
- Organism compatibility: [list organisms covered]
- Expected N-relevance: [number of nitrogen terms at this level]

### Pick 2 — ...
### Pick 3 — ...

## Organism-ontology compatibility matrix

| Ontology | MED4 | <CTX-org-1> | <CTX-org-2> |
|----------|------|-------------|-------------|
| ...      | ✓    | ✓           | ✗ (no annot)|
```

- [ ] **Step 5: Identify key-pathway term_ids per selected ontology**

For each selected ontology, look up canonical N-response term IDs via `search_ontology`. Write `$ANALYSIS_DIR/exploration/key_pathways.csv`:

```csv
ontology,term_id,term_name,expected_direction,canonical_gene_marker,rationale
cyanorak_role,cyanorak.role:E.4,nitrogen metabolism,up,glnA,canonical N-limitation response pathway
cyanorak_role,cyanorak.role:J.2,photosystem II,down,psbA,N-limitation reduces photosynthesis
...
```

- [ ] **Step 6: Sanity-check key-pathway term_ids**

For each key-pathway term_id, call `genes_by_ontology(term_id=..., organism="MED4", limit=None)` and verify canonical gene markers appear in the result. Example:

```python
from multiomics_explorer import genes_by_ontology
r = genes_by_ontology(term_id="cyanorak.role:E.4", organism="MED4", limit=None)
assert any("glnA" in g.get("gene_name", "") for g in r["results"]), "glnA should be in E.4"
```

If canonical markers are missing, the term_id is wrong — **stop** and redo ontology/level selection before proceeding.

- [ ] **Step 7: Update DATA_MANIFEST.md**

Add rows for `landscape_<org>.csv` and `nitrogen_ontology_search.csv`.

- [ ] **Step 8: Commit 1 (do-phase)**

```bash
# Collect the landscape_*.csv paths into an array — shell globs in git add are
# fragile (literal pattern if nothing matches, and "git add" complains).
mapfile -t LANDSCAPE_FILES < <(ls "$ANALYSIS_DIR/data/"landscape_*.csv 2>/dev/null || true)
if [ ${#LANDSCAPE_FILES[@]} -eq 0 ]; then
    echo "No landscape_*.csv files found — did 02_ontology_landscape.py run?" >&2
    exit 1
fi

git add "$ANALYSIS_DIR/scripts/02_ontology_landscape.py" \
        "${LANDSCAPE_FILES[@]}" \
        "$ANALYSIS_DIR/data/nitrogen_ontology_search.csv" \
        "$ANALYSIS_DIR/data/DATA_MANIFEST.md" \
        "$ANALYSIS_DIR/ontology_selection.md" \
        "$ANALYSIS_DIR/exploration/key_pathways.csv" \
        "$ANALYSIS_DIR/logs/step1b.log"
git commit -m "step 1b do: ontology landscape + selection + key-pathway panel"
```

---

## Task 4 — Step 1b show / explore / decide

**Spec reference:** §5 Step 1b show/explore/decide

- [ ] **Step 1: Present show tables in chat**

- Ranked landscape per organism (top 10 by relevance_rank).
- Genome coverage gaps per selected ontology × organism.
- Top populated terms per pick.
- Key-pathway term_ids + their member counts + canonical-marker check result.

- [ ] **Step 2: Append show section to notebook entry**

Append to `$ANALYSIS_DIR/exploration/notebook.md`.

- [ ] **Step 3: Interactive walkthrough with researcher**

- Are the selected ontologies complementary (coverage vs. relevance)?
- Are the key pathways sensible? Any missing?
- Organism-coverage gaps acceptable?

Capture Q → data → finding → impact to notebook.

- [ ] **Step 4: Researcher decides**

Decision: ontology set locked / redo. Append to notebook.

- [ ] **Step 5: Commit 2 (decide-phase)**

```bash
git add "$ANALYSIS_DIR/exploration/notebook.md"
git commit -m "step 1b decide: ontology set + key-pathway panel locked"
```

---

## Task 5 — Step 2 do: pathway enrichment run (split by background)

**Spec reference:** §5 Step 2 do

**Files:**
- Create: `$ANALYSIS_DIR/scripts/03_run_enrichment.py`
- Create: `$ANALYSIS_DIR/data/enrichment_all.csv`
- Create: `$ANALYSIS_DIR/data/enrichment_results.pkl`
- Create: `$ANALYSIS_DIR/exploration/qc/step2_key_pathway_heatmap.png` (diagnostic — not publication)
- Create: `$ANALYSIS_DIR/logs/step2.log`

- [ ] **Step 1: Write the enrichment script**

Write `$ANALYSIS_DIR/scripts/03_run_enrichment.py`:

```python
"""Run pathway_enrichment per (organism, ontology, background_used).

Splits experiments by table_scope (all_detected_genes -> 'table_scope' background;
restricted scope -> 'organism' background) so Fisher's 2x2 is valid in each call.
"""
from __future__ import annotations

import logging
import pickle
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
from multiomics_explorer import pathway_enrichment

ANALYSIS_DIR = Path(__file__).resolve().parent.parent

RESTRICTED_SCOPES = {"significant_only", "significant_any_timepoint", "top_n", "filtered_subset"}


def run_one(organism, ontology, level, experiment_ids, background, tree=None):
    kwargs = dict(
        organism=organism,
        experiment_ids=experiment_ids,
        ontology=ontology,
        level=level,
        background=background,
        direction="both",
        significant_only=True,
    )
    if tree is not None:
        kwargs["tree"] = tree
    return pathway_enrichment(**kwargs)


def main() -> int:
    logging.basicConfig(
        filename=str(ANALYSIS_DIR / "logs" / "step2.log"),
        level=logging.INFO,
        format="%(asctime)s [%(levelname)s] %(message)s",
    )
    log = logging.getLogger("step2")

    classified = pd.read_csv(ANALYSIS_DIR / "data" / "experiments_classified.csv")
    log.info(f"Loaded {len(classified)} classified experiments")

    # Read selected ontologies from ontology_selection.md manually,
    # OR encode here (researcher updates during Step 1b):
    ONTOLOGIES = [
        # {"ontology": "cyanorak_role", "level": 1, "tree": None},
        # {"ontology": "kegg", "level": 2, "tree": None},
        # {"ontology": "brite", "level": 1, "tree": "transporters"},
    ]
    if not ONTOLOGIES:
        log.error("ONTOLOGIES is empty — update from Step 1b selection.")
        return 1

    results: dict[tuple, object] = {}
    all_rows: list[pd.DataFrame] = []

    # Dedupe to one row per experiment before splitting — experiments_to_dataframe
    # returns one row per (experiment × timepoint), so the class / table_scope / id
    # columns repeat across timepoints.
    classified_unique = classified.drop_duplicates("experiment_id")

    # PICKLE STAGE-1 PROBE — run a single small enrichment call first and
    # verify single-object round-trip + .explain() BEFORE the main loop.
    # Rationale (spec §5 Step 2): fail fast, don't discover pickle breakage
    # after every call has completed. Stage 2 (dict round-trip) happens after
    # the full results dict is built.
    probe_org_row = classified_unique.iloc[0]
    probe_ont = ONTOLOGIES[0]
    try:
        probe = run_one(
            organism=probe_org_row["organism_name"],
            ontology=probe_ont["ontology"],
            level=probe_ont["level"],
            experiment_ids=[probe_org_row["experiment_id"]],
            background="organism",  # safest default for a single-experiment probe
            tree=probe_ont.get("tree"),
        )
    except Exception as e:
        log.error(f"Pickle probe enrichment call failed: {e} — cannot verify pickle round-trip")
        return 1

    probe_pkl = ANALYSIS_DIR / "data" / "_probe_single.pkl"
    try:
        with open(probe_pkl, "wb") as f:
            pickle.dump(probe, f)
        with open(probe_pkl, "rb") as f:
            loaded_probe = pickle.load(f)
        if not loaded_probe.results.empty:
            first_row = loaded_probe.results.iloc[0]
            _ = loaded_probe.explain(first_row["cluster"], first_row["term_id"])
            log.info("Pickle stage 1 probe (single object): OK (.explain verified)")
        else:
            log.info("Pickle stage 1 probe (single object): OK (probe returned empty results, .explain skipped)")
    except Exception as e:
        log.error(f"Pickle stage 1 probe FAILED: {e} — see spec §8 risk 7 for fallbacks")
        probe_pkl.unlink(missing_ok=True)
        return 2
    finally:
        probe_pkl.unlink(missing_ok=True)

    for org, org_group in classified_unique.groupby("organism_name"):
        # Split by table_scope.
        is_detected = org_group["table_scope"] == "all_detected_genes"
        detected_ids = org_group.loc[is_detected, "experiment_id"].tolist()
        restricted_ids = org_group.loc[~is_detected, "experiment_id"].tolist()
        log.info(f"Organism {org}: {len(detected_ids)} detected, {len(restricted_ids)} restricted (experiment-level counts)")

        for ont in ONTOLOGIES:
            for bg_label, exp_ids, bg_value in [
                ("table_scope", detected_ids, "table_scope"),
                ("organism", restricted_ids, "organism"),
            ]:
                if not exp_ids:
                    continue
                try:
                    result = run_one(
                        organism=org,
                        ontology=ont["ontology"],
                        level=ont["level"],
                        experiment_ids=exp_ids,
                        background=bg_value,
                        tree=ont.get("tree"),
                    )
                except Exception as e:
                    log.error(f"pathway_enrichment failed for ({org}, {ont['ontology']}, {bg_label}): {e}")
                    continue
                key = (org, ont["ontology"], bg_label)
                results[key] = result
                df = result.results.copy()
                df["organism"] = org
                df["ontology"] = ont["ontology"]
                df["background_used"] = bg_label
                all_rows.append(df)
                log.info(f"Ran {key}: {len(df)} rows, n_significant={result.generate_summary().get('n_significant', 0)}")

    combined = pd.concat(all_rows, ignore_index=True) if all_rows else pd.DataFrame()
    combined.to_csv(ANALYSIS_DIR / "data" / "enrichment_all.csv", index=False)
    log.info(f"Wrote {len(combined)} total rows to enrichment_all.csv")

    # PICKLE STAGE-2 — full-dict round-trip. Stage-1 probe was done before the
    # main loop. See spec §8 risk 7 (empirical status: both stages validated on
    # 2026-04-19; escalation paths are precautionary).
    pkl_path = ANALYSIS_DIR / "data" / "enrichment_results.pkl"
    with open(pkl_path, "wb") as f:
        pickle.dump(results, f)
    with open(pkl_path, "rb") as f:
        loaded_dict = pickle.load(f)
    try:
        verified = False
        for key, inst in loaded_dict.items():
            if not inst.results.empty:
                first_row = inst.results.iloc[0]
                _ = inst.explain(first_row["cluster"], first_row["term_id"])
                log.info(f"Pickle stage 2 (dict): OK (.explain verified on {key})")
                verified = True
                break
        if not verified:
            log.warning("Pickle stage 2: dict loaded but all values empty; .explain not exercised")
    except Exception as e:
        log.error(f"Pickle stage 2 FAILED: {e} — see spec §8 risk 7 for fallbacks")
        pkl_path.unlink(missing_ok=True)
        return 3

    # Diagnostic heatmap of key pathways × clusters (pre-figure).
    key_paths = pd.read_csv(ANALYSIS_DIR / "exploration" / "key_pathways.csv")
    key_term_ids = set(key_paths["term_id"])
    diag = combined[combined["term_id"].isin(key_term_ids)].copy()
    if len(diag) > 0:
        pivot = diag.pivot_table(
            index="term_name",
            columns=["organism", "experiment_id", "timepoint", "direction"],
            values="signed_score",
            aggfunc="first",
        )
        # Cap at +/-10 for visualization
        pivot = pivot.clip(lower=-10, upper=10)
        fig, ax = plt.subplots(figsize=(max(12, 0.3 * pivot.shape[1]), max(4, 0.3 * pivot.shape[0])))
        sns.heatmap(pivot, center=0, cmap="RdBu_r", vmin=-10, vmax=10, ax=ax, cbar_kws={"label": "signed_score (capped)"})
        ax.set_title("Step 2 diagnostic: key pathways × clusters")
        plt.tight_layout()
        qc_path = ANALYSIS_DIR / "exploration" / "qc" / "step2_key_pathway_heatmap.png"
        plt.savefig(qc_path, dpi=150, bbox_inches="tight")
        plt.close(fig)
        log.info(f"Wrote diagnostic heatmap to {qc_path}")

    # Summary counts.
    print("\n=== Cluster/test counts per (org, ontology, background_used) ===")
    print(combined.groupby(["organism", "ontology", "background_used"])
          .agg(n_tests=("term_id", "count"),
               n_significant=("p_adjust", lambda x: (x < 0.05).sum()),
               n_clusters=("cluster", "nunique")).to_string())

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
```

- [ ] **Step 2: Fill in ONTOLOGIES list from Step 1b selection**

- [ ] **Step 3: Run the script**

```bash
uv run "$ANALYSIS_DIR/scripts/03_run_enrichment.py"
```

Expected:
- `data/enrichment_all.csv` with all rows.
- `data/enrichment_results.pkl` after stage-2 pickle round-trip passes.
- `exploration/qc/step2_key_pathway_heatmap.png` diagnostic figure.
- Log shows `Pickle round-trip stage 1: OK` and `stage 2: OK`.
- Exit code 0.

If pickle round-trip fails (exit code 2 or 3), STOP and apply spec §8 risk 7 fallbacks in order: split-to-subdir, dill/cloudpickle, constituent-state, inline re-run.

- [ ] **Step 4: Update DATA_MANIFEST.md**

Add rows for `enrichment_all.csv` and `enrichment_results.pkl`.

- [ ] **Step 5: Commit 1 (do-phase)**

```bash
git add "$ANALYSIS_DIR/scripts/03_run_enrichment.py" \
        "$ANALYSIS_DIR/data/enrichment_all.csv" \
        "$ANALYSIS_DIR/data/enrichment_results.pkl" \
        "$ANALYSIS_DIR/data/DATA_MANIFEST.md" \
        "$ANALYSIS_DIR/exploration/qc/step2_key_pathway_heatmap.png" \
        "$ANALYSIS_DIR/logs/step2.log"
git commit -m "step 2 do: pathway enrichment split by table_scope"
```

---

## Task 6 — Step 2 show / explore / decide

**Spec reference:** §5 Step 2 show/explore/decide

- [ ] **Step 1: Present show tables + diagnostic figure in chat**

- Per-(org, ontology, background_used) cluster/test/significance counts (from script stdout).
- `exploration/qc/step2_key_pathway_heatmap.png` embedded in notebook / inspected in chat.
- Sanity check: are R clusters showing up-enrichment on N-metabolism, down-enrichment on photosynthesis?

- [ ] **Step 2: Write explore scripts as needed**

For each R cluster × key pathway combination, drill into contributing genes via `.explain()`. Example script `$ANALYSIS_DIR/scripts/explore_step2_r_cluster_markers.py`:

```python
"""Drill-down: do canonical marker genes appear in key pathway enrichments?"""
from __future__ import annotations

import pickle
import sys
from pathlib import Path

ANALYSIS_DIR = Path(__file__).resolve().parent.parent

with open(ANALYSIS_DIR / "data" / "enrichment_results.pkl", "rb") as f:
    results = pickle.load(f)

# Pick an R cluster and the N-metabolism term.
# Replace with actual R cluster + term_id from Step 1b key_pathways.
target_key = ("Prochlorococcus MED4", "cyanorak_role", "table_scope")
r_cluster = "<experiment_id>|<timepoint_label>|up"   # fill in
# Cluster naming: "{experiment_id}|{timepoint}|{direction}". Single-timepoint
# experiments (no timepoint) serialize as "{experiment_id}|NA|{direction}".
# Verify the exact label by inspecting result.results["cluster"].unique().
term_id = "cyanorak.role:E.4"

if target_key not in results:
    print(f"target_key {target_key} not in results. Available keys:", file=sys.stderr)
    for k in results:
        print(f"  {k}", file=sys.stderr)
    sys.exit(1)

result = results[target_key]

available_clusters = set(result.results["cluster"].unique()) if not result.results.empty else set()
if r_cluster not in available_clusters:
    print(f"cluster {r_cluster!r} not found in result. Available clusters for {target_key}:", file=sys.stderr)
    for c in sorted(available_clusters):
        print(f"  {c}", file=sys.stderr)
    sys.exit(1)

available_terms = set(result.results["term_id"].unique()) if not result.results.empty else set()
if term_id not in available_terms:
    print(f"term_id {term_id!r} not present in this cluster's enrichment results.", file=sys.stderr)
    print(f"(The term may exist in the ontology but not pass min_gene_set_size for this cluster.)", file=sys.stderr)
    sys.exit(1)

explanation = result.explain(r_cluster, term_id)
print(explanation._repr_markdown_())
```

Commit the explore script with its output (terminal transcript or saved markdown).

- [ ] **Step 3: Append show + exploration to notebook**

- [ ] **Step 4: Interactive walkthrough with researcher**

Q → data → finding → impact capture. Typical questions:
- R cluster N-metabolism up-enrichment: does `glnA` appear? `amt`? `cynA`?
- R cluster photosynthesis down-enrichment: `rbcL`, `psbA`?
- NC clusters: any accidental enrichment on key pathways?

- [ ] **Step 5: Researcher decides**

- [ ] **Step 6: Commit 2 (decide-phase)**

```bash
git add "$ANALYSIS_DIR/exploration/notebook.md" \
        "$ANALYSIS_DIR/scripts/"explore_step2_*.py
git commit -m "step 2 decide: enrichment QC passed"
```

---

## Task 7 — Step 3 do: reference signature derivation (MED4-only)

**Spec reference:** §5 Step 3 do

**Files:**
- Create: `$ANALYSIS_DIR/scripts/signature.py` (shared primitive — reused by Task 10 LOO-R)
- Create: `$ANALYSIS_DIR/scripts/04_derive_signature.py` (runner — I/O + logging)
- Create: `$ANALYSIS_DIR/data/reference_signature.csv`
- Create: `$ANALYSIS_DIR/data/signature_dropped.csv`
- Create: `$ANALYSIS_DIR/logs/step3.log`

- [ ] **Step 1a: Write the shared signature primitive module**

Write `$ANALYSIS_DIR/scripts/signature.py`:

```python
"""Shared signature-derivation primitive.

Pure functions — no I/O, no logging. Imported by 04_derive_signature.py (main
run) and 05_compute_scores.py (LOO-R stability re-derivation). Kept as a
single-file module, not a package, to match spec §9's "no _utils/ package"
discipline while still avoiding cross-script importlib hacks.

Promotion candidate: if this primitive proves reusable across analyses, a
future `multiomics_explorer.analysis.signature` module could export
`derive_for_ontology` + the threshold constants. Defer that decision until a
second use case appears.
"""
from __future__ import annotations

import json
from collections import defaultdict

import pandas as pd

TIMEPOINT_HOURS_CUTOFF = 3.0
CORE_SUPPORT = 3
FALLBACK_SUPPORT = 2
NOTABLE_SIGNED_SCORE = 3.0


def derive_for_ontology(
    df: pd.DataFrame, support_threshold: int
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Return (signature_df, dropped_df) for one ontology.

    Caller is responsible for upstream filtering (restricting to R clusters of
    one ontology, applying any temporal filter). Required columns on df:
    term_id, term_name, cluster, direction, experiment_id, p_adjust,
    signed_score.
    """
    sig_rows = []
    dropped_rows = []
    for (term_id, term_name), grp in df.groupby(["term_id", "term_name"]):
        sig = grp[grp["p_adjust"] < 0.05]
        n_up = (sig["direction"] == "up").sum()
        n_down = (sig["direction"] == "down").sum()
        up_clusters = sig.loc[sig["direction"] == "up", "cluster"].tolist()
        down_clusters = sig.loc[sig["direction"] == "down", "cluster"].tolist()
        max_abs = grp["signed_score"].abs().max() if len(grp) else 0.0

        per_exp = defaultdict(lambda: {"up": 0, "down": 0})
        for _, r in sig.iterrows():
            per_exp[r["experiment_id"]][r["direction"]] += 1
        per_exp_json = json.dumps(dict(per_exp))

        is_up = n_up >= support_threshold and n_down < support_threshold
        is_down = n_down >= support_threshold and n_up < support_threshold
        is_bidirectional = n_up >= support_threshold and n_down >= support_threshold

        if is_up:
            sig_rows.append({
                "term_id": term_id, "term_name": term_name, "direction": "up",
                "n_clusters_supporting": n_up, "n_up": n_up, "n_down": n_down,
                "contributing_clusters": "|".join(up_clusters),
                "per_experiment_breakdown": per_exp_json,
            })
        elif is_down:
            sig_rows.append({
                "term_id": term_id, "term_name": term_name, "direction": "down",
                "n_clusters_supporting": n_down, "n_up": n_up, "n_down": n_down,
                "contributing_clusters": "|".join(down_clusters),
                "per_experiment_breakdown": per_exp_json,
            })
        elif is_bidirectional:
            dropped_rows.append({
                "term_id": term_id, "term_name": term_name,
                "n_up": n_up, "n_down": n_down,
                "drop_reason": "bidirectional",
                "max_signed_score": max_abs,
                "contributing_clusters_up": "|".join(up_clusters),
                "contributing_clusters_down": "|".join(down_clusters),
                "per_experiment_breakdown": per_exp_json,
            })
        elif max(n_up, n_down) < support_threshold and (
            max_abs >= NOTABLE_SIGNED_SCORE or (n_up + n_down) >= 2
        ):
            dropped_rows.append({
                "term_id": term_id, "term_name": term_name,
                "n_up": n_up, "n_down": n_down,
                "drop_reason": "below_threshold_notable",
                "max_signed_score": max_abs,
                "contributing_clusters_up": "|".join(up_clusters),
                "contributing_clusters_down": "|".join(down_clusters),
                "per_experiment_breakdown": per_exp_json,
            })

    return pd.DataFrame(sig_rows), pd.DataFrame(dropped_rows)
```

- [ ] **Step 1b: Write the runner**

Write `$ANALYSIS_DIR/scripts/04_derive_signature.py`:

```python
"""Derive MED4-only reference N-limitation signature.

I/O + logging runner. Pure derivation logic lives in signature.py.
"""
from __future__ import annotations

import logging
from pathlib import Path

import pandas as pd

from signature import (
    CORE_SUPPORT,
    FALLBACK_SUPPORT,
    TIMEPOINT_HOURS_CUTOFF,
    derive_for_ontology,
)

ANALYSIS_DIR = Path(__file__).resolve().parent.parent


MED4_ORGANISM_SUBSTRING = "MED4"  # organism_name filter (substring, case-insensitive)


def _append_fallback_to_decisions(ontology: str, core_size: int, analysis_dir: Path) -> None:
    """Append a fallback-rule entry to decisions.md (spec §5 Step 3)."""
    decisions = analysis_dir / "decisions.md"
    from datetime import datetime
    stamp = datetime.now().strftime("%Y-%m-%d %H:%M")
    section = (
        f"\n## D-auto — Signature derivation fallback ({ontology})\n\n"
        f"**Applied:** {stamp}\n\n"
        f"Ontology `{ontology}` produced {core_size} signature pathways under the "
        f"core rule (≥{CORE_SUPPORT} same-direction R clusters). Because {core_size} < 5, "
        f"the fallback threshold (≥{FALLBACK_SUPPORT}) was applied. The resulting "
        f"signature is weaker — document in `caveats.md`.\n"
    )
    with decisions.open("a") as f:
        f.write(section)


def main() -> int:
    logging.basicConfig(
        filename=str(ANALYSIS_DIR / "logs" / "step3.log"),
        level=logging.INFO,
        format="%(asctime)s [%(levelname)s] %(message)s",
    )
    log = logging.getLogger("step3")

    enrichment = pd.read_csv(ANALYSIS_DIR / "data" / "enrichment_all.csv")
    classified = pd.read_csv(ANALYSIS_DIR / "data" / "experiments_classified.csv")

    # Restrict to MED4 R clusters — match MED4 by substring so we're robust to
    # display-name variants ("Prochlorococcus MED4", "MED4", case drift, etc.).
    organism_mask = classified["organism_name"].fillna("").str.contains(
        MED4_ORGANISM_SUBSTRING, case=False, regex=False
    )
    r_exp_ids = set(
        classified.loc[(classified["class"] == "R") & organism_mask, "experiment_id"]
    )
    if not r_exp_ids:
        log.error(
            f"No MED4 R experiments found (filter: class=R AND organism_name contains "
            f"{MED4_ORGANISM_SUBSTRING!r}). Check experiments_classified.csv classifications "
            f"and organism_name spelling."
        )
        return 1

    r_clusters = enrichment[enrichment["experiment_id"].isin(r_exp_ids)].copy()
    log.info(
        f"MED4 R experiments: {len(r_exp_ids)}, "
        f"clusters: {r_clusters['cluster'].nunique()}"
    )
    if r_clusters.empty:
        log.error(
            f"MED4 R experiments ({sorted(r_exp_ids)}) classified but produced zero "
            f"enrichment rows. Verify enrichment_all.csv covers these experiments."
        )
        return 1

    # Temporal filter: drop clusters with timepoint_hours < 3 (keep NaN / missing).
    if "timepoint_hours" in r_clusters.columns:
        hours_col = r_clusters["timepoint_hours"]
        mask_keep = hours_col.isna() | (hours_col >= TIMEPOINT_HOURS_CUTOFF)
        dropped_early = int((~mask_keep).sum())
        r_clusters = r_clusters[mask_keep].copy()
        log.info(
            f"Temporal filter (<{TIMEPOINT_HOURS_CUTOFF}h): "
            f"dropped {dropped_early} early-cluster rows"
        )
    else:
        log.warning(
            "timepoint_hours column missing from enrichment_all.csv; temporal filter skipped. "
            "Verify 03_run_enrichment.py passes cluster metadata through."
        )

    if r_clusters.empty:
        log.error("After temporal filter, r_clusters is empty — cannot derive signature.")
        return 1

    # Per-ontology derivation.
    all_sig = []
    all_dropped = []
    fallback_ontologies: list[tuple[str, int]] = []  # (ontology, core_size)
    for ontology, ont_df in r_clusters.groupby("ontology"):
        sig_df, drop_df = derive_for_ontology(ont_df, CORE_SUPPORT)
        rule = "core"
        if len(sig_df) < 5:
            core_size = len(sig_df)
            log.warning(
                f"Ontology {ontology}: {core_size} < 5 signature pathways "
                f"under core rule. Applying fallback (>={FALLBACK_SUPPORT})."
            )
            sig_df, drop_df = derive_for_ontology(ont_df, FALLBACK_SUPPORT)
            rule = "fallback"
            fallback_ontologies.append((ontology, core_size))
        sig_df["ontology"] = ontology
        sig_df["rule_applied"] = rule
        drop_df["ontology"] = ontology
        all_sig.append(sig_df)
        all_dropped.append(drop_df)
        log.info(
            f"Ontology {ontology}: signature size = {len(sig_df)} ({rule}), "
            f"dropped = {len(drop_df)}"
        )

    signature_df = pd.concat(all_sig, ignore_index=True) if all_sig else pd.DataFrame()
    dropped = pd.concat(all_dropped, ignore_index=True) if all_dropped else pd.DataFrame()

    if signature_df.empty:
        log.error(
            "Derived signature is empty across all ontologies — cannot proceed. "
            "Either R set is too small, or no pathway meets the support threshold."
        )
        return 1

    signature_df.to_csv(ANALYSIS_DIR / "data" / "reference_signature.csv", index=False)
    dropped.to_csv(ANALYSIS_DIR / "data" / "signature_dropped.csv", index=False)
    log.info(f"Wrote {len(signature_df)} signature, {len(dropped)} dropped")

    # Record any fallback applications to decisions.md (spec §5 Step 3).
    for ont, core_size in fallback_ontologies:
        _append_fallback_to_decisions(ont, core_size, ANALYSIS_DIR)
        log.info(f"Appended fallback-rule entry for ontology={ont} to decisions.md")

    print("\n=== Signature sizes per ontology ===")
    print(signature_df.groupby(["ontology", "rule_applied", "direction"]).size().to_string())
    print("\n=== Dropped sizes per ontology ===")
    print(dropped.groupby(["ontology", "drop_reason"]).size().to_string())

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
```

**If `decisions.md` doesn't exist yet:** the runner's `open("a")` creates it on first fallback write. If no fallback applies, the file is never touched here and gets created in Task 9 (pre-registration). Either order is safe; the append-mode open never truncates.

**Import note:** `uv run scripts/04_derive_signature.py` puts `scripts/` on `sys.path`, so `from signature import ...` resolves to `scripts/signature.py` without any packaging. No `__init__.py`; no `sys.path` munging.

- [ ] **Step 2: Run the script**

```bash
uv run "$ANALYSIS_DIR/scripts/04_derive_signature.py"
```

Expected: `data/reference_signature.csv`, `data/signature_dropped.csv`, stdout shows per-ontology sizes.

- [ ] **Step 3: Update DATA_MANIFEST.md**

- [ ] **Step 4: Commit 1 (do-phase)**

```bash
# If the runner appended fallback entries to decisions.md, include it in this commit.
DECISIONS_ARG=()
if [ -f "$ANALYSIS_DIR/decisions.md" ] && git status --short "$ANALYSIS_DIR/decisions.md" | grep -q .; then
    DECISIONS_ARG=("$ANALYSIS_DIR/decisions.md")
fi

git add "$ANALYSIS_DIR/scripts/signature.py" \
        "$ANALYSIS_DIR/scripts/04_derive_signature.py" \
        "$ANALYSIS_DIR/data/reference_signature.csv" \
        "$ANALYSIS_DIR/data/signature_dropped.csv" \
        "$ANALYSIS_DIR/data/DATA_MANIFEST.md" \
        "$ANALYSIS_DIR/logs/step3.log" \
        "${DECISIONS_ARG[@]}"
git commit -m "step 3 do: MED4 reference signature with temporal + tie-breaker filters"
```

---

## Task 8 — Step 3 show / explore / decide

**Spec reference:** §5 Step 3 show/explore/decide

- [ ] **Step 1: Present show tables in chat**

- Signature size per ontology (with rule_applied core/fallback).
- Distribution of `n_clusters_supporting` per ontology (histogram or stats).
- Per-R-experiment breakdown for every signature pathway as a markdown table.
- Top-10 dropped pathways by `max_signed_score`.
- Key pathways — which are in the signature, which are dropped, which are silently omitted.

- [ ] **Step 2: Write explore scripts if needed (e.g., dominance check)**

Example: `$ANALYSIS_DIR/scripts/explore_step3_single_exp_dominance.py` — flag any signature pathway where per_experiment_breakdown shows all support from one experiment.

- [ ] **Step 2b: Within-ontology pathway-gene-overlap audit (REQUIRED — new sub-task from Step 2 decide)**

**Why.** LOO stability (§5 Step 4 M3) detects *fragility* (single-pathway/experiment dominance). It does NOT detect *redundancy* — correlated gene overlap across signature pathways that makes Layer A scoring double-count the same biological event. Step 2 exploration found the concern is real: within-KEGG, `ko00190 Oxphos ∩ ko00195 Photosynthesis` = 9 atp genes (the atpA-I operon); if both enter the signature the atp signal is 2× counted. See `gaps_and_friction.md` "within-ontology pathway redundancy is not covered by LOO" for rationale.

**What to do.** Write `$ANALYSIS_DIR/scripts/explore_step3_redundancy_audit.py`:

```python
"""Within-ontology pathway-gene-overlap audit for signature pathways.

For each ontology, compute pairwise Jaccard(term_i.genes, term_j.genes) on
MED4 gene membership (from result.term2gene) among signature pathways.
Flag: strict subset (one pathway ⊂ another) or Jaccard > 0.5.

Outputs:
    exploration/qc/step3_signature_redundancy_{ontology}.csv
    exploration/qc/step3_signature_redundancy_{ontology}.png  (heatmap)

Note: audit runs within ontology only — spec §5 Step 4 scores per-ontology,
so cross-ontology overlap (cyanorak vs kegg) doesn't inflate a single
score_A value. Only within-ontology overlap matters for scoring.
"""
from __future__ import annotations

import pickle
from itertools import combinations
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

ANALYSIS_DIR = Path(__file__).resolve().parent.parent
JACCARD_FLAG_THRESHOLD = 0.5


def main() -> None:
    signature = pd.read_csv(ANALYSIS_DIR / "data" / "reference_signature.csv")
    with open(ANALYSIS_DIR / "data" / "enrichment_results.pkl", "rb") as f:
        results = pickle.load(f)

    for ontology, sig_ont in signature.groupby("ontology"):
        # Find any EnrichmentResult for this ontology; term2gene is consistent
        # across background_used within the same ontology.
        ont_result = next(
            (r for (org, o, bg), r in results.items() if o == ontology),
            None,
        )
        if ont_result is None or len(sig_ont) < 2:
            continue
        t2g = ont_result.term2gene  # long DataFrame: locus_tag × term_id
        members = {
            tid: set(t2g.loc[t2g["term_id"] == tid, "locus_tag"])
            for tid in sig_ont["term_id"]
        }
        term_ids = list(members.keys())
        n = len(term_ids)

        # Pairwise Jaccard + subset detection.
        rows = []
        for a, b in combinations(term_ids, 2):
            A, B = members[a], members[b]
            if not A or not B:
                continue
            inter = A & B
            union = A | B
            jac = len(inter) / len(union) if union else 0.0
            subset = ("a_in_b" if A <= B else
                      "b_in_a" if B <= A else "none")
            flag = jac > JACCARD_FLAG_THRESHOLD or subset != "none"
            rows.append({
                "term_id_a": a, "n_a": len(A),
                "term_id_b": b, "n_b": len(B),
                "intersect": len(inter), "union": len(union),
                "jaccard": jac, "subset_relation": subset,
                "flag": flag,
            })
        pairs = pd.DataFrame(rows).sort_values("jaccard", ascending=False)
        pairs.to_csv(
            ANALYSIS_DIR / "exploration" / "qc"
            / f"step3_signature_redundancy_{ontology}.csv",
            index=False,
        )

        # Jaccard matrix heatmap.
        J = np.zeros((n, n))
        for i, a in enumerate(term_ids):
            for j, b in enumerate(term_ids):
                if i == j:
                    J[i, j] = 1.0
                    continue
                A, B = members[a], members[b]
                if not A or not B:
                    continue
                J[i, j] = len(A & B) / max(len(A | B), 1)
        fig, ax = plt.subplots(figsize=(max(6, 0.4 * n + 2),
                                         max(5, 0.4 * n + 1.5)))
        im = ax.imshow(J, cmap="OrRd", vmin=0, vmax=1, aspect="auto")
        labels = [t.split(":", 1)[-1] for t in term_ids]
        ax.set_xticks(range(n)); ax.set_xticklabels(labels, rotation=90, fontsize=8)
        ax.set_yticks(range(n)); ax.set_yticklabels(labels, fontsize=8)
        ax.axhline(y=-0.5, color="black", lw=0.5)
        for i in range(n):
            for j in range(n):
                if i != j and J[i, j] > 0:
                    ax.text(j, i, f"{J[i, j]:.2f}", ha="center", va="center",
                            fontsize=6, color="black" if J[i, j] < 0.7 else "white")
        fig.colorbar(im, ax=ax, label="Jaccard")
        ax.set_title(f"Step 3 signature pathway gene-overlap — {ontology}")
        fig.tight_layout()
        fig.savefig(
            ANALYSIS_DIR / "exploration" / "qc"
            / f"step3_signature_redundancy_{ontology}.png",
            dpi=140, bbox_inches="tight",
        )
        plt.close(fig)

        flagged = pairs[pairs["flag"]]
        print(f"\n=== {ontology}: {len(flagged)} flagged pair(s) ===")
        if not flagged.empty:
            print(flagged[[
                "term_id_a", "n_a", "term_id_b", "n_b",
                "intersect", "jaccard", "subset_relation",
            ]].to_string(index=False))


if __name__ == "__main__":
    main()
```

Run it, write the output to notebook exploration section as a markdown table, and have the researcher decide ONE of three options per ontology — document the decision in `decisions.md`:

| option | what | when appropriate |
|---|---|---|
| **(A) Audit-only** | Keep all signature pathways; report overlaps as an M4 stability observation in Step 4; add a `caveats.md` entry noting the overlap list. | No strict subset, Jaccard < 0.7 max. Default if nothing dramatic surfaces. |
| **(B) Post-filter in Step 3** | Within each overlap cluster (pathways with Jaccard > threshold or subset relation), keep only the most-enriched signature member; drop the others into `signature_dropped.csv` with `drop_reason="redundant_with_<term_id>"`. | Strict subset or Jaccard > 0.8; redundancy inflates score materially. |
| **(C) De-weight during scoring** | Carry all pathways into the signature but scale each pathway's Layer A contribution by `1 − max_jaccard_with_other_signature_members` so redundant contributions are dampened. Step 4 scoring rule adjustment. | Middle ground — redundancy is real but pathways carry somewhat different biology. More spec-level change; defer unless audit surfaces it. |

Choosing A is the default. If the researcher picks B, implement the filter in `04_derive_signature.py` in a redo commit (per step-protocol redo path — new commit, not amend). If C, the adjustment happens in `05_compute_scores.py` at Task 10; flag it in decisions.md now so Task 10 implements.

- [ ] **Step 3: Append show + exploration to notebook**

- [ ] **Step 4: Interactive walkthrough with researcher**

- Are all key pathways in the signature?
- Any fragile (single-experiment-dominated) signature members?
- Interesting biology in dropped pathways?

- [ ] **Step 5: Researcher decides (adjust threshold / proceed)**

- [ ] **Step 6: Commit 2 (decide-phase)**

```bash
git add "$ANALYSIS_DIR/exploration/notebook.md" \
        "$ANALYSIS_DIR/scripts/"explore_step3_*.py \
        "$ANALYSIS_DIR/exploration/qc/"step3_signature_redundancy_*.csv \
        "$ANALYSIS_DIR/exploration/qc/"step3_signature_redundancy_*.png \
        2>/dev/null || true
# If option B (post-filter) was chosen, also include the redo'd signature files:
# git add "$ANALYSIS_DIR/data/reference_signature.csv" "$ANALYSIS_DIR/data/signature_dropped.csv"
# If option C (de-weight) was chosen, decisions.md must already carry the note
# for Task 10 — include it here too.
git commit -m "step 3 decide: signature locked"
```

---

## Task 9 — Step 4 pre-registration (do-phase part 1)

**Spec reference:** §5 Step 4 do (pre-registration)

- [ ] **Step 1: Create/update decisions.md with pre-registration section**

Write or append to `$ANALYSIS_DIR/decisions.md`:

**Note on numbering:** `decisions.md` was seeded at Step 2 decide with D1–D4 (temporal filter, heatmap design, AA-anchor falsification, NC calibration exclusion — see commit after Step 2 decide). This Task 9 entry appends as **D5**, not D1. Do not overwrite D1–D4; the file is append-only.

```markdown
## D5 — Pre-registration of expected T condition outcomes

**Written:** 2026-04-19 HH:MM (before Step 4 scoring)

For each Weissberg T cluster, record the expected outcome BEFORE computing
scores. This prevents post-hoc storytelling. Matched-background calibration
is used (scores compared only within their (ontology, background_used) group).

### Weissberg axenic RNA-seq
- Table scope: [all_detected_genes / significant_only / ...]
- Background used: [table_scope / organism]
- Expected: [...]
- Rationale: [...]

### Weissberg axenic proteomics
...

### Weissberg coculture RNA-seq
...

### Weissberg coculture proteomics
...
```

Fill in each condition's expectations based on the signature from Step 3, B1 findings, and biological reasoning. Commit BEFORE the scoring script runs.

- [ ] **Step 2: Commit pre-registration (will be bundled with Task 10's Commit 1)**

Do NOT commit yet — this will be combined with the scoring script commit in Task 10. Save `decisions.md` with the pre-registration content but leave uncommitted until Task 10 Step 5.

---

## Task 10 — Step 4 do part 2: scoring + stability

**Spec reference:** §5 Step 4 do (scoring)

**Files:**
- Create: `$ANALYSIS_DIR/scripts/05_compute_scores.py`
- Create: `$ANALYSIS_DIR/results/scores_all.csv`
- Create: `$ANALYSIS_DIR/results/score_summary.csv`
- Create: `$ANALYSIS_DIR/results/loo_signature.csv`
- Create: `$ANALYSIS_DIR/results/loo_r_experiments.csv`
- Create: `$ANALYSIS_DIR/logs/step4.log`

- [ ] **Step 1: Write the scoring script**

Write `$ANALYSIS_DIR/scripts/05_compute_scores.py`:

```python
"""Compute Layer A scores with sign-weighted alignment + SCORE_CAP + stability checks.

Stability check 2 (LOO R experiments) re-derives the signature per exclusion,
re-scores all clusters, re-computes calibration, and re-classifies T clusters —
so the flag is a *classification flip*, not just a score delta (spec §5 Step 4).
"""
from __future__ import annotations

import logging
from pathlib import Path

import numpy as np
import pandas as pd

from signature import (
    CORE_SUPPORT,
    FALLBACK_SUPPORT,
    TIMEPOINT_HOURS_CUTOFF,
    derive_for_ontology,
)

ANALYSIS_DIR = Path(__file__).resolve().parent.parent

SCORE_CAP = 10.0
NC_KEY_PATHWAY_EXCLUSION_PADJ = 1e-3  # Step 2 decide, decision #4 (option b)


def capped_signed_score(s: float) -> float:
    if pd.isna(s):
        return 0.0
    return float(np.sign(s)) * min(abs(s), SCORE_CAP)


def compute_score(cluster_rows: pd.DataFrame, signature: pd.DataFrame) -> float:
    """Layer A score (M2 = beta) for one cluster × ontology."""
    contributions = []
    sig_map = dict(zip(signature["term_id"], signature["direction"]))
    for term_id, ref_dir in sig_map.items():
        match = cluster_rows[cluster_rows["term_id"] == term_id]
        if match.empty:
            continue  # not tested in this cluster
        raw = match.iloc[0]["signed_score"]
        capped = capped_signed_score(raw)
        sign_ref = 1.0 if ref_dir == "up" else -1.0
        contributions.append(sign_ref * capped)
    return float(np.mean(contributions)) if contributions else np.nan


def score_all_clusters(
    enrichment_df: pd.DataFrame,
    signature_df: pd.DataFrame,
    class_map: dict[str, str],
) -> pd.DataFrame:
    """Score every (cluster, ontology, background_used) against signature_df."""
    rows = []
    for (cluster, ontology, bg), grp in enrichment_df.groupby(
        ["cluster", "ontology", "background_used"]
    ):
        ont_sig = signature_df[signature_df["ontology"] == ontology]
        if ont_sig.empty:
            continue
        score = compute_score(grp, ont_sig)
        first = grp.iloc[0]
        rows.append({
            "cluster": cluster,
            "ontology": ontology,
            "background_used": bg,
            "score": score,
            "organism": first["organism"],
            "experiment_id": first["experiment_id"],
            "timepoint": first.get("timepoint"),
            "timepoint_hours": first.get("timepoint_hours"),
            "direction": first.get("direction"),
            "class": class_map.get(first["experiment_id"]),
        })
    return pd.DataFrame(rows)


def disqualified_nc_clusters(
    enrichment_df: pd.DataFrame,
    key_panel: pd.DataFrame,
    nc_clusters_df: pd.DataFrame,
) -> set[tuple[str, str, str]]:
    """Return set of (cluster, ontology, background_used) tuples to exclude
    from NC calibration per Step 2 decide, decision #4 (option b).

    An NC cluster is disqualified for (ontology, bg) if it shows significant
    enrichment at padj < NC_KEY_PATHWAY_EXCLUSION_PADJ on ANY term in the a priori
    key-pathway panel for that ontology. Per-cluster granularity (up and down
    clusters at the same (exp, tp) are evaluated separately).
    """
    disqualified: set[tuple[str, str, str]] = set()
    for ontology, kp in key_panel.groupby("ontology"):
        key_ids = set(kp["term_id"])
        # Only consider enrichment rows in NC clusters × this ontology's key panel.
        nc_match = enrichment_df[
            enrichment_df["cluster"].isin(nc_clusters_df["cluster"])
            & (enrichment_df["ontology"] == ontology)
            & enrichment_df["term_id"].isin(key_ids)
            & (enrichment_df["p_adjust"] < NC_KEY_PATHWAY_EXCLUSION_PADJ)
        ]
        for (cluster, bg), _ in nc_match.groupby(["cluster", "background_used"]):
            disqualified.add((cluster, ontology, bg))
    return disqualified


def compute_calibration(
    scores_df: pd.DataFrame,
    nc_exclusions: set[tuple[str, str, str]] | None = None,
) -> pd.DataFrame:
    """NC noise floor + PC peak per (ontology, background_used).

    nc_exclusions: optional set of (cluster, ontology, background_used) tuples
    to drop from NC calibration (per decision #4). When None, no exclusion.
    """
    nc_exclusions = nc_exclusions or set()
    nc = scores_df[scores_df["class"] == "NC"].copy()
    # Drop disqualified NC clusters.
    if nc_exclusions:
        nc["_excl_key"] = list(zip(nc["cluster"], nc["ontology"], nc["background_used"]))
        nc = nc[~nc["_excl_key"].isin(nc_exclusions)].drop(columns=["_excl_key"])
    pc = scores_df[scores_df["class"] == "PC"]
    rows = []
    groups = set(zip(nc["ontology"], nc["background_used"])) | set(
        zip(pc["ontology"], pc["background_used"])
    )
    for ont, bg in sorted(groups):
        nc_grp = nc[(nc["ontology"] == ont) & (nc["background_used"] == bg)]
        pc_grp = pc[(pc["ontology"] == ont) & (pc["background_used"] == bg)]
        rows.append({
            "ontology": ont, "background_used": bg,
            "nc_mean": nc_grp["score"].mean() if len(nc_grp) else np.nan,
            "nc_std": nc_grp["score"].std(ddof=1) if len(nc_grp) >= 2 else np.nan,
            "nc_n": len(nc_grp),
            "pc_mean": pc_grp["score"].mean() if len(pc_grp) else np.nan,
            "pc_std": pc_grp["score"].std(ddof=1) if len(pc_grp) >= 2 else np.nan,
            "pc_n": len(pc_grp),
        })
    return pd.DataFrame(rows)


def classify(score: float, calib_df: pd.DataFrame, ontology: str, bg: str) -> str:
    """Map (score, matched calibration) → label per spec §5 Step 4."""
    c = calib_df[(calib_df["ontology"] == ontology) & (calib_df["background_used"] == bg)]
    if c.empty:
        return "no_matched_calibration"
    row = c.iloc[0]
    if pd.isna(row["nc_std"]) or row["nc_std"] == 0:
        return "insufficient_nc"
    if pd.isna(score):
        return "no_signature_coverage"
    if score >= row["nc_mean"] + 2 * row["nc_std"]:
        if (
            not pd.isna(row["pc_mean"])
            and not pd.isna(row["pc_std"])
            and score >= row["pc_mean"] - 2 * row["pc_std"]
        ):
            return "pc_like"
        return "detectable"
    if abs(score - row["nc_mean"]) < 2 * row["nc_std"]:
        return "no_signal"
    return "below_nc"


def rederive_signature_loo(
    enrichment_df: pd.DataFrame,
    r_exp_ids: set[str],
    excluded: str,
) -> pd.DataFrame:
    """Re-run Step 3 derivation logic excluding one R experiment."""
    r_subset = enrichment_df[
        enrichment_df["experiment_id"].isin(r_exp_ids - {excluded})
    ].copy()
    if "timepoint_hours" in r_subset.columns:
        hours = r_subset["timepoint_hours"]
        r_subset = r_subset[hours.isna() | (hours >= TIMEPOINT_HOURS_CUTOFF)]
    parts = []
    for ont, ont_df in r_subset.groupby("ontology"):
        s, _ = derive_for_ontology(ont_df, CORE_SUPPORT)
        if len(s) < 5:
            s, _ = derive_for_ontology(ont_df, FALLBACK_SUPPORT)
        s["ontology"] = ont
        parts.append(s)
    return pd.concat(parts, ignore_index=True) if parts else pd.DataFrame()


def main() -> int:
    logging.basicConfig(
        filename=str(ANALYSIS_DIR / "logs" / "step4.log"),
        level=logging.INFO,
        format="%(asctime)s [%(levelname)s] %(message)s",
    )
    log = logging.getLogger("step4")

    enrichment = pd.read_csv(ANALYSIS_DIR / "data" / "enrichment_all.csv")
    classified = pd.read_csv(ANALYSIS_DIR / "data" / "experiments_classified.csv")
    signature = pd.read_csv(ANALYSIS_DIR / "data" / "reference_signature.csv")
    key_panel = pd.read_csv(ANALYSIS_DIR / "exploration" / "key_pathways.csv")

    # Dedupe to one row per experiment — classified has one row per (experiment × timepoint).
    classified_unique = classified.drop_duplicates("experiment_id")
    class_map = classified_unique.set_index("experiment_id")["class"].to_dict()

    # Main scoring + calibration.
    scores = score_all_clusters(enrichment, signature, class_map)
    scores.to_csv(ANALYSIS_DIR / "results" / "scores_all.csv", index=False)
    log.info(f"Wrote {len(scores)} cluster scores")

    # NC calibration exclusion (Step 2 decide, decision #4 option b):
    # NC clusters with padj < NC_KEY_PATHWAY_EXCLUSION_PADJ on any a priori
    # key-pathway anchor are disqualified from the (ontology, bg) NC group.
    nc_clusters = scores[scores["class"] == "NC"][["cluster"]].drop_duplicates()
    nc_exclusions = disqualified_nc_clusters(enrichment, key_panel, nc_clusters)
    for (cluster, ontology, bg) in sorted(nc_exclusions):
        log.info(
            f"NC exclusion: cluster={cluster} ontology={ontology} "
            f"background_used={bg} (decision #4, padj<{NC_KEY_PATHWAY_EXCLUSION_PADJ})"
        )

    calib_df = compute_calibration(scores, nc_exclusions=nc_exclusions)
    log.info(
        f"Calibration: {len(calib_df)} (ontology, background_used) groups; "
        f"NC exclusions applied: {len(nc_exclusions)}"
    )
    # Persist the exclusion list as a results artifact for transparency.
    if nc_exclusions:
        pd.DataFrame(sorted(nc_exclusions),
                      columns=["cluster", "ontology", "background_used"]).to_csv(
            ANALYSIS_DIR / "results" / "nc_calibration_exclusions.csv",
            index=False,
        )

    # T condition summary with classification.
    t_scores = scores[scores["class"] == "T"].copy()
    t_summary = []
    for _, row in t_scores.iterrows():
        t_summary.append({
            **row.to_dict(),
            "classification": classify(row["score"], calib_df, row["ontology"], row["background_used"]),
        })
    summary = pd.DataFrame(t_summary)
    summary.to_csv(ANALYSIS_DIR / "results" / "score_summary.csv", index=False)

    # Build a lookup of original classification keyed by (cluster, ontology, background_used).
    orig_class_lookup = {
        (r["cluster"], r["ontology"], r["background_used"]): r["classification"]
        for r in t_summary
    }

    # Stability check 1: LOO signature pathways per T cluster.
    loo_sig_rows = []
    for _, trow in t_scores.iterrows():
        ont_sig = signature[signature["ontology"] == trow["ontology"]]
        cluster_rows = enrichment[
            (enrichment["cluster"] == trow["cluster"])
            & (enrichment["ontology"] == trow["ontology"])
        ]
        orig_score = trow["score"]
        for _, sigrow in ont_sig.iterrows():
            reduced = ont_sig[ont_sig["term_id"] != sigrow["term_id"]]
            new_score = compute_score(cluster_rows, reduced)
            flip = (np.sign(orig_score) != np.sign(new_score)) if not pd.isna(new_score) else False
            drop = (abs(new_score) < 0.5 * abs(orig_score)) if not pd.isna(new_score) else False
            loo_sig_rows.append({
                "cluster": trow["cluster"], "ontology": trow["ontology"],
                "removed_term_id": sigrow["term_id"], "removed_term_name": sigrow["term_name"],
                "orig_score": orig_score, "new_score": new_score,
                "flag_sign_flip": flip, "flag_large_drop": drop,
            })
    loo_sig = pd.DataFrame(loo_sig_rows)
    loo_sig.to_csv(ANALYSIS_DIR / "results" / "loo_signature.csv", index=False)

    # Stability check 2: LOO R experiments (classification-flip per spec §5 Step 4).
    # Per exclusion: re-derive signature, re-score ALL clusters, re-compute
    # calibration, re-classify T. Flag classification_flip if label changes.
    r_exp_ids = {eid for eid, c in class_map.items() if c == "R"}
    loo_r_rows = []
    for excl in sorted(r_exp_ids):
        new_sig = rederive_signature_loo(enrichment, r_exp_ids, excl)
        if new_sig.empty:
            log.warning(f"LOO-R excl={excl}: empty signature; skipping classification")
            for _, trow in t_scores.iterrows():
                orig_class = orig_class_lookup[(trow["cluster"], trow["ontology"], trow["background_used"])]
                loo_r_rows.append({
                    "excluded_experiment": excl,
                    "cluster": trow["cluster"], "ontology": trow["ontology"],
                    "background_used": trow["background_used"],
                    "orig_score": trow["score"], "new_score": np.nan,
                    "orig_classification": orig_class,
                    "new_classification": "empty_signature",
                    "classification_flip": True,
                    "new_signature_size": 0,
                })
            continue

        new_scores = score_all_clusters(enrichment, new_sig, class_map)
        # Apply the same NC exclusion under LOO — the disqualification is a
        # property of the NC cluster, not of the signature.
        new_calib = compute_calibration(new_scores, nc_exclusions=nc_exclusions)
        new_sig_size_by_ont = new_sig.groupby("ontology").size().to_dict()

        for _, trow in t_scores.iterrows():
            key = (trow["cluster"], trow["ontology"], trow["background_used"])
            orig_class = orig_class_lookup[key]
            match = new_scores[
                (new_scores["cluster"] == trow["cluster"])
                & (new_scores["ontology"] == trow["ontology"])
                & (new_scores["background_used"] == trow["background_used"])
            ]
            if match.empty:
                new_score = np.nan
                new_class = "no_signature_coverage"
            else:
                new_score = match.iloc[0]["score"]
                new_class = classify(new_score, new_calib, trow["ontology"], trow["background_used"])
            loo_r_rows.append({
                "excluded_experiment": excl,
                "cluster": trow["cluster"], "ontology": trow["ontology"],
                "background_used": trow["background_used"],
                "orig_score": trow["score"], "new_score": new_score,
                "orig_classification": orig_class,
                "new_classification": new_class,
                "classification_flip": orig_class != new_class,
                "new_signature_size": int(new_sig_size_by_ont.get(trow["ontology"], 0)),
            })
    loo_r = pd.DataFrame(loo_r_rows)
    loo_r.to_csv(ANALYSIS_DIR / "results" / "loo_r_experiments.csv", index=False)

    n_flips = int(loo_r["classification_flip"].sum()) if not loo_r.empty else 0
    log.info(f"LOO-R: {len(loo_r)} (excl × T cluster × ontology) rows, {n_flips} classification flips")

    # Print summary.
    print("\n=== T condition scores with classification ===")
    print(summary[[
        "experiment_id", "timepoint", "direction", "ontology", "background_used",
        "score", "classification",
    ]].to_string(index=False))
    print(f"\n=== LOO-R classification flips: {n_flips} / {len(loo_r)} ===")
    if n_flips:
        print(loo_r[loo_r["classification_flip"]][[
            "excluded_experiment", "cluster", "ontology", "background_used",
            "orig_classification", "new_classification",
        ]].to_string(index=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
```

**Classification flip semantics (per spec §5 Step 4):** the LOO-R flag fires when removing any single R experiment changes a T cluster's threshold label (`detectable` / `no_signal` / `pc_like` / `below_nc` / `insufficient_nc` / `no_matched_calibration`). This requires re-running the full pipeline (derive → score all → calibrate → classify) per exclusion — `rederive_signature_loo` + `score_all_clusters` + `compute_calibration` + `classify`. Score-only flags (sign flip, >50% drop) are kept as Stability check 1 (LOO signature pathways), which is a different failure mode.

- [ ] **Step 2: Run the script**

```bash
uv run "$ANALYSIS_DIR/scripts/05_compute_scores.py"
```

Expected: `results/scores_all.csv`, `score_summary.csv`, `loo_signature.csv`, `loo_r_experiments.csv`. Stdout shows T summary with classification per ontology per background.

- [ ] **Step 3: Update RESULTS_MANIFEST.md**

- [ ] **Step 4: Check that pre-registered expectations match observed outcomes**

Compare decisions.md pre-registration against scores. Document discrepancies in notebook (decisions.md is ALREADY written; don't edit post-hoc).

- [ ] **Step 5: Commit 1 (do-phase, bundling pre-registration + scoring)**

**Before committing — verify `decisions.md` is ready to stage.** It was written in Task 9 Step 1 and deliberately left uncommitted so it lands in a single atomic commit alongside the scoring outputs (spec §5 Step 4). Confirm:

```bash
git status --short -- "$ANALYSIS_DIR/decisions.md"
# Expected: either "?? ...decisions.md" (untracked) or " M ...decisions.md" (modified).
# If empty, the pre-registration was never written — go back to Task 9.
```

Then commit:

```bash
git add "$ANALYSIS_DIR/scripts/05_compute_scores.py" \
        "$ANALYSIS_DIR/results/scores_all.csv" \
        "$ANALYSIS_DIR/results/score_summary.csv" \
        "$ANALYSIS_DIR/results/loo_signature.csv" \
        "$ANALYSIS_DIR/results/loo_r_experiments.csv" \
        "$ANALYSIS_DIR/results/RESULTS_MANIFEST.md" \
        "$ANALYSIS_DIR/decisions.md" \
        "$ANALYSIS_DIR/logs/step4.log"
# Conditional: nc_calibration_exclusions.csv is written only if any NC cluster
# gets disqualified. Include if present.
if [ -f "$ANALYSIS_DIR/results/nc_calibration_exclusions.csv" ]; then
    git add "$ANALYSIS_DIR/results/nc_calibration_exclusions.csv"
fi
git commit -m "step 4 do: pre-registration + scores + stability checks"
```

After commit, verify `decisions.md` was included: `git log -1 --stat | grep decisions.md` must return a line; if not, the commit is malformed and must be redone.

---

## Task 11 — Step 4 show / explore / decide

**Spec reference:** §5 Step 4 show/explore/decide

- [ ] **Step 1: Present score summary in chat**

- Per T condition × ontology × background: score + classification.
- NC calibration per (ontology, background).
- LOO signature pathway flags (any score flips or >50% drops).
- LOO R experiment stability (does any single exclusion flip classification?).
- Cross-ontology agreement per T condition.

- [ ] **Step 2: Write explore scripts for contribution decomposition**

`$ANALYSIS_DIR/scripts/explore_step4_contribution_decomposition.py` — per T cluster, show top-5 contributing pathways by absolute magnitude, with `.explain()` drill-down on the top contributors.

- [ ] **Step 3: Append show + exploration to notebook**

- [ ] **Step 4: Interactive walkthrough with researcher**

Per-T chat capture:
- Does the score match the pre-registered expectation?
- Top contributing pathways — canonical drivers?
- Any LOO instability?
- Cross-ontology disagreement?

- [ ] **Step 5: Researcher decides (scores final / redo)**

- [ ] **Step 6: Commit 2 (decide-phase)**

```bash
git add "$ANALYSIS_DIR/exploration/notebook.md" \
        "$ANALYSIS_DIR/scripts/"explore_step4_*.py 2>/dev/null || true
git commit -m "step 4 decide: scores approved"
```

---

## Task 12 — Step 5 do: figures + write-up

**Spec reference:** §5 Step 5 do

**Files:**
- Create: `$ANALYSIS_DIR/scripts/06_make_figures.py`
- Create: `$ANALYSIS_DIR/results/fig1_heatmap.png`, `.pdf`
- Create: `$ANALYSIS_DIR/results/fig2_trajectories.png`, `.pdf`
- Create: `$ANALYSIS_DIR/methods.md`
- Create: `$ANALYSIS_DIR/caveats.md`
- Create: `$ANALYSIS_DIR/gaps_and_friction.md`
- Create: `$ANALYSIS_DIR/api_coverage.md`
- Create: `$ANALYSIS_DIR/README.md`
- Create: `$ANALYSIS_DIR/references.md`
- Create: `$ANALYSIS_DIR/logs/step5.log`

- [ ] **Step 1: Write the figures script (skeleton)**

Write `$ANALYSIS_DIR/scripts/06_make_figures.py`. **This is a starting skeleton.** After the first run, inspect both figures and iterate on sizes/labels/highlights with the researcher in the show phase — expect 2–3 iterations before final. The structural logic (row selection policy, column ordering, color cap) is spec-determined and should not drift; everything marked `# TWEAK:` is cosmetic and will change based on actual data.

Spec refs: Fig 1 = spec §5 Step 5 "unified signed-enrichment heatmap"; Fig 2 = spec §5 Step 5 "score-trajectory lineplots". `SCORE_CAP = 10` matches the scoring cap.

```python
"""Build Fig 1 (unified signed-enrichment heatmap) and Fig 2 (score trajectories).

Starting skeleton — cosmetic choices (fontsize, row density, highlight styling)
are marked `# TWEAK:` and will be iterated in Step 5 show/explore/decide.
Structural choices (axis semantics, row selection, column ordering, color cap)
are spec-locked and should not drift.
"""
from __future__ import annotations

import logging
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

ANALYSIS_DIR = Path(__file__).resolve().parent.parent

SCORE_CAP = 10.0  # matches scoring cap — color scale saturates here
CLASS_ORDER = ["T", "R", "PC", "NC", "CTX"]
FIG1_ROW_CAP = 80  # TWEAK: may need to drop to ~60 for legibility, or switch to multi-panel


def load_inputs() -> dict:
    enrichment = pd.read_csv(ANALYSIS_DIR / "data" / "enrichment_all.csv")
    signature = pd.read_csv(ANALYSIS_DIR / "data" / "reference_signature.csv")
    key_paths = pd.read_csv(ANALYSIS_DIR / "exploration" / "key_pathways.csv")
    scores = pd.read_csv(ANALYSIS_DIR / "results" / "scores_all.csv")
    summary = pd.read_csv(ANALYSIS_DIR / "results" / "score_summary.csv")
    classified = pd.read_csv(ANALYSIS_DIR / "data" / "experiments_classified.csv")
    class_map = dict(zip(classified["experiment_id"], classified["class"]))
    enrichment["class"] = enrichment["experiment_id"].map(class_map)
    return dict(
        enrichment=enrichment, signature=signature, key_paths=key_paths,
        scores=scores, summary=summary, classified=classified, class_map=class_map,
    )


def select_fig1_rows(inputs: dict) -> pd.DataFrame:
    """Row selection policy (spec-locked): signature ∪ key-pathway panel ∪ any
    pathway sig-enriched in ≥1 T cluster."""
    enrichment = inputs["enrichment"]
    signature = inputs["signature"]
    key_paths = inputs["key_paths"]

    sig_ids = set(signature["term_id"])
    key_ids = set(key_paths["term_id"])
    t_sig_ids = set(
        enrichment[(enrichment["class"] == "T") & (enrichment["p_adjust"] < 0.05)]["term_id"]
    )
    visible = sig_ids | key_ids | t_sig_ids

    rows = (
        enrichment[enrichment["term_id"].isin(visible)]
        .drop_duplicates(["term_id", "term_name", "ontology"])
        [["term_id", "term_name", "ontology"]]
    )
    rows["in_signature"] = rows["term_id"].isin(sig_ids)
    rows["in_key_panel"] = rows["term_id"].isin(key_ids)
    return rows


def order_columns(enrichment: pd.DataFrame) -> list[str]:
    """Column ordering: T → R → PC → NC → CTX; within class by experiment, timepoint."""
    cluster_df = (
        enrichment.drop_duplicates("cluster")
        [["cluster", "class", "experiment_id", "timepoint", "timepoint_hours", "direction", "organism"]]
    )
    cluster_df["class_rank"] = cluster_df["class"].map({c: i for i, c in enumerate(CLASS_ORDER)}).fillna(99)
    cluster_df = cluster_df.sort_values(
        ["class_rank", "experiment_id", "timepoint_hours", "timepoint", "direction"]
    )
    return cluster_df["cluster"].tolist()


def build_fig1_pivot(enrichment: pd.DataFrame, visible_rows: pd.DataFrame,
                     column_order: list[str]) -> pd.DataFrame:
    sub = enrichment[enrichment["term_id"].isin(visible_rows["term_id"])]
    pivot = sub.pivot_table(
        index=["ontology", "term_name", "term_id"],
        columns="cluster",
        values="signed_score",
        aggfunc="first",
    )
    # Preserve column order.
    pivot = pivot.reindex(columns=[c for c in column_order if c in pivot.columns])
    # Saturate at ±SCORE_CAP for the visualization (raw values preserved in CSV).
    pivot = pivot.clip(lower=-SCORE_CAP, upper=SCORE_CAP)
    return pivot


def make_fig1(inputs: dict) -> None:
    visible = select_fig1_rows(inputs)
    column_order = order_columns(inputs["enrichment"])
    pivot = build_fig1_pivot(inputs["enrichment"], visible, column_order)

    if len(pivot) > FIG1_ROW_CAP:
        # TWEAK: if overflow, switch to per-ontology multi-panel. For now, warn + truncate.
        logging.warning(
            f"Fig 1: {len(pivot)} rows exceeds cap {FIG1_ROW_CAP}; "
            f"keeping top by max |signed_score|. Consider multi-panel."
        )
        rank = pivot.abs().max(axis=1).sort_values(ascending=False)
        pivot = pivot.loc[rank.head(FIG1_ROW_CAP).index]

    # TWEAK: figsize proportional to matrix shape; adjust after first run.
    fig, ax = plt.subplots(figsize=(max(14, 0.25 * pivot.shape[1]),
                                    max(8, 0.22 * pivot.shape[0])))
    sns.heatmap(
        pivot, center=0, cmap="RdBu_r", vmin=-SCORE_CAP, vmax=SCORE_CAP,
        ax=ax, cbar_kws={"label": f"signed_score (saturated at ±{SCORE_CAP:.0f})"},
    )
    # TWEAK: signature row highlighting. Simplest form — bold + left annotation.
    # After first run, may switch to shaded band or separate annotation column.
    sig_ids = set(inputs["signature"]["term_id"])
    for i, (_, _, term_id) in enumerate(pivot.index):
        if term_id in sig_ids:
            ax.get_yticklabels()[i].set_fontweight("bold")
    # TWEAK: MED4 column highlighting — bold/colored x-tick labels.
    med4_clusters = set(
        inputs["enrichment"]
        .query("organism == 'Prochlorococcus MED4'")["cluster"]
        .unique()
    )
    for i, col in enumerate(pivot.columns):
        if col in med4_clusters:
            ax.get_xticklabels()[i].set_color("#0b5394")
    ax.set_title("Fig 1 — pathway enrichment across classes (signed score)")
    plt.setp(ax.get_xticklabels(), rotation=90, fontsize=7)  # TWEAK: rotation / fontsize
    plt.setp(ax.get_yticklabels(), fontsize=8)               # TWEAK
    plt.tight_layout()
    out = ANALYSIS_DIR / "results" / "fig1_heatmap"
    plt.savefig(out.with_suffix(".png"), dpi=200, bbox_inches="tight")
    plt.savefig(out.with_suffix(".pdf"), bbox_inches="tight")
    plt.close(fig)


def make_fig2(inputs: dict) -> None:
    scores = inputs["scores"].copy()
    classified = inputs["classified"]

    # Precompute NC noise-floor reference lines per (ontology, background_used).
    nc = scores[scores["class"] == "NC"]
    calib = (
        nc.groupby(["ontology", "background_used"])["score"]
        .agg(mean="mean", std="std").reset_index()
    )

    experiments = sorted(
        scores["experiment_id"].dropna().unique(),
        key=lambda eid: (
            CLASS_ORDER.index(inputs["class_map"].get(eid, "CTX"))
            if inputs["class_map"].get(eid, "CTX") in CLASS_ORDER else 99,
            eid,
        ),
    )
    n = len(experiments)
    # TWEAK: panel grid shape. 4 cols is a decent default for 10–20 experiments.
    ncols = 4
    nrows = (n + ncols - 1) // ncols
    fig, axes = plt.subplots(nrows, ncols,
                             figsize=(4 * ncols, 3 * nrows),  # TWEAK
                             squeeze=False)
    ontologies = sorted(scores["ontology"].dropna().unique())
    palette = dict(zip(ontologies, sns.color_palette("tab10", len(ontologies))))

    for ax_idx, eid in enumerate(experiments):
        ax = axes[ax_idx // ncols][ax_idx % ncols]
        sub = scores[scores["experiment_id"] == eid].copy()
        cls = inputs["class_map"].get(eid, "CTX")
        for ont, g in sub.groupby("ontology"):
            g = g.sort_values("timepoint_hours", na_position="first")
            x = g["timepoint_hours"].fillna(0).to_numpy() if g["timepoint_hours"].notna().any() \
                else np.arange(len(g))
            ax.plot(x, g["score"], marker="o", color=palette.get(ont), label=ont)

        # NC reference band per ontology.
        for _, row in calib.iterrows():
            if pd.isna(row["mean"]) or pd.isna(row["std"]):
                continue
            ax.axhspan(row["mean"] - 2 * row["std"], row["mean"] + 2 * row["std"],
                       alpha=0.08, color=palette.get(row["ontology"], "grey"))
        ax.axhline(0, color="grey", linewidth=0.5)
        title = f"{eid[:40]}…\n[{cls}]"  # TWEAK: truncate length
        ax.set_title(title, fontsize=8)  # TWEAK
        ax.tick_params(labelsize=7)
        # TWEAK: frame MED4 T panels specifically.
        if cls == "T":
            for spine in ax.spines.values():
                spine.set_edgecolor("#0b5394")
                spine.set_linewidth(1.5)

    # Hide empty subplots.
    for j in range(n, nrows * ncols):
        axes[j // ncols][j % ncols].axis("off")

    handles = [plt.Line2D([0], [0], color=palette[o], marker="o", label=o) for o in ontologies]
    fig.legend(handles=handles, loc="upper right", bbox_to_anchor=(1.0, 1.0))
    fig.suptitle("Fig 2 — Layer A score trajectories by experiment", fontsize=12)
    plt.tight_layout(rect=[0, 0, 0.95, 0.97])
    out = ANALYSIS_DIR / "results" / "fig2_trajectories"
    plt.savefig(out.with_suffix(".png"), dpi=200, bbox_inches="tight")
    plt.savefig(out.with_suffix(".pdf"), bbox_inches="tight")
    plt.close(fig)


def main() -> int:
    logging.basicConfig(
        filename=str(ANALYSIS_DIR / "logs" / "step5.log"),
        level=logging.INFO,
        format="%(asctime)s [%(levelname)s] %(message)s",
    )
    inputs = load_inputs()
    make_fig1(inputs)
    make_fig2(inputs)
    logging.info("Figures written to results/ (png + pdf)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
```

**Iteration protocol (spec §5 Step 5 show/explore/decide):** first run produces `fig1_heatmap.png` and `fig2_trajectories.png` at reasonable-but-unpolished quality. Show in chat; researcher flags which `# TWEAK:` spots need adjustment (row cap, panel grid, fontsize, highlight styling). Commit each iteration with its notebook entry. Do not treat the first-run output as final.

- [ ] **Step 2: Run the figures script**

```bash
uv run "$ANALYSIS_DIR/scripts/06_make_figures.py"
```

Expected: two PNG + two PDF/SVG in `results/`.

- [ ] **Step 3: Write `methods.md`**

Per artifacts-guide §methods.md, required sections:
- Research question (operationalized).
- Data scope (organisms, DOIs, experiment IDs, inclusions/exclusions, KG version from `kg_schema()`).
- Gene selection (how identified via search_ontology/genes_by_ontology, key pathway panel).
- Statistical methods (Fisher's ORA via pathway_enrichment, BH per cluster, signed_score cap at SCORE_CAP=10, sign-weighted alignment score M2, NC noise floor calibration per (ontology, background_used)).
- Results summary (per-T-condition scores + classifications).
- Limitations (pulled from caveats.md).

- [ ] **Step 4: Write `caveats.md`**

Per artifacts-guide §caveats.md. Include:
- Cross-background incomparability (fold_enrichment magnitudes differ, signed_score still works).
- SCORE_CAP effect on extreme-enrichment pathways.
- Temporal filter exclusions.
- Any fragile signature members (single-experiment dominance from Step 3).
- Any T conditions with no matched calibration.
- B1 caveats that carry over (catch-all categories).

- [ ] **Step 5: Write `gaps_and_friction.md`**

Per artifacts-guide §gaps_and_friction.md format. Four sections: KG data bugs, KG gaps, MCP friction, Skill/methodology friction. Plus Process retrospective (what worked / didn't / proposed changes).

**Pre-seeded items** (carry these over from pre-execution review + Step 1a execution; expand with anything new encountered during execution). Source: `docs/superpowers/specs/2026-04-18-research-methodology-v3-improvements-from-b2.md` §4.

### MCP friction

- **`list_experiments` lacks an `experiment_ids` parameter.** Once the researcher has classified a concrete set of IDs, there's no direct way to fetch metadata only for those IDs — the workaround is a full-landscape pull + local filter. Proposed fix: add `experiment_ids: list[str] | None = None` to `list_experiments`, consistent with `pathway_enrichment`, `ontology_landscape`, and `gene_overview`. (v3 meta-doc §4.5 A1)
- **`pathway_enrichment` "Package import equivalent" doc says it returns a dict, but it actually returns an `EnrichmentResult` object.** Caused real confusion during pre-execution review (led to a retracted blocker-bug report about missing columns). Propose: fix the MCP doc to say "returns `EnrichmentResult`; call `.to_envelope(...)` for the MCP dict shape." (v3 meta-doc §4.5 A2)
- **Cluster naming convention `{experiment_id}|{timepoint}|{direction}` (NaN → `"NA"`) is not prominent in the docs.** Every `.explain()` drill-down depends on it. Propose: move to top of Response format section. (v3 meta-doc §4.5 A3)
- **`ontology_landscape` default `limit=10` silently truncates.** Surfaces as "landscape looks thin" only if you know to look for `truncated=True`. Propose: change default to `None` or highlight truncation more prominently. (v3 meta-doc §4.5 A4)

### Skill / methodology friction

- **Plan-review surfaced hygiene issues at ~4:1 ratio to analysis-specific bugs.** The hygiene patterns (F1–F4 in v3 meta-doc §4.1) aren't captured as explicit rules anywhere in the skill — propose a dedicated `research-code-hygiene.md` with the four failure modes and their corrections.
- **Empirical probe rule missing from step-protocol.** Spec §8 risk 7 had four escalation paths for pickle failure; a 10-minute empirical probe ruled all of them out. When any risk in a spec has ≥3 escalation paths, the plan must probe empirically before planning contingencies. (v3 meta-doc §4.2)
- **Spec-to-plan parameterization drift.** Spec §5.0 hardcoded a dated notebook filename that went stale when the plan ran a day later. Specs should reference parameters (`<target_organism>`, `<notebook_filename>`), not concrete values. (v3 meta-doc §4.3)
- **Figures iterate within a step; don't spec cosmetics upfront.** Already C19 — reinforced during review. Figure code belongs in skeleton + `# TWEAK:` form, not polished inline. (v3 meta-doc §4.4)

### Process retrospective — plan-review effectiveness

- Pre-execution plan review caught 20+ issues in one multi-pass session: one blocker (broken organism heuristic), several high-severity silent-wrong-behavior patterns (duplicated experiment_ids, dict-zip on multi-row frames, hardcoded organism strings, empty-intermediate blindness), and a handful of medium brittleness items (ANALYSIS_DIR env lifetime, shell globs in git add, nested markdown fences).
- ~80% of issues caught were generalizable hygiene patterns, not analysis-specific. Proposes: a plan-review sub-skill or a code-reviewer subagent primed with the v3 meta-doc §4.6 hygiene checklist could mechanize most of this work.
- Track during B2 execution: how many new hygiene issues surface that the §4.6 checklist missed? Any that show up is a candidate new rule.

- [ ] **Step 6: Write `api_coverage.md`**

Per spec §10. Table of MCP tools + Python primitives used, for what, where each worked or failed. Compared against B1's friction list (look at `analyses/2026-04-09-1713-pathway_enrichment_b1/gaps_and_friction.md`).

**Pre-seeded gaps to include** (from pre-execution review + execution). Source: v3 meta-doc §4.5.

| Tool / Function | Step used | Purpose | Worked well | Friction / gaps |
|---|---|---|---|---|
| `list_experiments` | 1a | Pull metadata for the classified experiment IDs | Filters (organism, treatment_type, search_text) work; single unfiltered `limit=None` call returns the full landscape cheaply | **No `experiment_ids` filter.** Cannot pass a list of specific IDs to fetch only those rows. Workaround: pull unfiltered and filter locally by `experiment_id`. Not costly today (~76 experiments), but as the KG grows this becomes wasteful. Proposed: add an `experiment_ids: list[str] \| None` parameter mirroring the filter shape of `pathway_enrichment` / `ontology_landscape`. (§4.5 A1) |
| `pathway_enrichment` | 2 | Run Fisher ORA per (org, ontology, background) | Returns `EnrichmentResult` with `.results` DataFrame, `.explain()`, `.generate_summary()`, `.to_envelope()`. Cluster metadata (experiment_id, timepoint, timepoint_hours, direction, etc.) is merged into `.results` by the Python wrapper — usable directly for downstream filtering. | **Doc says "returns dict" in "Package import equivalent" block** — actually returns `EnrichmentResult`. Caused a retracted blocker-bug report during plan review. (§4.5 A2)<br>**Cluster naming convention not prominent in response docs.** `{experiment_id}\|{timepoint}\|{direction}`, NaN → `"NA"` — documented only in a common-mistakes bullet. (§4.5 A3) |
| `ontology_landscape` | 1b | Rank (ontology, level) by enrichment suitability | Per-organism ranking is clear; `experiment_ids` filter works; `relevance_rank` column stable. | **Default `limit=10` silently truncates.** Surfaces only via `truncated=True` in envelope. (§4.5 A4) |
| `search_ontology` | 1b | Find nitrogen-related terms across ontologies | Returns a clean flat list, `to_dataframe`-friendly. | [Fill in during execution.] |
| `genes_by_ontology` | 1b (sanity check) | Verify canonical marker genes in selected ontology term_ids | Used in the key-pathways sanity loop; returns `term_id → locus_tag` mapping. | [Fill in during execution.] |
| `EnrichmentResult.explain` | 2, 3, 4 | Drill into per-pathway gene overlap | Works on loaded-from-pickle instances (empirically verified on 2026-04-19). | [Fill in during execution — any cases where cluster or term_id isn't found cleanly.] |
| Standard pickle round-trip | 2 | Persist `EnrichmentResult` dict for downstream drill-down | Works cleanly on 2-key heterogeneous dict at B2 scale (~0.4 MB/obj, ~20 MB full dict). | Spec §8 risk 7 contingencies are precautionary only; empirically validated (§4.2). |
| `fisher_ora` (direct Python primitive) | (not used in B2) | Would enable custom-term2gene ORA | Available via `multiomics_explorer.analysis.enrichment`. | Not exercised in B2 — B2 uses the DE-wired `pathway_enrichment` wrapper. Left for a future analysis. |

- [ ] **Step 7: Write `README.md`**

Top-level analysis README:
- Question + one-paragraph summary.
- Key findings bullet list.
- Links to methods, decisions, caveats, gaps, notebook.
- Reproduction instructions (`uv run scripts/01_...` through `06_...`).
- Artifact layout.
- Per-T-condition score summary table (from `score_summary.csv`).

- [ ] **Step 8: Write `references.md`**

Bibliographic — Weissberg 2025 DOI, Tolonen and other reference papers (pull DOIs from experiments_classified.csv + `list_publications`), methodology papers (clusterProfiler, etc.), tool/KG citations.

- [ ] **Step 9: Update RESULTS_MANIFEST.md with figure rows**

- [ ] **Step 10: Commit 1 (do-phase)**

```bash
git add "$ANALYSIS_DIR/scripts/06_make_figures.py" \
        "$ANALYSIS_DIR/results/fig1_heatmap.png" "$ANALYSIS_DIR/results/fig1_heatmap.pdf" \
        "$ANALYSIS_DIR/results/fig2_trajectories.png" "$ANALYSIS_DIR/results/fig2_trajectories.pdf" \
        "$ANALYSIS_DIR/results/RESULTS_MANIFEST.md" \
        "$ANALYSIS_DIR/methods.md" "$ANALYSIS_DIR/caveats.md" \
        "$ANALYSIS_DIR/gaps_and_friction.md" "$ANALYSIS_DIR/api_coverage.md" \
        "$ANALYSIS_DIR/README.md" "$ANALYSIS_DIR/references.md" \
        "$ANALYSIS_DIR/logs/step5.log"
git commit -m "step 5 do: figures and write-up"
```

---

## Task 13 — Step 5 show / explore / decide

**Spec reference:** §5 Step 5 show/explore/decide

- [ ] **Step 1: Present final figures and README summary table in chat**

- [ ] **Step 2: Append show section to notebook**

- [ ] **Step 3: Interactive walkthrough with researcher**

- MED4 distinguishable in Fig 1?
- Trajectories readable in Fig 2?
- Story clarity: does the combination of figures + summary answer the research question?
- Anything missing from methods/caveats?

- [ ] **Step 4: Researcher decides: publish**

- [ ] **Step 5: Commit 2 (decide-phase)**

```bash
git add "$ANALYSIS_DIR/exploration/notebook.md"
git commit -m "step 5 decide: analysis published"
```

---

## Task 14 — Final wrap-up + meta-doc merge

- [ ] **Step 1: Verify all manifests are complete**

Every file in `data/` and `results/` should have a row in the corresponding manifest. Check and add any missing rows.

- [ ] **Step 2: Merge analysis-local `gaps_and_friction.md` findings into the meta-doc**

Open `docs/superpowers/specs/2026-04-18-research-methodology-v3-improvements-from-b2.md` and append findings from `$ANALYSIS_DIR/gaps_and_friction.md` under the appropriate sections:

- **§2 (content gaps)** — new C-numbered items for missing methodology rules.
- **§3 (meta-process gaps)** — new M-lettered items for framing/drift issues.
- **§4 (B2 plan-review conclusions)** — add any new execution-time hygiene patterns that weren't caught by the §4.6 checklist. Update §4.7 retrospective with the "how did §4.6 hold up" observations.
- **§4.5 (API gaps)** — any new MCP / explorer API friction encountered during execution.

Note: §4 was seeded before execution from the plan review. The pre-seeded items (A1–A4 API gaps, F1–F4 code hygiene, empirical probe rule, parameterization drift) are already in the meta-doc — the analysis-local `gaps_and_friction.md` should focus on NEW friction encountered, not re-documenting the pre-seeded items. Treat the pre-seed as the baseline; record the delta.

- [ ] **Step 3: Commit meta-doc update**

```bash
git add docs/superpowers/specs/2026-04-18-research-methodology-v3-improvements-from-b2.md
git commit -m "meta: merge B2 execution friction findings into skill-v3 proposals"
```

- [ ] **Step 4: Tell the researcher the analysis is complete**

Report: analysis directory path, final commit hash, location of figures and summary table, any open items or follow-ups to discuss.

**Do not push to `origin` without explicit researcher approval.** `git push origin main` affects shared state; the researcher decides when the branch is ready to share. If they approve, push; otherwise leave the commits local.

---

## Self-review checklist (for the plan writer, not the executor)

- [ ] Every spec section has at least one task implementing it.
- [ ] No placeholders, no "TBD", no "similar to Task N" — every step has complete content.
- [ ] Method signatures are consistent (e.g., `compute_score()` signature stable across tasks).
- [ ] File paths are exact; exploration scripts use `scripts/explore_*.py` convention.
- [ ] Each step is bite-sized where possible; research-interactive steps are larger but still list concrete sub-actions.
- [ ] Every step-protocol phase (do/show/explore/decide) is a separate task — not collapsed.
- [ ] Three hard gates referenced (step boundary, manifest currency, chat-capture).
- [ ] Pickle round-trip + fallbacks addressed (§8 risk 7).
- [ ] Cross-background calibration addressed.
- [ ] SCORE_CAP documented.
- [ ] Temporal filter documented.
- [ ] Bidirectional tie-breaker documented.
- [ ] Per-experiment breakdown included.
- [ ] Dropped-notable log included.
- [ ] Pre-registration bundled with scoring commit.
- [ ] Meta-doc merge at end.
