# N-Limitation Signature Analysis v2 — Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Rebuild the N-limitation signature analysis with proper notebook discipline, producing a clean `sig_utils` package and interactive walkthrough at every step.

**Architecture:** Seven pipeline steps (discover → extract references → build signature → extract targets → design metric → score → interpret), each with do→show→explore→decide. Reusable methodology lives in `sig_utils/` (extraction, signature construction, scoring); scripts handle orchestration and diagnostics. One markdown notebook logs all steps; optional Jupyter companion at the end.

**Tech Stack:** Python 3.11+, pandas, numpy, scipy, matplotlib, multiomics_explorer Python API, uv

**Spec:** `docs/superpowers/specs/2026-04-07-n-limitation-signature-v2-design.md`

---

## File Structure

```
analyses/2026-04-08-HHMM-n_limitation_signature_v2/
├── exploration/
│   └── 2026-04-08-notebook.md
├── data/
│   └── DATA_MANIFEST.md
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
├── results/
│   └── RESULTS_MANIFEST.md
├── superpowers/
│   ├── spec.md
│   ├── plan.md
│   └── brainstorm-log.md
├── README.md
├── methods.md
├── decisions.md
├── caveats.md
├── gaps_and_friction.md
└── references.md
```

**Note on `HHMM`:** The scaffolding task fills this with the actual time at creation.

**Note on the research notebook:** The notebook (`exploration/2026-04-08-notebook.md`) is the primary deliverable — it's where each step's QC, exploration, and decisions are recorded. Every task that produces data includes a notebook entry. The **explore** and **decide** phases happen interactively in chat with the researcher, then get captured in the notebook.

**Note on marker genes:** Five genes are traced through every step: glnA (PMM0920, up), cynA (PMM0370, up), rbcL (PMM0550, down), atpD (PMM1452, down), PMM0030 (up). Edge-case tracers added as encountered.

---

### Task 1: Scaffold the analysis directory

**Files:**
- Create: all directories and stub files in the file structure above

- [ ] **Step 1: Create the directory tree**

```bash
TIMESTAMP=$(date +%H%M)
ANALYSIS_DIR="analyses/2026-04-08-${TIMESTAMP}-n_limitation_signature_v2"

mkdir -p "${ANALYSIS_DIR}"/{exploration,data,scripts,sig_utils/tests,logs,results,superpowers}

touch "${ANALYSIS_DIR}/sig_utils/__init__.py"
```

- [ ] **Step 2: Create stub documents**

Create these files with their headers only (content filled during analysis):

`README.md`:
```markdown
# N-Limitation Signature Analysis v2

**Question:** Can we quantify nitrogen limitation in Prochlorococcus MED4 molecularly, using a gene signature from independent reference studies?

**Approach:** Methodological rebuild of v1 with proper notebook discipline, single primary metric (rank score), positive/negative controls, and package-ready sig_utils.

**Spec:** `superpowers/spec.md`
**Plan:** `superpowers/plan.md`
**Brainstorm:** `superpowers/brainstorm-log.md`
**Notebook:** `exploration/2026-04-08-notebook.md`
**Predecessor:** `analyses/2026-04-06-1432-n_limitation_signature/`
```

`methods.md`:
```markdown
# N-Limitation Signature Analysis v2: Methods

*Living document — updated incrementally as findings become established.*
```

`decisions.md`:
```markdown
# Decision Log

Design decisions with rationale — WHY the analysis was done this way, not what it does (that's in methods.md).
```

`caveats.md`:
```markdown
# Caveats for Interpretation

Things a reader of these results needs to know before drawing conclusions.
```

`gaps_and_friction.md`:
```markdown
# Gaps and friction points
```

`references.md`:
```markdown
# References and Citations

## Publications
*Populated during Step 1 (KG discovery).*

## Software
- multiomics_explorer Python API
- Neo4j knowledge graph
- Python packages: pandas, numpy, scipy, matplotlib
```

`data/DATA_MANIFEST.md`:
```markdown
# Data Manifest

All files produced by extraction scripts from the multiomics KG.
Run scripts from `multiomics_research` root with `uv run`.
```

`results/RESULTS_MANIFEST.md`:
```markdown
# Results Manifest

All files produced by scoring, plotting, and analysis scripts.
```

`exploration/2026-04-08-notebook.md`:
```markdown
# Research Notebook: N-Limitation Signature v2

Chronological, append-only lab notebook. Each step follows do→show→explore→decide.

---
```

- [ ] **Step 3: Copy superpowers products into the analysis**

```bash
cp docs/superpowers/specs/2026-04-07-n-limitation-signature-v2-design.md "${ANALYSIS_DIR}/superpowers/spec.md"
cp docs/superpowers/plans/2026-04-08-n-limitation-signature-v2.md "${ANALYSIS_DIR}/superpowers/plan.md"
```

- [ ] **Step 4: Write the brainstorm log**

Create `superpowers/brainstorm-log.md` capturing the Q&A from the brainstorming session:

```markdown
# Brainstorming Log: N-Limitation Signature v2

**Date:** 2026-04-07 to 2026-04-08
**Participants:** Researcher + Claude

## Context

The original analysis (v1) at `analyses/2026-04-06-1432-n_limitation_signature/` completed Approach A but:
- Methodology skill loaded too late (retrofitted)
- Results are a black box — researcher can't follow or defend the method
- Process shortcuts (skipped reviews, untested plan code)

Goals for v2:
1. Clean method suitable for packaging
2. Deep understanding of the signature scoring approach
3. Correct methodology from the start

## Decisions

### Q1: Redo or walk-through?
**Options:** (A) Walk through existing scripts/outputs, fixing as we go. (B) Rebuild from scratch, reusing concepts.
**Decision:** B — rebuild from scratch. Cleaner for packaging, avoids anchoring to potentially wrong decisions.

### Q2: Scope — which parts?
**Options:** All four parts (extraction, signature, scoring, Alteromonas) or just 1-3.
**Decision:** Steps 1-3 only (reference extraction, signature building, scoring). Alteromonas N-recycling is a separate analysis.

### Q3: Study selection approach?
**Options:** (A) Keep same Tolonen+Read pair. (B) Start from KG discovery, let selection emerge. (C) Other.
**Decision:** B — data-driven from KG. Part of understanding the method from the ground up.

### Q4: Scoring metrics?
**Options:** (A) One metric, understood thoroughly. (B) Carry forward all three. (C) Decide after seeing data.
**Decision:** A — rank score as sole primary metric. One understood metric > three black boxes.

### Q5: How interactive?
**Options:** (A) Deep — trace 3-5 genes at each step. (B) Moderate — QC + deep-dive on surprises. (C) Light.
**Decision:** A — deep. This is how we break the black box.

### Q6: Where does reusable logic live?
**Options:** (A) Inside analysis dir, productize later. (B) Shared location from the start.
**Decision:** A — analysis-first, per the skill's phase 1 rule.

### Q7: Signature application in utils?
**Decision:** Yes — `apply_signature()` is reusable methodology, not script glue. Returns the inspectable intermediate (which genes matched, concordance, values).

### Q8: Per-step scripts with logging?
**Decision:** Yes. Scripts write diagnostic logs to `logs/`. Explore phase split: script diagnostics (reproducible) + chat reasoning (notebook).

### Q9: Where do logs sit?
**Decision:** `logs/` inside the analysis directory, alongside data and results.

### Q10: Companion Jupyter notebooks?
**Decision:** One optional notebook after step 6 (not per-step). Decided at that point whether it's worth creating.

### Q11: Table scope handling?
**Decision:** Capture `table_scope` from DE response envelope, add as column on every row in CSVs. Critical for interpreting absent genes.

### Q12: Toy tests as scripts?
**Decision:** Yes — saved as `sig_utils/tests/test_*.py` to seed the test suite for productization.

### Q13: Shared extraction utility?
**Decision:** `sig_utils/extraction.py` wraps DE extraction logic. All download scripts use it.

### Q14: Permutation testing scope?
**Decision:** Part of scoring, not limited to specific experiments. `score_with_significance()` wraps rank_score + permutation into one call.

### Q15: Superpowers products in analysis dir?
**Decision:** Yes — spec, plan, and brainstorm log copied into `superpowers/` so decision history is self-contained.

### Q16: Positive/negative controls?
**Decision:** Every discovered experiment gets classified. Positive controls (N-stress) must score high, negative controls (non-N or early timepoints) must score low. Gate: controls must separate before scoring targets.
```

- [ ] **Step 5: Commit scaffold**

```bash
git add "${ANALYSIS_DIR}"
git commit -m "scaffold: analysis directory for N-limitation signature v2"
```

---

### Task 2: Build sig_utils/extraction.py — DE extraction utility

**Files:**
- Create: `sig_utils/extraction.py`
- Create: `sig_utils/io.py`

- [ ] **Step 1: Write extraction.py**

```python
"""Reusable DE data extraction from the multiomics KG.

Wraps differential_expression_by_gene() to produce DataFrames in a
consistent schema with table_scope on every row. All extraction scripts
use this module.
"""

import pandas as pd
from multiomics_explorer import differential_expression_by_gene
from multiomics_explorer.analysis import to_dataframe


def extract_de(
    experiment_ids: list[str],
    organism: str = "MED4",
    verbose: bool = True,
    significant_only: bool = False,
) -> tuple[pd.DataFrame, dict]:
    """Extract DE data for experiments, returning DataFrame + envelope metadata.

    Calls differential_expression_by_gene with verbose=True, limit=None.
    Captures table_scope from the envelope and ensures it's on every row
    (verbose=True already includes it per-row, but we verify).

    Args:
        experiment_ids: List of experiment IDs to extract.
        organism: Organism name (fuzzy matched).
        verbose: Include product, gene_category, etc. Default True.
        significant_only: If True, only return significant rows.

    Returns:
        (df, envelope) where:
        - df: DataFrame with all DE rows, table_scope as a column
        - envelope: dict with metadata (total_matching, experiments, etc.)
    """
    result = differential_expression_by_gene(
        experiment_ids=experiment_ids,
        organism=organism,
        verbose=verbose,
        significant_only=significant_only,
        limit=None,
    )

    assert not result["truncated"], (
        f"Results truncated! total_matching={result['total_matching']}, "
        f"returned={result['returned']}"
    )

    df = to_dataframe(result)

    # Verify table_scope is present (verbose=True puts it on rows)
    if "table_scope" not in df.columns and result.get("experiments"):
        # Fallback: get from envelope and broadcast
        scope_map = {
            exp["experiment_id"]: exp["table_scope"]
            for exp in result["experiments"]
        }
        df["table_scope"] = df["experiment_id"].map(scope_map)

    # Ensure rank columns are nullable int
    for col in ["rank", "rank_up", "rank_down"]:
        if col in df.columns:
            df[col] = pd.to_numeric(df[col], errors="coerce").astype("Int64")

    # Build envelope summary (without the full results)
    envelope = {
        k: v for k, v in result.items()
        if k not in ("results",)
    }

    return df, envelope


def extraction_summary(df: pd.DataFrame, envelope: dict) -> str:
    """Produce a human-readable summary of an extraction for logging.

    Args:
        df: DataFrame from extract_de.
        envelope: Envelope dict from extract_de.

    Returns:
        Multi-line summary string.
    """
    lines = [
        f"Rows: {len(df)}",
        f"Genes: {df['locus_tag'].nunique()}",
        f"Experiments: {envelope.get('experiment_count', '?')}",
    ]

    if "timepoint" in df.columns:
        lines.append(f"Timepoints: {sorted(df['timepoint'].unique())}")

    if "expression_status" in df.columns:
        lines.append(f"Status: {df['expression_status'].value_counts().to_dict()}")

    if "table_scope" in df.columns:
        lines.append(f"Table scope: {df['table_scope'].unique().tolist()}")

    if "omics_type" in df.columns:
        lines.append(f"Omics type: {df['omics_type'].unique().tolist()}")

    return "\n".join(lines)


def check_marker_genes(
    df: pd.DataFrame,
    marker_locus_tags: list[str],
) -> pd.DataFrame:
    """Extract rows for marker genes from a DE DataFrame.

    Returns a filtered DataFrame with only the marker gene rows,
    useful for QC and exploration logging.

    Args:
        df: DE DataFrame.
        marker_locus_tags: List of locus tags to check.

    Returns:
        DataFrame filtered to marker genes, sorted by locus_tag and timepoint.
    """
    marker_df = df[df["locus_tag"].isin(marker_locus_tags)].copy()

    sort_cols = ["locus_tag"]
    if "timepoint_order" in marker_df.columns:
        sort_cols.append("timepoint_order")
    elif "timepoint" in marker_df.columns:
        sort_cols.append("timepoint")

    return marker_df.sort_values(sort_cols).reset_index(drop=True)
```

- [ ] **Step 2: Write io.py**

```python
"""I/O helpers for signature analysis.

Load/save DE data, signatures, and score tables with consistent column types.
"""

import pandas as pd
from pathlib import Path

# Base path for data/ and results/ directories
ANALYSIS_DIR = Path(__file__).parent.parent
DATA_DIR = ANALYSIS_DIR / "data"
RESULTS_DIR = ANALYSIS_DIR / "results"
LOGS_DIR = ANALYSIS_DIR / "logs"


def load_de_csv(path: str | Path) -> pd.DataFrame:
    """Load a DE extract CSV.

    Ensures consistent types for rank columns (nullable int)
    and expression_status (string).

    Args:
        path: Path to CSV file (absolute or relative to DATA_DIR).
    """
    path = Path(path)
    if not path.is_absolute():
        path = DATA_DIR / path

    df = pd.read_csv(path)

    for col in ["rank", "rank_up", "rank_down"]:
        if col in df.columns:
            df[col] = pd.to_numeric(df[col], errors="coerce").astype("Int64")

    return df


def load_signature_csv(path: str | Path) -> pd.DataFrame:
    """Load a signature definition CSV.

    Expected columns: locus_tag, gene_name, direction, signature_type,
    plus rank columns.

    Args:
        path: Path to CSV file (absolute or relative to DATA_DIR).
    """
    path = Path(path)
    if not path.is_absolute():
        path = DATA_DIR / path

    df = pd.read_csv(path)

    for col in df.columns:
        if "rank" in col:
            df[col] = pd.to_numeric(df[col], errors="coerce").astype("Int64")

    return df


def save_csv(df: pd.DataFrame, path: str | Path, subdir: str = "data") -> Path:
    """Save a DataFrame to a CSV in the analysis directory.

    Args:
        df: DataFrame to save.
        path: Filename or path. If relative, placed under subdir.
        subdir: "data" or "results". Ignored if path is absolute.

    Returns:
        The absolute path where the file was saved.
    """
    path = Path(path)
    if not path.is_absolute():
        base = DATA_DIR if subdir == "data" else RESULTS_DIR
        base.mkdir(exist_ok=True)
        path = base / path

    df.to_csv(path, index=False)
    return path


def write_log(content: str, filename: str) -> Path:
    """Write diagnostic log content to logs/.

    Args:
        content: Log text to write.
        filename: Log filename (e.g., '01_discover_experiments.log').

    Returns:
        The absolute path where the log was saved.
    """
    LOGS_DIR.mkdir(exist_ok=True)
    path = LOGS_DIR / filename
    path.write_text(content)
    return path
```

- [ ] **Step 3: Commit**

```bash
git add sig_utils/extraction.py sig_utils/io.py
git commit -m "feat: sig_utils extraction and I/O utilities"
```

---

### Task 3: Script 01 — KG discovery

**Files:**
- Create: `scripts/01_discover_experiments.py`

This step is **interactive**. The script queries the KG and produces a scoping table. The researcher and Claude walk through it together, classify experiments, and decide which are references/controls/targets. The notebook entry is written after the interactive session.

- [ ] **Step 1: Write the discovery script**

```python
"""Step 1: Discover what MED4 N-limitation data exists in the KG.

Queries list_experiments and list_publications for MED4, produces a
scoping table with experiment metadata for classification.

Outputs:
    data/experiment_scoping.csv
    logs/01_discover_experiments.log

Run from multiomics_research root:
    uv run analyses/ANALYSIS_DIR/scripts/01_discover_experiments.py
"""

import sys
from pathlib import Path
from datetime import datetime

ANALYSIS_DIR = Path(__file__).parent.parent
sys.path.insert(0, str(ANALYSIS_DIR))

import pandas as pd
from multiomics_explorer import list_experiments, list_publications, resolve_gene
from multiomics_explorer.analysis import to_dataframe, experiments_to_dataframe
from sig_utils.io import save_csv, write_log


# Marker genes to confirm in the KG
MARKER_GENES = {
    "glnA": "PMM0920",
    "cynA": "PMM0370",
    "rbcL": "PMM0550",
    "atpD": "PMM1452",
    "PMM0030": "PMM0030",
}


def main():
    log_lines = [f"Step 1: KG Discovery — {datetime.now().isoformat()}", ""]

    # 1. List all MED4 experiments
    print("Querying MED4 experiments...")
    exp_result = list_experiments(organism="MED4", verbose=True, limit=None)
    exp_df = to_dataframe(exp_result)
    exp_tp_df = experiments_to_dataframe(exp_result)

    log_lines.append(f"Total experiments: {len(exp_df)}")
    log_lines.append(f"Total experiment x timepoints: {len(exp_tp_df)}")
    log_lines.append("")

    # 2. List publications
    print("Querying MED4 publications...")
    pub_result = list_publications(organism="MED4", limit=None)
    pub_df = to_dataframe(pub_result)

    log_lines.append(f"Publications: {len(pub_df)}")
    for _, row in pub_df.iterrows():
        log_lines.append(f"  {row.get('doi', '?')}: {row.get('title', '?')}")
    log_lines.append("")

    # 3. Build scoping table
    print("Building scoping table...")

    # Select columns for scoping
    scoping_cols = [
        "experiment_id", "experiment_name", "omics_type",
        "treatment_type", "table_scope", "table_scope_detail",
        "is_time_course", "background_factors", "coculture_partner",
    ]
    # Keep only columns that exist
    scoping_cols = [c for c in scoping_cols if c in exp_df.columns]
    scoping_df = exp_df[scoping_cols].copy()

    # Add a classification column (to be filled interactively)
    scoping_df["classification"] = ""
    scoping_df["notes"] = ""

    save_csv(scoping_df, "experiment_scoping.csv")

    log_lines.append("Scoping table columns: " + ", ".join(scoping_df.columns.tolist()))
    log_lines.append(f"Experiments by treatment_type:")
    if "treatment_type" in scoping_df.columns:
        for tt, count in scoping_df["treatment_type"].value_counts().items():
            log_lines.append(f"  {tt}: {count}")
    log_lines.append("")
    if "table_scope" in scoping_df.columns:
        log_lines.append(f"Experiments by table_scope:")
        for ts, count in scoping_df["table_scope"].value_counts().items():
            log_lines.append(f"  {ts}: {count}")
        log_lines.append("")

    # 4. Confirm marker genes
    print("Confirming marker genes...")
    log_lines.append("Marker gene confirmation:")
    for gene_name, locus_tag in MARKER_GENES.items():
        result = resolve_gene(search_term=locus_tag, organism="MED4")
        found = len(result.get("results", [])) > 0
        status = "FOUND" if found else "NOT FOUND"
        log_lines.append(f"  {gene_name} ({locus_tag}): {status}")
        if found:
            gene = result["results"][0]
            log_lines.append(f"    product: {gene.get('product', '?')}")
    log_lines.append("")

    # 5. Print summary for interactive exploration
    print("\n=== SCOPING TABLE ===")
    print(scoping_df.to_string(index=False))
    print(f"\nTotal: {len(scoping_df)} experiments")
    print("\nClassification needed: positive_control / negative_control / target / irrelevant")
    print("Discuss with researcher to fill the 'classification' column.")

    # Write log
    log_lines.append("--- Full scoping table ---")
    log_lines.append(scoping_df.to_string(index=False))
    write_log("\n".join(log_lines), "01_discover_experiments.log")

    print(f"\nOutputs: data/experiment_scoping.csv, logs/01_discover_experiments.log")


if __name__ == "__main__":
    main()
```

- [ ] **Step 2: Run the script**

```bash
uv run analyses/ANALYSIS_DIR/scripts/01_discover_experiments.py
```

(Replace `ANALYSIS_DIR` with the actual directory name.)

- [ ] **Step 3: Interactive exploration (RESEARCHER CHECKPOINT)**

This is the **explore→decide** phase. In chat:
- Walk through the scoping table together
- Classify each experiment: positive_control / negative_control / target / irrelevant
- Discuss which experiments to use as references for signature building
- Confirm marker gene locus tags
- Add any edge-case tracer genes identified from the experiment list

After the interactive session, update `data/experiment_scoping.csv` with classifications and write the notebook entry.

- [ ] **Step 4: Write notebook entry**

Append to `exploration/2026-04-08-notebook.md`:

```markdown
---

## YYYY-MM-DD HH:MM — Step 1: KG discovery

### Command
```bash
uv run scripts/01_discover_experiments.py
```

### Outputs
- [experiment_scoping.csv](../data/experiment_scoping.csv) — N experiments, classifications
- [01_discover_experiments.log](../logs/01_discover_experiments.log) — full query results

### QC
- Total experiments found: N
- Treatment types: ...
- Table scopes: ...
- All 5 marker genes confirmed in KG

### Exploration
- [Filled during interactive session — experiment classifications, rationale]

### Decision
- References: [list]
- Positive controls: [list]
- Negative controls: [list]
- Targets: [list]
```

- [ ] **Step 5: Update DATA_MANIFEST.md and commit**

```bash
git add data/experiment_scoping.csv logs/01_discover_experiments.log exploration/ data/DATA_MANIFEST.md
git commit -m "feat: Step 1 — KG discovery and experiment classification"
```

---

### Task 4: Script 02 — Reference and control DE extraction

**Files:**
- Create: `scripts/02_extract_reference_de.py`

- [ ] **Step 1: Write the extraction script**

```python
"""Step 2: Extract DE data for reference and control experiments.

Uses sig_utils/extraction.py for consistent schema.
Experiment IDs come from the classified scoping table (Step 1).

Outputs:
    data/de_*.csv — one per experiment
    logs/02_extract_reference_de.log

Run from multiomics_research root:
    uv run analyses/ANALYSIS_DIR/scripts/02_extract_reference_de.py [--explore]
"""

import sys
import argparse
from pathlib import Path
from datetime import datetime

ANALYSIS_DIR = Path(__file__).parent.parent
sys.path.insert(0, str(ANALYSIS_DIR))

import pandas as pd
from sig_utils.extraction import extract_de, extraction_summary, check_marker_genes
from sig_utils.io import save_csv, write_log, load_de_csv, DATA_DIR

# Marker genes for exploration
MARKERS = ["PMM0920", "PMM0370", "PMM0550", "PMM1452", "PMM0030"]

# Experiment registry — filled after Step 1 classification
# Format: (experiment_id, output_filename, label, role)
EXPERIMENTS = [
    # === REFERENCES (for signature building) ===
    # Populated after Step 1 interactive session, e.g.:
    # ("10.1038/msb4100087_nitrogen_...", "de_ref_tolonen_ndep.csv", "Tolonen 2006 N-dep", "reference"),
    # ("10.1038/ismej.2017.88_nitrogen_...", "de_ref_read_ndep.csv", "Read 2017 N-dep", "reference"),

    # === POSITIVE CONTROLS ===
    # Other N-stress experiments not used for signature building

    # === NEGATIVE CONTROLS ===
    # Non-N stresses, alternative N sources, early timepoints
    # ("10.1038/msb4100087_growth_medium_cyanate_...", "de_ctrl_tolonen_cyanate.csv", "Tolonen cyanate", "negative_control"),
    # ("10.1038/msb4100087_growth_medium_urea_...", "de_ctrl_tolonen_urea.csv", "Tolonen urea", "negative_control"),
]


def main():
    parser = argparse.ArgumentParser(description="Extract reference and control DE data")
    parser.add_argument("--explore", action="store_true", help="Print marker gene traces")
    args = parser.parse_args()

    log_lines = [f"Step 2: Reference DE Extraction — {datetime.now().isoformat()}", ""]

    for exp_id, filename, label, role in EXPERIMENTS:
        print(f"\nExtracting {label} ({role})...")
        df, envelope = extract_de(experiment_ids=[exp_id])

        summary = extraction_summary(df, envelope)
        print(summary)
        log_lines.append(f"=== {label} ({role}) ===")
        log_lines.append(summary)

        # Check for metadata issues
        if "timepoint" in df.columns:
            single_tp = df[df["timepoint"].isin(["single", None, ""])]
            if len(single_tp) > 0:
                warning = f"  WARNING: {len(single_tp)} rows with timepoint='single' or null"
                print(warning)
                log_lines.append(warning)

        # Save
        save_csv(df, filename)
        log_lines.append(f"Saved: data/{filename}")

        # Marker gene exploration
        if args.explore:
            marker_df = check_marker_genes(df, MARKERS)
            if len(marker_df) > 0:
                print(f"\n  Marker genes in {label}:")
                cols = ["locus_tag", "gene_name", "timepoint", "log2fc",
                        "expression_status", "rank_up", "rank_down"]
                cols = [c for c in cols if c in marker_df.columns]
                print(marker_df[cols].to_string(index=False))

            log_lines.append(f"\nMarker genes:")
            log_lines.append(marker_df[cols].to_string(index=False) if len(marker_df) > 0 else "  No marker genes found")

        log_lines.append("")

    write_log("\n".join(log_lines), "02_extract_reference_de.log")
    print(f"\nDone. See logs/02_extract_reference_de.log")


if __name__ == "__main__":
    main()
```

- [ ] **Step 2: Fill in EXPERIMENTS list after Step 1 classification**

Replace the placeholder comments in the `EXPERIMENTS` list with actual experiment IDs and filenames from the Step 1 scoping table. This happens during the interactive session.

- [ ] **Step 3: Run the script**

```bash
uv run analyses/ANALYSIS_DIR/scripts/02_extract_reference_de.py --explore
```

- [ ] **Step 4: Interactive exploration (RESEARCHER CHECKPOINT)**

In chat:
- Review marker genes in each reference dataset
- Check: are glnA and cynA significant_up? rbcL and atpD significant_down? What ranks?
- Check table_scope per experiment
- Flag any metadata issues (single timepoints, missing genes)
- Identify edge-case tracer genes (e.g., gene significant in one study but absent in another)

- [ ] **Step 5: Write notebook entry and update manifests**

Append notebook entry (same template as Task 3 step 4). Update `data/DATA_MANIFEST.md` with extracted files.

- [ ] **Step 6: Commit**

```bash
git add scripts/02_extract_reference_de.py data/ logs/ exploration/ 
git commit -m "feat: Step 2 — reference and control DE extraction"
```

---

### Task 5: Build sig_utils/signature.py — signature construction

**Files:**
- Create: `sig_utils/signature.py`
- Create: `sig_utils/tests/test_signature.py`

This task has a design→implement→toy-test→apply sub-cycle.

- [ ] **Step 1: Design the signature building logic (RESEARCHER CHECKPOINT)**

Before writing code, present the logic with worked examples in chat:

1. **`summarize_per_gene(de_df)`** — for a gene significant_up at 3 timepoints and significant_down at 1, majority direction = up. Best directional rank = min(rank_up across up-timepoints). Show with concrete numbers.

2. **`intersect_references(study_a_summary, study_b_summary)`** — gene X is up in both → core. Gene Y is up in A, down in B → discordant, excluded. Gene Z is up in A, absent from B → extended. Show the merge logic.

3. **`classify_non_core()`** — explain the categories: `study_a_only_study_b_ns` (present but not significant in B), `study_a_only_study_b_absent` (not in B's dataset at all), `study_b_only`.

Log the worked examples in the notebook.

- [ ] **Step 2: Write test_signature.py with toy data**

```python
"""Toy data tests for signature construction.

Hand-calculated expected values for every operation.
These serve as the seed for a real test suite during productization.
"""

import pandas as pd
import pytest
import sys
from pathlib import Path

# Add sig_utils to path
sys.path.insert(0, str(Path(__file__).parent.parent.parent))

from sig_utils.signature import summarize_per_gene, intersect_references, classify_non_core


# === Toy DE data ===

def make_toy_study_a():
    """Study A: 5 genes, 2 timepoints."""
    return pd.DataFrame([
        # Gene G1: up at both timepoints, rank_up 2 at tp1, rank_up 1 at tp2
        {"locus_tag": "G1", "gene_name": "geneA", "timepoint": "tp1",
         "log2fc": 2.0, "expression_status": "significant_up",
         "rank": 3, "rank_up": 2, "rank_down": pd.NA},
        {"locus_tag": "G1", "gene_name": "geneA", "timepoint": "tp2",
         "log2fc": 3.0, "expression_status": "significant_up",
         "rank": 1, "rank_up": 1, "rank_down": pd.NA},
        # Gene G2: down at tp1, not significant at tp2
        {"locus_tag": "G2", "gene_name": "geneB", "timepoint": "tp1",
         "log2fc": -1.5, "expression_status": "significant_down",
         "rank": 5, "rank_up": pd.NA, "rank_down": 3},
        {"locus_tag": "G2", "gene_name": "geneB", "timepoint": "tp2",
         "log2fc": -0.3, "expression_status": "not_significant",
         "rank": 10, "rank_up": pd.NA, "rank_down": pd.NA},
        # Gene G3: up at tp1 only (only in study A)
        {"locus_tag": "G3", "gene_name": "geneC", "timepoint": "tp1",
         "log2fc": 1.0, "expression_status": "significant_up",
         "rank": 7, "rank_up": 4, "rank_down": pd.NA},
        {"locus_tag": "G3", "gene_name": "geneC", "timepoint": "tp2",
         "log2fc": 0.1, "expression_status": "not_significant",
         "rank": 20, "rank_up": pd.NA, "rank_down": pd.NA},
        # Gene G4: mixed direction — up at tp1, down at tp2 (majority = tie, resolved arbitrarily)
        {"locus_tag": "G4", "gene_name": "geneD", "timepoint": "tp1",
         "log2fc": 1.2, "expression_status": "significant_up",
         "rank": 6, "rank_up": 3, "rank_down": pd.NA},
        {"locus_tag": "G4", "gene_name": "geneD", "timepoint": "tp2",
         "log2fc": -1.8, "expression_status": "significant_down",
         "rank": 4, "rank_up": pd.NA, "rank_down": 2},
        # Gene G5: not significant at either timepoint
        {"locus_tag": "G5", "gene_name": "geneE", "timepoint": "tp1",
         "log2fc": 0.05, "expression_status": "not_significant",
         "rank": 50, "rank_up": pd.NA, "rank_down": pd.NA},
        {"locus_tag": "G5", "gene_name": "geneE", "timepoint": "tp2",
         "log2fc": -0.1, "expression_status": "not_significant",
         "rank": 45, "rank_up": pd.NA, "rank_down": pd.NA},
    ])


def make_toy_study_b():
    """Study B: 4 genes, 1 timepoint. G3 absent from dataset."""
    return pd.DataFrame([
        # Gene G1: up, rank_up 1 — concordant with study A
        {"locus_tag": "G1", "gene_name": "geneA", "timepoint": "tp1",
         "log2fc": 4.0, "expression_status": "significant_up",
         "rank": 1, "rank_up": 1, "rank_down": pd.NA},
        # Gene G2: down, rank_down 2 — concordant with study A
        {"locus_tag": "G2", "gene_name": "geneB", "timepoint": "tp1",
         "log2fc": -2.0, "expression_status": "significant_down",
         "rank": 3, "rank_up": pd.NA, "rank_down": 2},
        # Gene G4: up — DISCORDANT if A says down (depends on A's majority)
        {"locus_tag": "G4", "gene_name": "geneD", "timepoint": "tp1",
         "log2fc": 1.5, "expression_status": "significant_up",
         "rank": 4, "rank_up": 2, "rank_down": pd.NA},
        # Gene G5: not significant (present but not DE)
        {"locus_tag": "G5", "gene_name": "geneE", "timepoint": "tp1",
         "log2fc": 0.1, "expression_status": "not_significant",
         "rank": 20, "rank_up": pd.NA, "rank_down": pd.NA},
        # Gene G6: down — only in study B
        {"locus_tag": "G6", "gene_name": "geneF", "timepoint": "tp1",
         "log2fc": -3.0, "expression_status": "significant_down",
         "rank": 2, "rank_up": pd.NA, "rank_down": 1},
    ])


class TestSummarizePerGene:
    def test_majority_direction(self):
        """G1 is up at both timepoints → direction = up."""
        df = make_toy_study_a()
        result = summarize_per_gene(df)
        g1 = result[result["locus_tag"] == "G1"].iloc[0]
        assert g1["direction"] == "up"

    def test_best_directional_rank(self):
        """G1 has rank_up=2 at tp1, rank_up=1 at tp2 → best_dir_rank = 1."""
        df = make_toy_study_a()
        result = summarize_per_gene(df)
        g1 = result[result["locus_tag"] == "G1"].iloc[0]
        assert g1["best_dir_rank"] == 1

    def test_single_timepoint_significant(self):
        """G2 is significant_down at tp1 only → direction = down."""
        df = make_toy_study_a()
        result = summarize_per_gene(df)
        g2 = result[result["locus_tag"] == "G2"].iloc[0]
        assert g2["direction"] == "down"

    def test_not_significant_excluded(self):
        """G5 is not significant at any timepoint → excluded from summary."""
        df = make_toy_study_a()
        result = summarize_per_gene(df)
        assert "G5" not in result["locus_tag"].values

    def test_gene_count(self):
        """4 genes significant in at least one timepoint (G1, G2, G3, G4)."""
        df = make_toy_study_a()
        result = summarize_per_gene(df)
        assert len(result) == 4


class TestIntersectReferences:
    def test_concordant_core(self):
        """G1 (up in both) and G2 (down in both) → core."""
        a_summary = summarize_per_gene(make_toy_study_a())
        b_summary = summarize_per_gene(make_toy_study_b())
        core, extended = intersect_references(a_summary, b_summary)
        core_tags = core["locus_tag"].tolist()
        assert "G1" in core_tags
        assert "G2" in core_tags

    def test_core_direction(self):
        """Core genes have correct direction."""
        a_summary = summarize_per_gene(make_toy_study_a())
        b_summary = summarize_per_gene(make_toy_study_b())
        core, _ = intersect_references(a_summary, b_summary)
        g1 = core[core["locus_tag"] == "G1"].iloc[0]
        assert g1["direction"] == "up"
        g2 = core[core["locus_tag"] == "G2"].iloc[0]
        assert g2["direction"] == "down"

    def test_b_only_extended(self):
        """G6 is only in study B → extended, study_b_only."""
        a_summary = summarize_per_gene(make_toy_study_a())
        b_summary = summarize_per_gene(make_toy_study_b())
        _, extended = intersect_references(a_summary, b_summary)
        g6 = extended[extended["locus_tag"] == "G6"]
        assert len(g6) == 1
        assert "study_b_only" in g6.iloc[0]["signature_type"]


class TestClassifyNonCore:
    def test_a_only_b_absent(self):
        """G3 is in study A but absent from study B dataset → a_only_b_absent."""
        a_summary = summarize_per_gene(make_toy_study_a())
        b_summary = summarize_per_gene(make_toy_study_b())
        study_b_all_tags = {"G1", "G2", "G4", "G5", "G6"}  # G3 absent
        _, extended = intersect_references(
            a_summary, b_summary,
            study_b_all_locus_tags=study_b_all_tags,
        )
        g3 = extended[extended["locus_tag"] == "G3"]
        assert len(g3) == 1
        assert "absent" in g3.iloc[0]["signature_type"]
```

- [ ] **Step 3: Run tests to verify they fail**

```bash
cd analyses/ANALYSIS_DIR && uv run python -m pytest sig_utils/tests/test_signature.py -v
```

Expected: FAIL — `sig_utils.signature` module doesn't exist yet.

- [ ] **Step 4: Write signature.py**

```python
"""Signature building: summarize DE per gene, intersect references, classify.

All functions take DataFrames and return DataFrames. No KG dependency.
"""

import pandas as pd


def summarize_per_gene(de_df: pd.DataFrame) -> pd.DataFrame:
    """Summarize a multi-timepoint DE DataFrame to one row per gene.

    For each gene:
    1. Filter to significant rows only
    2. Determine majority direction across timepoints
    3. Pick the timepoint with best directional rank in majority direction

    Args:
        de_df: DataFrame with columns: locus_tag, gene_name, timepoint,
            log2fc, expression_status, rank, rank_up, rank_down.
            May also include product, gene_category (preserved if present).

    Returns:
        DataFrame with one row per gene: locus_tag, gene_name, direction,
        peak_timepoint, best_dir_rank, best_global_rank.
        Plus product, gene_category if present in input.
    """
    result_cols = [
        "locus_tag", "gene_name", "direction",
        "peak_timepoint", "best_dir_rank", "best_global_rank",
    ]
    # Preserve metadata columns if present
    meta_cols = [c for c in ["product", "gene_category"] if c in de_df.columns]

    sig = de_df[de_df["expression_status"].isin(["significant_up", "significant_down"])].copy()
    if sig.empty:
        return pd.DataFrame(columns=result_cols + meta_cols)

    sig["direction"] = sig["expression_status"].map({
        "significant_up": "up", "significant_down": "down"
    })
    sig["dir_rank"] = sig.apply(
        lambda r: r["rank_up"] if r["direction"] == "up" else r["rank_down"], axis=1
    )

    # Majority direction per gene
    direction_counts = (
        sig.groupby("locus_tag")["direction"]
        .value_counts()
        .unstack(fill_value=0)
    )
    majority_direction = direction_counts.idxmax(axis=1).rename("majority_direction")

    # Filter to majority direction rows
    sig = sig.merge(majority_direction, on="locus_tag")
    sig_majority = sig[sig["direction"] == sig["majority_direction"]]

    # Peak timepoint: best directional rank in majority direction
    peak = sig_majority.loc[sig_majority.groupby("locus_tag")["dir_rank"].idxmin()]

    # Best ranks
    best_dir = sig_majority.groupby("locus_tag")["dir_rank"].min().rename("best_dir_rank")
    best_global = sig.groupby("locus_tag")["rank"].min().rename("best_global_rank")

    result = peak[["locus_tag", "gene_name", "direction", "timepoint"] + meta_cols].copy()
    result = result.rename(columns={"timepoint": "peak_timepoint"})
    result["direction"] = result["locus_tag"].map(majority_direction)
    result = result.merge(best_dir, on="locus_tag")
    result = result.merge(best_global, on="locus_tag")

    return result.reset_index(drop=True)


def intersect_references(
    study_a: pd.DataFrame,
    study_b: pd.DataFrame,
    study_a_name: str = "study_a",
    study_b_name: str = "study_b",
    study_b_all_locus_tags: set[str] | None = None,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Intersect two per-gene DE summaries into core and extended signatures.

    Core: genes significant in both studies with concordant direction.
    Extended: genes significant in only one study.
    Discordant genes (different direction in each study) are excluded from both.

    Args:
        study_a: Output of summarize_per_gene for study A.
        study_b: Output of summarize_per_gene for study B.
        study_a_name: Label for study A (used in column prefixes).
        study_b_name: Label for study B (used in column prefixes).
        study_b_all_locus_tags: All locus tags in study B's dataset
            (not just significant). Used to distinguish "absent from
            dataset" vs "present but not significant" for study-A-only genes.

    Returns:
        (core_df, extended_df) — each sorted by cross_study_best_dir_rank.
    """
    pa, pb = study_a_name, study_b_name

    # Prefix study-specific columns
    a = study_a.rename(columns={
        "direction": "direction_a",
        "gene_name": "gene_name_a",
        "peak_timepoint": f"{pa}_peak_timepoint",
        "best_dir_rank": f"{pa}_best_dir_rank",
        "best_global_rank": f"{pa}_best_global_rank",
    })
    b = study_b.rename(columns={
        "direction": "direction_b",
        "gene_name": "gene_name_b",
        "peak_timepoint": f"{pb}_peak_timepoint",
        "best_dir_rank": f"{pb}_best_dir_rank",
        "best_global_rank": f"{pb}_best_global_rank",
    })

    # Also rename metadata columns if present
    for meta in ["product", "gene_category"]:
        if meta in a.columns:
            a = a.rename(columns={meta: f"{meta}_a"})
        if meta in b.columns:
            b = b.rename(columns={meta: f"{meta}_b"})

    merged = a.merge(b, on="locus_tag", how="outer")

    # Resolve gene_name: prefer A, fall back to B
    for meta in ["gene_name", "product", "gene_category"]:
        col_a, col_b = f"{meta}_a", f"{meta}_b"
        if col_a in merged.columns and col_b in merged.columns:
            merged[meta] = merged[col_a].fillna(merged[col_b])
            merged = merged.drop(columns=[col_a, col_b])

    # CORE: both present, concordant direction
    both = merged.dropna(subset=["direction_a", "direction_b"])
    concordant = both[both["direction_a"] == both["direction_b"]].copy()
    concordant["direction"] = concordant["direction_a"]
    concordant["signature_type"] = "core"
    concordant["cross_study_best_dir_rank"] = concordant[
        [f"{pa}_best_dir_rank", f"{pb}_best_dir_rank"]
    ].min(axis=1)

    # EXTENDED: in one study only (discordant excluded)
    a_only_mask = merged["direction_b"].isna() & merged["direction_a"].notna()
    b_only_mask = merged["direction_a"].isna() & merged["direction_b"].notna()

    a_only = merged[a_only_mask].copy()
    a_only["direction"] = a_only["direction_a"]
    if study_b_all_locus_tags is not None:
        a_only["signature_type"] = a_only["locus_tag"].apply(
            lambda lt: f"{pa}_only_{pb}_ns" if lt in study_b_all_locus_tags
            else f"{pa}_only_{pb}_absent"
        )
    else:
        a_only["signature_type"] = f"{pa}_only"

    b_only = merged[b_only_mask].copy()
    b_only["direction"] = b_only["direction_b"]
    b_only["signature_type"] = f"{pb}_only"

    extended = pd.concat([a_only, b_only], ignore_index=True)
    extended["cross_study_best_dir_rank"] = extended.apply(
        lambda r: r[f"{pa}_best_dir_rank"] if pd.notna(r.get(f"{pa}_best_dir_rank"))
        else r.get(f"{pb}_best_dir_rank"),
        axis=1,
    )

    # Assemble output columns
    cols = [
        "locus_tag", "gene_name", "direction", "signature_type",
        f"{pa}_peak_timepoint", f"{pa}_best_dir_rank", f"{pa}_best_global_rank",
        f"{pb}_peak_timepoint", f"{pb}_best_dir_rank", f"{pb}_best_global_rank",
        "cross_study_best_dir_rank",
    ]
    # Add metadata columns if present
    for meta in ["product", "gene_category"]:
        if meta in concordant.columns:
            cols.append(meta)

    for df in [concordant, extended]:
        for c in cols:
            if c not in df.columns:
                df[c] = pd.NA

    core_df = concordant[cols].sort_values("cross_study_best_dir_rank").reset_index(drop=True)
    extended_df = extended[cols].sort_values("cross_study_best_dir_rank").reset_index(drop=True)

    return core_df, extended_df
```

- [ ] **Step 5: Run tests to verify they pass**

```bash
cd analyses/ANALYSIS_DIR && uv run python -m pytest sig_utils/tests/test_signature.py -v
```

Expected: All tests PASS.

- [ ] **Step 6: Commit**

```bash
git add sig_utils/signature.py sig_utils/tests/test_signature.py
git commit -m "feat: sig_utils signature construction with toy tests"
```

---

### Task 6: Script 03 — Build signature from real data

**Files:**
- Create: `scripts/03_build_signature.py`

- [ ] **Step 1: Write the build script**

```python
"""Step 3: Build N-limitation gene signature from reference studies.

Uses sig_utils/signature.py functions on real DE data from Step 2.

Outputs:
    data/core_signature.csv
    data/extended_signature.csv
    logs/03_build_signature.log

Run from multiomics_research root:
    uv run analyses/ANALYSIS_DIR/scripts/03_build_signature.py [--explore]
"""

import sys
import argparse
from pathlib import Path
from datetime import datetime

ANALYSIS_DIR = Path(__file__).parent.parent
sys.path.insert(0, str(ANALYSIS_DIR))

import pandas as pd
from sig_utils.signature import summarize_per_gene, intersect_references
from sig_utils.extraction import check_marker_genes
from sig_utils.io import load_de_csv, save_csv, write_log

MARKERS = ["PMM0920", "PMM0370", "PMM0550", "PMM1452", "PMM0030"]

# --- FILL AFTER STEP 1/2 ---
STUDY_A_FILE = "de_ref_tolonen_ndep.csv"  # placeholder
STUDY_B_FILE = "de_ref_read_ndep.csv"     # placeholder
STUDY_A_NAME = "tolonen"
STUDY_B_NAME = "read"
# Timepoints to exclude from study A before summarization
STUDY_A_EXCLUDE_TIMEPOINTS = []  # e.g., ["0h", "3h"] — decided in Step 2
# ---


def main():
    parser = argparse.ArgumentParser(description="Build N-limitation signature")
    parser.add_argument("--explore", action="store_true", help="Print marker gene traces")
    args = parser.parse_args()

    log_lines = [f"Step 3: Build Signature — {datetime.now().isoformat()}", ""]

    # Load reference DE data
    print("Loading reference DE data...")
    de_a = load_de_csv(STUDY_A_FILE)
    de_b = load_de_csv(STUDY_B_FILE)

    log_lines.append(f"Study A ({STUDY_A_NAME}): {len(de_a)} rows, {de_a['locus_tag'].nunique()} genes")
    log_lines.append(f"Study B ({STUDY_B_NAME}): {len(de_b)} rows, {de_b['locus_tag'].nunique()} genes")

    # Exclude timepoints if configured
    if STUDY_A_EXCLUDE_TIMEPOINTS:
        before = len(de_a)
        de_a = de_a[~de_a["timepoint"].isin(STUDY_A_EXCLUDE_TIMEPOINTS)]
        log_lines.append(f"Excluded {STUDY_A_NAME} timepoints {STUDY_A_EXCLUDE_TIMEPOINTS}: {before} → {len(de_a)} rows")

    # Filter funnel
    print("\n=== Filter funnel ===")

    # Summarize per gene
    summary_a = summarize_per_gene(de_a)
    summary_b = summarize_per_gene(de_b)
    log_lines.append(f"\n--- Filter funnel ---")
    log_lines.append(f"{STUDY_A_NAME} significant genes: {len(summary_a)}")
    log_lines.append(f"  up: {(summary_a['direction'] == 'up').sum()}, down: {(summary_a['direction'] == 'down').sum()}")
    log_lines.append(f"{STUDY_B_NAME} significant genes: {len(summary_b)}")
    log_lines.append(f"  up: {(summary_b['direction'] == 'up').sum()}, down: {(summary_b['direction'] == 'down').sum()}")

    print(f"{STUDY_A_NAME}: {len(summary_a)} significant genes (up={( summary_a['direction'] == 'up').sum()}, down={(summary_a['direction'] == 'down').sum()})")
    print(f"{STUDY_B_NAME}: {len(summary_b)} significant genes (up={(summary_b['direction'] == 'up').sum()}, down={(summary_b['direction'] == 'down').sum()})")

    # All locus tags in study B (for classifying A-only genes)
    study_b_all_tags = set(de_b["locus_tag"].unique())

    # Intersect
    core, extended = intersect_references(
        summary_a, summary_b,
        study_a_name=STUDY_A_NAME,
        study_b_name=STUDY_B_NAME,
        study_b_all_locus_tags=study_b_all_tags,
    )

    print(f"\nCore signature: {len(core)} genes (up={(core['direction'] == 'up').sum()}, down={(core['direction'] == 'down').sum()})")
    print(f"Extended: {len(extended)} genes")
    if "signature_type" in extended.columns:
        print(f"  By type: {extended['signature_type'].value_counts().to_dict()}")

    log_lines.append(f"\nCore: {len(core)} genes (up={(core['direction'] == 'up').sum()}, down={(core['direction'] == 'down').sum()})")
    log_lines.append(f"Extended: {len(extended)} genes")
    if "signature_type" in extended.columns:
        log_lines.append(f"  {extended['signature_type'].value_counts().to_dict()}")

    # Save
    save_csv(core, "core_signature.csv")
    save_csv(extended, "extended_signature.csv")

    # Marker gene traces
    if args.explore:
        print("\n=== Marker gene traces ===")
        for tag in MARKERS:
            in_core = core[core["locus_tag"] == tag]
            in_ext = extended[extended["locus_tag"] == tag]
            if len(in_core) > 0:
                row = in_core.iloc[0]
                print(f"  {tag} ({row.get('gene_name', '?')}): CORE, direction={row['direction']}, "
                      f"cross_rank={row['cross_study_best_dir_rank']}")
                log_lines.append(f"  {tag}: CORE, dir={row['direction']}, rank={row['cross_study_best_dir_rank']}")
            elif len(in_ext) > 0:
                row = in_ext.iloc[0]
                print(f"  {tag} ({row.get('gene_name', '?')}): EXTENDED ({row['signature_type']}), "
                      f"direction={row['direction']}")
                log_lines.append(f"  {tag}: EXTENDED ({row['signature_type']})")
            else:
                # Check if in either study summary
                in_a = summary_a[summary_a["locus_tag"] == tag]
                in_b = summary_b[summary_b["locus_tag"] == tag]
                if len(in_a) > 0 and len(in_b) > 0:
                    print(f"  {tag}: DISCORDANT — {STUDY_A_NAME}={in_a.iloc[0]['direction']}, "
                          f"{STUDY_B_NAME}={in_b.iloc[0]['direction']}")
                    log_lines.append(f"  {tag}: DISCORDANT")
                else:
                    print(f"  {tag}: NOT IN SIGNATURE (A={len(in_a)>0}, B={len(in_b)>0})")
                    log_lines.append(f"  {tag}: NOT IN SIGNATURE")

    write_log("\n".join(log_lines), "03_build_signature.log")
    print(f"\nOutputs: data/core_signature.csv, data/extended_signature.csv")


if __name__ == "__main__":
    main()
```

- [ ] **Step 2: Run the script**

```bash
uv run analyses/ANALYSIS_DIR/scripts/03_build_signature.py --explore
```

- [ ] **Step 3: Interactive exploration (RESEARCHER CHECKPOINT)**

In chat:
- Review filter funnel: do the numbers make biological sense?
- Trace all 5 marker genes: are they in core? correct direction? reasonable rank?
- Identify edge-case genes: discordant genes, genes at the core/extended boundary
- Check the top 10 core genes by rank: are they recognizable N-limitation biology?
- Compare with v1 (198 core genes) — similar count? What's different?

- [ ] **Step 4: Write notebook entry, update manifests, commit**

```bash
git add scripts/03_build_signature.py data/core_signature.csv data/extended_signature.csv logs/ exploration/ data/DATA_MANIFEST.md
git commit -m "feat: Step 3 — build N-limitation signature from reference data"
```

---

### Task 7: Script 04 — Target DE extraction

**Files:**
- Create: `scripts/04_extract_target_de.py`

- [ ] **Step 1: Write the target extraction script**

```python
"""Step 4: Extract Weissberg 2025 MED4 DE data for scoring.

Uses the same sig_utils/extraction.py utility as Step 2.

Outputs:
    data/de_weissberg_*.csv — one per experiment/platform/condition
    logs/04_extract_target_de.log

Run from multiomics_research root:
    uv run analyses/ANALYSIS_DIR/scripts/04_extract_target_de.py [--explore]
"""

import sys
import argparse
from pathlib import Path
from datetime import datetime

ANALYSIS_DIR = Path(__file__).parent.parent
sys.path.insert(0, str(ANALYSIS_DIR))

import pandas as pd
from sig_utils.extraction import extract_de, extraction_summary, check_marker_genes
from sig_utils.io import save_csv, write_log, load_signature_csv

MARKERS = ["PMM0920", "PMM0370", "PMM0550", "PMM1452", "PMM0030"]

# --- FILL AFTER STEP 1 ---
# Format: (experiment_id, output_filename, label)
WEISSBERG_EXPERIMENTS = [
    # e.g.:
    # ("exp_id_rnaseq_axenic", "de_weissberg_rnaseq_axenic.csv", "Weissberg RNA-seq axenic"),
    # ("exp_id_rnaseq_coculture", "de_weissberg_rnaseq_coculture.csv", "Weissberg RNA-seq coculture"),
    # ("exp_id_proteomics_axenic", "de_weissberg_proteomics_axenic.csv", "Weissberg proteomics axenic"),
    # ("exp_id_proteomics_coculture", "de_weissberg_proteomics_coculture.csv", "Weissberg proteomics coculture"),
]


def main():
    parser = argparse.ArgumentParser(description="Extract Weissberg 2025 DE data")
    parser.add_argument("--explore", action="store_true", help="Print marker gene traces and coverage")
    args = parser.parse_args()

    log_lines = [f"Step 4: Target DE Extraction — {datetime.now().isoformat()}", ""]

    # Load signature for coverage check
    try:
        core_sig = load_signature_csv("core_signature.csv")
        sig_tags = set(core_sig["locus_tag"].unique())
        log_lines.append(f"Core signature loaded: {len(sig_tags)} genes")
    except FileNotFoundError:
        sig_tags = set()
        log_lines.append("WARNING: core_signature.csv not found, skipping coverage check")

    for exp_id, filename, label in WEISSBERG_EXPERIMENTS:
        print(f"\nExtracting {label}...")
        df, envelope = extract_de(experiment_ids=[exp_id])

        summary = extraction_summary(df, envelope)
        print(summary)
        log_lines.append(f"=== {label} ===")
        log_lines.append(summary)

        # Signature coverage
        if sig_tags:
            exp_tags = set(df["locus_tag"].unique())
            overlap = sig_tags & exp_tags
            missing = sig_tags - exp_tags
            log_lines.append(f"Signature coverage: {len(overlap)}/{len(sig_tags)} genes detected")
            if missing:
                log_lines.append(f"Missing signature genes: {len(missing)}")
            print(f"  Signature coverage: {len(overlap)}/{len(sig_tags)} ({len(missing)} missing)")

        # Check for single timepoint issues
        if "timepoint" in df.columns:
            tps = df["timepoint"].unique()
            if "single" in tps or any(pd.isna(t) for t in tps):
                warning = f"  WARNING: contains single/null timepoints: {tps}"
                print(warning)
                log_lines.append(warning)

        save_csv(df, filename)
        log_lines.append(f"Saved: data/{filename}")

        if args.explore:
            marker_df = check_marker_genes(df, MARKERS)
            if len(marker_df) > 0:
                print(f"\n  Marker genes in {label}:")
                cols = ["locus_tag", "gene_name", "timepoint", "log2fc",
                        "expression_status", "rank_up", "rank_down"]
                cols = [c for c in cols if c in marker_df.columns]
                print(marker_df[cols].to_string(index=False))
                log_lines.append(f"\nMarker genes:")
                log_lines.append(marker_df[cols].to_string(index=False))

        log_lines.append("")

    write_log("\n".join(log_lines), "04_extract_target_de.log")
    print(f"\nDone. See logs/04_extract_target_de.log")


if __name__ == "__main__":
    main()
```

- [ ] **Step 2: Fill WEISSBERG_EXPERIMENTS after Step 1**

- [ ] **Step 3: Run the script**

```bash
uv run analyses/ANALYSIS_DIR/scripts/04_extract_target_de.py --explore
```

- [ ] **Step 4: Interactive exploration (RESEARCHER CHECKPOINT)**

In chat:
- Check marker genes in axenic vs coculture: does glnA look different?
- Check platform coverage: how many signature genes in RNA-seq vs proteomics?
- Flag any timepoint issues
- Compare gene counts with v1 (RNA-seq 1849, proteomics 1424)

- [ ] **Step 5: Write notebook entry, update manifests, commit**

```bash
git add scripts/04_extract_target_de.py data/de_weissberg_*.csv logs/ exploration/ data/DATA_MANIFEST.md
git commit -m "feat: Step 4 — extract Weissberg 2025 target DE data"
```

---

### Task 8: Build sig_utils/scoring.py — metric design and verification

**Files:**
- Create: `sig_utils/scoring.py`
- Create: `sig_utils/tests/test_scoring.py`

- [ ] **Step 1: Design the rank score metric on paper (RESEARCHER CHECKPOINT)**

Present in chat with worked example:

**Rank score formula:**
Given a signature S (genes with expected directions) and a target DE dataset D:

For each signature gene i present in D:
- `concordance_i` = +1 if DE direction matches signature direction, -1 if opposite, 0 if not significant
- `dir_rank_i` = directional rank among significant genes in that direction (rank_up or rank_down)
- `normalized_rank_i` = 1 - (dir_rank_i / total_genes_in_experiment) if significant, else 0
- `contribution_i` = concordance_i * normalized_rank_i

`rank_score = mean(contribution_i)` over all signature genes present in D.

**Worked example with 4 genes:**

| Gene | Sig direction | DE status | dir_rank | total_genes | normalized_rank | concordance | contribution |
|------|--------------|-----------|----------|-------------|-----------------|-------------|-------------|
| G1 | up | significant_up | 3 | 1000 | 0.997 | +1 | +0.997 |
| G2 | down | significant_down | 10 | 1000 | 0.990 | +1 | +0.990 |
| G3 | up | significant_down | 5 | 1000 | 0.995 | -1 | -0.995 |
| G4 | down | not_significant | — | 1000 | 0 | 0 | 0 |

rank_score = mean(0.997, 0.990, -0.995, 0) = **0.248**

Log this worked example in the notebook.

- [ ] **Step 2: Write test_scoring.py with toy data**

```python
"""Toy data tests for scoring and permutation.

Hand-calculated expected values for rank score and permutation test.
"""

import pandas as pd
import numpy as np
import pytest
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent.parent))

from sig_utils.scoring import apply_signature, rank_score, permutation_test, score_with_significance


def make_toy_signature():
    """4-gene signature: 2 up, 2 down."""
    return pd.DataFrame([
        {"locus_tag": "G1", "gene_name": "up1", "direction": "up"},
        {"locus_tag": "G2", "gene_name": "down1", "direction": "down"},
        {"locus_tag": "G3", "gene_name": "up2", "direction": "up"},
        {"locus_tag": "G4", "gene_name": "down2", "direction": "down"},
    ])


def make_toy_de_all_concordant():
    """All 4 genes significant in expected direction."""
    return pd.DataFrame([
        {"locus_tag": "G1", "gene_name": "up1", "log2fc": 3.0,
         "expression_status": "significant_up", "rank": 2, "rank_up": 1, "rank_down": pd.NA,
         "timepoint": "tp1"},
        {"locus_tag": "G2", "gene_name": "down1", "log2fc": -2.5,
         "expression_status": "significant_down", "rank": 3, "rank_up": pd.NA, "rank_down": 1,
         "timepoint": "tp1"},
        {"locus_tag": "G3", "gene_name": "up2", "log2fc": 1.5,
         "expression_status": "significant_up", "rank": 5, "rank_up": 2, "rank_down": pd.NA,
         "timepoint": "tp1"},
        {"locus_tag": "G4", "gene_name": "down2", "log2fc": -1.0,
         "expression_status": "significant_down", "rank": 8, "rank_up": pd.NA, "rank_down": 3,
         "timepoint": "tp1"},
        # Background genes (not in signature)
        *[{"locus_tag": f"BG{i}", "gene_name": f"bg{i}", "log2fc": 0.1,
           "expression_status": "not_significant", "rank": 20+i, "rank_up": pd.NA, "rank_down": pd.NA,
           "timepoint": "tp1"} for i in range(96)],
    ])


def make_toy_de_mixed():
    """G1 concordant, G2 concordant, G3 reversed, G4 not significant."""
    return pd.DataFrame([
        {"locus_tag": "G1", "gene_name": "up1", "log2fc": 2.0,
         "expression_status": "significant_up", "rank": 3, "rank_up": 2, "rank_down": pd.NA,
         "timepoint": "tp1"},
        {"locus_tag": "G2", "gene_name": "down1", "log2fc": -1.5,
         "expression_status": "significant_down", "rank": 5, "rank_up": pd.NA, "rank_down": 3,
         "timepoint": "tp1"},
        {"locus_tag": "G3", "gene_name": "up2", "log2fc": -1.0,
         "expression_status": "significant_down", "rank": 8, "rank_up": pd.NA, "rank_down": 5,
         "timepoint": "tp1"},
        {"locus_tag": "G4", "gene_name": "down2", "log2fc": 0.2,
         "expression_status": "not_significant", "rank": 50, "rank_up": pd.NA, "rank_down": pd.NA,
         "timepoint": "tp1"},
        *[{"locus_tag": f"BG{i}", "gene_name": f"bg{i}", "log2fc": 0.1,
           "expression_status": "not_significant", "rank": 20+i, "rank_up": pd.NA, "rank_down": pd.NA,
           "timepoint": "tp1"} for i in range(96)],
    ])


def make_toy_de_gene_absent():
    """G1 present and concordant, G2 absent from dataset entirely."""
    rows = [
        {"locus_tag": "G1", "gene_name": "up1", "log2fc": 2.0,
         "expression_status": "significant_up", "rank": 3, "rank_up": 2, "rank_down": pd.NA,
         "timepoint": "tp1"},
        {"locus_tag": "G3", "gene_name": "up2", "log2fc": 1.0,
         "expression_status": "significant_up", "rank": 5, "rank_up": 3, "rank_down": pd.NA,
         "timepoint": "tp1"},
        {"locus_tag": "G4", "gene_name": "down2", "log2fc": -0.5,
         "expression_status": "not_significant", "rank": 40, "rank_up": pd.NA, "rank_down": pd.NA,
         "timepoint": "tp1"},
    ]
    rows.extend([
        {"locus_tag": f"BG{i}", "gene_name": f"bg{i}", "log2fc": 0.1,
         "expression_status": "not_significant", "rank": 20+i, "rank_up": pd.NA, "rank_down": pd.NA,
         "timepoint": "tp1"} for i in range(97)
    ])
    return pd.DataFrame(rows)


class TestApplySignature:
    def test_all_present(self):
        """All 4 signature genes found in DE data."""
        sig = make_toy_signature()
        de = make_toy_de_all_concordant()
        applied = apply_signature(sig, de)
        assert len(applied) == 4
        assert set(applied["locus_tag"]) == {"G1", "G2", "G3", "G4"}

    def test_concordance_correct(self):
        """Mixed DE: G1 concordant (+1), G3 reversed (-1), G4 not sig (0)."""
        sig = make_toy_signature()
        de = make_toy_de_mixed()
        applied = apply_signature(sig, de)
        g1 = applied[applied["locus_tag"] == "G1"].iloc[0]
        g3 = applied[applied["locus_tag"] == "G3"].iloc[0]
        g4 = applied[applied["locus_tag"] == "G4"].iloc[0]
        assert g1["concordance"] == 1
        assert g3["concordance"] == -1
        assert g4["concordance"] == 0

    def test_absent_gene(self):
        """G2 absent from dataset → row present with NaN DE values."""
        sig = make_toy_signature()
        de = make_toy_de_gene_absent()
        applied = apply_signature(sig, de)
        g2 = applied[applied["locus_tag"] == "G2"]
        assert len(g2) == 1
        assert pd.isna(g2.iloc[0]["log2fc"])


class TestRankScore:
    def test_all_concordant_positive(self):
        """All concordant → score > 0."""
        sig = make_toy_signature()
        de = make_toy_de_all_concordant()
        applied = apply_signature(sig, de)
        score = rank_score(applied, total_genes=100)
        assert score > 0

    def test_mixed_lower_than_concordant(self):
        """Mixed concordance → lower score than all concordant."""
        sig = make_toy_signature()
        de_conc = make_toy_de_all_concordant()
        de_mixed = make_toy_de_mixed()
        applied_conc = apply_signature(sig, de_conc)
        applied_mixed = apply_signature(sig, de_mixed)
        score_conc = rank_score(applied_conc, total_genes=100)
        score_mixed = rank_score(applied_mixed, total_genes=100)
        assert score_conc > score_mixed


class TestPermutationTest:
    def test_strong_signal_low_p(self):
        """All concordant with strong ranks → p should be low."""
        sig = make_toy_signature()
        de = make_toy_de_all_concordant()
        result = permutation_test(sig, de, n_perms=100, seed=42)
        assert result["empirical_p"] < 0.1

    def test_returns_expected_keys(self):
        sig = make_toy_signature()
        de = make_toy_de_all_concordant()
        result = permutation_test(sig, de, n_perms=50, seed=42)
        assert "observed" in result
        assert "empirical_p" in result
        assert "n_permutations" in result


class TestScoreWithSignificance:
    def test_returns_score_and_p(self):
        """Wrapper returns both score and p-value."""
        sig = make_toy_signature()
        de = make_toy_de_all_concordant()
        result = score_with_significance(sig, de, n_perms=50)
        assert "rank_score" in result
        assert "empirical_p" in result
```

- [ ] **Step 3: Run tests to verify they fail**

```bash
cd analyses/ANALYSIS_DIR && uv run python -m pytest sig_utils/tests/test_scoring.py -v
```

Expected: FAIL — `sig_utils.scoring` doesn't exist yet.

- [ ] **Step 4: Write scoring.py**

```python
"""Signature scoring: apply signature to DE data, compute rank score, permutation test.

All functions take DataFrames and return DataFrames or dicts. No KG dependency.
"""

import numpy as np
import pandas as pd


def apply_signature(
    signature_df: pd.DataFrame,
    de_df: pd.DataFrame,
) -> pd.DataFrame:
    """Apply a signature to a DE dataset, producing the scored intermediate.

    For each signature gene, looks it up in the DE data and computes
    concordance (does the DE direction match the signature direction?).

    Args:
        signature_df: Signature with at least locus_tag, direction columns.
        de_df: DE data for one experiment x timepoint. Must have locus_tag,
            expression_status, log2fc, rank, rank_up, rank_down.

    Returns:
        DataFrame with one row per signature gene, columns:
        locus_tag, direction (from signature), gene_name, log2fc,
        expression_status, rank_up, rank_down, de_direction,
        concordance (+1, -1, or 0), dir_rank, normalized_rank, contribution.
        Genes absent from DE get NaN for DE columns and concordance=0.
    """
    de_cols = ["locus_tag", "log2fc", "expression_status", "rank", "rank_up", "rank_down"]
    if "gene_name" in de_df.columns:
        de_cols.append("gene_name")

    sig_cols = ["locus_tag", "direction"]
    if "gene_name" in signature_df.columns:
        sig_cols.append("gene_name")

    merged = signature_df[sig_cols].merge(
        de_df[de_cols].drop_duplicates(subset=["locus_tag"]),
        on="locus_tag",
        how="left",
        suffixes=("_sig", "_de"),
    )

    # Resolve gene_name if both present
    if "gene_name_sig" in merged.columns:
        merged["gene_name"] = merged["gene_name_sig"].fillna(merged.get("gene_name_de"))
        merged = merged.drop(columns=[c for c in ["gene_name_sig", "gene_name_de"] if c in merged.columns])

    # DE direction
    merged["de_direction"] = merged["expression_status"].map({
        "significant_up": "up", "significant_down": "down",
    })

    # Concordance
    merged["concordance"] = np.where(
        merged["de_direction"] == merged["direction"], 1,
        np.where(
            merged["de_direction"].notna() & (merged["de_direction"] != merged["direction"]), -1,
            0
        )
    )

    # Directional rank
    merged["dir_rank"] = np.where(
        merged["expression_status"] == "significant_up", merged["rank_up"],
        np.where(
            merged["expression_status"] == "significant_down", merged["rank_down"],
            np.nan
        )
    )

    return merged


def rank_score(
    applied_df: pd.DataFrame,
    total_genes: int,
) -> float:
    """Compute rank score from an applied signature DataFrame.

    rank_score = mean(concordance_i * normalized_rank_i)
    where normalized_rank_i = 1 - (dir_rank_i / total_genes) if significant, else 0.

    Args:
        applied_df: Output of apply_signature.
        total_genes: Total genes in the experiment (for normalization).

    Returns:
        Float score. Positive = N-limitation signal, negative = reversal, 0 = no signal.
    """
    if len(applied_df) == 0:
        return np.nan

    normalized_rank = np.where(
        pd.notna(applied_df["dir_rank"]),
        1 - (applied_df["dir_rank"].astype(float) / total_genes),
        0,
    )

    contribution = applied_df["concordance"].values * normalized_rank
    return float(np.mean(contribution))


def permutation_test(
    signature_df: pd.DataFrame,
    de_df: pd.DataFrame,
    n_perms: int = 1000,
    seed: int = 42,
) -> dict:
    """Permutation test for rank score significance.

    Generates n_perms random gene sets of the same size as the signature,
    computes rank_score for each, returns empirical p-value.

    Args:
        signature_df: Signature with locus_tag, direction.
        de_df: Full DE data for one condition x timepoint (all genes, not just signature).
        n_perms: Number of permutations.
        seed: Random seed.

    Returns:
        dict with: observed, empirical_p, n_permutations, n_signature_genes.
    """
    # Compute observed score
    applied = apply_signature(signature_df, de_df)
    total_genes = de_df["locus_tag"].nunique()
    observed = rank_score(applied, total_genes)
    n_sig = len(signature_df)

    if np.isnan(observed) or n_sig < 30:
        return {
            "observed": observed,
            "empirical_p": np.nan,
            "n_permutations": 0,
            "n_signature_genes": n_sig,
        }

    rng = np.random.default_rng(seed)
    all_tags = de_df["locus_tag"].unique()
    directions = signature_df["direction"].values

    null_scores = np.empty(n_perms)
    for i in range(n_perms):
        random_tags = rng.choice(all_tags, size=n_sig, replace=False)
        random_sig = pd.DataFrame({
            "locus_tag": random_tags,
            "direction": directions,  # keep same direction distribution
        })
        random_applied = apply_signature(random_sig, de_df)
        null_scores[i] = rank_score(random_applied, total_genes)

    empirical_p = float((np.abs(null_scores) >= np.abs(observed)).mean())

    return {
        "observed": float(observed),
        "empirical_p": empirical_p,
        "n_permutations": n_perms,
        "n_signature_genes": n_sig,
    }


def score_with_significance(
    signature_df: pd.DataFrame,
    de_df: pd.DataFrame,
    n_perms: int = 1000,
    seed: int = 42,
) -> dict:
    """Score a signature against DE data with permutation p-value.

    Main entry point: wraps apply_signature + rank_score + permutation_test.

    Args:
        signature_df: Signature with locus_tag, direction.
        de_df: Full DE data for one condition x timepoint.
        n_perms: Number of permutations.
        seed: Random seed.

    Returns:
        dict with: rank_score, empirical_p, n_permutations, n_signature_genes,
        n_concordant, n_reversed, n_not_significant, n_absent.
    """
    applied = apply_signature(signature_df, de_df)
    total_genes = de_df["locus_tag"].nunique()
    score = rank_score(applied, total_genes)

    perm_result = permutation_test(signature_df, de_df, n_perms=n_perms, seed=seed)

    n_concordant = int((applied["concordance"] == 1).sum())
    n_reversed = int((applied["concordance"] == -1).sum())
    n_not_sig = int((applied["concordance"] == 0).sum() & applied["log2fc"].notna())
    n_absent = int(applied["log2fc"].isna().sum())

    return {
        "rank_score": score,
        "empirical_p": perm_result["empirical_p"],
        "n_permutations": perm_result["n_permutations"],
        "n_signature_genes": len(signature_df),
        "n_concordant": n_concordant,
        "n_reversed": n_reversed,
        "n_not_significant": n_not_sig,
        "n_absent": n_absent,
    }
```

- [ ] **Step 5: Run tests to verify they pass**

```bash
cd analyses/ANALYSIS_DIR && uv run python -m pytest sig_utils/tests/test_scoring.py -v
```

Expected: All tests PASS.

- [ ] **Step 6: Interactive verification (RESEARCHER CHECKPOINT)**

Walk through the worked example from Step 1 with the actual code. Verify edge cases:
- Gene absent from dataset
- Gene not significant
- Gene in opposite direction
- All concordant vs all discordant

- [ ] **Step 7: Commit**

```bash
git add sig_utils/scoring.py sig_utils/tests/test_scoring.py
git commit -m "feat: sig_utils scoring with rank score, permutation, and toy tests"
```

---

### Task 9: Script 05 — Score all experiments

**Files:**
- Create: `scripts/05_score_experiments.py`

- [ ] **Step 1: Write the scoring script**

```python
"""Step 6: Score all experiments (references, controls, targets) against signature.

Uses sig_utils/scoring.py for all computation.

Outputs:
    data/applied_*.csv — signature applied to each experiment
    results/scores_all.csv — all scores with role classification
    logs/06_score_experiments.log

Run from multiomics_research root:
    uv run analyses/ANALYSIS_DIR/scripts/05_score_experiments.py [--explore]
"""

import sys
import argparse
from pathlib import Path
from datetime import datetime

ANALYSIS_DIR = Path(__file__).parent.parent
sys.path.insert(0, str(ANALYSIS_DIR))

import pandas as pd
import numpy as np
from sig_utils.scoring import apply_signature, rank_score, score_with_significance
from sig_utils.extraction import check_marker_genes
from sig_utils.io import load_de_csv, load_signature_csv, save_csv, write_log

MARKERS = ["PMM0920", "PMM0370", "PMM0550", "PMM1452", "PMM0030"]

# --- FILL AFTER STEPS 1-4 ---
# Format: (de_filename, label, role, skip_timepoints)
# role: "reference", "positive_control", "negative_control", "target"
ALL_EXPERIMENTS = [
    # References
    # ("de_ref_tolonen_ndep.csv", "Tolonen N-dep", "reference", {"0h", "3h"}),
    # ("de_ref_read_ndep.csv", "Read N-dep", "reference", set()),

    # Negative controls
    # ("de_ctrl_tolonen_cyanate.csv", "Tolonen cyanate", "negative_control", set()),
    # ("de_ctrl_tolonen_urea.csv", "Tolonen urea", "negative_control", set()),

    # Targets
    # ("de_weissberg_rnaseq_axenic.csv", "Weissberg RNA-seq axenic", "target", set()),
    # ("de_weissberg_rnaseq_coculture.csv", "Weissberg RNA-seq coculture", "target", set()),
    # ("de_weissberg_proteomics_axenic.csv", "Weissberg proteomics axenic", "target", set()),
    # ("de_weissberg_proteomics_coculture.csv", "Weissberg proteomics coculture", "target", set()),
]

N_PERMS = 1000
SEED = 42


def score_experiment(
    signature_df: pd.DataFrame,
    de_df: pd.DataFrame,
    label: str,
    role: str,
    skip_timepoints: set,
    n_perms: int = N_PERMS,
) -> tuple[list[dict], list[pd.DataFrame]]:
    """Score all timepoints in one experiment.

    Returns:
        (score_rows, applied_dfs) — score rows for the results table,
        and applied DataFrames for saving.
    """
    score_rows = []
    applied_dfs = []

    if "timepoint" in de_df.columns and de_df["timepoint"].nunique() > 1:
        for tp, tp_df in de_df.groupby("timepoint"):
            if tp in skip_timepoints:
                continue

            result = score_with_significance(signature_df, tp_df, n_perms=n_perms, seed=SEED)
            applied = apply_signature(signature_df, tp_df)
            total_genes = tp_df["locus_tag"].nunique()

            result.update({
                "label": label,
                "role": role,
                "timepoint": tp,
                "timepoint_hours": tp_df["timepoint_hours"].iloc[0] if "timepoint_hours" in tp_df.columns else np.nan,
                "total_genes_in_experiment": total_genes,
            })
            score_rows.append(result)
            applied_dfs.append(applied.assign(timepoint=tp, label=label))
    else:
        result = score_with_significance(signature_df, de_df, n_perms=n_perms, seed=SEED)
        applied = apply_signature(signature_df, de_df)
        total_genes = de_df["locus_tag"].nunique()

        result.update({
            "label": label,
            "role": role,
            "timepoint": "single",
            "timepoint_hours": np.nan,
            "total_genes_in_experiment": total_genes,
        })
        score_rows.append(result)
        applied_dfs.append(applied.assign(timepoint="single", label=label))

    return score_rows, applied_dfs


def main():
    parser = argparse.ArgumentParser(description="Score all experiments against signature")
    parser.add_argument("--explore", action="store_true", help="Print marker gene traces")
    args = parser.parse_args()

    log_lines = [f"Step 6: Score All Experiments — {datetime.now().isoformat()}", ""]

    # Load signature
    core = load_signature_csv("core_signature.csv")
    print(f"Core signature: {len(core)} genes")
    log_lines.append(f"Core signature: {len(core)} genes")

    all_score_rows = []
    all_applied = []

    for de_filename, label, role, skip_tps in ALL_EXPERIMENTS:
        print(f"\nScoring {label} ({role})...")
        de = load_de_csv(de_filename)

        score_rows, applied_dfs = score_experiment(core, de, label, role, skip_tps)
        all_score_rows.extend(score_rows)
        all_applied.extend(applied_dfs)

        for row in score_rows:
            print(f"  {row['timepoint']}: rank_score={row['rank_score']:.3f}, p={row['empirical_p']:.3f}, "
                  f"concordant={row['n_concordant']}, reversed={row['n_reversed']}")
            log_lines.append(f"  {label} {row['timepoint']}: score={row['rank_score']:.3f}, "
                             f"p={row['empirical_p']:.3f}")

        # Save applied subset
        applied_filename = f"applied_{de_filename}"
        combined_applied = pd.concat(applied_dfs, ignore_index=True)
        save_csv(combined_applied, applied_filename)

        if args.explore:
            marker_applied = combined_applied[combined_applied["locus_tag"].isin(MARKERS)]
            if len(marker_applied) > 0:
                print(f"\n  Marker genes in {label}:")
                cols = ["locus_tag", "direction", "timepoint", "log2fc",
                        "concordance", "dir_rank"]
                cols = [c for c in cols if c in marker_applied.columns]
                print(marker_applied[cols].to_string(index=False))
                log_lines.append(f"\nMarker genes:")
                log_lines.append(marker_applied[cols].to_string(index=False))

    # Save all scores
    scores_df = pd.DataFrame(all_score_rows)
    save_csv(scores_df, "scores_all.csv", subdir="results")

    # Control separation summary
    print("\n=== Control separation ===")
    log_lines.append("\n=== Control separation ===")
    for role_name in ["reference", "positive_control", "negative_control", "target"]:
        role_scores = scores_df[scores_df["role"] == role_name]
        if len(role_scores) > 0:
            mean_s = role_scores["rank_score"].mean()
            line = f"  {role_name}: n={len(role_scores)}, mean_rank_score={mean_s:.3f}, range=[{role_scores['rank_score'].min():.3f}, {role_scores['rank_score'].max():.3f}]"
            print(line)
            log_lines.append(line)

    write_log("\n".join(log_lines), "05_score_experiments.log")
    print(f"\nOutputs: results/scores_all.csv, data/applied_*.csv")


if __name__ == "__main__":
    main()
```

- [ ] **Step 2: Fill ALL_EXPERIMENTS after Steps 1-4**

- [ ] **Step 3: Run the script**

```bash
uv run analyses/ANALYSIS_DIR/scripts/05_score_experiments.py --explore
```

- [ ] **Step 4: Interactive exploration (RESEARCHER CHECKPOINT)**

In chat:
- Full results table: do positive controls score high? Negative controls low?
- **Control separation gate:** Is the gap between positive and negative controls clear? If not, the metric or signature needs adjustment before proceeding.
- Trace marker genes through scoring for specific conditions (axenic day 14, coculture day 18)
- Compare RNA-seq vs proteomics scores for the same conditions
- Identify surprises: any negative controls scoring high? Any targets unexpected?
- Compare with v1 results

- [ ] **Step 5: Write notebook entry, update manifests, commit**

```bash
git add scripts/05_score_experiments.py results/scores_all.csv data/applied_*.csv logs/ exploration/ results/RESULTS_MANIFEST.md
git commit -m "feat: Step 6 — score all experiments with rank score and permutation tests"
```

---

### Task 10: Script 06 — Plot results

**Files:**
- Create: `scripts/06_plot_results.py`

- [ ] **Step 1: Write the plotting script**

```python
"""Step 7: Visualize scoring results.

Produces:
    results/trajectory_rnaseq.png — rank score over time, RNA-seq
    results/trajectory_proteomics.png — rank score over time, proteomics
    results/control_separation.png — positive vs negative controls

Run from multiomics_research root:
    uv run analyses/ANALYSIS_DIR/scripts/06_plot_results.py
"""

import sys
from pathlib import Path

ANALYSIS_DIR = Path(__file__).parent.parent
sys.path.insert(0, str(ANALYSIS_DIR))

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

from sig_utils.io import RESULTS_DIR

RESULTS_DIR.mkdir(exist_ok=True)


def load_scores():
    """Load scores_all.csv from results/."""
    return pd.read_csv(RESULTS_DIR / "scores_all.csv")


def plot_trajectories(scores_df: pd.DataFrame):
    """Plot rank score trajectories for RNA-seq and proteomics."""

    for platform_label, platform_filter in [("rnaseq", "rnaseq"), ("proteomics", "proteomics")]:
        # Filter to target experiments for this platform
        # The exact filtering depends on how labels are structured after Step 1
        # This is a template — adjust column names after seeing actual data
        target = scores_df[
            (scores_df["role"] == "target") &
            (scores_df["label"].str.contains(platform_label, case=False))
        ].copy()

        if target.empty:
            print(f"  No {platform_label} target data, skipping trajectory plot")
            continue

        # Reference baselines
        ref = scores_df[scores_df["role"].isin(["reference", "positive_control"])].copy()

        fig, ax = plt.subplots(1, 1, figsize=(10, 6))

        # Plot reference range as a band
        if not ref.empty:
            ref_min = ref["rank_score"].min()
            ref_max = ref["rank_score"].max()
            ax.axhspan(ref_min, ref_max, alpha=0.1, color="green", label="Reference range")

        # Plot target trajectories
        for label_val in target["label"].unique():
            subset = target[target["label"] == label_val].sort_values("timepoint_hours")
            if subset["timepoint_hours"].notna().any():
                ax.plot(subset["timepoint_hours"], subset["rank_score"],
                        marker="o", label=label_val)
            else:
                # Single timepoint — plot as horizontal marker
                ax.axhline(y=subset["rank_score"].iloc[0],
                           linestyle="--", alpha=0.7, label=label_val)

        ax.set_xlabel("Time (hours)")
        ax.set_ylabel("Rank Score")
        ax.set_title(f"N-Limitation Signature Score — {platform_label.upper()}")
        ax.axhline(y=0, color="gray", linestyle=":", alpha=0.5)
        ax.legend(bbox_to_anchor=(1.05, 1), loc="upper left")
        fig.tight_layout()

        outpath = RESULTS_DIR / f"trajectory_{platform_label}.png"
        fig.savefig(outpath, dpi=300, bbox_inches="tight")
        plt.close(fig)
        print(f"  Saved {outpath}")


def plot_control_separation(scores_df: pd.DataFrame):
    """Plot positive vs negative control scores."""
    fig, ax = plt.subplots(1, 1, figsize=(8, 6))

    colors = {
        "reference": "green",
        "positive_control": "darkgreen",
        "negative_control": "red",
        "target": "blue",
    }

    for role in ["reference", "positive_control", "negative_control", "target"]:
        subset = scores_df[scores_df["role"] == role]
        if subset.empty:
            continue
        ax.scatter(
            [role] * len(subset),
            subset["rank_score"],
            c=colors.get(role, "gray"),
            s=80,
            alpha=0.7,
            label=role,
        )

    ax.set_ylabel("Rank Score")
    ax.set_title("Control Separation: Rank Score by Experiment Role")
    ax.axhline(y=0, color="gray", linestyle=":", alpha=0.5)
    fig.tight_layout()

    outpath = RESULTS_DIR / f"control_separation.png"
    fig.savefig(outpath, dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved {outpath}")


def main():
    print("Loading scores...")
    scores = load_scores()
    print(f"  {len(scores)} score rows loaded")

    print("\nPlotting trajectories...")
    plot_trajectories(scores)

    print("\nPlotting control separation...")
    plot_control_separation(scores)

    print("\nDone.")


if __name__ == "__main__":
    main()
```

- [ ] **Step 2: Run the script**

```bash
uv run analyses/ANALYSIS_DIR/scripts/06_plot_results.py
```

- [ ] **Step 3: Interactive exploration (RESEARCHER CHECKPOINT)**

In chat:
- Review trajectory figures: does the pattern match expectations?
- Review control separation: clean gap between positive and negative?
- Discuss biological interpretation
- RNA/protein discordance: does it show up in v2 as it did in v1?
- Flag anything unexpected

- [ ] **Step 4: Decide on companion Jupyter notebook**

Ask the researcher: "We have all the data and figures. Would a companion Jupyter notebook for interactive exploration be useful, or is the markdown notebook sufficient?"

If yes, create `exploration/notebooks/signature_explorer.ipynb`.

- [ ] **Step 5: Write notebook entry, update manifests, commit**

```bash
git add scripts/06_plot_results.py results/*.png logs/ exploration/ results/RESULTS_MANIFEST.md
git commit -m "feat: Step 7 — trajectory and control separation plots"
```

---

### Task 11: Finalize documentation

**Files:**
- Modify: `methods.md`, `decisions.md`, `caveats.md`, `gaps_and_friction.md`, `references.md`, `README.md`

- [ ] **Step 1: Update methods.md with final results**

Fill in all sections: research question, data scope (with experiment IDs and DOIs from Step 1), gene selection (filter funnel from Step 3), statistical methods (rank score formula, permutation test), results summary, limitations.

- [ ] **Step 2: Update decisions.md with choices made during analysis**

Add entries for each decision from the interactive sessions: study selection, timepoint exclusions, direction handling, scoring metric, permutation parameters.

- [ ] **Step 3: Update caveats.md**

Add all interpretation caveats discovered during the analysis: platform coverage, acute vs chronic, self-scoring circularity, table_scope filtering, etc.

- [ ] **Step 4: Update gaps_and_friction.md**

Log any KG bugs, gaps, MCP friction, or methodology issues discovered during this rebuild.

- [ ] **Step 5: Update README.md with key findings and file index**

- [ ] **Step 6: Update references.md with DOIs and versions**

- [ ] **Step 7: Assess sig_utils for package readiness**

Add a section to the notebook: which functions are package-ready? What needs refinement? What's the proposed API for `multiomics_explorer.analysis`?

- [ ] **Step 8: Final commit**

```bash
git add .
git commit -m "docs: finalize analysis documentation and package assessment"
```

---

## Self-Review Checklist

**Spec coverage:**
- [x] Step 1: KG discovery → Task 3
- [x] Step 2: Reference DE extraction → Task 4 (uses extraction utility from Task 2)
- [x] Step 3: Signature building → Task 5 (sig_utils/signature.py) + Task 6 (script)
- [x] Step 4: Target DE extraction → Task 7
- [x] Step 5: Metric design and verification → Task 8 (sig_utils/scoring.py)
- [x] Step 6: Score all experiments → Task 9
- [x] Step 7: Interpret and document → Task 10 + Task 11
- [x] sig_utils/extraction.py → Task 2
- [x] sig_utils/io.py → Task 2
- [x] sig_utils/tests/ → Task 5 (test_signature.py) + Task 8 (test_scoring.py)
- [x] Superpowers products in analysis dir → Task 1
- [x] Brainstorm log → Task 1
- [x] Companion Jupyter notebook (optional) → Task 10
- [x] Positive/negative control validation gate → Task 9 step 4
- [x] table_scope from envelope → Task 2 extraction.py
- [x] Marker gene traces at every step → all scripts have --explore flag

**Placeholder scan:** No TBDs. Experiment IDs are intentionally left as comments (filled interactively after Step 1) — this is a feature, not a placeholder.

**Type consistency:** `apply_signature` returns DataFrame with `concordance`, `dir_rank` columns — used consistently in scoring.py and test_scoring.py. `score_with_significance` returns dict with `rank_score`, `empirical_p` — used consistently in test and script.
