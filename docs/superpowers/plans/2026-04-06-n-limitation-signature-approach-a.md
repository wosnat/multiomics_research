# N-Limitation Signature Analysis (Approach A) Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Build an N-limitation gene signature from reference studies and score Weissberg 2025 axenic vs coculture experiments against it to quantify nitrogen limitation severity over time.

**Architecture:** Extract DE data from KG via Python API → build core/extended signatures by intersection → score each Weissberg condition/timepoint with reusable `sig_utils` metrics → visualize trajectories with reference baselines → annotate by pathway → check Alteromonas N-recycling support.

**Tech Stack:** Python 3, pandas, scipy (Spearman), matplotlib, multiomics_explorer Python API

**Spec:** `docs/superpowers/specs/2026-04-06-n-limitation-signature-analysis-design.md`

---

## File Structure

```
analyses/2026-04-06-1432-n_limitation_signature/
├── sig_utils/
│   ├── __init__.py
│   ├── io.py                  # load_de_csv, save_scores_csv, load_signature_csv
│   ├── signature.py           # intersect_de_lists, assign_direction, assign_ranks, build_core_signature
│   └── metrics.py             # hit_rate, mean_signed_log2fc, mean_signed_normalized_rank, rank_correlation
├── scripts/
│   ├── extract_reference_de.py
│   ├── extract_weissberg_de.py
│   ├── build_signature.py
│   ├── score_signature.py
│   ├── plot_trajectories.py
│   └── explore_alteromonas.py
├── data/                      # created by extraction scripts
├── results/                   # created by scoring/plotting scripts
├── exploration/
├── superpowers/               # spec snapshot
├── methods.md
├── gaps_and_friction.md
├── references.md
└── README.md
```

All scripts run from `multiomics_research` root via `uv run analyses/2026-04-06-1432-n_limitation_signature/scripts/<script>.py`.

---

### Task 1: Scaffold the analysis directory and logging files

**Files:**
- Create: `analyses/2026-04-06-1432-n_limitation_signature/` (full tree)

- [ ] **Step 1: Create directory structure**

```bash
mkdir -p analyses/2026-04-06-1432-n_limitation_signature/{sig_utils,scripts,data,results,exploration,superpowers}
```

- [ ] **Step 2: Copy spec and plan snapshots**

```bash
cp docs/superpowers/specs/2026-04-06-n-limitation-signature-analysis-design.md \
   analyses/2026-04-06-1432-n_limitation_signature/superpowers/
cp docs/superpowers/plans/2026-04-06-n-limitation-signature-approach-a.md \
   analyses/2026-04-06-1432-n_limitation_signature/superpowers/
```

- [ ] **Step 3: Create methods.md skeleton**

Write to `analyses/2026-04-06-1432-n_limitation_signature/methods.md`:

```markdown
# N-Limitation Signature Analysis: Methods

## Research question

Can we quantify the degree of nitrogen limitation in Prochlorococcus MED4
molecularly, and show that coculture with Alteromonas alleviates it over
time, using a gene signature derived from independent reference studies?

## Data scope

### Reference studies
- Tolonen et al. 2006 (DOI: 10.1038/msb4100087) — MED4 N-deprivation time course (microarray, 0-48h)
- Read et al. 2017 (DOI: 10.1038/ismej.2017.88) — MED4 N-depleted vs N-replete (RNA-seq, 3-24h)

### Target study
- Weissberg et al. 2025 (DOI: 10.1101/2025.11.24.690089) — MED4 in PRO99-lowN, axenic vs coculture with Alteromonas HOT1A3 (RNA-seq + proteomics, d14-89)

### Inclusions/exclusions
- Tolonen timepoints 0h and 3h excluded (nearly zero DE genes)
- Read dataset filtered to top 50% of genes by expression level (table_scope: filtered_subset)
- Weissberg d60+89 combined timepoint included in data but excluded from trajectory plots

## Gene selection

(Updated as analysis proceeds)

## Statistical methods

(Updated as analysis proceeds)

## Results summary

(Updated as analysis proceeds)

## Limitations

(Updated as analysis proceeds)
```

- [ ] **Step 4: Create gaps_and_friction.md skeleton**

Write to `analyses/2026-04-06-1432-n_limitation_signature/gaps_and_friction.md`:

```markdown
# Gaps and friction points

## KG data bugs

## KG gaps

## MCP friction

## Skill/methodology friction

## Process retrospective

### What worked

### What didn't work

### Proposed changes

**To the skill:**

**To the MCP/KG:**
```

- [ ] **Step 5: Create references.md**

Write to `analyses/2026-04-06-1432-n_limitation_signature/references.md`:

```markdown
# References

## Data publications
- Tolonen et al. 2006. Global gene expression of Prochlorococcus ecotypes in response to changes in nitrogen availability. Molecular Systems Biology. DOI: 10.1038/msb4100087
- Read et al. 2017. Nitrogen cost minimization is promoted by structural changes in the transcriptome of N-deprived Prochlorococcus cells. The ISME Journal. DOI: 10.1038/ismej.2017.88
- Weissberg et al. 2025. Transcriptomic and Proteomic Analysis Reveals Nitrogen Recycling as a Core Mechanism for Prochlorococcus Prolonged Survival. bioRxiv. DOI: 10.1101/2025.11.24.690089

## Tools and software
- multiomics_explorer Python API (version: check with `kg_schema()`)
- KG build date: (fill in from `kg_schema()` output)
- Python: (fill in from `sys.version`)
- pandas, scipy, matplotlib: (fill in versions)
```

- [ ] **Step 6: Create README.md skeleton**

Write to `analyses/2026-04-06-1432-n_limitation_signature/README.md`:

```markdown
# N-Limitation Signature Analysis

**Question:** Is Prochlorococcus nitrogen-limited in Weissberg 2025's axenic
and coculture experiments, and how does the severity differ over time?

**Approach:** Reference signature scoring (Approach A of 3).

## Key findings

(Updated as analysis proceeds)

## File index

- `methods.md` — publication-ready methods
- `gaps_and_friction.md` — issues log
- `references.md` — citations and versions
- `data/` — staged KG extracts
- `scripts/` — extraction and analysis scripts
- `sig_utils/` — reusable scoring utilities
- `results/` — output tables and figures
- `exploration/` — exploration logs

## Exploration logs

(Links added as iterations proceed)
```

- [ ] **Step 7: Commit**

```bash
git add analyses/2026-04-06-1432-n_limitation_signature/
git commit -m "scaffold: analysis directory for N-limitation signature"
```

---

### Task 2: Build `sig_utils/io.py`

**Files:**
- Create: `analyses/2026-04-06-1432-n_limitation_signature/sig_utils/__init__.py`
- Create: `analyses/2026-04-06-1432-n_limitation_signature/sig_utils/io.py`
- Test: inline `assert` checks (no separate test dir — analysis project)

- [ ] **Step 1: Create `__init__.py`**

Write to `analyses/2026-04-06-1432-n_limitation_signature/sig_utils/__init__.py`:

```python
"""Reusable building blocks for signature-based scoring."""
```

- [ ] **Step 2: Write `io.py`**

Write to `analyses/2026-04-06-1432-n_limitation_signature/sig_utils/io.py`:

```python
"""I/O helpers for signature analysis.

Load/save DE data and signature CSVs with consistent column types.
"""

import pandas as pd
from pathlib import Path

# Base path for data/ and results/ directories
ANALYSIS_DIR = Path(__file__).parent.parent
DATA_DIR = ANALYSIS_DIR / "data"
RESULTS_DIR = ANALYSIS_DIR / "results"


def load_de_csv(filename: str) -> pd.DataFrame:
    """Load a DE extract CSV from data/.

    Ensures consistent types for rank columns (nullable int)
    and expression_status (string).
    """
    path = DATA_DIR / filename
    df = pd.read_csv(path)

    # rank_up and rank_down are nullable int (None when not significant in that direction)
    for col in ["rank", "rank_up", "rank_down"]:
        if col in df.columns:
            df[col] = pd.to_numeric(df[col], errors="coerce").astype("Int64")

    return df


def load_signature_csv(filename: str) -> pd.DataFrame:
    """Load a signature definition CSV from data/.

    Expected columns: locus_tag, gene_name, direction, signature_type,
    plus rank columns.
    """
    path = DATA_DIR / filename
    df = pd.read_csv(path)

    for col in df.columns:
        if "rank" in col:
            df[col] = pd.to_numeric(df[col], errors="coerce").astype("Int64")

    return df


def save_scores_csv(df: pd.DataFrame, filename: str) -> Path:
    """Save a scores DataFrame to results/."""
    RESULTS_DIR.mkdir(exist_ok=True)
    path = RESULTS_DIR / filename
    df.to_csv(path, index=False)
    return path


def save_data_csv(df: pd.DataFrame, filename: str) -> Path:
    """Save a DataFrame to data/."""
    DATA_DIR.mkdir(exist_ok=True)
    path = DATA_DIR / filename
    df.to_csv(path, index=False)
    return path
```

- [ ] **Step 3: Verify io module imports**

```bash
uv run python -c "
from analyses.sig_utils_import_test import test
" 2>&1
```

Actually, since the analysis dir isn't a package on sys.path, verify with a direct import:

```bash
uv run python -c "
import sys; sys.path.insert(0, 'analyses/2026-04-06-1432-n_limitation_signature')
from sig_utils.io import load_de_csv, save_data_csv, ANALYSIS_DIR, DATA_DIR
print('ANALYSIS_DIR:', ANALYSIS_DIR)
print('DATA_DIR:', DATA_DIR)
assert DATA_DIR.name == 'data'
print('OK')
"
```

Expected: prints paths and `OK`.

- [ ] **Step 4: Commit**

```bash
git add analyses/2026-04-06-1432-n_limitation_signature/sig_utils/
git commit -m "feat: sig_utils/io.py — I/O helpers for signature analysis"
```

---

### Task 3: Build `sig_utils/signature.py`

**Files:**
- Create: `analyses/2026-04-06-1432-n_limitation_signature/sig_utils/signature.py`

- [ ] **Step 1: Write `signature.py`**

Write to `analyses/2026-04-06-1432-n_limitation_signature/sig_utils/signature.py`:

```python
"""Signature building: intersect DE lists, assign directions and ranks."""

import pandas as pd


def summarize_de_per_gene(de_df: pd.DataFrame) -> pd.DataFrame:
    """Summarize a multi-timepoint DE DataFrame to one row per gene.

    For each gene, picks the timepoint with max |log2fc| among significant rows.
    Records: direction, peak log2fc, peak timepoint, best directional rank,
    best global rank.

    Args:
        de_df: DataFrame with columns: locus_tag, gene_name, timepoint,
            log2fc, expression_status, rank, rank_up, rank_down.

    Returns:
        DataFrame with one row per gene: locus_tag, gene_name, direction,
        peak_log2fc, peak_timepoint_by_rank (primary),
        peak_timepoint_by_fc (secondary), best_dir_rank, best_global_rank.
    """
    sig = de_df[de_df["expression_status"].isin(["significant_up", "significant_down"])].copy()
    if sig.empty:
        return pd.DataFrame(columns=[
            "locus_tag", "gene_name", "direction", "peak_log2fc",
            "peak_timepoint_by_rank", "peak_timepoint_by_fc",
            "best_dir_rank", "best_global_rank",
        ])

    sig["abs_log2fc"] = sig["log2fc"].abs()
    sig["direction"] = sig["expression_status"].map({
        "significant_up": "up", "significant_down": "down"
    })
    sig["dir_rank"] = sig.apply(
        lambda r: r["rank_up"] if r["direction"] == "up" else r["rank_down"], axis=1
    )

    # For each gene, determine majority direction across timepoints
    direction_counts = sig.groupby("locus_tag")["direction"].value_counts().unstack(fill_value=0)
    majority_direction = direction_counts.idxmax(axis=1).rename("majority_direction")

    # Filter to majority direction
    sig = sig.merge(majority_direction, on="locus_tag")
    sig_majority = sig[sig["direction"] == sig["majority_direction"]]

    # Peak by rank (primary): timepoint where gene has best directional rank
    peak_by_rank = sig_majority.loc[sig_majority.groupby("locus_tag")["dir_rank"].idxmin()]

    # Peak by log2fc (secondary): timepoint where gene has max |log2fc|
    peak_by_fc = sig_majority.loc[sig_majority.groupby("locus_tag")["abs_log2fc"].idxmax()]

    # Best directional rank across all timepoints in majority direction
    best_dir = sig_majority.groupby("locus_tag")["dir_rank"].min().rename("best_dir_rank")

    # Best global rank across all significant timepoints
    best_global = sig.groupby("locus_tag")["rank"].min().rename("best_global_rank")

    # Use rank-based peak as primary, record both peak timepoints
    result = peak_by_rank[["locus_tag", "gene_name", "direction", "log2fc", "timepoint"]].copy()
    result = result.rename(columns={"log2fc": "peak_log2fc", "timepoint": "peak_timepoint_by_rank"})
    result["direction"] = result["locus_tag"].map(majority_direction)

    # Add fc-based peak timepoint for reference
    fc_peaks = peak_by_fc[["locus_tag", "timepoint"]].rename(columns={"timepoint": "peak_timepoint_by_fc"})
    result = result.merge(fc_peaks, on="locus_tag")

    result = result.merge(best_dir, on="locus_tag")
    result = result.merge(best_global, on="locus_tag")

    return result.reset_index(drop=True)


def intersect_de_lists(
    study_a: pd.DataFrame,
    study_b: pd.DataFrame,
    study_b_all_locus_tags: set[str] | None = None,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Intersect two per-gene DE summaries into core and extended signatures.

    Core: genes in both studies with concordant direction.
    Extended: genes in only one study, tagged by source.

    Args:
        study_a: Output of summarize_de_per_gene for study A (Tolonen).
        study_b: Output of summarize_de_per_gene for study B (Read).
        study_b_all_locus_tags: Set of all locus tags present in study B's
            dataset (not just significant). If None, cannot distinguish
            'absent from dataset' vs 'present but not significant'.

    Returns:
        (core_df, extended_df) — each with columns:
            locus_tag, gene_name, direction,
            study_a_peak_log2fc, study_a_peak_timepoint_by_rank, study_a_peak_timepoint_by_fc, study_a_best_dir_rank, study_a_best_global_rank,
            study_b_peak_log2fc, study_b_peak_timepoint_by_rank, study_b_peak_timepoint_by_fc, study_b_best_dir_rank, study_b_best_global_rank,
            cross_study_best_dir_rank,
            signature_type ('core', 'tolonen_only_read_absent', 'tolonen_only_read_ns', 'read_only')
    """
    a = study_a.add_prefix("study_a_").rename(columns={
        "study_a_locus_tag": "locus_tag", "study_a_gene_name": "gene_name",
        "study_a_direction": "direction_a"
    })
    b = study_b.add_prefix("study_b_").rename(columns={
        "study_b_locus_tag": "locus_tag", "study_b_gene_name": "gene_name",
        "study_b_direction": "direction_b"
    })

    merged = a.merge(b, on="locus_tag", how="outer", suffixes=("", "_b"))
    # Use gene_name from whichever study has it
    merged["gene_name"] = merged["gene_name"].fillna(merged.get("gene_name_b"))
    if "gene_name_b" in merged.columns:
        merged = merged.drop(columns=["gene_name_b"])

    # Core: both present, concordant direction
    both = merged.dropna(subset=["direction_a", "direction_b"])
    concordant = both[both["direction_a"] == both["direction_b"]].copy()
    concordant["direction"] = concordant["direction_a"]
    concordant["signature_type"] = "core"

    # Cross-study best directional rank
    concordant["cross_study_best_dir_rank"] = concordant[
        ["study_a_best_dir_rank", "study_b_best_dir_rank"]
    ].min(axis=1)

    # Extended: in one study only
    a_only_mask = merged["direction_b"].isna() & merged["direction_a"].notna()
    b_only_mask = merged["direction_a"].isna() & merged["direction_b"].notna()
    discordant_mask = (
        merged["direction_a"].notna() & merged["direction_b"].notna()
        & (merged["direction_a"] != merged["direction_b"])
    )

    a_only = merged[a_only_mask].copy()
    a_only["direction"] = a_only["direction_a"]
    if study_b_all_locus_tags is not None:
        a_only["signature_type"] = a_only["locus_tag"].apply(
            lambda lt: "tolonen_only_read_ns" if lt in study_b_all_locus_tags
            else "tolonen_only_read_absent"
        )
    else:
        a_only["signature_type"] = "tolonen_only"

    b_only = merged[b_only_mask].copy()
    b_only["direction"] = b_only["direction_b"]
    b_only["signature_type"] = "read_only"

    extended = pd.concat([a_only, b_only], ignore_index=True)
    extended["cross_study_best_dir_rank"] = extended.apply(
        lambda r: r["study_a_best_dir_rank"] if pd.notna(r["study_a_best_dir_rank"])
        else r["study_b_best_dir_rank"],
        axis=1,
    )

    # Select and order columns
    cols = [
        "locus_tag", "gene_name", "direction", "signature_type",
        "study_a_peak_log2fc", "study_a_peak_timepoint_by_rank", "study_a_peak_timepoint_by_fc",
        "study_a_best_dir_rank", "study_a_best_global_rank",
        "study_b_peak_log2fc", "study_b_peak_timepoint_by_rank", "study_b_peak_timepoint_by_fc",
        "study_b_best_dir_rank", "study_b_best_global_rank",
        "cross_study_best_dir_rank",
    ]
    for df in [concordant, extended]:
        for c in cols:
            if c not in df.columns:
                df[c] = pd.NA

    core_df = concordant[cols].sort_values("cross_study_best_dir_rank").reset_index(drop=True)
    extended_df = extended[cols].sort_values("cross_study_best_dir_rank").reset_index(drop=True)

    return core_df, extended_df
```

- [ ] **Step 2: Verify with a smoke test**

```bash
uv run python -c "
import sys; sys.path.insert(0, 'analyses/2026-04-06-1432-n_limitation_signature')
import pandas as pd
from sig_utils.signature import summarize_de_per_gene, intersect_de_lists

# Fake DE data
a_data = pd.DataFrame({
    'locus_tag': ['PMM0920', 'PMM0920', 'PMM0246', 'PMM0100'],
    'gene_name': ['glnA', 'glnA', 'ntcA', 'fake'],
    'timepoint': ['6h', '12h', '6h', '6h'],
    'log2fc': [5.0, 4.0, 3.0, -2.0],
    'expression_status': ['significant_up', 'significant_up', 'significant_up', 'significant_down'],
    'rank': [5, 9, 7, 10],
    'rank_up': [5, 4, 6, None],
    'rank_down': [None, None, None, 3],
})
b_data = pd.DataFrame({
    'locus_tag': ['PMM0920', 'PMM0246', 'PMM0200'],
    'gene_name': ['glnA', 'ntcA', 'other'],
    'timepoint': ['12h', '12h', '12h'],
    'log2fc': [3.0, 2.0, -1.5],
    'expression_status': ['significant_up', 'significant_up', 'significant_down'],
    'rank': [3, 10, 20],
    'rank_up': [2, 5, None],
    'rank_down': [None, None, 8],
})

sa = summarize_de_per_gene(a_data)
sb = summarize_de_per_gene(b_data)
print('Study A summary:'); print(sa.to_string())
print('Study B summary:'); print(sb.to_string())

core, ext = intersect_de_lists(sa, sb, study_b_all_locus_tags={'PMM0920','PMM0246','PMM0200'})
print('Core:'); print(core[['locus_tag','direction','signature_type','cross_study_best_dir_rank']].to_string())
print('Extended:'); print(ext[['locus_tag','direction','signature_type','cross_study_best_dir_rank']].to_string())
assert len(core) == 2  # glnA and ntcA
assert set(core['locus_tag']) == {'PMM0920', 'PMM0246'}
assert len(ext) == 2  # PMM0100 (tolonen_only_read_absent) and PMM0200 (read_only)
print('OK')
"
```

Expected: prints summaries and `OK`.

- [ ] **Step 3: Commit**

```bash
git add analyses/2026-04-06-1432-n_limitation_signature/sig_utils/signature.py
git commit -m "feat: sig_utils/signature.py — intersect DE lists, build core/extended signatures"
```

---

### Task 4: Build `sig_utils/metrics.py`

**Files:**
- Create: `analyses/2026-04-06-1432-n_limitation_signature/sig_utils/metrics.py`

- [ ] **Step 1: Write `metrics.py`**

Write to `analyses/2026-04-06-1432-n_limitation_signature/sig_utils/metrics.py`:

```python
"""Scoring metrics for signature activation.

All functions take a DE DataFrame (genes x timepoint) and a signature
DataFrame, and return scores. No KG dependency — pure computation on
staged data.
"""

import numpy as np
import pandas as pd
from scipy import stats


def _merge_signature_with_de(
    signature_df: pd.DataFrame,
    de_df: pd.DataFrame,
) -> pd.DataFrame:
    """Merge signature genes with DE data for one experiment x timepoint.

    Returns merged DataFrame with signature direction and DE results.
    Signature genes not found in DE get NaN for DE columns.
    """
    return signature_df[["locus_tag", "direction"]].merge(
        de_df[["locus_tag", "log2fc", "expression_status", "rank", "rank_up", "rank_down"]],
        on="locus_tag",
        how="left",
    )


def hit_rate(
    signature_df: pd.DataFrame,
    de_df: pd.DataFrame,
) -> dict:
    """Fraction of signature genes significant in expected direction.

    Returns:
        dict with keys: concordant_hits, reversed_hits, not_significant,
        total, hit_rate_concordant, hit_rate_reversed.
    """
    m = _merge_signature_with_de(signature_df, de_df)
    total = len(m)
    if total == 0:
        return {
            "concordant_hits": 0, "reversed_hits": 0, "not_significant": 0,
            "total": 0, "hit_rate_concordant": np.nan, "hit_rate_reversed": np.nan,
        }

    m["de_direction"] = m["expression_status"].map({
        "significant_up": "up", "significant_down": "down",
    })

    concordant = ((m["de_direction"] == m["direction"]) & m["de_direction"].notna()).sum()
    reversed_ = (
        (m["de_direction"] != m["direction"])
        & m["de_direction"].notna()
        & m["direction"].notna()
    ).sum()
    not_sig = total - concordant - reversed_

    return {
        "concordant_hits": int(concordant),
        "reversed_hits": int(reversed_),
        "not_significant": int(not_sig),
        "total": int(total),
        "hit_rate_concordant": concordant / total,
        "hit_rate_reversed": reversed_ / total,
    }


def mean_signed_log2fc(
    signature_df: pd.DataFrame,
    de_df: pd.DataFrame,
) -> dict:
    """Mean of sign_i * log2FC_weissberg_i for detected signature genes.

    sign_i = +1 if reference direction is up, -1 if down.

    Returns:
        dict with keys: score, n_detected, n_total.
    """
    m = _merge_signature_with_de(signature_df, de_df)
    detected = m.dropna(subset=["log2fc"])

    if len(detected) == 0:
        return {"score": np.nan, "n_detected": 0, "n_total": len(m)}

    sign = detected["direction"].map({"up": 1, "down": -1})
    score = (sign * detected["log2fc"]).mean()

    return {"score": float(score), "n_detected": len(detected), "n_total": len(m)}


def mean_signed_normalized_rank(
    signature_df: pd.DataFrame,
    de_df: pd.DataFrame,
    total_genes: int,
) -> dict:
    """Mean of concordance_i * normalized_rank_i for detected signature genes.

    concordance_i = +1 if Weissberg direction matches reference, -1 if opposite, 0 if not significant.
    normalized_rank_i = 1 - (directional_rank_i / total_genes).

    Returns:
        dict with keys: score, n_concordant, n_reversed, n_not_significant, n_total.
    """
    m = _merge_signature_with_de(signature_df, de_df)
    total = len(m)

    if total == 0:
        return {
            "score": np.nan, "n_concordant": 0, "n_reversed": 0,
            "n_not_significant": 0, "n_total": 0,
        }

    m["de_direction"] = m["expression_status"].map({
        "significant_up": "up", "significant_down": "down",
    })

    # Directional rank: rank_up if significant_up, rank_down if significant_down
    m["dir_rank"] = np.where(
        m["expression_status"] == "significant_up", m["rank_up"],
        np.where(m["expression_status"] == "significant_down", m["rank_down"], np.nan)
    )

    m["concordance"] = np.where(
        m["de_direction"] == m["direction"], 1,
        np.where(m["de_direction"].notna() & (m["de_direction"] != m["direction"]), -1, 0)
    )

    m["normalized_rank"] = np.where(
        m["dir_rank"].notna(),
        1 - (m["dir_rank"].astype(float) / total_genes),
        0,
    )

    m["contribution"] = m["concordance"] * m["normalized_rank"]
    score = m["contribution"].mean()

    n_conc = (m["concordance"] == 1).sum()
    n_rev = (m["concordance"] == -1).sum()
    n_ns = (m["concordance"] == 0).sum()

    return {
        "score": float(score),
        "n_concordant": int(n_conc),
        "n_reversed": int(n_rev),
        "n_not_significant": int(n_ns),
        "n_total": int(total),
    }


def rank_correlation(
    signature_df: pd.DataFrame,
    de_ref_df: pd.DataFrame,
    de_target_df: pd.DataFrame,
) -> dict:
    """Spearman rho between reference and target directional ranks.

    Only includes genes significant in the expected direction in BOTH
    reference and target.

    Args:
        signature_df: Signature with locus_tag, direction.
        de_ref_df: Reference DE data (one timepoint).
        de_target_df: Target DE data (one timepoint).

    Returns:
        dict with keys: rho, pvalue, n_genes.
    """
    ref = _merge_signature_with_de(signature_df, de_ref_df).rename(
        columns={"rank_up": "ref_rank_up", "rank_down": "ref_rank_down",
                 "expression_status": "ref_status"}
    )
    target = _merge_signature_with_de(signature_df, de_target_df).rename(
        columns={"rank_up": "tgt_rank_up", "rank_down": "tgt_rank_down",
                 "expression_status": "tgt_status"}
    )

    merged = ref[["locus_tag", "direction", "ref_rank_up", "ref_rank_down", "ref_status"]].merge(
        target[["locus_tag", "tgt_rank_up", "tgt_rank_down", "tgt_status"]],
        on="locus_tag",
    )

    # Get directional ranks for genes concordant in both
    def get_dir_rank(row, prefix):
        if row["direction"] == "up":
            return row[f"{prefix}_rank_up"]
        else:
            return row[f"{prefix}_rank_down"]

    merged["ref_dir_rank"] = merged.apply(lambda r: get_dir_rank(r, "ref"), axis=1)
    merged["tgt_dir_rank"] = merged.apply(lambda r: get_dir_rank(r, "tgt"), axis=1)

    # Filter: significant in expected direction in both
    ref_concordant = merged.apply(
        lambda r: (r["direction"] == "up" and r["ref_status"] == "significant_up")
        or (r["direction"] == "down" and r["ref_status"] == "significant_down"),
        axis=1,
    )
    tgt_concordant = merged.apply(
        lambda r: (r["direction"] == "up" and r["tgt_status"] == "significant_up")
        or (r["direction"] == "down" and r["tgt_status"] == "significant_down"),
        axis=1,
    )

    valid = merged[ref_concordant & tgt_concordant].dropna(
        subset=["ref_dir_rank", "tgt_dir_rank"]
    )

    if len(valid) < 3:
        return {"rho": np.nan, "pvalue": np.nan, "n_genes": len(valid)}

    rho, pval = stats.spearmanr(valid["ref_dir_rank"], valid["tgt_dir_rank"])

    return {"rho": float(rho), "pvalue": float(pval), "n_genes": len(valid)}


def permutation_test_mean_signed_log2fc(
    signature_df: pd.DataFrame,
    de_df: pd.DataFrame,
    all_genes_log2fc: pd.Series,
    n_permutations: int = 1000,
    seed: int = 42,
) -> dict:
    """Permutation test for mean signed log2FC.

    Shuffles gene labels n_permutations times, recomputes mean signed
    log2FC each time, returns empirical p-value.

    Args:
        signature_df: Core/extended signature.
        de_df: DE data for one condition x timepoint.
        all_genes_log2fc: Series of log2fc for all genes in the experiment
            (background), indexed by locus_tag.
        n_permutations: Number of permutations.
        seed: Random seed for reproducibility.

    Returns:
        dict with keys: observed, empirical_p, n_permutations, n_signature_genes.
    """
    observed_result = mean_signed_log2fc(signature_df, de_df)
    observed = observed_result["score"]
    n_sig = len(signature_df)

    if np.isnan(observed) or n_sig < 30:
        return {
            "observed": observed,
            "empirical_p": np.nan,
            "n_permutations": 0,
            "n_signature_genes": n_sig,
        }

    rng = np.random.default_rng(seed)
    all_tags = all_genes_log2fc.index.tolist()
    directions = signature_df["direction"].values

    null_scores = np.empty(n_permutations)
    for i in range(n_permutations):
        random_tags = rng.choice(all_tags, size=n_sig, replace=False)
        random_fc = all_genes_log2fc.reindex(random_tags).values
        signs = np.where(directions == "up", 1, -1)
        null_scores[i] = np.nanmean(signs * random_fc)

    empirical_p = (np.abs(null_scores) >= np.abs(observed)).mean()

    return {
        "observed": float(observed),
        "empirical_p": float(empirical_p),
        "n_permutations": n_permutations,
        "n_signature_genes": n_sig,
    }
```

- [ ] **Step 2: Verify with a smoke test**

```bash
uv run python -c "
import sys; sys.path.insert(0, 'analyses/2026-04-06-1432-n_limitation_signature')
import pandas as pd
import numpy as np
from sig_utils.metrics import hit_rate, mean_signed_log2fc, mean_signed_normalized_rank

sig = pd.DataFrame({
    'locus_tag': ['PMM0920', 'PMM0246', 'PMM0100'],
    'direction': ['up', 'up', 'down'],
})
de = pd.DataFrame({
    'locus_tag': ['PMM0920', 'PMM0246', 'PMM0100'],
    'log2fc': [3.0, 2.0, -1.5],
    'expression_status': ['significant_up', 'significant_up', 'significant_down'],
    'rank': [5, 10, 20],
    'rank_up': [3, 5, None],
    'rank_down': [None, None, 8],
})

hr = hit_rate(sig, de)
print('hit_rate:', hr)
assert hr['concordant_hits'] == 3
assert hr['reversed_hits'] == 0

ms = mean_signed_log2fc(sig, de)
print('mean_signed_log2fc:', ms)
# (+1*3.0 + +1*2.0 + -1*-1.5) / 3 = (3+2+1.5)/3 = 2.167
assert abs(ms['score'] - 2.1667) < 0.01

mr = mean_signed_normalized_rank(sig, de, total_genes=1849)
print('mean_signed_normalized_rank:', mr)
assert mr['n_concordant'] == 3
assert mr['score'] > 0

print('OK')
"
```

Expected: prints metrics and `OK`.

- [ ] **Step 3: Commit**

```bash
git add analyses/2026-04-06-1432-n_limitation_signature/sig_utils/metrics.py
git commit -m "feat: sig_utils/metrics.py — hit_rate, mean_signed_log2fc, rank_score, rank_correlation"
```

---

### Task 5: Write `extract_reference_de.py`

**Files:**
- Create: `analyses/2026-04-06-1432-n_limitation_signature/scripts/extract_reference_de.py`

- [ ] **Step 1: Write the extraction script**

Write to `analyses/2026-04-06-1432-n_limitation_signature/scripts/extract_reference_de.py`:

```python
"""Extract DE data from reference nitrogen experiments (Tolonen 2006, Read 2017).

Outputs:
    data/de_reference_tolonen_ndep.csv — all DE rows for Tolonen N-deprivation
    data/de_reference_read_ndep.csv — all DE rows for Read N-depleted
    data/de_reference_tolonen_cyanate.csv — Tolonen cyanate (supplementary)
    data/de_reference_tolonen_urea.csv — Tolonen urea (supplementary)

Run from multiomics_research root:
    uv run analyses/2026-04-06-1432-n_limitation_signature/scripts/extract_reference_de.py
"""

import sys
from pathlib import Path

# Add analysis dir to path for sig_utils
ANALYSIS_DIR = Path(__file__).parent.parent
sys.path.insert(0, str(ANALYSIS_DIR))

from multiomics_explorer import differential_expression_by_gene
from multiomics_explorer.analysis import to_dataframe
from sig_utils.io import save_data_csv

# Experiment IDs
TOLONEN_NDEP = "10.1038/msb4100087_nitrogen_nitrogen_deprivation_med4_med4_microarray"
TOLONEN_CYANATE = "10.1038/msb4100087_growth_medium_growth_on_cyanate_as_med4_microarray"
TOLONEN_UREA = "10.1038/msb4100087_growth_medium_growth_on_urea_as_med4_microarray"
READ_NDEP = "10.1038/ismej.2017.88_nitrogen_stress_ndepleted_pro99_medium_med4_rnaseq"


def extract_experiment(experiment_id: str, label: str) -> None:
    """Extract all DE rows for an experiment and save to CSV."""
    print(f"Extracting {label}...")

    # First check total size
    result = differential_expression_by_gene(
        experiment_ids=[experiment_id],
        organism="MED4",
        summary=True,
    )
    total = result["total_matching"]
    print(f"  Total rows: {total}")

    # Extract all rows
    result = differential_expression_by_gene(
        experiment_ids=[experiment_id],
        organism="MED4",
        verbose=True,
        limit=None,
    )
    assert not result["truncated"], f"Results truncated for {label}!"

    df = to_dataframe(result)
    print(f"  Extracted {len(df)} rows, {df['locus_tag'].nunique()} genes")
    print(f"  Timepoints: {sorted(df['timepoint'].unique())}")
    print(f"  Status: {df['expression_status'].value_counts().to_dict()}")

    return df


def main():
    # Tolonen N-deprivation (main reference)
    df_tolonen = extract_experiment(TOLONEN_NDEP, "Tolonen 2006 N-deprivation")
    save_data_csv(df_tolonen, "de_reference_tolonen_ndep.csv")

    # Read N-depleted (second reference)
    df_read = extract_experiment(READ_NDEP, "Read 2017 N-depleted")
    save_data_csv(df_read, "de_reference_read_ndep.csv")

    # Supplementary: Tolonen cyanate and urea
    df_cyan = extract_experiment(TOLONEN_CYANATE, "Tolonen 2006 cyanate")
    save_data_csv(df_cyan, "de_reference_tolonen_cyanate.csv")

    df_urea = extract_experiment(TOLONEN_UREA, "Tolonen 2006 urea")
    save_data_csv(df_urea, "de_reference_tolonen_urea.csv")

    print("\nDone. Files saved to data/")


if __name__ == "__main__":
    main()
```

- [ ] **Step 2: Run the extraction**

```bash
uv run analyses/2026-04-06-1432-n_limitation_signature/scripts/extract_reference_de.py
```

Expected: prints row counts and timepoints for each experiment, saves 4 CSVs to `data/`.

- [ ] **Step 3: Verify outputs**

```bash
uv run python -c "
import pandas as pd
for f in ['de_reference_tolonen_ndep.csv', 'de_reference_read_ndep.csv']:
    df = pd.read_csv(f'analyses/2026-04-06-1432-n_limitation_signature/data/{f}')
    print(f'{f}: {len(df)} rows, {df[\"locus_tag\"].nunique()} genes')
    sig = df[df['expression_status'].isin(['significant_up','significant_down'])]
    print(f'  Significant: {len(sig)} rows, {sig[\"locus_tag\"].nunique()} unique genes')
"
```

- [ ] **Step 4: Commit**

```bash
git add analyses/2026-04-06-1432-n_limitation_signature/scripts/extract_reference_de.py
git add analyses/2026-04-06-1432-n_limitation_signature/data/de_reference_*.csv
git commit -m "feat: extract reference DE data (Tolonen 2006, Read 2017)"
```

---

### Task 6: Write `extract_weissberg_de.py`

**Files:**
- Create: `analyses/2026-04-06-1432-n_limitation_signature/scripts/extract_weissberg_de.py`

- [ ] **Step 1: Write the extraction script**

Write to `analyses/2026-04-06-1432-n_limitation_signature/scripts/extract_weissberg_de.py`:

```python
"""Extract DE data from Weissberg 2025 MED4 experiments.

Outputs:
    data/de_weissberg_med4_rnaseq_axenic.csv
    data/de_weissberg_med4_rnaseq_coculture.csv
    data/de_weissberg_med4_proteomics_axenic.csv
    data/de_weissberg_med4_proteomics_coculture.csv
    data/de_weissberg_med4_all.csv — combined for convenience

Run from multiomics_research root:
    uv run analyses/2026-04-06-1432-n_limitation_signature/scripts/extract_weissberg_de.py
"""

import sys
from pathlib import Path

ANALYSIS_DIR = Path(__file__).parent.parent
sys.path.insert(0, str(ANALYSIS_DIR))

import pandas as pd
from multiomics_explorer import differential_expression_by_gene
from multiomics_explorer.analysis import to_dataframe
from sig_utils.io import save_data_csv

WEISSBERG_DOI = "10.1101/2025.11.24.690089"

EXPERIMENTS = {
    "de_weissberg_med4_rnaseq_axenic.csv": (
        f"{WEISSBERG_DOI}_growth_state_pro99lown_nutrient_starvation_med4_rnaseq_axenic"
    ),
    "de_weissberg_med4_rnaseq_coculture.csv": (
        f"{WEISSBERG_DOI}_growth_state_pro99lown_nutrient_starvation_med4_rnaseq_coculture"
    ),
    "de_weissberg_med4_proteomics_axenic.csv": (
        f"{WEISSBERG_DOI}_growth_state_pro99lown_nutrient_starvation_med4_proteomics_axenic"
    ),
    "de_weissberg_med4_proteomics_coculture.csv": (
        f"{WEISSBERG_DOI}_growth_state_pro99lown_nutrient_starvation_med4_proteomics_coculture"
    ),
}


def main():
    all_dfs = []

    for filename, experiment_id in EXPERIMENTS.items():
        print(f"Extracting {filename}...")

        result = differential_expression_by_gene(
            experiment_ids=[experiment_id],
            organism="MED4",
            verbose=True,
            limit=None,
        )
        assert not result["truncated"], f"Results truncated for {filename}!"

        df = to_dataframe(result)
        print(f"  {len(df)} rows, {df['locus_tag'].nunique()} genes")
        if "timepoint" in df.columns:
            print(f"  Timepoints: {sorted(df['timepoint'].dropna().unique())}")

        save_data_csv(df, filename)
        df["source_file"] = filename
        all_dfs.append(df)

    # Combined file
    combined = pd.concat(all_dfs, ignore_index=True)
    save_data_csv(combined, "de_weissberg_med4_all.csv")
    print(f"\nCombined: {len(combined)} rows")
    print("Done. Files saved to data/")


if __name__ == "__main__":
    main()
```

- [ ] **Step 2: Run the extraction**

```bash
uv run analyses/2026-04-06-1432-n_limitation_signature/scripts/extract_weissberg_de.py
```

- [ ] **Step 3: Verify outputs**

```bash
uv run python -c "
import pandas as pd
for f in ['de_weissberg_med4_rnaseq_axenic.csv','de_weissberg_med4_rnaseq_coculture.csv',
          'de_weissberg_med4_proteomics_axenic.csv','de_weissberg_med4_proteomics_coculture.csv']:
    df = pd.read_csv(f'analyses/2026-04-06-1432-n_limitation_signature/data/{f}')
    print(f'{f}: {len(df)} rows, {df[\"locus_tag\"].nunique()} genes')
"
```

- [ ] **Step 4: Commit**

```bash
git add analyses/2026-04-06-1432-n_limitation_signature/scripts/extract_weissberg_de.py
git add analyses/2026-04-06-1432-n_limitation_signature/data/de_weissberg_*.csv
git commit -m "feat: extract Weissberg 2025 MED4 DE data (RNA-seq + proteomics)"
```

---

### Task 7: Write `build_signature.py`

**Files:**
- Create: `analyses/2026-04-06-1432-n_limitation_signature/scripts/build_signature.py`

- [ ] **Step 1: Write the signature building script**

Write to `analyses/2026-04-06-1432-n_limitation_signature/scripts/build_signature.py`:

```python
"""Build core and extended N-limitation signatures from reference DE data.

Reads: data/de_reference_tolonen_ndep.csv, data/de_reference_read_ndep.csv
Outputs: data/core_signature_genes.csv, data/extended_signature_genes.csv

Run from multiomics_research root:
    uv run analyses/2026-04-06-1432-n_limitation_signature/scripts/build_signature.py
"""

import sys
from pathlib import Path

ANALYSIS_DIR = Path(__file__).parent.parent
sys.path.insert(0, str(ANALYSIS_DIR))

from sig_utils.io import load_de_csv, save_data_csv
from sig_utils.signature import summarize_de_per_gene, intersect_de_lists


def main():
    # Load reference DE data
    tolonen = load_de_csv("de_reference_tolonen_ndep.csv")
    read = load_de_csv("de_reference_read_ndep.csv")

    print(f"Tolonen: {len(tolonen)} rows, {tolonen['locus_tag'].nunique()} genes")
    print(f"Read: {len(read)} rows, {read['locus_tag'].nunique()} genes")

    # Filter Tolonen to informative timepoints (skip 0h, 3h)
    tolonen_filtered = tolonen[~tolonen["timepoint"].isin(["0h", "3h"])]
    print(f"Tolonen after filtering 0h/3h: {len(tolonen_filtered)} rows")

    # Summarize to one row per gene (peak timepoint)
    tolonen_summary = summarize_de_per_gene(tolonen_filtered)
    read_summary = summarize_de_per_gene(read)

    print(f"\nTolonen significant genes: {len(tolonen_summary)}")
    print(f"  Up: {(tolonen_summary['direction'] == 'up').sum()}")
    print(f"  Down: {(tolonen_summary['direction'] == 'down').sum()}")

    print(f"Read significant genes: {len(read_summary)}")
    print(f"  Up: {(read_summary['direction'] == 'up').sum()}")
    print(f"  Down: {(read_summary['direction'] == 'down').sum()}")

    # Get all locus tags in Read dataset (not just significant) for tagging
    read_all_tags = set(read["locus_tag"].unique())
    print(f"Read total locus tags (including non-significant): {len(read_all_tags)}")

    # Build signatures
    core, extended = intersect_de_lists(
        tolonen_summary, read_summary,
        study_b_all_locus_tags=read_all_tags,
    )

    print(f"\n=== CORE SIGNATURE ===")
    print(f"Total genes: {len(core)}")
    print(f"  Up: {(core['direction'] == 'up').sum()}")
    print(f"  Down: {(core['direction'] == 'down').sum()}")

    if len(core) < 30:
        print("  WARNING: Core signature < 30 genes. Permutation test will be underpowered.")
        print("  Extended signature results should be given more weight.")

    print(f"\n=== EXTENDED SIGNATURE ===")
    print(f"Total genes: {len(extended)}")
    print(f"  By type: {extended['signature_type'].value_counts().to_dict()}")

    # Save
    save_data_csv(core, "core_signature_genes.csv")
    save_data_csv(extended, "extended_signature_genes.csv")

    print("\nDone. Saved to data/core_signature_genes.csv and data/extended_signature_genes.csv")


if __name__ == "__main__":
    main()
```

- [ ] **Step 2: Run the script**

```bash
uv run analyses/2026-04-06-1432-n_limitation_signature/scripts/build_signature.py
```

Expected: prints gene counts for core and extended signatures.

- [ ] **Step 3: Inspect outputs**

```bash
uv run python -c "
import pandas as pd
core = pd.read_csv('analyses/2026-04-06-1432-n_limitation_signature/data/core_signature_genes.csv')
ext = pd.read_csv('analyses/2026-04-06-1432-n_limitation_signature/data/extended_signature_genes.csv')
print('Core signature:')
print(core[['locus_tag','gene_name','direction','cross_study_best_dir_rank']].head(20).to_string())
print(f'\nTotal core: {len(core)}, Extended: {len(ext)}')
print(f'Extended types: {ext[\"signature_type\"].value_counts().to_dict()}')
"
```

- [ ] **Step 4: Start exploration log**

Write to `analyses/2026-04-06-1432-n_limitation_signature/exploration/2026-04-06-signature-building.md`:

```markdown
# 2026-04-06: Signature building

## Question
How many MED4 genes form a reproducible N-limitation signature across
Tolonen 2006 (microarray) and Read 2017 (RNA-seq)?
Type: exploratory

## Approach
Extract all DE genes from both studies, summarize per gene (peak timepoint),
intersect with concordant direction requirement.

## Findings
(Fill in after running build_signature.py — record actual counts)

- [KG] Tolonen N-deprivation: X significant genes across 6-48h timepoints
- [KG] Read N-depleted: X significant genes across 3-24h timepoints
- [KG] Core signature (concordant intersection): X genes (Y up, Z down)
- [KG] Extended signature: X genes (breakdown by type)

## Assessment
(Fill in: established / preliminary / speculative)

## Gaps and friction
(Fill in any issues encountered)

## Next
Score core and extended signatures against Weissberg 2025 experiments.
```

- [ ] **Step 5: Commit**

```bash
git add analyses/2026-04-06-1432-n_limitation_signature/scripts/build_signature.py
git add analyses/2026-04-06-1432-n_limitation_signature/data/core_signature_genes.csv
git add analyses/2026-04-06-1432-n_limitation_signature/data/extended_signature_genes.csv
git add analyses/2026-04-06-1432-n_limitation_signature/exploration/
git commit -m "feat: build core + extended N-limitation signatures from Tolonen/Read"
```

---

### Task 8: Write `score_signature.py`

**Files:**
- Create: `analyses/2026-04-06-1432-n_limitation_signature/scripts/score_signature.py`

- [ ] **Step 1: Write the scoring script**

Write to `analyses/2026-04-06-1432-n_limitation_signature/scripts/score_signature.py`:

```python
"""Score Weissberg and reference experiments against N-limitation signatures.

Reads: data/core_signature_genes.csv, data/extended_signature_genes.csv,
       data/de_weissberg_med4_*.csv, data/de_reference_*.csv
Outputs: results/signature_scores_core.csv, results/signature_scores_extended.csv,
         results/reference_baseline_scores.csv, results/core_vs_extended_comparison.csv

Run from multiomics_research root:
    uv run analyses/2026-04-06-1432-n_limitation_signature/scripts/score_signature.py
"""

import sys
from pathlib import Path

ANALYSIS_DIR = Path(__file__).parent.parent
sys.path.insert(0, str(ANALYSIS_DIR))

import pandas as pd
import numpy as np
from sig_utils.io import load_de_csv, load_signature_csv, save_scores_csv
from sig_utils.metrics import (
    hit_rate,
    mean_signed_log2fc,
    mean_signed_normalized_rank,
    permutation_test_mean_signed_log2fc,
)

# Experiment files and their metadata
WEISSBERG_EXPERIMENTS = [
    ("de_weissberg_med4_rnaseq_axenic.csv", "rnaseq", "axenic"),
    ("de_weissberg_med4_rnaseq_coculture.csv", "rnaseq", "coculture"),
    ("de_weissberg_med4_proteomics_axenic.csv", "proteomics", "axenic"),
    ("de_weissberg_med4_proteomics_coculture.csv", "proteomics", "coculture"),
]

REFERENCE_EXPERIMENTS = [
    ("de_reference_tolonen_ndep.csv", "microarray", "tolonen_ndep"),
    ("de_reference_read_ndep.csv", "rnaseq", "read_ndep"),
]

# Timepoints to skip
SKIP_TIMEPOINTS_TOLONEN = {"0h", "3h"}


def score_one_timepoint(
    signature_df: pd.DataFrame,
    de_timepoint_df: pd.DataFrame,
    total_genes: int,
) -> dict:
    """Compute all metrics for one signature x one timepoint."""
    hr = hit_rate(signature_df, de_timepoint_df)
    ms = mean_signed_log2fc(signature_df, de_timepoint_df)
    mr = mean_signed_normalized_rank(signature_df, de_timepoint_df, total_genes)

    return {
        "hit_rate_concordant": hr["hit_rate_concordant"],
        "hit_rate_reversed": hr["hit_rate_reversed"],
        "concordant_hits": hr["concordant_hits"],
        "reversed_hits": hr["reversed_hits"],
        "not_significant": hr["not_significant"],
        "mean_signed_log2fc": ms["score"],
        "n_detected": ms["n_detected"],
        "rank_score": mr["score"],
        "n_concordant_rank": mr["n_concordant"],
        "n_reversed_rank": mr["n_reversed"],
    }


def score_experiment(
    signature_df: pd.DataFrame,
    de_df: pd.DataFrame,
    platform: str,
    condition: str,
    study: str,
    skip_timepoints: set | None = None,
) -> list[dict]:
    """Score all timepoints in an experiment."""
    rows = []

    if "timepoint" in de_df.columns and de_df["timepoint"].nunique() > 1:
        for tp, tp_df in de_df.groupby("timepoint"):
            if skip_timepoints and tp in skip_timepoints:
                continue
            total_genes = tp_df["locus_tag"].nunique()
            scores = score_one_timepoint(signature_df, tp_df, total_genes)
            scores.update({
                "study": study,
                "platform": platform,
                "condition": condition,
                "timepoint": tp,
                "timepoint_hours": tp_df["timepoint_hours"].iloc[0] if "timepoint_hours" in tp_df.columns else np.nan,
                "total_genes_in_experiment": total_genes,
                "signature_genes_total": len(signature_df),
            })
            rows.append(scores)
    else:
        total_genes = de_df["locus_tag"].nunique()
        scores = score_one_timepoint(signature_df, de_df, total_genes)
        scores.update({
            "study": study,
            "platform": platform,
            "condition": condition,
            "timepoint": "single",
            "timepoint_hours": np.nan,
            "total_genes_in_experiment": total_genes,
            "signature_genes_total": len(signature_df),
        })
        rows.append(scores)

    return rows


def main():
    # Load signatures
    core = load_signature_csv("core_signature_genes.csv")
    extended = load_signature_csv("extended_signature_genes.csv")
    print(f"Core signature: {len(core)} genes")
    print(f"Extended signature: {len(extended)} genes")

    for sig_name, sig_df, output_file in [
        ("core", core, "signature_scores_core.csv"),
        ("extended", extended, "signature_scores_extended.csv"),
    ]:
        print(f"\n=== Scoring {sig_name} signature ===")
        all_rows = []

        # Score Weissberg experiments
        for filename, platform, condition in WEISSBERG_EXPERIMENTS:
            print(f"  {filename}...")
            de = load_de_csv(filename)
            rows = score_experiment(sig_df, de, platform, condition, "weissberg_2025")
            all_rows.extend(rows)

        # Score reference experiments (baselines)
        for filename, platform, study in REFERENCE_EXPERIMENTS:
            print(f"  {filename} (baseline)...")
            de = load_de_csv(filename)
            skip = SKIP_TIMEPOINTS_TOLONEN if "tolonen" in filename else None
            rows = score_experiment(sig_df, de, platform, "reference", study, skip)
            all_rows.extend(rows)

        scores_df = pd.DataFrame(all_rows)
        save_scores_csv(scores_df, output_file)
        print(f"  Saved {len(scores_df)} rows to results/{output_file}")

    # Permutation tests for core signature (Weissberg only)
    print("\n=== Permutation tests (core, mean_signed_log2fc) ===")
    perm_rows = []
    for filename, platform, condition in WEISSBERG_EXPERIMENTS:
        de = load_de_csv(filename)
        all_fc = de.set_index("locus_tag")["log2fc"]

        if "timepoint" in de.columns and de["timepoint"].nunique() > 1:
            for tp, tp_df in de.groupby("timepoint"):
                tp_fc = tp_df.set_index("locus_tag")["log2fc"]
                result = permutation_test_mean_signed_log2fc(core, tp_df, tp_fc)
                result.update({"platform": platform, "condition": condition, "timepoint": tp})
                perm_rows.append(result)
        else:
            result = permutation_test_mean_signed_log2fc(core, de, all_fc)
            result.update({"platform": platform, "condition": condition, "timepoint": "single"})
            perm_rows.append(result)

    perm_df = pd.DataFrame(perm_rows)
    save_scores_csv(perm_df, "permutation_tests.csv")
    print(f"  Saved {len(perm_df)} permutation tests")

    # Core vs extended comparison
    print("\n=== Core vs extended comparison ===")
    core_scores = pd.read_csv(ANALYSIS_DIR / "results" / "signature_scores_core.csv")
    ext_scores = pd.read_csv(ANALYSIS_DIR / "results" / "signature_scores_extended.csv")

    core_scores = core_scores.rename(columns={
        "hit_rate_concordant": "core_hit_rate",
        "mean_signed_log2fc": "core_mean_signed_log2fc",
        "rank_score": "core_rank_score",
    })
    ext_scores = ext_scores.rename(columns={
        "hit_rate_concordant": "extended_hit_rate",
        "mean_signed_log2fc": "extended_mean_signed_log2fc",
        "rank_score": "extended_rank_score",
    })

    merge_keys = ["study", "platform", "condition", "timepoint"]
    comparison = core_scores[merge_keys + ["core_hit_rate", "core_mean_signed_log2fc", "core_rank_score"]].merge(
        ext_scores[merge_keys + ["extended_hit_rate", "extended_mean_signed_log2fc", "extended_rank_score"]],
        on=merge_keys,
    )
    save_scores_csv(comparison, "core_vs_extended_comparison.csv")
    print(f"  Saved {len(comparison)} comparison rows")

    print("\nDone.")


if __name__ == "__main__":
    main()
```

- [ ] **Step 2: Run the scoring**

```bash
uv run analyses/2026-04-06-1432-n_limitation_signature/scripts/score_signature.py
```

- [ ] **Step 3: Verify outputs**

```bash
uv run python -c "
import pandas as pd
for f in ['signature_scores_core.csv', 'signature_scores_extended.csv',
          'reference_baseline_scores.csv', 'core_vs_extended_comparison.csv',
          'permutation_tests.csv']:
    path = f'analyses/2026-04-06-1432-n_limitation_signature/results/{f}'
    try:
        df = pd.read_csv(path)
        print(f'{f}: {len(df)} rows')
    except FileNotFoundError:
        print(f'{f}: NOT FOUND')
"
```

- [ ] **Step 4: Commit**

```bash
git add analyses/2026-04-06-1432-n_limitation_signature/scripts/score_signature.py
git add analyses/2026-04-06-1432-n_limitation_signature/results/
git commit -m "feat: score core + extended signatures against Weissberg and reference experiments"
```

---

### Task 9: Write `plot_trajectories.py`

**Files:**
- Create: `analyses/2026-04-06-1432-n_limitation_signature/scripts/plot_trajectories.py`

- [ ] **Step 1: Write the plotting script**

Write to `analyses/2026-04-06-1432-n_limitation_signature/scripts/plot_trajectories.py`:

```python
"""Plot activation score trajectories for axenic vs coculture.

Reads: results/signature_scores_core.csv, results/signature_scores_extended.csv
Outputs: results/trajectory_rnaseq.png, results/trajectory_proteomics.png

Run from multiomics_research root:
    uv run analyses/2026-04-06-1432-n_limitation_signature/scripts/plot_trajectories.py
"""

import sys
from pathlib import Path

ANALYSIS_DIR = Path(__file__).parent.parent
sys.path.insert(0, str(ANALYSIS_DIR))

import pandas as pd
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from sig_utils.io import RESULTS_DIR


def load_weissberg_scores(filename: str) -> pd.DataFrame:
    """Load scores and filter to Weissberg experiments with numeric timepoints."""
    df = pd.read_csv(RESULTS_DIR / filename)
    df = df[df["study"] == "weissberg_2025"].copy()
    # Exclude combined timepoints (d60+89)
    df = df[~df["timepoint"].str.contains(r"\+", na=False)]
    # Parse day number from timepoint strings like "day 18", "day 31"
    df["day"] = df["timepoint"].str.extract(r"(\d+)").astype(float)
    return df


def load_reference_scores(filename: str) -> pd.DataFrame:
    """Load reference baseline scores."""
    df = pd.read_csv(RESULTS_DIR / filename)
    df = df[df["study"] != "weissberg_2025"].copy()
    return df


METRICS = [
    ("hit_rate_concordant", "Hit rate (concordant)", 0, 1),
    ("mean_signed_log2fc", "Mean signed log2FC", None, None),
    ("rank_score", "Rank score", None, None),
]


def plot_trajectories_for_platform(
    platform: str,
    core_df: pd.DataFrame,
    extended_df: pd.DataFrame,
    ref_df: pd.DataFrame,
    output_path: Path,
):
    """Plot 3-panel figure (one per metric) for a platform."""
    core_plat = core_df[core_df["platform"] == platform]
    ext_plat = extended_df[extended_df["platform"] == platform]
    ref_plat = ref_df[ref_df["study"].str.contains("tolonen" if platform == "microarray" else "")]

    fig, axes = plt.subplots(1, 3, figsize=(18, 5))
    fig.suptitle(f"N-Limitation Signature Activation — {platform.upper()}", fontsize=14)

    for ax, (metric, label, ymin, ymax) in zip(axes, METRICS):
        # Weissberg core
        for condition, style, color in [("axenic", "--", "red"), ("coculture", "-", "blue")]:
            subset = core_plat[core_plat["condition"] == condition].sort_values("day")
            if not subset.empty:
                ax.plot(subset["day"], subset[metric], style, color=color,
                        marker="o", label=f"{condition} (core)", linewidth=2)

        # Weissberg extended (lighter)
        for condition, style, color in [("axenic", "--", "salmon"), ("coculture", "-", "lightblue")]:
            subset = ext_plat[ext_plat["condition"] == condition].sort_values("day")
            if not subset.empty:
                ax.plot(subset["day"], subset[metric], style, color=color,
                        marker="s", alpha=0.6, label=f"{condition} (extended)", linewidth=1)

        # Reference baselines as horizontal bands
        for study_label, ref_color in [("tolonen_ndep", "gray"), ("read_ndep", "darkgray")]:
            ref_study = ref_plat[ref_plat["study"] == study_label]
            if not ref_study.empty:
                ref_min = ref_study[metric].min()
                ref_max = ref_study[metric].max()
                ax.axhspan(ref_min, ref_max, alpha=0.15, color=ref_color,
                          label=f"{study_label} range")

        ax.set_xlabel("Day")
        ax.set_ylabel(label)
        if ymin is not None:
            ax.set_ylim(ymin, ymax)
        ax.legend(fontsize=8, loc="best")
        ax.grid(True, alpha=0.3)

    plt.tight_layout()
    plt.savefig(output_path, dpi=300, bbox_inches="tight")
    plt.close()
    print(f"Saved: {output_path}")


def main():
    core_w = load_weissberg_scores("signature_scores_core.csv")
    ext_w = load_weissberg_scores("signature_scores_extended.csv")
    core_all = pd.read_csv(RESULTS_DIR / "signature_scores_core.csv")
    ref = core_all[core_all["study"] != "weissberg_2025"]

    for platform in ["rnaseq", "proteomics"]:
        output = RESULTS_DIR / f"trajectory_{platform}.png"
        plot_trajectories_for_platform(platform, core_w, ext_w, ref, output)

    print("Done.")


if __name__ == "__main__":
    main()
```

- [ ] **Step 2: Run the plotting**

```bash
uv run analyses/2026-04-06-1432-n_limitation_signature/scripts/plot_trajectories.py
```

Expected: creates `trajectory_rnaseq.png` and `trajectory_proteomics.png` in results/.

- [ ] **Step 3: Commit**

```bash
git add analyses/2026-04-06-1432-n_limitation_signature/scripts/plot_trajectories.py
git add analyses/2026-04-06-1432-n_limitation_signature/results/trajectory_*.png
git commit -m "feat: trajectory plots for N-limitation signature scores"
```

---

### Task 10: Write `explore_alteromonas.py` and pathway summary

**Files:**
- Create: `analyses/2026-04-06-1432-n_limitation_signature/scripts/explore_alteromonas.py`

- [ ] **Step 1: Write Alteromonas exploration script**

Write to `analyses/2026-04-06-1432-n_limitation_signature/scripts/explore_alteromonas.py`:

```python
"""Targeted check: Alteromonas HOT1A3 N-recycling genes in Weissberg coculture.

Step 7 of Approach A: query DE status of N-recycling genes in HOT1A3
coculture time course experiments.

Outputs: data/de_weissberg_hot1a3_nrecycling.csv, results/alteromonas_nrecycling.csv

Run from multiomics_research root:
    uv run analyses/2026-04-06-1432-n_limitation_signature/scripts/explore_alteromonas.py
"""

import sys
from pathlib import Path

ANALYSIS_DIR = Path(__file__).parent.parent
sys.path.insert(0, str(ANALYSIS_DIR))

import pandas as pd
from multiomics_explorer import genes_by_function, differential_expression_by_gene
from multiomics_explorer.analysis import to_dataframe
from sig_utils.io import save_data_csv, save_scores_csv

WEISSBERG_DOI = "10.1101/2025.11.24.690089"

HOT1A3_COCULTURE_EXPERIMENTS = [
    f"{WEISSBERG_DOI}_growth_state_pro99lown_nutrient_starvation_hot1a3_rnaseq_coculture",
    f"{WEISSBERG_DOI}_growth_state_pro99lown_nutrient_starvation_hot1a3_proteomics_coculture",
]

# N-recycling search terms to find relevant HOT1A3 genes
NRECYCLING_SEARCHES = [
    "ammonium AND transport",
    "amino acid AND degradation",
    "peptidase OR protease",
    "deaminase",
    "urease",
    "glutamate AND dehydrogenase",
    "organic nitrogen",
]


def find_nrecycling_genes() -> list[str]:
    """Search for N-recycling genes in HOT1A3."""
    all_tags = set()

    for search in NRECYCLING_SEARCHES:
        result = genes_by_function(
            search_text=search,
            organism="HOT1A3",
            limit=None,
        )
        tags = {r["locus_tag"] for r in result["results"]}
        print(f"  '{search}': {len(tags)} genes")
        all_tags.update(tags)

    print(f"  Total unique N-recycling locus tags: {len(all_tags)}")
    return sorted(all_tags)


def main():
    print("Finding HOT1A3 N-recycling genes...")
    tags = find_nrecycling_genes()

    if not tags:
        print("No N-recycling genes found! Check search terms.")
        return

    print(f"\nExtracting DE for {len(tags)} genes across coculture experiments...")
    all_dfs = []

    for exp_id in HOT1A3_COCULTURE_EXPERIMENTS:
        result = differential_expression_by_gene(
            locus_tags=tags,
            experiment_ids=[exp_id],
            verbose=True,
            limit=None,
        )
        if result["results"]:
            df = to_dataframe(result)
            all_dfs.append(df)
            print(f"  {exp_id.split('_')[-1]}: {len(df)} rows")

    if not all_dfs:
        print("No DE data found for N-recycling genes.")
        return

    combined = pd.concat(all_dfs, ignore_index=True)
    save_data_csv(combined, "de_weissberg_hot1a3_nrecycling.csv")

    # Summary: which genes are upregulated at which timepoints?
    sig = combined[combined["expression_status"].isin(["significant_up", "significant_down"])]
    summary = sig.groupby(["locus_tag", "gene_name", "product", "timepoint", "omics_type"]).agg(
        log2fc=("log2fc", "first"),
        expression_status=("expression_status", "first"),
        rank=("rank", "first"),
        rank_up=("rank_up", "first"),
        rank_down=("rank_down", "first"),
    ).reset_index()

    save_scores_csv(summary, "alteromonas_nrecycling.csv")
    print(f"\nSignificant N-recycling genes: {summary['locus_tag'].nunique()}")
    print(f"Summary saved to results/alteromonas_nrecycling.csv")
    print("\nSignificant genes per timepoint:")
    print(summary.groupby(["timepoint", "expression_status"]).size().unstack(fill_value=0))


if __name__ == "__main__":
    main()
```

- [ ] **Step 2: Run the script**

```bash
uv run analyses/2026-04-06-1432-n_limitation_signature/scripts/explore_alteromonas.py
```

- [ ] **Step 3: Commit**

```bash
git add analyses/2026-04-06-1432-n_limitation_signature/scripts/explore_alteromonas.py
git add analyses/2026-04-06-1432-n_limitation_signature/data/de_weissberg_hot1a3_nrecycling.csv
git add analyses/2026-04-06-1432-n_limitation_signature/results/alteromonas_nrecycling.csv
git commit -m "feat: Alteromonas HOT1A3 N-recycling gene check"
```

---

### Task 11: Pathway annotation and methods.md update

**Files:**
- Modify: `analyses/2026-04-06-1432-n_limitation_signature/methods.md`
- Modify: `analyses/2026-04-06-1432-n_limitation_signature/exploration/2026-04-06-signature-building.md`

This task is manual/interpretive — it uses the outputs from Tasks 5-10.

- [ ] **Step 1: Generate pathway summary from core signature**

```bash
uv run python -c "
import sys; sys.path.insert(0, 'analyses/2026-04-06-1432-n_limitation_signature')
import pandas as pd
from sig_utils.io import load_signature_csv, RESULTS_DIR

core = load_signature_csv('core_signature_genes.csv')

# Use gene_category from the DE data (verbose output includes it)
de = pd.read_csv('analyses/2026-04-06-1432-n_limitation_signature/data/de_reference_tolonen_ndep.csv')
categories = de[['locus_tag', 'gene_category']].drop_duplicates()
core_cat = core.merge(categories, on='locus_tag', how='left')

summary = core_cat.groupby(['gene_category', 'direction']).agg(
    count=('locus_tag', 'nunique'),
    genes=('locus_tag', lambda x: ', '.join(sorted(x.unique()))),
    gene_names=('gene_name', lambda x: ', '.join(str(n) for n in sorted(x.dropna().unique()))),
).reset_index()

summary.to_csv(RESULTS_DIR / 'pathway_summary.csv', index=False)
print(summary.to_string())
"
```

- [ ] **Step 2: Update exploration log with actual findings**

Read the outputs from all previous scripts and fill in the `## Findings` and `## Assessment` sections of the exploration log with actual numbers. Tag each finding with `[KG]`, `[interpretation]`, or `[gap]`.

- [ ] **Step 3: Update methods.md with actual gene counts and thresholds**

Fill in the Gene selection and Statistical methods sections based on actual results.

- [ ] **Step 4: Commit**

```bash
git add analyses/2026-04-06-1432-n_limitation_signature/results/pathway_summary.csv
git add analyses/2026-04-06-1432-n_limitation_signature/exploration/
git add analyses/2026-04-06-1432-n_limitation_signature/methods.md
git commit -m "docs: pathway summary, exploration log, and methods update"
```

---

### Task 12: Final review and README update

- [ ] **Step 1: Review all outputs exist**

```bash
ls -la analyses/2026-04-06-1432-n_limitation_signature/data/
ls -la analyses/2026-04-06-1432-n_limitation_signature/results/
```

Expected files:
- `data/`: 8+ CSVs (reference DE, Weissberg DE, signatures)
- `results/`: 7+ files (scores, comparison, permutation, pathway, trajectories, Alteromonas)

- [ ] **Step 2: Update README with key findings**

Update `analyses/2026-04-06-1432-n_limitation_signature/README.md` with:
- Core signature size and composition
- Key trajectory observations (axenic vs coculture divergence)
- Alteromonas N-recycling findings

- [ ] **Step 3: Update gaps_and_friction.md with any issues encountered**

- [ ] **Step 4: Final commit**

```bash
git add analyses/2026-04-06-1432-n_limitation_signature/
git commit -m "docs: complete Approach A analysis — README, gaps, plan snapshot"
```
