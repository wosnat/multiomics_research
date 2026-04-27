"""Permutation test for every (axis × panel × experiment × TP) cell.

For each cell we have a fixed gene set with calibrated directions and an observed
axis_score. We permute *which genes count as "the axis"*: K random draws of N
genes from the genome (excluding the actual axis genes), each draw assigned the
same multiset of directions as the actual axis (shuffled). Empirical p-value =
fraction of permuted scores that match-or-exceed the observed (one-tailed in the
direction of the observed score; cells with score = 0 get p = 1.0).

Multiple-testing correction across all cells: Benjamini-Hochberg FDR.

Inputs:
  - 5_analyze/data/all_axes_scores.csv         (108 cell scores)
  - 5_analyze/data/panel_definitions.csv       (axis × panel × locus × direction)
  - 5_analyze/data/genome_de_<short_id>.csv    (4 files; pulled in step 5)

Outputs:
  - data/permutation_pvalues.csv               one row per cell, observed score + p_raw + p_bh
  - data/01_permutation_test.log               progress log

Run:
  uv run analyses/2026-04-27-1117-prochlorococcus_stress_axenic_vs_coculture/6_evaluate/scripts/01_permutation_test.py
"""

from __future__ import annotations

from pathlib import Path
import sys
import time

import numpy as np
import pandas as pd

ANALYSIS_DIR = Path(__file__).resolve().parents[2]

DATA_5 = ANALYSIS_DIR / "5_analyze" / "data"
SCORES_CSV = DATA_5 / "all_axes_scores.csv"
PANEL_CSV = DATA_5 / "panel_definitions.csv"

OUT_DIR = Path(__file__).resolve().parents[1] / "data"
OUT_DIR.mkdir(parents=True, exist_ok=True)
LOG_PATH = OUT_DIR / "01_permutation_test.log"

K_PERMUTATIONS = 10_000
SEED = 42


def log(msg: str, *, fh) -> None:
    print(msg, flush=True)
    fh.write(msg + "\n")
    fh.flush()


def load_genome_de(experiment_short: str) -> pd.DataFrame:
    path = DATA_5 / f"genome_de_{experiment_short}.csv"
    df = pd.read_csv(path)
    df["log2fc"] = pd.to_numeric(df["log2fc"], errors="coerce")
    return df


def axis_score_array(
    log2fc_arr: np.ndarray,        # shape (n_genes,) — log2FC of all genes at this TP
    axis_idx: np.ndarray,          # shape (N,) — indices into log2fc_arr for the axis
    direction_arr: np.ndarray,     # shape (N,) — ±1 directions
) -> float:
    """Vectorized axis_stress_score for a single (axis-set, log2FC vector) pair.

    Equivalent to stress_score.axis_stress_score but operates on numpy arrays
    for speed. NaN handling matches: NaN log2FC genes are excluded from both
    axis and background.
    """
    axis_lfc = log2fc_arr[axis_idx]
    axis_dir = direction_arr
    valid_axis = ~np.isnan(axis_lfc)
    if valid_axis.sum() < 1:
        return np.nan
    axis_signed = axis_lfc[valid_axis] * axis_dir[valid_axis]
    axis_mean = axis_signed.mean()

    bg_mask = np.ones_like(log2fc_arr, dtype=bool)
    bg_mask[axis_idx] = False
    bg_lfc = log2fc_arr[bg_mask]
    bg_lfc = bg_lfc[~np.isnan(bg_lfc)]
    if len(bg_lfc) < 2:
        return np.nan

    bg_mean = bg_lfc.mean()
    bg_sd = bg_lfc.std(ddof=0)
    if bg_sd == 0:
        return np.nan

    return (axis_mean - bg_mean) / bg_sd


def permutation_pvalue_for_cell(
    log2fc_arr: np.ndarray,        # all genes' log2FC at this TP
    actual_axis_idx: np.ndarray,   # actual axis gene indices
    direction_arr: np.ndarray,     # actual directions for those axis genes
    observed_score: float,
    rng: np.random.Generator,
    K: int = K_PERMUTATIONS,
) -> tuple[float, int, float, float]:
    """Return (p_value, n_exceed, perm_mean, perm_sd).

    One-tailed test in the direction of the observed score.
    """
    if np.isnan(observed_score):
        return float("nan"), 0, float("nan"), float("nan")

    n_axis = len(actual_axis_idx)
    n_total = len(log2fc_arr)
    eligible_mask = np.ones(n_total, dtype=bool)
    eligible_mask[actual_axis_idx] = False
    eligible_idx = np.where(eligible_mask)[0]

    perm_scores = np.empty(K, dtype=np.float64)
    for k in range(K):
        # Draw N random genes from non-axis pool
        random_axis_idx = rng.choice(eligible_idx, size=n_axis, replace=False)
        # Shuffle the direction multiset and assign
        permuted_directions = direction_arr.copy()
        rng.shuffle(permuted_directions)
        perm_scores[k] = axis_score_array(log2fc_arr, random_axis_idx, permuted_directions)

    perm_scores_finite = perm_scores[~np.isnan(perm_scores)]
    if observed_score >= 0:
        n_exceed = int((perm_scores_finite >= observed_score).sum())
    else:
        n_exceed = int((perm_scores_finite <= observed_score).sum())
    K_eff = len(perm_scores_finite)
    p_value = (n_exceed + 1) / (K_eff + 1)

    return (
        p_value,
        n_exceed,
        float(perm_scores_finite.mean()) if K_eff else float("nan"),
        float(perm_scores_finite.std(ddof=0)) if K_eff else float("nan"),
    )


def bh_fdr(pvalues: np.ndarray) -> np.ndarray:
    """Benjamini-Hochberg FDR. Returns adjusted p-values."""
    pvalues = np.asarray(pvalues, dtype=float)
    n = len(pvalues)
    order = np.argsort(pvalues)
    ranked = np.arange(1, n + 1)
    adj = pvalues[order] * n / ranked
    # Enforce monotonicity from the largest down
    adj = np.minimum.accumulate(adj[::-1])[::-1]
    out = np.empty(n, dtype=float)
    out[order] = np.minimum(adj, 1.0)
    return out


def main() -> None:
    with LOG_PATH.open("w") as fh:
        rng = np.random.default_rng(SEED)
        scores = pd.read_csv(SCORES_CSV)
        panels = pd.read_csv(PANEL_CSV)

        log(f"Permutation test: K = {K_PERMUTATIONS}, seed = {SEED}", fh=fh)
        log(f"Scoring {len(scores)} cells", fh=fh)

        # Cache genome DE per experiment_short
        genome_cache: dict[str, pd.DataFrame] = {}
        for s in scores["experiment_short"].unique():
            genome_cache[s] = load_genome_de(s)
            log(f"  loaded genome DE for {s}: {len(genome_cache[s])} rows", fh=fh)

        results = []
        t_start = time.time()

        for i, row in scores.iterrows():
            axis = row["axis"]
            panel_kind = row["panel_kind"]
            short = row["experiment_short"]
            tp = row["timepoint"]
            observed_score = row["axis_score"]

            # Get the panel's gene list and directions
            panel = panels[(panels["axis"] == axis) & (panels["panel_kind"] == panel_kind)].copy()
            if panel.empty:
                continue
            axis_loci = panel["locus_tag"].tolist()
            axis_dirs = panel.set_index("locus_tag").loc[axis_loci, "direction"].astype(int).to_numpy()

            # Get genome DE at this TP
            genome = genome_cache[short]
            if pd.isna(tp) or tp == "single":
                tp_de = genome[genome["timepoint"].isna() | (genome["timepoint"].astype(str) == "nan")]
                if tp_de.empty:
                    # for axenic-RNA single contrast, no per-TP filter
                    tp_de = genome.copy()
            else:
                tp_de = genome[genome["timepoint"] == tp].copy()
            if tp_de.empty:
                log(f"  WARN: no DE rows for {axis} {panel_kind} {short} tp={tp}", fh=fh)
                continue

            # Build the log2fc array indexed by locus_tag
            tp_de = tp_de.drop_duplicates(subset=["locus_tag"]).reset_index(drop=True)
            locus_to_idx = {g: i for i, g in enumerate(tp_de["locus_tag"])}
            log2fc_arr = tp_de["log2fc"].to_numpy(dtype=float)

            # Get actual axis indices; some axis genes may be missing from this TP's DE
            actual_axis_idx = np.array(
                [locus_to_idx[g] for g in axis_loci if g in locus_to_idx], dtype=int
            )
            present_dirs = np.array(
                [d for g, d in zip(axis_loci, axis_dirs) if g in locus_to_idx], dtype=int
            )
            if len(actual_axis_idx) < 1:
                log(f"  SKIP: no axis genes present in DE for {axis}/{panel_kind}/{short}/{tp}", fh=fh)
                continue

            p_val, n_exceed, perm_mean, perm_sd = permutation_pvalue_for_cell(
                log2fc_arr=log2fc_arr,
                actual_axis_idx=actual_axis_idx,
                direction_arr=present_dirs,
                observed_score=float(observed_score) if pd.notna(observed_score) else float("nan"),
                rng=rng,
                K=K_PERMUTATIONS,
            )

            results.append({
                "axis": axis,
                "panel_kind": panel_kind,
                "experiment_short": short,
                "condition": row["condition"],
                "omics": row["omics"],
                "timepoint": tp,
                "timepoint_hours": row["timepoint_hours"],
                "growth_phase": row["growth_phase"],
                "n_axis_with_data": int(len(actual_axis_idx)),
                "axis_score": float(observed_score) if pd.notna(observed_score) else float("nan"),
                "permutation_p_raw": p_val,
                "permutation_n_exceed": n_exceed,
                "permutation_mean": perm_mean,
                "permutation_sd": perm_sd,
            })

            if (i + 1) % 10 == 0 or (i + 1) == len(scores):
                elapsed = time.time() - t_start
                log(f"  [{i+1}/{len(scores)}] elapsed {elapsed:.1f}s — last: {axis} {panel_kind} {short} tp={tp} → p={p_val:.4f}", fh=fh)

        out_df = pd.DataFrame(results)

        # BH-FDR correction across all cells with a valid p-value
        valid_mask = out_df["permutation_p_raw"].notna()
        out_df["permutation_p_bh"] = float("nan")
        if valid_mask.any():
            out_df.loc[valid_mask, "permutation_p_bh"] = bh_fdr(out_df.loc[valid_mask, "permutation_p_raw"].to_numpy())

        out_df["sig_05"] = (out_df["permutation_p_bh"] < 0.05) & out_df["permutation_p_bh"].notna()
        out_df["sig_01"] = (out_df["permutation_p_bh"] < 0.01) & out_df["permutation_p_bh"].notna()

        out_path = OUT_DIR / "permutation_pvalues.csv"
        out_df.to_csv(out_path, index=False)
        log(f"\nWrote {len(out_df)} rows to {out_path}", fh=fh)

        # Summary
        log("\nSummary (p_bh < 0.05):", fh=fh)
        sig = out_df[out_df["sig_05"]]
        log(f"  {len(sig)}/{len(out_df)} cells significant after BH-FDR correction", fh=fh)
        for _, r in sig.sort_values("permutation_p_bh").iterrows():
            log(
                f"  {r['axis']:<12} {r['panel_kind']:<10} {r['experiment_short']:<25} "
                f"tp={str(r['timepoint']):<12} z={r['axis_score']:+.2f}  p_bh={r['permutation_p_bh']:.4g}",
                fh=fh,
            )

        log(f"\nTotal elapsed: {time.time() - t_start:.1f}s", fh=fh)


if __name__ == "__main__":
    main()
