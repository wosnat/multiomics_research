"""Distribution of signed_score across all Step 2 enrichment results.

Asks: is the spec's ±10 cap (SCORE_CAP) right for this dataset, or is it
clipping genuine signal?

Outputs:
    exploration/qc/step2_score_distribution.png  — multi-panel view
    stdout — quantile and tail summary
"""
from __future__ import annotations

import os
from pathlib import Path

os.environ.setdefault("MPLBACKEND", "Agg")

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

ANALYSIS_DIR = Path(__file__).resolve().parent.parent
SCORE_CAP = 10.0

CLASS_COLORS = {
    "T": "#2c7fb8",
    "R": "#d73027",
    "PC": "#fdae61",
    "NC": "#888888",
    "CTX": "#4daf4a",
}


def main() -> None:
    df = pd.read_csv(ANALYSIS_DIR / "data" / "enrichment_all.csv")
    classified = pd.read_csv(ANALYSIS_DIR / "data" / "experiments_classified.csv")
    exp_to_class = dict(zip(classified["experiment_id"], classified["class"]))
    df["class"] = df["experiment_id"].map(exp_to_class)

    # Clean / focus.
    sig = df[df["p_adjust"] < 0.05].copy()
    all_s = df["signed_score"].replace([np.inf, -np.inf], np.nan).dropna()
    sig_s = sig["signed_score"].replace([np.inf, -np.inf], np.nan).dropna()

    print("=== Score distribution ===")
    print(f"All rows: {len(df)}, finite signed_score: {len(all_s)}")
    print(f"Significant rows (padj<0.05): {len(sig)}, finite: {len(sig_s)}")

    def q(s, p):
        return float(np.quantile(np.abs(s), p))

    quants = [0.5, 0.75, 0.9, 0.95, 0.99, 0.995, 0.999, 1.0]
    print("\n|signed_score| quantiles (significant rows):")
    for p in quants:
        print(f"  p={p:5.3f}: {q(sig_s, p):8.2f}")

    print("\n|signed_score| quantiles (all rows):")
    for p in quants:
        print(f"  p={p:5.3f}: {q(all_s, p):8.2f}")

    thresholds = [3, 5, 7, 10, 15, 20, 30]
    print("\n# rows exceeding |signed_score| thresholds (significant only):")
    for t in thresholds:
        above = int((np.abs(sig_s) > t).sum())
        frac = above / max(len(sig_s), 1)
        print(f"  |s| > {t:2d}: {above:5d}  ({100*frac:5.2f}%)")

    # Top 15 significant rows by |signed_score|.
    top = sig.assign(abs_s=sig["signed_score"].abs()).sort_values(
        "abs_s", ascending=False
    ).head(15)
    print("\nTop-15 largest |signed_score| hits (significant rows):")
    cols = ["class", "ontology", "term_id", "term_name", "direction",
            "experiment_id", "timepoint", "signed_score", "p_adjust"]
    with pd.option_context("display.max_colwidth", 50, "display.width", 200):
        print(top[cols].to_string(index=False))

    # ---- Figure ----
    fig, axes = plt.subplots(2, 2, figsize=(12, 8))

    # Panel 1: histogram of signed_score, all rows, log-count.
    ax = axes[0, 0]
    ax.hist(all_s, bins=80, color="#cccccc", edgecolor="white", label="all")
    ax.hist(sig_s, bins=80, color="#4a77aa", edgecolor="white", alpha=0.75,
            label="padj<0.05")
    ax.axvline(SCORE_CAP, color="red", ls="--", lw=1, label=f"±{int(SCORE_CAP)} cap")
    ax.axvline(-SCORE_CAP, color="red", ls="--", lw=1)
    ax.set_yscale("log")
    ax.set_xlabel("signed_score")
    ax.set_ylabel("count (log)")
    ax.set_title("Panel A — Histogram (all vs significant)")
    ax.legend(fontsize=8)

    # Panel 2: histogram of significant rows by class, |signed_score|.
    ax = axes[0, 1]
    for cls, color in CLASS_COLORS.items():
        sub = sig.loc[sig["class"] == cls, "signed_score"].abs()
        sub = sub.replace([np.inf], np.nan).dropna()
        if len(sub) == 0:
            continue
        ax.hist(sub, bins=50, alpha=0.55, color=color,
                label=f"{cls} (n={len(sub)})",
                edgecolor="white")
    ax.axvline(SCORE_CAP, color="red", ls="--", lw=1, label=f"{int(SCORE_CAP)} cap")
    ax.set_yscale("log")
    ax.set_xlabel("|signed_score|")
    ax.set_ylabel("count (log)")
    ax.set_title("Panel B — |signed_score| by class (significant only)")
    ax.legend(fontsize=8)

    # Panel 3: CDF of |signed_score|.
    ax = axes[1, 0]
    abs_sig = np.sort(np.abs(sig_s.values))
    cdf = np.arange(1, len(abs_sig) + 1) / len(abs_sig)
    ax.plot(abs_sig, cdf, color="#4a77aa")
    for p in [0.5, 0.9, 0.95, 0.99]:
        val = q(sig_s, p)
        ax.axvline(val, color="#888", ls=":", lw=0.7)
        ax.text(val, p, f" p{int(p*100)}={val:.1f}", fontsize=7, va="center")
    ax.axvline(SCORE_CAP, color="red", ls="--", lw=1, label=f"{int(SCORE_CAP)} cap")
    ax.set_xlabel("|signed_score|")
    ax.set_ylabel("CDF (significant rows)")
    ax.set_title("Panel C — CDF of |signed_score| (significant)")
    ax.legend(fontsize=8)
    ax.grid(alpha=0.3)

    # Panel 4: sorted bar of top-50 largest |signed_score|.
    ax = axes[1, 1]
    top50 = sig.assign(abs_s=sig["signed_score"].abs()).sort_values(
        "abs_s", ascending=False
    ).head(50)
    colors = [CLASS_COLORS.get(c, "#888") for c in top50["class"]]
    ax.barh(range(len(top50)), top50["abs_s"], color=colors)
    ax.axvline(SCORE_CAP, color="red", ls="--", lw=1)
    # Short row labels
    labels = [f"{r['term_id'].split(':', 1)[-1]}|{r['class']}"
              for _, r in top50.iterrows()]
    ax.set_yticks(range(len(top50)))
    ax.set_yticklabels(labels, fontsize=6)
    ax.invert_yaxis()
    ax.set_xlabel("|signed_score|")
    ax.set_title("Panel D — Top 50 largest hits")

    plt.tight_layout()
    outpath = ANALYSIS_DIR / "exploration" / "qc" / "step2_score_distribution.png"
    fig.savefig(outpath, dpi=140, bbox_inches="tight")
    plt.close(fig)
    print(f"\nwrote {outpath}")


if __name__ == "__main__":
    main()
