"""
QC plots for the candidate-set filter funnel and FC-stratified facet availability.

Inputs (default: ../data/):
    sig_up_de_rows.csv        all sig_up DE rows (for full-sig_up FC distribution background)
    sig_up_405_overview.csv   gene_overview for all 405 sig_up
    candidate_set.csv         sig_up x annotation_quality <= 1, sorted by log2fc desc

Outputs (default: ../figures/, ../data/):
    figures/qc_annotation_quality_dist.png       AQ distribution within sig_up + within candidate set
    figures/qc_log2fc_by_annotation_quality.png  log2fc histogram: full sig_up (n=405) background + candidate set split by AQ
    figures/qc_facet_availability_by_fc.png      facet-availability counts in top-quartile FC vs rest of candidate set
    data/qc_facet_availability_by_fc.csv         tabular form of the above
"""

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data-dir", type=Path, default=Path(__file__).resolve().parent.parent / "data")
    parser.add_argument("--fig-dir", type=Path, default=Path(__file__).resolve().parent.parent / "figures")
    args = parser.parse_args()
    args.fig_dir.mkdir(parents=True, exist_ok=True)

    sig_up = pd.read_csv(args.data_dir / "sig_up_405_overview.csv")
    cand = pd.read_csv(args.data_dir / "candidate_set.csv")
    de_rows = pd.read_csv(args.data_dir / "sig_up_de_rows.csv")  # for full-sig_up log2fc background

    # Figure 1: annotation_quality distribution among 405 sig_up (highlights candidate cut at <= 1)
    fig, ax = plt.subplots(figsize=(7, 4))
    aq_counts = sig_up["annotation_quality"].value_counts(dropna=False).sort_index()
    bars = ax.bar(aq_counts.index.astype(str), aq_counts.values, color=["#d62728", "#ff7f0e", "#2ca02c", "#1f77b4"])
    for bar, n in zip(bars, aq_counts.values):
        ax.text(bar.get_x() + bar.get_width() / 2, bar.get_height(), f"{int(n)}", ha="center", va="bottom")
    ax.set_xlabel("annotation_quality (0=hypothetical, 1=has description, 2=named product, 3=well-annotated)")
    ax.set_ylabel("# genes")
    ax.set_title(f"annotation_quality distribution among {len(sig_up)} sig_up genes\n"
                 f"candidate set (<=1) = {(sig_up['annotation_quality'] <= 1).sum()}")
    ax.axvline(x=1.5, color="black", linestyle="--", alpha=0.5)
    ax.text(1.5, max(aq_counts.values) * 0.9, "candidate cut\n(<=1)", ha="center", fontsize=9)
    fig.tight_layout()
    fig.savefig(args.fig_dir / "qc_annotation_quality_dist.png", dpi=300)
    plt.close(fig)

    # Figure 2: log2fc histogram with full sig_up (n=405) as background, candidate set split by AQ on top.
    # The full sig_up distribution shows where the candidate set sits in the broader response.
    fig, ax = plt.subplots(figsize=(9, 4.5))
    sig_up_log2fc = de_rows["log2fc"]
    upper = max(sig_up_log2fc.max(), cand["log2fc"].max()) + 0.5
    bins = np.linspace(0.5, upper, 30)
    # Background: full 405 sig_up
    ax.hist(sig_up_log2fc, bins=bins, color="#cccccc", alpha=0.9,
            label=f"all sig_up (n={len(sig_up_log2fc)})", zorder=1)
    # Foreground: candidate set, split by AQ
    aq0 = cand.loc[cand["annotation_quality"] == 0, "log2fc"]
    aq1 = cand.loc[cand["annotation_quality"] == 1, "log2fc"]
    ax.hist(aq1, bins=bins, color="#ff7f0e", alpha=0.85,
            label=f"candidate AQ=1 (has desc, no named product) n={len(aq1)}", zorder=2)
    ax.hist(aq0, bins=bins, color="#d62728", alpha=0.85,
            label=f"candidate AQ=0 (pure hypothetical) n={len(aq0)}", zorder=3)
    ax.set_xlabel("log2fc")
    ax.set_ylabel("# genes")
    ax.set_title(
        f"log2fc distribution: full sig_up (background) vs candidate set (foreground)\n"
        f"sig_up n={len(sig_up_log2fc)} max={sig_up_log2fc.max():.2f} median={sig_up_log2fc.median():.2f} | "
        f"candidate n={len(cand)} max={cand['log2fc'].max():.2f} median={cand['log2fc'].median():.2f}"
    )
    ax.legend(loc="upper right")
    fig.tight_layout()
    fig.savefig(args.fig_dir / "qc_log2fc_by_annotation_quality.png", dpi=300)
    plt.close(fig)

    # Figure 3 + table: facet availability stratified by FC bucket
    # Top quartile = top 25% by log2fc within candidate set; rest = bottom 75%
    fc_q3 = cand["log2fc"].quantile(0.75)
    cand = cand.assign(fc_bucket=np.where(cand["log2fc"] >= fc_q3, "top_quartile", "rest"))

    facets = {
        "annotation_types": ("has_any_ontology", lambda s: s.fillna("").astype(str).str.len() > 0),
        "expression_edge_count": ("has_expression", lambda s: s.fillna(0).astype(float) > 0),
        "closest_ortholog_group_size": ("has_orthologs", lambda s: s.fillna(0).astype(float) > 0),
        "cluster_membership_count": ("has_clusters", lambda s: s.fillna(0).astype(float) > 0),
        "derived_metric_count": ("has_derived_metrics", lambda s: s.fillna(0).astype(float) > 0),
    }
    rows = []
    for col, (label, test) in facets.items():
        present = test(cand[col])
        for bucket in ["top_quartile", "rest"]:
            mask = cand["fc_bucket"] == bucket
            rows.append({
                "facet": label,
                "fc_bucket": bucket,
                "n_with": int(present[mask].sum()),
                "n_total": int(mask.sum()),
                "fraction_with": float(present[mask].sum() / mask.sum()) if mask.sum() else float("nan"),
            })
    avail = pd.DataFrame(rows)
    avail.to_csv(args.data_dir / "qc_facet_availability_by_fc.csv", index=False)

    # Stacked bar plot
    fig, ax = plt.subplots(figsize=(9, 4))
    pivot = avail.pivot(index="facet", columns="fc_bucket", values="fraction_with")
    pivot = pivot.reindex(columns=["top_quartile", "rest"])
    x = np.arange(len(pivot))
    width = 0.4
    ax.bar(x - width / 2, pivot["top_quartile"].values, width, label=f"top quartile FC (log2fc >= {fc_q3:.2f})", color="#1f77b4")
    ax.bar(x + width / 2, pivot["rest"].values, width, label="rest", color="#aec7e8")
    ax.set_xticks(x)
    ax.set_xticklabels(pivot.index, rotation=20, ha="right")
    ax.set_ylabel("fraction with data")
    ax.set_ylim(0, 1.05)
    ax.set_title(f"Facet availability in candidate set by FC bucket (n_top={int((cand['fc_bucket']=='top_quartile').sum())}, n_rest={int((cand['fc_bucket']=='rest').sum())})")
    ax.legend(loc="lower left")
    for i, facet in enumerate(pivot.index):
        for j, bucket in enumerate(["top_quartile", "rest"]):
            offset = -width / 2 if j == 0 else width / 2
            v = pivot.loc[facet, bucket]
            ax.text(i + offset, v + 0.02, f"{v:.0%}", ha="center", fontsize=8)
    fig.tight_layout()
    fig.savefig(args.fig_dir / "qc_facet_availability_by_fc.png", dpi=300)
    plt.close(fig)

    print("Wrote figures to", args.fig_dir)
    print("Wrote table to", args.data_dir / "qc_facet_availability_by_fc.csv")
    print("\nFC bucket boundaries: top_quartile = log2fc >=", round(fc_q3, 3))
    print("\nFacet availability table:")
    print(avail.to_string(index=False))


if __name__ == "__main__":
    main()
