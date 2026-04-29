"""
Aggregate views over the 116-candidate dossier corpus.

Produces:
  - Per-axis fill-rate sanity check (re-confirms step-2/3 numbers)
  - Cluster pivot: per-cluster candidate count + curated description, sorted desc
  - Treatment pivot: per-treatment candidate response count + breadth distribution
  - Floor genes table (F1 no-ontology + F3 RefSeq-only-no-cluster-or-homolog)
  - Broad-stress responders table (breadth=4 candidates)
  - Top-FC dossier signal summary (top 10 by log2fc, with cluster + ortholog probe takeaways)
  - Ortholog probe headline: cumulative `has_*` counts across the candidate set's groups

Inputs (default ../data/):
    dossiers.json
    group_probe_cache.json

Outputs (default ../data/, ../figures/):
    data/agg_axis_fill.csv
    data/agg_clusters_pivot.csv
    data/agg_treatments_pivot.csv
    data/agg_floor_genes.csv
    data/agg_broad_responders.csv
    data/agg_top_fc_signals.csv
    data/agg_ortholog_probe_summary.csv
    figures/agg_clusters_pivot.png
    figures/agg_treatments_pivot.png
    figures/agg_breadth_distribution.png
    figures/agg_axis_fill.png
"""

import argparse
import json
import logging
from collections import Counter
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

    log_path = args.data_dir / "02_aggregate_views.log"
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s %(levelname)s %(message)s",
        handlers=[logging.FileHandler(log_path, mode="w"), logging.StreamHandler()],
    )
    log = logging.getLogger(__name__)

    cards = json.loads((args.data_dir / "dossiers.json").read_text())
    n = len(cards)
    log.info("loaded %d candidate dossiers", n)

    # =================================================================
    # 1. Per-axis fill rate
    # =================================================================
    axis_fill = []
    fields = [
        ("clusters", lambda c: bool(c.get("clusters"))),
        ("ortholog_groups", lambda c: bool(c.get("ortholog_groups"))),
        ("ontology_terms", lambda c: bool(c.get("ontology"))),
        ("derived_metrics", lambda c: bool(c.get("derived_metrics"))),
    ]
    for axis, has in fields:
        n_with = sum(1 for c in cards if has(c))
        axis_fill.append({"axis": axis, "n_with": n_with, "n_total": n, "fraction": round(n_with / n, 3)})
    pd.DataFrame(axis_fill).to_csv(args.data_dir / "agg_axis_fill.csv", index=False)

    # Bar figure
    af_df = pd.DataFrame(axis_fill)
    fig, ax = plt.subplots(figsize=(7, 4))
    bars = ax.bar(af_df["axis"], af_df["n_with"], color="#1f77b4")
    for b, v in zip(bars, af_df["n_with"]):
        ax.text(b.get_x() + b.get_width() / 2, b.get_height(), f"{v}", ha="center", va="bottom")
    ax.axhline(n, color="grey", linestyle="--", alpha=0.5, label=f"n={n} candidates")
    ax.set_ylabel("# candidates with data")
    ax.set_title(f"Per-axis fill rate across candidate set (n={n})")
    ax.legend()
    fig.tight_layout()
    fig.savefig(args.fig_dir / "agg_axis_fill.png", dpi=300)
    plt.close(fig)
    log.info("wrote agg_axis_fill (CSV + PNG)")

    # =================================================================
    # 2. Cluster pivot
    # =================================================================
    cluster_rows = []
    for c in cards:
        for cl in c.get("clusters") or []:
            cluster_rows.append({
                "locus_tag": c["locus_tag"],
                "cluster_id": cl.get("cluster_id"),
                "cluster_name": cl.get("cluster_name"),
                "cluster_type": cl.get("cluster_type"),
                "member_count": cl.get("member_count"),
                "cluster_functional_description": cl.get("cluster_functional_description"),
                "cluster_expression_dynamics": cl.get("cluster_expression_dynamics"),
                "cluster_temporal_pattern": cl.get("cluster_temporal_pattern"),
                "analysis_id": cl.get("analysis_id"),
                "analysis_name": cl.get("analysis_name"),
            })
    cdf = pd.DataFrame(cluster_rows)

    if not cdf.empty:
        pivot = cdf.groupby(
            ["cluster_id", "cluster_name", "cluster_type", "analysis_name", "member_count"]
        ).agg(
            n_candidates=("locus_tag", "nunique"),
            cluster_functional_description=("cluster_functional_description", "first"),
            cluster_expression_dynamics=("cluster_expression_dynamics", "first"),
            cluster_temporal_pattern=("cluster_temporal_pattern", "first"),
        ).reset_index().sort_values("n_candidates", ascending=False)
        pivot["candidate_share_of_cluster"] = (pivot["n_candidates"] / pivot["member_count"]).round(3)
        pivot.to_csv(args.data_dir / "agg_clusters_pivot.csv", index=False)
        log.info("wrote agg_clusters_pivot (%d unique clusters touched)", len(pivot))

        # Top-N cluster figure
        top = pivot.head(15).sort_values("n_candidates")
        fig, ax = plt.subplots(figsize=(11, 0.45 * len(top) + 1))
        y = np.arange(len(top))
        bars = ax.barh(y, top["n_candidates"].values, color="#1f77b4")
        labels = [
            f"[{r.cluster_type}] {r.cluster_name}\n  ({int(r.member_count)} members; {r.analysis_name})"
            for r in top.itertuples()
        ]
        ax.set_yticks(y)
        ax.set_yticklabels(labels, fontsize=8)
        for b, v, total in zip(bars, top["n_candidates"].values, top["member_count"].values):
            ax.text(v + 0.5, b.get_y() + b.get_height() / 2,
                    f"{v}/{int(total)} ({v / total:.0%})", va="center", fontsize=8)
        ax.set_xlabel("# candidates in cluster (n_candidates / cluster_member_count)")
        ax.set_title(f"Top-15 clusters by candidate count (n_candidates total = {n})")
        fig.tight_layout()
        fig.savefig(args.fig_dir / "agg_clusters_pivot.png", dpi=300)
        plt.close(fig)
        log.info("wrote figures/agg_clusters_pivot.png")

    # =================================================================
    # 3. Treatment pivot + breadth distribution
    # =================================================================
    treatment_rows = []
    breadth = []
    for c in cards:
        resp = c.get("response") or {}
        full = resp.get("full") or {}
        groups_responded = full.get("groups_responded") or []
        breadth.append(len(groups_responded))
        per_t = full.get("per_treatment") or {}
        for t, ts in per_t.items():
            treatment_rows.append({
                "locus_tag": c["locus_tag"],
                "treatment": t,
                "experiments_total": ts.get("experiments_total", 0),
                "experiments_tested": ts.get("experiments_tested", 0),
                "experiments_up": ts.get("experiments_up", 0),
                "experiments_down": ts.get("experiments_down", 0),
                "responded": t in groups_responded,
            })
    tdf = pd.DataFrame(treatment_rows)
    if not tdf.empty:
        treatment_pivot = tdf.groupby("treatment").agg(
            n_candidates_responded=("responded", "sum"),
            n_candidates_with_treatment_data=("experiments_tested", lambda s: (s > 0).sum()),
        ).reset_index()
        treatment_pivot["fraction_responded"] = (
            treatment_pivot["n_candidates_responded"] / n
        ).round(3)
        treatment_pivot = treatment_pivot.sort_values("n_candidates_responded", ascending=False)
        treatment_pivot.to_csv(args.data_dir / "agg_treatments_pivot.csv", index=False)
        log.info("wrote agg_treatments_pivot")

        # Treatment figure
        fig, ax = plt.subplots(figsize=(8, 4))
        bars = ax.bar(
            treatment_pivot["treatment"],
            treatment_pivot["n_candidates_responded"],
            color="#2ca02c",
        )
        for b, v in zip(bars, treatment_pivot["n_candidates_responded"]):
            ax.text(b.get_x() + b.get_width() / 2, b.get_height(), f"{int(v)}",
                    ha="center", va="bottom", fontsize=9)
        ax.set_ylabel("# candidates responded")
        ax.set_title(f"Per-treatment candidate response count (n={n})")
        ax.tick_params(axis="x", rotation=20)
        fig.tight_layout()
        fig.savefig(args.fig_dir / "agg_treatments_pivot.png", dpi=300)
        plt.close(fig)

    # Breadth distribution
    breadth_counts = Counter(breadth)
    breadth_rows = sorted(
        [{"n_treatments_responded": k, "n_candidates": v} for k, v in breadth_counts.items()],
        key=lambda r: r["n_treatments_responded"],
    )
    pd.DataFrame(breadth_rows).to_csv(args.data_dir / "agg_breadth_distribution.csv", index=False)
    fig, ax = plt.subplots(figsize=(7, 4))
    bars = ax.bar(
        [r["n_treatments_responded"] for r in breadth_rows],
        [r["n_candidates"] for r in breadth_rows],
        color="#ff7f0e",
    )
    for b, v in zip(bars, [r["n_candidates"] for r in breadth_rows]):
        ax.text(b.get_x() + b.get_width() / 2, b.get_height(), f"{v}",
                ha="center", va="bottom")
    ax.set_xlabel("# treatments where candidate responded")
    ax.set_ylabel("# candidates")
    ax.set_title(f"Cross-treatment response breadth distribution (n={n})")
    fig.tight_layout()
    fig.savefig(args.fig_dir / "agg_breadth_distribution.png", dpi=300)
    plt.close(fig)
    log.info("wrote agg_breadth_distribution (CSV + PNG)")

    # =================================================================
    # 4. Floor genes (F1 no-ontology, F3 no-cluster + no-homolog)
    # =================================================================
    floor_rows = []
    for c in cards:
        n_clusters = len(c.get("clusters") or [])
        n_groups = len(c.get("ortholog_groups") or [])
        n_ontology = len(c.get("ontology") or [])
        n_dms = len(c.get("derived_metrics") or [])
        f1 = n_ontology == 0
        f3 = n_clusters == 0 and n_groups == 0  # F3 = no cluster AND no homolog (RefSeq-only floor)
        if f1 or f3:
            ident = c.get("identity") or {}
            de = c.get("de_evidence") or {}
            floor_rows.append({
                "locus_tag": c["locus_tag"],
                "gene_name": ident.get("gene_name"),
                "product": ident.get("product"),
                "annotation_quality": ident.get("annotation_quality"),
                "log2fc": de.get("log2fc"),
                "n_clusters": n_clusters,
                "n_ortholog_groups": n_groups,
                "n_ontology": n_ontology,
                "n_dms": n_dms,
                "is_F1_no_ontology": f1,
                "is_F3_no_cluster_no_homolog": f3,
            })
    fdf = pd.DataFrame(floor_rows).sort_values("log2fc", ascending=False)
    fdf.to_csv(args.data_dir / "agg_floor_genes.csv", index=False)
    log.info("wrote agg_floor_genes (%d genes; F1=%d, F3=%d)",
             len(fdf), fdf["is_F1_no_ontology"].sum(), fdf["is_F3_no_cluster_no_homolog"].sum())

    # =================================================================
    # 5. Broad-stress responders (breadth >= 3)
    # =================================================================
    broad_rows = []
    for c in cards:
        groups_responded = (c.get("response") or {}).get("full", {}).get("groups_responded") or []
        if len(groups_responded) >= 3:
            ident = c.get("identity") or {}
            de = c.get("de_evidence") or {}
            broad_rows.append({
                "locus_tag": c["locus_tag"],
                "gene_name": ident.get("gene_name"),
                "product": ident.get("product"),
                "annotation_quality": ident.get("annotation_quality"),
                "log2fc": de.get("log2fc"),
                "breadth": len(groups_responded),
                "groups_responded": ", ".join(groups_responded),
            })
    bdf = pd.DataFrame(broad_rows).sort_values(["breadth", "log2fc"], ascending=[False, False])
    bdf.to_csv(args.data_dir / "agg_broad_responders.csv", index=False)
    log.info("wrote agg_broad_responders (%d genes with breadth>=3; breadth=4: %d)",
             len(bdf), (bdf["breadth"] == 4).sum())

    # =================================================================
    # 6. Top-FC signals
    # =================================================================
    top_rows = []
    for c in sorted(cards, key=lambda x: -((x.get("de_evidence") or {}).get("log2fc") or 0))[:10]:
        ident = c.get("identity") or {}
        de = c.get("de_evidence") or {}
        cluster_anchors = []
        for cl in c.get("clusters") or []:
            desc = cl.get("cluster_functional_description") or ""
            if desc and desc.strip() != "N/A":
                cluster_anchors.append(f"[{cl['cluster_type']}] {cl['cluster_name']}: {desc}")
        n_groups = len(c.get("ortholog_groups") or [])
        n_ontology = len(c.get("ontology") or [])
        n_dms = len(c.get("derived_metrics") or [])
        groups_responded = (c.get("response") or {}).get("full", {}).get("groups_responded") or []
        top_rows.append({
            "locus_tag": c["locus_tag"],
            "gene_name": ident.get("gene_name"),
            "product": ident.get("product"),
            "annotation_quality": ident.get("annotation_quality"),
            "log2fc": de.get("log2fc"),
            "rank_up": de.get("rank_up"),
            "n_clusters": len(c.get("clusters") or []),
            "cluster_with_curated_desc": cluster_anchors[0] if cluster_anchors else "",
            "n_ortholog_groups": n_groups,
            "n_ontology": n_ontology,
            "n_dms": n_dms,
            "breadth": len(groups_responded),
            "groups_responded": ", ".join(groups_responded),
        })
    pd.DataFrame(top_rows).to_csv(args.data_dir / "agg_top_fc_signals.csv", index=False)
    log.info("wrote agg_top_fc_signals (top 10 by log2fc)")

    # =================================================================
    # 7. Ortholog probe summary across the candidate set
    # =================================================================
    probe_path = args.data_dir / "group_probe_cache.json"
    if probe_path.exists():
        cache_data = json.loads(probe_path.read_text())
        rows = []
        for group_id, probe in cache_data.items():
            if probe is None:
                rows.append({"group_id": group_id, "is_singleton": True,
                             "n_members_excl_candidate": 0})
                continue
            rows.append({
                "group_id": group_id,
                "is_singleton": False,
                "n_members_excl_candidate": probe.get("n_members_excl_candidate", 0),
                "has_expression": probe.get("has_expression", 0),
                "has_significant_expression": probe.get("has_significant_expression", 0),
                "has_clusters": probe.get("has_clusters", 0),
                "has_derived_metrics": probe.get("has_derived_metrics", 0),
            })
        pdf = pd.DataFrame(rows)
        pdf.to_csv(args.data_dir / "agg_ortholog_probe_summary.csv", index=False)
        # Headline counts
        non_singleton = pdf[~pdf["is_singleton"]]
        log.info("ortholog probe summary: %d total probes, %d singletons; "
                 "non-singleton groups have median has_expression=%s, has_clusters=%s, has_derived_metrics=%s",
                 len(pdf), pdf["is_singleton"].sum(),
                 non_singleton["has_expression"].median() if len(non_singleton) else None,
                 non_singleton["has_clusters"].median() if len(non_singleton) else None,
                 non_singleton["has_derived_metrics"].median() if len(non_singleton) else None)

    log.info("done")


if __name__ == "__main__":
    main()
