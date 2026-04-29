"""
Extract the candidate gene set: significant_up × annotation_quality ≤ 1
from the locked MED4 axenic RNA-seq experiment (Weissberg 2025).

Inputs:
    --experiment-id   KG experiment id (default: locked MED4 axenic RNA-seq)
    --annotation-cut  max annotation_quality to include (default: 1)
    --out-dir         output directory (default: ../data)

Outputs (under --out-dir):
    sig_up_de_rows.csv         All sig_up DE rows for the experiment
    sig_up_405_overview.csv    gene_overview for every sig_up locus_tag
    candidate_set.csv          Sig_up × annotation_quality ≤ cut, sorted by log2fc desc
    01_extract_candidate_set.log   Filter funnel + per-facet availability counts
"""

import argparse
import logging
from pathlib import Path

import pandas as pd

from multiomics_explorer import differential_expression_by_gene, gene_overview
from multiomics_explorer.analysis import to_dataframe

LOCKED_EXPERIMENT_ID = (
    "10.1101/2025.11.24.690089_growth_state_pro99lown_nutrient_starvation_med4_rnaseq_axenic"
)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--experiment-id", default=LOCKED_EXPERIMENT_ID)
    parser.add_argument("--annotation-cut", type=int, default=1)
    parser.add_argument(
        "--out-dir",
        type=Path,
        default=Path(__file__).resolve().parent.parent / "data",
    )
    args = parser.parse_args()
    args.out_dir.mkdir(parents=True, exist_ok=True)

    log_path = args.out_dir / "01_extract_candidate_set.log"
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s %(levelname)s %(message)s",
        handlers=[logging.FileHandler(log_path, mode="w"), logging.StreamHandler()],
    )
    log = logging.getLogger(__name__)
    log.info("experiment_id=%s annotation_cut=%d", args.experiment_id, args.annotation_cut)

    # 1. Pull all sig_up DE rows for the locked experiment.
    de_result = differential_expression_by_gene(
        experiment_ids=[args.experiment_id],
        direction="up",
        significant_only=True,
        limit=None,
        verbose=True,
    )
    log.info(
        "DE call: matching_genes=%d total_matching=%d truncated=%s",
        de_result["matching_genes"],
        de_result["total_matching"],
        de_result["truncated"],
    )
    if de_result["truncated"]:
        raise RuntimeError("DE call was truncated despite limit=None")

    de_df = to_dataframe(de_result)
    de_path = args.out_dir / "sig_up_de_rows.csv"
    de_df.to_csv(de_path, index=False)
    log.info("wrote %s (%d rows, %d cols)", de_path, len(de_df), len(de_df.columns))

    locus_tags = sorted(de_df["locus_tag"].unique().tolist())
    log.info("distinct sig_up locus_tags: %d", len(locus_tags))
    if len(locus_tags) != de_result["matching_genes"]:
        raise RuntimeError(
            f"locus_tag set size {len(locus_tags)} != matching_genes {de_result['matching_genes']}"
        )

    # 2. Pull gene_overview for every sig_up locus_tag.
    ov_result = gene_overview(locus_tags=locus_tags, limit=None, verbose=True)
    log.info(
        "gene_overview: total_matching=%d returned=%d truncated=%s not_found=%s",
        ov_result["total_matching"],
        ov_result["returned"],
        ov_result["truncated"],
        ov_result["not_found"],
    )
    if ov_result["truncated"]:
        raise RuntimeError("gene_overview was truncated despite limit=None")

    ov_df = to_dataframe(ov_result)
    ov_path = args.out_dir / "sig_up_405_overview.csv"
    ov_df.to_csv(ov_path, index=False)
    log.info("wrote %s (%d rows, %d cols)", ov_path, len(ov_df), len(ov_df.columns))

    # 3. Distribution of annotation_quality across sig_up.
    aq_dist = ov_df["annotation_quality"].value_counts(dropna=False).sort_index()
    log.info("annotation_quality distribution among %d sig_up genes:", len(ov_df))
    for q, n in aq_dist.items():
        log.info("    annotation_quality=%s: %d", q, n)

    # 4. Build candidate set = sig_up × annotation_quality ≤ cut.
    de_keep = de_df[["locus_tag", "log2fc", "padj", "rank", "rank_up", "product", "gene_category"]]
    de_keep = de_keep.rename(columns={"product": "de_product", "gene_category": "de_gene_category"})
    cand_df = ov_df.merge(de_keep, on="locus_tag", how="inner")
    cand_df = cand_df[cand_df["annotation_quality"].notna()]
    cand_df = cand_df[cand_df["annotation_quality"] <= args.annotation_cut]
    cand_df = cand_df.sort_values("log2fc", ascending=False).reset_index(drop=True)

    cand_path = args.out_dir / "candidate_set.csv"
    cand_df.to_csv(cand_path, index=False)
    log.info(
        "candidate set (annotation_quality <= %d): %d genes — %s",
        args.annotation_cut,
        len(cand_df),
        cand_path,
    )

    # 5. Filter funnel.
    log.info("=" * 60)
    log.info("FILTER FUNNEL")
    log.info(
        "1849 detected (KG list_experiments) -> %d sig_up -> %d candidate (annotation_quality <= %d)",
        len(locus_tags),
        len(cand_df),
        args.annotation_cut,
    )

    # 6. Per-facet availability counts across the candidate set.
    log.info("=" * 60)
    log.info("CANDIDATE-SET FACET AVAILABILITY (from gene_overview):")
    facet_cols = {
        "annotation_types": "has any ontology",
        "expression_edge_count": "has expression data (>0)",
        "significant_up_count": "is sig_up in any other study (>0)",
        "closest_ortholog_group_size": "has ortholog group (>0)",
        "cluster_membership_count": "has cluster membership (>0)",
        "derived_metric_count": "has derived metrics (>0)",
    }
    for col, label in facet_cols.items():
        if col == "annotation_types":
            present = cand_df[col].fillna("").astype(str).str.len() > 0
        else:
            present = cand_df[col].fillna(0).astype(float) > 0
        log.info("    %-30s %s: %d / %d", col, label, present.sum(), len(cand_df))

    # 7. log2fc distribution within candidate set.
    log.info("=" * 60)
    log.info("CANDIDATE-SET log2fc summary:")
    log.info("    min=%.3f max=%.3f median=%.3f mean=%.3f",
             cand_df["log2fc"].min(),
             cand_df["log2fc"].max(),
             cand_df["log2fc"].median(),
             cand_df["log2fc"].mean())

    log.info("=" * 60)
    log.info("done.")


if __name__ == "__main__":
    main()
