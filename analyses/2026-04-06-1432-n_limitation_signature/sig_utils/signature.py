"""Signature building: intersect DE lists, assign directions and ranks."""

import pandas as pd


def summarize_de_per_gene(de_df: pd.DataFrame) -> pd.DataFrame:
    """Summarize a multi-timepoint DE DataFrame to one row per gene.

    For each gene, determines majority direction across significant timepoints,
    then picks the timepoint with best directional rank in that direction.

    Args:
        de_df: DataFrame with columns: locus_tag, gene_name, timepoint,
            log2fc, expression_status, rank, rank_up, rank_down.

    Returns:
        DataFrame with one row per gene: locus_tag, gene_name, direction,
        peak_timepoint, best_dir_rank, best_global_rank.
    """
    empty = pd.DataFrame(columns=[
        "locus_tag", "gene_name", "direction",
        "peak_timepoint", "best_dir_rank", "best_global_rank",
    ])

    sig = de_df[de_df["expression_status"].isin(["significant_up", "significant_down"])].copy()
    if sig.empty:
        return empty

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

    # Peak timepoint: where gene has best directional rank
    peak = sig_majority.loc[sig_majority.groupby("locus_tag")["dir_rank"].idxmin()]

    # Best directional rank across all timepoints in majority direction
    best_dir = sig_majority.groupby("locus_tag")["dir_rank"].min().rename("best_dir_rank")

    # Best global rank across all significant timepoints
    best_global = sig.groupby("locus_tag")["rank"].min().rename("best_global_rank")

    result = peak[["locus_tag", "gene_name", "direction", "timepoint"]].copy()
    result = result.rename(columns={"timepoint": "peak_timepoint"})
    result["direction"] = result["locus_tag"].map(majority_direction)
    result = result.merge(best_dir, on="locus_tag")
    result = result.merge(best_global, on="locus_tag")

    return result.reset_index(drop=True)


def intersect_de_lists(
    study_a: pd.DataFrame,
    study_b: pd.DataFrame,
    study_a_name: str = "study_a",
    study_b_name: str = "study_b",
    study_b_all_locus_tags: set[str] | None = None,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Intersect two per-gene DE summaries into core and extended signatures.

    Core: genes in both studies with concordant direction.
    Extended: genes in only one study, tagged by source.

    Args:
        study_a: Output of summarize_de_per_gene for study A.
        study_b: Output of summarize_de_per_gene for study B.
        study_a_name: Label for study A (e.g. 'tolonen'). Used in column prefixes.
        study_b_name: Label for study B (e.g. 'read'). Used in column prefixes.
        study_b_all_locus_tags: Set of all locus tags present in study B's
            dataset (not just significant). If None, cannot distinguish
            'absent from dataset' vs 'present but not significant'.

    Returns:
        (core_df, extended_df)
    """
    pa, pb = study_a_name, study_b_name

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

    merged = a.merge(b, on="locus_tag", how="outer")
    merged["gene_name"] = merged["gene_name_a"].fillna(merged["gene_name_b"])
    merged = merged.drop(columns=["gene_name_a", "gene_name_b"])

    # Core: both present, concordant direction
    both = merged.dropna(subset=["direction_a", "direction_b"])
    concordant = both[both["direction_a"] == both["direction_b"]].copy()
    concordant["direction"] = concordant["direction_a"]
    concordant["signature_type"] = "core"

    concordant["cross_study_best_dir_rank"] = concordant[
        [f"{pa}_best_dir_rank", f"{pb}_best_dir_rank"]
    ].min(axis=1)

    # Extended: in one study only
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
        lambda r: r[f"{pa}_best_dir_rank"] if pd.notna(r[f"{pa}_best_dir_rank"])
        else r[f"{pb}_best_dir_rank"],
        axis=1,
    )

    # Select and order columns
    cols = [
        "locus_tag", "gene_name", "direction", "signature_type",
        f"{pa}_peak_timepoint",
        f"{pa}_best_dir_rank", f"{pa}_best_global_rank",
        f"{pb}_peak_timepoint",
        f"{pb}_best_dir_rank", f"{pb}_best_global_rank",
        "cross_study_best_dir_rank",
    ]
    for df in [concordant, extended]:
        for c in cols:
            if c not in df.columns:
                df[c] = pd.NA

    core_df = concordant[cols].sort_values("cross_study_best_dir_rank").reset_index(drop=True)
    extended_df = extended[cols].sort_values("cross_study_best_dir_rank").reset_index(drop=True)

    return core_df, extended_df
