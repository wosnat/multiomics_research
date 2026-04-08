"""Signature building: summarize DE per gene, intersect references, classify.

All functions take DataFrames and return DataFrames. No KG dependency.
"""

import pandas as pd
import numpy as np


def summarize_per_gene(de_df: pd.DataFrame) -> pd.DataFrame:
    """Summarize a multi-timepoint DE DataFrame to one row per gene.

    For each gene:
    1. Filter to significant rows only
    2. Determine majority direction across timepoints
    3. On tie: break by best directional rank (lower = stronger signal wins)
    4. Pick the timepoint with best directional rank in chosen direction

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
    meta_cols = [c for c in ["product", "gene_category"] if c in de_df.columns]

    sig = de_df[de_df["expression_status"].isin(["significant_up", "significant_down"])].copy()
    if sig.empty:
        return pd.DataFrame(columns=result_cols + meta_cols)

    sig["direction"] = sig["expression_status"].map({
        "significant_up": "up", "significant_down": "down"
    })
    sig["dir_rank"] = np.where(
        sig["direction"] == "up",
        sig["rank_up"],
        sig["rank_down"],
    )

    # Majority direction per gene, with tie-breaking by best rank
    direction_counts = (
        sig.groupby("locus_tag")["direction"]
        .value_counts()
        .unstack(fill_value=0)
    )

    # Best directional rank per gene per direction
    best_rank_per_dir = (
        sig.groupby(["locus_tag", "direction"])["dir_rank"]
        .min()
        .unstack(fill_value=np.inf)
    )

    # Determine majority direction
    majority = []
    for tag in direction_counts.index:
        counts = direction_counts.loc[tag]
        up_count = counts.get("up", 0)
        down_count = counts.get("down", 0)

        if up_count > down_count:
            majority.append("up")
        elif down_count > up_count:
            majority.append("down")
        else:
            # Tie: break by best rank (lower rank = stronger signal)
            up_rank = best_rank_per_dir.loc[tag].get("up", np.inf)
            down_rank = best_rank_per_dir.loc[tag].get("down", np.inf)
            majority.append("down" if down_rank < up_rank else "up")

    majority_direction = pd.Series(majority, index=direction_counts.index, name="majority_direction")

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
    study_a_all_locus_tags: set[str] | None = None,
    study_b_all_locus_tags: set[str] | None = None,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """Intersect two per-gene DE summaries into core and extended signatures.

    Core: genes significant in both studies with concordant direction.
    Extended: genes significant in only one study.
    Discordant genes (different direction in each study) are excluded from both.

    Args:
        study_a: Output of summarize_per_gene for study A.
        study_b: Output of summarize_per_gene for study B.
        study_a_name: Label for study A (used in column prefixes).
        study_b_name: Label for study B (used in column prefixes).
        study_a_all_locus_tags: All locus tags in study A's dataset
            (not just significant). Used to classify study-B-only genes
            as "absent from dataset" vs "present but not significant".
        study_b_all_locus_tags: All locus tags in study B's dataset
            (not just significant). Used to classify study-A-only genes
            as "absent from dataset" vs "present but not significant".

    Returns:
        (core_df, extended_df, discordant_df) — core and extended sorted
        by cross_study_best_dir_rank; discordant contains genes significant
        in both studies but with opposite directions.
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

    # Rename metadata columns if present
    for meta in ["product", "gene_category"]:
        if meta in a.columns:
            a = a.rename(columns={meta: f"{meta}_a"})
        if meta in b.columns:
            b = b.rename(columns={meta: f"{meta}_b"})

    merged = a.merge(b, on="locus_tag", how="outer")

    # Resolve gene_name and metadata: prefer A, fall back to B
    for meta in ["gene_name", "product", "gene_category"]:
        col_a, col_b = f"{meta}_a", f"{meta}_b"
        if col_a in merged.columns and col_b in merged.columns:
            merged[meta] = merged[col_a].fillna(merged[col_b])
            merged = merged.drop(columns=[col_a, col_b])

    # Genes present in both studies
    both = merged.dropna(subset=["direction_a", "direction_b"])

    # DISCORDANT: both present, opposite direction
    discordant = both[both["direction_a"] != both["direction_b"]].copy()
    discordant["direction_a_label"] = discordant["direction_a"]
    discordant["direction_b_label"] = discordant["direction_b"]

    # CORE: both present, concordant direction
    concordant = both[both["direction_a"] == both["direction_b"]].copy()
    concordant["direction"] = concordant["direction_a"]
    concordant["signature_type"] = "core"
    concordant["cross_study_best_dir_rank"] = concordant[
        [f"{pa}_best_dir_rank", f"{pb}_best_dir_rank"]
    ].min(axis=1)

    # EXTENDED: in one study only (discordant excluded entirely)
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
    if study_a_all_locus_tags is not None:
        b_only["signature_type"] = b_only["locus_tag"].apply(
            lambda lt: f"{pb}_only_{pa}_ns" if lt in study_a_all_locus_tags
            else f"{pb}_only_{pa}_absent"
        )
    else:
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
    for meta in ["product", "gene_category"]:
        if meta in concordant.columns:
            cols.append(meta)

    for df in [concordant, extended]:
        for c in cols:
            if c not in df.columns:
                df[c] = pd.NA

    core_df = concordant[cols].sort_values("cross_study_best_dir_rank").reset_index(drop=True)
    extended_df = extended[cols].sort_values("cross_study_best_dir_rank").reset_index(drop=True)

    # Discordant output
    disc_cols = [
        "locus_tag", "gene_name", "direction_a_label", "direction_b_label",
        f"{pa}_best_dir_rank", f"{pb}_best_dir_rank",
    ]
    for meta in ["product", "gene_category"]:
        if meta in discordant.columns:
            disc_cols.append(meta)
    for c in disc_cols:
        if c not in discordant.columns:
            discordant[c] = pd.NA
    discordant_df = discordant[disc_cols].reset_index(drop=True)

    return core_df, extended_df, discordant_df
