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
        concordance (+1, -1, or 0), dir_rank.
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
    de_df: pd.DataFrame,
) -> float:
    """Compute rank score from an applied signature DataFrame.

    rank_score = mean(concordance_i * normalized_rank_i)
    where normalized_rank_i = 1 - (dir_rank_i / n_significant_in_direction).

    Normalizing by n_significant (not total genes) gives meaningful spread:
    rank 10 of 300 significant_up genes → 0.967, rank 200 of 300 → 0.333.

    Only includes genes that are present in the dataset (non-NaN
    expression_status). Genes absent from the dataset are excluded.

    Args:
        applied_df: Output of apply_signature.
        de_df: Full DE data for this condition x timepoint (needed to count
            significant genes per direction).

    Returns:
        Float score. Positive = N-limitation signal, negative = reversal, 0 = no signal.
        NaN if no genes are matched.
    """
    # Exclude genes absent from the dataset
    matched = applied_df[applied_df["expression_status"].notna()]

    if len(matched) == 0:
        return np.nan

    # Count significant genes per direction in the full experiment
    n_sig_up = (de_df["expression_status"] == "significant_up").sum()
    n_sig_down = (de_df["expression_status"] == "significant_down").sum()

    # Normalize by n_significant in the gene's DE direction
    n_sig_for_gene = np.where(
        matched["expression_status"] == "significant_up", n_sig_up,
        np.where(
            matched["expression_status"] == "significant_down", n_sig_down,
            1,  # avoid division by zero for not_significant (will be multiplied by 0)
        )
    )

    normalized_rank = np.where(
        pd.notna(matched["dir_rank"]),
        1 - (matched["dir_rank"].astype(float) / n_sig_for_gene),
        0,
    )

    contribution = matched["concordance"].values * normalized_rank
    return float(np.mean(contribution))


def permutation_test(
    signature_df: pd.DataFrame,
    de_df: pd.DataFrame,
    n_perms: int = 1000,
    seed: int = 42,
) -> dict:
    """Permutation test for rank score significance.

    Generates n_perms random gene sets of the same size as the matched
    signature genes, computes rank_score for each, returns empirical p-value.

    Args:
        signature_df: Signature with locus_tag, direction.
        de_df: Full DE data for one condition x timepoint (all genes).
        n_perms: Number of permutations.
        seed: Random seed.

    Returns:
        dict with: observed, empirical_p, n_permutations, n_signature_genes.
    """
    applied = apply_signature(signature_df, de_df)
    observed = rank_score(applied, de_df)
    n_matched = applied["expression_status"].notna().sum()

    if np.isnan(observed) or n_matched < 30:
        return {
            "observed": observed,
            "empirical_p": np.nan,
            "n_permutations": 0,
            "n_signature_genes": len(signature_df),
        }

    rng = np.random.default_rng(seed)
    all_tags = de_df["locus_tag"].unique()
    directions = signature_df["direction"].values

    # Match the size to what's actually in the dataset
    n_draw = min(len(directions), len(all_tags))

    null_scores = np.empty(n_perms)
    for i in range(n_perms):
        random_tags = rng.choice(all_tags, size=n_draw, replace=False)
        # Recycle direction labels to match draw size
        random_dirs = directions[:n_draw] if n_draw <= len(directions) else np.resize(directions, n_draw)
        random_sig = pd.DataFrame({
            "locus_tag": random_tags,
            "direction": random_dirs,
        })
        random_applied = apply_signature(random_sig, de_df)
        null_scores[i] = rank_score(random_applied, de_df)

    empirical_p = float((np.abs(null_scores) >= np.abs(observed)).mean())

    return {
        "observed": float(observed),
        "empirical_p": empirical_p,
        "n_permutations": n_perms,
        "n_signature_genes": len(signature_df),
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
        dict with: rank_score, empirical_p, n_permutations,
        n_signature_genes, n_matched, n_absent,
        n_concordant, n_reversed, n_not_significant, hit_rate.
    """
    applied = apply_signature(signature_df, de_df)
    score = rank_score(applied, de_df)

    perm_result = permutation_test(signature_df, de_df, n_perms=n_perms, seed=seed)

    # Breakdown
    matched = applied["expression_status"].notna()
    n_matched = int(matched.sum())
    n_absent = int((~matched).sum())
    n_concordant = int((applied["concordance"] == 1).sum())
    n_reversed = int((applied["concordance"] == -1).sum())
    # Not significant = matched but concordance == 0
    n_not_significant = int(n_matched - n_concordant - n_reversed)
    hit_rate = n_concordant / n_matched if n_matched > 0 else np.nan

    return {
        "rank_score": score,
        "empirical_p": perm_result["empirical_p"],
        "n_permutations": perm_result["n_permutations"],
        "n_signature_genes": len(signature_df),
        "n_matched": n_matched,
        "n_absent": n_absent,
        "n_concordant": n_concordant,
        "n_reversed": n_reversed,
        "n_not_significant": n_not_significant,
        "hit_rate": hit_rate,
    }
