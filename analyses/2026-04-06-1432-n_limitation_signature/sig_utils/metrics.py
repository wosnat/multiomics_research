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
        pd.notna(m["dir_rank"]),
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
    ref = _merge_signature_with_de(signature_df, de_ref_df)
    target = _merge_signature_with_de(signature_df, de_target_df)

    # Rename to avoid collisions
    ref = ref.rename(columns={
        "rank_up": "ref_rank_up", "rank_down": "ref_rank_down",
        "expression_status": "ref_status", "log2fc": "ref_log2fc", "rank": "ref_rank",
    })
    target = target.rename(columns={
        "rank_up": "tgt_rank_up", "rank_down": "tgt_rank_down",
        "expression_status": "tgt_status", "log2fc": "tgt_log2fc", "rank": "tgt_rank",
    })

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
