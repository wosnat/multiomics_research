"""Step 4 scoring — D8 methodology.

Per (experiment × timepoint × ontology), combine the UP and DOWN gene-list
clusters into a single signature-engagement row. For each signature term,
take the max-absolute-value of the sign-matched contribution across the UP
and DOWN clusters (bigger-magnitude evidence wins, regardless of which
cluster carried it). Apply a significance threshold `|contrib| >= 1.301`
(corresponds to padj < 0.05) so subthreshold noise contributes 0. Mean the
max-abs contributions separately across up-expected and down-expected
signature terms to get `score_up` and `score_down`; `final_score` is the
n-weighted mean.

Wires D1 (temporal filter > 3.0), D4 (NC calibration exclusion padj<1e-3),
D6 (kegg redundancy sensitivity), D7 (3-category pre-registration:
coculture / axenic_dying / axenic_dead), D8 (this scoring methodology).

Per (category × ontology × omics × direction) means aggregate per-TP rows.
D7 formal check aggregates category means across both RNA and Prot
(omics-agnostic summary); per-omics breakdowns are reported for the Task 11
decide interpretation.
"""
from __future__ import annotations

import logging
from pathlib import Path

import numpy as np
import pandas as pd

from signature import (
    CORE_SUPPORT,
    FALLBACK_SUPPORT,
    TIMEPOINT_HOURS_CUTOFF,
    derive_for_ontology,
)

ANALYSIS_DIR = Path(__file__).resolve().parent.parent

SCORE_CAP = 10.0
SIGNIFICANCE_THRESHOLD = -np.log10(0.05)  # ≈ 1.301
NC_KEY_PATHWAY_EXCLUSION_PADJ = 1e-3  # D4

# Ontology-level term totals (MED4 landscape, Step 1b). Used only for the
# whole-ontology QA column (n_total_terms).
ONTOLOGY_TOTAL_TERMS: dict[str, int] = {
    "cyanorak_role": 69,  # cyanorak L1
    "kegg": 97,           # kegg L2
}

# D7 categories, renamed per Task 10 methodology walkthrough (D8).
CATEGORY_COCULTURE = "coculture"
CATEGORY_AXENIC_DYING = "axenic_dying"
CATEGORY_AXENIC_DEAD = "axenic_dead"

# D6 Part 1: kegg redundancy sensitivity variants.
KEGG_REDUNDANCY_VARIANTS: dict[str, set[str]] = {
    "full": set(),
    "no_ko00710": {"kegg.pathway:ko00710"},
    "no_ko00710_ko00195": {"kegg.pathway:ko00710", "kegg.pathway:ko00195"},
}

# D7 pre-registration claim identifiers (see decisions.md D7 for text).
PREREG_ORDERING_CLAIMS = [
    ("P1", "cyanorak_role", "up"),
    ("P2", "cyanorak_role", "down"),
    ("P3", "kegg", "up"),
    ("P4", "kegg", "down"),
]
PREREG_THRESHOLD_CLAIMS = [
    ("T1", "cyanorak_role", "up"),
    ("T2", "cyanorak_role", "down"),
    ("T3", "kegg", "up"),
    ("T4", "kegg", "down"),
]


# --------------------------------------------------------------------------
# Primitive: thresholded + capped + sign-matched contribution
# --------------------------------------------------------------------------

def thresholded_capped_contrib(signed_score: float, sign_ref: float) -> float:
    """Per-term, per-cluster contribution under D8 scoring rule.

    signed_score: raw -log10(p_adjust) × sign(cluster_direction) from enrichment.
    sign_ref:     +1 if signature expects "up", -1 if "down".

    Returns:
        sign_ref × capped_signed_score  if |capped| >= SIGNIFICANCE_THRESHOLD
        0                               otherwise (subthreshold is noise)
        0                               if signed_score is NaN (term not tested)
    """
    if pd.isna(signed_score):
        return 0.0
    capped = float(np.sign(signed_score)) * min(abs(signed_score), SCORE_CAP)
    if abs(capped) < SIGNIFICANCE_THRESHOLD:
        return 0.0
    return sign_ref * capped


# --------------------------------------------------------------------------
# Main primitive: compute_signature_score for one (exp × tp × ontology)
# --------------------------------------------------------------------------

def compute_signature_score(
    up_cluster_enrichment: pd.DataFrame,
    down_cluster_enrichment: pd.DataFrame,
    signature_ontology: pd.DataFrame,
    n_total_terms_ontology: int,
) -> dict:
    """Combine two clusters' enrichment evidence into a single signature-score row.

    For each signature term, compute per-cluster contribution
    (sign-matched + thresholded + capped), then take max-|abs| across the
    UP and DOWN clusters. Mean by direction → score_up, score_down.

    Args:
        up_cluster_enrichment: enrichment_all.csv rows for the UP cluster on
            this ontology. May be empty.
        down_cluster_enrichment: enrichment_all.csv rows for the DOWN cluster.
        signature_ontology: signature rows for THIS ontology only.
        n_total_terms_ontology: total ontology term count (69 cyan, 97 kegg).

    Returns a dict with all schema columns (without category/omics/tp/ontology —
    caller adds those).
    """
    up_contribs: list[float] = []
    down_contribs: list[float] = []
    n_up_sig = 0
    n_down_sig = 0

    for _, sigrow in signature_ontology.iterrows():
        term_id = sigrow["term_id"]
        ref_dir = sigrow["direction"]
        sign_ref = 1.0 if ref_dir == "up" else -1.0

        up_match = up_cluster_enrichment[up_cluster_enrichment["term_id"] == term_id]
        up_signed = up_match.iloc[0]["signed_score"] if not up_match.empty else np.nan
        up_c = thresholded_capped_contrib(up_signed, sign_ref)

        down_match = down_cluster_enrichment[down_cluster_enrichment["term_id"] == term_id]
        down_signed = down_match.iloc[0]["signed_score"] if not down_match.empty else np.nan
        down_c = thresholded_capped_contrib(down_signed, sign_ref)

        # Max-abs pick across clusters.
        max_abs = up_c if abs(up_c) >= abs(down_c) else down_c

        if ref_dir == "up":
            up_contribs.append(max_abs)
            if max_abs != 0:
                n_up_sig += 1
        else:
            down_contribs.append(max_abs)
            if max_abs != 0:
                n_down_sig += 1

    n_up_total = len(up_contribs)
    n_down_total = len(down_contribs)
    score_up = float(np.mean(up_contribs)) if up_contribs else np.nan
    score_down = float(np.mean(down_contribs)) if down_contribs else np.nan

    n_terms_sig = (n_up_total + n_down_total)
    if n_terms_sig > 0:
        final_score = (
            (n_up_total * (score_up if not pd.isna(score_up) else 0.0)
             + n_down_total * (score_down if not pd.isna(score_down) else 0.0))
            / n_terms_sig
        )
    else:
        final_score = np.nan

    # Whole-ontology QA: count ALL ontology terms significant at padj<0.05
    # across BOTH clusters (sum, can double-count a term that's sig in both).
    n_total_sig_up = (
        int((up_cluster_enrichment["p_adjust"] < 0.05).sum())
        if len(up_cluster_enrichment) else 0
    )
    n_total_sig_down = (
        int((down_cluster_enrichment["p_adjust"] < 0.05).sum())
        if len(down_cluster_enrichment) else 0
    )
    n_total_sig = n_total_sig_up + n_total_sig_down

    return {
        "score_up": score_up,
        "score_down": score_down,
        "final_score": final_score,
        "n_up_sig": n_up_sig,
        "n_down_sig": n_down_sig,
        "n_up_total": n_up_total,
        "n_down_total": n_down_total,
        "n_total_sig": n_total_sig,
        "n_total_terms": n_total_terms_ontology,
    }


# --------------------------------------------------------------------------
# Category + omics assignment
# --------------------------------------------------------------------------

def assign_category(experiment_id: str, timepoint: str | None, class_label: str) -> str:
    """D7 categories for T experiments; empty for non-T."""
    if class_label != "T":
        return ""
    eid = experiment_id.lower() if isinstance(experiment_id, str) else ""
    tp = timepoint if isinstance(timepoint, str) else ""
    if "coculture" in eid:
        return CATEGORY_COCULTURE
    if "axenic" in eid:
        if tp == "day 14":
            return CATEGORY_AXENIC_DYING
        if tp in ("day 31", "day 89"):
            return CATEGORY_AXENIC_DEAD
    return ""


def assign_omics(experiment_id: str) -> str:
    eid = experiment_id.lower() if isinstance(experiment_id, str) else ""
    if "proteomics" in eid:
        return "Prot"
    if "rnaseq" in eid:
        return "RNA"
    if "microarray" in eid:
        return "Microarray"
    return ""


# --------------------------------------------------------------------------
# Loop: score all (exp × tp × ontology) rows
# --------------------------------------------------------------------------

def score_all_exp_tp(
    enrichment_df: pd.DataFrame,
    signature_df: pd.DataFrame,
    class_map: dict[str, str],
) -> pd.DataFrame:
    """For every (experiment × timepoint × ontology), find UP + DOWN clusters
    and compute the signature-score row via compute_signature_score.

    Returns DataFrame with full schema: category, omics, timepoint, ontology,
    experiment_id, background_used, score_up, score_down, final_score, n_*.
    """
    rows = []
    # Group by (experiment, timepoint, ontology, background). One iteration
    # per atomic unit; within that, up_cluster + down_cluster subsets.
    group_cols = ["experiment_id", "timepoint", "ontology", "background_used"]
    for keys, grp in enrichment_df.groupby(group_cols, dropna=False):
        experiment_id, timepoint, ontology, background_used = keys
        up_rows = grp[grp["direction"] == "up"]
        down_rows = grp[grp["direction"] == "down"]

        ont_sig = signature_df[signature_df["ontology"] == ontology]
        if ont_sig.empty:
            continue

        n_total = ONTOLOGY_TOTAL_TERMS.get(ontology, 0)
        res = compute_signature_score(up_rows, down_rows, ont_sig, n_total)

        class_label = class_map.get(experiment_id, "")
        first = grp.iloc[0]
        rows.append({
            "category": assign_category(experiment_id, timepoint, class_label),
            "omics": assign_omics(experiment_id),
            "timepoint": timepoint,
            "timepoint_hours": first.get("timepoint_hours"),
            "ontology": ontology,
            "experiment_id": experiment_id,
            "class": class_label,
            "background_used": background_used,
            "organism": first.get("organism"),
            **res,
        })
    return pd.DataFrame(rows)


# --------------------------------------------------------------------------
# D4 NC calibration exclusion (unchanged from previous revision)
# --------------------------------------------------------------------------

def disqualified_nc_clusters(
    enrichment_df: pd.DataFrame,
    key_panel: pd.DataFrame,
    nc_clusters_df: pd.DataFrame,
) -> set[tuple[str, str, str]]:
    """Return (cluster, ontology, background_used) tuples to exclude from NC
    calibration per D4 (option b, padj<1e-3 criterion). Note: exclusions
    remain at per-cluster granularity since D4 was defined before the
    per-(exp × tp) aggregation; we collapse them to per-(exp × tp) below.
    """
    disqualified: set[tuple[str, str, str]] = set()
    for ontology, kp in key_panel.groupby("ontology"):
        key_ids = set(kp["term_id"])
        nc_match = enrichment_df[
            enrichment_df["cluster"].isin(nc_clusters_df["cluster"])
            & (enrichment_df["ontology"] == ontology)
            & enrichment_df["term_id"].isin(key_ids)
            & (enrichment_df["p_adjust"] < NC_KEY_PATHWAY_EXCLUSION_PADJ)
        ]
        for (cluster, bg), _ in nc_match.groupby(["cluster", "background_used"]):
            disqualified.add((cluster, ontology, bg))
    return disqualified


def disqualified_nc_exp_tp(
    nc_cluster_exclusions: set[tuple[str, str, str]],
    enrichment_df: pd.DataFrame,
) -> set[tuple[str, str, str, str]]:
    """Collapse per-cluster D4 exclusions to per-(exp × tp × ontology × bg)
    exclusions. If EITHER the up cluster or down cluster of an (exp, tp) is
    disqualified for an (ontology, bg), exclude that (exp, tp, ontology, bg)
    from NC calibration.
    """
    # Build cluster → (experiment_id, timepoint) map. Cluster is globally
    # unique per (experiment_id × timepoint × direction); we want exp×tp, so
    # drop_duplicates on cluster alone.
    cluster_idx = (
        enrichment_df[["cluster", "experiment_id", "timepoint"]]
        .drop_duplicates(subset=["cluster"])
        .set_index("cluster")
    )
    exp_tp_excl: set[tuple[str, str, str, str]] = set()
    for (cluster, ontology, bg) in nc_cluster_exclusions:
        if cluster in cluster_idx.index:
            row = cluster_idx.loc[cluster]
            exp_tp_excl.add((row["experiment_id"], row["timepoint"], ontology, bg))
    return exp_tp_excl


# --------------------------------------------------------------------------
# NC calibration + classification on `final_score`
# --------------------------------------------------------------------------

def compute_calibration(
    scores_df: pd.DataFrame,
    nc_exclusions: set[tuple[str, str, str, str]] | None = None,
    value_col: str = "final_score",
) -> pd.DataFrame:
    """NC noise floor + PC peak per (ontology, background_used) on
    `value_col` (final_score, score_up, or score_down). Applies D4-derived
    per-(exp × tp × ontology × bg) exclusions.
    """
    nc_exclusions = nc_exclusions or set()
    nc = scores_df[scores_df["class"] == "NC"].copy()
    if nc_exclusions:
        nc["_excl_key"] = list(zip(
            nc["experiment_id"], nc["timepoint"], nc["ontology"], nc["background_used"],
        ))
        nc = nc[~nc["_excl_key"].isin(nc_exclusions)].drop(columns=["_excl_key"])
    pc = scores_df[scores_df["class"] == "PC"]
    rows = []
    groups = set(zip(nc["ontology"], nc["background_used"])) | set(
        zip(pc["ontology"], pc["background_used"])
    )
    for ont, bg in sorted(groups):
        nc_grp = nc[(nc["ontology"] == ont) & (nc["background_used"] == bg)]
        pc_grp = pc[(pc["ontology"] == ont) & (pc["background_used"] == bg)]
        rows.append({
            "ontology": ont,
            "background_used": bg,
            "value_col": value_col,
            "nc_mean": nc_grp[value_col].mean() if len(nc_grp) else np.nan,
            "nc_std": nc_grp[value_col].std(ddof=1) if len(nc_grp) >= 2 else np.nan,
            "nc_n": len(nc_grp),
            "pc_mean": pc_grp[value_col].mean() if len(pc_grp) else np.nan,
            "pc_std": pc_grp[value_col].std(ddof=1) if len(pc_grp) >= 2 else np.nan,
            "pc_n": len(pc_grp),
        })
    return pd.DataFrame(rows)


def classify(score: float, calib_df: pd.DataFrame, ontology: str, bg: str) -> str:
    c = calib_df[(calib_df["ontology"] == ontology) & (calib_df["background_used"] == bg)]
    if c.empty:
        return "no_matched_calibration"
    row = c.iloc[0]
    if pd.isna(row["nc_std"]) or row["nc_std"] == 0:
        return "insufficient_nc"
    if pd.isna(score):
        return "no_signature_coverage"
    if score >= row["nc_mean"] + 2 * row["nc_std"]:
        if (
            not pd.isna(row["pc_mean"])
            and not pd.isna(row["pc_std"])
            and score >= row["pc_mean"] - 2 * row["pc_std"]
        ):
            return "pc_like"
        return "detectable"
    if abs(score - row["nc_mean"]) < 2 * row["nc_std"]:
        return "no_signal"
    return "below_nc"


# --------------------------------------------------------------------------
# D7 category-mean aggregation + pre-registration check
# --------------------------------------------------------------------------

def compute_category_means(scores_df: pd.DataFrame) -> pd.DataFrame:
    """Aggregate per-TP scores to per (category × omics × ontology × direction).
    Also computes omics-agnostic per (category × ontology × direction) for D7
    formal check (marked omics='ALL').
    """
    t = scores_df[scores_df["category"] != ""].copy()
    rows = []
    # Per-omics breakdown (primary interpretation).
    for (cat, omics, ont), grp in t.groupby(["category", "omics", "ontology"]):
        rows.append({
            "category": cat, "omics": omics, "ontology": ont, "direction": "up",
            "mean": grp["score_up"].mean(),
            "n_tp": len(grp),
        })
        rows.append({
            "category": cat, "omics": omics, "ontology": ont, "direction": "down",
            "mean": grp["score_down"].mean(),
            "n_tp": len(grp),
        })
        rows.append({
            "category": cat, "omics": omics, "ontology": ont, "direction": "final",
            "mean": grp["final_score"].mean(),
            "n_tp": len(grp),
        })
    # Omics-agnostic (D7 formal check level).
    for (cat, ont), grp in t.groupby(["category", "ontology"]):
        rows.append({
            "category": cat, "omics": "ALL", "ontology": ont, "direction": "up",
            "mean": grp["score_up"].mean(),
            "n_tp": len(grp),
        })
        rows.append({
            "category": cat, "omics": "ALL", "ontology": ont, "direction": "down",
            "mean": grp["score_down"].mean(),
            "n_tp": len(grp),
        })
        rows.append({
            "category": cat, "omics": "ALL", "ontology": ont, "direction": "final",
            "mean": grp["final_score"].mean(),
            "n_tp": len(grp),
        })
    return pd.DataFrame(rows)


def check_preregistration(
    category_means: pd.DataFrame,
    calib_up: pd.DataFrame,
    calib_down: pd.DataFrame,
) -> pd.DataFrame:
    """Evaluate D7 P1-P4 and T1-T4 against category means (omics='ALL')."""
    rows = []

    def cat_mean(cat, ont, direction):
        m = category_means[
            (category_means["category"] == cat)
            & (category_means["omics"] == "ALL")
            & (category_means["ontology"] == ont)
            & (category_means["direction"] == direction)
        ]
        return float(m.iloc[0]["mean"]) if not m.empty else np.nan

    def threshold(ont, direction, bg="table_scope"):
        calib = calib_up if direction == "up" else calib_down
        c = calib[(calib["ontology"] == ont) & (calib["background_used"] == bg)]
        if c.empty or pd.isna(c.iloc[0]["nc_std"]):
            return np.nan, np.nan, np.nan
        return (
            float(c.iloc[0]["nc_mean"]),
            float(c.iloc[0]["nc_std"]),
            float(c.iloc[0]["nc_mean"] + 2 * c.iloc[0]["nc_std"]),
        )

    # Ordering claims P1-P4: B (axenic_dying) > A (coculture)
    for claim_id, ont, direction in PREREG_ORDERING_CLAIMS:
        b_mean = cat_mean(CATEGORY_AXENIC_DYING, ont, direction)
        a_mean = cat_mean(CATEGORY_COCULTURE, ont, direction)
        holds = (b_mean > a_mean) if not (pd.isna(b_mean) or pd.isna(a_mean)) else False
        rows.append({
            "claim_id": claim_id,
            "claim_type": "ordering_B_gt_A",
            "ontology": ont,
            "direction": direction,
            "cat_B_mean": b_mean,
            "cat_A_mean": a_mean,
            "delta_B_minus_A": (b_mean - a_mean) if not (pd.isna(b_mean) or pd.isna(a_mean)) else np.nan,
            "nc_mean": np.nan, "nc_std": np.nan, "threshold": np.nan,
            "holds": holds,
        })

    # Threshold claims T1-T4: A > nc_mean + 2σ
    for claim_id, ont, direction in PREREG_THRESHOLD_CLAIMS:
        a_mean = cat_mean(CATEGORY_COCULTURE, ont, direction)
        nc_mean, nc_std, th = threshold(ont, direction)
        holds = (a_mean > th) if not (pd.isna(a_mean) or pd.isna(th)) else False
        rows.append({
            "claim_id": claim_id,
            "claim_type": "threshold_A_above_noise",
            "ontology": ont,
            "direction": direction,
            "cat_B_mean": np.nan,
            "cat_A_mean": a_mean,
            "delta_B_minus_A": np.nan,
            "nc_mean": nc_mean, "nc_std": nc_std, "threshold": th,
            "holds": holds,
        })

    return pd.DataFrame(rows)


# --------------------------------------------------------------------------
# D6 kegg redundancy sensitivity
# --------------------------------------------------------------------------

def kegg_redundancy_sensitivity(
    scores_df: pd.DataFrame,
    enrichment_df: pd.DataFrame,
    signature: pd.DataFrame,
) -> pd.DataFrame:
    """For each T × kegg (exp × tp), recompute signature scores with
    reduced kegg signatures: full / no_ko00710 / no_ko00710_ko00195.
    Report final_score per variant and flag changes > 50% relative.
    """
    rows = []
    t_kegg = scores_df[(scores_df["class"] == "T") & (scores_df["ontology"] == "kegg")]
    full_kegg_sig = signature[signature["ontology"] == "kegg"]
    for _, trow in t_kegg.iterrows():
        up_rows = enrichment_df[
            (enrichment_df["experiment_id"] == trow["experiment_id"])
            & (enrichment_df["timepoint"] == trow["timepoint"])
            & (enrichment_df["ontology"] == "kegg")
            & (enrichment_df["direction"] == "up")
        ]
        down_rows = enrichment_df[
            (enrichment_df["experiment_id"] == trow["experiment_id"])
            & (enrichment_df["timepoint"] == trow["timepoint"])
            & (enrichment_df["ontology"] == "kegg")
            & (enrichment_df["direction"] == "down")
        ]
        variant_scores: dict[str, float] = {}
        for vname, drop_ids in KEGG_REDUNDANCY_VARIANTS.items():
            sig_variant = full_kegg_sig[~full_kegg_sig["term_id"].isin(drop_ids)]
            res = compute_signature_score(
                up_rows, down_rows, sig_variant, ONTOLOGY_TOTAL_TERMS["kegg"],
            )
            variant_scores[vname] = res["final_score"]
        entry = {
            "experiment_id": trow["experiment_id"],
            "timepoint": trow["timepoint"],
            "category": trow["category"],
            "omics": trow["omics"],
            "background_used": trow["background_used"],
        }
        for vname, s in variant_scores.items():
            entry[f"final_score_{vname}"] = s
        # Relative change flags (>50% change treated as material shift).
        full = variant_scores["full"]
        if not pd.isna(full) and abs(full) > 1e-9:
            entry["rel_change_no_ko00710"] = (variant_scores["no_ko00710"] - full) / abs(full)
            entry["rel_change_no_ko00710_ko00195"] = (variant_scores["no_ko00710_ko00195"] - full) / abs(full)
        else:
            entry["rel_change_no_ko00710"] = np.nan
            entry["rel_change_no_ko00710_ko00195"] = np.nan
        entry["material_shift_no_ko00710"] = (
            abs(entry["rel_change_no_ko00710"]) > 0.5
            if not pd.isna(entry["rel_change_no_ko00710"]) else False
        )
        entry["material_shift_no_ko00710_ko00195"] = (
            abs(entry["rel_change_no_ko00710_ko00195"]) > 0.5
            if not pd.isna(entry["rel_change_no_ko00710_ko00195"]) else False
        )
        rows.append(entry)
    return pd.DataFrame(rows)


# --------------------------------------------------------------------------
# LOO-R re-derivation (D1: hours > cutoff strict)
# --------------------------------------------------------------------------

def rederive_signature_loo(
    enrichment_df: pd.DataFrame,
    r_exp_ids: set[str],
    excluded: str,
) -> pd.DataFrame:
    r_subset = enrichment_df[
        enrichment_df["experiment_id"].isin(r_exp_ids - {excluded})
    ].copy()
    if "timepoint_hours" in r_subset.columns:
        hours = r_subset["timepoint_hours"]
        r_subset = r_subset[hours.isna() | (hours > TIMEPOINT_HOURS_CUTOFF)]  # D1
    parts = []
    for ont, ont_df in r_subset.groupby("ontology"):
        s, _ = derive_for_ontology(ont_df, CORE_SUPPORT)
        if len(s) < 5:
            s, _ = derive_for_ontology(ont_df, FALLBACK_SUPPORT)
        s["ontology"] = ont
        parts.append(s)
    return pd.concat(parts, ignore_index=True) if parts else pd.DataFrame()


# --------------------------------------------------------------------------
# Main
# --------------------------------------------------------------------------

def main() -> int:
    logging.basicConfig(
        filename=str(ANALYSIS_DIR / "logs" / "step4.log"),
        level=logging.INFO,
        format="%(asctime)s [%(levelname)s] %(message)s",
    )
    log = logging.getLogger("step4")

    enrichment = pd.read_csv(ANALYSIS_DIR / "data" / "enrichment_all.csv")
    classified = pd.read_csv(ANALYSIS_DIR / "data" / "experiments_classified.csv")
    signature = pd.read_csv(ANALYSIS_DIR / "data" / "reference_signature.csv")
    key_panel = pd.read_csv(ANALYSIS_DIR / "exploration" / "key_pathways.csv")

    classified_unique = classified.drop_duplicates("experiment_id")
    class_map = classified_unique.set_index("experiment_id")["class"].to_dict()

    # ---- Per (exp × tp × ontology) scoring via compute_signature_score ----
    scores = score_all_exp_tp(enrichment, signature, class_map)
    scores.to_csv(ANALYSIS_DIR / "results" / "scores_all.csv", index=False)
    log.info(f"Wrote {len(scores)} per-(exp × tp × ontology) signature-score rows")

    # ---- D4 NC exclusions (cluster-granular → exp×tp granular) ----
    nc_clusters = (
        enrichment[enrichment["experiment_id"].isin(
            {e for e, c in class_map.items() if c == "NC"}
        )][["cluster"]]
        .drop_duplicates()
    )
    nc_cluster_exclusions = disqualified_nc_clusters(enrichment, key_panel, nc_clusters)
    nc_exp_tp_exclusions = disqualified_nc_exp_tp(nc_cluster_exclusions, enrichment)
    for (excl_key) in sorted(nc_cluster_exclusions):
        log.info(f"NC cluster exclusion (D4): {excl_key}")
    for excl_key in sorted(nc_exp_tp_exclusions):
        log.info(f"NC exp×tp exclusion (D4 collapsed): {excl_key}")
    if nc_exp_tp_exclusions:
        pd.DataFrame(
            sorted(nc_exp_tp_exclusions),
            columns=["experiment_id", "timepoint", "ontology", "background_used"],
        ).to_csv(
            ANALYSIS_DIR / "results" / "nc_calibration_exclusions.csv",
            index=False,
        )

    # ---- Calibration on final_score + per-direction ----
    calib_final = compute_calibration(scores, nc_exp_tp_exclusions, value_col="final_score")
    calib_up = compute_calibration(scores, nc_exp_tp_exclusions, value_col="score_up")
    calib_down = compute_calibration(scores, nc_exp_tp_exclusions, value_col="score_down")
    calibs = pd.concat([calib_final, calib_up, calib_down], ignore_index=True)
    calibs.to_csv(ANALYSIS_DIR / "results" / "nc_calibration.csv", index=False)
    log.info(
        f"Calibrations: final/up/down × {len(calib_final)} (ont×bg) groups; "
        f"NC exclusions applied: {len(nc_exp_tp_exclusions)}"
    )

    # ---- Classify T rows on final_score ----
    t_rows = scores[scores["class"] == "T"].copy()
    t_rows["classification"] = [
        classify(row["final_score"], calib_final, row["ontology"], row["background_used"])
        for _, row in t_rows.iterrows()
    ]
    t_rows.to_csv(ANALYSIS_DIR / "results" / "score_summary.csv", index=False)

    # ---- D7 category means + pre-reg check ----
    category_means = compute_category_means(scores)
    category_means.to_csv(
        ANALYSIS_DIR / "results" / "category_mean_scores.csv", index=False
    )

    prereg = check_preregistration(category_means, calib_up, calib_down)
    prereg.to_csv(ANALYSIS_DIR / "results" / "preregistration_check.csv", index=False)
    for _, pr in prereg.iterrows():
        log.info(
            f"Prereg {pr['claim_id']} ({pr['claim_type']}, {pr['ontology']}×{pr['direction']}): "
            f"holds={pr['holds']}"
        )

    # ---- D6 kegg redundancy sensitivity ----
    kegg_red = kegg_redundancy_sensitivity(t_rows, enrichment, signature)
    kegg_red.to_csv(
        ANALYSIS_DIR / "results" / "kegg_redundancy_sensitivity.csv", index=False
    )
    shifts_710 = int(kegg_red["material_shift_no_ko00710"].sum()) if not kegg_red.empty else 0
    shifts_710_195 = int(kegg_red["material_shift_no_ko00710_ko00195"].sum()) if not kegg_red.empty else 0
    log.info(
        f"D6 kegg redundancy: {len(kegg_red)} T-kegg (exp×tp); "
        f"material shifts >50%: no_ko00710={shifts_710}, no_ko00710_ko00195={shifts_710_195}"
    )

    # ---- LOO signature (per T × ontology × signature term removed) ----
    loo_sig_rows = []
    for _, trow in t_rows.iterrows():
        ont_sig = signature[signature["ontology"] == trow["ontology"]]
        up_rows = enrichment[
            (enrichment["experiment_id"] == trow["experiment_id"])
            & (enrichment["timepoint"] == trow["timepoint"])
            & (enrichment["ontology"] == trow["ontology"])
            & (enrichment["direction"] == "up")
        ]
        down_rows = enrichment[
            (enrichment["experiment_id"] == trow["experiment_id"])
            & (enrichment["timepoint"] == trow["timepoint"])
            & (enrichment["ontology"] == trow["ontology"])
            & (enrichment["direction"] == "down")
        ]
        orig = trow["final_score"]
        for _, sigrow in ont_sig.iterrows():
            reduced = ont_sig[ont_sig["term_id"] != sigrow["term_id"]]
            res = compute_signature_score(up_rows, down_rows, reduced, ONTOLOGY_TOTAL_TERMS[trow["ontology"]])
            new = res["final_score"]
            loo_sig_rows.append({
                "experiment_id": trow["experiment_id"],
                "timepoint": trow["timepoint"],
                "category": trow["category"],
                "omics": trow["omics"],
                "ontology": trow["ontology"],
                "background_used": trow["background_used"],
                "removed_term_id": sigrow["term_id"],
                "removed_term_name": sigrow["term_name"],
                "orig_final_score": orig,
                "new_final_score": new,
                "flag_sign_flip": (np.sign(orig) != np.sign(new)) if not pd.isna(new) else False,
                "flag_large_drop": (abs(new) < 0.5 * abs(orig)) if not pd.isna(new) and not pd.isna(orig) and abs(orig) > 1e-9 else False,
            })
    loo_sig = pd.DataFrame(loo_sig_rows)
    loo_sig.to_csv(ANALYSIS_DIR / "results" / "loo_signature.csv", index=False)
    log.info(f"LOO-signature: {len(loo_sig)} rows")

    # ---- LOO-R (classification flip) ----
    r_exp_ids = {eid for eid, c in class_map.items() if c == "R"}
    orig_class_lookup = {
        (r["experiment_id"], r["timepoint"], r["ontology"], r["background_used"]): r["classification"]
        for _, r in t_rows.iterrows()
    }

    loo_r_rows = []
    for excl in sorted(r_exp_ids):
        new_sig = rederive_signature_loo(enrichment, r_exp_ids, excl)
        if new_sig.empty:
            log.warning(f"LOO-R excl={excl}: empty signature")
            for _, trow in t_rows.iterrows():
                key = (trow["experiment_id"], trow["timepoint"], trow["ontology"], trow["background_used"])
                loo_r_rows.append({
                    "excluded_experiment": excl,
                    **{k: trow[k] for k in ("experiment_id", "timepoint", "category", "omics", "ontology", "background_used")},
                    "orig_final_score": trow["final_score"],
                    "new_final_score": np.nan,
                    "orig_classification": orig_class_lookup[key],
                    "new_classification": "empty_signature",
                    "classification_flip": True,
                    "new_signature_size": 0,
                })
            continue

        new_scores = score_all_exp_tp(enrichment, new_sig, class_map)
        new_calib = compute_calibration(new_scores, nc_exp_tp_exclusions, value_col="final_score")
        new_sig_size_by_ont = new_sig.groupby("ontology").size().to_dict()

        for _, trow in t_rows.iterrows():
            key = (trow["experiment_id"], trow["timepoint"], trow["ontology"], trow["background_used"])
            orig_class = orig_class_lookup[key]
            match = new_scores[
                (new_scores["experiment_id"] == trow["experiment_id"])
                & (new_scores["timepoint"] == trow["timepoint"])
                & (new_scores["ontology"] == trow["ontology"])
                & (new_scores["background_used"] == trow["background_used"])
            ]
            if match.empty:
                new_final = np.nan
                new_class = "no_signature_coverage"
            else:
                new_final = match.iloc[0]["final_score"]
                new_class = classify(new_final, new_calib, trow["ontology"], trow["background_used"])
            # D8-i raw-score stability flags.
            orig_final = trow["final_score"]
            if pd.isna(new_final) or pd.isna(orig_final):
                flag_sign_flip_raw = False
                flag_large_change_raw = True  # missing new score is itself a large change
                rel_change = np.nan
            elif abs(orig_final) < 1e-9:
                flag_sign_flip_raw = False
                flag_large_change_raw = abs(new_final) > 0.5  # orig ~ 0, any appreciable new is "change"
                rel_change = np.nan
            else:
                flag_sign_flip_raw = np.sign(orig_final) != np.sign(new_final)
                rel_change = (new_final - orig_final) / abs(orig_final)
                flag_large_change_raw = abs(rel_change) > 0.5
            loo_r_rows.append({
                "excluded_experiment": excl,
                **{k: trow[k] for k in ("experiment_id", "timepoint", "category", "omics", "ontology", "background_used")},
                "orig_final_score": orig_final,
                "new_final_score": new_final,
                "rel_change": rel_change,
                "flag_sign_flip_raw": flag_sign_flip_raw,
                "flag_large_change_raw": flag_large_change_raw,
                "orig_classification": orig_class,
                "new_classification": new_class,
                "classification_flip": orig_class != new_class,
                "new_signature_size": int(new_sig_size_by_ont.get(trow["ontology"], 0)),
            })
    loo_r = pd.DataFrame(loo_r_rows)
    loo_r.to_csv(ANALYSIS_DIR / "results" / "loo_r_experiments.csv", index=False)
    n_class_flips = int(loo_r["classification_flip"].sum()) if not loo_r.empty else 0
    n_raw_sign_flips = int(loo_r["flag_sign_flip_raw"].sum()) if not loo_r.empty else 0
    n_raw_large_changes = int(loo_r["flag_large_change_raw"].sum()) if not loo_r.empty else 0
    log.info(
        f"LOO-R: {len(loo_r)} rows — classification_flips={n_class_flips}, "
        f"raw_sign_flips={n_raw_sign_flips}, raw_large_changes={n_raw_large_changes}"
    )

    # ---- Stdout summary ----
    print("\n=== T scores (per exp × tp × ontology) — final_score + up/down ===")
    print(t_rows[[
        "category", "omics", "timepoint", "ontology", "experiment_id",
        "score_up", "score_down", "final_score",
        "n_up_sig", "n_down_sig", "n_up_total", "n_down_total",
        "n_total_sig", "n_total_terms", "classification",
    ]].sort_values(["category", "omics", "timepoint", "ontology"]).to_string(index=False))

    print("\n=== Category means (per-omics breakdown + D7 omics='ALL') ===")
    print(category_means.sort_values(["category", "omics", "ontology", "direction"]).to_string(index=False))

    print("\n=== D7 pre-reg check ===")
    print(prereg[[
        "claim_id", "claim_type", "ontology", "direction",
        "cat_B_mean", "cat_A_mean", "threshold", "holds",
    ]].to_string(index=False))

    print(f"\n=== D6 kegg redundancy: material shifts >50%: "
          f"no_ko00710={shifts_710}, no_ko00710_ko00195={shifts_710_195} / {len(kegg_red)} ===")

    print(f"\n=== LOO-R ({len(loo_r)} rows) ===")
    print(f"  classification flips: {n_class_flips}")
    print(f"  raw sign flips:       {n_raw_sign_flips}")
    print(f"  raw large changes (>50%): {n_raw_large_changes}")
    if n_raw_large_changes or n_raw_sign_flips:
        flagged = loo_r[loo_r["flag_large_change_raw"] | loo_r["flag_sign_flip_raw"]]
        print(flagged[[
            "excluded_experiment", "category", "omics", "timepoint", "ontology",
            "orig_final_score", "new_final_score", "rel_change",
            "flag_sign_flip_raw", "flag_large_change_raw", "new_signature_size",
        ]].sort_values(["category", "omics", "timepoint", "ontology"]).to_string(index=False))

    # LOO-signature summary
    n_sig_flips = int(loo_sig["flag_sign_flip"].sum()) if not loo_sig.empty else 0
    n_sig_drops = int(loo_sig["flag_large_drop"].sum()) if not loo_sig.empty else 0
    print(f"\n=== LOO-signature ({len(loo_sig)} rows) ===")
    print(f"  sign flips:          {n_sig_flips}")
    print(f"  large drops (>50%):  {n_sig_drops}")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
