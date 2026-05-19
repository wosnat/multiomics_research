"""Statistical tests for the cyano N-gene inventory analysis.

Three tests, each addressing a specific Discussion-relevant question:

  Test 1 — Mantel test on Jaccard vs taxonomic distance.
    Does the data-driven clustering correlate with taxonomy at the
    distance-matrix level? Tests the step-3 pre-registered
    "taxa-clustering" hypothesis. Permutation-based p-value.

  Test 2 — Per (clade-group × pathway) Fisher's exact test.
    For each of the 8 clade groups × 10 pathways, ask whether at least
    one gene from that pathway is enriched (or depleted) in that group
    relative to the rest of the cyano. BH-corrected p-values.

  Test 3 — Per-strain silhouette score on the Jaccard distance using
    clade-group labels. Negative or near-zero scores flag strains that
    fit their clade poorly (we expect SS120 and WH7803 to flag).

Inputs:
    2_kg_selection/data/05_inventory_matrix.csv
    2_kg_selection/data/02_strain_table.csv
    2_kg_selection/data/04_ortholog_groups.csv
Outputs:
    data/01_mantel_test.csv          Mantel r and p-value
    data/01_fisher_pathway_clade.csv per (clade_group, pathway) p, q, OR
    data/01_silhouette.csv           per-strain silhouette + clade label
    data/01_statistical_tests.log
"""
from __future__ import annotations

import logging
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.spatial.distance import pdist, squareform
from scipy.stats import fisher_exact, mannwhitneyu
from sklearn.metrics import silhouette_samples
from statsmodels.stats.multitest import multipletests

HERE = Path(__file__).resolve()
ANALYSIS_ROOT = HERE.parents[2]
STEP_DIR = HERE.parents[1]
DATA_DIR = STEP_DIR / "data"
INVENTORY_DIR = ANALYSIS_ROOT / "2_kg_selection" / "data"
LOG_PATH = DATA_DIR / "01_statistical_tests.log"

sys.path.insert(0, str(ANALYSIS_ROOT / "4_methods"))
from n_pathway_categories import CATEGORY_ORDER, categorize  # noqa: E402

# Clade-group definitions (mirror step-5 02_pathway_views.py)
CLADE_GROUPS = {
    "Pro HLI":           ["Prochlorococcus MED4", "Prochlorococcus RSP50"],
    "Pro HLII":          ["Prochlorococcus AS9601", "Prochlorococcus MIT9301",
                          "Prochlorococcus MIT9312"],
    "Pro LLI":           ["Prochlorococcus MIT0801", "Prochlorococcus NATL1A",
                          "Prochlorococcus NATL2A"],
    "Pro LLII":          ["Prochlorococcus marinus subsp. marinus CCMP1375 (SS120)"],
    "Pro LLIV":          ["Prochlorococcus MIT9303", "Prochlorococcus MIT9313"],
    "Marine Syn":        ["Synechococcus WH8102", "Synechococcus WH7803",
                          "Synechococcus CC9311", "Synechococcus sp. BL107"],
    "Non-marine Syn":    ["Synechococcus PCC 7002",
                          "Synechococcus elongatus PCC 7942",
                          "Synechococcus elongatus UTEX 2973"],
    "Thermosyn":         ["Thermosynechococcus vestitus BP-1"],
}


def strain_to_clade_group() -> dict[str, str]:
    out = {}
    for group, strains in CLADE_GROUPS.items():
        for s in strains:
            out[s] = group
    return out


def taxonomic_distance(strain_a: str, strain_b: str, strain_table: pd.DataFrame) -> float:
    """Simple taxonomic distance: 0 same clade_canonical, 1 same genus diff clade,
    2 different genus, 3 different genus AND different parent (Pro vs Syn/Thermosyn)."""
    if strain_a == strain_b:
        return 0.0
    a = strain_table.loc[strain_a]
    b = strain_table.loc[strain_b]
    if a["genus"] == b["genus"]:
        if pd.notna(a["clade_canonical"]) and a["clade_canonical"] == b["clade_canonical"]:
            return 0.0
        return 1.0
    # different genus
    pro = {"Prochlorococcus"}
    if {a["genus"], b["genus"]}.issubset(pro):
        return 2.0
    syn_like = {"Synechococcus", "Synechococcus elongatus", "Parasynechococcus",
                "Picosynechococcus", "Thermosynechococcus"}
    if a["genus"] in pro and b["genus"] in pro:
        return 2.0
    if a["genus"] in syn_like and b["genus"] in syn_like:
        return 2.0
    return 3.0


def mantel_test(d1_condensed: np.ndarray, d2_condensed: np.ndarray,
                permutations: int = 9999, rng: np.random.Generator | None = None) -> tuple[float, float]:
    """Mantel test using Pearson correlation; permutation p-value on r."""
    if rng is None:
        rng = np.random.default_rng(0)
    d1 = squareform(d1_condensed)
    d2 = squareform(d2_condensed)
    n = d1.shape[0]
    observed = np.corrcoef(d1_condensed, d2_condensed)[0, 1]
    perms = np.empty(permutations)
    for k in range(permutations):
        perm = rng.permutation(n)
        d1p = d1[perm][:, perm]
        d1p_cond = d1p[np.triu_indices(n, k=1)]
        perms[k] = np.corrcoef(d1p_cond, d2_condensed)[0, 1]
    p = (np.sum(np.abs(perms) >= np.abs(observed)) + 1) / (permutations + 1)
    return float(observed), float(p)


def per_pathway_fisher(binary: pd.DataFrame, groups: pd.DataFrame,
                       strain_groups: dict[str, str]) -> pd.DataFrame:
    """For each (clade_group, pathway), test whether strains in that group
    are over- or under-represented in having ≥1 gene from the pathway.

    Returns a DataFrame with one row per (clade_group, pathway): odds ratio,
    raw p-value, BH-corrected q-value, count summary.
    """
    cat_map = groups.set_index("group_id")["consensus_gene_name"].apply(categorize)
    # For each strain, does it have ≥1 gene from each pathway?
    pathway_present = pd.DataFrame(0, index=binary.index, columns=CATEGORY_ORDER, dtype=int)
    for cat in CATEGORY_ORDER:
        cols = [gid for gid in binary.columns if cat_map.get(gid, "other") == cat]
        if cols:
            pathway_present[cat] = (binary[cols].sum(axis=1) > 0).astype(int)
    rows = []
    cg_list = list(CLADE_GROUPS.keys())
    for cg in cg_list:
        strains_in = [s for s in CLADE_GROUPS[cg] if s in binary.index]
        strains_out = [s for s in binary.index if s not in strains_in]
        for cat in CATEGORY_ORDER:
            in_pos = int(pathway_present.loc[strains_in, cat].sum())
            in_neg = len(strains_in) - in_pos
            out_pos = int(pathway_present.loc[strains_out, cat].sum())
            out_neg = len(strains_out) - out_pos
            table = [[in_pos, in_neg], [out_pos, out_neg]]
            try:
                or_val, p = fisher_exact(table, alternative="two-sided")
            except Exception:
                or_val, p = float("nan"), 1.0
            rows.append({
                "clade_group": cg,
                "pathway": cat,
                "n_strains_in_group": len(strains_in),
                "in_present": in_pos, "in_absent": in_neg,
                "out_present": out_pos, "out_absent": out_neg,
                "odds_ratio": or_val,
                "p_value": p,
            })
    df = pd.DataFrame(rows)
    # BH correction across all 80 tests
    valid = df["p_value"].notna()
    _, qvals, _, _ = multipletests(df.loc[valid, "p_value"].values, method="fdr_bh")
    df.loc[valid, "q_value_bh"] = qvals
    return df


def main() -> None:
    DATA_DIR.mkdir(parents=True, exist_ok=True)
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s %(levelname)s %(message)s",
        handlers=[logging.FileHandler(LOG_PATH, mode="w"), logging.StreamHandler()],
    )

    matrix = pd.read_csv(INVENTORY_DIR / "05_inventory_matrix.csv", index_col=0)
    strain_table = pd.read_csv(INVENTORY_DIR / "02_strain_table.csv").set_index("strain")
    groups = pd.read_csv(INVENTORY_DIR / "04_ortholog_groups.csv")
    binary = (matrix > 0).astype(int)
    sg_map = strain_to_clade_group()

    # ---------- Test 1: Mantel ----------
    logging.info("===== Test 1: Mantel (Jaccard vs taxonomic distance) =====")
    jaccard_d = pdist(binary.values, metric="jaccard")
    n = binary.shape[0]
    tax_d = np.zeros((n, n))
    strains = binary.index.tolist()
    for i in range(n):
        for j in range(i + 1, n):
            d = taxonomic_distance(strains[i], strains[j], strain_table)
            tax_d[i, j] = tax_d[j, i] = d
    tax_d_cond = tax_d[np.triu_indices(n, k=1)]
    rng = np.random.default_rng(0)
    r, p = mantel_test(jaccard_d, tax_d_cond, permutations=9999, rng=rng)
    logging.info("Mantel r = %.4f, p (9999 perms, two-sided) = %.4g", r, p)
    pd.DataFrame([{
        "test": "mantel_jaccard_vs_taxonomy",
        "n_strains": n,
        "n_distances": len(jaccard_d),
        "permutations": 9999,
        "mantel_r": r,
        "p_value": p,
    }]).to_csv(DATA_DIR / "01_mantel_test.csv", index=False)

    # ---------- Test 2a: Fisher's exact per (clade_group, pathway) — broad scan ----------
    logging.info("===== Test 2a: Fisher's exact (clade group × pathway), BH × 80 =====")
    fisher_df = per_pathway_fisher(binary, groups, sg_map)
    fisher_df.to_csv(DATA_DIR / "01_fisher_pathway_clade.csv", index=False)
    sig = fisher_df[fisher_df["q_value_bh"] < 0.05].copy()
    logging.info("Significant (q < 0.05): %d of %d. (Conservative due to small clade-group N.)",
                 len(sig), len(fisher_df))
    # Report raw-p < 0.05 for transparency
    raw_sig = fisher_df[fisher_df["p_value"] < 0.05].copy().sort_values("p_value")
    logging.info("Raw p < 0.05 (uncorrected, %d contrasts):", len(raw_sig))
    for _, r_ in raw_sig.iterrows():
        direction = "ENRICHED" if r_["odds_ratio"] > 1 else "DEPLETED"
        logging.info(
            "  %-18s × %-18s %s in=%d/%d out=%d/%d OR=%.2f p=%.3g q=%.3g",
            r_["clade_group"], r_["pathway"], direction,
            int(r_["in_present"]), int(r_["in_present"] + r_["in_absent"]),
            int(r_["out_present"]), int(r_["out_present"] + r_["out_absent"]),
            r_["odds_ratio"], r_["p_value"], r_["q_value_bh"],
        )

    # ---------- Test 2b: pre-specified pairwise contrasts (small, focused) ----------
    logging.info("===== Test 2b: Pre-specified contrasts (Fisher's exact, BH × 5) =====")
    cat_map = groups.set_index("group_id")["consensus_gene_name"].apply(categorize)
    pathway_present = pd.DataFrame(0, index=binary.index, columns=CATEGORY_ORDER, dtype=int)
    for cat in CATEGORY_ORDER:
        cols = [gid for gid in binary.columns if cat_map.get(gid, "other") == cat]
        if cols:
            pathway_present[cat] = (binary[cols].sum(axis=1) > 0).astype(int)

    contrasts = [
        # (name, pathway, group_a_strains, group_b_strains)
        ("HLI Pro vs other Pro: cyanate",
         "cyanate", CLADE_GROUPS["Pro HLI"],
         [s for cg, ss in CLADE_GROUPS.items() if cg.startswith("Pro") and cg != "Pro HLI" for s in ss]),
        ("Syn (all) vs Pro: nitrate",
         "nitrate_nitrite",
         CLADE_GROUPS["Marine Syn"] + CLADE_GROUPS["Non-marine Syn"] + CLADE_GROUPS["Thermosyn"],
         [s for cg in ["Pro HLI", "Pro HLII", "Pro LLI", "Pro LLII", "Pro LLIV"] for s in CLADE_GROUPS[cg]]),
        ("Syn (all) vs Pro: Mo cofactor",
         "mo_cofactor",
         CLADE_GROUPS["Marine Syn"] + CLADE_GROUPS["Non-marine Syn"] + CLADE_GROUPS["Thermosyn"],
         [s for cg in ["Pro HLI", "Pro HLII", "Pro LLI", "Pro LLII", "Pro LLIV"] for s in CLADE_GROUPS[cg]]),
        ("SS120 vs other Pro: urea_uptake",
         "urea_uptake", CLADE_GROUPS["Pro LLII"],
         [s for cg in ["Pro HLI", "Pro HLII", "Pro LLI", "Pro LLIV"] for s in CLADE_GROUPS[cg]]),
        ("Thermosyn vs others: met biosynthesis",
         "met_biosynthesis", CLADE_GROUPS["Thermosyn"],
         [s for cg in CLADE_GROUPS if cg != "Thermosyn" for s in CLADE_GROUPS[cg]]),
    ]
    rows = []
    for name, pathway, group_a, group_b in contrasts:
        a_pos = int(pathway_present.loc[group_a, pathway].sum())
        a_neg = len(group_a) - a_pos
        b_pos = int(pathway_present.loc[group_b, pathway].sum())
        b_neg = len(group_b) - b_pos
        or_val, p = fisher_exact([[a_pos, a_neg], [b_pos, b_neg]], alternative="two-sided")
        rows.append({
            "contrast": name, "pathway": pathway,
            "a_present": a_pos, "a_total": len(group_a),
            "b_present": b_pos, "b_total": len(group_b),
            "odds_ratio": or_val, "p_value": p,
        })
    contrast_df = pd.DataFrame(rows)
    _, qvals, _, _ = multipletests(contrast_df["p_value"].values, method="fdr_bh")
    contrast_df["q_value_bh"] = qvals
    contrast_df.to_csv(DATA_DIR / "01_fisher_prespecified_contrasts.csv", index=False)
    for _, r_ in contrast_df.iterrows():
        marker = "*" if r_["q_value_bh"] < 0.05 else " "
        logging.info(
            "  %s %-50s a=%d/%d b=%d/%d OR=%.2f p=%.3g q=%.3g",
            marker, r_["contrast"],
            int(r_["a_present"]), int(r_["a_total"]),
            int(r_["b_present"]), int(r_["b_total"]),
            r_["odds_ratio"], r_["p_value"], r_["q_value_bh"],
        )

    # ---------- Test 2c: Mann-Whitney U on % pathway coverage ----------
    logging.info("===== Test 2c: Mann-Whitney U on %% pathway coverage =====")
    # Build the %-coverage matrix (strains × 10 pathways) — same shape as Figure 2 input
    pct_present = pd.DataFrame(0.0, index=binary.index, columns=CATEGORY_ORDER)
    for cat in CATEGORY_ORDER:
        cols = [gid for gid in binary.columns if cat_map.get(gid, "other") == cat]
        if cols:
            pct_present[cat] = binary[cols].sum(axis=1) / len(cols)
    mwu_rows = []
    mwu_contrasts = [
        ("Syn vs Pro: nitrate_nitrite",
         CLADE_GROUPS["Marine Syn"] + CLADE_GROUPS["Non-marine Syn"] + CLADE_GROUPS["Thermosyn"],
         [s for cg in ["Pro HLI", "Pro HLII", "Pro LLI", "Pro LLII", "Pro LLIV"] for s in CLADE_GROUPS[cg]],
         "nitrate_nitrite"),
        ("Syn vs Pro: mo_cofactor",
         CLADE_GROUPS["Marine Syn"] + CLADE_GROUPS["Non-marine Syn"] + CLADE_GROUPS["Thermosyn"],
         [s for cg in ["Pro HLI", "Pro HLII", "Pro LLI", "Pro LLII", "Pro LLIV"] for s in CLADE_GROUPS[cg]],
         "mo_cofactor"),
        ("Marine Syn vs Non-marine Syn: met_biosynthesis",
         CLADE_GROUPS["Marine Syn"],
         CLADE_GROUPS["Non-marine Syn"] + CLADE_GROUPS["Thermosyn"],
         "met_biosynthesis"),
        ("HL Pro vs LL Pro: nitrate_nitrite",
         CLADE_GROUPS["Pro HLI"] + CLADE_GROUPS["Pro HLII"],
         CLADE_GROUPS["Pro LLI"] + CLADE_GROUPS["Pro LLII"] + CLADE_GROUPS["Pro LLIV"],
         "nitrate_nitrite"),
        ("HL Pro vs LL Pro: cyanate",
         CLADE_GROUPS["Pro HLI"] + CLADE_GROUPS["Pro HLII"],
         CLADE_GROUPS["Pro LLI"] + CLADE_GROUPS["Pro LLII"] + CLADE_GROUPS["Pro LLIV"],
         "cyanate"),
    ]
    for name, a, b, pathway in mwu_contrasts:
        a_vals = pct_present.loc[a, pathway].values
        b_vals = pct_present.loc[b, pathway].values
        try:
            stat, p = mannwhitneyu(a_vals, b_vals, alternative="two-sided")
        except ValueError:
            stat, p = float("nan"), float("nan")
        mwu_rows.append({
            "contrast": name, "pathway": pathway,
            "a_n": len(a_vals), "a_median_pct": float(np.median(a_vals)),
            "b_n": len(b_vals), "b_median_pct": float(np.median(b_vals)),
            "U_statistic": stat, "p_value": p,
        })
    mwu_df = pd.DataFrame(mwu_rows)
    _, qvals, _, _ = multipletests(mwu_df["p_value"].fillna(1).values, method="fdr_bh")
    mwu_df["q_value_bh"] = qvals
    mwu_df.to_csv(DATA_DIR / "01_mannwhitney_contrasts.csv", index=False)
    for _, r_ in mwu_df.iterrows():
        marker = "*" if r_["q_value_bh"] < 0.05 else " "
        logging.info(
            "  %s %-52s a_med=%.2f (n=%d) b_med=%.2f (n=%d) U=%.1f p=%.3g q=%.3g",
            marker, r_["contrast"],
            r_["a_median_pct"], int(r_["a_n"]), r_["b_median_pct"], int(r_["b_n"]),
            r_["U_statistic"], r_["p_value"], r_["q_value_bh"],
        )

    # ---------- Test 3: Silhouette per strain ----------
    logging.info("===== Test 3: Silhouette per strain on Jaccard distance, by clade group =====")
    labels = np.array([sg_map[s] for s in strains])
    # silhouette_samples needs ≥2 labels with ≥2 members; Pro LLII and Thermosyn have n=1.
    # Aggregate those with their closest larger group for the silhouette computation:
    # Pro LLII (SS120) → "Pro LL" pool (i.e., all LL Pro)
    # Thermosyn (BP-1) → "Other cyano" pool with Non-marine Syn
    # We compute silhouette twice: (a) raw clade groups (SS120, BP-1 get NaN), (b) merged labels.
    # Report both.
    merged_label_map = {
        "Pro HLI": "Pro HL",
        "Pro HLII": "Pro HL",
        "Pro LLI": "Pro LL",
        "Pro LLII": "Pro LL",      # SS120 merged
        "Pro LLIV": "Pro LL",
        "Marine Syn": "Marine Syn",
        "Non-marine Syn": "Non-marine cyano",
        "Thermosyn": "Non-marine cyano",  # BP-1 merged
    }
    merged_labels = np.array([merged_label_map[l] for l in labels])
    jaccard_sq = squareform(jaccard_d)
    sil_merged = silhouette_samples(jaccard_sq, merged_labels, metric="precomputed")
    out = pd.DataFrame({
        "strain": strains,
        "clade_group": labels,
        "merged_label": merged_labels,
        "silhouette_merged": sil_merged,
    }).sort_values("silhouette_merged")
    out.to_csv(DATA_DIR / "01_silhouette.csv", index=False)
    logging.info("Silhouette (merged labels): mean=%.3f, min=%.3f (= worst fit), max=%.3f",
                 sil_merged.mean(), sil_merged.min(), sil_merged.max())
    logging.info("Worst-fitting strains (lowest silhouette = farthest from own group centroid):")
    for _, r_ in out.head(5).iterrows():
        logging.info("  %-55s s=%.3f  (%s → %s)",
                     r_["strain"], r_["silhouette_merged"],
                     r_["clade_group"], r_["merged_label"])
    logging.info("Best-fitting strains (highest silhouette):")
    for _, r_ in out.tail(5).iterrows():
        logging.info("  %-55s s=%.3f  (%s → %s)",
                     r_["strain"], r_["silhouette_merged"],
                     r_["clade_group"], r_["merged_label"])

    logging.info("Wrote 3 result CSVs.")


if __name__ == "__main__":
    main()
