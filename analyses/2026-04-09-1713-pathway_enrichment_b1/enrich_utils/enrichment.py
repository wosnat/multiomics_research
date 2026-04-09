"""Fisher's exact test + BH-FDR enrichment for pathway analysis.

Pure DataFrame in / DataFrame out. No KG dependency.

Public API:
    run_enrichment(de_df, pathway_defs, gene_universe, table_scope) -> DataFrame
    run_enrichment_all_timepoints(de_df, pathway_defs, table_scope) -> DataFrame

2×2 contingency table (over-representation, one-sided):

    DE (in direction)   | Not DE (in direction)
    -------------------------------------------------
    In pathway    |  a  |  c
    Not in pathway|  b  |  d

fisher_exact([[a, b], [c, d]], alternative='greater')

Note: scipy.stats.fisher_exact takes [[row1_col1, row1_col2], [row2_col1, row2_col2]].
Here we treat rows = pathway membership, cols = DE status.
  [[a, c], [b, d]] with alternative='greater' tests whether pathway genes are
  over-represented among DE genes.
"""

import math
import warnings

import numpy as np
import pandas as pd
from scipy.stats import fisher_exact
from statsmodels.stats.multitest import multipletests

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

TABLE_SCOPE_TO_TEST_TYPE = {
    "all_detected_genes": "vs_genome",
    "filtered_subset": "vs_filtered_genome",
    "significant_any_timepoint": "vs_all_responsive",
    "significant_only": "descriptive_only",
}

UNDERPOWERED_THRESHOLD = 5  # fewer than this many pathway genes in universe → underpowered

DIRECTIONS = ("up", "down", "combined")


# ---------------------------------------------------------------------------
# Core enrichment
# ---------------------------------------------------------------------------

def _compute_2x2(pathway_in_universe: set, de_genes: set, universe: set):
    """Return (a, b, c, d) for a single pathway × direction test.

    a = pathway ∩ DE
    b = (universe − pathway) ∩ DE
    c = pathway ∩ (universe − DE)
    d = (universe − pathway) ∩ (universe − DE)
    """
    pw = pathway_in_universe
    not_pw = universe - pw
    a = len(pw & de_genes)
    c = len(pw) - a
    b = len(not_pw & de_genes)
    d = len(not_pw) - b
    return a, b, c, d


def _fisher_row(a, b, c, d):
    """Run Fisher's exact (greater) and return (odds_ratio, p_value).

    Returns (NaN, NaN) if the table is degenerate (all zeros in a row/column).
    """
    if a + c == 0:
        # No pathway genes in universe — degenerate
        return math.nan, math.nan
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        odds_ratio, p_value = fisher_exact([[a, c], [b, d]], alternative="greater")
    return float(odds_ratio), float(p_value)


def _fold_enrichment(a, b, c, d):
    """Fold enrichment = observed / expected.

    observed fraction in pathway = a / (a + c)
    expected fraction (background) = (a + b) / (a + b + c + d)
    fold_enrichment = observed / expected
    """
    n_pathway = a + c
    n_universe = a + b + c + d
    n_de = a + b
    if n_pathway == 0 or n_universe == 0 or n_de == 0:
        return math.nan
    expected_rate = n_de / n_universe
    observed_rate = a / n_pathway
    return observed_rate / expected_rate if expected_rate > 0 else math.nan


def _expected_count(a, b, c, d):
    """Expected overlap under independence = (a+c) * (a+b) / N."""
    n_pathway = a + c
    n_de = a + b
    n_total = a + b + c + d
    if n_total == 0:
        return math.nan
    return (n_pathway * n_de) / n_total


def _enrich_one_pathway_direction(
    pathway_id: str,
    pathway_name: str,
    pathway_genes_total: int,
    pathway_in_universe: set,
    de_genes: set,
    universe: set,
    direction: str,
    test_type: str,
) -> dict:
    """Compute enrichment stats for one pathway × one direction."""
    n_in_universe = len(pathway_in_universe)
    coverage = n_in_universe / pathway_genes_total if pathway_genes_total > 0 else 0.0
    underpowered = n_in_universe < UNDERPOWERED_THRESHOLD

    if n_in_universe == 0:
        return {
            "pathway_id": pathway_id,
            "pathway_name": pathway_name,
            "direction": direction,
            "a": 0,
            "b": len(de_genes),
            "c": 0,
            "d": len(universe) - len(de_genes),
            "odds_ratio": math.nan,
            "p_value": math.nan,
            "fold_enrichment": math.nan,
            "expected": math.nan,
            "padj": math.nan,
            "test_type": test_type,
            "pathway_coverage": 0.0,
            "n_pathway_genes_in_universe": 0,
            "n_pathway_genes_total": pathway_genes_total,
            "underpowered": underpowered,
        }

    a, b, c, d = _compute_2x2(pathway_in_universe, de_genes, universe)
    odds_ratio, p_value = _fisher_row(a, b, c, d)
    fe = _fold_enrichment(a, b, c, d)
    expected = _expected_count(a, b, c, d)

    return {
        "pathway_id": pathway_id,
        "pathway_name": pathway_name,
        "direction": direction,
        "a": a,
        "b": b,
        "c": c,
        "d": d,
        "odds_ratio": odds_ratio,
        "p_value": p_value,
        "fold_enrichment": fe,
        "expected": expected,
        "padj": math.nan,  # filled in later by BH correction
        "test_type": test_type,
        "pathway_coverage": coverage,
        "n_pathway_genes_in_universe": n_in_universe,
        "n_pathway_genes_total": pathway_genes_total,
        "underpowered": underpowered,
    }


def _apply_bh_fdr(df: pd.DataFrame) -> pd.DataFrame:
    """Apply BH FDR correction within each direction group. Modifies padj in-place."""
    result = df.copy()
    for direction in DIRECTIONS:
        mask = result["direction"] == direction
        p_vals = result.loc[mask, "p_value"].values.astype(float)

        # Only correct non-NaN values
        valid_mask = ~np.isnan(p_vals)
        if valid_mask.sum() == 0:
            continue

        padj_vals = np.full_like(p_vals, fill_value=math.nan)
        _, corrected, _, _ = multipletests(p_vals[valid_mask], method="fdr_bh")
        padj_vals[valid_mask] = corrected
        result.loc[mask, "padj"] = padj_vals

    return result


def run_enrichment(
    de_df: pd.DataFrame,
    pathway_defs: pd.DataFrame,
    gene_universe: set,
    table_scope: str,
) -> pd.DataFrame:
    """Run Fisher's exact + BH-FDR enrichment for all pathways × directions.

    Parameters
    ----------
    de_df:
        DataFrame with columns: locus_tag, expression_status
        (values: significant_up, significant_down, not_significant).
        May include a 'timepoint' column (ignored here; handled by
        run_enrichment_all_timepoints).
    pathway_defs:
        DataFrame with columns: pathway_id, pathway_name, locus_tags (set), gene_count (int).
    gene_universe:
        Set of locus_tag strings defining the background.
    table_scope:
        One of: all_detected_genes, filtered_subset,
                significant_any_timepoint, significant_only.

    Returns
    -------
    DataFrame with one row per pathway × direction, including padj (BH-corrected
    within each direction group).
    """
    if table_scope not in TABLE_SCOPE_TO_TEST_TYPE:
        raise ValueError(
            f"Unknown table_scope: {table_scope!r}. "
            f"Valid values: {list(TABLE_SCOPE_TO_TEST_TYPE)}"
        )
    test_type = TABLE_SCOPE_TO_TEST_TYPE[table_scope]

    # Build direction → DE gene sets (intersected with universe)
    up_genes = set(de_df.loc[de_df["expression_status"] == "significant_up", "locus_tag"]) & gene_universe
    down_genes = set(de_df.loc[de_df["expression_status"] == "significant_down", "locus_tag"]) & gene_universe
    combined_genes = up_genes | down_genes

    direction_sets = {
        "up": up_genes,
        "down": down_genes,
        "combined": combined_genes,
    }

    rows = []
    for _, pw_row in pathway_defs.iterrows():
        pathway_id = pw_row["pathway_id"]
        pathway_name = pw_row["pathway_name"]
        pathway_genes_total = int(pw_row["gene_count"])
        pathway_in_universe = pw_row["locus_tags"] & gene_universe

        for direction, de_genes in direction_sets.items():
            row = _enrich_one_pathway_direction(
                pathway_id=pathway_id,
                pathway_name=pathway_name,
                pathway_genes_total=pathway_genes_total,
                pathway_in_universe=pathway_in_universe,
                de_genes=de_genes,
                universe=gene_universe,
                direction=direction,
                test_type=test_type,
            )
            rows.append(row)

    if not rows:
        return pd.DataFrame(columns=[
            "pathway_id", "pathway_name", "direction", "a", "b", "c", "d",
            "odds_ratio", "p_value", "fold_enrichment", "expected", "padj",
            "test_type", "pathway_coverage", "n_pathway_genes_in_universe",
            "n_pathway_genes_total", "underpowered",
        ])

    result_df = pd.DataFrame(rows)
    result_df = _apply_bh_fdr(result_df)
    return result_df.reset_index(drop=True)


def run_enrichment_all_timepoints(
    de_df: pd.DataFrame,
    pathway_defs: pd.DataFrame,
    table_scope: str,
) -> pd.DataFrame:
    """Run enrichment per timepoint and concatenate results.

    Requires de_df to have a 'timepoint' column.
    Gene universe for each timepoint = all locus_tags observed at that timepoint.

    Returns
    -------
    DataFrame with all columns from run_enrichment plus 'timepoint'.
    """
    if "timepoint" not in de_df.columns:
        raise ValueError("de_df must have a 'timepoint' column for run_enrichment_all_timepoints")

    timepoint_results = []
    for tp, tp_df in de_df.groupby("timepoint"):
        universe = set(tp_df["locus_tag"].dropna())
        result = run_enrichment(tp_df, pathway_defs, universe, table_scope)
        result["timepoint"] = tp
        timepoint_results.append(result)

    if not timepoint_results:
        return pd.DataFrame()

    return pd.concat(timepoint_results, ignore_index=True)
