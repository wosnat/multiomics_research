"""Tests for enrichment.py — Fisher's exact + FDR.

All 2×2 tables follow the convention:
    a = pathway genes that are DE (in the test set)
    b = non-pathway genes that are DE
    c = pathway genes that are NOT DE
    d = non-pathway genes that are NOT DE

One-sided Fisher's exact (over-representation): alternative='greater'.
"""

import math
import numpy as np
import pandas as pd
import pytest
from scipy.stats import fisher_exact

from enrich_utils.enrichment import run_enrichment, run_enrichment_all_timepoints


# ---------------------------------------------------------------------------
# Helpers to build minimal input DataFrames
# ---------------------------------------------------------------------------

def make_de_df(up_genes, down_genes, not_sig_genes, timepoint=None):
    """Build a de_df with expression_status and optional timepoint columns."""
    rows = []
    for g in up_genes:
        row = {"locus_tag": g, "expression_status": "significant_up"}
        if timepoint is not None:
            row["timepoint"] = timepoint
        rows.append(row)
    for g in down_genes:
        row = {"locus_tag": g, "expression_status": "significant_down"}
        if timepoint is not None:
            row["timepoint"] = timepoint
        rows.append(row)
    for g in not_sig_genes:
        row = {"locus_tag": g, "expression_status": "not_significant"}
        if timepoint is not None:
            row["timepoint"] = timepoint
        rows.append(row)
    return pd.DataFrame(rows)


def make_pathway(pathway_id, pathway_name, locus_tags):
    return {
        "pathway_id": pathway_id,
        "pathway_name": pathway_name,
        "locus_tags": set(locus_tags),
        "gene_count": len(locus_tags),
    }


def make_pathway_defs(pathways):
    return pd.DataFrame(pathways)


# ---------------------------------------------------------------------------
# Case 1: Strong up-enrichment
# ---------------------------------------------------------------------------

def test_strong_up_enrichment():
    """
    100-gene universe. 10-gene pathway. 30 DE up total. 8 of 10 pathway genes are DE up.

    2×2 table for 'up' direction:
        a = 8  (pathway ∩ DE-up)
        b = 22 (non-pathway ∩ DE-up)  = 30 - 8
        c = 2  (pathway ∩ not DE-up)  = 10 - 8
        d = 68 (non-pathway ∩ not DE-up) = 100 - 30 - (10 - 8) = 68

    Expected values computed with scipy.stats.fisher_exact([[8,22],[2,68]], alternative='greater').
    """
    pathway_genes = [f"gene{i}" for i in range(10)]
    up_genes = pathway_genes[:8] + [f"up{i}" for i in range(22)]  # 8 pathway + 22 non-pathway = 30 up
    not_sig = pathway_genes[8:] + [f"ns{i}" for i in range(68)]   # 2 pathway + 68 non-pathway = 70 not sig
    down_genes = []

    universe = set(pathway_genes + [f"up{i}" for i in range(22)] + [f"ns{i}" for i in range(68)])
    assert len(universe) == 100

    de_df = make_de_df(up_genes, down_genes, not_sig)
    pathway_defs = make_pathway_defs([make_pathway("PW001", "Test Pathway", pathway_genes)])

    result = run_enrichment(de_df, pathway_defs, universe, "all_detected_genes")

    up_row = result[(result["pathway_id"] == "PW001") & (result["direction"] == "up")]
    assert len(up_row) == 1, "Expected one row for PW001 / up"
    row = up_row.iloc[0]

    assert row["a"] == 8
    assert row["b"] == 22
    assert row["c"] == 2
    assert row["d"] == 68

    expected_or, expected_p = fisher_exact([[8, 22], [2, 68]], alternative="greater")
    assert math.isclose(row["odds_ratio"], expected_or, rel_tol=1e-6)
    assert math.isclose(row["p_value"], expected_p, rel_tol=1e-6)
    assert row["test_type"] == "vs_genome"
    assert row["p_value"] < 0.05


# ---------------------------------------------------------------------------
# Case 2: No enrichment
# ---------------------------------------------------------------------------

def test_no_enrichment():
    """
    100-gene universe. 10-gene pathway. 30 DE up. Only 1 of 10 pathway genes DE up.

    a=1, b=29, c=9, d=61. Should be not significant.
    """
    pathway_genes = [f"gene{i}" for i in range(10)]
    up_genes = pathway_genes[:1] + [f"up{i}" for i in range(29)]
    not_sig = pathway_genes[1:] + [f"ns{i}" for i in range(61)]

    universe = set(pathway_genes) | set(up_genes) | set(not_sig)
    assert len(universe) == 100

    de_df = make_de_df(up_genes, [], not_sig)
    pathway_defs = make_pathway_defs([make_pathway("PW002", "No Enrich Pathway", pathway_genes)])

    result = run_enrichment(de_df, pathway_defs, universe, "all_detected_genes")
    up_row = result[(result["pathway_id"] == "PW002") & (result["direction"] == "up")]
    assert len(up_row) == 1
    assert up_row.iloc[0]["p_value"] > 0.05


# ---------------------------------------------------------------------------
# Case 3: Down enrichment (mirror of case 1)
# ---------------------------------------------------------------------------

def test_strong_down_enrichment():
    """Same as case 1 but direction = down."""
    pathway_genes = [f"gene{i}" for i in range(10)]
    down_genes = pathway_genes[:8] + [f"dn{i}" for i in range(22)]
    not_sig = pathway_genes[8:] + [f"ns{i}" for i in range(68)]

    universe = set(pathway_genes) | set(down_genes) | set(not_sig)
    assert len(universe) == 100

    de_df = make_de_df([], down_genes, not_sig)
    pathway_defs = make_pathway_defs([make_pathway("PW003", "Down Pathway", pathway_genes)])

    result = run_enrichment(de_df, pathway_defs, universe, "all_detected_genes")
    down_row = result[(result["pathway_id"] == "PW003") & (result["direction"] == "down")]
    assert len(down_row) == 1
    row = down_row.iloc[0]

    assert row["a"] == 8
    assert row["b"] == 22
    assert row["c"] == 2
    assert row["d"] == 68

    expected_or, expected_p = fisher_exact([[8, 22], [2, 68]], alternative="greater")
    assert math.isclose(row["odds_ratio"], expected_or, rel_tol=1e-6)
    assert math.isclose(row["p_value"], expected_p, rel_tol=1e-6)
    assert row["p_value"] < 0.05


# ---------------------------------------------------------------------------
# Case 4: Combined test counts both up and down
# ---------------------------------------------------------------------------

def test_combined_direction():
    """
    10-gene pathway: 5 are DE up, 5 are DE down.
    Combined test: a = 10 (all pathway genes are DE in either direction).
    Separate up: a = 5.
    """
    pathway_genes = [f"gene{i}" for i in range(10)]
    up_genes = pathway_genes[:5] + [f"up{i}" for i in range(25)]   # 5 pathway + 25 non-pathway = 30 up
    down_genes = pathway_genes[5:] + [f"dn{i}" for i in range(25)] # 5 pathway + 25 non-pathway = 30 down
    not_sig = [f"ns{i}" for i in range(40)]

    universe = set(pathway_genes) | set(up_genes) | set(down_genes) | set(not_sig)
    assert len(universe) == 100

    de_df = make_de_df(up_genes, down_genes, not_sig)
    pathway_defs = make_pathway_defs([make_pathway("PW004", "Combined Pathway", pathway_genes)])

    result = run_enrichment(de_df, pathway_defs, universe, "all_detected_genes")

    combined_row = result[(result["pathway_id"] == "PW004") & (result["direction"] == "combined")]
    up_row = result[(result["pathway_id"] == "PW004") & (result["direction"] == "up")]

    assert len(combined_row) == 1
    assert combined_row.iloc[0]["a"] == 10

    assert len(up_row) == 1
    assert up_row.iloc[0]["a"] == 5


# ---------------------------------------------------------------------------
# Case 5: Empty pathway in universe (pathway genes disjoint from universe)
# ---------------------------------------------------------------------------

def test_empty_pathway_in_universe():
    """
    Pathway genes are NOT in the gene universe.
    pathway_coverage should be 0. p_value should be NaN (or skipped).
    """
    pathway_genes = ["phantom1", "phantom2", "phantom3"]
    universe_genes = [f"gene{i}" for i in range(50)]
    up_genes = [f"gene{i}" for i in range(10)]
    not_sig = [f"gene{i}" for i in range(10, 50)]

    de_df = make_de_df(up_genes, [], not_sig)
    pathway_defs = make_pathway_defs([make_pathway("PW005", "Ghost Pathway", pathway_genes)])
    universe = set(universe_genes)

    result = run_enrichment(de_df, pathway_defs, universe, "all_detected_genes")

    pw_rows = result[result["pathway_id"] == "PW005"]
    assert len(pw_rows) > 0, "Should still have rows for the pathway"
    for _, row in pw_rows.iterrows():
        assert row["pathway_coverage"] == 0.0
        assert math.isnan(row["p_value"]) or row["p_value"] is None or np.isnan(row["p_value"])


# ---------------------------------------------------------------------------
# Case 6: Underpowered (< 5 pathway genes in universe)
# ---------------------------------------------------------------------------

def test_underpowered_flag():
    """
    10-gene pathway but only 3 genes are in the universe.
    coverage = 0.3, underpowered = True.
    """
    pathway_genes = [f"gene{i}" for i in range(10)]
    in_universe = pathway_genes[:3]  # only 3 in universe
    other_universe = [f"other{i}" for i in range(97)]
    universe = set(in_universe) | set(other_universe)  # 100 genes, 3 from pathway

    up_genes = in_universe + [f"other{i}" for i in range(27)]  # 30 up
    not_sig = [f"other{i}" for i in range(27, 70)]

    de_df = make_de_df(up_genes, [], not_sig)
    pathway_defs = make_pathway_defs([make_pathway("PW006", "Underpowered", pathway_genes)])

    result = run_enrichment(de_df, pathway_defs, universe, "all_detected_genes")

    pw_rows = result[result["pathway_id"] == "PW006"]
    assert len(pw_rows) > 0
    for _, row in pw_rows.iterrows():
        assert math.isclose(row["pathway_coverage"], 3 / 10, rel_tol=1e-6)
        assert row["underpowered"] is True or row["underpowered"] == True


# ---------------------------------------------------------------------------
# Case 7: All pathway genes are DE — very significant
# ---------------------------------------------------------------------------

def test_all_pathway_genes_de():
    """
    10-gene pathway, all 10 are DE up. 30 total DE up out of 100 genes.
    a=10, c=0. Should have a very significant p-value.
    """
    pathway_genes = [f"gene{i}" for i in range(10)]
    up_genes = pathway_genes + [f"up{i}" for i in range(20)]  # 10 pathway + 20 non-pathway = 30 up
    not_sig = [f"ns{i}" for i in range(70)]

    universe = set(pathway_genes) | set(up_genes) | set(not_sig)
    assert len(universe) == 100

    de_df = make_de_df(up_genes, [], not_sig)
    pathway_defs = make_pathway_defs([make_pathway("PW007", "Full Hit Pathway", pathway_genes)])

    result = run_enrichment(de_df, pathway_defs, universe, "all_detected_genes")
    up_row = result[(result["pathway_id"] == "PW007") & (result["direction"] == "up")]
    assert len(up_row) == 1
    row = up_row.iloc[0]

    assert row["a"] == 10
    assert row["c"] == 0

    expected_or, expected_p = fisher_exact([[10, 20], [0, 70]], alternative="greater")
    assert math.isclose(row["p_value"], expected_p, rel_tol=1e-6)
    assert row["p_value"] < 0.001


# ---------------------------------------------------------------------------
# Case 8: FDR correction across multiple pathways
# ---------------------------------------------------------------------------

def test_fdr_correction():
    """
    20 pathways: 2 strongly enriched, 18 not enriched.
    padj >= p_value for all rows (BH correction never makes p smaller).
    Enriched pathways have padj < 0.05.
    For direction='up' group only (to keep the test focused).
    """
    n_total = 200
    universe = set(f"gene{i}" for i in range(n_total))

    # Build 20 pathways: first 2 strongly enriched (8/10 DE up), rest have 1/10
    pathways = []
    up_genes_set = set()
    not_sig_set = set()

    all_pathway_genes = []
    for pw_idx in range(20):
        pw_genes = [f"gene{pw_idx * 10 + j}" for j in range(10)]
        all_pathway_genes.extend(pw_genes)
        if pw_idx < 2:
            up_genes_set.update(pw_genes[:8])
            not_sig_set.update(pw_genes[8:])
        else:
            up_genes_set.add(pw_genes[0])
            not_sig_set.update(pw_genes[1:])
        pathways.append(make_pathway(f"PW{pw_idx:03d}", f"Pathway {pw_idx}", pw_genes))

    # Fill remaining universe with not-significant non-pathway genes
    used = set(all_pathway_genes)
    remaining = [g for g in universe if g not in used]
    not_sig_set.update(remaining)

    pathway_defs = make_pathway_defs(pathways)
    de_df = make_de_df(list(up_genes_set), [], list(not_sig_set))

    result = run_enrichment(de_df, pathway_defs, universe, "all_detected_genes")

    up_results = result[result["direction"] == "up"].copy()

    # padj >= p_value for all rows
    valid = up_results[~up_results["p_value"].isna()]
    assert (valid["padj"] >= valid["p_value"] - 1e-10).all(), \
        "padj should be >= p_value (BH never reduces p)"

    # Enriched pathways have padj < 0.05
    enriched = up_results[up_results["pathway_id"].isin(["PW000", "PW001"])]
    assert (enriched["padj"] < 0.05).all(), "Strongly enriched pathways should have padj < 0.05"

    # Most non-enriched pathways have padj > 0.05
    non_enriched = up_results[~up_results["pathway_id"].isin(["PW000", "PW001"])]
    assert (non_enriched["padj"] > 0.05).mean() > 0.8, \
        "Most non-enriched pathways should have padj > 0.05"


# ---------------------------------------------------------------------------
# Case 9: test_type mapping for all 4 table_scope values
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("table_scope,expected_test_type", [
    ("all_detected_genes", "vs_genome"),
    ("filtered_subset", "vs_filtered_genome"),
    ("significant_any_timepoint", "vs_all_responsive"),
    ("significant_only", "descriptive_only"),
])
def test_test_type_mapping(table_scope, expected_test_type):
    """Each table_scope maps to a specific test_type string in all result rows."""
    pathway_genes = [f"gene{i}" for i in range(5)]
    up_genes = pathway_genes[:2] + [f"up{i}" for i in range(10)]
    not_sig = pathway_genes[2:] + [f"ns{i}" for i in range(30)]
    universe = set(pathway_genes) | set(up_genes) | set(not_sig)

    de_df = make_de_df(up_genes, [], not_sig)
    pathway_defs = make_pathway_defs([make_pathway("PW_SCOPE", "Scope Test", pathway_genes)])

    result = run_enrichment(de_df, pathway_defs, universe, table_scope)
    assert len(result) > 0
    assert (result["test_type"] == expected_test_type).all(), \
        f"Expected test_type='{expected_test_type}' for table_scope='{table_scope}'"


# ---------------------------------------------------------------------------
# Case: run_enrichment_all_timepoints concatenates results per timepoint
# ---------------------------------------------------------------------------

def test_run_enrichment_all_timepoints():
    """
    Two timepoints (T1, T2). Should produce results for each timepoint.
    Output has a 'timepoint' column and both timepoints are represented.
    """
    pathway_genes = [f"gene{i}" for i in range(10)]
    up_genes_t1 = pathway_genes[:5] + [f"up_t1_{i}" for i in range(25)]
    not_sig_t1 = pathway_genes[5:] + [f"ns_t1_{i}" for i in range(65)]

    up_genes_t2 = pathway_genes[:8] + [f"up_t2_{i}" for i in range(22)]
    not_sig_t2 = pathway_genes[8:] + [f"ns_t2_{i}" for i in range(68)]

    df_t1 = make_de_df(up_genes_t1, [], not_sig_t1, timepoint="T1")
    df_t2 = make_de_df(up_genes_t2, [], not_sig_t2, timepoint="T2")
    de_df = pd.concat([df_t1, df_t2], ignore_index=True)

    pathway_defs = make_pathway_defs([make_pathway("PW_TP", "Timepoint Pathway", pathway_genes)])

    # Universe is all genes across both timepoints
    universe = (set(pathway_genes) | set(up_genes_t1) | set(not_sig_t1) |
                set(up_genes_t2) | set(not_sig_t2))

    result = run_enrichment_all_timepoints(de_df, pathway_defs, "all_detected_genes")

    assert "timepoint" in result.columns, "Output must have 'timepoint' column"
    assert set(result["timepoint"].unique()) == {"T1", "T2"}
    assert "pathway_id" in result.columns


# ---------------------------------------------------------------------------
# Signed enrichment score tests
# ---------------------------------------------------------------------------

def test_signed_score_up_only():
    """Pathway significant up only → positive score."""
    from enrich_utils.enrichment import signed_enrichment_score
    df = pd.DataFrame([
        {"pathway_id": "P1", "pathway_name": "Path1", "direction": "up", "padj": 0.001},
        {"pathway_id": "P1", "pathway_name": "Path1", "direction": "down", "padj": 0.5},
        {"pathway_id": "P1", "pathway_name": "Path1", "direction": "combined", "padj": 0.01},
    ])
    result = signed_enrichment_score(df)
    assert len(result) == 1
    assert result.iloc[0]["score"] == pytest.approx(-np.log10(0.001))
    assert result.iloc[0]["dominant_direction"] == "up"


def test_signed_score_down_only():
    """Pathway significant down only → negative score."""
    from enrich_utils.enrichment import signed_enrichment_score
    df = pd.DataFrame([
        {"pathway_id": "P1", "pathway_name": "Path1", "direction": "up", "padj": 0.8},
        {"pathway_id": "P1", "pathway_name": "Path1", "direction": "down", "padj": 0.001},
        {"pathway_id": "P1", "pathway_name": "Path1", "direction": "combined", "padj": 0.01},
    ])
    result = signed_enrichment_score(df)
    assert result.iloc[0]["score"] == pytest.approx(np.log10(0.001))  # negative
    assert result.iloc[0]["dominant_direction"] == "down"


def test_signed_score_both_sig_up_wins():
    """Both significant, up has lower padj → positive score from up."""
    from enrich_utils.enrichment import signed_enrichment_score
    df = pd.DataFrame([
        {"pathway_id": "P1", "pathway_name": "Path1", "direction": "up", "padj": 0.001},
        {"pathway_id": "P1", "pathway_name": "Path1", "direction": "down", "padj": 0.01},
        {"pathway_id": "P1", "pathway_name": "Path1", "direction": "combined", "padj": 0.001},
    ])
    result = signed_enrichment_score(df)
    assert result.iloc[0]["score"] > 0
    assert result.iloc[0]["dominant_direction"] == "up"


def test_signed_score_neither_sig():
    """Neither direction significant → score = 0."""
    from enrich_utils.enrichment import signed_enrichment_score
    df = pd.DataFrame([
        {"pathway_id": "P1", "pathway_name": "Path1", "direction": "up", "padj": 0.5},
        {"pathway_id": "P1", "pathway_name": "Path1", "direction": "down", "padj": 0.8},
        {"pathway_id": "P1", "pathway_name": "Path1", "direction": "combined", "padj": 0.3},
    ])
    result = signed_enrichment_score(df)
    assert result.iloc[0]["score"] == 0.0
    assert result.iloc[0]["dominant_direction"] == "ns"
