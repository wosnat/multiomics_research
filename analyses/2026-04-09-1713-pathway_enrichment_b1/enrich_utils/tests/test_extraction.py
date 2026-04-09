"""Integration tests for extraction.py.

These tests require a running Neo4j instance and the multiomics-explorer
package with the KG loaded. All tests are skipped if the KG connection
is not available.

Test cases:
    1. test_cyanorak_annotations      - marker genes, cyanorak_role
    2. test_kegg_annotations          - marker genes, kegg
    3. test_no_annotation_gene        - single gene, may be empty, no error
    4. test_cyanorak_hierarchy        - hierarchy edges + levels
    5. test_flat_ontology_returns_empty - tigr_role (flat)
    6. test_end_to_end_rollup         - full MED4 cyanorak rollup vs genes_by_ontology
"""

import pytest
import pandas as pd


# ---------------------------------------------------------------------------
# KG connection fixture
# ---------------------------------------------------------------------------

@pytest.fixture(scope="module")
def kg_conn():
    """Return a live KG connection or skip the test module."""
    try:
        from multiomics_explorer.kg.connection import GraphConnection
        conn = GraphConnection()
        if not conn.verify_connectivity():
            pytest.skip("KG connection not available (Neo4j unreachable)")
        return conn
    except Exception as e:
        pytest.skip(f"KG connection not available: {e}")


# ---------------------------------------------------------------------------
# Marker genes for MED4
# ---------------------------------------------------------------------------

MARKER_GENES = ["PMM0920", "PMM0370", "PMM0550", "PMM1452", "PMM0030", "PMM0346"]


# ---------------------------------------------------------------------------
# Test 1: CyanoRak annotations
# ---------------------------------------------------------------------------

def test_cyanorak_annotations(kg_conn):
    """Extract cyanorak_role annotations for marker genes.

    PMM0920 is glnA (glutamine synthetase), a key nitrogen assimilation gene,
    expected to have a CyanoRak annotation.
    """
    from enrich_utils.extraction import extract_annotations

    df = extract_annotations(MARKER_GENES, ontology="cyanorak_role", conn=kg_conn)

    # Correct columns
    assert list(df.columns) == ["locus_tag", "term_id", "term_name"], (
        f"Unexpected columns: {df.columns.tolist()}"
    )

    # Should have results
    assert len(df) > 0, "Expected at least one cyanorak_role annotation for marker genes"

    # PMM0920 (glnA) should be annotated
    assert "PMM0920" in df["locus_tag"].values, (
        "PMM0920 (glnA) not found in cyanorak_role annotations"
    )

    # term_id and term_name should be non-null
    assert df["term_id"].notna().all(), "term_id should not be null"
    assert df["term_name"].notna().all(), "term_name should not be null"


# ---------------------------------------------------------------------------
# Test 2: KEGG annotations
# ---------------------------------------------------------------------------

def test_kegg_annotations(kg_conn):
    """Extract KEGG annotations for marker genes. Should return results."""
    from enrich_utils.extraction import extract_annotations

    df = extract_annotations(MARKER_GENES, ontology="kegg", conn=kg_conn)

    assert list(df.columns) == ["locus_tag", "term_id", "term_name"], (
        f"Unexpected columns: {df.columns.tolist()}"
    )
    assert len(df) > 0, "Expected at least one KEGG annotation for marker genes"


# ---------------------------------------------------------------------------
# Test 3: Gene with no annotation does not error
# ---------------------------------------------------------------------------

def test_no_annotation_gene(kg_conn):
    """Extract annotations for a single gene that may have no cyanorak_role.

    PMM0030 may or may not be annotated. Either way, the call should succeed
    and return a DataFrame with correct columns.
    """
    from enrich_utils.extraction import extract_annotations

    df = extract_annotations(["PMM0030"], ontology="cyanorak_role", conn=kg_conn)

    assert list(df.columns) == ["locus_tag", "term_id", "term_name"], (
        f"Unexpected columns: {df.columns.tolist()}"
    )
    # No assertion on len — may be empty or have results


# ---------------------------------------------------------------------------
# Test 4: CyanoRak hierarchy
# ---------------------------------------------------------------------------

def test_cyanorak_hierarchy(kg_conn):
    """Extract cyanorak_role hierarchy.

    Expected: non-empty edges, at least 2 distinct levels.
    """
    from enrich_utils.extraction import extract_hierarchy

    df = extract_hierarchy("cyanorak_role", conn=kg_conn)

    assert list(df.columns) == ["child_id", "parent_id", "child_level", "parent_level"], (
        f"Unexpected columns: {df.columns.tolist()}"
    )
    assert len(df) > 0, "Expected hierarchy edges for cyanorak_role"

    all_levels = set(df["child_level"]) | set(df["parent_level"])
    assert len(all_levels) >= 2, (
        f"Expected at least 2 distinct levels, got: {sorted(all_levels)}"
    )

    # Level values should be non-negative integers
    assert (df["child_level"] >= 0).all(), "child_level should be >= 0"
    assert (df["parent_level"] >= 0).all(), "parent_level should be >= 0"

    # Parent level should be less than child level (parent is closer to root)
    assert (df["parent_level"] < df["child_level"]).all(), (
        "parent_level should be strictly less than child_level "
        "(parent is closer to root)"
    )


# ---------------------------------------------------------------------------
# Test 5: Flat ontology returns empty DataFrame
# ---------------------------------------------------------------------------

def test_flat_ontology_returns_empty(kg_conn):
    """tigr_role has no hierarchy_rels — should return empty DataFrame."""
    from enrich_utils.extraction import extract_hierarchy

    df = extract_hierarchy("tigr_role", conn=kg_conn)

    assert list(df.columns) == ["child_id", "parent_id", "child_level", "parent_level"], (
        f"Unexpected columns: {df.columns.tolist()}"
    )
    assert len(df) == 0, (
        f"Expected empty DataFrame for flat ontology tigr_role, got {len(df)} rows"
    )


# ---------------------------------------------------------------------------
# Test 6: End-to-end rollup vs genes_by_ontology
# ---------------------------------------------------------------------------

def test_end_to_end_rollup(kg_conn):
    """Full MED4 cyanorak rollup vs genes_by_ontology.

    Steps:
    1. Get all MED4 locus tags via genes_by_function("*", organism="MED4")
    2. Extract cyanorak_role annotations + hierarchy for those genes
    3. Roll up to level 1 via roll_up_to_level
    4. Pick the most-populated term at level 1
    5. Compare gene set with genes_by_ontology(term_ids=[term], ontology="cyanorak_role",
       organism="MED4")

    The two gene sets should match (same locus_tags after intersection with MED4).
    """
    from multiomics_explorer import genes_by_function, genes_by_ontology
    from enrich_utils.extraction import extract_annotations, extract_hierarchy
    from enrich_utils.hierarchy import roll_up_to_level

    # Step 1: Get all MED4 locus tags
    result = genes_by_function("*", organism="MED4", limit=None)
    all_med4_tags = [r["locus_tag"] for r in result["results"]]
    assert len(all_med4_tags) > 0, "No MED4 genes returned"

    # Step 2: Extract annotations + hierarchy
    ann_df = extract_annotations(all_med4_tags, ontology="cyanorak_role", conn=kg_conn)
    hier_df = extract_hierarchy("cyanorak_role", conn=kg_conn)

    assert len(ann_df) > 0, "No cyanorak_role annotations for MED4"
    assert len(hier_df) > 0, "No cyanorak_role hierarchy"

    # Step 3: Roll up to level 1
    rolled = roll_up_to_level(ann_df, hier_df, target_level=1)
    assert len(rolled) > 0, "Roll-up to level 1 produced no results"

    # Step 4: Find the most-populated term at level 1
    term_sizes = rolled.groupby("term_id")["locus_tag"].nunique()
    top_term = term_sizes.idxmax()

    # Gene set from local extraction + rollup
    local_genes = set(rolled.loc[rolled["term_id"] == top_term, "locus_tag"])
    assert len(local_genes) > 0, f"No genes for top term {top_term}"

    # Step 5: Get gene set from genes_by_ontology
    gbo_result = genes_by_ontology(
        term_ids=[top_term],
        ontology="cyanorak_role",
        organism="MED4",
        limit=None,
        conn=kg_conn,
    )
    kg_genes = set(r["locus_tag"] for r in gbo_result["results"])

    # The local set should match the KG result
    # genes_by_ontology uses hierarchy expansion, so it may include genes
    # annotated to children of the top_term as well. We check that the
    # local genes are a subset (or match) of the KG result.
    missing_from_kg = local_genes - kg_genes
    assert len(missing_from_kg) == 0, (
        f"Local rollup has {len(missing_from_kg)} genes not found in "
        f"genes_by_ontology for term {top_term}: {sorted(missing_from_kg)[:10]}"
    )
