"""Reusable DE data extraction from the multiomics KG.

Wraps differential_expression_by_gene() to produce DataFrames in a
consistent schema with table_scope on every row. All extraction scripts
use this module.
"""

import pandas as pd
from multiomics_explorer import differential_expression_by_gene
from multiomics_explorer.analysis import to_dataframe


def extract_de(
    experiment_ids: list[str],
    organism: str = "MED4",
    verbose: bool = True,
    significant_only: bool = False,
) -> tuple[pd.DataFrame, dict]:
    """Extract DE data for experiments, returning DataFrame + envelope metadata.

    Calls differential_expression_by_gene with verbose=True, limit=None.
    Captures table_scope from the envelope and ensures it's on every row
    (verbose=True already includes it per-row, but we verify).

    Args:
        experiment_ids: List of experiment IDs to extract.
        organism: Organism name (fuzzy matched).
        verbose: Include product, gene_category, etc. Default True.
        significant_only: If True, only return significant rows.

    Returns:
        (df, envelope) where:
        - df: DataFrame with all DE rows, table_scope as a column
        - envelope: dict with metadata (total_matching, experiments, etc.)
    """
    result = differential_expression_by_gene(
        experiment_ids=experiment_ids,
        organism=organism,
        verbose=verbose,
        significant_only=significant_only,
        limit=None,
    )

    assert not result["truncated"], (
        f"Results truncated! total_matching={result['total_matching']}, "
        f"returned={result['returned']}"
    )

    df = to_dataframe(result)

    # Verify table_scope is present (verbose=True puts it on rows)
    if "table_scope" not in df.columns and result.get("experiments"):
        # Fallback: get from envelope and broadcast
        scope_map = {
            exp["experiment_id"]: exp["table_scope"]
            for exp in result["experiments"]
        }
        df["table_scope"] = df["experiment_id"].map(scope_map)

    # Ensure rank columns are nullable int
    for col in ["rank", "rank_up", "rank_down"]:
        if col in df.columns:
            df[col] = pd.to_numeric(df[col], errors="coerce").astype("Int64")

    # Build envelope summary (without the full results)
    envelope = {
        k: v for k, v in result.items()
        if k != "results"
    }

    return df, envelope


def extraction_summary(df: pd.DataFrame, envelope: dict) -> str:
    """Produce a human-readable summary of an extraction for logging.

    Args:
        df: DataFrame from extract_de.
        envelope: Envelope dict from extract_de.

    Returns:
        Multi-line summary string.
    """
    lines = [
        f"Rows: {len(df)}",
        f"Genes: {df['locus_tag'].nunique()}",
        f"Experiments: {envelope.get('experiment_count', '?')}",
    ]

    if "timepoint" in df.columns:
        lines.append(f"Timepoints: {sorted(df['timepoint'].unique())}")

    if "expression_status" in df.columns:
        lines.append(f"Status: {df['expression_status'].value_counts().to_dict()}")

    if "table_scope" in df.columns:
        lines.append(f"Table scope: {df['table_scope'].unique().tolist()}")

    if "omics_type" in df.columns:
        lines.append(f"Omics type: {df['omics_type'].unique().tolist()}")

    return "\n".join(lines)


def check_marker_genes(
    df: pd.DataFrame,
    marker_locus_tags: list[str],
) -> pd.DataFrame:
    """Extract rows for marker genes from a DE DataFrame.

    Returns a filtered DataFrame with only the marker gene rows,
    useful for QC and exploration logging.

    Args:
        df: DE DataFrame.
        marker_locus_tags: List of locus tags to check.

    Returns:
        DataFrame filtered to marker genes, sorted by locus_tag and timepoint.
    """
    marker_df = df[df["locus_tag"].isin(marker_locus_tags)].copy()

    sort_cols = ["locus_tag"]
    if "timepoint_order" in marker_df.columns:
        sort_cols.append("timepoint_order")
    elif "timepoint" in marker_df.columns:
        sort_cols.append("timepoint")

    return marker_df.sort_values(sort_cols).reset_index(drop=True)
