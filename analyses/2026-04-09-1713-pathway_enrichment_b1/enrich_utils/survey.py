"""Annotation landscape profiling and pathway definition utilities.

Public API:
    survey_ontology(annotations_df, hierarchy_df, gene_universe) -> dict
    rank_ontologies(profiles) -> DataFrame
    build_pathway_definitions(annotations_df, hierarchy_df, level,
                              min_genes, term_names=None) -> DataFrame
    scope_pathways_to_universe(pathway_defs, gene_universe) -> DataFrame
"""

import numpy as np
import pandas as pd

from enrich_utils.hierarchy import roll_up_to_level


# ---------------------------------------------------------------------------
# survey_ontology
# ---------------------------------------------------------------------------

def survey_ontology(
    annotations_df: pd.DataFrame,
    hierarchy_df: pd.DataFrame,
    gene_universe: set,
) -> dict:
    """Profile annotation coverage and term-size distribution per level.

    Parameters
    ----------
    annotations_df:
        DataFrame with columns: locus_tag, term_id.
    hierarchy_df:
        DataFrame with columns: child_id, parent_id, child_level, parent_level.
    gene_universe:
        Set of locus_tag strings defining the study background.

    Returns
    -------
    dict with:
        coverage         - fraction of gene_universe with >= 1 annotation
        n_annotated      - count of annotated genes in universe
        n_unannotated    - count of unannotated genes in universe
        per_level        - list of dicts, one per level present in hierarchy_df:
            level                - integer level
            n_genes_at_level     - distinct genes reachable at this level
            gene_coverage        - n_genes_at_level / n_annotated (fraction of
                                   annotated genes that appear at this level)
            n_terms_with_genes   - number of terms that have >= 1 annotated gene
            min_genes            - minimum term gene count (across terms with >= 1)
            q1_genes             - 25th percentile
            median_genes         - 50th percentile
            q3_genes             - 75th percentile
            max_genes            - maximum term gene count
    """
    # Annotated genes restricted to universe
    annotated_in_universe = set(annotations_df["locus_tag"]) & gene_universe
    n_annotated = len(annotated_in_universe)
    n_unannotated = len(gene_universe) - n_annotated
    coverage = n_annotated / len(gene_universe) if gene_universe else 0.0

    # Determine all levels present in hierarchy
    levels = sorted(set(hierarchy_df["child_level"]) | set(hierarchy_df["parent_level"]))

    per_level = []
    for level in levels:
        rolled = roll_up_to_level(annotations_df, hierarchy_df, target_level=level)
        # Restrict to universe genes
        rolled_in_universe = rolled[rolled["locus_tag"].isin(gene_universe)]
        n_genes_at_level = rolled_in_universe["locus_tag"].nunique()
        gene_cov = n_genes_at_level / n_annotated if n_annotated > 0 else 0.0

        if len(rolled_in_universe) == 0:
            per_level.append({
                "level": level,
                "n_genes_at_level": 0,
                "gene_coverage": 0.0,
                "n_terms_with_genes": 0,
                "min_genes": None,
                "q1_genes": None,
                "median_genes": None,
                "q3_genes": None,
                "max_genes": None,
            })
            continue

        term_sizes = rolled_in_universe.groupby("term_id")["locus_tag"].nunique()
        sizes = term_sizes.values

        per_level.append({
            "level": level,
            "n_genes_at_level": n_genes_at_level,
            "gene_coverage": gene_cov,
            "n_terms_with_genes": len(sizes),
            "min_genes": int(sizes.min()),
            "q1_genes": float(np.percentile(sizes, 25)),
            "median_genes": float(np.median(sizes)),
            "q3_genes": float(np.percentile(sizes, 75)),
            "max_genes": int(sizes.max()),
        })

    return {
        "coverage": coverage,
        "n_annotated": n_annotated,
        "n_unannotated": n_unannotated,
        "per_level": per_level,
    }


# ---------------------------------------------------------------------------
# rank_ontologies
# ---------------------------------------------------------------------------

def rank_ontologies(profiles: dict) -> pd.DataFrame:
    """Rank ontologies by coverage + sweet-spot term size.

    Parameters
    ----------
    profiles:
        dict mapping ontology name -> profile dict (as returned by survey_ontology).
        At minimum each profile must have:
            coverage  - float 0-1
            per_level - list of per-level dicts with median_genes

    Returns
    -------
    DataFrame with columns: ontology, coverage, sweet_spot_score, rank.
    Sorted by rank ascending (best first).

    Scoring per level:
        A level qualifies if: median term size 5-50, gene coverage >= 0.5
        Best level = qualifying level with highest gene_coverage * coverage
        score = coverage * best_level_gene_coverage  (0 if no qualifying level)
    """
    rows = []
    for ontology, profile in profiles.items():
        cov = profile.get("coverage", 0.0)

        # Find the best qualifying level
        best_level = None
        best_level_gene_cov = 0.0
        best_level_median = 0.0
        best_level_n_terms = 0

        for lvl in profile.get("per_level", []):
            median = lvl.get("median_genes")
            gene_cov = lvl.get("gene_coverage", 0.0)
            if median is not None and 5 <= median <= 50 and gene_cov >= 0.5:
                if gene_cov > best_level_gene_cov:
                    best_level = lvl["level"]
                    best_level_gene_cov = gene_cov
                    best_level_median = median
                    best_level_n_terms = lvl.get("n_terms_with_genes", 0)

        # Flat ontologies with no hierarchy: treat leaf as only level
        if not profile.get("per_level") and cov > 0:
            best_level_gene_cov = 0.0  # no qualifying level

        score = cov * best_level_gene_cov

        rows.append({
            "ontology": ontology,
            "coverage": cov,
            "best_level": best_level,
            "best_level_gene_coverage": best_level_gene_cov,
            "best_level_median_genes": best_level_median,
            "best_level_n_terms": best_level_n_terms,
            "score": score,
        })

    if not rows:
        return pd.DataFrame(columns=["ontology", "coverage", "median_term_size",
                                     "sweet_spot_score", "score", "rank"])

    df = pd.DataFrame(rows)
    df = df.sort_values("score", ascending=False).reset_index(drop=True)
    df["rank"] = df.index + 1
    return df


# ---------------------------------------------------------------------------
# build_pathway_definitions
# ---------------------------------------------------------------------------

def build_pathway_definitions(
    annotations_df: pd.DataFrame,
    hierarchy_df: pd.DataFrame,
    level: int,
    min_genes: int,
    term_names: dict | None = None,
) -> pd.DataFrame:
    """Build pathway definitions at a given ontology level.

    Rolls gene annotations up to `level` via roll_up_to_level, then groups
    by term to produce one pathway per term.

    Parameters
    ----------
    annotations_df:
        DataFrame with columns: locus_tag, term_id.
    hierarchy_df:
        DataFrame with columns: child_id, parent_id, child_level, parent_level.
    level:
        Target ontology level.
    min_genes:
        Minimum number of genes a pathway must have to be included.
    term_names:
        Optional dict mapping term_id -> human-readable name.
        If None, pathway_name defaults to pathway_id.

    Returns
    -------
    DataFrame with columns:
        pathway_id    - term identifier at target level
        pathway_name  - human-readable name (from term_names or = pathway_id)
        locus_tags    - set of locus_tag strings
        gene_count    - len(locus_tags)
    Filtered to rows where gene_count >= min_genes.
    """
    rolled = roll_up_to_level(annotations_df, hierarchy_df, target_level=level)

    if len(rolled) == 0:
        return pd.DataFrame(columns=["pathway_id", "pathway_name", "locus_tags", "gene_count"])

    grouped = rolled.groupby("term_id")["locus_tag"].apply(set).reset_index()
    grouped.columns = ["pathway_id", "locus_tags"]
    grouped["gene_count"] = grouped["locus_tags"].apply(len)

    # Apply min_genes filter
    grouped = grouped[grouped["gene_count"] >= min_genes].copy()

    # Apply term names
    if term_names is not None:
        grouped["pathway_name"] = grouped["pathway_id"].map(
            lambda tid: term_names.get(tid, tid)
        )
    else:
        grouped["pathway_name"] = grouped["pathway_id"]

    # Reorder columns
    result = grouped[["pathway_id", "pathway_name", "locus_tags", "gene_count"]]
    return result.reset_index(drop=True)


# ---------------------------------------------------------------------------
# scope_pathways_to_universe
# ---------------------------------------------------------------------------

def scope_pathways_to_universe(
    pathway_defs: pd.DataFrame,
    gene_universe: set,
) -> pd.DataFrame:
    """Add universe-scoped columns to a pathway definitions DataFrame.

    Parameters
    ----------
    pathway_defs:
        DataFrame as returned by build_pathway_definitions.
        Must have columns: pathway_id, locus_tags, gene_count.
    gene_universe:
        Set of locus_tag strings defining the study background.

    Returns
    -------
    Copy of pathway_defs with three additional columns:
        n_in_universe         - count of pathway genes in universe
        coverage              - n_in_universe / gene_count
        locus_tags_in_universe - set of pathway genes in universe
    """
    result = pathway_defs.copy()

    result["locus_tags_in_universe"] = result["locus_tags"].apply(
        lambda tags: tags & gene_universe
    )
    result["n_in_universe"] = result["locus_tags_in_universe"].apply(len)
    result["coverage"] = result.apply(
        lambda row: row["n_in_universe"] / row["gene_count"] if row["gene_count"] > 0 else 0.0,
        axis=1,
    )

    return result.reset_index(drop=True)
