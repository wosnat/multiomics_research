"""DAG-aware gene annotation roll-up to a target ontology level.

Public API:
    roll_up_to_level(annotations_df, hierarchy_df, target_level) -> DataFrame

The hierarchy_df encodes a DAG where each row is a (child, parent) edge with
level numbers.  Roll-up propagates gene annotations UPWARD only: a gene
annotated to a term at level N is projected to all ancestors at target_level.
There is NO downward propagation.

Deduplication: if a gene reaches the same ancestor via multiple paths (DAG
convergence), it is recorded once per ancestor term.
"""

from collections import defaultdict

import pandas as pd


def roll_up_to_level(
    annotations_df: pd.DataFrame,
    hierarchy_df: pd.DataFrame,
    target_level: int,
) -> pd.DataFrame:
    """Roll gene annotations up through a DAG to a target level.

    Parameters
    ----------
    annotations_df:
        DataFrame with columns: locus_tag, term_id.
        A gene may appear in multiple rows (one per annotation).
    hierarchy_df:
        DataFrame with columns: child_id, parent_id, child_level, parent_level.
        Each row is one directed edge in the ontology DAG.
    target_level:
        Integer level to project annotations to.
        Lower numbers are closer to the root (level 0 = root).

    Returns
    -------
    DataFrame with columns: locus_tag, term_id — deduplicated at target_level.
    Only genes that have at least one annotation reachable at target_level are
    included (no downward propagation; annotations below target_level are not
    in the output).
    """
    # Build term-level lookup: term_id -> level
    # Include both child and parent columns to cover all terms
    term_level: dict[str, int] = {}
    for _, row in hierarchy_df.iterrows():
        term_level[row["child_id"]] = int(row["child_level"])
        term_level[row["parent_id"]] = int(row["parent_level"])

    # Build upward adjacency: term -> set of direct parents
    parents: dict[str, set] = defaultdict(set)
    for _, row in hierarchy_df.iterrows():
        parents[row["child_id"]].add(row["parent_id"])

    # Build ancestor map: term -> set of ancestors AT target_level
    # Uses memoization to avoid redundant traversal.
    _cache: dict[str, frozenset] = {}

    def ancestors_at_level(term: str) -> frozenset:
        if term in _cache:
            return _cache[term]

        own_level = term_level.get(term)

        # If this term IS at the target level, it maps to itself
        if own_level == target_level:
            _cache[term] = frozenset({term})
            return _cache[term]

        # If this term is BELOW the target level (higher number), walk up
        if own_level is not None and own_level > target_level:
            result: set = set()
            for parent in parents.get(term, set()):
                result |= ancestors_at_level(parent)
            _cache[term] = frozenset(result)
            return _cache[term]

        # If this term is ABOVE the target level (lower number, closer to root),
        # it has no representation at target_level — no downward propagation
        _cache[term] = frozenset()
        return _cache[term]

    # Project each annotation
    output_rows = []
    for _, row in annotations_df.iterrows():
        gene = row["locus_tag"]
        term = row["term_id"]
        for ancestor in ancestors_at_level(term):
            output_rows.append({"locus_tag": gene, "term_id": ancestor})

    if not output_rows:
        return pd.DataFrame(columns=["locus_tag", "term_id"])

    result_df = pd.DataFrame(output_rows).drop_duplicates()
    return result_df.reset_index(drop=True)
