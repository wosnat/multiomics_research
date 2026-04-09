"""KG data extraction for pathway enrichment analysis.

Public API:
    extract_annotations(locus_tags, ontology, conn=None) -> DataFrame
    extract_hierarchy(ontology, max_depth=4, conn=None) -> DataFrame

These functions are the bridge between the KG and the pure DataFrame
functions in hierarchy.py, survey.py, and enrichment.py.
"""

import pandas as pd

# ---------------------------------------------------------------------------
# Ontology config: edge names and node labels
# ---------------------------------------------------------------------------

ONTOLOGY_EDGES = {
    "go_bp": {
        "gene_rel": "Gene_involved_in_biological_process",
        "hierarchy_rels": [
            "Biological_process_is_a_biological_process",
            "Biological_process_part_of_biological_process",
        ],
        "node_label": "BiologicalProcess",
    },
    "go_mf": {
        "gene_rel": "Gene_enables_molecular_function",
        "hierarchy_rels": [
            "Molecular_function_is_a_molecular_function",
            "Molecular_function_part_of_molecular_function",
        ],
        "node_label": "MolecularFunction",
    },
    "go_cc": {
        "gene_rel": "Gene_located_in_cellular_component",
        "hierarchy_rels": [
            "Cellular_component_is_a_cellular_component",
            "Cellular_component_part_of_cellular_component",
        ],
        "node_label": "CellularComponent",
    },
    "kegg": {
        "gene_rel": "Gene_has_kegg_ko",
        "hierarchy_rels": ["Kegg_term_is_a_kegg_term"],
        "node_label": "KeggTerm",
    },
    "cyanorak_role": {
        "gene_rel": "Gene_has_cyanorak_role",
        "hierarchy_rels": ["Cyanorak_role_is_a_cyanorak_role"],
        "node_label": "CyanorakRole",
    },
    "tigr_role": {
        "gene_rel": "Gene_has_tigr_role",
        "hierarchy_rels": [],
        "node_label": "TigrRole",
    },
    "cog_category": {
        "gene_rel": "Gene_in_cog_category",
        "hierarchy_rels": [],
        "node_label": "CogFunctionalCategory",
    },
    "ec": {
        "gene_rel": "Gene_catalyzes_ec_number",
        "hierarchy_rels": ["Ec_number_is_a_ec_number"],
        "node_label": "EcNumber",
    },
    "pfam": {
        "gene_rel": "Gene_has_pfam",
        "hierarchy_rels": ["Pfam_in_pfam_clan"],
        "node_label": "Pfam",
    },
}

_ANNOTATION_COLS = ["locus_tag", "term_id", "term_name"]
_HIERARCHY_COLS = ["child_id", "parent_id", "child_level", "parent_level"]


# ---------------------------------------------------------------------------
# extract_annotations
# ---------------------------------------------------------------------------

def extract_annotations(
    locus_tags: list[str],
    ontology: str,
    conn=None,
) -> pd.DataFrame:
    """Extract gene × term leaf annotations from the KG.

    Parameters
    ----------
    locus_tags:
        List of locus tag strings to query.
    ontology:
        Ontology name (must be a key in ONTOLOGY_EDGES).
    conn:
        Optional KG connection. If None, a default connection is created.

    Returns
    -------
    DataFrame with columns: locus_tag, term_id, term_name.
    Empty DataFrame with correct columns if no annotations found.

    Notes
    -----
    Genes in `not_found` (unknown to the KG) are silently ignored.
    """
    from multiomics_explorer import gene_ontology_terms
    from multiomics_explorer.analysis import to_dataframe

    if not locus_tags:
        return pd.DataFrame(columns=_ANNOTATION_COLS)

    # Batch to avoid Neo4j memory limits on large gene lists
    batch_size = 500
    frames = []
    for i in range(0, len(locus_tags), batch_size):
        batch = locus_tags[i : i + batch_size]
        result = gene_ontology_terms(
            locus_tags=batch,
            ontology=ontology,
            limit=None,
            conn=conn,
        )
        if result.get("results"):
            frames.append(to_dataframe(result))

    if not frames:
        return pd.DataFrame(columns=_ANNOTATION_COLS)

    df = pd.concat(frames, ignore_index=True)
    return df[["locus_tag", "term_id", "term_name"]].drop_duplicates().reset_index(drop=True)


# ---------------------------------------------------------------------------
# extract_hierarchy
# ---------------------------------------------------------------------------

def extract_hierarchy(
    ontology: str,
    max_depth: int = 4,
    conn=None,
) -> pd.DataFrame:
    """Extract parent-child hierarchy edges from the KG.

    For flat ontologies (no hierarchy_rels), returns an empty DataFrame.

    For hierarchical ontologies, finds all direct parent-child edges and
    assigns level numbers (0 = root, higher numbers = deeper leaves).

    Parameters
    ----------
    ontology:
        Ontology name (must be a key in ONTOLOGY_EDGES).
    max_depth:
        Maximum depth for hierarchy traversal. Used as a safety bound for
        level assignment. Edges deeper than this are still included if the
        KG contains them, but level assignment caps at max_depth.
    conn:
        Optional KG connection. If None, a default connection is created.

    Returns
    -------
    DataFrame with columns: child_id, parent_id, child_level, parent_level.
    Empty DataFrame with correct columns if ontology is flat or has no hierarchy.

    Notes
    -----
    Level assignment is distance-from-root: nodes with no outgoing hierarchy
    edge are roots at level 0. Their children are at level 1, and so on.
    For KEGG, the `level` property on KeggTerm nodes is used directly.
    For CyanoRak, the level is derived from the dot-count in the `code` property.
    For all others, BFS from roots is used.
    """
    from multiomics_explorer import run_cypher

    config = ONTOLOGY_EDGES.get(ontology)
    if config is None:
        raise ValueError(
            f"Unknown ontology: {ontology!r}. Valid: {sorted(ONTOLOGY_EDGES)}"
        )

    hierarchy_rels = config["hierarchy_rels"]
    if not hierarchy_rels:
        # Flat ontology — no hierarchy
        return pd.DataFrame(columns=_HIERARCHY_COLS)

    node_label = config["node_label"]

    # Use ontology-specific level derivation where available
    if ontology == "kegg":
        return _extract_hierarchy_kegg(hierarchy_rels, node_label, conn)
    elif ontology == "cyanorak_role":
        return _extract_hierarchy_cyanorak(hierarchy_rels, node_label, conn)
    else:
        return _extract_hierarchy_bfs(hierarchy_rels, node_label, conn)


# ---------------------------------------------------------------------------
# Hierarchy extraction strategies
# ---------------------------------------------------------------------------

def _build_rel_pattern(hierarchy_rels: list[str]) -> str:
    """Build a Cypher relationship type pattern for one or more rel types."""
    if len(hierarchy_rels) == 1:
        return f"[:{hierarchy_rels[0]}]"
    types = "|".join(f":{r}" for r in hierarchy_rels)
    return f"[{types}]"


def _extract_hierarchy_kegg(
    hierarchy_rels: list[str],
    node_label: str,
    conn,
) -> pd.DataFrame:
    """KEGG: use the `level` property on KeggTerm nodes.

    KEGG level property values: "category", "subcategory", "pathway", "ko".
    Map to integers: category=0, subcategory=1, pathway=2, ko=3.
    """
    from multiomics_explorer import run_cypher

    rel_pattern = _build_rel_pattern(hierarchy_rels)
    query = (
        f"MATCH (child:{node_label})-{rel_pattern}->(parent:{node_label}) "
        f"RETURN child.id AS child_id, parent.id AS parent_id, "
        f"child.level AS child_level_str, parent.level AS parent_level_str"
    )

    result = run_cypher(query=query, limit=10000, conn=conn)
    rows = result.get("results", [])
    if not rows:
        return pd.DataFrame(columns=_HIERARCHY_COLS)

    _KEGG_LEVEL_MAP = {"category": 0, "subcategory": 1, "pathway": 2, "ko": 3}

    records = []
    for row in rows:
        child_level = _KEGG_LEVEL_MAP.get(row.get("child_level_str", ""), None)
        parent_level = _KEGG_LEVEL_MAP.get(row.get("parent_level_str", ""), None)
        if child_level is None or parent_level is None:
            continue
        records.append({
            "child_id": row["child_id"],
            "parent_id": row["parent_id"],
            "child_level": child_level,
            "parent_level": parent_level,
        })

    if not records:
        return pd.DataFrame(columns=_HIERARCHY_COLS)

    return pd.DataFrame(records).reset_index(drop=True)


def _extract_hierarchy_cyanorak(
    hierarchy_rels: list[str],
    node_label: str,
    conn,
) -> pd.DataFrame:
    """CyanoRak: derive level from dot count in `code` property.

    code "A" -> 0 dots -> level 0 (root)
    code "A.1" -> 1 dot -> level 1
    code "A.1.1" -> 2 dots -> level 2
    """
    from multiomics_explorer import run_cypher

    rel_pattern = _build_rel_pattern(hierarchy_rels)
    query = (
        f"MATCH (child:{node_label})-{rel_pattern}->(parent:{node_label}) "
        f"RETURN child.id AS child_id, parent.id AS parent_id, "
        f"child.code AS child_code, parent.code AS parent_code"
    )

    result = run_cypher(query=query, limit=10000, conn=conn)
    rows = result.get("results", [])
    if not rows:
        return pd.DataFrame(columns=_HIERARCHY_COLS)

    def code_to_level(code: str | None) -> int | None:
        if not code:
            return None
        return code.count(".")

    records = []
    for row in rows:
        child_level = code_to_level(row.get("child_code"))
        parent_level = code_to_level(row.get("parent_code"))
        if child_level is None or parent_level is None:
            continue
        records.append({
            "child_id": row["child_id"],
            "parent_id": row["parent_id"],
            "child_level": child_level,
            "parent_level": parent_level,
        })

    if not records:
        return pd.DataFrame(columns=_HIERARCHY_COLS)

    return pd.DataFrame(records).reset_index(drop=True)


def _extract_hierarchy_bfs(
    hierarchy_rels: list[str],
    node_label: str,
    conn,
) -> pd.DataFrame:
    """General BFS-based level assignment: roots are at level 0.

    Roots are nodes that have no outgoing hierarchy edge (no parent).
    BFS propagates levels downward from roots.
    """
    from multiomics_explorer import run_cypher

    rel_pattern = _build_rel_pattern(hierarchy_rels)

    # Fetch all direct edges
    query = (
        f"MATCH (child:{node_label})-{rel_pattern}->(parent:{node_label}) "
        f"RETURN child.id AS child_id, parent.id AS parent_id"
    )

    result = run_cypher(query=query, limit=10000, conn=conn)
    rows = result.get("results", [])
    if not rows:
        return pd.DataFrame(columns=_HIERARCHY_COLS)

    # Build adjacency
    children_of: dict[str, set[str]] = {}  # parent -> set of children
    parents_of: dict[str, set[str]] = {}   # child -> set of parents
    all_nodes: set[str] = set()

    for row in rows:
        child_id = row["child_id"]
        parent_id = row["parent_id"]
        all_nodes.add(child_id)
        all_nodes.add(parent_id)
        children_of.setdefault(parent_id, set()).add(child_id)
        parents_of.setdefault(child_id, set()).add(parent_id)

    # Roots: nodes that appear as parent but never as child
    # (i.e., no outgoing hierarchy edge from them)
    roots = {n for n in all_nodes if n not in parents_of}

    # BFS to assign levels
    from collections import deque

    node_level: dict[str, int] = {}
    queue: deque[tuple[str, int]] = deque()
    for root in roots:
        node_level[root] = 0
        queue.append((root, 0))

    while queue:
        node, level = queue.popleft()
        for child in children_of.get(node, set()):
            child_level = level + 1
            # If already visited, keep the minimum level (shallowest path)
            if child not in node_level or node_level[child] > child_level:
                node_level[child] = child_level
                queue.append((child, child_level))

    # Build output DataFrame
    records = []
    for row in rows:
        child_id = row["child_id"]
        parent_id = row["parent_id"]
        child_level = node_level.get(child_id)
        parent_level = node_level.get(parent_id)
        if child_level is None or parent_level is None:
            continue
        records.append({
            "child_id": child_id,
            "parent_id": parent_id,
            "child_level": child_level,
            "parent_level": parent_level,
        })

    if not records:
        return pd.DataFrame(columns=_HIERARCHY_COLS)

    return pd.DataFrame(records).reset_index(drop=True)
