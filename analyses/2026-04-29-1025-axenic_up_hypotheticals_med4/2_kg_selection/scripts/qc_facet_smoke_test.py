"""
Per-gene KG facet smoke test on the top-3-by-log2fc candidates.

Each per-gene MCP function is called ONCE with all 3 locus_tags as a batch
(not 3 separate calls per function). Duplicate fields across responses
(gene_name, etc.) are tolerated — the goal is to confirm each facet
endpoint returns informative content for high-FC candidates before the
step 4/5 dossier construction commits to that surface.

Inputs (default: ../data/):
    candidate_set.csv         sorted by log2fc desc (top-3 picked from here)

Outputs (default: ../data/):
    qc_facet_smoke_test.json          full JSON of every facet response
    qc_facet_smoke_test_summary.csv   3-gene x N-facet matrix of "rows returned + 1-line content"
    qc_facet_smoke_test.log           runtime log
"""

import argparse
import json
import logging
from pathlib import Path

import pandas as pd

from multiomics_explorer import (
    gene_overview,
    gene_details,
    gene_ontology_terms,
    gene_clusters_by_gene,
    gene_derived_metrics,
    gene_homologs,
    gene_response_profile,
)


def _content_sample(facet: str, response: dict, locus_tag: str) -> str:
    """Return a 1-line content sample for the (facet, gene) cell."""
    rows = [r for r in response.get("results", []) if r.get("locus_tag") == locus_tag]

    if facet == "gene_overview":
        if not rows:
            return ""
        r = rows[0]
        return f"AQ={r.get('annotation_quality')} product='{r.get('product')}' annotation_types={r.get('annotation_types')}"

    if facet == "gene_details":
        if not rows:
            return ""
        r = rows[0]
        # gene_details returns ALL Gene node properties; pick a few sparse fields
        sparse = {k: v for k, v in r.items()
                  if k in ("ec_numbers", "transporter_classification", "cazy_ids",
                           "catalytic_activities", "kegg_ids", "cog_categories",
                           "pfam_ids", "signal_peptide", "transmembrane_regions")
                  and v not in (None, [], "")}
        return f"sparse_fields={sparse}" if sparse else "no sparse fields populated"

    if facet == "gene_ontology_terms":
        if not rows:
            return "no ontology terms"
        ont_types = sorted({r.get("ontology_type") for r in rows if r.get("ontology_type")})
        sample_terms = [f"{r.get('ontology_type')}:{r.get('term_name')}" for r in rows[:3]]
        return f"{len(rows)} terms across {ont_types}; sample={sample_terms}"

    if facet == "gene_clusters_by_gene":
        if not rows:
            return "no cluster memberships"
        cluster_names = [r.get("cluster_name") for r in rows]
        return f"{len(rows)} memberships: {cluster_names}"

    if facet == "gene_derived_metrics":
        if not rows:
            return "no derived metrics"
        names = [r.get("name") for r in rows]
        return f"{len(rows)} DMs: {names}"

    if facet == "gene_homologs":
        if not rows:
            return "no homolog groups"
        groups = [(r.get("source"), r.get("group_id"), r.get("taxonomic_level"), r.get("consensus_gene_name")) for r in rows]
        return f"{len(rows)} groups: {groups}"

    if facet == "gene_response_profile":
        if not rows:
            return "no response profile"
        r = rows[0]
        return (f"groups_responded={r.get('groups_responded')} "
                f"groups_not_responded={r.get('groups_not_responded')}")

    return ""


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data-dir", type=Path, default=Path(__file__).resolve().parent.parent / "data")
    parser.add_argument("--top-n", type=int, default=3)
    parser.add_argument("--organism", default="MED4")
    args = parser.parse_args()

    log_path = args.data_dir / "qc_facet_smoke_test.log"
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s %(levelname)s %(message)s",
        handlers=[logging.FileHandler(log_path, mode="w"), logging.StreamHandler()],
    )
    log = logging.getLogger(__name__)

    cand = pd.read_csv(args.data_dir / "candidate_set.csv")
    cand = cand.sort_values("log2fc", ascending=False)
    top = cand.head(args.top_n)
    locus_tags = top["locus_tag"].tolist()
    log.info("Top %d by log2fc: %s", args.top_n,
             list(zip(top["locus_tag"], top["log2fc"], top["annotation_quality"], top["product"])))

    # One batch call per facet.
    responses: dict[str, dict] = {}
    log.info("Calling gene_overview...")
    responses["gene_overview"] = gene_overview(locus_tags=locus_tags, verbose=True)
    log.info("Calling gene_details...")
    responses["gene_details"] = gene_details(locus_tags=locus_tags)
    log.info("Calling gene_ontology_terms...")
    responses["gene_ontology_terms"] = gene_ontology_terms(locus_tags=locus_tags, organism=args.organism, limit=None)
    log.info("Calling gene_clusters_by_gene...")
    responses["gene_clusters_by_gene"] = gene_clusters_by_gene(locus_tags=locus_tags, verbose=True, limit=None)
    log.info("Calling gene_derived_metrics...")
    responses["gene_derived_metrics"] = gene_derived_metrics(locus_tags=locus_tags, verbose=True, limit=None)
    log.info("Calling gene_homologs...")
    responses["gene_homologs"] = gene_homologs(locus_tags=locus_tags, verbose=True, limit=None)
    log.info("Calling gene_response_profile...")
    responses["gene_response_profile"] = gene_response_profile(locus_tags=locus_tags)

    # Dump full JSON.
    out_json = args.data_dir / "qc_facet_smoke_test.json"
    with out_json.open("w") as fh:
        json.dump(responses, fh, indent=2, default=str)
    log.info("wrote %s", out_json)

    # Build summary matrix.
    matrix_rows = []
    for facet, response in responses.items():
        for lt in locus_tags:
            n_rows = sum(1 for r in response.get("results", []) if r.get("locus_tag") == lt)
            matrix_rows.append({
                "facet": facet,
                "locus_tag": lt,
                "n_rows_for_gene": n_rows,
                "content_sample": _content_sample(facet, response, lt),
            })
    matrix = pd.DataFrame(matrix_rows)
    out_csv = args.data_dir / "qc_facet_smoke_test_summary.csv"
    matrix.to_csv(out_csv, index=False)
    log.info("wrote %s", out_csv)

    log.info("=" * 60)
    log.info("SMOKE TEST SUMMARY")
    for facet in responses:
        sub = matrix[matrix["facet"] == facet]
        per_gene = ", ".join(f"{r.locus_tag}={r.n_rows_for_gene}" for r in sub.itertuples())
        log.info("    %-25s rows_per_gene: %s", facet, per_gene)


if __name__ == "__main__":
    main()
