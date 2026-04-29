"""
Empirical discovery of each dossier axis on the candidate set + anchor genes.

For every axis the dossier card needs to cover (per qc_facet_availability_by_fc.png),
call the relevant discovery tool in verbose mode and observe what content actually
comes back. Framing decisions for the dossier surface emerge from this discovery,
not from an assumed schema.

Anchors (3): PMM0958 (driving example A — rich), PMM1898 (driving example B — F1 floor),
PMM0246 (positive control — ntcA, well-annotated N-regulator).

Per axis we capture:
- Anchor-level: full verbose response per anchor (saved as JSON for human review)
- Aggregate view across the candidate set (n=116) where useful (envelopes, breakdowns)

Outputs (default ../data/):
    discover_anchors.json       Full verbose responses per anchor × axis
    discover_anchors_summary.csv  Per-anchor × axis: n_rows, content sample
    discover_clusters_aggregate.csv  Cluster axis: per-analysis_id × candidate-count + curated-description-coverage
    discover_homologs_aggregate.csv  Ortholog axis: per-source × candidates with groups; specificity-rank distribution
    discover_dms_aggregate.csv  Derived-metric axis: per-metric_type × candidate-count
    discover_response_aggregate.csv  Response axis: per-treatment_type × candidate-count + breadth distribution
    01_discover_dossier_axes.log     runtime log
"""

import argparse
import json
import logging
from collections import Counter
from pathlib import Path

import pandas as pd

from multiomics_explorer import (
    gene_clusters_by_gene,
    gene_derived_metrics,
    gene_homologs,
    gene_ontology_terms,
    gene_overview,
    gene_response_profile,
)

ANCHORS = ["PMM0958", "PMM1898", "PMM0246"]


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data-dir", type=Path, default=Path(__file__).resolve().parent.parent / "data")
    parser.add_argument(
        "--candidate-csv",
        type=Path,
        default=Path(__file__).resolve().parent.parent.parent / "2_kg_selection" / "data" / "candidate_set.csv",
    )
    parser.add_argument("--organism", default="MED4")
    args = parser.parse_args()
    args.data_dir.mkdir(parents=True, exist_ok=True)

    log_path = args.data_dir / "01_discover_dossier_axes.log"
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s %(levelname)s %(message)s",
        handlers=[logging.FileHandler(log_path, mode="w"), logging.StreamHandler()],
    )
    log = logging.getLogger(__name__)

    cand = pd.read_csv(args.candidate_csv)
    cand_lts = cand["locus_tag"].tolist()
    log.info("Anchors: %s", ANCHORS)
    log.info("Candidate set: %d genes", len(cand_lts))

    # ============================================================
    # Anchor-level discovery (one verbose call per axis with all 3 anchors)
    # ============================================================
    anchor_responses: dict[str, dict] = {}

    log.info("[anchors] gene_overview verbose...")
    anchor_responses["gene_overview"] = gene_overview(locus_tags=ANCHORS, verbose=True)

    log.info("[anchors] gene_ontology_terms (organism=%s)...", args.organism)
    anchor_responses["gene_ontology_terms"] = gene_ontology_terms(
        locus_tags=ANCHORS, organism=args.organism, limit=None
    )

    log.info("[anchors] gene_clusters_by_gene verbose...")
    anchor_responses["gene_clusters_by_gene"] = gene_clusters_by_gene(
        locus_tags=ANCHORS, verbose=True, limit=None
    )

    log.info("[anchors] gene_derived_metrics verbose...")
    anchor_responses["gene_derived_metrics"] = gene_derived_metrics(
        locus_tags=ANCHORS, verbose=True, limit=None
    )

    log.info("[anchors] gene_homologs verbose...")
    anchor_responses["gene_homologs"] = gene_homologs(
        locus_tags=ANCHORS, verbose=True, limit=None
    )

    log.info("[anchors] gene_response_profile (default = all MED4 experiments)...")
    anchor_responses["gene_response_profile"] = gene_response_profile(
        locus_tags=ANCHORS, group_by="treatment_type"
    )

    out_json = args.data_dir / "discover_anchors.json"
    with out_json.open("w") as fh:
        json.dump(anchor_responses, fh, indent=2, default=str)
    log.info("wrote %s", out_json)

    # Per-anchor × axis summary table
    summary_rows = []
    for axis, resp in anchor_responses.items():
        rows = resp.get("results", [])
        for lt in ANCHORS:
            sub = [r for r in rows if r.get("locus_tag") == lt]
            summary_rows.append({
                "axis": axis,
                "locus_tag": lt,
                "n_rows": len(sub),
            })
    pd.DataFrame(summary_rows).to_csv(args.data_dir / "discover_anchors_summary.csv", index=False)
    log.info("wrote anchors_summary")

    # ============================================================
    # Aggregate views across the candidate set (n=116)
    # ============================================================

    # --- Cluster axis aggregate: per-analysis_id, how many candidates, curated coverage ---
    log.info("[candidate_set] gene_clusters_by_gene verbose...")
    clusters_resp = gene_clusters_by_gene(locus_tags=cand_lts, verbose=True, limit=None)
    cluster_rows = clusters_resp.get("results", [])
    log.info("    cluster rows for candidate set: %d (envelope total=%d)",
             len(cluster_rows), clusters_resp.get("total_matching", 0))
    log.info("    genes_with_clusters=%d genes_without_clusters=%d",
             clusters_resp.get("genes_with_clusters", 0),
             clusters_resp.get("genes_without_clusters", 0))

    if cluster_rows:
        df_clu = pd.DataFrame(cluster_rows)
        # Per analysis_id
        agg = df_clu.groupby("analysis_id").agg(
            analysis_name=("analysis_name", "first"),
            cluster_type=("cluster_type", "first"),
            cluster_method=("cluster_method", "first"),
            n_candidate_rows=("locus_tag", "count"),
            n_distinct_candidates=("locus_tag", "nunique"),
            n_distinct_clusters=("cluster_id", "nunique"),
        ).reset_index().sort_values("n_candidate_rows", ascending=False)
        # Add: of the distinct clusters in this analysis touched, how many have curated cluster_functional_description != N/A
        cur_desc_present = df_clu.assign(
            curated=df_clu["cluster_functional_description"].fillna("N/A").str.strip().str.upper().ne("N/A")
            & df_clu["cluster_functional_description"].notna()
        ).groupby("analysis_id")["curated"].sum().rename("rows_with_curated_func_desc")
        agg = agg.merge(cur_desc_present, on="analysis_id", how="left")
        agg.to_csv(args.data_dir / "discover_clusters_aggregate.csv", index=False)
        log.info("wrote clusters_aggregate (per-analysis_id)")
        log.info("\n%s", agg.to_string(index=False))

    # --- Ortholog axis aggregate ---
    log.info("[candidate_set] gene_homologs verbose...")
    hom_resp = gene_homologs(locus_tags=cand_lts, verbose=True, limit=None)
    hom_rows = hom_resp.get("results", [])
    log.info("    homolog rows for candidate set: %d (envelope total=%d)",
             len(hom_rows), hom_resp.get("total_matching", 0))
    log.info("    by_source: %s", hom_resp.get("by_source"))
    log.info("    no_groups: %s", hom_resp.get("no_groups", []))

    if hom_rows:
        df_hom = pd.DataFrame(hom_rows)
        # Specificity-rank distribution per source
        spec_dist = df_hom.groupby(["source", "taxonomic_level", "specificity_rank"]).size().reset_index(name="n_rows")
        spec_dist["n_distinct_candidates"] = df_hom.groupby(
            ["source", "taxonomic_level", "specificity_rank"]
        )["locus_tag"].nunique().values
        spec_dist = spec_dist.sort_values(["source", "specificity_rank"]).reset_index(drop=True)
        spec_dist.to_csv(args.data_dir / "discover_homologs_aggregate.csv", index=False)
        log.info("wrote homologs_aggregate")
        log.info("\n%s", spec_dist.to_string(index=False))

    # --- Derived-metric axis aggregate ---
    log.info("[candidate_set] gene_derived_metrics verbose...")
    dm_resp = gene_derived_metrics(locus_tags=cand_lts, verbose=True, limit=None)
    dm_rows = dm_resp.get("results", [])
    log.info("    DM rows for candidate set: %d (envelope total=%d)",
             len(dm_rows), dm_resp.get("total_matching", 0))
    log.info("    by_value_kind: %s", dm_resp.get("by_value_kind"))
    log.info("    by_metric_type: %s", dm_resp.get("by_metric_type"))
    log.info("    genes_with_metrics=%d genes_without_metrics=%d",
             dm_resp.get("genes_with_metrics", 0),
             dm_resp.get("genes_without_metrics", 0))

    if dm_rows:
        df_dm = pd.DataFrame(dm_rows)
        agg_dm = df_dm.groupby(["metric_type", "value_kind", "name"]).agg(
            publication_doi=("publication_doi", "first"),
            compartment=("compartment", "first"),
            n_rows=("locus_tag", "count"),
            n_distinct_candidates=("locus_tag", "nunique"),
        ).reset_index().sort_values("n_rows", ascending=False)
        agg_dm.to_csv(args.data_dir / "discover_dms_aggregate.csv", index=False)
        log.info("wrote dms_aggregate")
        log.info("\n%s", agg_dm.to_string(index=False))

    # --- Response-profile axis aggregate ---
    log.info("[candidate_set] gene_response_profile (default — all MED4)...")
    rp_resp = gene_response_profile(
        locus_tags=cand_lts, group_by="treatment_type", limit=len(cand_lts)
    )
    rp_rows = rp_resp.get("results", [])
    log.info("    response_profile rows for candidate set: %d (genes_with_response=%d)",
             len(rp_rows), rp_resp.get("genes_with_response", 0))

    if rp_rows:
        df_rp = pd.DataFrame(rp_rows)
        # Per-treatment_type: how many candidates respond / not respond / not_known
        all_treatments = set()
        for r in rp_rows:
            for k in ("groups_responded", "groups_not_responded",
                      "groups_tested_not_responded", "groups_not_known"):
                vals = r.get(k) or []
                if isinstance(vals, str):
                    # to_dataframe-style "a | b" — but this is raw API, expect list
                    vals = [v.strip() for v in vals.split("|") if v.strip()]
                all_treatments.update(vals)

        per_treatment_rows = []
        for t in sorted(all_treatments):
            n_responded = sum(1 for r in rp_rows if t in (r.get("groups_responded") or []))
            n_not_resp = sum(1 for r in rp_rows if t in (r.get("groups_not_responded") or []))
            n_tested_nr = sum(1 for r in rp_rows if t in (r.get("groups_tested_not_responded") or []))
            n_not_known = sum(1 for r in rp_rows if t in (r.get("groups_not_known") or []))
            per_treatment_rows.append({
                "treatment_type": t,
                "n_candidates_responded": n_responded,
                "n_candidates_not_responded": n_not_resp,
                "n_candidates_tested_not_responded": n_tested_nr,
                "n_candidates_not_known": n_not_known,
            })
        df_per_t = pd.DataFrame(per_treatment_rows).sort_values(
            "n_candidates_responded", ascending=False
        )
        df_per_t.to_csv(args.data_dir / "discover_response_aggregate.csv", index=False)
        log.info("wrote response_aggregate")
        log.info("\n%s", df_per_t.to_string(index=False))

        # Breadth distribution: how many treatments does each candidate respond to?
        breadth = Counter(len(r.get("groups_responded") or []) for r in rp_rows)
        log.info("breadth (n_treatments_responded): %s", sorted(breadth.items()))

    # ============================================================
    # Per-axis envelope summary for notebook
    # ============================================================
    log.info("=" * 60)
    log.info("ENVELOPE SUMMARY (candidate set, n=%d)", len(cand_lts))
    log.info("clusters: %d gene-cluster rows, %d distinct clusters; %d genes_with / %d genes_without",
             clusters_resp.get("total_matching", 0),
             clusters_resp.get("total_clusters", 0),
             clusters_resp.get("genes_with_clusters", 0),
             clusters_resp.get("genes_without_clusters", 0))
    log.info("homologs: %d gene-group rows; sources=%s; no_groups=%s",
             hom_resp.get("total_matching", 0),
             hom_resp.get("by_source", []),
             hom_resp.get("no_groups", []))
    log.info("derived_metrics: %d gene-DM rows; by_value_kind=%s",
             dm_resp.get("total_matching", 0),
             dm_resp.get("by_value_kind", []))
    log.info("response_profile: %d candidates with response data",
             rp_resp.get("genes_with_response", 0))


if __name__ == "__main__":
    main()
