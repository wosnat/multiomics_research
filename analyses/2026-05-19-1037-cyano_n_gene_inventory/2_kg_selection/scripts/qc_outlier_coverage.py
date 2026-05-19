"""Investigate two strains flagged in step 1 for low Cyanorak N-gene coverage:

 - SS120 (LLII Prochlorococcus): 12 anchor genes vs 24-28 for other Prochlorococcus
 - WH7803 (marine Synechococcus): 27 anchor genes vs 37-43 for other marine Syn

Question: is the low count a real biological difference (e.g., reduced
gene complement), an annotation gap (genes present but not assigned to
E.4/D.1.3), or a different gene set entirely?

Approach: for each outlier strain, list its anchor genes; for each anchor
gene the *median peer* strain has but the outlier doesn't, check whether
that gene exists in the outlier's proteome under a different annotation
(via free-text search on gene name / product).

Inputs:
    data/01_anchor_genes.csv (built by 01_anchor_genes.py)
Outputs:
    data/qc_outlier_coverage.csv     gene-level table of presence/absence + free-text hits
    data/qc_outlier_coverage.log
"""
from __future__ import annotations

import logging
from pathlib import Path

import pandas as pd
from multiomics_explorer import genes_by_function

OUT_DIR = Path(__file__).resolve().parents[1] / "data"
LOG_PATH = OUT_DIR / "qc_outlier_coverage.log"

# Outlier → peers (same lineage subset) for the comparison.
COMPARISONS = [
    {
        "outlier": "Prochlorococcus marinus subsp. marinus CCMP1375 (SS120)",
        "outlier_short": "SS120",
        "peers": [  # other Prochlorococcus strains with Cyanorak anchor data
            "Prochlorococcus AS9601", "Prochlorococcus MED4", "Prochlorococcus MIT9301",
            "Prochlorococcus MIT9303", "Prochlorococcus MIT9312", "Prochlorococcus MIT9313",
            "Prochlorococcus NATL1A", "Prochlorococcus NATL2A",
        ],
    },
    {
        "outlier": "Synechococcus WH7803",
        "outlier_short": "WH7803",
        "peers": ["Synechococcus WH8102", "Synechococcus CC9311", "Synechococcus sp. BL107"],
    },
]


def gene_name_in_strain(gene_name: str, strain: str) -> bool:
    """Check if a gene name exists in the strain's proteome via genes_by_function."""
    if not gene_name or pd.isna(gene_name):
        return False
    try:
        result = genes_by_function(
            search_text=f'gene_name:"{gene_name}"',
            organism=strain,
            limit=5,
        )
    except Exception:
        # fall back to plain text search
        try:
            result = genes_by_function(search_text=gene_name, organism=strain, limit=5)
        except Exception:
            return False
    # Conservative: only count exact gene_name match
    for r in result.get("results", []):
        if (r.get("gene_name") or "").lower() == gene_name.lower():
            return True
    return False


def main() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s %(levelname)s %(message)s",
        handlers=[logging.FileHandler(LOG_PATH, mode="w"), logging.StreamHandler()],
    )

    anchor = pd.read_csv(OUT_DIR / "01_anchor_genes.csv")

    rows = []
    for comp in COMPARISONS:
        outlier = comp["outlier"]
        outlier_short = comp["outlier_short"]
        peers = comp["peers"]

        outlier_genes = set(anchor[anchor["strain"] == outlier]["gene_name"].dropna())
        # Peer "consensus" set: gene names present in at least half the peers
        peer_df = anchor[anchor["strain"].isin(peers) & anchor["gene_name"].notna()]
        peer_strain_count = peer_df.groupby("gene_name")["strain"].nunique()
        threshold = max(1, len(peers) // 2)
        peer_consensus = set(peer_strain_count[peer_strain_count >= threshold].index)

        missing = peer_consensus - outlier_genes
        logging.info("%s: %d anchor genes; peer-consensus (in >=%d of %d peers): %d; missing from outlier: %d",
                     outlier_short, len(outlier_genes), threshold, len(peers),
                     len(peer_consensus), len(missing))

        for gene_name in sorted(missing):
            n_peers_with = int(peer_strain_count.get(gene_name, 0))
            present_under_other_annotation = gene_name_in_strain(gene_name, outlier)
            rows.append({
                "outlier_strain_short": outlier_short,
                "outlier_strain": outlier,
                "gene_name": gene_name,
                "n_peers_with_gene": n_peers_with,
                "n_peers_total": len(peers),
                "present_in_outlier_proteome_by_name": present_under_other_annotation,
            })

    out = pd.DataFrame(rows)
    out.to_csv(OUT_DIR / "qc_outlier_coverage.csv", index=False)

    logging.info("---")
    for outlier_short, group in out.groupby("outlier_strain_short"):
        n_missing = len(group)
        n_recovered = int(group["present_in_outlier_proteome_by_name"].sum())
        logging.info("%s: %d genes missing from anchor; %d found in proteome via gene_name search (annotation gap), %d truly absent",
                     outlier_short, n_missing, n_recovered, n_missing - n_recovered)
        for _, row in group.iterrows():
            tag = "annotation gap" if row["present_in_outlier_proteome_by_name"] else "absent"
            logging.info("  %s: %s (%d/%d peers have it) — %s",
                         outlier_short, row["gene_name"], row["n_peers_with_gene"], row["n_peers_total"], tag)
    logging.info("Wrote %s", OUT_DIR / "qc_outlier_coverage.csv")


if __name__ == "__main__":
    main()
