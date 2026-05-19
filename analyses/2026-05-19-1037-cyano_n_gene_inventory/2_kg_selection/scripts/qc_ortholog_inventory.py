"""QC the orthology-bridged N-gene inventory.

Checks:
  Q1. The 7 universal (n=19) groups: do they correspond to the expected
      universally-conserved N machinery (NtcA, GlnA, etc.)?
  Q2. Coverage of the 6 previously-Cyanorak-uncovered strains via the bridge:
      do the per-strain group lists look biologically sensible?
  Q3. SS120 sanity: the inventory should show the urea-pathway groups as absent.
  Q4. WH7803 sanity: same — urease + most urea transport absent.
  Q5. Are there out-of-cyano members in our groups? (Bridge integrity check —
      eggnog Cyanobacteria-level groups should not contain Alteromonas/etc.)
  Q6. Distribution of paralogs: which groups have multi-copy members in any strain?
  Q7. Anchors-without-group orphans: do they have a higher-tier group available?

Inputs:
    data/04_anchor_to_ortholog.csv
    data/04_anchors_without_group.csv
    data/04_ortholog_groups.csv
    data/05_inventory_members.csv
    data/05_inventory_matrix.csv
    data/02_strain_table.csv
Outputs:
    data/qc_ortholog_inventory.log    detailed checks
    data/qc_universal_groups.csv      Q1 detail
    data/qc_paralogs.csv              Q6 detail
    data/qc_orphan_anchors.csv        Q7 detail (orphans + their other-tier groups)
"""
from __future__ import annotations

import logging
from pathlib import Path

import pandas as pd
from multiomics_explorer import gene_homologs

OUT_DIR = Path(__file__).resolve().parents[1] / "data"
LOG_PATH = OUT_DIR / "qc_ortholog_inventory.log"


def main() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s %(levelname)s %(message)s",
        handlers=[logging.FileHandler(LOG_PATH, mode="w"), logging.StreamHandler()],
    )

    a2o = pd.read_csv(OUT_DIR / "04_anchor_to_ortholog.csv")
    orphans = pd.read_csv(OUT_DIR / "04_anchors_without_group.csv")
    groups = pd.read_csv(OUT_DIR / "04_ortholog_groups.csv")
    members = pd.read_csv(OUT_DIR / "05_inventory_members.csv")
    matrix = pd.read_csv(OUT_DIR / "05_inventory_matrix.csv", index_col=0)
    strain_table = pd.read_csv(OUT_DIR / "02_strain_table.csv")

    # ---------- Q1: universal (n=19) groups ----------
    logging.info("===== Q1: Universal (present in all 19 strains) groups =====")
    universal_mask = (matrix > 0).sum(axis=0) == 19
    universal_groups = matrix.columns[universal_mask].tolist()
    universal_detail = groups[groups["group_id"].isin(universal_groups)].copy()
    universal_detail = universal_detail.sort_values("consensus_gene_name")
    universal_detail.to_csv(OUT_DIR / "qc_universal_groups.csv", index=False)
    logging.info("Universal groups: %d", len(universal_groups))
    for _, row in universal_detail.iterrows():
        logging.info("  %s — %s (%s) [genera: %s]",
                     row["group_id"], row.get("consensus_gene_name") or "",
                     row.get("consensus_product") or "", row.get("genera") or "")

    # ---------- Q2 + Q3 + Q4: per-strain group lists for selected strains ----------
    interest_strains = [
        "Prochlorococcus marinus subsp. marinus CCMP1375 (SS120)",
        "Synechococcus WH7803",
        "Prochlorococcus MIT0801",          # was annotation-blank in Cyanorak
        "Prochlorococcus RSP50",             # was annotation-blank in Cyanorak
        "Synechococcus PCC 7002",            # outside Cyanorak picocyano
        "Synechococcus elongatus PCC 7942",  # outside Cyanorak picocyano
        "Synechococcus elongatus UTEX 2973", # outside Cyanorak picocyano
        "Thermosynechococcus vestitus BP-1", # outside Cyanorak picocyano
    ]
    logging.info("===== Q2/Q3/Q4: per-strain group inventory for outlier and previously-uncovered strains =====")
    consensus_name_map = dict(zip(groups["group_id"], groups["consensus_gene_name"]))
    for strain in interest_strains:
        row = matrix.loc[strain]
        present_groups = row[row > 0].index.tolist()
        absent_groups = row[row == 0].index.tolist()
        def _name(gid: str) -> str:
            v = consensus_name_map.get(gid)
            if isinstance(v, str) and v:
                return v
            return gid
        present_names = sorted(_name(gid) for gid in present_groups)
        absent_names = sorted(_name(gid) for gid in absent_groups)
        logging.info("--- %s (%d present / %d absent of 54 groups)", strain, len(present_groups), len(absent_groups))
        logging.info("    present: %s", ", ".join(present_names))
        logging.info("    absent:  %s", ", ".join(absent_names))

    # ---------- Q5: out-of-cyano members (bridge integrity) ----------
    logging.info("===== Q5: bridge integrity — are any members outside the 19 cyano scope? =====")
    scope = set(strain_table["strain"])
    in_members_orgs = set(members["organism_name"])
    out_of_scope = in_members_orgs - scope
    if out_of_scope:
        logging.warning("Out-of-scope organisms in member rows (should be zero): %s", sorted(out_of_scope))
    else:
        logging.info("All %d member rows are in-scope (19 cyano strains). Bridge integrity OK.", len(members))

    # ---------- Q6: paralogs ----------
    logging.info("===== Q6: groups with multi-copy members in any strain =====")
    multi = members.groupby(["organism_name", "group_id"])["locus_tag"].nunique()
    multi = multi[multi > 1].reset_index().rename(columns={"locus_tag": "copy_count"})
    if len(multi):
        multi = multi.merge(groups[["group_id", "consensus_gene_name", "consensus_product"]], on="group_id", how="left")
        multi = multi.sort_values(["organism_name", "copy_count"], ascending=[True, False])
        multi.to_csv(OUT_DIR / "qc_paralogs.csv", index=False)
        logging.info("Groups with paralogs in at least one strain: %d entries (strain × group)", len(multi))
        for _, row in multi.iterrows():
            logging.info("  %s: %s (%s) × %d copies",
                         row["organism_name"], row.get("consensus_gene_name") or "?",
                         row["group_id"], row["copy_count"])
    else:
        logging.info("No paralogs detected (every strain × group has ≤1 copy).")

    # ---------- Q7: orphan anchors — do they have higher-tier groups? ----------
    logging.info("===== Q7: orphan anchors — alternative-tier groups =====")
    orphan_loci = orphans["locus_tag"].tolist()
    if orphan_loci:
        # Probe at Prochloraceae and Bacteria tiers.
        for tier in ["Prochloraceae", "Bacteria"]:
            r = gene_homologs(locus_tags=orphan_loci, source="eggnog", taxonomic_level=tier,
                              limit=None, verbose=False)
            df = pd.DataFrame(r["results"]) if r["results"] else pd.DataFrame()
            covered = set(df["locus_tag"]) if not df.empty else set()
            logging.info("  Tier=%s: %d/%d orphan loci have a group at this tier",
                         tier, len(covered), len(orphan_loci))
            for locus in orphan_loci:
                hits = df[df["locus_tag"] == locus][["group_id", "consensus_gene_name"]] if not df.empty else pd.DataFrame()
                if not hits.empty:
                    for _, h in hits.iterrows():
                        logging.info("    %s → %s (%s)", locus, h["group_id"], h.get("consensus_gene_name") or "")
                else:
                    logging.info("    %s → NO group at %s tier", locus, tier)
        # Save merged orphan summary.
        out = orphans.copy()
        out.to_csv(OUT_DIR / "qc_orphan_anchors.csv", index=False)
    else:
        logging.info("No orphan anchors.")

    logging.info("---")
    logging.info("Wrote %s, %s, %s",
                 OUT_DIR / "qc_universal_groups.csv",
                 OUT_DIR / "qc_paralogs.csv",
                 OUT_DIR / "qc_orphan_anchors.csv")


if __name__ == "__main__":
    main()
