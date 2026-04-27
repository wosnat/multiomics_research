"""QC: per-experiment gene detection sets and pairwise overlap.

For each of the 8 selected experiments, queries the set of locus_tags with at
least one DE row via `differential_expression_by_gene(experiment_ids=[exp_id])`.
All 8 experiments use `table_scope == "all_detected_genes"`, so this set equals
the experiment's distinct_gene_count and is constant across timepoints within
the experiment.

For each (organism × condition) pairing of RNA-seq and proteomics experiments,
reports:
  - |RNA only|, |Protein only|, |both| (locus_tag intersection)
  - The intersection is the universe of genes available for paired-FC discordance.

Inputs:
  - data/paired_observations.csv  (from 02_build_paired_observations.py)
Outputs:
  - data/qc_gene_coverage.csv     per (organism, condition) coverage row
  - data/qc_gene_coverage.log

Run from the multiomics_research repo root:
  .venv/bin/python analyses/2026-04-27-1638-proteome_transcriptome_discordance/2_kg_selection/scripts/qc_gene_coverage.py
"""

from __future__ import annotations

from pathlib import Path
import sys

import pandas as pd

from multiomics_explorer import differential_expression_by_gene

DATA_DIR = Path(__file__).resolve().parents[1] / "data"
PAIRED_CSV = DATA_DIR / "paired_observations.csv"
OUT_CSV = DATA_DIR / "qc_gene_coverage.csv"
LOG_PATH = DATA_DIR / "qc_gene_coverage.log"


def log(msg: str, *, fh) -> None:
    print(msg)
    fh.write(msg + "\n")


def get_locus_tags(experiment_id: str) -> set[str]:
    """All distinct locus_tags with at least one DE row for the experiment.

    `differential_expression_by_gene(experiment_ids=[exp_id])` returns one row
    per (gene × timepoint) for the experiment; we collapse to the locus_tag set.
    The envelope's `matching_genes` field carries the same distinct count as a
    cross-check.
    """
    result = differential_expression_by_gene(
        experiment_ids=[experiment_id], limit=None
    )
    if result.get("truncated"):
        sys.exit(
            f"ERROR: differential_expression_by_gene returned truncated results "
            f"for {experiment_id}; rerun with limit=None"
        )
    locus_tags = {r["locus_tag"] for r in result["results"] if r.get("locus_tag")}
    # Cross-check against envelope
    expected = result.get("matching_genes")
    if expected is not None and expected != len(locus_tags):
        print(
            f"  WARNING: envelope matching_genes={expected} != distinct locus_tags={len(locus_tags)} "
            f"for {experiment_id}"
        )
    return locus_tags


def main() -> None:
    if not PAIRED_CSV.exists():
        sys.exit(f"ERROR: {PAIRED_CSV} not found. Run 02_build_paired_observations.py first.")

    paired = pd.read_csv(PAIRED_CSV)

    with LOG_PATH.open("w") as fh:
        # One (organism, condition) pair → one (rnaseq_exp, proteomics_exp) — all
        # the per-TP rows in the paired_observations CSV share the same exp ids.
        pairs = (
            paired[["organism_name", "condition", "rnaseq_experiment_id", "proteomics_experiment_id"]]
            .drop_duplicates()
            .reset_index(drop=True)
        )
        log(f"Found {len(pairs)} (organism, condition) experiment pairs", fh=fh)

        cache: dict[str, set[str]] = {}
        out_rows = []
        for _, p in pairs.iterrows():
            rna_id = p["rnaseq_experiment_id"]
            prot_id = p["proteomics_experiment_id"]

            log("", fh=fh)
            log(f"=== {p['organism_name']} / {p['condition']} ===", fh=fh)
            log(f"  RNAseq exp:     {rna_id}", fh=fh)
            log(f"  Proteomics exp: {prot_id}", fh=fh)

            if rna_id not in cache:
                cache[rna_id] = get_locus_tags(rna_id)
            if prot_id not in cache:
                cache[prot_id] = get_locus_tags(prot_id)

            rna_set = cache[rna_id]
            prot_set = cache[prot_id]
            both = rna_set & prot_set
            rna_only = rna_set - prot_set
            prot_only = prot_set - rna_set

            log(f"  |RNAseq detected|     = {len(rna_set):>5}", fh=fh)
            log(f"  |Proteomics detected| = {len(prot_set):>5}", fh=fh)
            log(f"  |Both (paired pool)|  = {len(both):>5}", fh=fh)
            log(f"  |RNA only|            = {len(rna_only):>5}", fh=fh)
            log(f"  |Protein only|        = {len(prot_only):>5}", fh=fh)
            log(f"  Pct of proteins also in RNA: "
                f"{len(both)/len(prot_set):.1%}" if prot_set else "  (no proteins detected)",
                fh=fh)

            out_rows.append({
                "organism_name": p["organism_name"],
                "condition": p["condition"],
                "rnaseq_experiment_id": rna_id,
                "proteomics_experiment_id": prot_id,
                "rnaseq_gene_count": len(rna_set),
                "proteomics_gene_count": len(prot_set),
                "intersection_gene_count": len(both),
                "rna_only_gene_count": len(rna_only),
                "prot_only_gene_count": len(prot_only),
                "pct_proteins_in_rna": (
                    len(both) / len(prot_set) if prot_set else float("nan")
                ),
            })

        out_df = pd.DataFrame(out_rows)
        out_df.to_csv(OUT_CSV, index=False)

        # Cross-condition consistency check: for each organism, are the
        # proteomics gene sets the same in axenic vs coculture? Same for RNAseq?
        log("", fh=fh)
        log("=== Cross-condition consistency (within organism) ===", fh=fh)
        for org in pairs["organism_name"].unique():
            rna_ids = pairs[pairs["organism_name"] == org]["rnaseq_experiment_id"].unique()
            prot_ids = pairs[pairs["organism_name"] == org]["proteomics_experiment_id"].unique()

            rna_sets = [cache[i] for i in rna_ids]
            prot_sets = [cache[i] for i in prot_ids]

            log(f"  {org}:", fh=fh)
            if len(rna_sets) >= 2:
                rna_inter = set.intersection(*rna_sets)
                rna_union = set.union(*rna_sets)
                log(f"    RNAseq:     |union|={len(rna_union)}  |intersection|={len(rna_inter)}  "
                    f"|axen-only|={len(rna_sets[0]-rna_sets[1])}  "
                    f"|cocul-only|={len(rna_sets[1]-rna_sets[0])}", fh=fh)
            if len(prot_sets) >= 2:
                prot_inter = set.intersection(*prot_sets)
                prot_union = set.union(*prot_sets)
                log(f"    Proteomics: |union|={len(prot_union)}  |intersection|={len(prot_inter)}  "
                    f"|axen-only|={len(prot_sets[0]-prot_sets[1])}  "
                    f"|cocul-only|={len(prot_sets[1]-prot_sets[0])}", fh=fh)

        log("", fh=fh)
        log(f"Wrote {OUT_CSV}", fh=fh)
        log(f"Wrote {LOG_PATH}", fh=fh)


if __name__ == "__main__":
    main()
