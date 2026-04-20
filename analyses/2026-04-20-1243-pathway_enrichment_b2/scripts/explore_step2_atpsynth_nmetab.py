"""Drill-down on ATP synthase (J.1), Ox-phos (ko00190), and N-metab (E.4).

For each, for a representative R cluster:
- Which genes (locus_tag + gene_name) are in the DE foreground ∩ pathway
- Pathway size in MED4
- Fold enrichment, p_adjust

Plus cross-pathway gene overlap: J.1 vs E.4, J.1 vs ko00190, ko00190 vs E.4.
"""
from __future__ import annotations

import pickle
from pathlib import Path

ANALYSIS_DIR = Path(__file__).resolve().parent.parent


def main() -> None:
    with open(ANALYSIS_DIR / "data" / "enrichment_results.pkl", "rb") as f:
        results = pickle.load(f)

    cyano_ts = results[("Prochlorococcus MED4", "cyanorak_role", "table_scope")]
    kegg_ts = results[("Prochlorococcus MED4", "kegg", "table_scope")]

    # Pick Tolonen 12h for both down (ribosome/PS signal) and up (N-metab signal).
    tol = "10.1038/msb4100087_nitrogen_nitrogen_deprivation_med4_med4_microarray"
    c_down = f"{tol}|12h|down"
    c_up = f"{tol}|12h|up"

    def dump(result, cluster: str, term: str, label: str) -> set[str]:
        print(f"\n--- {label}: cluster={cluster} term={term} ---")
        if cluster not in set(result.results["cluster"].unique()):
            print(f"cluster not in result")
            return set()
        explanation = result.explain(cluster, term)
        try:
            print(explanation._repr_markdown_())
        except Exception:
            print(explanation)
        # Pull the overlap gene set from result.inputs.gene_sets + term2gene.
        # Alternatively: read back from overlap_genes if accessor exists.
        try:
            overlap = result.overlap_genes(cluster, term)
        except Exception:
            overlap = None
        if overlap is not None:
            tags = [g.locus_tag for g in overlap]
            return set(tags)
        return set()

    genes_j1_down = dump(cyano_ts, c_down, "cyanorak.role:J.1",
                          "ATP synthase DOWN (Tolonen 12h down)")
    genes_koOP_down = dump(kegg_ts, c_down, "kegg.pathway:ko00190",
                            "Oxidative phosphorylation DOWN (Tolonen 12h down)")
    genes_e4_up = dump(cyano_ts, c_up, "cyanorak.role:E.4",
                        "N-metabolism UP (Tolonen 12h up)")

    # Now pathway-pathway overlap — what are the pathway member sets (not
    # cluster-specific)?  Use term2gene on the EnrichmentResult.
    print("\n\n=== Pathway member-gene overlap (MED4 pathway membership, NOT "
          "cluster-specific) ===")

    def pathway_members(result, term: str) -> set[str]:
        t2g = result.term2gene  # long DataFrame (locus_tag × term_id rows)
        sub = t2g[t2g["term_id"] == term]
        return set(sub["locus_tag"])

    def pathway_members_with_names(result, term: str) -> list[tuple[str, str]]:
        t2g = result.term2gene
        sub = t2g[t2g["term_id"] == term]
        return list(zip(sub["locus_tag"], sub["gene_name"].fillna("")))

    j1 = pathway_members(cyano_ts, "cyanorak.role:J.1")
    e4 = pathway_members(cyano_ts, "cyanorak.role:E.4")
    koOP = pathway_members(kegg_ts, "kegg.pathway:ko00190")
    ko195 = pathway_members(kegg_ts, "kegg.pathway:ko00195")

    def fmt(pairs):
        return ", ".join(f"{n}({t})" if n else t for t, n in sorted(pairs))

    j1_p = pathway_members_with_names(cyano_ts, "cyanorak.role:J.1")
    e4_p = pathway_members_with_names(cyano_ts, "cyanorak.role:E.4")
    koOP_p = pathway_members_with_names(kegg_ts, "kegg.pathway:ko00190")
    ko195_p = pathway_members_with_names(kegg_ts, "kegg.pathway:ko00195")

    print(f"\nJ.1 ATP synthase ({len(j1)}): {fmt(j1_p)}")
    print(f"\nE.4 N metab ({len(e4)}): {fmt(e4_p)}")
    print(f"\nko00190 Ox phos ({len(koOP)}): {fmt(koOP_p)}")
    print(f"\nko00195 Photosynthesis ({len(ko195)}): {fmt(ko195_p)}")

    print("\n--- pairwise overlaps (locus_tag sets) ---")
    print(f"J.1 ∩ E.4: {sorted(j1 & e4)}")
    print(f"J.1 ∩ ko00190: {sorted(j1 & koOP)}")
    print(f"J.1 ∩ ko00195: {sorted(j1 & ko195)}")
    print(f"ko00190 ∩ E.4: {sorted(koOP & e4)}")
    print(f"ko00190 ∩ ko00195: {sorted(koOP & ko195)}")
    print(f"E.4 ∩ ko00195: {sorted(e4 & ko195)}")


if __name__ == "__main__":
    main()
