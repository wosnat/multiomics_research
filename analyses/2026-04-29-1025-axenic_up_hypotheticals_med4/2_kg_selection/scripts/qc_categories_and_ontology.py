"""
QC: gene_category distribution + candidate-set ontology landscape (term-level).

Three independent QC outputs:

1. gene_category cross-tab — sig_up (n=405) vs candidate AQ<=1 (n=116) vs
   non-candidate AQ>=2 (n=289), sorted by sig_up count desc. Also a bar
   figure highlighting the non-Unknown candidates.

2. Per-source ontology coverage on the candidate set (descriptive only —
   does NOT classify the source itself as informative or not). Saved as
   per-gene flags. The earlier source-level "hypclass-only" claim is
   retracted; term-level content evaluation belongs to step 4/5 dossier
   construction.

3. Candidate-set ontology landscape via gene_ontology_terms envelope —
   uses summary=True call. Reports by_term (term frequency) and
   by_ontology (term + gene counts per source) directly from the envelope,
   no per-row processing. Full term-list saved to CSV.

Inputs (default: ../data/):
    sig_up_405_overview.csv
    candidate_set.csv

Outputs (default: ../data/, ../figures/):
    data/qc_gene_category_dist.csv          gene_category cross-tab
    data/qc_ontology_per_source_flags.csv   per-candidate per-source presence flags
    data/qc_candidate_ontology_terms.csv    full gene × term rows for the candidate set (one row per gene × term)
    data/qc_candidate_ontology_terms_by_term.csv  by_term envelope: term × source × n_genes
    data/qc_candidate_ontology_by_source.csv      by_ontology envelope: source × n_terms × n_genes
    figures/qc_gene_category_dist.png       sig_up vs candidate breakdown by gene_category
    figures/qc_candidate_ontology_terms_top.png   top-N candidate-set terms by gene-frequency
"""

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from multiomics_explorer import gene_ontology_terms

ALL_ONTOLOGY_SOURCES = [
    "go_bp", "go_mf", "go_cc", "kegg", "ec", "cog_category",
    "cyanorak_role", "tigr_role", "pfam", "brite",
]


def types_set(s) -> set[str]:
    if pd.isna(s) or not s:
        return set()
    return {t.strip() for t in str(s).split("|") if t.strip()}


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data-dir", type=Path, default=Path(__file__).resolve().parent.parent / "data")
    parser.add_argument("--fig-dir", type=Path, default=Path(__file__).resolve().parent.parent / "figures")
    parser.add_argument("--organism", default="MED4")
    parser.add_argument("--top-terms", type=int, default=20)
    args = parser.parse_args()
    args.fig_dir.mkdir(parents=True, exist_ok=True)

    sig_up = pd.read_csv(args.data_dir / "sig_up_405_overview.csv")
    cand = pd.read_csv(args.data_dir / "candidate_set.csv")
    cand_lt = set(cand["locus_tag"])

    # --- 1. gene_category cross-tab + figure ---
    cat_full = sig_up["gene_category"].fillna("(null)").value_counts()
    sig_up_cand = sig_up[sig_up["locus_tag"].isin(cand_lt)]
    sig_up_noncand = sig_up[~sig_up["locus_tag"].isin(cand_lt)]
    cat_cand = sig_up_cand["gene_category"].fillna("(null)").value_counts()
    cat_noncand = sig_up_noncand["gene_category"].fillna("(null)").value_counts()
    cats = sorted(set(cat_full.index) | set(cat_cand.index) | set(cat_noncand.index),
                  key=lambda c: -cat_full.get(c, 0))
    cat_df = pd.DataFrame({
        "gene_category": cats,
        "sig_up_n405": [int(cat_full.get(c, 0)) for c in cats],
        "candidate_aq_le_1_n116": [int(cat_cand.get(c, 0)) for c in cats],
        "noncandidate_aq_ge_2_n289": [int(cat_noncand.get(c, 0)) for c in cats],
    })
    cat_df["candidate_share_of_sig_up"] = (cat_df["candidate_aq_le_1_n116"] / cat_df["sig_up_n405"]).round(3)
    cat_path = args.data_dir / "qc_gene_category_dist.csv"
    cat_df.to_csv(cat_path, index=False)
    print(f"wrote {cat_path}")

    # gene_category figure: stacked horizontal bars, candidate (AQ<=1) vs non-candidate (AQ>=2).
    # All categories shown, even zero-candidate ones, so the reader can see which functional categories contribute zero hypotheticals.
    fig, ax = plt.subplots(figsize=(9, 7))
    y = np.arange(len(cat_df))[::-1]  # top of plot = most-populated category
    cand_vals = cat_df["candidate_aq_le_1_n116"].values
    noncand_vals = cat_df["noncandidate_aq_ge_2_n289"].values
    ax.barh(y, cand_vals, color="#d62728", label=f"candidate AQ<=1 (n={len(cand)})")
    ax.barh(y, noncand_vals, left=cand_vals, color="#aec7e8", label="non-candidate AQ>=2")
    ax.set_yticks(y)
    ax.set_yticklabels(cat_df["gene_category"], fontsize=9)
    ax.set_xlabel("# sig_up genes in category")
    ax.set_title("Sig_up gene_category distribution: candidate (AQ<=1) vs non-candidate (AQ>=2)")
    for i, (yi, c, nc) in enumerate(zip(y, cand_vals, noncand_vals)):
        if c > 0:
            ax.text(c / 2, yi, f"{c}", ha="center", va="center", fontsize=8, color="white", fontweight="bold")
        if c + nc > 0:
            ax.text(c + nc + 1, yi, f"{c + nc}", ha="left", va="center", fontsize=8)
    ax.legend(loc="lower right")
    fig.tight_layout()
    fig.savefig(args.fig_dir / "qc_gene_category_dist.png", dpi=300)
    plt.close(fig)
    print(f"wrote {args.fig_dir / 'qc_gene_category_dist.png'}")

    # --- 2. per-source presence flags on candidate set (DESCRIPTIVE ONLY) ---
    # Earlier source-level "hypclass-only" classification is retracted.
    # Here we just record which sources are populated per gene so step 4/5
    # dossier construction can decide informativeness term-by-term.
    cand_flags = cand[["locus_tag", "log2fc", "annotation_quality", "annotation_types"]].copy()
    cand_flags["types"] = cand_flags["annotation_types"].apply(types_set)
    cand_flags["has_any_ontology_row"] = cand_flags["types"].apply(bool)
    for src in ALL_ONTOLOGY_SOURCES:
        cand_flags[f"has_{src}"] = cand_flags["types"].apply(lambda s, k=src: k in s)
    cand_flags = cand_flags.drop(columns=["types"])
    flags_path = args.data_dir / "qc_ontology_per_source_flags.csv"
    cand_flags.to_csv(flags_path, index=False)
    print(f"wrote {flags_path}")

    # --- 3. Candidate-set ontology landscape via gene_ontology_terms ---
    locus_tags = cand["locus_tag"].tolist()
    print(f"\nCalling gene_ontology_terms (full) on {len(locus_tags)} candidates...")
    full = gene_ontology_terms(
        locus_tags=locus_tags,
        organism=args.organism,
        limit=None,
    )
    rows_df = pd.DataFrame(full.get("results", []))
    if not rows_df.empty:
        rows_df.to_csv(args.data_dir / "qc_candidate_ontology_terms.csv", index=False)
        print(f"wrote {args.data_dir / 'qc_candidate_ontology_terms.csv'} ({len(rows_df)} gene × term rows)")
    # Use the same envelope for by_term / by_ontology — already populated regardless of summary mode.
    env = full
    print("envelope keys:", list(env.keys()))

    by_term_df = pd.DataFrame(env.get("by_term", []))
    if not by_term_df.empty:
        by_term_df = by_term_df.sort_values("count", ascending=False)
        by_term_df.to_csv(args.data_dir / "qc_candidate_ontology_terms_by_term.csv", index=False)
        print(f"wrote {args.data_dir / 'qc_candidate_ontology_terms_by_term.csv'} ({len(by_term_df)} unique terms)")

    by_ont_df = pd.DataFrame(env.get("by_ontology", []))
    if not by_ont_df.empty:
        by_ont_df = by_ont_df.sort_values("term_count", ascending=False)
        by_ont_df.to_csv(args.data_dir / "qc_candidate_ontology_by_source.csv", index=False)
        print(f"wrote {args.data_dir / 'qc_candidate_ontology_by_source.csv'}")

    # Top-terms figure
    if not by_term_df.empty:
        top = by_term_df.head(args.top_terms)
        fig, ax = plt.subplots(figsize=(11, 0.4 * len(top) + 1))
        y = np.arange(len(top))[::-1]
        # Color by ontology_type
        unique_ont = sorted(top["ontology_type"].unique())
        cmap = plt.get_cmap("tab10")
        color_for = {ont: cmap(i % 10) for i, ont in enumerate(unique_ont)}
        colors = [color_for[ont] for ont in top["ontology_type"]]
        ax.barh(y, top["count"].values, color=colors)
        labels = [f"[{r.ontology_type}] {r.term_name}" for r in top.itertuples()]
        ax.set_yticks(y)
        ax.set_yticklabels(labels, fontsize=8)
        ax.set_xlabel("# candidate genes carrying this term")
        ax.set_title(f"Top {len(top)} ontology terms in candidate set (n=116)\n"
                     f"total_matching={env['total_matching']} total_terms={env['total_terms']} "
                     f"terms_per_gene median={env['terms_per_gene_median']}")
        # Legend by ontology source
        from matplotlib.patches import Patch
        ax.legend(handles=[Patch(color=color_for[o], label=o) for o in unique_ont], loc="lower right")
        for yi, n in zip(y, top["count"].values):
            ax.text(n + 0.5, yi, f"{int(n)}", va="center", fontsize=8)
        fig.tight_layout()
        fig.savefig(args.fig_dir / "qc_candidate_ontology_terms_top.png", dpi=300)
        plt.close(fig)
        print(f"wrote {args.fig_dir / 'qc_candidate_ontology_terms_top.png'}")

    # --- Print summary ---
    n = len(cand)
    print(f"\n=== Candidate set ontology coverage (n={n}) ===")
    print(f"  Genes with at least one ontology row: {int(cand_flags['has_any_ontology_row'].sum())} ({cand_flags['has_any_ontology_row'].mean():.0%})")
    print(f"  Per-source presence:")
    for src in ALL_ONTOLOGY_SOURCES:
        col = f"has_{src}"
        n_with = int(cand_flags[col].sum())
        if n_with:
            print(f"    {src:14s}: {n_with:3d} ({n_with / n:.0%})")
    print(f"\n=== Envelope landscape ===")
    print(f"  total_matching (gene-term rows): {env['total_matching']}")
    print(f"  total_terms (distinct terms):    {env['total_terms']}")
    print(f"  total_genes (with terms):        {env['total_genes']}")
    print(f"  terms_per_gene min/median/max:   {env['terms_per_gene_min']} / {env['terms_per_gene_median']} / {env['terms_per_gene_max']}")
    if not by_term_df.empty:
        print(f"\n  Top {min(10, len(by_term_df))} terms in candidate set:")
        for r in by_term_df.head(10).itertuples():
            print(f"    [{r.ontology_type:14s}] {r.count:3d} genes — {r.term_name}")


if __name__ == "__main__":
    main()
