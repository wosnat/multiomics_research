"""Within-ontology pathway-gene-overlap audit for signature pathways.

For each ontology, compute pairwise Jaccard(term_i.genes, term_j.genes) on
MED4 gene membership (from result.term2gene) among signature pathways.
Flag: strict subset (one pathway ⊂ another) or Jaccard > 0.5.

Outputs:
    exploration/qc/step3_signature_redundancy_{ontology}.csv
    exploration/qc/step3_signature_redundancy_{ontology}.png  (heatmap)

Note: audit runs within ontology only — spec §5 Step 4 scores per-ontology,
so cross-ontology overlap (cyanorak vs kegg) doesn't inflate a single
score_A value. Only within-ontology overlap matters for scoring.
"""
from __future__ import annotations

import pickle
from itertools import combinations
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

ANALYSIS_DIR = Path(__file__).resolve().parent.parent
JACCARD_FLAG_THRESHOLD = 0.5


def main() -> None:
    signature = pd.read_csv(ANALYSIS_DIR / "data" / "reference_signature.csv")
    with open(ANALYSIS_DIR / "data" / "enrichment_results.pkl", "rb") as f:
        results = pickle.load(f)

    for ontology, sig_ont in signature.groupby("ontology"):
        # Find any EnrichmentResult for this ontology; term2gene is consistent
        # across background_used within the same ontology.
        ont_result = next(
            (r for (org, o, bg), r in results.items() if o == ontology),
            None,
        )
        if ont_result is None or len(sig_ont) < 2:
            continue
        t2g = ont_result.term2gene  # long DataFrame: locus_tag × term_id
        members = {
            tid: set(t2g.loc[t2g["term_id"] == tid, "locus_tag"])
            for tid in sig_ont["term_id"]
        }
        term_ids = list(members.keys())
        n = len(term_ids)

        # Pairwise Jaccard + subset detection.
        rows = []
        for a, b in combinations(term_ids, 2):
            A, B = members[a], members[b]
            if not A or not B:
                continue
            inter = A & B
            union = A | B
            jac = len(inter) / len(union) if union else 0.0
            subset = ("a_in_b" if A <= B else
                      "b_in_a" if B <= A else "none")
            flag = jac > JACCARD_FLAG_THRESHOLD or subset != "none"
            rows.append({
                "term_id_a": a, "n_a": len(A),
                "term_id_b": b, "n_b": len(B),
                "intersect": len(inter), "union": len(union),
                "jaccard": jac, "subset_relation": subset,
                "flag": flag,
            })
        pairs = pd.DataFrame(rows).sort_values("jaccard", ascending=False)
        pairs.to_csv(
            ANALYSIS_DIR / "exploration" / "qc"
            / f"step3_signature_redundancy_{ontology}.csv",
            index=False,
        )

        # Jaccard matrix heatmap.
        J = np.zeros((n, n))
        for i, a in enumerate(term_ids):
            for j, b in enumerate(term_ids):
                if i == j:
                    J[i, j] = 1.0
                    continue
                A, B = members[a], members[b]
                if not A or not B:
                    continue
                J[i, j] = len(A & B) / max(len(A | B), 1)
        fig, ax = plt.subplots(
            figsize=(max(6, 0.5 * n + 2), max(5, 0.5 * n + 1.5))
        )
        im = ax.imshow(J, cmap="OrRd", vmin=0, vmax=1, aspect="auto")
        labels = [t.split(":", 1)[-1] for t in term_ids]
        ax.set_xticks(range(n))
        ax.set_xticklabels(labels, rotation=90, fontsize=8)
        ax.set_yticks(range(n))
        ax.set_yticklabels(labels, fontsize=8)
        for i in range(n):
            for j in range(n):
                if i != j and J[i, j] > 0:
                    ax.text(
                        j, i, f"{J[i, j]:.2f}",
                        ha="center", va="center",
                        fontsize=6,
                        color="black" if J[i, j] < 0.7 else "white",
                    )
        fig.colorbar(im, ax=ax, label="Jaccard")
        ax.set_title(f"Step 3 signature pathway gene-overlap — {ontology}")
        fig.tight_layout()
        fig.savefig(
            ANALYSIS_DIR / "exploration" / "qc"
            / f"step3_signature_redundancy_{ontology}.png",
            dpi=140, bbox_inches="tight",
        )
        plt.close(fig)

        flagged = pairs[pairs["flag"]]
        print(f"\n=== {ontology}: {len(flagged)} flagged pair(s) (Jaccard>{JACCARD_FLAG_THRESHOLD} or subset) ===")
        if not flagged.empty:
            print(flagged[[
                "term_id_a", "n_a", "term_id_b", "n_b",
                "intersect", "jaccard", "subset_relation",
            ]].to_string(index=False))
        print(f"\n=== {ontology}: top-5 pairs by Jaccard (all) ===")
        print(pairs.head(5)[[
            "term_id_a", "n_a", "term_id_b", "n_b",
            "intersect", "jaccard", "subset_relation",
        ]].to_string(index=False))


if __name__ == "__main__":
    main()
