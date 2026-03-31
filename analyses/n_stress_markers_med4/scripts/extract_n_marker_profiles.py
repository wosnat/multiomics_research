#!/usr/bin/env python3
"""Extract differential expression data for candidate N-stress marker genes
across ALL MED4 experiments. Produces a per-gene × per-treatment summary.

Outputs:
  data/n_marker_de_all.csv       — full DE rows for all 21 genes across all experiments
  data/n_marker_specificity.csv  — per-gene summary: which treatments responded, direction, ranks
"""

import pandas as pd
from multiomics_explorer import differential_expression_by_gene, list_experiments

# --- Config ---
GENES = [
    "PMM0246",  # ntcA
    "PMM1463",  # glnB
    "PMM0393",  # pipX
    "PMM0920",  # glnA
    "PMM1512",  # glsF
    "PMM0263",  # amt1
    "PMM0370",  # cynA
    "PMM0371",  # cynB
    "PMM0372",  # cynD
    "PMM0373",  # cynS
    "PMM0970",  # urtA
    "PMM0971",  # urtB
    "PMM0972",  # urtC
    "PMM0973",  # urtD
    "PMM0974",  # urtE
    "PMM0963",  # ureC
    "PMM0964",  # ureB
    "PMM0965",  # ureA
    "PMM0958",  # hypothetical DUF1830
    "PMM1462",  # hypothetical (adj. glnB)
    "PMM0687",  # hypothetical
]

ORGANISM = "Prochlorococcus MED4"

# --- Step 1: Get all MED4 experiment IDs ---
print("Fetching experiment list...")
exp_data = list_experiments(organism="MED4", limit=100)
all_exp_ids = [e["experiment_id"] for e in exp_data["results"]]
print(f"  {len(all_exp_ids)} experiments total")

# Build experiment metadata lookup
exp_meta = {}
for e in exp_data["results"]:
    exp_meta[e["experiment_id"]] = {
        "treatment_type": e["treatment_type"],
        "omics_type": e["omics_type"],
        "publication_doi": e["publication_doi"],
        "table_scope": e["table_scope"],
        "is_time_course": e["is_time_course"],
    }

# --- Step 2: Extract DE data for all genes across all experiments ---
print("Extracting DE data for 21 genes across all experiments...")
de_data = differential_expression_by_gene(
    organism=ORGANISM,
    locus_tags=GENES,
    significant_only=False,
)
print(f"  {de_data['total_matching']} total rows returned")

df = pd.DataFrame(de_data["results"])
if len(df) == 0:
    print("ERROR: No DE data returned")
    exit(1)

# API already includes treatment_type, omics_type, publication_doi, table_scope

# Save full data
df.to_csv("data/n_marker_de_all.csv", index=False)
print(f"  Saved {len(df)} rows to data/n_marker_de_all.csv")

# --- Step 3: Build specificity summary ---
print("Building specificity summary...")

# For each gene × treatment_type: was gene tested? did it respond? direction?
summary_rows = []
for locus_tag in GENES:
    gene_df = df[df["locus_tag"] == locus_tag]
    gene_name = gene_df["gene_name"].iloc[0] if len(gene_df) > 0 and pd.notna(gene_df["gene_name"].iloc[0]) else ""

    for ttype in sorted(df["treatment_type"].unique()):
        tt_df = gene_df[gene_df["treatment_type"] == ttype]

        if len(tt_df) == 0:
            summary_rows.append({
                "locus_tag": locus_tag,
                "gene_name": gene_name,
                "gene_name_col": gene_name,
                "treatment_type": ttype,
                "status": "not_tested",
                "n_experiments": 0,
                "n_timepoints": 0,
                "n_sig_up": 0,
                "n_sig_down": 0,
                "max_log2fc": None,
                "min_log2fc": None,
                "median_log2fc": None,
                "best_rank_up": None,
                "best_rank_down": None,
            })
            continue

        n_exp = tt_df["experiment_id"].nunique()
        n_tp = len(tt_df)
        sig_df = tt_df[tt_df["expression_status"].isin(["significant_up", "significant_down"])]

        sig_up = sig_df[sig_df["log2fc"] > 0] if len(sig_df) > 0 else pd.DataFrame()
        sig_down = sig_df[sig_df["log2fc"] < 0] if len(sig_df) > 0 else pd.DataFrame()

        status = "not_responded"
        if len(sig_up) > 0 and len(sig_down) > 0:
            status = "mixed"
        elif len(sig_up) > 0:
            status = "up"
        elif len(sig_down) > 0:
            status = "down"

        summary_rows.append({
            "locus_tag": locus_tag,
            "gene_name": gene_name,
            "treatment_type": ttype,
            "status": status,
            "n_experiments": n_exp,
            "n_timepoints": n_tp,
            "n_sig_up": len(sig_up),
            "n_sig_down": len(sig_down),
            "max_log2fc": tt_df["log2fc"].max() if len(tt_df) > 0 else None,
            "min_log2fc": tt_df["log2fc"].min() if len(tt_df) > 0 else None,
            "median_log2fc": tt_df["log2fc"].median() if len(tt_df) > 0 else None,
            "best_rank_up": int(sig_up["rank"].min()) if len(sig_up) > 0 and "rank" in sig_up.columns and sig_up["rank"].notna().any() else None,
            "best_rank_down": int(sig_down["rank"].min()) if len(sig_down) > 0 and "rank" in sig_down.columns and sig_down["rank"].notna().any() else None,
        })

spec_df = pd.DataFrame(summary_rows)
spec_df.to_csv("data/n_marker_specificity.csv", index=False)
print(f"  Saved {len(spec_df)} rows to data/n_marker_specificity.csv")

# --- Step 4: Print summary table ---
print("\n=== N-STRESS MARKER SPECIFICITY ===\n")
treatment_types = sorted(df["treatment_type"].unique())
header = f"{'Gene':<12} {'Name':<6} " + " ".join(f"{t[:8]:>8}" for t in treatment_types) + "  Specificity"
print(header)
print("-" * len(header))

for locus_tag in GENES:
    gene_spec = spec_df[spec_df["locus_tag"] == locus_tag]
    gene_name = gene_spec["gene_name"].iloc[0] if len(gene_spec) > 0 else ""

    cells = []
    responded_types = []
    not_tested_types = []
    for ttype in treatment_types:
        row = gene_spec[gene_spec["treatment_type"] == ttype]
        if len(row) == 0 or row.iloc[0]["status"] == "not_tested":
            cells.append("   ?   ")
            not_tested_types.append(ttype)
        elif row.iloc[0]["status"] == "not_responded":
            cells.append("   -   ")
        elif row.iloc[0]["status"] == "up":
            cells.append(f"  UP({row.iloc[0]['n_sig_up']:>2})")
            responded_types.append(ttype)
        elif row.iloc[0]["status"] == "down":
            cells.append(f"  DN({row.iloc[0]['n_sig_down']:>2})")
            responded_types.append(ttype)
        elif row.iloc[0]["status"] == "mixed":
            cells.append(" MIX   ")
            responded_types.append(ttype)

    # Classify specificity
    n_only = all(t in ["nitrogen_stress", "coculture"] for t in responded_types)
    if len(responded_types) == 0:
        specificity = "none"
    elif n_only and len(not_tested_types) <= 2:
        specificity = "N-specific*"
    elif n_only:
        specificity = f"N-only(untested:{len(not_tested_types)})"
    elif len(responded_types) <= 2:
        specificity = "narrow"
    else:
        specificity = "broad"

    print(f"{locus_tag:<12} {str(gene_name):<6} " + " ".join(f"{c:>8}" for c in cells) + f"  {specificity}")

print("\nKey: UP(n)=sig up in n timepoints, DN(n)=sig down, MIX=both, -=tested not sig, ?=not tested")
print("N-specific*: responds only to N-stress/coculture AND tested in ≥6 treatment types")
print(f"N-only(untested:N): responds only to N/coculture but untested in N other types")
