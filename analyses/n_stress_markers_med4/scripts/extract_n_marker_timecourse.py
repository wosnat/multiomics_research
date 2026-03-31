#!/usr/bin/env python3
"""Extract temporal profiles for Tier 1-2 N-stress marker genes
across all nitrogen_stress time course experiments.

Outputs:
  data/n_marker_timecourse.csv — per-gene, per-experiment, per-timepoint DE data
"""

import pandas as pd
from multiomics_explorer import differential_expression_by_gene

# Tier 1: N-specific
# Tier 2: N+coculture (biologically coherent)
TIER1 = ["PMM0965", "PMM0964", "PMM1463", "PMM1512", "PMM1462"]  # ureA, ureB, glnB, glsF, hyp
TIER2 = ["PMM0246", "PMM0371", "PMM0372", "PMM0393", "PMM0687", "PMM0973"]  # ntcA, cynB, cynD, pipX, hyp, urtD
ALL_MARKERS = TIER1 + TIER2

# N-stress time course experiment IDs
N_TIMECOURSE_EXPS = [
    "10.1038/msb4100087_nitrogen_stress_nitrogen_deprivation_med4_med4_microarray",  # Tolonen 0-48h
    "10.1038/ismej.2017.88_nitrogen_stress_ndepleted_pro99_medium_med4_rnaseq",      # Read 3-24h
    "10.1101/2025.11.24.690089_growth_state_pro99lown_nutrient_starvation_med4_rnaseq_axenic",     # Weissberg axenic RNA-seq
    "10.1101/2025.11.24.690089_growth_state_pro99lown_nutrient_starvation_med4_rnaseq_coculture",  # Weissberg coculture RNA-seq
    "10.1101/2025.11.24.690089_growth_state_pro99lown_nutrient_starvation_med4_proteomics_axenic", # Weissberg axenic proteomics
    "10.1101/2025.11.24.690089_growth_state_pro99lown_nutrient_starvation_med4_proteomics_coculture", # Weissberg coculture proteomics
]

# Also get coculture experiments for mirror comparison
COCULTURE_EXPS = [
    "10.1038/ismej.2016.70_coculture_alteromonas_hot1a3_med4_rnaseq",
    "10.1101/2025.11.24.690089_coculture_alteromonas_hot1a3_med4_rnaseq",
]

ALL_EXPS = N_TIMECOURSE_EXPS + COCULTURE_EXPS

print("Extracting DE data for 11 Tier 1-2 markers across N-stress + coculture experiments...")
de_data = differential_expression_by_gene(
    organism="Prochlorococcus MED4",
    locus_tags=ALL_MARKERS,
    experiment_ids=ALL_EXPS,
    significant_only=False,
)
print(f"  {de_data['total_matching']} rows returned")

df = pd.DataFrame(de_data["results"])

# Add tier label
df["tier"] = df["locus_tag"].apply(lambda x: "Tier1_N-specific" if x in TIER1 else "Tier2_N+coculture")

# Save
df.to_csv("data/n_marker_timecourse.csv", index=False)
print(f"  Saved {len(df)} rows to data/n_marker_timecourse.csv")

# Print summary: gene × timepoint for Tolonen time course
print("\n=== TOLONEN TIME COURSE (microarray, 0-48h) ===\n")
tolonen = df[df["experiment_id"].str.contains("msb4100087") & df["experiment_id"].str.contains("nitrogen")]
if len(tolonen) > 0:
    pivot = tolonen.pivot_table(
        index=["tier", "locus_tag", "gene_name"],
        columns="timepoint",
        values="log2fc",
        aggfunc="first"
    )
    # Reorder columns by timepoint_hours
    tp_order = tolonen.drop_duplicates("timepoint").sort_values("timepoint_hours")["timepoint"].tolist()
    pivot = pivot[[c for c in tp_order if c in pivot.columns]]
    print(pivot.to_string())

# Print coculture comparison
print("\n=== COCULTURE vs AXENIC (Weissberg RNA-seq) ===\n")
cocult = df[df["treatment_type"] == "coculture"]
if len(cocult) > 0:
    for _, row in cocult.iterrows():
        sig = "***" if row.get("expression_status", "") in ["significant_up", "significant_down"] else ""
        print(f"  {row['locus_tag']:10s} {str(row.get('gene_name','')):6s} log2FC={row['log2fc']:+.2f} {sig}  ({row['experiment_id'][:30]}...)")

# Print Weissberg proteomics coculture time course
print("\n=== WEISSBERG PROTEOMICS COCULTURE TIME COURSE ===\n")
prot_cocult = df[df["experiment_id"].str.contains("proteomics_coculture")]
if len(prot_cocult) > 0:
    pivot = prot_cocult.pivot_table(
        index=["tier", "locus_tag", "gene_name"],
        columns="timepoint",
        values="log2fc",
        aggfunc="first"
    )
    tp_order = prot_cocult.drop_duplicates("timepoint").sort_values("timepoint_hours")["timepoint"].tolist()
    pivot = pivot[[c for c in tp_order if c in pivot.columns]]
    print(pivot.to_string())
