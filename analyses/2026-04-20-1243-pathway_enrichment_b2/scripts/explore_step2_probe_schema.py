"""Inspect EnrichmentResult API to find the right accessors."""
from __future__ import annotations

import pickle
from pathlib import Path

ANALYSIS_DIR = Path(__file__).resolve().parent.parent

with open(ANALYSIS_DIR / "data" / "enrichment_results.pkl", "rb") as f:
    results = pickle.load(f)

r = results[("Prochlorococcus MED4", "cyanorak_role", "table_scope")]

print("=== type ===")
print(type(r))
print("\n=== public attrs ===")
for a in dir(r):
    if not a.startswith("_"):
        print(f"  {a}")

print("\n=== inputs attrs ===")
for a in dir(r.inputs):
    if not a.startswith("_"):
        v = getattr(r.inputs, a)
        desc = f" — type={type(v).__name__}" if not callable(v) else " (callable)"
        print(f"  {a}{desc}")

print("\n=== sample clusters ===")
clusters = sorted(r.results["cluster"].unique())
for c in clusters[:10]:
    print(f"  {c}")
print(f"  ... ({len(clusters)} total)")

print("\n=== any Tolonen clusters ===")
for c in clusters:
    if "msb4100087" in c:
        print(f"  {c}")
