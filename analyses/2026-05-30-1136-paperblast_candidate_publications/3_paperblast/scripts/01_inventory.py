#!/usr/bin/env python
"""
01_inventory.py — map every saved PaperBLAST .html to one of the 32 seeds.

Filenames are unreliable (PaperBLAST(N).html, typos, PB_ prefixes), so each file
is identified by its in-page query line ("PaperBLAST Hits for ...") cross-checked
against the seed table. The query line carries, variously: the locus_tag, the
RefSeq protein_id, a UniProt accession, and/or the gene name (GN=) + product.

Matching strategy (first hit wins):
  1. any seed locus_tag appears in query text or filename
  2. any seed protein_id appears in query text or filename
  3. GN=<gene> from query matches a seed gene_name (only if unambiguous)

Outputs (../data/):
  01_inventory.csv   file -> matched seed (pool, sel_rank, locus_tag) + query + n_hits
  01_coverage.csv    per-seed: covered? by which file(s); duplicates flagged

Run:  env -u VIRTUAL_ENV uv run python 3_paperblast/scripts/01_inventory.py
(from repo root)
"""
from __future__ import annotations

import html
import re
from pathlib import Path

import pandas as pd

BASE = Path(__file__).resolve().parents[1]          # 3_paperblast/
RESULTS = BASE / "results"
DATA = BASE / "data"
SEEDS = BASE / "seeds_to_blast.csv"


def flatten(t: str) -> str:
    return re.sub(r"\s+", " ", html.unescape(re.sub(r"<[^>]+>", " ", t)))


def query_line(flat: str) -> str:
    m = re.search(r"PaperBLAST Hits for (.*?)\(\d+ a\.a", flat)
    if m:
        return m.group(1).strip()
    m = re.search(r"PaperBLAST Hits for (.{0,120})", flat)
    return m.group(1).strip() if m else ""


def n_hits(flat: str) -> int | None:
    m = re.search(r"Found (\d+) similar", flat)
    return int(m.group(1)) if m else None


def main() -> int:
    DATA.mkdir(exist_ok=True)
    seeds = pd.read_csv(SEEDS)
    loci = {str(r.locus_tag): r for r in seeds.itertuples()}
    pids = {str(r.protein_id): r for r in seeds.itertuples()}
    # gene_name map (may be absent); only used as last resort, must be unique
    seeds_gn = pd.read_csv(BASE.parent / "2_kg_selection" / "data" / "05_selected_papers.csv")
    gn_counts = seeds_gn["rep_gene_name"].value_counts()
    gn_unique = {g: lt for g, lt in zip(seeds_gn["rep_gene_name"], seeds_gn["rep_locus_tag"])
                 if pd.notna(g) and gn_counts[g] == 1}

    rows = []
    for f in sorted(RESULTS.glob("*.html")):
        t = f.read_text(encoding="utf-8", errors="replace")
        flat = flatten(t)
        q = query_line(flat)
        hay = q + " " + f.name  # search query text + filename
        match, how = None, None
        # locus_tag as plain substring (distinctive enough; handles PB-PMM1028, etc.)
        for lt, r in sorted(loci.items(), key=lambda kv: -len(kv[0])):
            if lt in hay:
                match, how = r, "locus_tag"
                break
        if match is None:
            for pid, r in pids.items():
                if pid in hay:
                    match, how = r, "protein_id"
                    break
        if match is None:
            gn = re.search(r"GN=([A-Za-z0-9_]+)", q)
            if gn and gn.group(1) in gn_unique:
                lt = gn_unique[gn.group(1)]
                match, how = loci[lt], "gene_name"
        rows.append({
            "file": f.name,
            "matched_pool": getattr(match, "pool", None),
            "matched_sel_rank": getattr(match, "sel_rank", None),
            "matched_locus_tag": getattr(match, "locus_tag", None),
            "matched_by": how,
            "n_hits": n_hits(flat),
            "query": q[:110],
        })
    inv = pd.DataFrame(rows)
    inv.to_csv(DATA / "01_inventory.csv", index=False)

    # coverage per seed
    cov_rows = []
    for r in seeds.itertuples():
        hits = inv[inv["matched_locus_tag"] == r.locus_tag]
        cov_rows.append({
            "pool": r.pool, "sel_rank": r.sel_rank, "locus_tag": r.locus_tag,
            "protein_id": r.protein_id,
            "n_files": len(hits),
            "files": "; ".join(hits["file"].tolist()),
            "status": ("MISSING" if len(hits) == 0 else
                       "ok" if len(hits) == 1 else "DUPLICATE"),
        })
    cov = pd.DataFrame(cov_rows).sort_values(["pool", "sel_rank"])
    cov.to_csv(DATA / "01_coverage.csv", index=False)

    unmatched = inv[inv["matched_locus_tag"].isna()]
    print(f"files scanned: {len(inv)}")
    print(f"matched to a seed: {inv['matched_locus_tag'].notna().sum()}")
    print(f"unmatched files (junk/format-test/other): {len(unmatched)}")
    print(f"seeds covered: {(cov['status']!='MISSING').sum()} / {len(cov)}")
    print(f"  MISSING: {cov[cov['status']=='MISSING']['locus_tag'].tolist()}")
    print(f"  DUPLICATE: {cov[cov['status']=='DUPLICATE']['locus_tag'].tolist()}")
    print(f"\nunmatched files: {unmatched['file'].tolist()}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
