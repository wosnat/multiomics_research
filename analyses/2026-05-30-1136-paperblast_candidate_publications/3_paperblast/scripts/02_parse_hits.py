#!/usr/bin/env python
"""
02_parse_hits.py — parse saved PaperBLAST result pages into (hit, paper) rows,
classify each hit's homolog organism by the relevance criteria, flag MAGs.

Relevance criteria — POOL-AWARE (locked 2026-05-30):
  pro / syn pools  → keep a hit only if the homolog is a cyanobacterium.
  alt / other_hetero pools → keep a hit if the homolog is in the KG heterotroph
        phylum (Proteobacteria) — both KG genera and other Proteobacteria.
  AND in all pools the homolog must come from a proper isolate genome, not a
  MAG/SAG/metagenome.
Rationale: heterotroph hits are noise for a cyano seed and vice-versa, so the
keep-set tracks the seed's pool. Everything else (plants, fungi, metazoa,
archaea, Firmicutes, and cross-lifestyle hits) is dropped but logged to
data/02_dropped_genera.csv for review.

Hit organism is tagged `kg_genus` (a genus already in the KG) vs `other_proteo` /
`cyano` so the keep-set can be tightened to KG-genera-only later if wanted.

Inputs:  results/*.html  (skips _quarantine/), seeds_to_blast.csv, inventory
Outputs (data/):
  02_hits.csv          one row per (seed, homolog hit) with organism/identity/keep flag
  02_papers.csv        one row per (seed, hit, paper) for KEPT hits — pre-dedup
  02_dropped_genera.csv  genera dropped by the filter, with counts (review)
  02_parse.log

Run:  env -u VIRTUAL_ENV uv run python 3_paperblast/scripts/02_parse_hits.py  (repo root)
"""
from __future__ import annotations

import html
import re
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from pb_io import load_pages  # noqa: E402

BASE = Path(__file__).resolve().parents[1]
RESULTS = BASE / "results"
DATA = BASE / "data"
SEEDS = BASE / "seeds_to_blast.csv"
CYANORAK_DIR = Path("/home/osnat/github/multiomics_biocypher_kg/data")
CYANORAK_CSVS = ["Cyanorak  Organism Table  prochlorococcus.csv",
                 "Cyanorak  Organism Table  synechococcus.csv"]


def load_cyanorak_tokens() -> set[str]:
    """Strain codes from the Cyanorak organism tables (Pro + marine Syn/Cyanobium).
    Used to tag whether a hit organism is IN the KG's marine-picocyano scope
    (vs merely 'some cyanobacterium' like Synechocystis / freshwater elongatus)."""
    import csv
    toks = set()
    for fn in CYANORAK_CSVS:
        p = CYANORAK_DIR / fn
        if not p.exists():
            continue
        for r in csv.DictReader(p.open(encoding="utf-8-sig")):
            nm = (r.get("Name") or "").strip()
            if "_" in nm:
                toks.add(nm.split("_", 1)[1])
    return {t for t in toks if t}


CYANORAK_TOKENS = load_cyanorak_tokens()


def in_cyanorak(organism: str) -> bool:
    """True if the hit organism's strain matches a Cyanorak strain code.
    Normalizes 'MIT 9313'->'MIT9313', 'WH 8102'->'WH8102' before matching."""
    norm = re.sub(r"(?<=[A-Za-z])\s+(?=\d)", "", organism)  # join 'MIT 9313'
    norm = norm.replace(" ", "")
    for tok in CYANORAK_TOKENS:
        if tok and tok.replace(" ", "") in norm:
            return True
    return False

# --- taxonomy keep-sets (curated; [interpretation], reviewable via dropped CSV) ---
CYANO = {
    "Prochlorococcus", "Synechococcus", "Parasynechococcus", "Picosynechococcus",
    "Synechocystis", "Cyanobium", "Nostoc", "Anabaena", "Cyanothece",
    "Thermosynechococcus", "Acaryochloris", "Microcystis", "Crocosphaera",
    "Leptolyngbya", "Trichormus", "Phormidium", "Gloeothece", "Halothece",
    "Dolichospermum", "Halomicronema", "Arthrospira", "Limnospira",
}
# Proteobacteria genera (the KG heterotrophs' phylum). KG_GENERA are those already
# represented in the KG; OTHER_PROTEO are phylum-mates seen in the hit set.
KG_GENERA = {"Alteromonas", "Marinobacter", "Pseudomonas", "Shewanella", "Ruegeria"}
OTHER_PROTEO = {
    "Escherichia", "Salmonella", "Vibrio", "Aliivibrio", "Citrobacter",
    "Proteus", "Klebsiella", "Enterobacter", "Yersinia", "Erwinia", "Pantoea",
    "Photobacterium", "Serratia", "Acinetobacter", "Pseudonocardia",
    "Rhodopseudomonas", "Bradyrhizobium", "Nitrobacter", "Rhizobium",
    "Sinorhizobium", "Agrobacterium", "Cupriavidus", "Ralstonia", "Burkholderia",
    "Polynucleobacter", "Paracidovorax", "Xanthomonas", "Xylella", "Neisseria",
    "Bordetella", "Helicobacter", "Campylobacter", "Sulfurospirillum",
    "Geobacter", "Desulfovibrio", "Thiobacillus", "Thioalkalivibrio",
    "Thiocapsa", "Allochromatium", "Halothiobacillus", "Hydrogenovibrio",
    "Thiomicrospira", "Acidithiobacillus", "Acidihalobacter", "Ferrovum",
    "Sideroxydans", "Mariprofundus", "Nitrosomonas", "Methylococcus",
    "Methylocaldum", "Thiobacter", "Legionella", "Francisella", "Haemophilus",
    "Pasteurella", "Actinobacillus", "Rickettsia", "Brucella", "Bartonella",
    "Dinoroseobacter", "Asaia", "Acetobacter", "Gluconobacter", "Rhodospirillum",
    "Xanthobacter", "Sodalis", "Baumannia", "Candidatus", "Thioglobus",
    "Gaiella", "Paracoccus", "Sphingomonas", "Caulobacter",
}
MAG_MARKERS = ("MAG", "metagenome", "uncultured", " SAG", "bin ", "Candidatus")


def flatten(t: str) -> str:
    return re.sub(r"\s+", " ", html.unescape(re.sub(r"<[^>]+>", " ", t)))


def query_id(flat: str, seeds: pd.DataFrame, fname: str) -> str | None:
    m = re.search(r"PaperBLAST Hits for (.{0,140})", flat)
    q = (m.group(1) if m else "") + " " + fname
    for lt in sorted(seeds["locus_tag"].astype(str), key=len, reverse=True):
        if lt in q:
            return lt
    for pid in seeds["protein_id"].astype(str):
        if pid in q or re.sub(r"\.\d+$", "", pid) in q:  # tolerate dropped ".1" version
            return seeds.loc[seeds["protein_id"] == pid, "locus_tag"].iloc[0]
    return None


def tier_of(genus: str) -> str:
    if genus in CYANO:
        return "cyano"
    if genus in KG_GENERA:
        return "kg_genus"
    if genus in OTHER_PROTEO:
        return "other_proteo"
    return "out"


def keep_for_pool(tier: str, pool: str, identity: int, min_identity_hetero: int) -> bool:
    """Pool-aware relevance filter:
      pro / syn          → cyanobacterial hits only (any identity).
      alt / other_hetero → KG heterotroph genera only AND identity >=
                           min_identity_hetero. The broad `other_proteo` tier is
                           dropped (it readmitted clinical/plant-pathogen
                           Proteobacteria like P. aeruginosa / P. syringae); the
                           identity floor further removes distant clinical
                           paralogs within the KG genera, keeping close marine
                           relatives. (locked 2026-06-03)"""
    if pool in ("pro", "syn"):
        return tier == "cyano"
    if pool in ("alt", "other_hetero"):
        return tier == "kg_genus" and identity >= min_identity_hetero
    return False


def parse_file(raw: str):
    """Yield hit dicts: organism, genus, identity, coverage, papers[]. `raw` is
    the page HTML (already decoded — mhtml handled upstream by pb_io)."""
    # split on each characterized-hit anchor: "...from <i>ORG</i><BR><a ...>NN% identity, NN% coverage</a><UL>...papers...</UL>"
    # iterate identity anchors; for each, look back for the nearest 'from <i>..</i>'
    hits = []
    for m in re.finditer(
        r"from <i>([^<]+)</i>.*?>(\d+)% identity,\s*(\d+)% coverage</a>(.*?)(?=from <i>|<h3|</body|$)",
        raw, flags=re.S,
    ):
        org = html.unescape(m.group(1)).strip()
        identity = int(m.group(2))
        coverage = int(m.group(3))
        tail = m.group(4)
        papers = []
        for pm in re.finditer(
            r'<a href="(http[^"]*(?:pmc/articles/PMC\d+|pubmed/\d+|doi\.org/[^"]+))"[^>]*>(.*?)</a>'
            r'<br\s*/?>\s*<small>(.*?)</small>', tail, flags=re.S,
        ):
            url = pm.group(1)
            title = html.unescape(re.sub(r"<[^>]+>", "", pm.group(2))).strip()
            meta = html.unescape(re.sub(r"<[^>]+>", " ", pm.group(3))).strip()
            doi = None
            md = re.search(r"doi\.org/([^\s\"<>]+)", url)
            if md:
                doi = md.group(1)
            pmcid = re.search(r"PMC(\d+)", url)
            pmid = re.search(r"pubmed/(\d+)", url)
            year = re.search(r"\b(19|20)\d{2}\b", meta)
            papers.append({
                "title": title, "meta": re.sub(r"\s+", " ", meta)[:120],
                "doi": doi, "pmcid": pmcid.group(0) if pmcid else None,
                "pmid": pmid.group(1) if pmid else None,
                "year": year.group(0) if year else None,
            })
        hits.append({"organism": org, "genus": org.split()[0].strip("[]"),
                     "identity": identity, "coverage": coverage, "papers": papers})
    return hits


def main() -> int:
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument("--min-identity-hetero", type=int, default=50,
                    help="min %%identity for alt/other_hetero KG-genus hits (drops distant clinical paralogs)")
    args = ap.parse_args()

    DATA.mkdir(exist_ok=True)
    seeds = pd.read_csv(SEEDS)
    logf = open(DATA / "02_parse.log", "w")
    log = lambda m: (print(m), logf.write(m + "\n"))
    log(f"== parse_hits: min_identity_hetero={args.min_identity_hetero} ==")

    hit_rows, paper_rows = [], []
    dropped = {}
    for f, raw in load_pages(RESULTS):
        flat = flatten(raw)
        seed = query_id(flat, seeds, f.name)
        if seed is None:
            log(f"  ? {f.name}: no seed match — skipped")
            continue
        pool = seeds.loc[seeds["locus_tag"] == seed, "pool"].iloc[0]
        hits = parse_file(raw)
        for h in hits:
            is_mag = any(mk.lower() in h["organism"].lower() for mk in MAG_MARKERS)
            tier = tier_of(h["genus"])
            keep_tax = keep_for_pool(tier, pool, h["identity"], args.min_identity_hetero)
            keep = keep_tax and not is_mag
            if not keep_tax:
                dropped[h["genus"]] = dropped.get(h["genus"], 0) + 1
            cyanorak = in_cyanorak(h["organism"])
            hit_rows.append({
                "seed": seed, "pool": pool, "organism": h["organism"],
                "genus": h["genus"], "identity": h["identity"],
                "coverage": h["coverage"], "tier": tier, "in_cyanorak": cyanorak,
                "is_mag": is_mag, "keep": keep, "n_papers": len(h["papers"]),
            })
            if keep:
                for p in h["papers"]:
                    paper_rows.append({
                        "seed": seed, "pool": pool, "homolog_organism": h["organism"],
                        "homolog_genus": h["genus"], "tier": tier,
                        "in_cyanorak": cyanorak,
                        "identity": h["identity"], "coverage": h["coverage"],
                        **p,
                    })
        log(f"  ✓ {f.name} [{seed}/{pool}]: {len(hits)} hits, "
            f"{sum(r['keep'] for r in hit_rows if r['seed']==seed)} kept")

    hits_df = pd.DataFrame(hit_rows)
    papers_df = pd.DataFrame(paper_rows)
    hits_df.to_csv(DATA / "02_hits.csv", index=False)
    papers_df.to_csv(DATA / "02_papers.csv", index=False)
    drop_df = (pd.DataFrame(sorted(dropped.items()), columns=["genus", "hit_count"])
               .sort_values("hit_count", ascending=False))
    drop_df.to_csv(DATA / "02_dropped_genera.csv", index=False)

    log(f"\n=== summary ===")
    log(f"files parsed: {hits_df['seed'].nunique()} seeds")
    log(f"total hits: {len(hits_df)} | kept: {int(hits_df['keep'].sum())} | "
        f"dropped (taxonomy): {int((~hits_df['keep']).sum())}")
    log(f"  of dropped: {int(hits_df['is_mag'].sum())} were MAG-flagged")
    log(f"kept-hit tier mix:\n{hits_df[hits_df['keep']].groupby('tier').size().to_string()}")
    log(f"\npaper rows (kept hits, pre-dedup): {len(papers_df)}")
    if len(papers_df):
        log(f"  unique DOIs: {papers_df['doi'].nunique()} | unique PMCIDs: {papers_df['pmcid'].nunique()}")
    log(f"\ntop dropped genera: {drop_df.head(12).to_dict('records')}")
    logf.close()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
