#!/usr/bin/env python
"""
03_dedup_triage.py — turn parsed kept-hit paper rows into the candidate shortlist.

Pipeline:
  1. collapse 02_papers.csv to distinct papers (key: pmcid > pmid > doi > title).
  2. resolve PMCID/PMID -> DOI via NCBI id-converter (batched, cached to
     data/03_idconv_cache.json) so we can dedup against the KG.
  3. drop papers already in the KG (data/kg_dois.txt, 43 DOIs).
  4. classify each survivor: omics-dataset vs characterization (title/journal
     keywords) — a heuristic tag, not a final call.
  5. score by: #distinct seeds it was hit by (breadth) and best identity.
  6. write the ranked candidate shortlist.

Outputs (data/):
  03_candidates.csv       ranked distinct NEW papers (the shortlist)
  03_in_kg.csv            papers dropped because already in the KG
  03_idconv_cache.json    PMCID/PMID -> {doi,pmid,...} cache
  03_dedup_triage.log

Run:  env -u VIRTUAL_ENV uv run python 3_paperblast/scripts/03_dedup_triage.py  (repo root)
"""
from __future__ import annotations

import json
import os
import re
import time
import urllib.parse
import urllib.request
from pathlib import Path

import pandas as pd
from dotenv import load_dotenv

BASE = Path(__file__).resolve().parents[1]
DATA = BASE / "data"
CACHE = DATA / "03_idconv_cache.json"
EMAIL = "osnat.weissberg@gmail.com"
TOOL = "multiomics_kg"

# NCBI_API_KEY from the repo-root .env (raises NCBI limit 3->10 req/s).
load_dotenv(BASE.parents[1] / ".env")
NCBI_KEY = os.environ.get("NCBI_API_KEY")

# KG-fit = produces gene- or metabolite-level data the KG can ingest. Per the
# researcher: omics PLUS comparative-genomics / bioinformatics count, not just
# wet-lab omics. Three positive families + a negative (single-protein) family.
OMICS_KW = re.compile(
    r"\b(transcriptom|proteom|metabolom|lipidom|exoproteom|metaproteom|"
    r"RNA-?seq|mRNA-?seq|microarray|deep sequencing|expression profil|"
    r"gene expression|global (gene|transcription)|multi-?omic|"
    r"differential expression|expression dynamics|expression pattern)\b", re.I)
GENOMICS_KW = re.compile(
    r"\b(comparative genom|pan-?genom|genome-?wide|genomic analysis|"
    r"bioinformatic|in silico|phylogenom|genome sequenc|whole genome|"
    r"reconstruct|regulon|metabolic (model|network|reconstruction)|"
    r"genome annotation)\b", re.I)
CHARACTERIZATION_KW = re.compile(
    r"\b(crystal structure|3d structure|mutant|knockout|deletion strain|"
    r"kinetics|enzymatic|biochemical char|purif|catalytic|reaction mechanism|"
    r"substrate specificity|crystalliz)\b", re.I)


def kg_fit(title: str, meta: str):
    """Return (kind, score 0..3). Higher = better fit to the KG's gene/metabolite
    data model. omics=3, comparative-genomics/bioinformatics=2, mixed/ambiguous=1,
    pure single-protein characterization or other=0."""
    s = f"{title} {meta}"
    omics = bool(OMICS_KW.search(s))
    genomics = bool(GENOMICS_KW.search(s))
    char = bool(CHARACTERIZATION_KW.search(s))
    if omics:
        return ("omics", 3)
    if genomics:
        return ("comparative_genomics", 2)
    if char:
        return ("characterization", 0)
    return ("other", 1)


def idconv(ids: list[str], cache: dict, log, batch_size: int = 100) -> dict:
    """Resolve PMCIDs/PMIDs to records via NCBI id-converter, with on-disk cache.

    Uses NCBI_API_KEY when available (limit 10 req/s, else 3). Backoff-retry on
    HTTP 429. Cache written after every batch so a mid-run failure never loses
    progress."""
    spacing = 0.15 if NCBI_KEY else 1.0
    todo = [i for i in ids if i and i not in cache]
    for k in range(0, len(todo), batch_size):
        batch = todo[k:k + batch_size]
        params = {"ids": ",".join(batch), "format": "json", "tool": TOOL, "email": EMAIL}
        if NCBI_KEY:
            params["api_key"] = NCBI_KEY
        url = f"https://www.ncbi.nlm.nih.gov/pmc/utils/idconv/v1.0/?{urllib.parse.urlencode(params)}"
        for attempt in range(5):
            try:
                req = urllib.request.Request(url, headers={"User-Agent": f"{TOOL}/1.0"})
                with urllib.request.urlopen(req, timeout=30) as r:
                    d = json.load(r)
                for rec in d.get("records", []):
                    key = rec.get("pmcid") or rec.get("pmid")
                    if key:
                        cache[key] = {"doi": (rec.get("doi") or "").lower() or None,
                                      "pmid": rec.get("pmid"), "pmcid": rec.get("pmcid")}
                log(f"  idconv batch {k//batch_size+1}: +{len(batch)} ids")
                break
            except Exception as e:
                wait = 2 ** attempt
                log(f"  idconv batch {k//batch_size+1} attempt {attempt+1} failed ({e!r}); retry in {wait}s")
                time.sleep(wait)
        CACHE.write_text(json.dumps(cache))
        time.sleep(spacing)
    return cache


def main() -> int:
    log_f = open(DATA / "03_dedup_triage.log", "w")
    log = lambda m: (print(m), log_f.write(m + "\n"))

    papers = pd.read_csv(DATA / "02_papers.csv")
    log(f"input paper rows (kept hits): {len(papers)}")

    # distinct-paper key
    def pkey(r):
        return r.get("pmcid") or (f"PMID{r['pmid']}" if pd.notna(r.get("pmid")) else None) \
            or r.get("doi") or ("T:" + str(r["title"]).lower()[:60])
    papers["pkey"] = papers.apply(pkey, axis=1)

    # aggregate to one row per distinct paper, tracking which seeds/pools hit it
    agg = papers.groupby("pkey").agg(
        title=("title", "first"),
        meta=("meta", "first"),
        year=("year", "first"),
        doi_raw=("doi", "first"),
        pmcid=("pmcid", "first"),
        pmid=("pmid", "first"),
        n_seeds=("seed", "nunique"),
        seeds=("seed", lambda s: ";".join(sorted(set(s)))),
        pools=("pool", lambda s: ";".join(sorted(set(s)))),
        best_identity=("identity", "max"),
        in_cyanorak=("in_cyanorak", "max"),   # any hit in Cyanorak scope
        homolog_genera=("homolog_genus", lambda s: ";".join(sorted(set(s)))),
    ).reset_index()
    log(f"distinct papers: {len(agg)}")

    # resolve ids -> doi. Dedup BEFORE calling NCBI: only DOI-less distinct
    # papers need resolution, and each contributes exactly ONE id (prefer PMCID,
    # else PMID). This avoids sending both ids per paper and re-resolving papers
    # whose DOI we already scraped from the HTML.
    def needs_id(r):
        if pd.notna(r["doi_raw"]) and r["doi_raw"]:
            return None
        if pd.notna(r["pmcid"]) and r["pmcid"]:
            return r["pmcid"]
        if pd.notna(r["pmid"]):
            return f"{int(r['pmid'])}"
        return None
    want = sorted({x for x in agg.apply(needs_id, axis=1) if x})
    log(f"distinct papers needing DOI resolution: {len(want)} "
        f"(of {len(agg)}; {int(agg['doi_raw'].notna().sum())} already had a DOI)")
    cache = json.loads(CACHE.read_text()) if CACHE.exists() else {}
    cache = idconv(want, cache, log)

    def resolve_doi(r):
        if pd.notna(r["doi_raw"]) and r["doi_raw"]:
            return str(r["doi_raw"]).lower()
        for key in (r["pmcid"], f"{int(r['pmid'])}" if pd.notna(r["pmid"]) else None):
            if key and key in cache and cache[key].get("doi"):
                return cache[key]["doi"]
        return None
    agg["doi"] = agg.apply(resolve_doi, axis=1)

    # dedup vs KG
    kg = set(Path(DATA / "kg_dois.txt").read_text().split())
    agg["in_kg"] = agg["doi"].apply(lambda d: bool(d) and d in kg)
    in_kg = agg[agg["in_kg"]].copy()
    new = agg[~agg["in_kg"]].copy()
    log(f"already in KG: {len(in_kg)} | new candidates: {len(new)}")

    # triage: KG-fit kind + score
    kf = new.apply(lambda r: kg_fit(str(r["title"]), str(r["meta"])), axis=1)
    new["kind"] = [k for k, _ in kf]
    new["kg_fit_score"] = [s for _, s in kf]

    # priority tier (lower = higher priority):
    #   1 = in-Cyanorak homolog (directly about a KG marine-picocyano organism)
    #   2 = pro/syn pool, cyano homolog but NOT in Cyanorak (other cyanobacterium)
    #   3 = alt/other_hetero pool (heterotroph — lower priority per researcher)
    def priority(r):
        if r["in_cyanorak"]:
            return 1
        pools = set(str(r["pools"]).split(";"))
        if pools & {"pro", "syn"}:
            return 2
        return 3
    new["priority"] = new.apply(priority, axis=1)

    # rank within: priority, then KG-fit, then seed-breadth, then identity
    new = new.sort_values(
        ["priority", "kg_fit_score", "n_seeds", "best_identity"],
        ascending=[True, False, False, False])

    cols = ["priority", "in_cyanorak", "kind", "kg_fit_score", "n_seeds", "seeds",
            "pools", "best_identity", "year", "title", "doi", "pmcid", "pmid",
            "homolog_genera", "meta"]
    new[cols].to_csv(DATA / "03_candidates.csv", index=False)
    # split files by priority so the high-value set is separate
    new[new["priority"] == 1][cols].to_csv(DATA / "03_candidates_p1_cyanorak.csv", index=False)
    new[new["priority"] == 2][cols].to_csv(DATA / "03_candidates_p2_othercyano.csv", index=False)
    new[new["priority"] == 3][cols].to_csv(DATA / "03_candidates_p3_heterotroph.csv", index=False)
    in_kg[["title", "doi", "year", "n_seeds", "seeds", "in_cyanorak"]].to_csv(
        DATA / "03_in_kg.csv", index=False)

    log(f"\n=== candidate shortlist: {len(new)} new papers ===")
    log(f"by priority tier:\n{new['priority'].value_counts().sort_index().to_string()}")
    log(f"  p1 = in-Cyanorak | p2 = other-cyano (pro/syn) | p3 = heterotroph")
    log(f"by kind:\n{new['kind'].value_counts().to_string()}")
    log(f"in_cyanorak: {int(new['in_cyanorak'].sum())} | "
        f"omics+comparative-genomics (score>=2): {int((new['kg_fit_score']>=2).sum())}")
    log(f"multi-seed (hit by >1 seed): {(new['n_seeds']>1).sum()}")
    log(f"resolved DOI: {new['doi'].notna().sum()}/{len(new)}")
    log_f.close()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
