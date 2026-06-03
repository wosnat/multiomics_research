# Step 3 — PaperBLAST run + candidate triage (pro pool)

## Context

Step 2 produced 32 seed genes (8 per pool) ranked for PaperBLAST hit-likelihood.
This step is the hybrid workflow's middle + back end: the **researcher runs
PaperBLAST offline** (the one non-scriptable link — PaperBLAST is Cloudflare-gated,
see `gaps_and_friction.md`) and saves one HTML per seed; **Claude parses, filters,
dedups, and triages** the saved pages into a candidate-publication shortlist.

This commit covers the **pro pool (8 seeds)** end-to-end, validating the whole
pipeline before the other three pools are run.

> **Update (pro + syn, 16 seeds).** The syn pool was added. Two pipeline changes:
> (a) `pb_io.py` shared loader now reads both `.html` and `.mhtml` (Single-File)
> saves and scans `results/<pool>/` subfolders recursively; (b) pro files moved
> into `results/prochlorococcus/`. Combined run: 1418 hits → 157 kept → 164
> distinct papers → 10 already in KG → **154 new candidates** (P1 in-Cyanorak = 19,
> KG-fit ≥2 = 11). The syn P1 set is more in-scope than pro's (WH8102 phosphorus
> proteomics, coastal *Synechococcus* comparisons, genome-wide DNA-damage
> transcriptional response) — less Synechocystis spillover. The per-pool numbers
> and "two axes pull apart" finding below are the original pro-only figures; the
> final combined write-up lands after alt + other_hetero are run.

> **Update (+ alt, 24 seeds).** Adding the alt pool exposed that the "Proteobacteria
> = same phylum" keep-rule was too loose for heterotroph seeds — it readmitted
> clinical/plant-pathogen Proteobacteria (P. aeruginosa, P. syringae) and ballooned
> the candidate list to 1073 (915 P3). Two fixes (`02_parse_hits.py`): the
> heterotroph pools now keep **KG genera only** (dropping the broad `other_proteo`
> tier) **and** apply an **identity floor** (`--min-identity-hetero`, default 50%)
> — a close homolog is more likely the actual marine relative, not a distant
> clinical paralog. This cut alt to 36 candidates and the total to **190** (P1=19,
> P2=135, P3=36; KG-fit ≥2 = 17). Residual alt noise is *Pseudomonas*-specific:
> P. putida is a KG genus but its literature is dominated by clinical/soil
> studies, so even high-identity hits surface off-domain papers. The genuinely
> marine alt candidates (~6) are Shewanella oneidensis + Ruegeria pomeroyi
> omics/genomics. Inventory also gained version-suffix tolerance
> (WP_014948722 ~ WP_014948722.1) since some saves drop the ".1". Two alt seeds
> (slyA, fecA) returned "No hits to characterized proteins" in PaperBLAST but
> still list uncharacterized-homolog papers.

## What I did

Three scripts, run from repo root with `env -u VIRTUAL_ENV uv run python <script>`:

1. **`01_inventory.py`** — map each saved `results/*.html` to its seed by the
   in-page query line (filenames are unreliable), cross-checked against
   `seeds_to_blast.csv`. Flags missing / duplicate / unmatched files.
2. **`02_parse_hits.py`** — parse each page into (hit → paper) rows. Extracts per
   homolog hit: organism, %identity, %coverage, and its linked papers (title,
   PMCID/PMID/DOI, year). Applies the **pool-aware relevance filter** and tags
   each hit `in_cyanorak`.
3. **`03_dedup_triage.py`** — collapse to distinct papers; resolve PMCID/PMID→DOI
   via the NCBI id-converter (API key from `.env`, on-disk cache, dedup-before-call
   so one id per paper); drop papers already in the KG's 43 publications; score
   KG-fit; assign priority tier; write split candidate files.

### Decisions baked into the scripts

- **Save format: "Webpage, HTML Only", one file per seed.** Compared against
  Complete / Single-File (MHTML); HTML-only is smallest and fully parseable
  (query line + hits + identities + PMCID/DOI links all present). The 3 formats
  carry identical text content.
- **Pool-aware relevance filter.** pro/syn seeds keep only **cyanobacterial**
  homolog hits; alt/other_hetero seeds keep **Proteobacteria** (the KG
  heterotrophs' phylum). Cross-lifestyle hits are noise for the seed's pool.
  All pools additionally drop MAG/SAG/metagenome sources (proper genomes only).
- **`in_cyanorak` flag.** A hit is "in Cyanorak" if its strain matches a strain
  code from the Cyanorak organism tables (`prochlorococcus.csv` /
  `synechococcus.csv` in `multiomics_biocypher_kg/data`). This separates "an
  actual KG marine-picocyano organism" from "some other cyanobacterium"
  (Synechocystis, Anabaena, freshwater elongatus — cyano but **not** Cyanorak).
- **KG-fit score (researcher-defined).** omics = 3; comparative-genomics /
  bioinformatics = 2; other = 1; single-protein characterization = 0. Broadened
  beyond wet-lab omics to anything matching the KG's gene/metabolite data model.
- **Priority tiers (lower = higher):** 1 = in-Cyanorak homolog; 2 = pro/syn pool,
  cyano but not Cyanorak; 3 = heterotroph pool (lower priority per researcher).
  Combined file sorted priority-first, then KG-fit, then seed-breadth, then
  identity; split into per-tier files.

## Results

### Coverage (pro pool)
8 of 8 pro seeds saved and matched. One missave caught and fixed by the inventory
(`PMT1742.html` was byte-identical to `PMT9312_0549.html` — rplN re-saved
correctly). The ntcA format-test page (not a seed) is in `results/_quarantine/`.

### Filter funnel (pro, 8 files)
| stage | count |
|---|---|
| raw homolog hits | 624 |
| kept (cyanobacterial, non-MAG) | 88 |
| dropped (plants, Gram-positives, gut/clinical Proteobacteria, MAGs) | 536 |
| paper rows from kept hits | 189 |
| distinct papers | 102 |
| already in KG (dropped) | 8 |
| **new candidates** | **94** |

### Candidate breakdown (94 new)
| priority | meaning | n |
|---|---|---|
| 1 | in-Cyanorak homolog | 12 |
| 2 | other cyanobacterium (pro/syn) | 82 |
| 3 | heterotroph | 0 (pools not yet run) |

| KG-fit kind | n |
|---|---|
| omics | 5 |
| comparative_genomics | 4 |
| other | 84 |
| characterization | 1 |

9 of 94 score KG-fit ≥2 (omics + comparative-genomics). 19 hit by >1 seed.
All 94 DOIs resolved; the 8 KG-matches verified correct (Tolonen N, salinity,
glucose proteomics, etc.).

### Key finding — the two ranking axes pull apart for pro
The in-scope papers (P1, in-Cyanorak) are mostly **not** ingestible datasets
(ncRNA studies, the MED4 genome paper, ProPortal, regulon inference); the genuine
omics/genomics datasets (KG-fit ≥2) are about **other** cyanobacteria
(Synechocystis, Crocosphaera, Anabaena, Leptolyngbya) surfaced because the Pro
seed genes are conserved in them. This matches the step-1 expectation that
Prochlorococcus — heavily studied and already well-represented in the KG — would
yield few genuinely new in-scope datasets. The heterotroph pools (alt,
other_hetero) are the likelier source of novel ingestible papers.

## Surprises

- **PaperBLAST links papers by PMCID, not DOI** (102 PMCIDs vs 3 DOIs in the
  scraped HTML). DOI resolution via NCBI is required for KG dedup; built with
  cache + backoff (NCBI throttles to 3 req/s without a key, 10 with).
- **"Same phylum" was too loose.** The first pass kept clinical/gut Proteobacteria
  (P. aeruginosa, E. coli, Salmonella) for cyano seeds. Fixed with the pool-aware
  filter (cyano-only for pro/syn).
- **Cyanorak ≠ "cyanobacteria".** Cyanorak is marine picocyano only — Synechocystis
  and freshwater Synechococcus elongatus are cyano but out of Cyanorak scope. The
  `in_cyanorak` flag uses the strain tables, not a genus match, to get this right.

## Advance rationale

The pro pool is processed end-to-end and the pipeline (inventory → parse+filter →
dedup+triage) is validated and reproducible. The candidate files are the
deliverable shape. Ready to run the other 3 pools through the same scripts (no
code changes expected — the pool-aware filter already handles heterotrophs) and
re-run for a combined 32-seed shortlist.

---

## Decide-gate checklist

- **Outputs produced.**
  - `scripts/01_inventory.py`, `02_parse_hits.py`, `03_dedup_triage.py`
  - `seeds_to_blast.csv` (32 seeds: pool, sel_rank, locus_tag, protein_id)
  - `results/*.html` (8 pro PaperBLAST saves; `_quarantine/` holds the ntcA format-test)
  - `data/01_inventory.csv`, `01_coverage.csv`, `02_hits.csv`, `02_papers.csv`, `02_dropped_genera.csv`, `kg_dois.txt`
  - `data/03_candidates.csv` (94) + `03_candidates_p1_cyanorak.csv` (12) / `03_candidates_p2_othercyano.csv` (82); `03_in_kg.csv` (8)
  - `03_idconv_cache.json` gitignored (reproducible from NCBI API)
  - Run order: `01_inventory.py` → `02_parse_hits.py` → `03_dedup_triage.py`.
- **Results presented.** Filter funnel, priority/KG-fit breakdowns, and the P1 + high-fit P2 lists shown inline above (same as in chat).
- **QC gate.**
  - Inventory caught + fixed the PMT1742 missave → 8/8 pro seeds covered.
  - Pool-aware filter verified to drop clinical/plant/Gram-positive hits for cyano seeds.
  - `in_cyanorak` verified against Cyanorak strain tables (excludes Synechocystis / freshwater elongatus).
  - KG dedup verified: 8 dropped DOIs cross-checked as genuinely in the KG.
  - DOI resolution: 102/102 PMCIDs resolved.
- **Decisions made this step.** HTML-only save format; pool-aware relevance filter; in_cyanorak flag from strain tables; KG-fit score (omics 3 / comp-genomics 2 / other 1 / characterization 0); priority tiers + split files; priority-first sort.
- **Advance rationale.** Pipeline validated on pro end-to-end; candidate files produced; ready for the remaining 3 pools.
