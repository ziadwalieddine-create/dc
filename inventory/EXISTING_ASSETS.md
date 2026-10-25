# Existing assets: inventory and verified status

Generated 2026-10-25T09:58:08Z at commit `ddafb14d` by `inventory/audit_existing.py` (regenerate rather than edit).

Repository visibility at audit time: `{"private": false, "visibility": "public"}`. While the repository is public, only openly redistributable data and metadata may be committed.

## Branches

| Branch | Tip | Files | Role |
|---|---|---|---|
| `origin/ccr-b13fd75a-ra9gfy` | `7aef816` | 98 | working/PR branch |
| `origin/data/mavedb-scores` | `7cfbd17` | 1981 | orphan data branch: the 1,980 MaveDB score/count tables + licence notice (not to be merged) |
| `origin/land/dc-extended-pack` | `3547d0a` | 90 | legacy branch from the original dating-project repository; unrelated, left untouched |
| `origin/main` | `ddafb14` | 98 | default branch; holds scrape/ (metadata, scripts, manifests) and the new catalogue |

## Assets

| Asset | Path | What it is | Files | Size | Status |
|---|---|---|---|---|---|
| `nhgri-gwas` | `scrape/nhgri-gwas` | glossary entry + 4 linked entries (HTML, text, audio, images) | 43 | 10.7 MB | **downloaded_verified** |
| `nhgri-1000-genomes` | `scrape/nhgri-1000-genomes` | archived project page, 2 tutorial pages, 11 slide decks (PDF + text), 1 staff profile | 33 | 36.2 MB | **downloaded_verified** |
| `ensembl-human-genomes` | `scrape/ensembl-human-genomes` | genome/assembly registry for human (563 records) + details for GRCh38.p14, GRCh37.p13, T2T-CHM13v2.0 | 9 | 0.7 MB | **downloaded_verified** |
| `mavedb-human` | `scrape/mavedb-human` | 1,402 published human-target score sets: metadata + 1,980 score/count tables (tables on branch data/mavedb-scores) | 12 | 37.7 MB | **downloaded_verified** |

## Verification detail

### nhgri-gwas
- Source: https://www.genome.gov/genetics-glossary/Genome-Wide-Association-Studies
- Evidence value: background/terminology only (no variant-level data)
- Checks: `{"json_valid": true, "related_entries": 4, "media_assets": 31, "media_magic_failures": []}`

### nhgri-1000-genomes
- Source: https://www.genome.gov/27528684/1000-genomes-project
- Evidence value: documentation of the 1000 Genomes project (no genotype data)
- Checks: `{"json_valid": true, "tutorial_talks": 11, "slide_files_missing": [], "media_magic_failures": []}`

### ensembl-human-genomes
- Source: https://www.ensembl.org/genome-selector/search?species_taxonomy_id=9606
- Evidence value: assembly/annotation registry; no sequences or annotations downloaded
- Checks: `{"records": 563, "total_hits": 563, "csv_rows": 563, "consistent": true}`

### mavedb-human
- Source: https://www.mavedb.org/search (API https://api.mavedb.org/api/v1)
- Evidence value: direct experimental variant effects (deep mutational scanning and related MAVE assays); mapped to protein/nucleotide variants per score set
- Checks: `{"verify_scores_exit": 0, "verify_scores_output": ["score sets 1402 | files 1980 | score rows 5170093 | count rows 1544328", "known tmp: accession rows: 19100", "PROBLEMS: 0"], "scores_present_locally": true}`

## Preservation guard

`inventory/baseline_checksums.sha256` holds the sha256 of every git-tracked file under `scrape/` at the baseline commit. Run `python3 inventory/audit_existing.py --check` at any time: it fails if any existing file was changed or removed. New work is added outside `scrape/` (in `catalogue/`, `inventory/`, `acquisition/`).
