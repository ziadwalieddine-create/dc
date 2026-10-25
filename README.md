# dc: human genome reference data

Scrapes of public human-genome resources, kept for later analysis. Each folder under `scrape/` has its own README with sources, method and caveats.

| Folder | What it holds |
|---|---|
| `scrape/nhgri-gwas/` | NHGRI Talking Glossary entry for GWAS plus its four linked entries (text, HTML, audio, images) |
| `scrape/nhgri-1000-genomes/` | NHGRI archived 1000 Genomes Project page, tutorial pages and 11 tutorial slide decks (PDF and text) |
| `scrape/ensembl-human-genomes/` | Ensembl genome selector for human: 563 genome/assembly records and details for GRCh38.p14, GRCh37.p13 and T2T-CHM13v2.0 |
| `scrape/mavedb-human/` | MaveDB: 1,402 published score sets with a human target (1,320 by organism filter plus 82 defined by GRCh38 accessions), full metadata, and a resumable downloader and verifier for the 1,980 score and count tables (the tables are stored on the `data/mavedb-scores` branch, git-ignored on `main`) |

The dating-thread engine that used to live here (dc-extended) moved to `ziadwalieddine-create/DC-extended-`; its history up to the split is preserved there and in this repository's earlier commits.

## Reference-library extension

The collections above are the foundation; the work below extends them without altering them (`python3 inventory/audit_existing.py --check` proves `scrape/` is unchanged apart from one logged correction).

| Path | What it holds |
|---|---|
| `catalogue/` | Evidence-backed catalogue of external human-genomics resources: one record per source, generated `MASTER_INVENTORY.md`, `COVERAGE_MATRIX.md`, `GAP_MATRIX.md`, `REUSE_PLAN.md`, `ACQUIRED_FILES.md`, `ATTRIBUTION.md` |
| `acquisition/` | Reproducible downloader with source-checksum verification, per-protein and remote-tabix (gnomAD) retrieval, validators, compatibility (join) tests with unit tests for the helpers, and an additive set of 52 superseded MaveDB score sets |
| `inventory/` | Audit of the existing assets, a checksum guard for `scrape/`, and a repository policy check (`check_repo_policy.py`: file sizes, personal-genome file types, credential patterns, live visibility) |
| `AUDIT_LOG.md` | What was added, changed, corrected and left unresolved |

Downloaded external data live in `data/external/` and are git-ignored. This repository is public (confirmed live by `inventory/check_repo_policy.py`); licences differ per source (see the catalogue) and restricted or personal material must not be committed.
