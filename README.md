# dc: human genome reference data

Scrapes of public human-genome resources, kept for later analysis. Each folder under `scrape/` has its own README with sources, method and caveats.

| Folder | What it holds |
|---|---|
| `scrape/nhgri-gwas/` | NHGRI Talking Glossary entry for GWAS plus its four linked entries (text, HTML, audio, images) |
| `scrape/nhgri-1000-genomes/` | NHGRI archived 1000 Genomes Project page, tutorial pages and 11 tutorial slide decks (PDF and text) |
| `scrape/ensembl-human-genomes/` | Ensembl genome selector for human: 563 genome/assembly records and details for GRCh38.p14, GRCh37.p13 and T2T-CHM13v2.0 |
| `scrape/mavedb-human/` | MaveDB: 1,402 published score sets with a human target (1,320 by organism filter plus 82 defined by GRCh38 accessions), full metadata, and a resumable downloader and verifier for the 1,980 score and count tables (the tables are stored on the `data/mavedb-scores` branch, git-ignored on `main`) |

The dating-thread engine that used to live here (dc-extended) moved to `ziadwalieddine-create/DC-extended-`; its history up to the split is preserved there and in this repository's earlier commits.
