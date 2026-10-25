# dc: human genome reference data

Scrapes of public human-genome resources, kept for later analysis. Each folder under `scrape/` has its own README with sources, method and caveats.

| Folder | What it holds |
|---|---|
| `scrape/nhgri-gwas/` | NHGRI Talking Glossary entry for GWAS plus its four linked entries (text, HTML, audio, images) |
| `scrape/nhgri-1000-genomes/` | NHGRI archived 1000 Genomes Project page, tutorial pages and 11 tutorial slide decks (PDF and text) |
| `scrape/ensembl-human-genomes/` | Ensembl genome selector for human: 563 genome/assembly records and details for GRCh38.p14, GRCh37.p13 and T2T-CHM13v2.0 |
| `scrape/mavedb-human/` | MaveDB: 1,320 published score sets with a human target, full metadata, and a resumable downloader for the 1,871 score and count tables (data itself is git-ignored) |

The dating-thread engine that used to live here (dc-extended) moved to `ziadwalieddine-create/DC-extended-`; its history up to the split is preserved there and in this repository's earlier commits.
