# MaveDB superseded score sets: data notice

`data_extensions/mavedb_superseded/` holds the score and count tables of 52 published MaveDB score sets that the public search omits because they are *superseded* versions (51 have a human or accession target). They are additive: the main collection on this branch (`scrape/mavedb-human/scores/`) is unchanged.

- Source: MaveDB Bulk Download v4 (Zenodo record 20840937, dump generated 2026-06-24, md5 `2cb818f9add61a378c30766d9b4b45ce`). The dump contains only CC0-licensed records, and all 52 of these sets are CC0.
- Files are the dump's `csv/<urn>.scores.csv` / `.counts.csv`, gzip-compressed without other change.
- Row counts equal `numVariants` from `GET /api/v1/score-sets/{urn}` for every table; SHA-256 values are in `acquisition/mavedb_superseded/superseded_manifest.jsonl` on `main`.
- These are older versions of experiments that have newer versions; the index on `main` (`acquisition/mavedb_superseded/superseded_index.csv`) names the superseding score set for each. Do not pool them with their successors.
- Restore: `git archive origin/data/mavedb-scores data_extensions/mavedb_superseded | tar -x`.
