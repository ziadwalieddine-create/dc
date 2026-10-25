# MaveDB: human score sets

Source page: https://www.mavedb.org/search (fetched 2026-10-25; MaveDB API v2026.3.2)

The page is a JavaScript app (empty HTML shell), so this was scraped from the API it uses (https://api.mavedb.org/api/v1, OpenAPI spec at `/openapi.json`). "Just the human genome" is applied as two selections of **published** score sets, 1,402 in total:

| Selection | How | Sets |
|---|---|---|
| `organism_filter` | The page's own filter: `POST /score-sets/search` with `{"published": true, "targetOrganismNames": ["Homo sapiens"]}` (NCBI taxonomy 9606) | 1,320 |
| `accession_target` | Human genes whose target is a GRCh38/transcript accession (`NM_`, `ENST`, `NC_`) rather than a sequence; found with `targetAccessions` set to all 27 accessions from `GET /statistics/target/accession/accession`. The organism filter does not match these (BRCA1, BRCA2, TP53, BAP1, PALB2, RAD51C/D and others). | 82 |

Together: 1,110 experiments, 712 distinct target names, 5,170,093 variants.

**Coverage check:** all 2,938 published score sets that the search API lists were enumerated and tested. Every one with a *Homo sapiens* sequence target (1,320) or an accession target (82) is in this dataset, and nothing here is missing from that list. Cross-check against `GET /experiments/`: its experiments reference 2,957 distinct score-set URNs; 2,938 are the searchable ones, and the other 19 (17 `tmp:<uuid>` placeholders and `urn:mavedb:00001208-a-1`, `urn:mavedb:00001209-a-1`) return 404 from `GET /score-sets/{urn}`, so every score set the public API can return is accounted for. The statistics endpoints report larger totals (2,995 published sets, 1,617 human targets); those most likely include draft or otherwise non-public records, but that was not confirmed.

## Files
- `score_sets.json`: the 1,320 search records for the organism filter (title, description, experiment, publications, license, target genes)
- `score_sets_accession_targets.json`: the 82 search records for the accession selection, with the accession list used
- `score_set_details.jsonl`, `score_set_details_accession_targets.jsonl`: full `GET /score-sets/{urn}` records for each selection, one JSON object per line (adds abstract, methods, dataset columns, DOI identifiers, mapping and processing state, score distribution and calibrations)
- `score_sets_index.csv`: flat index of all 1,402 sets (URN, selection, title, published date, variant count, target genes, accessions, assemblies, Ensembl/UniProt/RefSeq IDs, experiment, publications, license, states, score columns). Rebuilt by `build_index.py`.
- `taxonomy_homo_sapiens.json`: the taxonomy record used for the filter
- `download_scores.py`, `verify_scores.py`, `build_index.py`: see below
- `scores_manifest.jsonl`: one record per downloaded table (URN, kind, file name, row count, column names, gzipped size, sha256 of the gzip file)

## Score tables
`download_scores.py` fetches, for each of the 1,402 score sets, `GET /score-sets/{urn}/scores` and, for the 578 sets with count columns, `GET /score-sets/{urn}/counts`. Output goes to `scores/` as `urn_mavedb_<id>.<scores|counts>.csv.gz`: 1,980 files, about 237 MB gzipped, 5,170,093 score rows and 1,544,328 count rows. `scores/` is git-ignored and the container it was built in is ephemeral, so rerun `python3 download_scores.py` to recreate it (resumable; tens of minutes). `python3 verify_scores.py` re-checks everything (exit status 0 means no problems).

- Verification covers: the metadata files agree on the same set of score sets; the manifest lists exactly the expected files; every file's sha256, header and row count match the manifest and the score set's `numVariants`; all rows have the header's width. Last run: 0 problems.
- Whole-table requests for very large sets time out at the gateway, so sets over 100k variants (only `urn:mavedb:00000053-a-1`, 648,022 variants) are fetched in 50k-row pages. Its 648,022 accessions are all unique.
- Early failures (502/504/429) were transient gateway errors; the script retries with backoff.
- Columns differ per score set (for example `score`, `sd`, `se` or assay-specific columns); see the manifest and `score_sets_index.csv`.

## Caveats
- One score set, `urn:mavedb:00001239-a-6` (Insulin Binding Assay, 19,100 variants), has accessions of the form `tmp:<uuid>#N` instead of `<urn>#N`. This is how the API returns it (re-queried to confirm), not a download fault; all other rows use their own score set's URN.
- The search API's pagination order is unstable (a single pass returns duplicates and misses records). Sets were collected by unioning passes at several page sizes until the count matched `numScoreSets`. Compare against `numScoreSets` when refreshing.
- Only published score sets are included; private records need authentication.
- 322 sets have `processingState: incomplete`, and mapping state is `complete` for only 888 of 1,402 (321 incomplete, 181 failed, 9 queued, 1 processing, 2 none); check these fields before analysis.
- Each score set carries its own license (index column `license`): 1,340 CC0, 45 CC BY 4.0, 11 CC BY-SA 4.0, 6 CC BY-NC-SA 4.0. Check it before reuse.
