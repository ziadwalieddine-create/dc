# MaveDB: human-target score sets

Source page: https://www.mavedb.org/search (fetched 2026-10-25; MaveDB API v2026.3.2)

The page is a JavaScript app (empty HTML shell), so this was scraped from the API it uses (https://api.mavedb.org/api/v1, OpenAPI spec at /openapi.json). "Just the human genome" is applied as the page's own organism filter: published score sets with a target organism of *Homo sapiens* (NCBI taxonomy 9606).

- Search: `POST /score-sets/search` with `{"published": true, "targetOrganismNames": ["Homo sapiens"], "limit": N, "offset": M}` returned `numScoreSets = 1320`.
- Detail: `GET /score-sets/{urn}` for each of the 1320 URNs.
- Of 2,995 published score sets overall, 1,320 have a human target (1,029 experiments, 700 distinct target names, 4,989,344 variants in total).

## Files
- `score_sets.json`: the 1320 search records (title, description, experiment, publications, license, target genes)
- `score_set_details.jsonl`: the 1320 full detail records, one JSON object per line (adds abstract, methods, dataset columns, DOI identifiers, mapping and processing state, score distribution and calibrations)
- `score_sets_index.csv`: flat index of the above (URN, title, published date, variant count, target genes, Ensembl/UniProt/RefSeq IDs, experiment, publications, license, states, score columns)
- `taxonomy_homo_sapiens.json`: the taxonomy record used for the filter

## Caveats
- **No variant-level score tables are included.** Those are per score set (`GET /score-sets/{urn}/scores`, `/counts`, `/variants/data`) and total about 5 million variants, which is too large to commit here without a decision. The index gives the URNs, variant counts and column names needed to fetch them.
- The search API's pagination order is unstable (a single pass returns duplicates and misses records). The set was completed by unioning passes at several page sizes until all 1320 unique URNs were present. Compare against `numScoreSets` when refreshing.
- Only published score sets are included; private records need authentication.
- 322 sets have `processingState: incomplete`, and mapping state is `complete` for only 883; check these fields before analysis.
- Each score set carries its own license (index column `license`): 1289 CC0, 14 CC BY 4.0, 11 CC BY-SA 4.0, 6 CC BY-NC-SA 4.0. Check it before reuse.
