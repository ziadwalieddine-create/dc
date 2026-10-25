# MaveDB score tables: data notice

This branch (`data/mavedb-scores`) holds the downloaded MaveDB score and count tables for the human-target score sets listed in `scrape/mavedb-human/score_sets_index.csv` on `main`. It is not meant to be merged into `main`.

## Contents
`scrape/mavedb-human/scores/` contains 1980 gzipped CSV files, one per score set and kind (`<urn>.scores.csv.gz`, and `<urn>.counts.csv.gz` where the score set has count columns). They are the unmodified output of `GET https://api.mavedb.org/api/v1/score-sets/{urn}/scores` and `/counts`, fetched 2026-10-25 and compressed with gzip. Paged downloads of the largest set were concatenated without altering rows.

## Restore and verify
From a checkout of `main`:

```
git fetch origin data/mavedb-scores
git archive origin/data/mavedb-scores scrape/mavedb-human/scores | tar -x
python3 scrape/mavedb-human/verify_scores.py   # expects 0 problems
```

`git archive` extracts the files without staging them, and `scrape/mavedb-human/scores/` is git-ignored on `main`, so your working tree stays clean.

## Source and licences
All data originate from MaveDB (https://www.mavedb.org) and its contributors. Each score set carries its own licence; the licence and publication for every set are in `scrape/mavedb-human/score_sets_index.csv` and the full records (creators, contributors, DOIs) in `score_set_details*.jsonl` on `main`.

| Licence | Score sets |
|---|---|
| CC0 | 1340 |
| CC BY 4.0 | 45 |
| CC BY-SA 4.0 | 11 |
| CC BY-NC-SA 4.0 | 6 |

Terms that apply to redistribution:
- **CC0**: no restrictions.
- **CC BY 4.0**: attribute the creators of each score set (see the publication and creator fields) and indicate that the files were converted to gzip.
- **CC BY-SA 4.0**: as above, and any adaptation must be shared under the same licence.
- **CC BY-NC-SA 4.0**: as above, and for non-commercial use only.

Score sets that are not CC0 (62):

| Score set | Licence | Title | Publications |
|---|---|---|---|
| `urn:mavedb:00000658-0-1` | CC BY 4.0 | DDX3X SGE scores for clinical use | PubMed:38057330 |
| `urn:mavedb:00000662-0-1` | CC BY 4.0 | BAP1 SGE | PubMed:38969833 |
| `urn:mavedb:00000662-a-2` | CC BY 4.0 | BAP1 SGE Exon 1 | n/a |
| `urn:mavedb:00000662-b-1` | CC BY 4.0 | BAP1 SGE Exon 2 | n/a |
| `urn:mavedb:00000662-c-1` | CC BY 4.0 | BAP1 SGE Exon 3 | n/a |
| `urn:mavedb:00000662-d-1` | CC BY 4.0 | BAP1 SGE Exon 4 | n/a |
| `urn:mavedb:00000662-e-1` | CC BY 4.0 | BAP1 SGE Exon 5 | n/a |
| `urn:mavedb:00000662-f-1` | CC BY 4.0 | BAP1 SGE Exon 6 | n/a |
| `urn:mavedb:00000662-g-1` | CC BY 4.0 | BAP1 SGE Exon 7 | n/a |
| `urn:mavedb:00000662-h-1` | CC BY 4.0 | BAP1 SGE Exon 8 | n/a |
| `urn:mavedb:00000662-i-1` | CC BY 4.0 | BAP1 SGE Exon 11 | n/a |
| `urn:mavedb:00000662-j-1` | CC BY 4.0 | BAP1 SGE Exon 9 | n/a |
| `urn:mavedb:00000662-k-1` | CC BY 4.0 | BAP1 SGE Exon 10 | n/a |
| `urn:mavedb:00000662-l-1` | CC BY 4.0 | BAP1 SGE Exon 12 | n/a |
| `urn:mavedb:00000662-m-1` | CC BY 4.0 | BAP1 SGE Exon 13 | n/a |
| `urn:mavedb:00000662-n-1` | CC BY 4.0 | BAP1 SGE Exon 14 | n/a |
| `urn:mavedb:00000662-o-1` | CC BY 4.0 | BAP1 SGE Exon 15 | n/a |
| `urn:mavedb:00000662-p-1` | CC BY 4.0 | BAP1 SGE Exon 16 | n/a |
| `urn:mavedb:00000662-q-1` | CC BY 4.0 | BAP1 SGE Exon 17 | n/a |
| `urn:mavedb:00000662-r-1` | CC BY 4.0 | BAP1 SGE additional variant scores | PubMed:38969833 |
| `urn:mavedb:00000665-a-1` | CC BY 4.0 | pax6_le9_geneticin | PubMed:38849565 |
| `urn:mavedb:00000673-0-1` | CC BY 4.0 | RAD51C SGE | PubMed:39299233 |
| `urn:mavedb:00000673-a-2` | CC BY 4.0 | RAD51C SGE Exon 1 | n/a |
| `urn:mavedb:00000673-b-1` | CC BY 4.0 | RAD51C SGE Exon 6 | n/a |
| `urn:mavedb:00000673-c-1` | CC BY 4.0 | RAD51C SGE Exon 2 | n/a |
| `urn:mavedb:00000673-d-1` | CC BY 4.0 | RAD51C SGE Exon 3 | n/a |
| `urn:mavedb:00000673-e-1` | CC BY 4.0 | RAD51C SGE Exon 4 | n/a |
| `urn:mavedb:00000673-f-1` | CC BY 4.0 | RAD51C SGE Exon 5 | n/a |
| `urn:mavedb:00000673-g-1` | CC BY 4.0 | RAD51C SGE Exon 7 | n/a |
| `urn:mavedb:00000673-h-1` | CC BY 4.0 | RAD51C SGE Exon 8 | n/a |
| `urn:mavedb:00000673-i-1` | CC BY 4.0 | RAD51C SGE Exon 9 | n/a |
| `urn:mavedb:00000675-a-1` | CC BY 4.0 | Complete VHL SGE data | PubMed:38969834 |
| `urn:mavedb:00001197-a-5` | CC BY 4.0 | FKRP variants effects on alpha-DG glycosylation | bioRxiv:2023.07.12.548370 |
| `urn:mavedb:00001203-a-1` | CC BY 4.0 | Human CHEK2 CDS SNVs functionally classified as damaging or tolerated based on complementation of yeast rad53 | n/a |
| `urn:mavedb:00001203-a-2` | CC BY 4.0 | Rad53 complementation scores for coding SNVs of human CHEK2 | n/a |
| `urn:mavedb:00001204-a-1` | CC BY 4.0 | EGFR gefitinib PC9 cells | n/a |
| `urn:mavedb:00001204-a-2` | CC BY 4.0 | MAP2K1 dabrafenib and cetuximab HT-29 cells | n/a |
| `urn:mavedb:00001204-a-3` | CC BY 4.0 | MAP2K1 trametinib HT-29 cells | n/a |
| `urn:mavedb:00001204-a-4` | CC BY 4.0 | Dabrafenib and Cetuximab HT-29 cells, base editing z-scores (predicted consequence) | n/a |
| `urn:mavedb:00001204-a-5` | CC BY 4.0 | Dabrafenib and Cetuximab HT-29 cells, base editing z-scores (predicted consequence) | n/a |
| `urn:mavedb:00001204-a-6` | CC BY 4.0 | Trametinib HT-29 cells, base editing z-scores (predicted consequence) | n/a |
| `urn:mavedb:00001206-a-1` | CC BY 4.0 | MAP2K1 trametinib HT-29 cells | n/a |
| `urn:mavedb:00001206-a-2` | CC BY 4.0 | MAP2K1 trametinib HT-29 cells | n/a |
| `urn:mavedb:00001254-a-1` | CC BY 4.0 | LARGE1 variants effects on alpha-DG glycosylation | bioRxiv:2023.07.12.548370 |
| `urn:mavedb:00001283-a-1` | CC BY 4.0 | SGCA variants effects on α-sarcoglycan surface expression | n/a |
| `urn:mavedb:00000107-a-1` | CC BY-NC-SA 4.0 | PSAT1 200 Solid Growth Assay | PubMed:32077105 |
| `urn:mavedb:00000107-b-1` | CC BY-NC-SA 4.0 | A SNV-accessible mutational scan of human PSAT1 using functional complementation in yeast. | n/a |
| `urn:mavedb:00000112-a-1` | CC BY-NC-SA 4.0 | SNV-accessible missense variant scores for OTC | PubMed:37146589 |
| `urn:mavedb:00001207-a-1` | CC BY-NC-SA 4.0 | GPR68 mutational scan of cAMP signaling at pH 6.5 | Crossref:10.1101/2024.04.17.590000 |
| `urn:mavedb:00001207-a-2` | CC BY-NC-SA 4.0 | GPR68 mutational scan of surface eppression | Crossref:10.1101/2024.04.17.590000 |
| `urn:mavedb:00001207-a-3` | CC BY-NC-SA 4.0 | GPR68 cAMP signaling at pH 5.5 | Crossref:10.1101/2024.04.17.590000 |
| `urn:mavedb:00001256-a-1` | CC BY-SA 4.0 | FGFR2 activation | n/a |
| `urn:mavedb:00001256-a-10` | CC BY-SA 4.0 | FGFR1 pemigatinib resistance | n/a |
| `urn:mavedb:00001256-a-11` | CC BY-SA 4.0 | FGFR4 activation | n/a |
| `urn:mavedb:00001256-a-2` | CC BY-SA 4.0 | FGFR3 futibatinib resistance | n/a |
| `urn:mavedb:00001256-a-3` | CC BY-SA 4.0 | FGFR3 activation | n/a |
| `urn:mavedb:00001256-a-4` | CC BY-SA 4.0 | FGFR1 futibatinib resistance | n/a |
| `urn:mavedb:00001256-a-5` | CC BY-SA 4.0 | FGFR3 pemigatinib resistance | n/a |
| `urn:mavedb:00001256-a-6` | CC BY-SA 4.0 | FGFR4 futibatinib resistance | n/a |
| `urn:mavedb:00001256-a-7` | CC BY-SA 4.0 | FGFR2 futibatinib resistance | n/a |
| `urn:mavedb:00001256-a-8` | CC BY-SA 4.0 | FGFR4 pemigatinib resistance | n/a |
| `urn:mavedb:00001256-a-9` | CC BY-SA 4.0 | FGFR2 pemigatinib resistance | n/a |
