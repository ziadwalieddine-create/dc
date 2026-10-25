# NHGRI Talking Glossary: GWAS (scrape)

Source: https://www.genome.gov/genetics-glossary/Genome-Wide-Association-Studies (page updated 2026-10-24; fetched 2026-10-25)

- `page.html`: raw HTML of the entry
- `content.txt`: extracted main-content text (definition, narration, narrator, related terms)
- `gwas.json`: structured record
- `related_*.html`: raw HTML of the four linked glossary entries (Risk, Genomic Variation, Genome, Polygenic Risk Score)
- `related_*.txt`: extracted main-content text of each linked entry
- `gwas.json`: also holds each linked entry's definition, narration, narrator, related terms and asset names under `related`
- `assets/`: definition and narration MP3s, narrator portraits, hero and thumbnail images for all five entries (flat, shared across entries)

Linked entries are scraped one level deep only; terms they in turn link to (e.g. SNPs, Gene Pool, Genetics) are not included.

Site navigation, footer and subscription form are boilerplate and are retained only in the raw HTML.
