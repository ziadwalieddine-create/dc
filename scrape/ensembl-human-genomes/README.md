# Ensembl genome selector: human (taxonomy ID 9606)

Source page: https://www.ensembl.org/genome-selector/search?species_taxonomy_id=9606 (fetched 2026-10-25)

The page is a JavaScript app (the HTML is an empty shell), so this was scraped from the same APIs the page calls:
- `GET https://www.ensembl.org/api/search/genomes/v3?species_taxonomy_id=9606&per_page=100&page=N` (the results table; 563 hits over 6 pages)
- `GET https://www.ensembl.org/api/search/genomes/v3/download?species_taxonomy_id=9606&format=csv` (the page's own "download full table" export)
- `GET https://www.ensembl.org/api/metadata/{releases?current_only=true, popular_species, genome_counts, genome/{id}/{details,stats,ftplinks,example_objects,explain}}`

## Files
- `genomes.json`: all 563 human genome records (assembly, annotation provider/method, coding gene count, contig N50, variation/regulation flags, release, FTP URL)
- `genomes.csv`: the page's CSV export of the same table (563 rows)
- `meta_releases.json`, `meta_popular_species.json`, `meta_genome_counts.json`: site-level metadata
- `genome_details/`: details, stats, FTP links and example objects for the three principal assemblies: GRCh38.p14 (the only reference genome, GENCODE annotation), GRCh37.p13, and T2T-CHM13v2.0

## Contents summary
563 records, all in integrated release 2026-07: GRCh38.p14, GRCh37.p13, T2T-CHM13v2.0, and 560 haplotype assemblies of individual samples (440 have `hprc` in the assembly name; the rest are other sample assemblies such as HG002, HG005 and PAN027). 96 records have variation data and 2 have regulation data; 562 are annotated by Ensembl and 1 (GRCh38.p14) by GENCODE.

## Not included
No sequence, annotation (GTF/GFF) or variation files were downloaded. `ftp_url` / `ftplinks` give their locations on ftp.ebi.ac.uk.
