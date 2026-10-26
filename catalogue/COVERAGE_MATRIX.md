# Coverage matrix: evidence dimensions by lifecycle stage

*Generated 2026-10-26 15:25 UTC by `catalogue/build_catalogue.py` from `catalogue/sources/*.json`, `acquisition/targets.json`, `acquisition/manifest.jsonl`, `acquisition/validation_report.json` and `acquisition/compat_report.json`. Edit the sources, not this file.*

`✓` done · `~` partial · `·` not done. Stages: **cat**alogued, **acc**essible from this environment, **acq**uired and checksum-verified, **val**idated against the source description, **map**ped to dc assets in a demonstrated join, **use**able (acquired + validated + mapped + licence verified). Stages beyond *acc* are evaluated per dimension, using only the downloaded files that support that dimension. `integrated` (a documented, tested query layer that returns records by variant, gene, protein or region) is defined in `vocab.json` and credited to no source, because no such layer exists yet.

Counts are numbers of sources reaching the stage **for that dimension**. A dimension is not complete because a source is catalogued, or even because some of its files are downloaded.

| Dimension | Priority | Sources | cat | acc | acq | val | map | use |
|---|---|---|---|---|---|---|---|---|
| D01 Direct experimental variant-effect measurements | highest | 4 | 4 | 4 | 3 | 3 | 3 | 3 |
| D02 Protein structure, stability, activity and interactions | high | 8 | 8 | 8 | 5 | 5 | 5 | 5 |
| D03 Transcription, splicing, regulatory elements, tissue-specific expression | high | 7 | 7 | 7 | 4 | 4 | 4 | 4 |
| D04 Knockout, knockdown, CRISPRi/a, overexpression, knock-in and rescue experiments | high | 5 | 5 | 4 | 3 | 3 | 3 | 3 |
| D05 Human natural knockouts, constraint and dosage sensitivity | highest | 5 | 5 | 5 | 4 | 4 | 4 | 4 |
| D06 Genetic association, fine-mapping, LD and colocalisation | high | 7 | 7 | 7 | 3 | 3 | 2 | 2 |
| D07 Functional haplotypes and pharmacogenomics | medium | 5 | 5 | 5 | 4 | 4 | 4 | 4 |
| D08 Structural variation and copy-number variation | high | 7 | 7 | 6 | 4 | 4 | 3 | 3 |
| D09 Repeat expansions and short tandem repeats (STR/VNTR) | medium | 2 | 2 | 2 | 1 | 1 | 1 | 1 |
| D10 Pathways, molecular interactions, epistasis and compensation | medium | 5 | 5 | 5 | 3 | 3 | 3 | 3 |
| D11 Clinical variant interpretation and disease association | highest | 7 | 7 | 5 | 4 | 4 | 4 | 4 |
| D12 Population allele frequencies and haplotype structure | highest | 4 | 4 | 4 | 1 | 1 | 1 | 1 |
| D13 Integrated variant annotation and predictors | high | 6 | 6 | 6 | 2 | 2 | 2 | 2 |
| D14 Reference genomes, transcripts and gene models | highest | 4 | 4 | 4 | 4 | 4 | 4 | 4 |
| D15 Resource catalogues and classifications | low | 2 | 2 | 2 | 0 | 0 | 0 | 0 |
| D16 Background documentation and terminology | low | 1 | 1 | 1 | 1 | 1 | 0 | 0 |
| D17 Complex loci, paralogs and region reliability | high | 5 | 5 | 5 | 3 | 3 | 3 | 3 |

### D01 Direct experimental variant-effect measurements

Multiplexed/deep mutational scanning, MPRA and other assays that measure the effect of specific variants.

| Source | cat | acc | acq | val | map | use | Files supporting this dimension | Notes |
|---|---|---|---|---|---|---|---|---|
| `ensembl_vep` | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | 2/2 downloaded | FTP listings and the plugin documentation page are reachable (the page returned HTTP 504 in an earli  |
| `mavedb` | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | 4/4 downloaded | API open; Zenodo reachable with an identifying user agent (generic curl user agent was refused earli existing assets mavedb-human: downloaded_verified |
| `mpravardb` | ✓ | ✓ | · | · | · | · | none targeted | Home page read; the data are behind an interactive application.  |
| `proteingym` | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | 4/4 downloaded | Plain HTTPS downloads; the server publishes no checksum, so sizes are pinned and zip integrity is te  |

### D02 Protein structure, stability, activity and interactions

Protein sequence/function annotation, domains, structural models and stability predictors.

| Source | cat | acc | acq | val | map | use | Files supporting this dimension | Notes |
|---|---|---|---|---|---|---|---|---|
| `alphafold_db` | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | 1/1 downloaded | API returns model metadata; per-protein files for the MaveDB proteins are fetched by acquisition/fet  |
| `alphamissense` | ✓ | ✓ | · | · | · | · | none targeted | GCS bucket metadata and MD5 verified.  |
| `dbnsfp` | ✓ | ✓ | · | · | · | · | none targeted | Site and licence pages readable; bulk download gated by registration.  |
| `ensembl_vep` | ✓ | ✓ | · | · | · | · | none targeted | FTP listings and the plugin documentation page are reachable (the page returned HTTP 504 in an earli  |
| `hpa` | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | 1/1 downloaded | Plain HTTPS downloads; no checksum is published, so size and SHA-256 are pinned. The URLs always ser  |
| `interpro` | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | 1/1 downloaded | API returns JSON; per-protein matches for the MaveDB proteins are fetched by acquisition/fetch_per_p  |
| `pdbe_sifts` | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | 2/2 downloaded | Flat files listed with dates in the FTP index; the files are replaced each week, so the manifest rec  |
| `uniprot` | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | 1/1 downloaded | REST total reviewed human = 20,431.  |

### D03 Transcription, splicing, regulatory elements, tissue-specific expression

Expression by tissue, regulatory element maps, splicing effect predictions and sequence-to-function models.

| Source | cat | acc | acq | val | map | use | Files supporting this dimension | Notes |
|---|---|---|---|---|---|---|---|---|
| `encode` | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | 8/8 downloaded | API and SCREEN download links resolved by rendering the downloads page.  |
| `ensembl_vep` | ✓ | ✓ | · | · | · | · | none targeted | FTP listings and the plugin documentation page are reachable (the page returned HTTP 504 in an earli  |
| `eqtl_catalogue` | ✓ | ✓ | · | · | · | · | none targeted | FTP directory listing and licence page readable; the old REST API is gone (its replacement is docume  |
| `gtex` | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | 2/2 downloaded | API v2 dataset list and bucket listing verified.  |
| `hpa` | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | 2/2 downloaded | Plain HTTPS downloads; no checksum is published, so size and SHA-256 are pinned. The URLs always ser  |
| `mpravardb` | ✓ | ✓ | · | · | · | · | none targeted | Home page read; the data are behind an interactive application.  |
| `opentargets` | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | 1/1 downloaded | FTP, GraphQL and docs reachable. A focused subset of 27 Parquet datasets (8.5 GiB) is acquired and v  |

### D04 Knockout, knockdown, CRISPRi/a, overexpression, knock-in and rescue experiments

Genome-scale perturbation screens and model-organism knockout phenotypes.

| Source | cat | acc | acq | val | map | use | Files supporting this dimension | Notes |
|---|---|---|---|---|---|---|---|---|
| `alliance_genome` | ✓ | ✓ | · | · | · | · | none targeted | The release information API and the public bucket listing answered; the website pages are not readab  |
| `biogrid_orcs` | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | 1/1 downloaded | Download directory readable.  |
| `depmap` | ✓ | · | · | · | · | · | none targeted | Automated access met a human-verification page; not bypassed.  |
| `impc` | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | 4/4 downloaded | FTP result files for release 24.0; the server publishes no checksums, so sizes and SHA-256 are pinne  |
| `opentargets` | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | 5/5 downloaded | FTP, GraphQL and docs reachable. A focused subset of 27 Parquet datasets (8.5 GiB) is acquired and v  |

### D05 Human natural knockouts, constraint and dosage sensitivity

Loss-of-function tolerance, observed/expected metrics and curated haploinsufficiency/triplosensitivity.

| Source | cat | acc | acq | val | map | use | Files supporting this dimension | Notes |
|---|---|---|---|---|---|---|---|---|
| `clingen` | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | 3/3 downloaded | FTP and CSV export reachable.  |
| `ensembl_vep` | ✓ | ✓ | · | · | · | · | none targeted | FTP listings and the plugin documentation page are reachable (the page returned HTTP 504 in an earli  |
| `gnomad_constraint` | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | 2/2 downloaded | S3 object size verified by HEAD.  |
| `impc` | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | 4/4 downloaded | FTP result files for release 24.0; the server publishes no checksums, so sizes and SHA-256 are pinne  |
| `opentargets` | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | 1/1 downloaded | FTP, GraphQL and docs reachable. A focused subset of 27 Parquet datasets (8.5 GiB) is acquired and v  |

### D06 Genetic association, fine-mapping, LD and colocalisation

GWAS associations, credible sets, molecular QTLs, colocalisation and locus-to-gene evidence.

| Source | cat | acc | acq | val | map | use | Files supporting this dimension | Notes |
|---|---|---|---|---|---|---|---|---|
| `eqtl_catalogue` | ✓ | ✓ | · | · | · | · | none targeted | FTP directory listing and licence page readable; the old REST API is gone (its replacement is docume  |
| `gnomad_ld` | ✓ | ✓ | · | · | · | · | none targeted | S3 prefix listing read.  |
| `gtex` | ✓ | ✓ | · | · | · | · | none targeted | API v2 dataset list and bucket listing verified.  |
| `gwas_catalog` | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | 4/4 downloaded | FTP listing verified; web REST API root URL guessed incorrectly (404) and not needed.  |
| `ldlink` | ✓ | ✓ | · | · | · | · | none targeted | Home page read; API access needs a token.  |
| `opentargets` | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | 4/4 downloaded | FTP, GraphQL and docs reachable. A focused subset of 27 Parquet datasets (8.5 GiB) is acquired and v  |
| `pgs_catalog` | ✓ | ✓ | ✓ | ✓ | · | · | 9/9 downloaded | FTP listing, bulk metadata files and the harmonised file for the first score were read; the web page  |

### D07 Functional haplotypes and pharmacogenomics

Star-allele/haplotype definitions, gene-drug annotations and guidelines.

| Source | cat | acc | acq | val | map | use | Files supporting this dimension | Notes |
|---|---|---|---|---|---|---|---|---|
| `cpic` | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | 10/10 downloaded | API returns JSON; a plain GET returns every row. The CPIC licence page needs JavaScript and its cont  |
| `opentargets` | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | 1/1 downloaded | FTP, GraphQL and docs reachable. A focused subset of 27 Parquet datasets (8.5 GiB) is acquired and v  |
| `pharmcat` | ✓ | ✓ | · | · | · | · | none targeted | Documentation home page read.  |
| `pharmgkb` | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | 9/9 downloaded | HEAD sizes and last-modified verified.  |
| `pharmvar` | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | 1/1 downloaded | Terms at /terms-and-conditions (footer link). The Download Complete Database button calls get-downlo  |

### D08 Structural variation and copy-number variation

Population SV/CNV catalogues, dosage-sensitive regions and clinically curated CNVs.

| Source | cat | acc | acq | val | map | use | Files supporting this dimension | Notes |
|---|---|---|---|---|---|---|---|---|
| `1000genomes` | ✓ | ✓ | · | · | · | · | none targeted | FTP and portal reachable.  |
| `clingen` | ✓ | ✓ | ✓ | ✓ | · | · | 2/2 downloaded | FTP and CSV export reachable.  |
| `dbvar` | ✓ | ✓ | · | · | · | · | none targeted | Only directory listings read.  |
| `decipher` | ✓ | · | · | · | · | · | none targeted | Access Denied.  |
| `dgv` | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | 1/1 downloaded | Download links extracted from the downloads page HTML.  |
| `gnomad_sv_cnv` | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | 2/2 downloaded | S3 listings give exact sizes; both files downloaded, checksum-verified and validated.  |
| `hgsvc3` | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | 4/4 downloaded | FTP files with a MANIFEST of sizes and MD5s; both VCFs verified against it.  |

### D09 Repeat expansions and short tandem repeats (STR/VNTR)

Pathogenic repeat loci, reference tandem-repeat catalogues and population repeat distributions.

| Source | cat | acc | acq | val | map | use | Files supporting this dimension | Notes |
|---|---|---|---|---|---|---|---|---|
| `strchive` | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | 3/3 downloaded | Site pages rendered; download links extracted.  |
| `trexplorer` | ✓ | ✓ | · | · | · | · | none targeted | Web pages and the public bucket listing were read; no data file was downloaded.  |

### D10 Pathways, molecular interactions, epistasis and compensation

Curated pathways, protein/genetic interactions and network context.

| Source | cat | acc | acq | val | map | use | Files supporting this dimension | Notes |
|---|---|---|---|---|---|---|---|---|
| `biogrid_interactions` | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | 1/1 downloaded | Archive generated on request: no size or checksum is published, so the manifest pins size and SHA-25  |
| `intact` | ✓ | ✓ | · | · | · | · | none targeted | FTP release listing read.  |
| `opentargets` | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | 2/2 downloaded | FTP, GraphQL and docs reachable. A focused subset of 27 Parquet datasets (8.5 GiB) is acquired and v  |
| `reactome` | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | 3/3 downloaded | Version endpoint returned 97.  |
| `string_db` | ✓ | ✓ | · | · | · | · | none targeted | Download page read (file names and sizes); the licence page was not.  |

### D11 Clinical variant interpretation and disease association

Clinically classified variants, gene-disease validity and phenotype links.

| Source | cat | acc | acq | val | map | use | Files supporting this dimension | Notes |
|---|---|---|---|---|---|---|---|---|
| `clingen` | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | 8/8 downloaded | FTP and CSV export reachable.  |
| `clinvar` | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | 6/6 downloaded | VCF 200,670,148 bytes with published MD5; the VCV XML stream is checked against the MD5 NCBI publish  |
| `decipher` | ✓ | · | · | · | · | · | none targeted | Access Denied.  |
| `mitomap` | ✓ | · | · | · | · | · | none targeted | HTTP 403 to the automated probe; recorded and not circumvented.  |
| `opentargets` | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | 8/8 downloaded | FTP, GraphQL and docs reachable. A focused subset of 27 Parquet datasets (8.5 GiB) is acquired and v  |
| `pharmgkb` | ✓ | ✓ | · | · | · | · | none targeted | HEAD sizes and last-modified verified.  |
| `uniprot` | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | 1/1 downloaded | REST total reviewed human = 20,431.  |

### D12 Population allele frequencies and haplotype structure

Allele counts/frequencies across ancestries, phased haplotypes and reference panels.

| Source | cat | acc | acq | val | map | use | Files supporting this dimension | Notes |
|---|---|---|---|---|---|---|---|---|
| `1000genomes` | ✓ | ✓ | · | · | · | · | none targeted | FTP and portal reachable.  |
| `dbnsfp` | ✓ | ✓ | · | · | · | · | none targeted | Site and licence pages readable; bulk download gated by registration.  |
| `gnomad_ld` | ✓ | ✓ | · | · | · | · | none targeted | S3 prefix listing read.  |
| `gnomad_variants` | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | 3/3 downloaded | Remote tabix over the proxy works (index opens in about 1 s; a region read takes about 0.5 s). Slice  |

### D13 Integrated variant annotation and predictors

Tools and resources that aggregate predictors and annotations per variant.

| Source | cat | acc | acq | val | map | use | Files supporting this dimension | Notes |
|---|---|---|---|---|---|---|---|---|
| `alphamissense` | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | 2/2 downloaded | GCS bucket metadata and MD5 verified.  |
| `cadd` | ✓ | ✓ | · | · | · | · | none targeted | The download page answers a plain curl request but the probe collector's agent receives an empty rep  |
| `dbnsfp` | ✓ | ✓ | · | · | · | · | none targeted | Site and licence pages readable; bulk download gated by registration.  |
| `ensembl_vep` | ✓ | ✓ | · | · | · | · | none targeted | FTP listings and the plugin documentation page are reachable (the page returned HTTP 504 in an earli  |
| `gnomad_variants` | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | 2/2 downloaded | Remote tabix over the proxy works (index opens in about 1 s; a region read takes about 0.5 s). Slice  |
| `opentargets` | ✓ | ✓ | · | · | · | · | none targeted | FTP, GraphQL and docs reachable. A focused subset of 27 Parquet datasets (8.5 GiB) is acquired and v  |

### D14 Reference genomes, transcripts and gene models

Assemblies, transcript sets and identifier mapping required to place every evidence type on one coordinate system.

| Source | cat | acc | acq | val | map | use | Files supporting this dimension | Notes |
|---|---|---|---|---|---|---|---|---|
| `ensembl_genomes` | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | 5/5 downloaded | REST and FTP reachable. The web docs for VEP plugins returned HTTP 504 at probe time. The FASTA and  existing assets ensembl-human-genomes: downloaded_verified |
| `gencode_mane` | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | 2/2 downloaded | Directory listings verified.  |
| `opentargets` | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | 5/5 downloaded | FTP, GraphQL and docs reachable. A focused subset of 27 Parquet datasets (8.5 GiB) is acquired and v  |
| `refseq` | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | 6/6 downloaded | FTP files with MD5 checksums published in md5checksums.txt; all six files verified against them.  |

### D15 Resource catalogues and classifications

Registries of databases, standards and tools used to discover and classify resources.

| Source | cat | acc | acq | val | map | use | Files supporting this dimension | Notes |
|---|---|---|---|---|---|---|---|---|
| `biotools` | ✓ | ✓ | · | · | · | · | none targeted | OpenEBench reachable; bio.tools API challenge page returned HTTP 200 with no data.  |
| `fairsharing` | ✓ | ✓ | · | · | · | · | none targeted | Web pages rendered.  |

### D16 Background documentation and terminology

Educational and documentation material that explains concepts and projects but contains no variant-level evidence.

| Source | cat | acc | acq | val | map | use | Files supporting this dimension | Notes |
|---|---|---|---|---|---|---|---|---|
| `nhgri_education` | ✓ | ✓ | ✓ | ✓ | · | · | none targeted | Preserved unchanged; see inventory/EXISTING_ASSETS.md. existing assets nhgri-gwas, nhgri-1000-genomes: downloaded_verified |

### D17 Complex loci, paralogs and region reliability

Context that decides whether a variant call and its interpretation are reliable or confounded: segmental duplications and paralogs, gene-conversion and pseudogene-rich loci, HLA, KIR and immunoglobulin loci, collapsed or falsely duplicated assembly regions, low-mappability and tandem-repeat regions, alternate and patch sequences.

| Source | cat | acc | acq | val | map | use | Files supporting this dimension | Notes |
|---|---|---|---|---|---|---|---|---|
| `giab_stratifications` | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | 27/27 downloaded | FTP directory with an MD5 list; all downloaded files verified against it. The set README says the fi  |
| `hgsvc3` | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | 4/4 downloaded | FTP files with a MANIFEST of sizes and MD5s; both VCFs verified against it.  |
| `ipd_imgt_hla` | ✓ | ✓ | · | · | · | · | none targeted | FTP listing, LICENCE.md, md5 list and allele list head were read; nothing was downloaded beyond thos  |
| `ipd_kir` | ✓ | ✓ | · | · | · | · | none targeted | FTP listing, LICENCE.md and the allele list head were read; nothing larger was downloaded.  |
| `refseq` | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | 3/3 downloaded | FTP files with MD5 checksums published in md5checksums.txt; all six files verified against them.  |
