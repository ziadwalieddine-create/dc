# Coverage matrix: evidence dimensions by lifecycle stage

*Generated 2026-10-25 10:38 UTC by `catalogue/build_catalogue.py` from `catalogue/sources/*.json`, `acquisition/targets.json`, `acquisition/manifest.jsonl`, `acquisition/validation_report.json` and `acquisition/compat_report.json`. Edit the sources, not this file.*

`✓` done · `~` partial · `·` not done. Stages: **cat**alogued, **acc**essible from this environment, **acq**uired and checksum-verified, **val**idated against the source description, **map**ped to dc assets in a demonstrated join, **use**able (acquired + validated + mapped + licence verified). Stages beyond *acc* are evaluated per dimension, using only the downloaded files that support that dimension.

Counts are numbers of sources reaching the stage **for that dimension**. A dimension is not complete because a source is catalogued, or even because some of its files are downloaded.

| Dimension | Priority | Sources | cat | acc | acq | val | map | use |
|---|---|---|---|---|---|---|---|---|
| D01 Direct experimental variant-effect measurements | highest | 2 | 2 | 2 | 2 | 2 | 2 | 2 |
| D02 Protein structure, stability, activity and interactions | high | 6 | 6 | 6 | 1 | 1 | 1 | 1 |
| D03 Transcription, splicing, regulatory elements, tissue-specific expression | high | 5 | 5 | 5 | 2 | 2 | 2 | 1 |
| D04 Knockout, knockdown, CRISPRi/a, overexpression, knock-in and rescue experiments | high | 4 | 4 | 3 | 2 | 2 | 1 | 0 |
| D05 Human natural knockouts, constraint and dosage sensitivity | highest | 5 | 5 | 5 | 3 | 3 | 2 | 2 |
| D06 Genetic association, fine-mapping, LD and colocalisation | high | 4 | 4 | 4 | 1 | 1 | 1 | 1 |
| D07 Functional haplotypes and pharmacogenomics | medium | 4 | 4 | 4 | 3 | 3 | 2 | 2 |
| D08 Structural variation and copy-number variation | high | 6 | 6 | 5 | 3 | 3 | 1 | 1 |
| D09 Repeat expansions and short tandem repeats (STR/VNTR) | medium | 1 | 1 | 1 | 1 | 1 | 1 | 1 |
| D10 Pathways, molecular interactions, epistasis and compensation | medium | 2 | 2 | 2 | 1 | 1 | 1 | 1 |
| D11 Clinical variant interpretation and disease association | highest | 5 | 5 | 4 | 2 | 2 | 2 | 2 |
| D12 Population allele frequencies and haplotype structure | highest | 3 | 3 | 3 | 0 | 0 | 0 | 0 |
| D13 Integrated variant annotation and predictors | high | 4 | 4 | 4 | 1 | 1 | 1 | 1 |
| D14 Reference genomes, transcripts and gene models | highest | 3 | 3 | 3 | 3 | 3 | 1 | 0 |
| D15 Resource catalogues and classifications | low | 2 | 2 | 2 | 0 | 0 | 0 | 0 |
| D16 Background documentation and terminology | low | 1 | 1 | 1 | 0 | 0 | 0 | 0 |

### D01 Direct experimental variant-effect measurements

Multiplexed/deep mutational scanning, MPRA and other assays that measure the effect of specific variants.

| Source | cat | acc | acq | val | map | use | Files supporting this dimension | Notes |
|---|---|---|---|---|---|---|---|---|
| `ensembl_vep` | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | 2/2 downloaded | FTP listings reachable; MaveDB CHECKSUMS file lists BSD sum 22727 / 547461 KiB.  |
| `mavedb` | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | 3/3 downloaded | API open; Zenodo reachable with an identifying user agent (generic curl user agent was refused earli existing asset mavedb-human: downloaded_verified |

### D02 Protein structure, stability, activity and interactions

Protein sequence/function annotation, domains, structural models and stability predictors.

| Source | cat | acc | acq | val | map | use | Files supporting this dimension | Notes |
|---|---|---|---|---|---|---|---|---|
| `alphafold_db` | ✓ | ✓ | · | · | · | · | none targeted | API returned model metadata.  |
| `alphamissense` | ✓ | ✓ | · | · | · | · | none targeted | GCS bucket metadata and MD5 verified.  |
| `dbnsfp` | ✓ | ✓ | · | · | · | · | none targeted | Site and licence pages readable; bulk download gated by registration.  |
| `ensembl_vep` | ✓ | ✓ | · | · | · | · | none targeted | FTP listings reachable; MaveDB CHECKSUMS file lists BSD sum 22727 / 547461 KiB.  |
| `interpro` | ✓ | ✓ | · | · | · | · | none targeted | API root and entry count returned JSON.  |
| `uniprot` | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | 1/1 downloaded | REST total reviewed human = 20,431.  |

### D03 Transcription, splicing, regulatory elements, tissue-specific expression

Expression by tissue, regulatory element maps, splicing effect predictions and sequence-to-function models.

| Source | cat | acc | acq | val | map | use | Files supporting this dimension | Notes |
|---|---|---|---|---|---|---|---|---|
| `encode` | ✓ | ✓ | ✓ | ✓ | ✓ | · | 8/8 downloaded | API and SCREEN download links resolved by rendering the downloads page.  |
| `ensembl_vep` | ✓ | ✓ | · | · | · | · | none targeted | FTP listings reachable; MaveDB CHECKSUMS file lists BSD sum 22727 / 547461 KiB.  |
| `eqtl_catalogue` | ✓ | ✓ | · | · | · | · | none targeted | FTP listing verified; API retired.  |
| `gtex` | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | 2/2 downloaded | API v2 dataset list and bucket listing verified.  |
| `opentargets` | ✓ | ✓ | · | · | · | · | 0/1 downloaded | FTP, GraphQL and docs reachable. A focused subset of 26 Parquet datasets (about 8.5 GiB) is acquired  |

### D04 Knockout, knockdown, CRISPRi/a, overexpression, knock-in and rescue experiments

Genome-scale perturbation screens and model-organism knockout phenotypes.

| Source | cat | acc | acq | val | map | use | Files supporting this dimension | Notes |
|---|---|---|---|---|---|---|---|---|
| `biogrid_orcs` | ✓ | ✓ | ✓ | ✓ | ✓ | · | 1/1 downloaded | Download directory readable.  |
| `depmap` | ✓ | · | · | · | · | · | none targeted | Automated access met a human-verification page; not bypassed.  |
| `impc` | ✓ | ✓ | · | · | · | · | none targeted | Release page and FTP listing readable.  |
| `opentargets` | ✓ | ✓ | ✓ | ✓ | · | · | 5/5 downloaded | FTP, GraphQL and docs reachable. A focused subset of 26 Parquet datasets (about 8.5 GiB) is acquired  |

### D05 Human natural knockouts, constraint and dosage sensitivity

Loss-of-function tolerance, observed/expected metrics and curated haploinsufficiency/triplosensitivity.

| Source | cat | acc | acq | val | map | use | Files supporting this dimension | Notes |
|---|---|---|---|---|---|---|---|---|
| `clingen` | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | 3/3 downloaded | FTP and CSV export reachable.  |
| `ensembl_vep` | ✓ | ✓ | · | · | · | · | none targeted | FTP listings reachable; MaveDB CHECKSUMS file lists BSD sum 22727 / 547461 KiB.  |
| `gnomad_constraint` | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | 1/1 downloaded | S3 object size verified by HEAD.  |
| `impc` | ✓ | ✓ | · | · | · | · | none targeted | Release page and FTP listing readable.  |
| `opentargets` | ✓ | ✓ | ✓ | ✓ | · | · | 1/1 downloaded | FTP, GraphQL and docs reachable. A focused subset of 26 Parquet datasets (about 8.5 GiB) is acquired  |

### D06 Genetic association, fine-mapping, LD and colocalisation

GWAS associations, credible sets, molecular QTLs, colocalisation and locus-to-gene evidence.

| Source | cat | acc | acq | val | map | use | Files supporting this dimension | Notes |
|---|---|---|---|---|---|---|---|---|
| `eqtl_catalogue` | ✓ | ✓ | · | · | · | · | none targeted | FTP listing verified; API retired.  |
| `gtex` | ✓ | ✓ | · | · | · | · | none targeted | API v2 dataset list and bucket listing verified.  |
| `gwas_catalog` | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | 1/1 downloaded | FTP listing verified; web REST API root URL guessed incorrectly (404) and not needed.  |
| `opentargets` | ✓ | ✓ | ~ | · | · | · | 2/4 downloaded | FTP, GraphQL and docs reachable. A focused subset of 26 Parquet datasets (about 8.5 GiB) is acquired  |

### D07 Functional haplotypes and pharmacogenomics

Star-allele/haplotype definitions, gene-drug annotations and guidelines.

| Source | cat | acc | acq | val | map | use | Files supporting this dimension | Notes |
|---|---|---|---|---|---|---|---|---|
| `cpic` | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | 2/2 downloaded | API returns JSON.  |
| `opentargets` | ✓ | ✓ | ✓ | ✓ | · | · | 1/1 downloaded | FTP, GraphQL and docs reachable. A focused subset of 26 Parquet datasets (about 8.5 GiB) is acquired  |
| `pharmgkb` | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | 5/5 downloaded | HEAD sizes and last-modified verified.  |
| `pharmvar` | ✓ | ✓ | · | · | · | · | none targeted | Home and download pages readable; terms page not found.  |

### D08 Structural variation and copy-number variation

Population SV/CNV catalogues, dosage-sensitive regions and clinically curated CNVs.

| Source | cat | acc | acq | val | map | use | Files supporting this dimension | Notes |
|---|---|---|---|---|---|---|---|---|
| `1000genomes` | ✓ | ✓ | · | · | · | · | none targeted | FTP and portal reachable.  |
| `clingen` | ✓ | ✓ | ✓ | ✓ | · | · | 2/2 downloaded | FTP and CSV export reachable.  |
| `dbvar` | ✓ | ✓ | · | · | · | · | none targeted | Only directory listings read.  |
| `decipher` | ✓ | · | · | · | · | · | none targeted | Access Denied.  |
| `dgv` | ✓ | ✓ | ✓ | ✓ | · | · | 1/1 downloaded | Download links extracted from the downloads page HTML.  |
| `gnomad_sv_cnv` | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | 1/1 downloaded | S3 listings give exact sizes.  |

### D09 Repeat expansions and short tandem repeats (STR/VNTR)

Pathogenic repeat loci, reference tandem-repeat catalogues and population repeat distributions.

| Source | cat | acc | acq | val | map | use | Files supporting this dimension | Notes |
|---|---|---|---|---|---|---|---|---|
| `strchive` | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | 3/3 downloaded | Site pages rendered; download links extracted.  |

### D10 Pathways, molecular interactions, epistasis and compensation

Curated pathways, protein/genetic interactions and network context.

| Source | cat | acc | acq | val | map | use | Files supporting this dimension | Notes |
|---|---|---|---|---|---|---|---|---|
| `opentargets` | ✓ | ✓ | · | · | · | · | 0/2 downloaded | FTP, GraphQL and docs reachable. A focused subset of 26 Parquet datasets (about 8.5 GiB) is acquired  |
| `reactome` | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | 3/3 downloaded | Version endpoint returned 97.  |

### D11 Clinical variant interpretation and disease association

Clinically classified variants, gene-disease validity and phenotype links.

| Source | cat | acc | acq | val | map | use | Files supporting this dimension | Notes |
|---|---|---|---|---|---|---|---|---|
| `clingen` | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | 1/1 downloaded | FTP and CSV export reachable.  |
| `clinvar` | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | 3/3 downloaded | VCF 200,670,148 bytes with published MD5.  |
| `decipher` | ✓ | · | · | · | · | · | none targeted | Access Denied.  |
| `opentargets` | ✓ | ✓ | ~ | · | · | · | 7/8 downloaded | FTP, GraphQL and docs reachable. A focused subset of 26 Parquet datasets (about 8.5 GiB) is acquired  |
| `pharmgkb` | ✓ | ✓ | · | · | · | · | none targeted | HEAD sizes and last-modified verified.  |

### D12 Population allele frequencies and haplotype structure

Allele counts/frequencies across ancestries, phased haplotypes and reference panels.

| Source | cat | acc | acq | val | map | use | Files supporting this dimension | Notes |
|---|---|---|---|---|---|---|---|---|
| `1000genomes` | ✓ | ✓ | · | · | · | · | none targeted | FTP and portal reachable.  |
| `dbnsfp` | ✓ | ✓ | · | · | · | · | none targeted | Site and licence pages readable; bulk download gated by registration.  |
| `gnomad_variants` | ✓ | ✓ | · | · | · | · | none targeted | Downloads page rendered; S3 listing reachable.  |

### D13 Integrated variant annotation and predictors

Tools and resources that aggregate predictors and annotations per variant.

| Source | cat | acc | acq | val | map | use | Files supporting this dimension | Notes |
|---|---|---|---|---|---|---|---|---|
| `alphamissense` | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | 2/2 downloaded | GCS bucket metadata and MD5 verified.  |
| `dbnsfp` | ✓ | ✓ | · | · | · | · | none targeted | Site and licence pages readable; bulk download gated by registration.  |
| `ensembl_vep` | ✓ | ✓ | · | · | · | · | none targeted | FTP listings reachable; MaveDB CHECKSUMS file lists BSD sum 22727 / 547461 KiB.  |
| `opentargets` | ✓ | ✓ | · | · | · | · | none targeted | FTP, GraphQL and docs reachable. A focused subset of 26 Parquet datasets (about 8.5 GiB) is acquired  |

### D14 Reference genomes, transcripts and gene models

Assemblies, transcript sets and identifier mapping required to place every evidence type on one coordinate system.

| Source | cat | acc | acq | val | map | use | Files supporting this dimension | Notes |
|---|---|---|---|---|---|---|---|---|
| `ensembl_genomes` | ✓ | ✓ | ✓ | ✓ | · | · | none targeted | REST and FTP reachable. The web docs for VEP plugins returned HTTP 504 at probe time. existing asset ensembl-human-genomes: downloaded_verified |
| `gencode_mane` | ✓ | ✓ | ✓ | ✓ | ✓ | · | 2/2 downloaded | Directory listings verified.  |
| `opentargets` | ✓ | ✓ | ✓ | ✓ | · | · | 5/5 downloaded | FTP, GraphQL and docs reachable. A focused subset of 26 Parquet datasets (about 8.5 GiB) is acquired  |

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
| `nhgri_education` | ✓ | ✓ | · | · | · | · | none targeted | Preserved unchanged; see inventory/EXISTING_ASSETS.md.  |
