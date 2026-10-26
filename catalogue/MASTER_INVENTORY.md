# Master inventory of reusable genomic resources

*Generated 2026-10-26 02:55 UTC by `catalogue/build_catalogue.py` from `catalogue/sources/*.json`, `acquisition/targets.json`, `acquisition/manifest.jsonl`, `acquisition/validation_report.json` and `acquisition/compat_report.json`. Edit the sources, not this file.*

`✓` done · `~` partial · `·` not done. Stages: **cat**alogued, **acc**essible from this environment, **acq**uired and checksum-verified, **val**idated against the source description, **map**ped to dc assets in a demonstrated join, **use**able (acquired + validated + mapped + licence verified). Stages beyond *acc* are evaluated per dimension, using only the downloaded files that support that dimension. `integrated` (a documented, tested query layer that returns records by variant, gene, protein or region) is defined in `vocab.json` and credited to no source, because no such layer exists yet.

53 sources catalogued. A source appears under its primary evidence dimension. Columns *acq/val/map* here cover all of the source's downloaded files; "Usable for" lists the dimensions for which the source is actually usable. See `COVERAGE_MATRIX.md` for the per-dimension view.

## D01 Direct experimental variant-effect measurements

| Source | Version / date | Access · licence | Variants · assemblies | cat | acc | acq | val | map | Usable for | Decision |
|---|---|---|---|---|---|---|---|---|---|---|
| **MaveDB (Multiplexed Assays of Variant Effect database)** (`mavedb`) | API 2026.3.2; bulk dump v4 (Zenodo 20840937) · 2026-06-25 (dump); API queried 2026-10-25 | open · Per record: CC0 (1,340), CC BY 4.0 (45), CC BY-SA 4.0 (11), CC BY-NC-S | snv, protein_level, indel · GRCh38, none/not_applicable | ✓ | ✓ | ✓ | ✓ | ✓ | D01 | reuse_existing |
| **MPRAVarDB: massively parallel reporter assay variant database** (`mpravardb`) | Not stated on the home page (unverified) · unknown | unverified · Not stated on the pages read. (unverified) | snv, indel, region_level · none/not_applicable | ✓ | ✓ | · | · | · | — | defer |
| **ProteinGym v1.3: deep mutational scanning benchmark assays (substitutions and indels)** (`proteingym`) | v1.3 (files dated 2025-04-13); the repository README lists r · 2025-04-13 | open_with_terms · MIT licence for the ProteinGym project (LICENSE of the OATML-Markslab/ | protein_level, indel · none/not_applicable | ✓ | ✓ | ✓ | ✓ | ✓ | D01 | adopt_curated |

## D02 Protein structure, stability, activity and interactions

| Source | Version / date | Access · licence | Variants · assemblies | cat | acc | acq | val | map | Usable for | Decision |
|---|---|---|---|---|---|---|---|---|---|---|
| **AlphaFold Protein Structure Database** (`alphafold_db`) | Model version 6 (P04637 model created 2025-08-01); human pro · 2025-08 | open · CC BY 4.0, academic and commercial use (FAQ) | protein_level · none/not_applicable | ✓ | ✓ | ✓ | ✓ | ✓ | D02 | retrieve_source |
| **InterPro** (`interpro`) | InterPro 110.0 (interpro-version response header of the API, · 2026-10-25 (query date) | open · CC0 1.0 Universal for InterPro, Pfam, PRINTS and SFLD downloadable dat | protein_level · none/not_applicable | ✓ | ✓ | ✓ | ✓ | ✓ | D02 | retrieve_source |
| **PDBe SIFTS: PDB chain to UniProt residue mappings** (`pdbe_sifts`) | SIFTS flat files of 2026-10-18 (PDB release 42.26, UniProt 2 · 2026-10-18 | open · EMBL-EBI Terms of Use: freely available, no additional restriction on  | protein_level · none/not_applicable | ✓ | ✓ | ✓ | ✓ | ✓ | D02 | retrieve_source |
| **UniProtKB (Swiss-Prot and TrEMBL)** (`uniprot`) | Release 2026_03 (150,006,383 entries; Swiss-Prot 575,748) · 2026-09-02 | open · CC BY 4.0 for copyrightable parts | protein_level, gene_level · none/not_applicable | ✓ | ✓ | ✓ | ✓ | ✓ | D02, D11 | retrieve_source |

## D03 Transcription, splicing, regulatory elements, tissue-specific expression

| Source | Version / date | Access · licence | Variants · assemblies | cat | acc | acq | val | map | Usable for | Decision |
|---|---|---|---|---|---|---|---|---|---|---|
| **ENCODE and SCREEN candidate cis-regulatory elements** (`encode`) | Registry V4 cCREs (Nature 2026 citation); ENCODE portal quer · 2026 | open · ENCODE Data Use Policy: 'all data produced will be available for unres | region_level · GRCh38, GRCh37 | ✓ | ✓ | ✓ | ✓ | ✓ | D03 | retrieve_source |
| **GTEx (Genotype-Tissue Expression)** (`gtex`) | v11 (946 subjects, 54 tissues, 50 eQTL tissues, GRCh38, GENC · 2025-08-22 (v11 RNA-seq files) | open · Open-access files free to download; raw sequence and full donor metada | snv, indel, gene_level · GRCh38, GRCh37 | ✓ | ✓ | ✓ | ✓ | ✓ | D03 | retrieve_source |
| **Human Protein Atlas v25.1: tissue and cell-type RNA expression and protein subcellular location** (`hpa`) | Human Protein Atlas version 25.1 (the download page names th · 2025-12-12 | open_with_terms · Creative Commons Attribution 4.0 International for all copyrightable p | gene_level · none/not_applicable | ✓ | ✓ | ✓ | ✓ | ✓ | D03, D02 | retrieve_source |

## D04 Knockout, knockdown, CRISPRi/a, overexpression, knock-in and rescue experiments

| Source | Version / date | Access · licence | Variants · assemblies | cat | acc | acq | val | map | Usable for | Decision |
|---|---|---|---|---|---|---|---|---|---|---|
| **Alliance of Genome Resources (model-organism and human gene knowledge)** (`alliance_genome`) | Release 9.1.0 (API release info dated 2026-09-02) · 2026-09-02 | unverified · Not retrieved: the Alliance website pages are JavaScript-only and one  (unverified) | gene_level · none/not_applicable | ✓ | ✓ | · | · | · | — | defer |
| **BioGRID ORCS (Open Repository of CRISPR Screens)** (`biogrid_orcs`) | Index version 2.0.18 · 2026 | open · MIT License (stated on the BioGRID download repository page for BioGRI | gene_level · none/not_applicable | ✓ | ✓ | ✓ | ✓ | ✓ | D04 | retrieve_source |
| **DepMap (Cancer Dependency Map)** (`depmap`) | Not retrieved: the portal presented a human-verification pag (unverified) · unknown | blocked_for_automation · Terms page not readable (verification wall); commonly CC BY 4.0 (unverified) | gene_level · none/not_applicable | ✓ | · | · | · | · | — | adopt_curated |
| **IMPC (International Mouse Phenotyping Consortium)** (`impc`) | Data Release 24.0 (GRCm39; OpenStats 1.20.0) · 2026-03-16 | open · Creative Commons Attribution 4.0 (site footer statement) | gene_level · other | ✓ | ✓ | ✓ | ✓ | ✓ | D04, D05 | adopt_curated |

## D05 Human natural knockouts, constraint and dosage sensitivity

| Source | Version / date | Access · licence | Variants · assemblies | cat | acc | acq | val | map | Usable for | Decision |
|---|---|---|---|---|---|---|---|---|---|---|
| **ClinGen (dosage sensitivity, gene-disease validity, region curation)** (`clingen`) | Files regenerated daily (gene curation list dated 2026-10-25 · 2026-10-25 | open · ClinGen terms: curated content released openly, free to capture and us | gene_level, region_level, cnv · GRCh37, GRCh38 | ✓ | ✓ | ✓ | ✓ | ✓ | D05, D11 | retrieve_source |
| **gnomAD constraint metrics (LoF/missense) and pext** (`gnomad_constraint`) | v4.1.2 constraint metrics (TSV.bgz 574.5 MB, 2026-09-23); v4 · 2026-09-23 | open · CC0 for primary data (policies page) | gene_level · GRCh38 | ✓ | ✓ | ✓ | ✓ | ✓ | D05 | retrieve_source |

## D06 Genetic association, fine-mapping, LD and colocalisation

| Source | Version / date | Access · licence | Variants · assemblies | cat | acc | acq | val | map | Usable for | Decision |
|---|---|---|---|---|---|---|---|---|---|---|
| **eQTL Catalogue** (`eqtl_catalogue`) | Release 7 and pre-release r8_beta (FTP directory dated 2026- · 2026-01-28 | open · Data under Creative Commons Attribution 4.0 International; code under  | snv, indel · GRCh38 | ✓ | ✓ | · | · | · | — | adopt_curated |
| **gnomAD v2.1.1 linkage disequilibrium matrices** (`gnomad_ld`) | gnomAD v2.1.1 LD matrices (S3 prefix release/2.1.1/ld/, 39 e · unknown | open · gnomAD data are released under CC0 1.0 (as listed for gnomAD in the Op (unverified) | snv · GRCh37 | ✓ | ✓ | · | · | · | — | defer |
| **NHGRI-EBI GWAS Catalog** (`gwas_catalog`) | Release files dated 2026-10-22 · 2026-10-22 | open · EMBL-EBI general terms of use; summary statistics CC0 unless otherwise | snv, indel · GRCh38 | ✓ | ✓ | ✓ | ✓ | ✓ | D06 | adopt_curated |
| **LDlink (NCI): linkage disequilibrium tools for 1000 Genomes populations** (`ldlink`) | LDlink 7.4.0 (as named on the home page) · unknown | registration_required · Not retrieved: usage conditions of the API (token, rate limits) were n (unverified) | snv, indel · GRCh37, GRCh38 | ✓ | ✓ | · | · | · | — | do_not_acquire |

## D07 Functional haplotypes and pharmacogenomics

| Source | Version / date | Access · licence | Variants · assemblies | cat | acc | acq | val | map | Usable for | Decision |
|---|---|---|---|---|---|---|---|---|---|---|
| **CPIC (Clinical Pharmacogenetics Implementation Consortium)** (`cpic`) | Public API tables queried 2026-10-25 (PostgREST 12.0.2); the · 2026-10-25 | open_with_terms · The CPIC licence URL resolves to the ClinPGx data-use text (CC BY-SA 4 | gene_level, snv, indel, cnv, mtdna · none/not_applicable | ✓ | ✓ | ✓ | ✓ | ✓ | D07 | retrieve_source |
| **PharmCAT (Pharmacogenomics Clinical Annotation Tool)** (`pharmcat`) | PharmCAT v3.4.0 (as named on the documentation home page) · unknown | unverified · Not retrieved: the documentation pages were read but the software lice (unverified) | snv, indel, cnv, gene_level · none/not_applicable | ✓ | ✓ | · | · | · | — | do_not_acquire |
| **ClinPGx / PharmGKB** (`pharmgkb`) | Archives dated 2026-10-05 · 2026-10-05 | open_with_terms · CC BY-SA 4.0 (data usage policy) | snv, indel, gene_level · none/not_applicable | ✓ | ✓ | ✓ | ✓ | ✓ | D07 | retrieve_source |
| **PharmVar (Pharmacogene Variation Consortium)** (`pharmvar`) | PharmVar 6.2.29 (site footer "Last Updated Sep 22 2026"; VCF · 2026-09-22 | open_with_terms · CC BY-SA 4.0 for the PharmVar database content, with the user agreeing | snv, indel, gene_level · GRCh37, GRCh38, none/not_applicable | ✓ | ✓ | ✓ | ✓ | ✓ | D07 | retrieve_source |

## D08 Structural variation and copy-number variation

| Source | Version / date | Access · licence | Variants · assemblies | cat | acc | acq | val | map | Usable for | Decision |
|---|---|---|---|---|---|---|---|---|---|---|
| **dbVar (NCBI database of genomic structural variation)** (`dbvar`) | FTP root lists ISCA/, archive/, clingen/, data/, sandbox/; d · 2021-01-21 | open · NCBI data usage policy: NCBI places no restrictions on the use or dist | sv, cnv · GRCh37, GRCh38 | ✓ | ✓ | · | · | · | — | defer |
| **DECIPHER** (`decipher`) | Not retrieved (unverified) · unknown | restricted_agreement · Automated access to the downloads page returned {"error":"Access Denie (unverified) | cnv, snv, indel · GRCh37, GRCh38 | ✓ | · | · | · | · | — | do_not_acquire |
| **DGV (Database of Genomic Variants)** (`dgv`) | Latest release 2025-12-01 (GRCh38 variants file 364,191,405  · 2025-12-01 | open_with_terms · TCAG Notice and Disclaimer of Database Use (data provided 'AS IS'; app | cnv, sv · GRCh37, GRCh38 | ✓ | ✓ | ✓ | ✓ | ✓ | D08 | retrieve_source |
| **gnomAD structural variants and copy-number variants** (`gnomad_sv_cnv`) | SVs v4.1.0, CNVs v4.1.0 (site header); also SVs v2.1 · 2024-04-25 (CNV VCF) | open · CC0 for primary data (policies page) | sv, cnv · GRCh38 | ✓ | ✓ | ✓ | ✓ | ✓ | D08 | retrieve_source |
| **HGSVC3 release 1.0: assembly-based structural variants from 65 phased genomes (GRCh38)** (`hgsvc3`) | Variant calls release 1.0 (directory dated 2024-09-18; PAV 2 · 2024-09-18 | open_with_terms · Data reuse statement of the project (README_HGSVC3_datareuse_statement | sv, complex_locus · GRCh38 | ✓ | ✓ | ✓ | ✓ | ✓ | D08, D17 | retrieve_source |

## D09 Repeat expansions and short tandem repeats (STR/VNTR)

| Source | Version / date | Access · licence | Variants · assemblies | cat | acc | acq | val | map | Usable for | Decision |
|---|---|---|---|---|---|---|---|---|---|---|
| **STRchive and related tandem-repeat catalogues** (`strchive`) | Site v2.28.0; 82 disease loci · 2026 | open · CC BY 4.0 (site footer); research use only | str · GRCh37, GRCh38, T2T-CHM13 | ✓ | ✓ | ✓ | ✓ | ✓ | D09 | retrieve_source |
| **TRExplorer: genome-wide tandem-repeat catalogue** (`trexplorer`) | Catalogue versions v1.0, v2.0 and v2.1 in the public bucket  · unknown | unverified · Not found: neither the about page, the FAQ nor the bucket listing stat (unverified) | str, region_level · GRCh38, T2T-CHM13 | ✓ | ✓ | · | · | · | — | defer |

## D10 Pathways, molecular interactions, epistasis and compensation

| Source | Version / date | Access · licence | Variants · assemblies | cat | acc | acq | val | map | Usable for | Decision |
|---|---|---|---|---|---|---|---|---|---|---|
| **BioGRID protein and genetic interactions (release 5.0.262, human file)** (`biogrid_interactions`) | BioGRID release 5.0.262 (read from the member names of the o · unknown (generated on request) | open · MIT License: "100% freely available to both academic and commercial us | gene_level, protein_level · none/not_applicable | ✓ | ✓ | ✓ | ✓ | ✓ | D10 | retrieve_source |
| **IntAct molecular interaction database** (`intact`) | Current release directory listed 2026-01-14 (all.zip 15 GB) · 2026-01-14 | open_with_terms · EMBL-EBI terms of use, as listed for IntAct in the Open Targets data-s (unverified) | gene_level, protein_level · none/not_applicable | ✓ | ✓ | · | · | · | — | defer |
| **Reactome** (`reactome`) | Version 97 (mapping files dated 2026-06-19/20) · 2026-06 | open · Data licensed CC0; illustrations CC BY 4.0 (license page) | gene_level, protein_level · none/not_applicable | ✓ | ✓ | ✓ | ✓ | ✓ | D10 | retrieve_source |
| **STRING protein-protein association networks** (`string_db`) | STRING v12.5 (file names of the download page) · unknown | open_with_terms · CC BY 4.0 as listed for STRING in the Open Targets data-source licence (unverified) | gene_level, protein_level · none/not_applicable | ✓ | ✓ | · | · | · | — | defer |

## D11 Clinical variant interpretation and disease association

| Source | Version / date | Access · licence | Variants · assemblies | cat | acc | acq | val | map | Usable for | Decision |
|---|---|---|---|---|---|---|---|---|---|---|
| **ClinVar** (`clinvar`) | VCF fileDate 2026-10-19 · 2026-10-19 | open · NCBI public domain policy; submitters may claim rights in some content | snv, indel, sv, cnv · GRCh37, GRCh38 | ✓ | ✓ | ✓ | ✓ | ✓ | D11 | retrieve_source |
| **MITOMAP: human mitochondrial genome database** (`mitomap`) | Not retrieved: the site answered HTTP 403 to the automated r (unverified) · unknown | blocked_for_automation · Not retrieved. (unverified) | mtdna · none/not_applicable | ✓ | · | · | · | · | — | defer |

## D12 Population allele frequencies and haplotype structure

| Source | Version / date | Access · licence | Variants · assemblies | cat | acc | acq | val | map | Usable for | Decision |
|---|---|---|---|---|---|---|---|---|---|---|
| **1000 Genomes Project / IGSR (International Genome Sample Resource)** (`1000genomes`) | FTP data_collections include 1000G_2504_high_coverage (3,202 · 2020-05-25 (sample index file) | open · IGSR provides open data; terms of use are in the IGSR data disclaimer  | snv, indel, sv · GRCh37, GRCh38 | ✓ | ✓ | ✓ | ✓ | ✓ | — | extend_integration |
| **gnomAD v4.1 (exomes and genomes) allele frequencies** (`gnomad_variants`) | v4.1 (4.1.2 for constraint); site shows v4.1.2, v3.1.2, v2.1 · 2026-09-23 (4.1.2 constraint) | open · Primary exome/genome data released under CC0; some annotations carry s | snv, indel · GRCh38, GRCh37 | ✓ | ✓ | ✓ | ✓ | ✓ | D12, D13 | retrieve_source |

## D13 Integrated variant annotation and predictors

| Source | Version / date | Access · licence | Variants · assemblies | cat | acc | acq | val | map | Usable for | Decision |
|---|---|---|---|---|---|---|---|---|---|---|
| **AlphaMissense** (`alphamissense`) | Static release; GCS objects dated 2023-08-03 · 2023-08-03 | open_with_terms · Conflicting statements, applied conservatively as CC BY-NC-SA 4.0 for  | snv, protein_level · GRCh37, GRCh38 | ✓ | ✓ | ✓ | ✓ | ✓ | D13 | retrieve_source |
| **CADD: combined annotation dependent depletion** (`cadd`) | CADD v1.7 (the download page names v1.7 as the current relea · unknown | open_with_terms · "CADD scores are freely available for all non-commercial applications" | snv, indel · GRCh37, GRCh38 | ✓ | ✓ | · | · | · | — | defer |
| **dbNSFP** (`dbnsfp`) | August 2026 release (GENCODE 50 / Ensembl 116); column docs  · 2026-08 | registration_required · Academic branch CC BY-NC-ND 4.0 (about 50 GB); commercial branch under | snv · GRCh38, GRCh37 | ✓ | ✓ | · | · | · | — | adopt_curated |
| **Ensembl VEP and release-116 plugin data** (`ensembl_vep`) | Ensembl release 116 (REST); plugin data directories dated 20 · 2026-08-07 (CHECKSUMS) | open_with_terms · Ensembl: no restrictions on data (verified); plugin resources carry up | snv, indel, sv, gene_level · GRCh37, GRCh38 | ✓ | ✓ | ✓ | ✓ | ✓ | D01 | adopt_curated |
| **Open Targets Platform** (`opentargets`) | Data 26.09 (API 26.9.0) · 2026-09 | open · CC0 1.0 for Platform data (docs page); a separate data-source licensin | snv, indel, gene_level, region_level, protein_level · GRCh38 | ✓ | ✓ | ✓ | ✓ | ✓ | D06, D04, D10, D11, D05, D07, D03, D14 | adopt_curated |

## D14 Reference genomes, transcripts and gene models

| Source | Version / date | Access · licence | Variants · assemblies | cat | acc | acq | val | map | Usable for | Decision |
|---|---|---|---|---|---|---|---|---|---|---|
| **Ensembl genome registry, assemblies and gene annotation (human)** (`ensembl_genomes`) | Ensembl release 116 (REST); integrated release 2026-07 in th · 2026-07 | open · Ensembl imposes no restrictions on access to or use of the data (discl | gene_level, region_level, complex_locus · GRCh37, GRCh38, T2T-CHM13, other | ✓ | ✓ | ✓ | ✓ | ✓ | D14 | retrieve_source |
| **GENCODE and MANE transcript sets** (`gencode_mane`) | GENCODE v50 (files dated 2026-06-09); MANE v1.5 (2025-12-04) · 2026-06-09 | open · GENCODE: 'all the GENCODE project data is open access' (GENCODE data-a | gene_level, region_level · GRCh37, GRCh38 | ✓ | ✓ | ✓ | ✓ | ✓ | D14 | retrieve_source |
| **NCBI RefSeq: GRCh38.p14 accession map, alternate-locus placement and gene annotation (release RS_2025_08)** (`refseq`) | GCF_000001405.40 (GRCh38.p14) with annotation release GCF_00 · 2025-08-06 | open · NCBI data usage policy: NCBI places no restrictions on the use or dist | gene_level, region_level, complex_locus · GRCh38 | ✓ | ✓ | ✓ | ✓ | ✓ | D14, D17 | retrieve_source |

## D15 Resource catalogues and classifications

| Source | Version / date | Access · licence | Variants · assemblies | cat | acc | acq | val | map | Usable for | Decision |
|---|---|---|---|---|---|---|---|---|---|---|
| **bio.tools and OpenEBench tool registry** (`biotools`) | OpenEBench snapshot: 180,237 records; newest record timestam · 2026-10-23 | open · OpenEBench records carry CC BY 4.0 (per-record @license) | n/a · none/not_applicable | ✓ | ✓ | · | · | · | — | do_not_acquire |
| **FAIRsharing** (`fairsharing`) | Live registry · 2026-10-25 | registration_required · API terms not retrieved; the API requires registration (unverified) | n/a · none/not_applicable | ✓ | ✓ | · | · | · | — | do_not_acquire |

## D16 Background documentation and terminology

| Source | Version / date | Access · licence | Variants · assemblies | cat | acc | acq | val | map | Usable for | Decision |
|---|---|---|---|---|---|---|---|---|---|---|
| **NHGRI educational pages (Talking Glossary, 1000 Genomes Project page)** (`nhgri_education`) | Glossary entries updated 2026-10-24; 1000 Genomes page archi · 2026-10-24 | open · NHGRI website content is in the public domain unless otherwise indicat | n/a · none/not_applicable | ✓ | ✓ | ✓ | ✓ | · | — | reuse_existing |

## D17 Complex loci, paralogs and region reliability

| Source | Version / date | Access · licence | Variants · assemblies | cat | acc | acq | val | map | Usable for | Decision |
|---|---|---|---|---|---|---|---|---|---|---|
| **GIAB genome stratifications v3.6 (GRCh38): difficult and complex regions** (`giab_stratifications`) | v3.6 (directory listed 2025-04-25; GRCh38 files dated 2025-0 · 2025-04-25 | open · "Publicly released data are freely available for reuse without embargo | region_level, complex_locus · GRCh38 | ✓ | ✓ | ✓ | ✓ | ✓ | D17 | retrieve_source |
| **IPD-IMGT/HLA Database (HLA allele sequences and nomenclature)** (`ipd_imgt_hla`) | IPD-IMGT/HLA 3.66.0 (header of Allelelist.txt dated 2026-10- · 2026-10-14 | open_with_terms · Creative Commons Attribution-NoDerivs: free to copy, distribute, displ | complex_locus, region_level · none/not_applicable | ✓ | ✓ | · | · | · | — | do_not_acquire |
| **IPD-KIR Database (killer-cell immunoglobulin-like receptor allele sequences)** (`ipd_kir`) | IPD-KIR 2.15.0 (Allelelist.txt header, dated 2025-12-10) · 2025-12-10 | open_with_terms · Creative Commons Attribution-NoDerivs: free to copy, distribute, displ | complex_locus, region_level · none/not_applicable | ✓ | ✓ | · | · | · | — | do_not_acquire |

## Existing dc assets (preserved, unchanged)

See `inventory/EXISTING_ASSETS.md`. They are referenced here through `overlap_with_existing`.
