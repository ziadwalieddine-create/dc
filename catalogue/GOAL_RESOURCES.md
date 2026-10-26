# Resources named in the project goal

*Generated 2026-10-26 02:55 UTC by `catalogue/build_catalogue.py` from `catalogue/sources/*.json`, `acquisition/targets.json`, `acquisition/manifest.jsonl`, `acquisition/validation_report.json` and `acquisition/compat_report.json`. Edit the sources, not this file.*

`✓` done · `~` partial · `·` not done. Stages: **cat**alogued, **acc**essible from this environment, **acq**uired and checksum-verified, **val**idated against the source description, **map**ped to dc assets in a demonstrated join, **use**able (acquired + validated + mapped + licence verified). Stages beyond *acc* are evaluated per dimension, using only the downloaded files that support that dimension. `integrated` (a documented, tested query layer that returns records by variant, gene, protein or region) is defined in `vocab.json` and credited to no source, because no such layer exists yet.

Every resource or category that the project goal names, with the catalogue records that cover it. The state of each is computed in GOAL_RESOURCES.md from the lifecycle of its records; nothing here asserts that a resource is usable.

*Best stage* is the furthest lifecycle stage that any of the listed records reaches over all of its files; it is computed, not asserted. A resource with no record is a gap.

## Reference and annotation

| Resource | Records (best stage · decision) | Best stage | Note |
|---|---|---|---|
| Ensembl VEP and its plugins | `ensembl_vep` (usable · adopt_curated) | usable | Plugin catalogue (74 plugins) and release-116 plugin data listings are recorded; MaveDB and LOEUF plugin files are acquired; SpliceAI, SpliceVault, Enformer and the other plugin data sets are listed with sizes and not mirrored. |
| GENCODE | `gencode_mane` (usable · retrieve_source)<br>`ensembl_genomes` (usable · retrieve_source) | usable | The Ensembl release-116 GTF (the GENCODE gene set for human) is mirrored with its checksum. |
| MANE | `gencode_mane` (usable · retrieve_source) | usable | MANE v1.5 summary and GTF acquired; transcripts on alternate and patch sequences resolve through the RefSeq report and MANE Select agrees with the RefSeq annotation for 99.8% of transcripts. |
| RefSeq | `refseq` (usable · retrieve_source) | usable | Assembly report, regions, alternate-locus placement, feature counts and the GFF3 annotation of GRCh38.p14 (release RS_2025_08) acquired with published MD5s. |
| dbNSFP | `dbnsfp` (accessible · adopt_curated) | accessible | Registration is required for the data; its predictors reach the library separately (AlphaMissense, CADD scores in gnomAD INFO). |

## Functional assays

| Resource | Records (best stage · decision) | Best stage | Note |
|---|---|---|---|
| MaveDB | `mavedb` (usable · reuse_existing) | usable | Existing verified asset: 1,402 human-target score sets and their tables. |
| ProteinGym | `proteingym` (usable · adopt_curated) | usable | 217 substitution and 66 indel DMS assays acquired and checked row by row; overlap with MaveDB experiments not yet matched. |
| MAVE and deep mutational scanning data | `mavedb` (usable · reuse_existing)<br>`proteingym` (usable · adopt_curated) | usable | Two independent routes; each experiment must be counted once. |
| MPRA and reporter-assay variant effects | `mavedb` (usable · reuse_existing)<br>`mpravardb` (accessible · defer) | usable | MaveDB holds MPRA score sets; MPRAVarDB is catalogued and deferred until its licence and genome build are known. |
| Knock-in and rescue experiments | — | no record | No catalogued source: such experiments live in the literature and in model-organism allele and construct records (the Alliance bucket has ALLELE and CONSTRUCT directories, not fetched). |

## Population and association

| Resource | Records (best stage · decision) | Best stage | Note |
|---|---|---|---|
| gnomAD | `gnomad_variants` (usable · retrieve_source)<br>`gnomad_constraint` (usable · retrieve_source)<br>`gnomad_sv_cnv` (usable · retrieve_source)<br>`gnomad_ld` (accessible · defer) | usable | Allele-frequency slices at MaveDB-mapped regions, constraint metrics and the SV/CNV callsets acquired; LD matrices (GRCh37, large) deferred. |
| 1000 Genomes | `1000genomes` (usable · extend_integration) | usable | Sample metadata and release files catalogued; the project page and tutorials are the existing NHGRI asset. |
| GWAS Catalog | `gwas_catalog` (usable · adopt_curated) | usable | Associations and studies tables acquired. |
| Open Targets | `opentargets` (usable · adopt_curated) | usable | 27 of 61 datasets of release 26.09 acquired and validated. |
| Linkage disequilibrium | `ldlink` (accessible · do_not_acquire)<br>`gnomad_ld` (accessible · defer)<br>`1000genomes` (usable · extend_integration) | usable | LDlink is an on-demand service and is not acquired; gnomAD LD matrices are deferred; LD can be computed from phased 1000 Genomes genotypes. |
| Fine-mapping and colocalisation | `eqtl_catalogue` (accessible · adopt_curated)<br>`opentargets` (usable · adopt_curated) | usable | Open Targets credible sets and the eQTL Catalogue fine-mapping are the sources; the large colocalisation dataset is deferred. |

## Regulation

| Resource | Records (best stage · decision) | Best stage | Note |
|---|---|---|---|
| GTEx | `gtex` (usable · retrieve_source) | usable | v11 gene median TPM and the dataset API acquired. |
| ENCODE and candidate cis-regulatory elements | `encode` (usable · retrieve_source) | usable | SCREEN cCRE class BED files acquired. |
| eQTL Catalogue | `eqtl_catalogue` (accessible · adopt_curated) | accessible | Catalogued; summary statistics and credible sets are not acquired. |
| eQTL and sQTL | `gtex` (usable · retrieve_source)<br>`eqtl_catalogue` (accessible · adopt_curated) | usable | GTEx and eQTL Catalogue are catalogued for this; association files are not yet acquired. |
| Chromatin accessibility and binding | `encode` (usable · retrieve_source) | usable | cCRE classification only; signal tracks are not acquired. |
| Tissue and cell-type specificity | `gtex` (usable · retrieve_source)<br>`hpa` (usable · retrieve_source) | usable | GTEx median TPM and HPA consensus tissue and single-cell-type tables acquired. |

## Perturbation

| Resource | Records (best stage · decision) | Best stage | Note |
|---|---|---|---|
| DepMap | `depmap` (catalogued · adopt_curated) | catalogued | Blocked for automation (human-verification page); not circumvented. |
| IMPC | `impc` (usable · adopt_curated) | usable | Catalogued; data files are not acquired. |
| BioGRID ORCS | `biogrid_orcs` (usable · retrieve_source) | usable | Human screens acquired. |
| Natural human knockouts | `gnomad_constraint` (usable · retrieve_source) | usable | Constraint metrics and gnomAD LoF curation acquired. |
| Knockout and knockdown screens | `depmap` (catalogued · adopt_curated)<br>`biogrid_orcs` (usable · retrieve_source)<br>`impc` (usable · adopt_curated)<br>`alliance_genome` (accessible · defer) | usable | BioGRID ORCS acquired; others catalogued. |
| CRISPR interference and activation | `biogrid_orcs` (usable · retrieve_source) | usable | BioGRID ORCS records screen types; the share of interference and activation screens was not counted. |
| Overexpression and gene dosage | `clingen` (usable · retrieve_source)<br>`biogrid_interactions` (usable · retrieve_source) | usable | ClinGen dosage curation and BioGRID dosage-rescue and dosage-lethality interactions acquired. |

## Clinical and pharmacogenomics

| Resource | Records (best stage · decision) | Best stage | Note |
|---|---|---|---|
| ClinVar | `clinvar` (usable · retrieve_source) | usable | GRCh38 VCF acquired; records on RefSeq NT_/NW_ contigs resolve through the assembly report. |
| ClinGen | `clingen` (usable · retrieve_source) | usable | Dosage, recurrent-region and gene-validity files acquired. |
| PharmGKB (ClinPGx) | `pharmgkb` (usable · retrieve_source) | usable | Annotation archives acquired. |
| PharmVar | `pharmvar` (usable · retrieve_source) | usable | Complete database 6.2.29 acquired (CC BY-SA 4.0, research use). |
| CPIC | `cpic` (usable · retrieve_source) | usable | Gene and guideline lists and eight database API tables (allele definitions, locations, function, phenotypes, diplotypes, recommendations) acquired. |
| PharmCAT | `pharmcat` (accessible · do_not_acquire) | accessible | A tool applying the CPIC and ClinPGx tables; catalogued, not acquired. |

## Proteins and pathways

| Resource | Records (best stage · decision) | Best stage | Note |
|---|---|---|---|
| UniProt | `uniprot` (usable · retrieve_source) | usable | Reviewed human entries and the humsavar variant index acquired. |
| InterPro | `interpro` (usable · retrieve_source) | usable | Per-protein annotation of MaveDB proteins acquired. |
| AlphaFold DB | `alphafold_db` (usable · retrieve_source) | usable | Canonical models of MaveDB proteins acquired. |
| Reactome | `reactome` (usable · retrieve_source) | usable | Pathway and mapping tables acquired. |
| Interaction databases | `biogrid_interactions` (usable · retrieve_source)<br>`intact` (accessible · defer)<br>`string_db` (accessible · defer)<br>`opentargets` (usable · adopt_curated) | usable | BioGRID (physical and genetic) acquired; IntAct, STRING and Reactome interactions reach the library through the Open Targets interaction dataset. |
| PDB residue mappings | `pdbe_sifts` (usable · retrieve_source) | usable | SIFTS mappings acquired. |

## Complex variants

| Resource | Records (best stage · decision) | Best stage | Note |
|---|---|---|---|
| Copy-number and structural variation | `gnomad_sv_cnv` (usable · retrieve_source)<br>`dgv` (usable · retrieve_source)<br>`dbvar` (accessible · defer)<br>`clingen` (usable · retrieve_source)<br>`hgsvc3` (usable · retrieve_source)<br>`1000genomes` (usable · extend_integration)<br>`decipher` (catalogued · do_not_acquire) | usable | gnomAD-SV/CNV, DGV, ClinGen regions and HGSVC3 acquired; dbVar deferred; DECIPHER is restricted. |
| Short tandem repeats, VNTRs and repeat expansions | `strchive` (usable · retrieve_source)<br>`trexplorer` (accessible · defer) | usable | STRchive disease loci acquired; the genome-wide TRExplorer catalogue is deferred until its licence is known. |
| Complex pharmacogene loci | `pharmvar` (usable · retrieve_source)<br>`cpic` (usable · retrieve_source)<br>`ipd_imgt_hla` (accessible · do_not_acquire)<br>`ipd_kir` (accessible · do_not_acquire)<br>`pharmcat` (accessible · do_not_acquire) | usable | PharmVar and CPIC acquired; HLA and KIR allele databases are not mirrored (no-derivatives licence, mirroring discouraged). |
| Callability and region reliability | `giab_stratifications` (usable · retrieve_source)<br>`refseq` (usable · retrieve_source)<br>`hgsvc3` (usable · retrieve_source) | usable | GIAB difficult-region stratifications, alternate-locus placement and HGSVC3 repeat flags acquired. |
| Mitochondrial DNA variants | `mitomap` (catalogued · defer)<br>`cpic` (usable · retrieve_source) | usable | MITOMAP refuses automated access; CPIC supplies MT-RNR1 alleles only. |

## Discovery

| Resource | Records (best stage · decision) | Best stage | Note |
|---|---|---|---|
| FAIRsharing | `fairsharing` (accessible · do_not_acquire) | accessible | Registration required; not acquired. |
| bio.tools and OpenEBench | `biotools` (accessible · do_not_acquire) | accessible | API refused automated clients; OpenEBench counts recorded. |

48 resources named; 41 have at least one record that is usable for something (acquired, validated, mapped and with a verified licence).

## Not yet usable

| Group | Resource | Best stage reached |
|---|---|---|
| Reference and annotation | dbNSFP | accessible |
| Functional assays | Knock-in and rescue experiments | no record |
| Regulation | eQTL Catalogue | accessible |
| Perturbation | DepMap | catalogued |
| Clinical and pharmacogenomics | PharmCAT | accessible |
| Discovery | FAIRsharing | accessible |
| Discovery | bio.tools and OpenEBench | accessible |
