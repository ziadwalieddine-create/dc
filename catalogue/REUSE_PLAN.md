# Source-by-source reuse, acquisition and integration plan

*Generated 2026-10-25 12:30 UTC by `catalogue/build_catalogue.py` from `catalogue/sources/*.json`, `acquisition/targets.json`, `acquisition/manifest.jsonl`, `acquisition/validation_report.json` and `acquisition/compat_report.json`. Edit the sources, not this file.*

Decision hierarchy: reuse existing verified material → adopt existing curated datasets → extend existing integrations → retrieve uncovered source data → build original methods only where justified.

## Reuse what is already verified (3)

### `ensembl_genomes`: Ensembl genome registry, assemblies and gene annotation (human)
- **Why:** Registry already collected and verified; fetch FASTA/GTF only when an analysis needs them, preferring GENCODE/MANE files already staged for acquisition.
- **Status (all files):** catalogued✓ accessible✓ acquired✓ validated✓ mapped✓ · usable for: D14 · existing assets ensembl-human-genomes: downloaded_verified
- **Access:** open (Ensembl imposes no restrictions on access to or use of the data (disclaimer page))

### `mavedb`: MaveDB (Multiplexed Assays of Variant Effect database)
- **Why:** Metadata and 1,980 score/count tables are already downloaded and verified (scrape/mavedb-human, data/mavedb-scores). The Ensembl coordinate-mapped file and the Zenodo v4 dump are being compared against it before any change.
- **Status (all files):** catalogued✓ accessible✓ acquired✓ validated✓ mapped✓ · usable for: D01 · existing assets mavedb-human: downloaded_verified
- **Access:** open (Per record: CC0 (1,340), CC BY 4.0 (45), CC BY-SA 4.0 (11), CC BY-NC-SA 4.0 (6) among our 1,402 human sets; Zenodo bulk )
- **Targets:** `mavedb_zenodo_dump_v4`, `ensembl_vep_mavedb_variants`, `ensembl_vep_mavedb_variants_tbi`

### `nhgri_education`: NHGRI educational pages (Talking Glossary, 1000 Genomes Project page)
- **Why:** Already collected, verified and preserved under scrape/. No further acquisition.
- **Status (all files):** catalogued✓ accessible✓ acquired✓ validated✓ mapped· · usable for: nothing yet · existing assets nhgri-gwas, nhgri-1000-genomes: downloaded_verified
- **Access:** open (NHGRI website content is in the public domain unless otherwise indicated (NHGRI Copyright Policy); material with a credi)

## Adopt an existing curated integration (7)

### `dbnsfp`: dbNSFP
- **Why:** Best single consolidated source of missense predictors and frequencies. Not acquired: it needs the user to register with an institutional email.
- **Status (all files):** catalogued✓ accessible✓ acquired· validated· mapped· · usable for: nothing yet
- **Access:** registration_required (Academic branch CC BY-NC-ND 4.0 (about 50 GB); commercial branch under paid licence)

### `depmap`: DepMap (Cancer Dependency Map)
- **Why:** Open Targets 26.09 carries a target_essentiality (DepMap) dataset and CRISPR screen evidence; use that. Direct download needs the user to fetch from the DepMap downloads page manually.
- **Status (all files):** catalogued✓ accessible· acquired· validated· mapped· · usable for: nothing yet
- **Access:** blocked_for_automation (Terms page not readable (verification wall); commonly CC BY 4.0)
- **Unverified:** version, licence, update_frequency

### `ensembl_vep`: Ensembl VEP and release-116 plugin data
- **Why:** Ensembl already integrates MaveDB, LOEUF, SpliceAI, SpliceVault, ProtVar, mutfunc and Enformer for VEP; reuse these rather than re-deriving them. The MaveDB file is acquired first to measure overlap with the API-built set.
- **Status (all files):** catalogued✓ accessible✓ acquired✓ validated✓ mapped✓ · usable for: D01
- **Access:** open_with_terms (Ensembl: no restrictions on data (verified); plugin resources carry upstream terms (e.g. gnomAD notes SpliceAI annotatio)
- **Targets:** `ensembl_vep_mavedb_variants`, `ensembl_vep_mavedb_variants_tbi`
- **Unverified:** complete plugin list, per-plugin licences

### `eqtl_catalogue`: eQTL Catalogue
- **Why:** Open Targets already consumes eQTL Catalogue credible sets for colocalisation; defer direct download of the very large sumstats until a locus-level analysis requires them.
- **Status (all files):** catalogued✓ accessible✓ acquired· validated· mapped· · usable for: nothing yet
- **Access:** open (Data under Creative Commons Attribution 4.0 International; code under Apache 2.0 (eQTL Catalogue License page))

### `gwas_catalog`: NHGRI-EBI GWAS Catalog
- **Why:** Open Targets integrates GWAS Catalog studies into credible sets; the flat associations table is small and fetched for direct lookups. The flat studies.tsv was fetched first and discarded: it has no study accession column; the associations file carries STUDY ACCESSION.
- **Status (all files):** catalogued✓ accessible✓ acquired✓ validated✓ mapped✓ · usable for: D06
- **Access:** open (EMBL-EBI general terms of use; summary statistics CC0 unless otherwise stated)
- **Targets:** `gwas_catalog_associations`
- **Unverified:** assembly

### `impc`: IMPC (International Mouse Phenotyping Consortium)
- **Why:** Open Targets integrates IMPC as evidence_impc and mouse_phenotype with human ortholog mapping; fetch source tables only if gene-level mouse detail is needed.
- **Status (all files):** catalogued✓ accessible✓ acquired· validated· mapped· · usable for: nothing yet
- **Access:** open (Creative Commons Attribution 4.0 (site footer statement))

### `opentargets`: Open Targets Platform
- **Why:** Already integrates GWAS credible sets, colocalisation, DepMap, CRISPR screens, IMPC, ClinGen, ClinVar (EVA) and pharmacogenomics. Use it as the primary route for those dimensions instead of scraping each source; keep the originals only where Open Targets omits information.
- **Status (all files):** catalogued✓ accessible✓ acquired✓ validated✓ mapped✓ · usable for: D06, D04, D10, D11, D05, D07, D03, D14
- **Access:** open (CC0 1.0 for Platform data (docs page); a separate data-source licensing list exists on the same page)
- **Targets:** `opentargets_26_09_biosample`, `opentargets_26_09_credible_set`, `opentargets_26_09_disease`, `opentargets_26_09_disease_hpo`, `opentargets_26_09_disease_phenotype`, `opentargets_26_09_enhancer_to_gene`, `opentargets_26_09_evidence_crispr`, `opentargets_26_09_evidence_crispr_screen`, `opentargets_26_09_evidence_gene2phenotype`, `opentargets_26_09_evidence_gene_burden`, `opentargets_26_09_evidence_genomics_england`, `opentargets_26_09_evidence_impc`, `opentargets_26_09_evidence_orphanet`, `opentargets_26_09_evidence_uniprot_literature`, `opentargets_26_09_evidence_uniprot_variants`, `opentargets_26_09_go`, `opentargets_26_09_homology`, `opentargets_26_09_interaction`, `opentargets_26_09_interaction_evidence`, `opentargets_26_09_l2g_prediction`, `opentargets_26_09_mouse_phenotype`, `opentargets_26_09_pharmacogenomics`, `opentargets_26_09_so`, `opentargets_26_09_study`, `opentargets_26_09_target`, `opentargets_26_09_target_essentiality`, `opentargets_26_09_variant`
- **Deferred (15.43 GiB):** colocalisation: Largest dataset (200 Parquet parts). Fetch when colocalisation of GWAS and molecular-QTL signals is needed for specific loci; it has no coordinate index, so it is read by full scan.
- **Deferred (1.58 GiB):** baseline_expression: Tissue and cell-type expression; GTEx median TPM is already acquired directly.
- **Deferred (9.69 GiB):** association_* (6 datasets: by datasource/datatype/overall, direct and indirect): Derived target-disease association scores; recomputable from the evidence datasets and not variant-level.
- **Deferred (7.48 GiB):** evidence_europepmc, literature_entity_lut, literature_vector: Text-mining evidence; weak provenance for variant-level interpretation.
- **Deferred (0.446 GiB):** evidence_eva, evidence_eva_somatic, evidence_clingen, evidence_gwas_credible_sets, evidence_reactome, reactome, evidence_expression_atlas: Re-expressions of ClinVar, ClinGen, GWAS credible sets and Reactome, which are acquired directly (credible_set is acquired in full).
- **Deferred (0.128 GiB):** cancer and drug datasets (evidence_cancer_biomarkers, evidence_cancer_gene_census, evidence_intogen, evidence_clinical_precedence, drug_molecule, drug_mechanism_of_action, drug_warning, clinical_indication, clinical_report, clinical_target, chemical_probe, target_tractability, target_safety_event, target_prioritisation, openfda_significant_adverse_drug_reactions): Somatic-cancer and drug-development content; outside germline genome interpretation for now.
- **Deferred (0.028 GiB):** transcript: Ensembl transcript table (29 MiB); the MANE summary already bridges Ensembl, RefSeq and HGNC identifiers for the joins in use. Fetch when transcript-level consequences are needed.

## Extend an existing integration (1)

### `1000genomes`: 1000 Genomes Project / IGSR (International Genome Sample Resource)
- **Why:** The existing NHGRI page gives project documentation only. Acquire the 3,202-sample population table now (links genotypes to ancestry); fetch chromosome-scale genotype files only on demand, since gnomAD/dbNSFP already carry 1000 Genomes frequencies.
- **Status (all files):** catalogued✓ accessible✓ acquired✓ validated✓ mapped· · usable for: nothing yet
- **Access:** open (IGSR provides open data; terms of use are in the IGSR data disclaimer (not parsed))
- **Targets:** `igsr_1000g_samples`

## Retrieve the source data (18)

### `alphafold_db`: AlphaFold Protein Structure Database
- **Why:** Per-protein retrieval for the UniProt proteins our MaveDB score sets are about (models, per-residue pLDDT and prediction metadata through the API); the 5.18 GB human proteome bundle is not needed for that set.
- **Status (all files):** catalogued✓ accessible✓ acquired✓ validated✓ mapped✓ · usable for: D02
- **Access:** open (CC BY 4.0, academic and commercial use (FAQ))
- **Targets:** `alphafold_db_mavedb_proteins`

### `alphamissense`: AlphaMissense
- **Why:** Also embedded in dbNSFP, which needs registration; the open file lets us test predictor integration now. Both the gene-level table and the 643 MB variant table (71M canonical-transcript missense SNVs) are fetched, checksum-verified and structurally validated. Licence CC BY-NC-SA 4.0 (non-commercial, share-alike), so the data stay in the git-ignored data/ area.
- **Status (all files):** catalogued✓ accessible✓ acquired✓ validated✓ mapped✓ · usable for: D13
- **Access:** open_with_terms (CC BY-NC-SA 4.0 (Zenodo record))
- **Targets:** `alphamissense_gene_hg38`, `alphamissense_hg38`

### `biogrid_orcs`: BioGRID ORCS (Open Repository of CRISPR Screens)
- **Why:** Open Targets lists evidence_crispr/evidence_crispr_screen, but its source mix is not verified to equal ORCS; fetch the human bundle for gene-level screen detail.
- **Status (all files):** catalogued✓ accessible✓ acquired✓ validated✓ mapped✓ · usable for: D04
- **Access:** open (MIT License (stated on the BioGRID download repository page for BioGRID-ORCS; data provided without warranty))
- **Targets:** `biogrid_orcs_human`

### `clingen`: ClinGen (dosage sensitivity, gene-disease validity, region curation)
- **Why:** Small, openly licensed and directly addresses dosage sensitivity and CNV interpretation; Open Targets carries only evidence_clingen gene-disease links.
- **Status (all files):** catalogued✓ accessible✓ acquired✓ validated✓ mapped✓ · usable for: D05, D11
- **Access:** open (ClinGen terms: curated content released openly, free to capture and use)
- **Targets:** `clingen_gene_curation_list`, `clingen_region_curation_list`, `clingen_haploinsufficiency_gene`, `clingen_triplosensitivity_gene`, `clingen_recurrent_cnv_v2_1`, `clingen_gene_validity`

### `clinvar`: ClinVar
- **Why:** Source of record for clinical classifications; Open Targets only carries a subset (EVA evidence), so the VCF is retrieved directly.
- **Status (all files):** catalogued✓ accessible✓ acquired✓ validated✓ mapped✓ · usable for: D11
- **Access:** open (NCBI public domain policy; submitters may claim rights in some content (NCBI policies page))
- **Targets:** `clinvar_vcf_grch38`, `clinvar_readme`, `clinvar_readme_vcf`

### `cpic`: CPIC (Clinical Pharmacogenetics Implementation Consortium)
- **Why:** Small API tables; combine with ClinPGx and PharmVar.
- **Status (all files):** catalogued✓ accessible✓ acquired✓ validated✓ mapped✓ · usable for: D07
- **Access:** open_with_terms (The CPIC licence URL resolves to the ClinPGx data-use text (CC BY-SA 4.0))
- **Targets:** `cpic_genes`, `cpic_guidelines`
- **Unverified:** CPIC-specific licence terms

### `dgv`: DGV (Database of Genomic Variants)
- **Why:** Complements gnomAD-SV with older, study-level control CNVs and structural variants. The 364 MB GRCh38 release (2025-12-01) is fetched, size-verified (the site publishes no checksum) and structurally validated; used under the TCAG disclaimer with acknowledgement (cite MacDonald et al. 2013) and not redistributed.
- **Status (all files):** catalogued✓ accessible✓ acquired✓ validated✓ mapped✓ · usable for: D08
- **Access:** open_with_terms (TCAG Notice and Disclaimer of Database Use (data provided 'AS IS'; appropriate acknowledgement required; cite MacDonald )
- **Targets:** `dgv_grch38_variants_2025_12_01`

### `encode`: ENCODE and SCREEN candidate cis-regulatory elements
- **Why:** SCREEN already integrates ENCODE and partner data into a registry; use the class BED files (about 129 MB for non-overlapping classes) instead of parsing thousands of experiments.
- **Status (all files):** catalogued✓ accessible✓ acquired✓ validated✓ mapped✓ · usable for: D03
- **Access:** open (ENCODE Data Use Policy: 'all data produced will be available for unrestricted use immediately upon release to public dat)
- **Targets:** `screen_ccre_pls`, `screen_ccre_pels`, `screen_ccre_dels`, `screen_ccre_ca_ctcf`, `screen_ccre_ca_h3k4me3`, `screen_ccre_ca_tf`, `screen_ccre_ca`, `screen_ccre_tf`

### `gencode_mane`: GENCODE and MANE transcript sets
- **Why:** MANE files are small and fix the transcript choice needed to interpret HGVS across ClinVar, MaveDB and gnomAD; the GENCODE genome FASTA (806 to 856 MB) is deferred.
- **Status (all files):** catalogued✓ accessible✓ acquired✓ validated✓ mapped✓ · usable for: D14
- **Access:** open (GENCODE: 'all the GENCODE project data is open access' (GENCODE data-access page), under the EMBL-EBI Terms of Use (free)
- **Targets:** `mane_summary_v1_5`, `mane_gtf_v1_5`

### `gnomad_constraint`: gnomAD constraint metrics (LoF/missense) and pext
- **Why:** Primary source of record and moderate size; the Ensembl LOEUF file is a derivative of the same resource so is not downloaded separately.
- **Status (all files):** catalogued✓ accessible✓ acquired✓ validated✓ mapped✓ · usable for: D05
- **Access:** open (CC0 for primary data (policies page))
- **Targets:** `gnomad_constraint_v4_1_2`

### `gnomad_sv_cnv`: gnomAD structural variants and copy-number variants
- **Why:** Both VCFs are fetched and validated: the exome CNV VCF (8.9 MB) and the genome SV sites VCF (1.74 GB; size from the S3 listing, MD5 from the downloads page). They complement the older control-cohort DGV catalogue and the ClinGen dosage regions.
- **Status (all files):** catalogued✓ accessible✓ acquired✓ validated✓ mapped✓ · usable for: D08
- **Access:** open (CC0 for primary data (policies page))
- **Targets:** `gnomad_exome_cnv_v4_1`, `gnomad_sv_v4_1_sites_vcf`

### `gnomad_variants`: gnomAD v4.1 (exomes and genomes) allele frequencies
- **Why:** Whole-chromosome VCFs (11-67 GiB each) are not mirrored. The files are bgzip + tabix on a public bucket, so regions are read by HTTP range requests: slices for the genomic positions that our assets name (MaveDB-mapped SNVs) are fetched and kept; the same primitive serves later lookups for variants of interest. Absence of a record in a slice means not observed or filtered out, not a frequency of zero.
- **Status (all files):** catalogued✓ accessible✓ acquired✓ validated✓ mapped✓ · usable for: D12, D13
- **Access:** open (Primary exome/genome data released under CC0; some annotations carry separate terms (SpliceAI is CC BY-NC 4.0))
- **Targets:** `gnomad_v4_1_exomes_regional`, `gnomad_v4_1_genomes_regional`, `gnomad_v4_1_joint_regional`

### `gtex`: GTEx (Genotype-Tissue Expression)
- **Why:** Median expression table is small and fetched now. Large eQTL tables are available but Open Targets and eQTL Catalogue provide fine-mapped, harmonised versions; pull only what is needed.
- **Status (all files):** catalogued✓ accessible✓ acquired✓ validated✓ mapped✓ · usable for: D03
- **Access:** open (Open-access files free to download; raw sequence and full donor metadata are protected (dbGaP/AnVIL))
- **Targets:** `gtex_v11_gene_median_tpm`, `gtex_datasets_api`

### `interpro`: InterPro
- **Why:** Per-protein retrieval of InterPro and member-database matches (with residue locations) for the UniProt proteins our MaveDB score sets are about; the bulk match files are not needed for that set.
- **Status (all files):** catalogued✓ accessible✓ acquired✓ validated✓ mapped✓ · usable for: D02
- **Access:** open (CC0 1.0 Universal for InterPro, Pfam, PRINTS and SFLD downloadable data; InterProScan under the Apache License (InterPro)
- **Targets:** `interpro_mavedb_proteins`

### `pharmgkb`: ClinPGx / PharmGKB
- **Why:** Small archives; Open Targets only carries a pharmacogenomics subset.
- **Status (all files):** catalogued✓ accessible✓ acquired✓ validated✓ mapped✓ · usable for: D07
- **Access:** open_with_terms (CC BY-SA 4.0 (data usage policy))
- **Targets:** `clinpgx_variantAnnotations`, `clinpgx_relationships`, `clinpgx_summaryAnnotations`, `clinpgx_guidelineAnnotations`, `clinpgx_drugLabels`

### `reactome`: Reactome
- **Why:** Small mapping files; Open Targets also carries a reactome dataset.
- **Status (all files):** catalogued✓ accessible✓ acquired✓ validated✓ mapped✓ · usable for: D10
- **Access:** open (Data licensed CC0; illustrations CC BY 4.0 (license page))
- **Targets:** `reactome_ensembl2reactome`, `reactome_uniprot2reactome_all_levels`, `reactome_pathways`

### `strchive`: STRchive and related tandem-repeat catalogues
- **Why:** Small, openly licensed disease-locus catalogue. Catalogue files are linked for download from the STRchive site itself (hosted as raw files on GitHub).
- **Status (all files):** catalogued✓ accessible✓ acquired✓ validated✓ mapped✓ · usable for: D09
- **Access:** open (CC BY 4.0 (site footer); research use only)
- **Targets:** `strchive_hg38_general`, `strchive_hg38_trgt`, `strchive_hg38_stranger_json`

### `uniprot`: UniProtKB (Swiss-Prot and TrEMBL)
- **Why:** Primary source; fetch the reviewed human table with identifier cross-references (needed to join protein-level MaveDB targets to genes and structures) and the humsavar missense variant index (a protein-level clinical annotation that joins to MaveDB variants through UniProt accession and amino-acid change).
- **Status (all files):** catalogued✓ accessible✓ acquired✓ validated✓ mapped✓ · usable for: D02, D11
- **Access:** open (CC BY 4.0 for copyrightable parts)
- **Targets:** `uniprot_human_reviewed`, `uniprot_humsavar`

## Deferred (2)

### `dbvar`: dbVar (NCBI database of genomic structural variation)
- **Why:** Catalogued only: ClinGen, gnomAD-SV and DGV cover the near-term needs; inspect dbVar contents and currency before any download.
- **Status (all files):** catalogued✓ accessible✓ acquired· validated· mapped· · usable for: nothing yet
- **Access:** open (NCBI data usage policy: NCBI places no restrictions on the use or distribution of the data in its databases; original su)
- **Unverified:** currency, contents

### `pharmvar`: PharmVar (Pharmacogene Variation Consortium)
- **Why:** Licence and bulk download route unverified. Use CPIC/ClinPGx tables meanwhile; revisit after reading the terms.
- **Status (all files):** catalogued✓ accessible✓ acquired· validated· mapped· · usable for: nothing yet
- **Access:** unverified (Terms and Conditions page not retrievable (HTTP 404 on the URLs tried))
- **Unverified:** licence, release date

## Not acquired by design (3)

### `biotools`: bio.tools and OpenEBench tool registry
- **Why:** Used for discovery only. Counts derived by catalogue/derive/openebench_counts.py from a 253 MB snapshot that is not stored.
- **Status (all files):** catalogued✓ accessible✓ acquired· validated· mapped· · usable for: nothing yet
- **Access:** open (OpenEBench records carry CC BY 4.0 (per-record @license))

### `decipher`: DECIPHER
- **Why:** Access-controlled and denied automated access; not circumvented. User would have to apply for access.
- **Status (all files):** catalogued✓ accessible· acquired· validated· mapped· · usable for: nothing yet
- **Access:** restricted_agreement (Automated access to the downloads page returned {"error":"Access Denied"}; data access normally requires an agreement)
- **Unverified:** version, licence, everything except the access denial

### `fairsharing`: FAIRsharing
- **Why:** A discovery and classification aid, not an evidence source. Use it to cross-check this catalogue for missing resources; programmatic harvesting would need the user to create an account.
- **Status (all files):** catalogued✓ accessible✓ acquired· validated· mapped· · usable for: nothing yet
- **Access:** registration_required (API terms not retrieved; the API requires registration)
- **Unverified:** API access terms, licence
