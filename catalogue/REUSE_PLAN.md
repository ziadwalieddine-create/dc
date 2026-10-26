# Source-by-source reuse, acquisition and integration plan

*Generated 2026-10-26 15:25 UTC by `catalogue/build_catalogue.py` from `catalogue/sources/*.json`, `acquisition/targets.json`, `acquisition/manifest.jsonl`, `acquisition/validation_report.json` and `acquisition/compat_report.json`. Edit the sources, not this file.*

Decision hierarchy: reuse existing verified material → adopt existing curated datasets → extend existing integrations → retrieve uncovered source data → build original methods only where justified.

## Reuse what is already verified (2)

### `mavedb`: MaveDB (Multiplexed Assays of Variant Effect database)
- **Why:** Metadata and 1,980 score/count tables are already downloaded and verified (scrape/mavedb-human, data/mavedb-scores). The Ensembl coordinate-mapped file and the Zenodo v4 dump were compared against them (acquisition/compat_tests.py) without any change to the existing assets; score tables of the sets that carry a post-mapped genomic HGVS are extracted from the dump byte for byte.
- **Status (all files):** catalogued✓ accessible✓ acquired✓ validated✓ mapped✓ · usable for: D01 · existing assets mavedb-human: downloaded_verified
- **Access:** open (Per record: CC0 (1,340), CC BY 4.0 (45), CC BY-SA 4.0 (11), CC BY-NC-SA 4.0 (6) among our 1,402 human sets; Zenodo bulk )
- **Targets:** `mavedb_zenodo_dump_v4`, `ensembl_vep_mavedb_variants`, `ensembl_vep_mavedb_variants_tbi`, `mavedb_scores_genomic_mapped`

### `nhgri_education`: NHGRI educational pages (Talking Glossary, 1000 Genomes Project page)
- **Why:** Already collected, verified and preserved under scrape/. No further acquisition.
- **Status (all files):** catalogued✓ accessible✓ acquired✓ validated✓ mapped· · usable for: nothing yet · existing assets nhgri-gwas, nhgri-1000-genomes: downloaded_verified
- **Access:** open (NHGRI website content is in the public domain unless otherwise indicated (NHGRI Copyright Policy); material with a credi)

## Adopt an existing curated integration (8)

### `dbnsfp`: dbNSFP
- **Why:** Best single consolidated source of missense predictors and frequencies. Not acquired: it needs the user to register with an institutional email.
- **Status (all files):** catalogued✓ accessible✓ acquired· validated· mapped· · usable for: nothing yet
- **Access:** registration_required (Academic branch CC BY-NC-ND 4.0 (about 50 GB); commercial branch under paid licence)

### `depmap`: DepMap (Cancer Dependency Map)
- **Why:** Open Targets 26.09 carries a target_essentiality (DepMap) dataset and CRISPR screen evidence; use that. Direct download needs the user to fetch from the DepMap downloads page manually.
- **Status (all files):** catalogued✓ accessible· acquired· validated· mapped· · usable for: nothing yet
- **Access:** blocked_for_automation (Terms page not readable (verification wall); commonly CC BY 4.0)
- **Unverified:** licence, scope details (RNAi screens, expression and mutation context), mapping_methods (join identifiers), quantitative (score definition, units, direction, uncertainty), update_frequency, version

### `ensembl_vep`: Ensembl VEP and release-116 plugin data
- **Why:** Ensembl already integrates MaveDB, LOEUF, SpliceAI, SpliceVault, ProtVar, mutfunc and Enformer for VEP; reuse these rather than re-deriving them. The MaveDB file is acquired first to measure overlap with the API-built set.
- **Status (all files):** catalogued✓ accessible✓ acquired✓ validated✓ mapped✓ · usable for: D01
- **Access:** open_with_terms (Ensembl: no restrictions on data (verified); plugin resources carry upstream terms (e.g. gnomAD notes SpliceAI annotatio)
- **Targets:** `ensembl_vep_mavedb_variants`, `ensembl_vep_mavedb_variants_tbi`
- **Unverified:** per-plugin data licences not stated on the plugin page

### `eqtl_catalogue`: eQTL Catalogue
- **Why:** Open Targets lists the eQTL Catalogue among its integrated data sources, and the acquired Open Targets credible_set dataset already holds 1,718,748 eQTL, 844,764 sQTL, 521,558 transcript-usage QTL and 123,630 single-cell eQTL credible-set rows with effect sizes. Which of its projects are eQTL Catalogue datasets was not matched, so the catalogue's own files (summary statistics and credible sets per dataset, now under a new data-access scheme with an r8_beta release) are not fetched until a gap in the Open Targets coverage is shown.
- **Status (all files):** catalogued✓ accessible✓ acquired· validated· mapped· · usable for: nothing yet
- **Access:** open (Data under Creative Commons Attribution 4.0 International; code under Apache 2.0 (eQTL Catalogue License page))

### `gwas_catalog`: NHGRI-EBI GWAS Catalog
- **Why:** Open Targets integrates GWAS Catalog studies into credible sets; the flat associations table is small and fetched for direct lookups. The flat studies.tsv was fetched first and discarded: it has no study accession column; the associations file carries STUDY ACCESSION. The studies and ancestry downloads of the same release, which do carry the accession, are retrieved for study and sample metadata, together with the trait to EFO mappings.
- **Status (all files):** catalogued✓ accessible✓ acquired✓ validated✓ mapped✓ · usable for: D06
- **Access:** open (EMBL-EBI general terms of use; summary statistics CC0 unless otherwise stated)
- **Targets:** `gwas_catalog_associations`, `gwas_catalog_studies`, `gwas_catalog_ancestries`, `gwas_catalog_efo_trait_mappings`
- **Unverified:** assembly

### `impc`: IMPC (International Mouse Phenotyping Consortium)
- **Why:** Open Targets integrates IMPC with human ortholog mapping (evidence_impc, mouse_phenotype), which is the route for orthology-mapped use; the three small result files are fetched in addition because they carry the allele, zygosity, sex and viability information that the integrated datasets drop, and viability is the mouse counterpart of human loss-of-function intolerance (tested against gnomAD LOEUF).
- **Status (all files):** catalogued✓ accessible✓ acquired✓ validated✓ mapped✓ · usable for: D04, D05
- **Access:** open (Creative Commons Attribution 4.0 (site footer statement))
- **Targets:** `impc_r24_readme`, `impc_r24_genotype_phenotype_assertions_all`, `impc_r24_gene_mp_association_1`, `impc_r24_viability`

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

### `proteingym`: ProteinGym v1.3: deep mutational scanning benchmark assays (substitutions and indels)
- **Why:** An established curated benchmark that already puts 283 published assays into one checked format; reusing it avoids re-curating each paper. Some assays are expected to be the same experiments as MaveDB score sets, so they are matched by publication and target before any combination, and each experiment is counted once.
- **Status (all files):** catalogued✓ accessible✓ acquired✓ validated✓ mapped✓ · usable for: D01
- **Access:** open_with_terms (MIT licence for the ProteinGym project (LICENSE of the OATML-Markslab/ProteinGym repository, read in the read-only clone)
- **Targets:** `proteingym_v1_3_dms_substitutions`, `proteingym_v1_3_dms_indels`, `proteingym_v1_3_dms_substitutions_reference`, `proteingym_v1_3_dms_indels_reference`
- **Unverified:** licences of the underlying published assay datasets, which assays are the same experiments as MaveDB score sets (not yet matched)

## Extend an existing integration (1)

### `1000genomes`: 1000 Genomes Project / IGSR (International Genome Sample Resource)
- **Why:** The existing NHGRI page gives project documentation only. Acquire the 3,202-sample population table now (links genotypes to ancestry); fetch chromosome-scale genotype files only on demand, since gnomAD/dbNSFP already carry 1000 Genomes frequencies.
- **Status (all files):** catalogued✓ accessible✓ acquired✓ validated✓ mapped✓ · usable for: nothing yet
- **Access:** open (IGSR provides open data; terms of use are in the IGSR data disclaimer (not parsed))
- **Targets:** `igsr_1000g_samples`

## Retrieve the source data (27)

### `alphafold_db`: AlphaFold Protein Structure Database
- **Why:** Per-protein retrieval for the UniProt proteins our MaveDB score sets are about (models, per-residue pLDDT and prediction metadata through the API); the 5.18 GB human proteome bundle is not needed for that set.
- **Status (all files):** catalogued✓ accessible✓ acquired✓ validated✓ mapped✓ · usable for: D02
- **Access:** open (CC BY 4.0, academic and commercial use (FAQ))
- **Targets:** `alphafold_db_mavedb_proteins`

### `alphamissense`: AlphaMissense
- **Why:** Also embedded in dbNSFP, which needs registration; the open file lets us test predictor integration now. Both the gene-level table and the 643 MB variant table (71M canonical-transcript missense SNVs) are fetched, checksum-verified and structurally validated. Licence CC BY-NC-SA 4.0 (non-commercial, share-alike), so the data stay in the git-ignored data/ area.
- **Status (all files):** catalogued✓ accessible✓ acquired✓ validated✓ mapped✓ · usable for: D13
- **Access:** open_with_terms (Conflicting statements, applied conservatively as CC BY-NC-SA 4.0 for the copy held. The bucket README.pdf (2024-03-13) )
- **Targets:** `alphamissense_gene_hg38`, `alphamissense_hg38`
- **Unverified:** training and calibration details (paper not retrievable), which licence statement governs the 2023 data files

### `biogrid_interactions`: BioGRID protein and genetic interactions (release 5.0.262, human file)
- **Why:** The only open resource in the catalogue that carries curated genetic interactions (synthetic lethality, dosage rescue, suppression) for epistasis and compensation, alongside physical interactions; MIT-licensed and a single 189 MB archive. Open Targets already integrates IntAct, STRING, Reactome and Signor interactions, so those are not mirrored separately.
- **Status (all files):** catalogued✓ accessible✓ acquired✓ validated✓ mapped✓ · usable for: D10
- **Access:** open (MIT License: "100% freely available to both academic and commercial users" and provided without warranty (BioGRID downlo)
- **Targets:** `biogrid_interactions_organism_tab3`
- **Unverified:** update frequency (monthly releases are stated on the BioGRID site, not retrieved here)

### `biogrid_orcs`: BioGRID ORCS (Open Repository of CRISPR Screens)
- **Why:** Open Targets lists evidence_crispr/evidence_crispr_screen, but its source mix is not verified to equal ORCS; fetch the human bundle for gene-level screen detail.
- **Status (all files):** catalogued✓ accessible✓ acquired✓ validated✓ mapped✓ · usable for: D04
- **Access:** open (MIT License (stated on the BioGRID download repository page for BioGRID-ORCS; data provided without warranty))
- **Targets:** `biogrid_orcs_human`

### `clingen`: ClinGen (dosage sensitivity, gene-disease validity, region curation)
- **Why:** Small, openly licensed and directly addresses dosage sensitivity and CNV interpretation; Open Targets carries only evidence_clingen gene-disease links.
- **Status (all files):** catalogued✓ accessible✓ acquired✓ validated✓ mapped✓ · usable for: D05, D11
- **Access:** open (ClinGen terms: curated content released openly, free to capture and use)
- **Targets:** `clingen_gene_curation_list`, `clingen_region_curation_list`, `clingen_haploinsufficiency_gene`, `clingen_triplosensitivity_gene`, `clingen_recurrent_cnv_v2_1`, `clingen_gene_validity`, `clingen_actionability_adult_summary_flat`, `clingen_actionability_adult_assertions`, `clingen_actionability_adult_scoring`, `clingen_actionability_pediatric_summary_flat`, `clingen_actionability_pediatric_assertions`, `clingen_actionability_pediatric_scoring`, `clingen_erepo_variant_classifications`

### `clinvar`: ClinVar
- **Why:** Source of record for clinical classifications; Open Targets only carries a subset (EVA evidence), so the VCF, the summary tables and the VCV release are retrieved directly. The VCV release adds what the VCF drops: individual submissions, structural-variant extents and composite records.
- **Status (all files):** catalogued✓ accessible✓ acquired✓ validated✓ mapped✓ · usable for: D11
- **Access:** open (NCBI public domain policy; submitters may claim rights in some content (NCBI policies page))
- **Targets:** `clinvar_vcf_grch38`, `clinvar_readme`, `clinvar_readme_vcf`, `clinvar_variant_summary`, `clinvar_submission_summary`, `clinvar_vcv_xml_weekly_2026_1018`

### `cpic`: CPIC (Clinical Pharmacogenetics Implementation Consortium)
- **Why:** Small API tables that give allele definitions, allele function with evidence strength and population frequency, phenotype translation and recommendations in one open, versioned source; combined with PharmVar nomenclature and ClinPGx annotations rather than rebuilt.
- **Status (all files):** catalogued✓ accessible✓ acquired✓ validated✓ mapped✓ · usable for: D07
- **Access:** open_with_terms (The CPIC licence URL resolves to the ClinPGx data-use text (CC BY-SA 4.0))
- **Targets:** `cpic_genes`, `cpic_guidelines`, `cpic_api_allele_definition`, `cpic_api_allele_location_value`, `cpic_api_sequence_location`, `cpic_api_allele`, `cpic_api_gene_result`, `cpic_api_gene_result_lookup`, `cpic_api_gene_result_diplotype`, `cpic_api_recommendation`
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

### `ensembl_genomes`: Ensembl genome registry, assemblies and gene annotation (human)
- **Why:** The GRCh38 toplevel FASTA (808 MB) and the release-116 GTF (141 MB) are the reference every coordinate file in the library is checked against, publish checksums, and have no curated alternative inside the library. They are fetched and verified against the Ensembl CHECKSUMS. For human the Ensembl README states that the GTF equals the GENCODE gene set, so GENCODE's own copy is not fetched. The genome registry (existing asset) is kept as is. The remote range reads of the same FASTA stay as an independent route and are tested against the local copy.
- **Status (all files):** catalogued✓ accessible✓ acquired✓ validated✓ mapped✓ · usable for: D14 · existing assets ensembl-human-genomes: downloaded_verified
- **Access:** open (Ensembl imposes no restrictions on access to or use of the data (disclaimer page))
- **Targets:** `ensembl_grch38_dna_toplevel`, `ensembl_grch38_dna_toplevel_fai`, `ensembl_grch38_dna_toplevel_gzi`, `ensembl_grch38_gtf_116`, `ensembl_grch38_reference_regional`

### `gencode_mane`: GENCODE and MANE transcript sets
- **Why:** MANE files are small and fix the transcript choice needed to interpret HGVS across ClinVar, MaveDB and gnomAD. GENCODE's own annotation and FASTA are not fetched: for human the Ensembl README states that the Ensembl GTF equals the GENCODE gene set, and the release-116 GTF and FASTA are held through the ensembl_genomes source (its genebuild was last updated 2025-11 and may differ from GENCODE v50, dated 2026-06-09).
- **Status (all files):** catalogued✓ accessible✓ acquired✓ validated✓ mapped✓ · usable for: D14
- **Access:** open (GENCODE: 'all the GENCODE project data is open access' (GENCODE data-access page), under the EMBL-EBI Terms of Use (free)
- **Targets:** `mane_summary_v1_5`, `mane_gtf_v1_5`

### `giab_stratifications`: GIAB genome stratifications v3.6 (GRCh38): difficult and complex regions
- **Why:** Before a call in a duplicated, polymorphic or low-mappability region is interpreted, the region has to be recognised; the maintained GA4GH/GIAB stratifications are the community standard for this and are small (23 files, a few MB each), so building own masks would repeat established work.
- **Status (all files):** catalogued✓ accessible✓ acquired✓ validated✓ mapped✓ · usable for: D17
- **Access:** open ("Publicly released data are freely available for reuse without embargo" (BACKGROUND file of the GRCh38 set, verified liv)
- **Targets:** `giab_v3_6_grch38_readme`, `giab_v3_6_grch38_background`, `giab_v3_6_grch38_stratification_list`, `giab_v3_6_grch38_md5s`, `giab_v3_6_grch38_alldifficultregions`, `giab_v3_6_grch38_alllowmapandsegdupregions`, `giab_v3_6_grch38_lowmappabilityall`, `giab_v3_6_grch38_segdups`, `giab_v3_6_grch38_segdups_gt10kb`, `giab_v3_6_grch38_mhc`, `giab_v3_6_grch38_kir`, `giab_v3_6_grch38_vdj`, `giab_v3_6_grch38_false_duplications_correct_copy`, `giab_v3_6_grch38_false_duplications_incorrect_copy`, `giab_v3_6_grch38_collapsed_duplication_fp_regions`, `giab_v3_6_grch38_population_cnv_fp_regions`, `giab_v3_6_grch38_l1h_gt500`, `giab_v3_6_grch38_gaps_slop15kb`, `giab_v3_6_grch38_allotherdifficultregions`, `giab_v3_6_grch38_cmrgv1_00_falselyduplicatedgenes`, `giab_v3_6_grch38_badpromoters`, `giab_v3_6_grch38_alltandemrepeats`, `giab_v3_6_grch38_satellites_slop5`, `giab_v3_6_grch38_chrx_par`, `giab_v3_6_grch38_chry_par`, `giab_v3_6_grch38_chrx_xtr`, `giab_v3_6_grch38_chrx_ampliconic`
- **Unverified:** how each individual stratification was derived (only the category summaries of the README were read)

### `gnomad_constraint`: gnomAD constraint metrics (LoF/missense) and pext
- **Why:** Primary source of record and moderate size; the Ensembl LOEUF file is a derivative of the same resource so is not downloaded separately.
- **Status (all files):** catalogued✓ accessible✓ acquired✓ validated✓ mapped✓ · usable for: D05
- **Access:** open (CC0 for primary data (policies page))
- **Targets:** `gnomad_constraint_v4_1_2`, `gnomad_v4_1_lof_curation_incomplete_penetrance`

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
- **Unverified:** coverage_limitations (donor ancestry composition and sampling conditions are from the GTEx publications, not retrieved)

### `hgsvc3`: HGSVC3 release 1.0: assembly-based structural variants from 65 phased genomes (GRCh38)
- **Why:** Long-read assembly-based SV calls resolve insertions, inversions and repeat-embedded events that short-read catalogues miss, with per-haplotype genotypes and built-in flags for segmental duplications and tandem repeats. The symbolic insertion/deletion and inversion call sets total 6.4 MB; the multi-hundred-megabyte SNV and indel releases are deferred.
- **Status (all files):** catalogued✓ accessible✓ acquired✓ validated✓ mapped✓ · usable for: D08, D17
- **Access:** open_with_terms (Data reuse statement of the project (README_HGSVC3_datareuse_statement.md, verified live): in the spirit of the Fort Lau)
- **Targets:** `hgsvc3_grch38_sv_insdel_sym_vcf`, `hgsvc3_grch38_sv_inv_sym_vcf`, `hgsvc3_variant_calls_readme`, `hgsvc3_variant_calls_manifest`
- **Unverified:** the calling and merging methods beyond the release README (PAV and SV-Pop are named; their parameters were not read)

### `hpa`: Human Protein Atlas v25.1: tissue and cell-type RNA expression and protein subcellular location
- **Why:** Tissue and cell-type specificity of gene expression decides whether a regulatory or loss-of-function mechanism is plausible in a tissue of interest, and subcellular location bears on protein-level mechanisms. The three tables are small (22 MB) and fill the cell-type-resolved and protein-localisation gaps that GTEx bulk tissue and Open Targets baseline expression leave.
- **Status (all files):** catalogued✓ accessible✓ acquired✓ validated✓ mapped✓ · usable for: D03, D02
- **Access:** open_with_terms (Creative Commons Attribution 4.0 International for all copyrightable parts of the Human Protein Atlas; data subject to t)
- **Targets:** `hpa_v25_1_rna_tissue_consensus`, `hpa_v25_1_rna_single_cell_type`, `hpa_v25_1_subcellular_location`
- **Unverified:** how the consensus expression levels are computed (the page lists the HPA, GTEx and consensus RNA datasets separately; the combination method was not read), licences of the external sources that the licence page lists

### `interpro`: InterPro
- **Why:** Per-protein retrieval of InterPro and member-database matches (with residue locations) for the UniProt proteins our MaveDB score sets are about; the bulk match files are not needed for that set.
- **Status (all files):** catalogued✓ accessible✓ acquired✓ validated✓ mapped✓ · usable for: D02
- **Access:** open (CC0 1.0 Universal for InterPro, Pfam, PRINTS and SFLD downloadable data; InterProScan under the Apache License (InterPro)
- **Targets:** `interpro_mavedb_proteins`

### `pdbe_sifts`: PDBe SIFTS: PDB chain to UniProt residue mappings
- **Why:** Experimental-structure coverage per UniProt residue is the complement to the AlphaFold DB predicted models and was an explicit gap in the library (PDB not catalogued). The two summary files are small (6.3 MB and 11.0 MB). Structure coordinate files are not mirrored.
- **Status (all files):** catalogued✓ accessible✓ acquired✓ validated✓ mapped✓ · usable for: D02
- **Access:** open (EMBL-EBI Terms of Use: freely available, no additional restriction on contributed data beyond the data owner's terms (te)
- **Targets:** `sifts_pdb_chain_uniprot`, `sifts_uniprot_segments_observed`
- **Unverified:** SIFTS-specific licence statement

### `pgs_catalog`: Polygenic Score (PGS) Catalog
- **Why:** The bulk metadata (scores, publications, traits, samples, cohorts, performance) is small and is the registry of models; the harmonised GRCh38 scoring files are retrieved for a documented selection and applied to a sample by a separate scoring step. Scores are not merged into the variant evidence index: a polygenic score is a model, not an assertion about a variant.
- **Status (all files):** catalogued✓ accessible✓ acquired✓ validated✓ mapped· · usable for: nothing yet
- **Access:** open_with_terms (EMBL-EBI terms of use by default; a scoring file's own licence field governs that score and can be more restrictive (for)
- **Targets:** `pgs_scores_list`, `pgs_metadata_scores`, `pgs_metadata_publications`, `pgs_metadata_efo_traits`, `pgs_metadata_performance_metrics`, `pgs_metadata_evaluation_sample_sets`, `pgs_metadata_score_development_samples`, `pgs_metadata_cohorts`, `pgs_scoring_files_grch38`

### `pharmgkb`: ClinPGx / PharmGKB
- **Why:** Small archives; Open Targets only carries a pharmacogenomics subset.
- **Status (all files):** catalogued✓ accessible✓ acquired✓ validated✓ mapped✓ · usable for: D07
- **Access:** open_with_terms (CC BY-SA 4.0 (data usage policy))
- **Targets:** `clinpgx_variantAnnotations`, `clinpgx_relationships`, `clinpgx_summaryAnnotations`, `clinpgx_guidelineAnnotations`, `clinpgx_drugLabels`, `clinpgx_clinicalAnnotations`, `clinpgx_clinicalVariants`, `clinpgx_variants`, `clinpgx_genes`

### `pharmvar`: PharmVar (Pharmacogene Variation Consortium)
- **Why:** The licence was readable after all (the earlier guessed URLs returned the application shell; the terms are linked from the footer of the rendered page), it permits research use with attribution and share-alike, and the complete database is one 193 MB archive that the site's own button serves. It closes the star-allele sequence gap that CPIC and ClinPGx tables leave.
- **Status (all files):** catalogued✓ accessible✓ acquired✓ validated✓ mapped✓ · usable for: D07
- **Access:** open_with_terms (CC BY-SA 4.0 for the PharmVar database content, with the user agreeing to use the data for research purposes and not to )
- **Targets:** `pharmvar_6_2_29`

### `reactome`: Reactome
- **Why:** Small mapping files; Open Targets also carries a reactome dataset.
- **Status (all files):** catalogued✓ accessible✓ acquired✓ validated✓ mapped✓ · usable for: D10
- **Access:** open (Data licensed CC0; illustrations CC BY 4.0 (license page))
- **Targets:** `reactome_ensembl2reactome`, `reactome_uniprot2reactome_all_levels`, `reactome_pathways`

### `refseq`: NCBI RefSeq: GRCh38.p14 accession map, alternate-locus placement and gene annotation (release RS_2025_08)
- **Why:** ClinVar reports HGVS on RefSeq transcripts and some records on RefSeq NT_/NW_ contigs; the assembly report is the authoritative map from those names to the Ensembl and UCSC names used elsewhere, and the placement files identify the alternate-locus and patch regions where one reference interval has competing representations. The 80 KB report, 45 KB regions file and 78 MB GFF are small relative to their use; sequence files are not needed because the Ensembl FASTA is mirrored.
- **Status (all files):** catalogued✓ accessible✓ acquired✓ validated✓ mapped✓ · usable for: D14, D17
- **Access:** open (NCBI data usage policy: NCBI places no restrictions on the use or distribution of the data in its databases; original su)
- **Targets:** `refseq_grch38p14_assembly_report`, `refseq_grch38p14_assembly_regions`, `refseq_grch38p14_alt_scaffold_placement`, `refseq_grch38p14_feature_count`, `refseq_annotation_release_readme`, `refseq_grch38p14_genomic_gff`
- **Unverified:** RefSeq-specific licence statement

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

## Deferred (9)

### `alliance_genome`: Alliance of Genome Resources (model-organism and human gene knowledge)
- **Why:** Mouse phenotypes (IMPC) and orthology (Open Targets) are already held; the Alliance adds zebrafish, rat, yeast, fly and worm models. Fetch specific phenotype files once the licence is read.
- **Status (all files):** catalogued✓ accessible✓ acquired· validated· mapped· · usable for: nothing yet
- **Access:** unverified (Not retrieved: the Alliance website pages are JavaScript-only and one probe was refused by the content delivery network;)
- **Unverified:** licence, file sizes of the human-relevant files

### `cadd`: CADD: combined annotation dependent depletion
- **Why:** CADD scores for observed variants are carried in gnomAD v4.1 INFO fields, and the full table is tens of gigabytes with non-commercial terms. Fetch only if a variant outside gnomAD needs a score, and read the terms for the intended use first.
- **Status (all files):** catalogued✓ accessible✓ acquired· validated· mapped· · usable for: nothing yet
- **Access:** open_with_terms ("CADD scores are freely available for all non-commercial applications"; commercial use needs a licence from the universi)
- **Unverified:** file sizes, training labels and annotation set of v1.7

### `dbvar`: dbVar (NCBI database of genomic structural variation)
- **Why:** Catalogued only: ClinGen, gnomAD-SV and DGV cover the near-term needs; inspect dbVar contents and currency before any download.
- **Status (all files):** catalogued✓ accessible✓ acquired· validated· mapped· · usable for: nothing yet
- **Access:** open (NCBI data usage policy: NCBI places no restrictions on the use or distribution of the data in its databases; original su)
- **Unverified:** currency, contents

### `gnomad_ld`: gnomAD v2.1.1 linkage disequilibrium matrices
- **Why:** LD for fine-mapping is wanted per region; the matrices are large and GRCh37-only, so fetching waits for a concrete locus and a liftover or a recomputation from phased GRCh38 genotypes.
- **Status (all files):** catalogued✓ accessible✓ acquired· validated· mapped· · usable for: nothing yet
- **Access:** open (gnomAD data are released under CC0 1.0 (as listed for gnomAD in the Open Targets data-source table); the gnomAD LD files)
- **Unverified:** sizes of the matrices, licence statement of the LD files, file format and sizes (the extensions suggest Hail matrices and tables; the gnomAD documentation was not read)

### `intact`: IntAct molecular interaction database
- **Why:** Already integrated in Open Targets (594,648 distinct human gene pairs), so a separate 15 GB archive adds mostly the same pairs. Revisit when interaction-level mutation annotations are needed.
- **Status (all files):** catalogued✓ accessible✓ acquired· validated· mapped· · usable for: nothing yet
- **Access:** open_with_terms (EMBL-EBI terms of use, as listed for IntAct in the Open Targets data-source licence table; the IntAct-specific licence s)
- **Unverified:** licence statement of IntAct itself, size of the human-only subset, quantitative (the definition of the confidence score was not read)

### `mitomap`: MITOMAP: human mitochondrial genome database
- **Why:** The mitochondrial class is a coverage gap. The resource is catalogued so that the gap is explicit; reading its terms and data access conditions by a person is the next step.
- **Status (all files):** catalogued✓ accessible· acquired· validated· mapped· · usable for: nothing yet
- **Access:** blocked_for_automation (Not retrieved.)
- **Unverified:** everything except the existence of the resource and its refusal of automated access, licence, version

### `mpravardb`: MPRAVarDB: massively parallel reporter assay variant database
- **Why:** A compilation of published MPRA studies, several of which are expected in MaveDB (which names MPRA among its assay types). Fetch the data only after the licence and genome build are confirmed and overlap with MaveDB MPRA score sets is matched by publication.
- **Status (all files):** catalogued✓ accessible✓ acquired· validated· mapped· · usable for: nothing yet
- **Access:** unverified (Not stated on the pages read.)
- **Unverified:** licence, version, genome build, scale and direction of the effect measure

### `string_db`: STRING protein-protein association networks
- **Why:** Already integrated in Open Targets (13.2 million rows, 6.5 million distinct gene pairs), and the full network is 124.6 GB. A human-only file could be fetched later if channel subscores are needed.
- **Status (all files):** catalogued✓ accessible✓ acquired· validated· mapped· · usable for: nothing yet
- **Access:** open_with_terms (CC BY 4.0 as listed for STRING in the Open Targets data-source licence table; the STRING licence page itself was not rea)
- **Unverified:** licence statement of STRING itself, size of the human-only files

### `trexplorer`: TRExplorer: genome-wide tandem-repeat catalogue
- **Why:** Wanted for loci outside the curated pathogenic set, but the licence is not stated and the useful annotated tables are about 1 GB. Ask for the terms, then fetch the 45 MB base BED and the extended-loci table first.
- **Status (all files):** catalogued✓ accessible✓ acquired· validated· mapped· · usable for: nothing yet
- **Access:** unverified (Not found: neither the about page, the FAQ nor the bucket listing states a data licence.)
- **Unverified:** licence, locus count of each version

## Not acquired by design (7)

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

### `ipd_imgt_hla`: IPD-IMGT/HLA Database (HLA allele sequences and nomenclature)
- **Why:** The licence forbids derivatives without permission and the maintainers strongly oppose mirroring. The HLA allele sequences are therefore used by reference (queried at their source when an HLA analysis is built), while the GIAB MHC region and the RefSeq HLA gene models are held for locating the region and recognising that calls there need allele-aware handling.
- **Status (all files):** catalogued✓ accessible✓ acquired· validated· mapped· · usable for: nothing yet
- **Access:** open_with_terms (Creative Commons Attribution-NoDerivs: free to copy, distribute, display and make commercial use of the databases with c)
- **Unverified:** number of alleles in the release (only the list head was read)

### `ipd_kir`: IPD-KIR Database (killer-cell immunoglobulin-like receptor allele sequences)
- **Why:** Same terms and request as IPD-IMGT/HLA: no derivatives without permission and mirroring is discouraged. KIR allele sequences are used by reference; the GIAB KIR region and the RefSeq KIR gene models locate the region here.
- **Status (all files):** catalogued✓ accessible✓ acquired· validated· mapped· · usable for: nothing yet
- **Access:** open_with_terms (Creative Commons Attribution-NoDerivs: free to copy, distribute, display and make commercial use with credit by citing t)
- **Unverified:** number of alleles in the release (only the list head was read)

### `ldlink`: LDlink (NCI): linkage disequilibrium tools for 1000 Genomes populations
- **Why:** An on-demand computation service, not a dataset. LD for the populations of interest can be computed from the 1000 Genomes genotypes, and gnomAD LD matrices (catalogued as gnomad_ld) are the bulk alternative; fine-mapped credible sets already come with the association sources.
- **Status (all files):** catalogued✓ accessible✓ acquired· validated· mapped· · usable for: nothing yet
- **Access:** registration_required (Not retrieved: usage conditions of the API (token, rate limits) were not read.)
- **Unverified:** licence (usage conditions of the API), release date

### `pharmcat`: PharmCAT (Pharmacogenomics Clinical Annotation Tool)
- **Why:** A tool that applies the allele definitions and guidelines already held through CPIC, PharmVar and ClinPGx; it carries no evidence that is not in those sources. It is the reference implementation to run, or to compare a future allele caller with, when a pharmacogenomic analysis is built.
- **Status (all files):** catalogued✓ accessible✓ acquired· validated· mapped· · usable for: nothing yet
- **Access:** unverified (Not retrieved: the documentation pages were read but the software licence file was not.)
- **Unverified:** licence, release date
