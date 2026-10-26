# Repository audit

Generated 2026-10-26T02:57:00Z at commit `27cb8cdb` by `inventory/repository_audit.py` (regenerate rather than edit). It records what exists, which lifecycle stage each dataset has reached, whether the checks pass, and whether the documentation says what the implementation does.

## Checks

| Check | Result | Detail |
|---|---|---|
| existing assets re-audited (NHGRI, Ensembl registry, MaveDB score tables) | PASS | {"nhgri-gwas": "downloaded_verified", "nhgri-1000-genomes": "downloaded_verified", "ensembl-human-genomes": "downloaded_verified", "mavedb-human": "downloaded_verified"} |
| preservation guard: tracked files under scrape/ unchanged apart from documented exceptions | PASS | baseline files: 97 | missing: 0 | changed: 0 | documented exceptions: 1 | added: 0 |
| unit tests (acquisition.test_helpers) | PASS | 17 tests |
| catalogue records validate against the vocabulary, evidence files and acquisition targets | PASS | 53 source records valid; 146 acquisition targets; 251 evidence probes |
| repository policy: size, personal-genome file types, credential patterns | PASS | visibility: public (ziadwalieddine-create/dc) |
| all scripts compile | PASS | 30 scripts |
| the claim checker detects a changed number in every claim (mutation self-test) | PASS | 31 claims |

## Documentation claims checked against the implementation: 31 pass, 0 fail, 0 skipped

| Document | Claim (pattern) | Result | Detail |
|---|---|---|---|
| `README.md` | `(\d+) genome/assembly records` | PASS | [563.0] |
| `README.md` | `MaveDB: ([\d,]+) published score sets with a human target \(([\d,]+) b` | PASS | [1402.0, 1320.0, 82.0] |
| `README.md` | `the ([\d,]+) score and count tables` | PASS | [1980.0] |
| `README.md` | `an additive set of (\d+) superseded MaveDB score sets` | PASS | [52.0] |
| `README.md` | `(\d+) tutorial slide decks` | PASS | [11.0] |
| `scrape/mavedb-human/README.md` | `([\d,]+) in total` | PASS | [1402.0] |
| `scrape/mavedb-human/README.md` | `\| `organism_filter` \|.*?\| ([\d,]+) \|` | PASS | [1320.0] |
| `scrape/mavedb-human/README.md` | `\| `accession_target` \|.*?\| (\d+) \|` | PASS | [82.0] |
| `scrape/mavedb-human/README.md` | `Together: ([\d,]+) experiments, ([\d,]+) distinct target names, ([\d,]` | PASS | [1110.0, 712.0, 5170093.0] |
| `scrape/mavedb-human/README.md` | `([\d,]+) files, about [\d]+ MB gzipped, ([\d,]+) score rows and ([\d,]` | PASS | [1980.0, 5170093.0, 1544328.0] |
| `scrape/mavedb-human/README.md` | `for the (\d+) sets with count columns` | PASS | [578.0] |
| `scrape/mavedb-human/README.md` | `(\d+) sets have `processingState: incomplete`` | PASS | [322.0] |
| `scrape/mavedb-human/README.md` | ``complete` for only (\d+) of ([\d,]+) \((\d+) incomplete, (\d+) failed` | PASS | [888.0, 1402.0, 321.0, 181.0, 9.0, 1.0, 2.0] |
| `scrape/mavedb-human/README.md` | `([\d,]+) CC0, (\d+) CC BY 4\.0, (\d+) CC BY-SA 4\.0, (\d+) CC BY-NC-SA` | PASS | [1340.0, 45.0, 11.0, 6.0] |
| `scrape/mavedb-human/README.md` | `only `urn:mavedb:00000053-a-1`, ([\d,]+) variants` | PASS | [648022.0] |
| `scrape/ensembl-human-genomes/README.md` | `^(\d+) records, all in integrated release 2026-07` | PASS | [563.0] |
| `scrape/ensembl-human-genomes/README.md` | `\((\d+) have `hprc` in the assembly name` | PASS | [440.0] |
| `scrape/ensembl-human-genomes/README.md` | `(\d+) records have variation data and (\d+) have regulation` | PASS | [96.0, 2.0] |
| `acquisition/README.md` | `\((\d[\d,]*) files, no differences\)` | PASS | [1297.0] |
| `acquisition/README.md` | `\((\d[\d,]*) of ([\d,]+), none elsewhere\)` | PASS | [3838.0, 71697556.0] |
| `acquisition/README.md` | `# (\d+) tests, no data needed` | PASS | [17.0] |
| `catalogue/MASTER_INVENTORY.md` | `(\d+) sources catalogued` | PASS | [53.0] |
| `latest:AUDIT_LOG.md` | ``validate\.py`: (\d+) of (\d+) targets pass` | PASS | [146.0, 146.0] |
| `latest:AUDIT_LOG.md` | `(\d+) manifest records` | PASS | [146.0] |
| `latest:AUDIT_LOG.md` | ``compat_tests\.py`: (\d+) tests, (\d+) PASS, (\d+) FAIL, (\d+) INFO` | PASS | [115.0, 101.0, 0.0, 14.0] |
| `latest:AUDIT_LOG.md` | `test_helpers`: (\d+) tests pass` | PASS | [17.0] |
| `catalogue/gaps.json` | `GTF: ([\d,]+) genes, ([\d,]+) transcripts, ([\d,]+) MANE Select` | PASS | [78941.0, 646577.0, 19299.0] |
| `catalogue/gaps.json` | `toplevel FASTA: (\d+) sequences` | PASS | [706.0] |
| `catalogue/gaps.json` | `\(([\d.]+)M PDB chain segments, ([\d,]+) PDB entries, ([\d,]+) UniProt` | PASS | [1.05, 246198.0, 76796.0] |
| `catalogue/gaps.json` | `\(([\d,]+) VCF files, ([\d,]+) variant records\)` | PASS | [7865.0, 35172.0] |
| `catalogue/gaps.json` | `with ([\d,]+) haplotype names` | PASS | [2458.0] |

## Branches

| Branch | Tip | Date | Tracked files | Ahead of / behind origin/main | Last commit |
|---|---|---|---|---|---|
| `ccr-b13fd75a-ra9gfy` | `27cb8cd` | 2026-10-25 | 459 | +6 / -0 | Add evidence profiles, mechanism matrix, repository audit and evidence-support checks |
| `data/mavedb-scores` | `1b58196` | 2026-10-25 | 2061 | +8 / -22 | Add 52 superseded MaveDB score sets (79 tables, CC0) as a data extension |
| `main` | `f388b5f` | 2026-09-30 | 90 | +0 / -19 | Merge pull request #1 from ziadwalieddine-create/land/dc-extended-pack |
| `origin/ccr-b13fd75a-ra9gfy` | `27cb8cd` | 2026-10-25 | 459 | +6 / -0 | Add evidence profiles, mechanism matrix, repository audit and evidence-support checks |
| `origin/data/mavedb-scores` | `1b58196` | 2026-10-25 | 2061 | +8 / -22 | Add 52 superseded MaveDB score sets (79 tables, CC0) as a data extension |
| `origin/land/dc-extended-pack` | `3547d0a` | 2026-09-30 | 90 | +0 / -20 | Land dc-extended skill pack |
| `origin/main` | `ddafb14` | 2026-10-25 | 98 | +0 / -0 | Merge pull request #6 from ziadwalieddine-create/ccr-b13fd75a-ra9gfy |

## Datasets (acquisition targets)

146 listed · 146 acquired (manifest ok and the file present with the recorded size) · 146 verified (validator passed) · 107 mapped (a passing join or consistency test used it) · 0 integrated (a documented, tested query layer returning records by variant, gene, protein or region; none is built in this phase).

| Target | Source record | Tier | Listed | Acquired | Verified | Mapped | Integrated | Size |
|---|---|---|---|---|---|---|---|---|
| `mavedb_zenodo_dump_v4` | `mavedb` | A | ✓ | ✓ | ✓ | ✓ | · | 1907.7 MB |
| `clinvar_vcf_grch38` | `clinvar` | A | ✓ | ✓ | ✓ | ✓ | · | 200.7 MB |
| `clinvar_readme` | `clinvar` | A | ✓ | ✓ | ✓ | · | · | 0.0 MB |
| `clinvar_readme_vcf` | `clinvar` | A | ✓ | ✓ | ✓ | · | · | 0.0 MB |
| `clingen_gene_curation_list` | `clingen` | A | ✓ | ✓ | ✓ | ✓ | · | 0.2 MB |
| `clingen_region_curation_list` | `clingen` | A | ✓ | ✓ | ✓ | · | · | 0.1 MB |
| `clingen_haploinsufficiency_gene` | `clingen` | A | ✓ | ✓ | ✓ | ✓ | · | 0.0 MB |
| `clingen_triplosensitivity_gene` | `clingen` | A | ✓ | ✓ | ✓ | ✓ | · | 0.0 MB |
| `clingen_recurrent_cnv_v2_1` | `clingen` | A | ✓ | ✓ | ✓ | · | · | 0.0 MB |
| `clingen_gene_validity` | `clingen` | A | ✓ | ✓ | ✓ | ✓ | · | 1.1 MB |
| `gnomad_constraint_v4_1_2` | `gnomad` | A | ✓ | ✓ | ✓ | ✓ | · | 574.5 MB |
| `gnomad_exome_cnv_v4_1` | `gnomad` | A | ✓ | ✓ | ✓ | ✓ | · | 8.9 MB |
| `gwas_catalog_associations` | `gwas_catalog` | A | ✓ | ✓ | ✓ | ✓ | · | 73.0 MB |
| `reactome_ensembl2reactome` | `reactome` | A | ✓ | ✓ | ✓ | ✓ | · | 183.3 MB |
| `reactome_uniprot2reactome_all_levels` | `reactome` | A | ✓ | ✓ | ✓ | ✓ | · | 117.8 MB |
| `reactome_pathways` | `reactome` | A | ✓ | ✓ | ✓ | · | · | 1.6 MB |
| `clinpgx_variantAnnotations` | `pharmgkb` | A | ✓ | ✓ | ✓ | · | · | 4.3 MB |
| `clinpgx_relationships` | `pharmgkb` | A | ✓ | ✓ | ✓ | ✓ | · | 4.3 MB |
| `clinpgx_summaryAnnotations` | `pharmgkb` | A | ✓ | ✓ | ✓ | · | · | 1.2 MB |
| `clinpgx_guidelineAnnotations` | `pharmgkb` | A | ✓ | ✓ | ✓ | · | · | 0.9 MB |
| `clinpgx_drugLabels` | `pharmgkb` | A | ✓ | ✓ | ✓ | · | · | 0.1 MB |
| `cpic_genes` | `cpic` | A | ✓ | ✓ | ✓ | ✓ | · | 0.1 MB |
| `cpic_guidelines` | `cpic` | A | ✓ | ✓ | ✓ | · | · | 0.0 MB |
| `gtex_v11_gene_median_tpm` | `gtex` | A | ✓ | ✓ | ✓ | ✓ | · | 10.1 MB |
| `gtex_datasets_api` | `gtex` | A | ✓ | ✓ | ✓ | · | · | 0.0 MB |
| `mane_summary_v1_5` | `gencode` | A | ✓ | ✓ | ✓ | ✓ | · | 1.1 MB |
| `mane_gtf_v1_5` | `gencode` | A | ✓ | ✓ | ✓ | ✓ | · | 8.6 MB |
| `igsr_1000g_samples` | `1000genomes` | A | ✓ | ✓ | ✓ | ✓ | · | 0.1 MB |
| `uniprot_human_reviewed` | `uniprot` | A | ✓ | ✓ | ✓ | ✓ | · | 18.3 MB |
| `biogrid_orcs_human` | `biogrid_orcs` | A | ✓ | ✓ | ✓ | ✓ | · | 752.7 MB |
| `alphamissense_gene_hg38` | `alphamissense` | A | ✓ | ✓ | ✓ | ✓ | · | 0.3 MB |
| `alphamissense_hg38` | `alphamissense` | B | ✓ | ✓ | ✓ | ✓ | · | 643.0 MB |
| `ensembl_vep_mavedb_variants` | `ensembl_vep` | A | ✓ | ✓ | ✓ | ✓ | · | 560.6 MB |
| `ensembl_vep_mavedb_variants_tbi` | `ensembl_vep` | A | ✓ | ✓ | ✓ | · | · | 0.0 MB |
| `strchive_hg38_general` | `strchive` | A | ✓ | ✓ | ✓ | ✓ | · | 0.0 MB |
| `strchive_hg38_trgt` | `strchive` | A | ✓ | ✓ | ✓ | · | · | 0.0 MB |
| `strchive_hg38_stranger_json` | `strchive` | A | ✓ | ✓ | ✓ | · | · | 0.0 MB |
| `dgv_grch38_variants_2025_12_01` | `dgv` | B | ✓ | ✓ | ✓ | ✓ | · | 364.2 MB |
| `screen_ccre_pls` | `encode` | A | ✓ | ✓ | ✓ | ✓ | · | 2.6 MB |
| `screen_ccre_pels` | `encode` | A | ✓ | ✓ | ✓ | ✓ | · | 13.7 MB |
| `screen_ccre_dels` | `encode` | A | ✓ | ✓ | ✓ | ✓ | · | 80.7 MB |
| `screen_ccre_ca_ctcf` | `encode` | A | ✓ | ✓ | ✓ | ✓ | · | 7.3 MB |
| `screen_ccre_ca_h3k4me3` | `encode` | A | ✓ | ✓ | ✓ | ✓ | · | 4.8 MB |
| `screen_ccre_ca_tf` | `encode` | A | ✓ | ✓ | ✓ | ✓ | · | 1.5 MB |
| `screen_ccre_ca` | `encode` | A | ✓ | ✓ | ✓ | ✓ | · | 13.0 MB |
| `screen_ccre_tf` | `encode` | A | ✓ | ✓ | ✓ | ✓ | · | 5.6 MB |
| `opentargets_26_09_target` | `opentargets` | A | ✓ | ✓ | ✓ | ✓ | · | 21.6 MB |
| `opentargets_26_09_go` | `opentargets` | A | ✓ | ✓ | ✓ | · | · | 2.9 MB |
| `opentargets_26_09_so` | `opentargets` | A | ✓ | ✓ | ✓ | · | · | 0.0 MB |
| `opentargets_26_09_homology` | `opentargets` | A | ✓ | ✓ | ✓ | ✓ | · | 24.5 MB |
| `opentargets_26_09_biosample` | `opentargets` | A | ✓ | ✓ | ✓ | · | · | 6.0 MB |
| `opentargets_26_09_disease` | `opentargets` | A | ✓ | ✓ | ✓ | ✓ | · | 6.3 MB |
| `opentargets_26_09_disease_hpo` | `opentargets` | A | ✓ | ✓ | ✓ | · | · | 1.7 MB |
| `opentargets_26_09_disease_phenotype` | `opentargets` | A | ✓ | ✓ | ✓ | · | · | 6.5 MB |
| `opentargets_26_09_mouse_phenotype` | `opentargets` | A | ✓ | ✓ | ✓ | ✓ | · | 6.8 MB |
| `opentargets_26_09_evidence_genomics_england` | `opentargets` | A | ✓ | ✓ | ✓ | ✓ | · | 2.7 MB |
| `opentargets_26_09_evidence_gene2phenotype` | `opentargets` | A | ✓ | ✓ | ✓ | ✓ | · | 0.4 MB |
| `opentargets_26_09_evidence_orphanet` | `opentargets` | A | ✓ | ✓ | ✓ | ✓ | · | 0.5 MB |
| `opentargets_26_09_evidence_uniprot_variants` | `opentargets` | A | ✓ | ✓ | ✓ | ✓ | · | 1.9 MB |
| `opentargets_26_09_pharmacogenomics` | `opentargets` | A | ✓ | ✓ | ✓ | ✓ | · | 3.4 MB |
| `opentargets_26_09_evidence_gene_burden` | `opentargets` | A | ✓ | ✓ | ✓ | ✓ | · | 3.9 MB |
| `opentargets_26_09_target_essentiality` | `opentargets` | A | ✓ | ✓ | ✓ | ✓ | · | 246.6 MB |
| `opentargets_26_09_evidence_crispr` | `opentargets` | A | ✓ | ✓ | ✓ | ✓ | · | 0.0 MB |
| `opentargets_26_09_evidence_crispr_screen` | `opentargets` | A | ✓ | ✓ | ✓ | ✓ | · | 0.9 MB |
| `opentargets_26_09_evidence_impc` | `opentargets` | A | ✓ | ✓ | ✓ | ✓ | · | 906.9 MB |
| `opentargets_26_09_study` | `opentargets` | A | ✓ | ✓ | ✓ | ✓ | · | 97.7 MB |
| `opentargets_26_09_credible_set` | `opentargets` | A | ✓ | ✓ | ✓ | ✓ | · | 2990.6 MB |
| `opentargets_26_09_l2g_prediction` | `opentargets` | A | ✓ | ✓ | ✓ | ✓ | · | 509.6 MB |
| `opentargets_26_09_variant` | `opentargets` | A | ✓ | ✓ | ✓ | ✓ | · | 1967.3 MB |
| `opentargets_26_09_enhancer_to_gene` | `opentargets` | A | ✓ | ✓ | ✓ | ✓ | · | 1770.1 MB |
| `opentargets_26_09_interaction` | `opentargets` | A | ✓ | ✓ | ✓ | ✓ | · | 251.9 MB |
| `opentargets_26_09_interaction_evidence` | `opentargets` | A | ✓ | ✓ | ✓ | ✓ | · | 344.1 MB |
| `opentargets_26_09_evidence_uniprot_literature` | `opentargets` | A | ✓ | ✓ | ✓ | ✓ | · | 0.5 MB |
| `gnomad_sv_v4_1_sites_vcf` | `gnomad` | B | ✓ | ✓ | ✓ | ✓ | · | 1743.5 MB |
| `alphafold_db_mavedb_proteins` | `alphafold_db` | A | ✓ | ✓ | ✓ | ✓ | · | 356.7 MB |
| `interpro_mavedb_proteins` | `interpro` | A | ✓ | ✓ | ✓ | ✓ | · | 8.0 MB |
| `gnomad_v4_1_exomes_regional` | `gnomad` | A | ✓ | ✓ | ✓ | ✓ | · | 7.9 MB |
| `gnomad_v4_1_genomes_regional` | `gnomad` | A | ✓ | ✓ | ✓ | ✓ | · | 2.4 MB |
| `uniprot_humsavar` | `uniprot` | A | ✓ | ✓ | ✓ | ✓ | · | 8.7 MB |
| `gnomad_v4_1_joint_regional` | `gnomad` | A | ✓ | ✓ | ✓ | ✓ | · | 4.3 MB |
| `ensembl_grch38_reference_regional` | `ensembl_genomes` | A | ✓ | ✓ | ✓ | ✓ | · | 0.1 MB |
| `ensembl_grch38_dna_toplevel` | `ensembl_genomes` | B | ✓ | ✓ | ✓ | ✓ | · | 847.7 MB |
| `ensembl_grch38_dna_toplevel_fai` | `ensembl_genomes` | B | ✓ | ✓ | ✓ | ✓ | · | 0.0 MB |
| `ensembl_grch38_dna_toplevel_gzi` | `ensembl_genomes` | B | ✓ | ✓ | ✓ | · | · | 0.8 MB |
| `ensembl_grch38_gtf_116` | `ensembl_genomes` | B | ✓ | ✓ | ✓ | ✓ | · | 141.1 MB |
| `sifts_pdb_chain_uniprot` | `pdbe_sifts` | A | ✓ | ✓ | ✓ | ✓ | · | 6.3 MB |
| `sifts_uniprot_segments_observed` | `pdbe_sifts` | A | ✓ | ✓ | ✓ | ✓ | · | 11.0 MB |
| `gnomad_v4_1_lof_curation_incomplete_penetrance` | `gnomad_constraint` | A | ✓ | ✓ | ✓ | ✓ | · | 0.1 MB |
| `pharmvar_6_2_29` | `pharmvar` | B | ✓ | ✓ | ✓ | ✓ | · | 193.1 MB |
| `refseq_grch38p14_assembly_report` | `refseq` | A | ✓ | ✓ | ✓ | ✓ | · | 0.1 MB |
| `refseq_grch38p14_assembly_regions` | `refseq` | A | ✓ | ✓ | ✓ | · | · | 0.0 MB |
| `refseq_grch38p14_alt_scaffold_placement` | `refseq` | A | ✓ | ✓ | ✓ | · | · | 0.1 MB |
| `refseq_grch38p14_feature_count` | `refseq` | A | ✓ | ✓ | ✓ | · | · | 0.0 MB |
| `refseq_annotation_release_readme` | `refseq` | A | ✓ | ✓ | ✓ | · | · | 0.0 MB |
| `refseq_grch38p14_genomic_gff` | `refseq` | B | ✓ | ✓ | ✓ | ✓ | · | 78.2 MB |
| `giab_v3_6_grch38_readme` | `giab_stratifications` | A | ✓ | ✓ | ✓ | · | · | 0.0 MB |
| `giab_v3_6_grch38_background` | `giab_stratifications` | A | ✓ | ✓ | ✓ | · | · | 0.0 MB |
| `giab_v3_6_grch38_stratification_list` | `giab_stratifications` | A | ✓ | ✓ | ✓ | · | · | 0.0 MB |
| `giab_v3_6_grch38_md5s` | `giab_stratifications` | A | ✓ | ✓ | ✓ | · | · | 0.0 MB |
| `giab_v3_6_grch38_alldifficultregions` | `giab_stratifications` | A | ✓ | ✓ | ✓ | ✓ | · | 31.5 MB |
| `giab_v3_6_grch38_alllowmapandsegdupregions` | `giab_stratifications` | A | ✓ | ✓ | ✓ | ✓ | · | 3.5 MB |
| `giab_v3_6_grch38_lowmappabilityall` | `giab_stratifications` | A | ✓ | ✓ | ✓ | ✓ | · | 4.5 MB |
| `giab_v3_6_grch38_segdups` | `giab_stratifications` | A | ✓ | ✓ | ✓ | ✓ | · | 0.1 MB |
| `giab_v3_6_grch38_segdups_gt10kb` | `giab_stratifications` | A | ✓ | ✓ | ✓ | ✓ | · | 0.0 MB |
| `giab_v3_6_grch38_mhc` | `giab_stratifications` | A | ✓ | ✓ | ✓ | ✓ | · | 0.0 MB |
| `giab_v3_6_grch38_kir` | `giab_stratifications` | A | ✓ | ✓ | ✓ | ✓ | · | 0.0 MB |
| `giab_v3_6_grch38_vdj` | `giab_stratifications` | A | ✓ | ✓ | ✓ | ✓ | · | 0.0 MB |
| `giab_v3_6_grch38_false_duplications_correct_copy` | `giab_stratifications` | A | ✓ | ✓ | ✓ | ✓ | · | 0.0 MB |
| `giab_v3_6_grch38_false_duplications_incorrect_copy` | `giab_stratifications` | A | ✓ | ✓ | ✓ | ✓ | · | 0.0 MB |
| `giab_v3_6_grch38_collapsed_duplication_fp_regions` | `giab_stratifications` | A | ✓ | ✓ | ✓ | ✓ | · | 0.0 MB |
| `giab_v3_6_grch38_population_cnv_fp_regions` | `giab_stratifications` | A | ✓ | ✓ | ✓ | ✓ | · | 0.0 MB |
| `giab_v3_6_grch38_l1h_gt500` | `giab_stratifications` | A | ✓ | ✓ | ✓ | ✓ | · | 0.0 MB |
| `giab_v3_6_grch38_gaps_slop15kb` | `giab_stratifications` | A | ✓ | ✓ | ✓ | ✓ | · | 0.0 MB |
| `giab_v3_6_grch38_allotherdifficultregions` | `giab_stratifications` | A | ✓ | ✓ | ✓ | ✓ | · | 0.3 MB |
| `giab_v3_6_grch38_cmrgv1_00_falselyduplicatedgenes` | `giab_stratifications` | A | ✓ | ✓ | ✓ | ✓ | · | 0.0 MB |
| `giab_v3_6_grch38_badpromoters` | `giab_stratifications` | A | ✓ | ✓ | ✓ | ✓ | · | 0.0 MB |
| `giab_v3_6_grch38_alltandemrepeats` | `giab_stratifications` | A | ✓ | ✓ | ✓ | ✓ | · | 9.2 MB |
| `giab_v3_6_grch38_satellites_slop5` | `giab_stratifications` | A | ✓ | ✓ | ✓ | ✓ | · | 0.0 MB |
| `giab_v3_6_grch38_chrx_par` | `giab_stratifications` | A | ✓ | ✓ | ✓ | ✓ | · | 0.0 MB |
| `giab_v3_6_grch38_chry_par` | `giab_stratifications` | A | ✓ | ✓ | ✓ | ✓ | · | 0.0 MB |
| `giab_v3_6_grch38_chrx_xtr` | `giab_stratifications` | A | ✓ | ✓ | ✓ | ✓ | · | 0.0 MB |
| `giab_v3_6_grch38_chrx_ampliconic` | `giab_stratifications` | A | ✓ | ✓ | ✓ | ✓ | · | 0.0 MB |
| `proteingym_v1_3_dms_substitutions` | `proteingym` | A | ✓ | ✓ | ✓ | ✓ | · | 43.0 MB |
| `proteingym_v1_3_dms_indels` | `proteingym` | A | ✓ | ✓ | ✓ | · | · | 6.7 MB |
| `proteingym_v1_3_dms_substitutions_reference` | `proteingym` | A | ✓ | ✓ | ✓ | ✓ | · | 0.2 MB |
| `proteingym_v1_3_dms_indels_reference` | `proteingym` | A | ✓ | ✓ | ✓ | ✓ | · | 0.0 MB |
| `hpa_v25_1_rna_tissue_consensus` | `hpa` | A | ✓ | ✓ | ✓ | ✓ | · | 5.3 MB |
| `hpa_v25_1_rna_single_cell_type` | `hpa` | A | ✓ | ✓ | ✓ | ✓ | · | 16.3 MB |
| `hpa_v25_1_subcellular_location` | `hpa` | A | ✓ | ✓ | ✓ | ✓ | · | 0.3 MB |
| `biogrid_interactions_organism_tab3` | `biogrid_interactions` | A | ✓ | ✓ | ✓ | ✓ | · | 188.9 MB |
| `hgsvc3_grch38_sv_insdel_sym_vcf` | `hgsvc3` | A | ✓ | ✓ | ✓ | ✓ | · | 6.4 MB |
| `hgsvc3_grch38_sv_inv_sym_vcf` | `hgsvc3` | A | ✓ | ✓ | ✓ | ✓ | · | 0.0 MB |
| `hgsvc3_variant_calls_readme` | `hgsvc3` | A | ✓ | ✓ | ✓ | · | · | 0.0 MB |
| `hgsvc3_variant_calls_manifest` | `hgsvc3` | A | ✓ | ✓ | ✓ | · | · | 0.0 MB |
| `cpic_api_allele_definition` | `cpic` | A | ✓ | ✓ | ✓ | ✓ | · | 0.2 MB |
| `cpic_api_allele_location_value` | `cpic` | A | ✓ | ✓ | ✓ | · | · | 0.3 MB |
| `cpic_api_sequence_location` | `cpic` | A | ✓ | ✓ | ✓ | ✓ | · | 0.3 MB |
| `cpic_api_allele` | `cpic` | A | ✓ | ✓ | ✓ | · | · | 1.1 MB |
| `cpic_api_gene_result` | `cpic` | A | ✓ | ✓ | ✓ | · | · | 0.1 MB |
| `cpic_api_gene_result_lookup` | `cpic` | A | ✓ | ✓ | ✓ | · | · | 0.1 MB |
| `cpic_api_gene_result_diplotype` | `cpic` | A | ✓ | ✓ | ✓ | · | · | 22.2 MB |
| `cpic_api_recommendation` | `cpic` | A | ✓ | ✓ | ✓ | · | · | 2.4 MB |
| `impc_r24_readme` | `impc` | A | ✓ | ✓ | ✓ | · | · | 0.0 MB |
| `impc_r24_genotype_phenotype_assertions_all` | `impc` | A | ✓ | ✓ | ✓ | ✓ | · | 5.0 MB |
| `impc_r24_gene_mp_association_1` | `impc` | A | ✓ | ✓ | ✓ | · | · | 0.2 MB |
| `impc_r24_viability` | `impc` | A | ✓ | ✓ | ✓ | ✓ | · | 0.5 MB |

## Manifests

| File | Records | What it records |
|---|---|---|
| `acquisition/manifest.jsonl` | 146 | 146 ok |
| `acquisition/manifest_files/alphafold_db_mavedb_proteins.tsv` | 1621 | per-file table |
| `acquisition/manifest_files/ensembl_grch38_reference_regional.tsv` | 3 | per-file table |
| `acquisition/manifest_files/gnomad_v4_1_exomes_regional.tsv` | 24 | per-file table |
| `acquisition/manifest_files/gnomad_v4_1_genomes_regional.tsv` | 24 | per-file table |
| `acquisition/manifest_files/gnomad_v4_1_joint_regional.tsv` | 24 | per-file table |
| `acquisition/manifest_files/interpro_mavedb_proteins.tsv` | 549 | per-file table |
| `acquisition/manifest_files/opentargets_26_09_biosample.tsv` | 1 | per-file table |
| `acquisition/manifest_files/opentargets_26_09_credible_set.tsv` | 200 | per-file table |
| `acquisition/manifest_files/opentargets_26_09_disease.tsv` | 1 | per-file table |
| `acquisition/manifest_files/opentargets_26_09_disease_hpo.tsv` | 1 | per-file table |
| `acquisition/manifest_files/opentargets_26_09_disease_phenotype.tsv` | 1 | per-file table |
| `acquisition/manifest_files/opentargets_26_09_enhancer_to_gene.tsv` | 30 | per-file table |
| `acquisition/manifest_files/opentargets_26_09_evidence_crispr.tsv` | 1 | per-file table |
| `acquisition/manifest_files/opentargets_26_09_evidence_crispr_screen.tsv` | 1 | per-file table |
| `acquisition/manifest_files/opentargets_26_09_evidence_gene2phenotype.tsv` | 1 | per-file table |
| `acquisition/manifest_files/opentargets_26_09_evidence_gene_burden.tsv` | 1 | per-file table |
| `acquisition/manifest_files/opentargets_26_09_evidence_genomics_england.tsv` | 1 | per-file table |
| `acquisition/manifest_files/opentargets_26_09_evidence_impc.tsv` | 25 | per-file table |
| `acquisition/manifest_files/opentargets_26_09_evidence_orphanet.tsv` | 1 | per-file table |
| `acquisition/manifest_files/opentargets_26_09_evidence_uniprot_literature.tsv` | 1 | per-file table |
| `acquisition/manifest_files/opentargets_26_09_evidence_uniprot_variants.tsv` | 1 | per-file table |
| `acquisition/manifest_files/opentargets_26_09_go.tsv` | 1 | per-file table |
| `acquisition/manifest_files/opentargets_26_09_homology.tsv` | 2 | per-file table |
| `acquisition/manifest_files/opentargets_26_09_interaction.tsv` | 2 | per-file table |
| `acquisition/manifest_files/opentargets_26_09_interaction_evidence.tsv` | 6 | per-file table |
| `acquisition/manifest_files/opentargets_26_09_l2g_prediction.tsv` | 200 | per-file table |
| `acquisition/manifest_files/opentargets_26_09_mouse_phenotype.tsv` | 1 | per-file table |
| `acquisition/manifest_files/opentargets_26_09_pharmacogenomics.tsv` | 1 | per-file table |
| `acquisition/manifest_files/opentargets_26_09_so.tsv` | 1 | per-file table |
| `acquisition/manifest_files/opentargets_26_09_study.tsv` | 1 | per-file table |
| `acquisition/manifest_files/opentargets_26_09_target.tsv` | 3 | per-file table |
| `acquisition/manifest_files/opentargets_26_09_target_essentiality.tsv` | 4 | per-file table |
| `acquisition/manifest_files/opentargets_26_09_variant.tsv` | 25 | per-file table |
| `scrape/mavedb-human/scores_manifest.jsonl` | 1980 | sha256 and row count of each score and count table |
| `inventory/baseline_checksums.sha256` | 97 | sha256 of every tracked file under scrape/ at the baseline |
| `catalogue/evidence/evidence_index.jsonl` | 295 | append-only record of every probe (URL, time, status, sha256) |

## Documents

| File | Title | Lines | Last commit |
|---|---|---|---|
| `AUDIT_LOG.md` | Audit log | 200 | 2026-10-25 |
| `README.md` | dc: human genome reference data | 25 | 2026-10-25 |
| `acquisition/README.md` | Acquisition, validation and compatibility tests | 45 | 2026-10-25 |
| `catalogue/ACQUIRED_FILES.md` | Acquired files register | 154 | 2026-10-25 |
| `catalogue/ATTRIBUTION.md` | Attribution and licence notes for the data we hold | 42 | 2026-10-25 |
| `catalogue/COVERAGE_MATRIX.md` | Coverage matrix: evidence dimensions by lifecycle stage | 229 | 2026-10-25 |
| `catalogue/EVIDENCE_PROFILES.md` | Evidence profiles | 856 | 2026-10-25 |
| `catalogue/GAP_MATRIX.md` | Gap matrix | 66 | 2026-10-25 |
| `catalogue/MASTER_INVENTORY.md` | Master inventory of reusable genomic resources | 149 | 2026-10-25 |
| `catalogue/MECHANISM_MATRIX.md` | Mechanism matrix: variant types by mechanism axis | 322 | 2026-10-25 |
| `catalogue/README.md` | Evidence-source catalogue | 47 | 2026-10-25 |
| `catalogue/REUSE_PLAN.md` | Source-by-source reuse, acquisition and integration plan | 350 | 2026-10-25 |
| `inventory/EXISTING_ASSETS.md` | Existing assets: inventory and verified status | 49 | 2026-10-25 |
| `inventory/REPOSITORY_AUDIT.md` | Repository audit | 251 | 2026-10-25 |
| `scrape/ensembl-human-genomes/README.md` | Ensembl genome selector: human (taxonomy ID 9606) | 20 | 2026-10-25 |
| `scrape/mavedb-human/README.md` | MaveDB: human score sets | 45 | 2026-10-25 |
| `scrape/nhgri-1000-genomes/README.md` | NHGRI: 1000 Genomes Project (archived page scrape) | 22 | 2026-10-25 |
| `scrape/nhgri-gwas/README.md` | NHGRI Talking Glossary: GWAS (scrape) | 17 | 2026-10-25 |

## Scripts

| File | Lines | Compiles | Purpose |
|---|---|---|---|
| `acquisition/acquire.py` | 325 | True | Download external resources listed in targets.json, verify them and record a manifest. |
| `acquisition/compat_opentargets.py` | 272 | True | Joins between the acquired Open Targets 26.09 Parquet datasets and the MaveDB human set (called from compat_te |
| `acquisition/compat_reference.py` | 405 | True | Tests that tie the mirrored GRCh38 reference sequence and the Ensembl release-116 gene annotation to the other |
| `acquisition/compat_tests.py` | 859 | True | Compatibility tests: do the acquired resources join to the existing dc assets (MaveDB human set)? |
| `acquisition/extract_dump.py` | 57 | True | Extract, from the MaveDB bulk dump v4, the members that the compatibility tests and the superseded-set builder |
| `acquisition/fetch_gnomad_regional.py` | 213 | True | Slice the gnomAD v4.1 exome and genome sites VCFs to the positions of MaveDB-mapped SNVs, by remote tabix. |
| `acquisition/fetch_per_protein.py` | 178 | True | Fetch per-protein resources for the UniProt proteins that the MaveDB human score sets are about. |
| `acquisition/fetch_reference_regional.py` | 155 | True | Read the GRCh38 reference sequence for the regions that contain MaveDB-mapped SNV positions, by remote random  |
| `acquisition/intervals.py` | 20 | True | Exact stabbing queries over half-open intervals, used by every genomic-overlap test in compat_tests.py. |
| `acquisition/mavedb_snvs.py` | 77 | True | Genomic SNV positions of the MaveDB human set, from the two mapping products that exist (see compat_tests.py). |
| `acquisition/mavedb_superseded/build_superseded.py` | 84 | True | Add the MaveDB score sets that the public search hides because they are superseded versions. |
| `acquisition/test_helpers.py` | 165 | True | Unit tests for the small functions that every overlap and join number in this directory rests on. |
| `acquisition/validate.py` | 1758 | True | Validate downloaded external files against each source's own description. |
| `catalogue/build_catalogue.py` | 572 | True | Validate catalogue/sources/*.json and generate the inventory, coverage matrix, gap matrix and reuse plan. |
| `catalogue/derive/evidence_support.py` | 222 | True | Find statements in the source records that none of their stored evidence supports. |
| `catalogue/derive/gnomad_sv_header.py` | 38 | True | Record the INFO and FILTER definitions of the gnomAD v4.1 SV sites VCF, read from the header of the downloaded |
| `catalogue/derive/gnomad_vcf_header.py` | 44 | True | Record the INFO/FILTER definitions of the gnomAD v4.1 exome and genome sites VCFs, read from the VCF headers b |
| `catalogue/derive/openebench_counts.py` | 69 | True | Count registry tools relevant to human variant interpretation, from the OpenEBench monitor snapshot. |
| `catalogue/derive/opentargets_datasets.py` | 40 | True | Extract the dataset inventory (name, category, one-line description) from the rendered Open Targets downloads  |
| `catalogue/derive/opentargets_sizes.py` | 55 | True | Measure the size of every Open Targets output dataset from the EBI FTP directory indexes (read-only listing re |
| `catalogue/derive/pdf_text.py` | 39 | True | Store the text of a PDF as evidence (the probe collector decodes bodies as text and cannot read PDFs). |
| `catalogue/derive/vep_plugins.py` | 85 | True | List the Ensembl VEP plugins from the documentation page and keep what the page says about their data and lice |
| `catalogue/probe.py` | 128 | True | Collect provenance evidence for catalogue entries. |
| `catalogue/render.js` | 15 | n/a |  |
| `inventory/audit_existing.py` | 231 | True | Audit the existing dc assets under scrape/ and guard them against accidental change. |
| `inventory/check_repo_policy.py` | 99 | True | Repository policy check: nothing that must stay out of git is tracked, and the visibility is stated. |
| `inventory/repository_audit.py` | 416 | True | Repository audit: what exists, what state it is in, and whether the documentation matches the implementation. |
| `scrape/mavedb-human/build_index.py` | 41 | True | Build score_sets_index.csv from both detail files (organism-filter and accession-target selections). |
| `scrape/mavedb-human/download_scores.py` | 112 | True | Download MaveDB score (and count) tables for the human-target score sets. |
| `scrape/mavedb-human/verify_scores.py` | 84 | True | Verify scores/ against the manifest and the metadata files. |
