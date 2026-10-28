# Repository audit

Generated 2026-10-28T05:20:56Z at commit `fb4dc78e` by `inventory/repository_audit.py` (regenerate rather than edit). It records what exists, which lifecycle stage each dataset has reached, whether the checks pass, and whether the documentation says what the implementation does.

## Checks

| Check | Result | Detail |
|---|---|---|
| existing assets re-audited (NHGRI, Ensembl registry, MaveDB score tables) | PASS | {"nhgri-gwas": "downloaded_verified", "nhgri-1000-genomes": "downloaded_verified", "ensembl-human-genomes": "downloaded_verified", "mavedb-human": "downloaded_verified"} |
| preservation guard: tracked files under scrape/ unchanged apart from documented exceptions | PASS | baseline files: 97 | missing: 0 | changed: 0 | documented exceptions: 1 | added: 0 |
| unit tests (acquisition.test_helpers) | PASS | 17 tests |
| unit tests (variant_index) | PASS | 557 tests |
| catalogue records validate against the vocabulary, evidence files and acquisition targets | PASS | 54 source records valid; 173 acquisition targets; 277 evidence probes |
| repository policy: size, personal-genome file types, credential patterns | PASS | visibility: public (ziadwalieddine-create/dc) |
| all scripts compile | PASS | 149 scripts |
| the claim checker detects a changed number in every claim (mutation self-test) | PASS | 42 claims |

## Documentation claims checked against the implementation: 42 pass, 0 fail, 0 skipped

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
| `catalogue/MASTER_INVENTORY.md` | `(\d+) sources catalogued` | PASS | [54.0] |
| `latest:AUDIT_LOG.md` | ``validate\.py`: (\d+) of (\d+) targets pass` | PASS | [173.0, 173.0] |
| `latest:AUDIT_LOG.md` | `(\d+) manifest records` | PASS | [173.0] |
| `latest:AUDIT_LOG.md` | ``compat_tests\.py`: (\d+) tests, (\d+) PASS, (\d+) FAIL, (\d+) INFO` | PASS | [115.0, 101.0, 0.0, 14.0] |
| `latest:AUDIT_LOG.md` | `test_helpers`: (\d+) tests pass` | PASS | [17.0] |
| `latest:AUDIT_LOG.md` | ``variant_index`: (\d+) unit tests pass` | PASS | [557.0] |
| `variant_index/STATUS.md` | `(\d+) unit tests pass` | PASS | [557.0] |
| `variant_index/STATUS.md` | `Requirements: (\d+), of which (\d+) verified, (\d+) partial, (\d+) una` | PASS | [44.0, 30.0, 12.0, 1.0, 1.0] |
| `variant_index/STATUS.md` | `\| 60x \((\d+)\) \| (\d+) of (\d+) \(([\d.]+)%\) \| (\d+) of (\d+) \((` | PASS | [59.0, 781.0, 800.0, 97.6, 390.0, 400.0, 97.5, 996.0, 1000.0, 0.0] |
| `variant_index/STATUS.md` | `\| 30x \((\d+)\) \| (\d+) of (\d+) \(([\d.]+)%\) \| (\d+) of (\d+) \((` | PASS | [29.0, 769.0, 800.0, 96.1, 387.0, 400.0, 96.8, 994.0, 1000.0, 0.0] |
| `variant_index/STATUS.md` | `\| 15x \((\d+)\) \| (\d+) of (\d+) \(([\d.]+)%\) \| (\d+) of (\d+) \((` | PASS | [15.0, 696.0, 800.0, 87.0, 379.0, 400.0, 94.8, 895.0, 1000.0, 0.0] |
| `variant_index/STATUS.md` | `([\d,]+) traits named by effect-asserting evidence roll up into (\d+) ` | PASS | [10306.0, 25.0, 203.0] |
| `variant_index/STATUS.md` | `of the ([\d,]+) ids of the HG002 run: ([\d,]+) by their own term, (\d+` | PASS | [10306.0, 9943.0, 160.0, 8.0] |
| `variant_index/STATUS.md` | `([\d,]+) records re-opened in the acquired files, ([\d,]+) digests re-` | PASS | [24195946.0, 24195946.0] |
| `variant_index/STATUS.md` | `of the ([\d,]+) high-value loci that the benchmark callset leaves unre` | PASS | [4982.0, 1283.0, 218.0, 3405.0] |
| `variant_index/STATUS.md` | `([\d,]+) evidence rows from (\d+) source releases, ([\d,]+) variant en` | PASS | [24195946.0, 72.0, 6857645.0] |
| `catalogue/gaps.json` | `GTF: ([\d,]+) genes, ([\d,]+) transcripts, ([\d,]+) MANE Select` | PASS | [78941.0, 646577.0, 19299.0] |
| `catalogue/gaps.json` | `toplevel FASTA: (\d+) sequences` | PASS | [706.0] |
| `catalogue/gaps.json` | `\(([\d.]+)M PDB chain segments, ([\d,]+) PDB entries, ([\d,]+) UniProt` | PASS | [1.05, 246198.0, 76796.0] |
| `catalogue/gaps.json` | `\(([\d,]+) VCF files, ([\d,]+) variant records\)` | PASS | [7865.0, 35172.0] |
| `catalogue/gaps.json` | `with ([\d,]+) haplotype names` | PASS | [2458.0] |

## Branches

| Branch | Tip | Date | Tracked files | Ahead of / behind origin/main | Last commit |
|---|---|---|---|---|---|
| `ccr-b13fd75a-ra9gfy` | `fb4dc78` | 2026-10-28 | 676 | +88 / -0 | Every evidence row re-derived from its source; phase-6 results, audits and documents |
| `data/mavedb-scores` | `1b58196` | 2026-10-25 | 2061 | +8 / -22 | Add 52 superseded MaveDB score sets (79 tables, CC0) as a data extension |
| `main` | `f388b5f` | 2026-09-30 | 90 | +0 / -19 | Merge pull request #1 from ziadwalieddine-create/land/dc-extended-pack |
| `origin/ccr-b13fd75a-ra9gfy` | `fb4dc78` | 2026-10-28 | 676 | +88 / -0 | Every evidence row re-derived from its source; phase-6 results, audits and documents |
| `origin/data/mavedb-scores` | `1b58196` | 2026-10-25 | 2061 | +8 / -22 | Add 52 superseded MaveDB score sets (79 tables, CC0) as a data extension |
| `origin/land/dc-extended-pack` | `3547d0a` | 2026-09-30 | 90 | +0 / -20 | Land dc-extended skill pack |
| `origin/main` | `ddafb14` | 2026-10-25 | 98 | +0 / -0 | Merge pull request #6 from ziadwalieddine-create/ccr-b13fd75a-ra9gfy |

## Datasets (acquisition targets)

173 listed · 173 acquired (manifest ok and the file present with the recorded size) · 173 verified (validator passed) · 122 mapped (a passing join or consistency test used it, or the built index holds an evidence row of it mapped to a normalised entity) · 72 integrated (the built index holds a release of the target; `data/index` is git-ignored, so a checkout without it reports none).

| Target | Source record | Tier | Listed | Acquired | Verified | Mapped | Integrated | Size |
|---|---|---|---|---|---|---|---|---|
| `mavedb_zenodo_dump_v4` | `mavedb` | A | ✓ | ✓ | ✓ | ✓ | ✓ | 1907.7 MB |
| `clinvar_vcf_grch38` | `clinvar` | A | ✓ | ✓ | ✓ | ✓ | ✓ | 200.7 MB |
| `clinvar_readme` | `clinvar` | A | ✓ | ✓ | ✓ | · | · | 0.0 MB |
| `clinvar_readme_vcf` | `clinvar` | A | ✓ | ✓ | ✓ | · | · | 0.0 MB |
| `clingen_gene_curation_list` | `clingen` | A | ✓ | ✓ | ✓ | ✓ | ✓ | 0.2 MB |
| `clingen_region_curation_list` | `clingen` | A | ✓ | ✓ | ✓ | ✓ | ✓ | 0.1 MB |
| `clingen_haploinsufficiency_gene` | `clingen` | A | ✓ | ✓ | ✓ | ✓ | · | 0.0 MB |
| `clingen_triplosensitivity_gene` | `clingen` | A | ✓ | ✓ | ✓ | ✓ | · | 0.0 MB |
| `clingen_recurrent_cnv_v2_1` | `clingen` | A | ✓ | ✓ | ✓ | · | ✓ | 0.0 MB |
| `clingen_gene_validity` | `clingen` | A | ✓ | ✓ | ✓ | ✓ | ✓ | 1.1 MB |
| `gnomad_constraint_v4_1_2` | `gnomad` | A | ✓ | ✓ | ✓ | ✓ | ✓ | 574.5 MB |
| `gnomad_exome_cnv_v4_1` | `gnomad` | A | ✓ | ✓ | ✓ | ✓ | ✓ | 8.9 MB |
| `gwas_catalog_associations` | `gwas_catalog` | A | ✓ | ✓ | ✓ | ✓ | ✓ | 73.0 MB |
| `reactome_ensembl2reactome` | `reactome` | A | ✓ | ✓ | ✓ | ✓ | ✓ | 183.3 MB |
| `reactome_uniprot2reactome_all_levels` | `reactome` | A | ✓ | ✓ | ✓ | ✓ | · | 117.8 MB |
| `reactome_pathways` | `reactome` | A | ✓ | ✓ | ✓ | · | · | 1.6 MB |
| `clinpgx_variantAnnotations` | `pharmgkb` | A | ✓ | ✓ | ✓ | ✓ | ✓ | 4.3 MB |
| `clinpgx_relationships` | `pharmgkb` | A | ✓ | ✓ | ✓ | ✓ | · | 4.3 MB |
| `clinpgx_summaryAnnotations` | `pharmgkb` | A | ✓ | ✓ | ✓ | ✓ | ✓ | 1.2 MB |
| `clinpgx_guidelineAnnotations` | `pharmgkb` | A | ✓ | ✓ | ✓ | ✓ | ✓ | 0.9 MB |
| `clinpgx_drugLabels` | `pharmgkb` | A | ✓ | ✓ | ✓ | ✓ | ✓ | 0.1 MB |
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
| `ensembl_vep_mavedb_variants` | `ensembl_vep` | A | ✓ | ✓ | ✓ | ✓ | ✓ | 560.6 MB |
| `ensembl_vep_mavedb_variants_tbi` | `ensembl_vep` | A | ✓ | ✓ | ✓ | · | · | 0.0 MB |
| `strchive_hg38_general` | `strchive` | A | ✓ | ✓ | ✓ | ✓ | ✓ | 0.0 MB |
| `strchive_hg38_trgt` | `strchive` | A | ✓ | ✓ | ✓ | · | ✓ | 0.0 MB |
| `strchive_hg38_stranger_json` | `strchive` | A | ✓ | ✓ | ✓ | ✓ | ✓ | 0.0 MB |
| `dgv_grch38_variants_2025_12_01` | `dgv` | B | ✓ | ✓ | ✓ | ✓ | ✓ | 364.2 MB |
| `screen_ccre_pls` | `encode` | A | ✓ | ✓ | ✓ | ✓ | · | 2.6 MB |
| `screen_ccre_pels` | `encode` | A | ✓ | ✓ | ✓ | ✓ | · | 13.7 MB |
| `screen_ccre_dels` | `encode` | A | ✓ | ✓ | ✓ | ✓ | · | 80.7 MB |
| `screen_ccre_ca_ctcf` | `encode` | A | ✓ | ✓ | ✓ | ✓ | · | 7.3 MB |
| `screen_ccre_ca_h3k4me3` | `encode` | A | ✓ | ✓ | ✓ | ✓ | · | 4.8 MB |
| `screen_ccre_ca_tf` | `encode` | A | ✓ | ✓ | ✓ | ✓ | · | 1.5 MB |
| `screen_ccre_ca` | `encode` | A | ✓ | ✓ | ✓ | ✓ | · | 13.0 MB |
| `screen_ccre_tf` | `encode` | A | ✓ | ✓ | ✓ | ✓ | · | 5.6 MB |
| `opentargets_26_09_target` | `opentargets` | A | ✓ | ✓ | ✓ | ✓ | ✓ | 21.6 MB |
| `opentargets_26_09_go` | `opentargets` | A | ✓ | ✓ | ✓ | · | · | 2.9 MB |
| `opentargets_26_09_so` | `opentargets` | A | ✓ | ✓ | ✓ | · | · | 0.0 MB |
| `opentargets_26_09_homology` | `opentargets` | A | ✓ | ✓ | ✓ | ✓ | · | 24.5 MB |
| `opentargets_26_09_biosample` | `opentargets` | A | ✓ | ✓ | ✓ | · | ✓ | 6.0 MB |
| `opentargets_26_09_disease` | `opentargets` | A | ✓ | ✓ | ✓ | ✓ | · | 6.3 MB |
| `opentargets_26_09_disease_hpo` | `opentargets` | A | ✓ | ✓ | ✓ | · | · | 1.7 MB |
| `opentargets_26_09_disease_phenotype` | `opentargets` | A | ✓ | ✓ | ✓ | · | · | 6.5 MB |
| `opentargets_26_09_mouse_phenotype` | `opentargets` | A | ✓ | ✓ | ✓ | ✓ | · | 6.8 MB |
| `opentargets_26_09_evidence_genomics_england` | `opentargets` | A | ✓ | ✓ | ✓ | ✓ | · | 2.7 MB |
| `opentargets_26_09_evidence_gene2phenotype` | `opentargets` | A | ✓ | ✓ | ✓ | ✓ | · | 0.4 MB |
| `opentargets_26_09_evidence_orphanet` | `opentargets` | A | ✓ | ✓ | ✓ | ✓ | · | 0.5 MB |
| `opentargets_26_09_evidence_uniprot_variants` | `opentargets` | A | ✓ | ✓ | ✓ | ✓ | · | 1.9 MB |
| `opentargets_26_09_pharmacogenomics` | `opentargets` | A | ✓ | ✓ | ✓ | ✓ | · | 3.4 MB |
| `opentargets_26_09_evidence_gene_burden` | `opentargets` | A | ✓ | ✓ | ✓ | ✓ | ✓ | 3.9 MB |
| `opentargets_26_09_target_essentiality` | `opentargets` | A | ✓ | ✓ | ✓ | ✓ | ✓ | 246.6 MB |
| `opentargets_26_09_evidence_crispr` | `opentargets` | A | ✓ | ✓ | ✓ | ✓ | ✓ | 0.0 MB |
| `opentargets_26_09_evidence_crispr_screen` | `opentargets` | A | ✓ | ✓ | ✓ | ✓ | ✓ | 0.9 MB |
| `opentargets_26_09_evidence_impc` | `opentargets` | A | ✓ | ✓ | ✓ | ✓ | ✓ | 906.9 MB |
| `opentargets_26_09_study` | `opentargets` | A | ✓ | ✓ | ✓ | ✓ | ✓ | 97.7 MB |
| `opentargets_26_09_credible_set` | `opentargets` | A | ✓ | ✓ | ✓ | ✓ | ✓ | 2990.6 MB |
| `opentargets_26_09_l2g_prediction` | `opentargets` | A | ✓ | ✓ | ✓ | ✓ | ✓ | 509.6 MB |
| `opentargets_26_09_variant` | `opentargets` | A | ✓ | ✓ | ✓ | ✓ | · | 1967.3 MB |
| `opentargets_26_09_enhancer_to_gene` | `opentargets` | A | ✓ | ✓ | ✓ | ✓ | · | 1770.1 MB |
| `opentargets_26_09_interaction` | `opentargets` | A | ✓ | ✓ | ✓ | ✓ | · | 251.9 MB |
| `opentargets_26_09_interaction_evidence` | `opentargets` | A | ✓ | ✓ | ✓ | ✓ | · | 344.1 MB |
| `opentargets_26_09_evidence_uniprot_literature` | `opentargets` | A | ✓ | ✓ | ✓ | ✓ | · | 0.5 MB |
| `gnomad_sv_v4_1_sites_vcf` | `gnomad` | B | ✓ | ✓ | ✓ | ✓ | ✓ | 1743.5 MB |
| `alphafold_db_mavedb_proteins` | `alphafold_db` | A | ✓ | ✓ | ✓ | ✓ | · | 356.7 MB |
| `interpro_mavedb_proteins` | `interpro` | A | ✓ | ✓ | ✓ | ✓ | · | 8.0 MB |
| `gnomad_v4_1_exomes_regional` | `gnomad` | A | ✓ | ✓ | ✓ | ✓ | ✓ | 7.9 MB |
| `gnomad_v4_1_genomes_regional` | `gnomad` | A | ✓ | ✓ | ✓ | ✓ | · | 2.4 MB |
| `uniprot_humsavar` | `uniprot` | A | ✓ | ✓ | ✓ | ✓ | · | 8.7 MB |
| `gnomad_v4_1_joint_regional` | `gnomad` | A | ✓ | ✓ | ✓ | ✓ | ✓ | 4.3 MB |
| `ensembl_grch38_reference_regional` | `ensembl_genomes` | A | ✓ | ✓ | ✓ | ✓ | · | 0.1 MB |
| `ensembl_grch38_dna_toplevel` | `ensembl_genomes` | B | ✓ | ✓ | ✓ | ✓ | · | 847.7 MB |
| `ensembl_grch38_dna_toplevel_fai` | `ensembl_genomes` | B | ✓ | ✓ | ✓ | ✓ | · | 0.0 MB |
| `ensembl_grch38_dna_toplevel_gzi` | `ensembl_genomes` | B | ✓ | ✓ | ✓ | · | · | 0.8 MB |
| `ensembl_grch38_gtf_116` | `ensembl_genomes` | B | ✓ | ✓ | ✓ | ✓ | · | 141.1 MB |
| `sifts_pdb_chain_uniprot` | `pdbe_sifts` | A | ✓ | ✓ | ✓ | ✓ | · | 6.3 MB |
| `sifts_uniprot_segments_observed` | `pdbe_sifts` | A | ✓ | ✓ | ✓ | ✓ | · | 11.0 MB |
| `gnomad_v4_1_lof_curation_incomplete_penetrance` | `gnomad_constraint` | A | ✓ | ✓ | ✓ | ✓ | ✓ | 0.1 MB |
| `pharmvar_6_2_29` | `pharmvar` | B | ✓ | ✓ | ✓ | ✓ | ✓ | 193.1 MB |
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
| `giab_v3_6_grch38_alldifficultregions` | `giab_stratifications` | A | ✓ | ✓ | ✓ | ✓ | ✓ | 31.5 MB |
| `giab_v3_6_grch38_alllowmapandsegdupregions` | `giab_stratifications` | A | ✓ | ✓ | ✓ | ✓ | ✓ | 3.5 MB |
| `giab_v3_6_grch38_lowmappabilityall` | `giab_stratifications` | A | ✓ | ✓ | ✓ | ✓ | ✓ | 4.5 MB |
| `giab_v3_6_grch38_segdups` | `giab_stratifications` | A | ✓ | ✓ | ✓ | ✓ | ✓ | 0.1 MB |
| `giab_v3_6_grch38_segdups_gt10kb` | `giab_stratifications` | A | ✓ | ✓ | ✓ | ✓ | ✓ | 0.0 MB |
| `giab_v3_6_grch38_mhc` | `giab_stratifications` | A | ✓ | ✓ | ✓ | ✓ | ✓ | 0.0 MB |
| `giab_v3_6_grch38_kir` | `giab_stratifications` | A | ✓ | ✓ | ✓ | ✓ | ✓ | 0.0 MB |
| `giab_v3_6_grch38_vdj` | `giab_stratifications` | A | ✓ | ✓ | ✓ | ✓ | ✓ | 0.0 MB |
| `giab_v3_6_grch38_false_duplications_correct_copy` | `giab_stratifications` | A | ✓ | ✓ | ✓ | ✓ | ✓ | 0.0 MB |
| `giab_v3_6_grch38_false_duplications_incorrect_copy` | `giab_stratifications` | A | ✓ | ✓ | ✓ | ✓ | ✓ | 0.0 MB |
| `giab_v3_6_grch38_collapsed_duplication_fp_regions` | `giab_stratifications` | A | ✓ | ✓ | ✓ | ✓ | ✓ | 0.0 MB |
| `giab_v3_6_grch38_population_cnv_fp_regions` | `giab_stratifications` | A | ✓ | ✓ | ✓ | ✓ | ✓ | 0.0 MB |
| `giab_v3_6_grch38_l1h_gt500` | `giab_stratifications` | A | ✓ | ✓ | ✓ | ✓ | ✓ | 0.0 MB |
| `giab_v3_6_grch38_gaps_slop15kb` | `giab_stratifications` | A | ✓ | ✓ | ✓ | ✓ | ✓ | 0.0 MB |
| `giab_v3_6_grch38_allotherdifficultregions` | `giab_stratifications` | A | ✓ | ✓ | ✓ | ✓ | ✓ | 0.3 MB |
| `giab_v3_6_grch38_cmrgv1_00_falselyduplicatedgenes` | `giab_stratifications` | A | ✓ | ✓ | ✓ | ✓ | ✓ | 0.0 MB |
| `giab_v3_6_grch38_badpromoters` | `giab_stratifications` | A | ✓ | ✓ | ✓ | ✓ | ✓ | 0.0 MB |
| `giab_v3_6_grch38_alltandemrepeats` | `giab_stratifications` | A | ✓ | ✓ | ✓ | ✓ | ✓ | 9.2 MB |
| `giab_v3_6_grch38_satellites_slop5` | `giab_stratifications` | A | ✓ | ✓ | ✓ | ✓ | ✓ | 0.0 MB |
| `giab_v3_6_grch38_chrx_par` | `giab_stratifications` | A | ✓ | ✓ | ✓ | ✓ | ✓ | 0.0 MB |
| `giab_v3_6_grch38_chry_par` | `giab_stratifications` | A | ✓ | ✓ | ✓ | ✓ | ✓ | 0.0 MB |
| `giab_v3_6_grch38_chrx_xtr` | `giab_stratifications` | A | ✓ | ✓ | ✓ | ✓ | ✓ | 0.0 MB |
| `giab_v3_6_grch38_chrx_ampliconic` | `giab_stratifications` | A | ✓ | ✓ | ✓ | ✓ | ✓ | 0.0 MB |
| `proteingym_v1_3_dms_substitutions` | `proteingym` | A | ✓ | ✓ | ✓ | ✓ | · | 43.0 MB |
| `proteingym_v1_3_dms_indels` | `proteingym` | A | ✓ | ✓ | ✓ | · | · | 6.7 MB |
| `proteingym_v1_3_dms_substitutions_reference` | `proteingym` | A | ✓ | ✓ | ✓ | ✓ | · | 0.2 MB |
| `proteingym_v1_3_dms_indels_reference` | `proteingym` | A | ✓ | ✓ | ✓ | ✓ | · | 0.0 MB |
| `hpa_v25_1_rna_tissue_consensus` | `hpa` | A | ✓ | ✓ | ✓ | ✓ | · | 5.3 MB |
| `hpa_v25_1_rna_single_cell_type` | `hpa` | A | ✓ | ✓ | ✓ | ✓ | · | 16.3 MB |
| `hpa_v25_1_subcellular_location` | `hpa` | A | ✓ | ✓ | ✓ | ✓ | · | 0.3 MB |
| `biogrid_interactions_organism_tab3` | `biogrid_interactions` | A | ✓ | ✓ | ✓ | ✓ | · | 188.9 MB |
| `hgsvc3_grch38_sv_insdel_sym_vcf` | `hgsvc3` | A | ✓ | ✓ | ✓ | ✓ | ✓ | 6.4 MB |
| `hgsvc3_grch38_sv_inv_sym_vcf` | `hgsvc3` | A | ✓ | ✓ | ✓ | ✓ | ✓ | 0.0 MB |
| `hgsvc3_variant_calls_readme` | `hgsvc3` | A | ✓ | ✓ | ✓ | · | · | 0.0 MB |
| `hgsvc3_variant_calls_manifest` | `hgsvc3` | A | ✓ | ✓ | ✓ | · | · | 0.0 MB |
| `cpic_api_allele_definition` | `cpic` | A | ✓ | ✓ | ✓ | ✓ | ✓ | 0.2 MB |
| `cpic_api_allele_location_value` | `cpic` | A | ✓ | ✓ | ✓ | · | · | 0.3 MB |
| `cpic_api_sequence_location` | `cpic` | A | ✓ | ✓ | ✓ | ✓ | · | 0.3 MB |
| `cpic_api_allele` | `cpic` | A | ✓ | ✓ | ✓ | ✓ | ✓ | 1.1 MB |
| `cpic_api_gene_result` | `cpic` | A | ✓ | ✓ | ✓ | · | · | 0.1 MB |
| `cpic_api_gene_result_lookup` | `cpic` | A | ✓ | ✓ | ✓ | · | · | 0.1 MB |
| `cpic_api_gene_result_diplotype` | `cpic` | A | ✓ | ✓ | ✓ | · | · | 22.2 MB |
| `cpic_api_recommendation` | `cpic` | A | ✓ | ✓ | ✓ | ✓ | ✓ | 2.4 MB |
| `impc_r24_readme` | `impc` | A | ✓ | ✓ | ✓ | · | · | 0.0 MB |
| `impc_r24_genotype_phenotype_assertions_all` | `impc` | A | ✓ | ✓ | ✓ | ✓ | · | 5.0 MB |
| `impc_r24_gene_mp_association_1` | `impc` | A | ✓ | ✓ | ✓ | · | · | 0.2 MB |
| `impc_r24_viability` | `impc` | A | ✓ | ✓ | ✓ | ✓ | · | 0.5 MB |
| `clingen_actionability_adult_summary_flat` | `clingen` | A | ✓ | ✓ | ✓ | · | · | 0.1 MB |
| `clingen_actionability_adult_assertions` | `clingen` | A | ✓ | ✓ | ✓ | ✓ | ✓ | 1.1 MB |
| `clingen_actionability_adult_scoring` | `clingen` | A | ✓ | ✓ | ✓ | ✓ | ✓ | 1.7 MB |
| `clingen_actionability_pediatric_summary_flat` | `clingen` | A | ✓ | ✓ | ✓ | · | · | 0.1 MB |
| `clingen_actionability_pediatric_assertions` | `clingen` | A | ✓ | ✓ | ✓ | ✓ | ✓ | 1.2 MB |
| `clingen_actionability_pediatric_scoring` | `clingen` | A | ✓ | ✓ | ✓ | ✓ | ✓ | 1.7 MB |
| `clingen_erepo_variant_classifications` | `clingen` | A | ✓ | ✓ | ✓ | ✓ | ✓ | 32.8 MB |
| `clinpgx_clinicalAnnotations` | `pharmgkb` | A | ✓ | ✓ | ✓ | · | · | 1.2 MB |
| `clinpgx_clinicalVariants` | `pharmgkb` | A | ✓ | ✓ | ✓ | · | · | 0.1 MB |
| `clinpgx_variants` | `pharmgkb` | A | ✓ | ✓ | ✓ | · | ✓ | 0.9 MB |
| `clinpgx_genes` | `pharmgkb` | A | ✓ | ✓ | ✓ | · | · | 2.9 MB |
| `gwas_catalog_studies` | `gwas_catalog` | A | ✓ | ✓ | ✓ | · | ✓ | 130.3 MB |
| `gwas_catalog_ancestries` | `gwas_catalog` | A | ✓ | ✓ | ✓ | · | ✓ | 53.8 MB |
| `gwas_catalog_efo_trait_mappings` | `gwas_catalog` | A | ✓ | ✓ | ✓ | · | · | 19.5 MB |
| `pgs_scores_list` | `pgs_catalog` | A | ✓ | ✓ | ✓ | · | · | 0.1 MB |
| `pgs_metadata_scores` | `pgs_catalog` | A | ✓ | ✓ | ✓ | · | · | 4.3 MB |
| `pgs_metadata_publications` | `pgs_catalog` | A | ✓ | ✓ | ✓ | · | · | 0.5 MB |
| `pgs_metadata_efo_traits` | `pgs_catalog` | A | ✓ | ✓ | ✓ | · | · | 0.2 MB |
| `pgs_metadata_performance_metrics` | `pgs_catalog` | A | ✓ | ✓ | ✓ | · | · | 5.5 MB |
| `pgs_metadata_evaluation_sample_sets` | `pgs_catalog` | A | ✓ | ✓ | ✓ | · | · | 1.7 MB |
| `pgs_metadata_score_development_samples` | `pgs_catalog` | A | ✓ | ✓ | ✓ | · | · | 2.7 MB |
| `pgs_metadata_cohorts` | `pgs_catalog` | A | ✓ | ✓ | ✓ | · | · | 0.1 MB |
| `clinvar_variant_summary` | `clinvar` | A | ✓ | ✓ | ✓ | · | · | 456.5 MB |
| `clinvar_submission_summary` | `clinvar` | A | ✓ | ✓ | ✓ | · | · | 404.8 MB |
| `pgs_scoring_files_grch38` | `pgs_catalog` | A | ✓ | ✓ | ✓ | · | · | 851.8 MB |
| `mavedb_scores_genomic_mapped` | `mavedb` | A | ✓ | ✓ | ✓ | ✓ | ✓ | 28.8 MB |
| `clinvar_vcv_xml_weekly_2026_1018` | `clinvar` | A | ✓ | ✓ | ✓ | ✓ | ✓ | 6071.4 MB |

## Manifests

| File | Records | What it records |
|---|---|---|
| `acquisition/manifest.jsonl` | 173 | 173 ok |
| `acquisition/manifest_files/alphafold_db_mavedb_proteins.tsv` | 1621 | per-file table |
| `acquisition/manifest_files/ensembl_grch38_reference_regional.tsv` | 3 | per-file table |
| `acquisition/manifest_files/gnomad_v4_1_exomes_regional.tsv` | 24 | per-file table |
| `acquisition/manifest_files/gnomad_v4_1_genomes_regional.tsv` | 24 | per-file table |
| `acquisition/manifest_files/gnomad_v4_1_joint_regional.tsv` | 24 | per-file table |
| `acquisition/manifest_files/interpro_mavedb_proteins.tsv` | 549 | per-file table |
| `acquisition/manifest_files/mavedb_scores_genomic_mapped.tsv` | 200 | per-file table |
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
| `acquisition/manifest_files/pgs_scoring_files_grch38.tsv` | 3467 | per-file table |
| `scrape/mavedb-human/scores_manifest.jsonl` | 1980 | sha256 and row count of each score and count table |
| `inventory/baseline_checksums.sha256` | 97 | sha256 of every tracked file under scrape/ at the baseline |
| `catalogue/evidence/evidence_index.jsonl` | 321 | append-only record of every probe (URL, time, status, sha256) |

## Documents

| File | Title | Lines | Last commit |
|---|---|---|---|
| `AUDIT_LOG.md` | Audit log | 320 | 2026-10-28 |
| `README.md` | dc: human genome reference data | 26 | 2026-10-26 |
| `acquisition/README.md` | Acquisition, validation and compatibility tests | 46 | 2026-10-26 |
| `catalogue/ACQUIRED_FILES.md` | Acquired files register | 181 | 2026-10-26 |
| `catalogue/ATTRIBUTION.md` | Attribution and licence notes for the data we hold | 43 | 2026-10-26 |
| `catalogue/COVERAGE_MATRIX.md` | Coverage matrix: evidence dimensions by lifecycle stage | 230 | 2026-10-26 |
| `catalogue/EVIDENCE_PROFILES.md` | Evidence profiles | 871 | 2026-10-26 |
| `catalogue/GAP_MATRIX.md` | Gap matrix | 66 | 2026-10-26 |
| `catalogue/GOAL_RESOURCES.md` | Resources named in the project goal | 116 | 2026-10-26 |
| `catalogue/MASTER_INVENTORY.md` | Master inventory of reusable genomic resources | 150 | 2026-10-26 |
| `catalogue/MECHANISM_MATRIX.md` | Mechanism matrix: variant types by mechanism axis | 323 | 2026-10-26 |
| `catalogue/README.md` | Evidence-source catalogue | 47 | 2026-10-26 |
| `catalogue/REUSE_PLAN.md` | Source-by-source reuse, acquisition and integration plan | 356 | 2026-10-26 |
| `inventory/EXISTING_ASSETS.md` | Existing assets: inventory and verified status | 49 | 2026-10-28 |
| `inventory/PHASE_AUDIT.md` | Phase gap audit | 38 | 2026-10-28 |
| `inventory/REPOSITORY_AUDIT.md` | Repository audit | 483 | 2026-10-28 |
| `scrape/ensembl-human-genomes/README.md` | Ensembl genome selector: human (taxonomy ID 9606) | 20 | 2026-10-25 |
| `scrape/mavedb-human/README.md` | MaveDB: human score sets | 45 | 2026-10-25 |
| `scrape/nhgri-1000-genomes/README.md` | NHGRI: 1000 Genomes Project (archived page scrape) | 22 | 2026-10-25 |
| `scrape/nhgri-gwas/README.md` | NHGRI Talking Glossary: GWAS (scrape) | 17 | 2026-10-25 |
| `variant_index/BUILD_REPORT.md` | Build report | 611 | 2026-10-28 |
| `variant_index/README.md` | High-Value Human Variant Reference Index (`vindex`) | 79 | 2026-10-28 |
| `variant_index/REQUIREMENTS.md` | Requirements and verification | 178 | 2026-10-28 |
| `variant_index/SCHEMA.md` | Index schema (version 0.4.0) | 297 | 2026-10-27 |
| `variant_index/STATUS.md` | Status of the High-Value Human Variant Reference Index | 133 | 2026-10-28 |

## Scripts

| File | Lines | Compiles | Purpose |
|---|---|---|---|
| `acquisition/acquire.py` | 325 | True | Download external resources listed in targets.json, verify them and record a manifest. |
| `acquisition/compat_opentargets.py` | 272 | True | Joins between the acquired Open Targets 26.09 Parquet datasets and the MaveDB human set (called from compat_te |
| `acquisition/compat_reference.py` | 405 | True | Tests that tie the mirrored GRCh38 reference sequence and the Ensembl release-116 gene annotation to the other |
| `acquisition/compat_resources.py` | 865 | True | Tests that tie the RefSeq accession map and annotation, the GIAB difficult-region stratifications and ProteinG |
| `acquisition/compat_tests.py` | 859 | True | Compatibility tests: do the acquired resources join to the existing dc assets (MaveDB human set)? |
| `acquisition/extract_dump.py` | 57 | True | Extract, from the MaveDB bulk dump v4, the members that the compatibility tests and the superseded-set builder |
| `acquisition/extract_mavedb_scores.py` | 95 | True | Extract the score tables of the MaveDB score sets that have genomic (GRCh38) mappings from the bulk dump. |
| `acquisition/fetch_gnomad_regional.py` | 213 | True | Slice the gnomAD v4.1 exome and genome sites VCFs to the positions of MaveDB-mapped SNVs, by remote tabix. |
| `acquisition/fetch_per_protein.py` | 178 | True | Fetch per-protein resources for the UniProt proteins that the MaveDB human score sets are about. |
| `acquisition/fetch_pgs_scores.py` | 199 | True | Fetch PGS Catalog harmonised GRCh38 scoring files for a documented, reproducible selection. |
| `acquisition/fetch_reference_regional.py` | 155 | True | Read the GRCh38 reference sequence for the regions that contain MaveDB-mapped SNV positions, by remote random  |
| `acquisition/intervals.py` | 20 | True | Exact stabbing queries over half-open intervals, used by every genomic-overlap test in compat_tests.py. |
| `acquisition/mavedb_snvs.py` | 77 | True | Genomic SNV positions of the MaveDB human set, from the two mapping products that exist (see compat_tests.py). |
| `acquisition/mavedb_superseded/build_superseded.py` | 84 | True | Add the MaveDB score sets that the public search hides because they are superseded versions. |
| `acquisition/stream_clinvar_vcv.py` | 357 | True | Stream the ClinVar VCV XML full release into Parquet tables without storing the 6 GB archive. |
| `acquisition/test_helpers.py` | 165 | True | Unit tests for the small functions that every overlap and join number in this directory rests on. |
| `acquisition/validate.py` | 2186 | True | Validate downloaded external files against each source's own description. |
| `catalogue/build_catalogue.py` | 572 | True | Validate catalogue/sources/*.json and generate the inventory, coverage matrix, gap matrix and reuse plan. |
| `catalogue/derive/evidence_support.py` | 222 | True | Find statements in the source records that none of their stored evidence supports. |
| `catalogue/derive/gnomad_sv_header.py` | 38 | True | Record the INFO and FILTER definitions of the gnomAD v4.1 SV sites VCF, read from the header of the downloaded |
| `catalogue/derive/gnomad_vcf_header.py` | 44 | True | Record the INFO/FILTER definitions of the gnomAD v4.1 exome and genome sites VCFs, read from the VCF headers b |
| `catalogue/derive/local_file.py` | 50 | True | Store a local file (for example a licence read from a sanctioned read-only clone) as evidence, with its proven |
| `catalogue/derive/openebench_counts.py` | 69 | True | Count registry tools relevant to human variant interpretation, from the OpenEBench monitor snapshot. |
| `catalogue/derive/opentargets_datasets.py` | 40 | True | Extract the dataset inventory (name, category, one-line description) from the rendered Open Targets downloads  |
| `catalogue/derive/opentargets_sizes.py` | 55 | True | Measure the size of every Open Targets output dataset from the EBI FTP directory indexes (read-only listing re |
| `catalogue/derive/ot_credible_set_counts.py` | 41 | True | Count the acquired Open Targets credible sets and studies by type and project. |
| `catalogue/derive/ot_interaction_counts.py` | 32 | True | Count the rows and distinct gene pairs per source database in the acquired Open Targets interaction dataset. |
| `catalogue/derive/pdf_text.py` | 39 | True | Store the text of a PDF as evidence (the probe collector decodes bodies as text and cannot read PDFs). |
| `catalogue/derive/vep_plugins.py` | 85 | True | List the Ensembl VEP plugins from the documentation page and keep what the page says about their data and lice |
| `catalogue/probe.py` | 128 | True | Collect provenance evidence for catalogue entries. |
| `catalogue/render.js` | 15 | n/a |  |
| `inventory/audit_existing.py` | 231 | True | Audit the existing dc assets under scrape/ and guard them against accidental change. |
| `inventory/check_repo_policy.py` | 99 | True | Repository policy check: nothing that must stay out of git is tracked, and the visibility is stated. |
| `inventory/gap_audit.py` | 552 | True | Adversarial gap audit of the phase requirements: every check either passes, fails, is skipped for a stated mis |
| `inventory/repository_audit.py` | 550 | True | Repository audit: what exists, what state it is in, and whether the documentation matches the implementation. |
| `scrape/mavedb-human/build_index.py` | 41 | True | Build score_sets_index.csv from both detail files (organism-filter and accession-target selections). |
| `scrape/mavedb-human/download_scores.py` | 112 | True | Download MaveDB score (and count) tables for the human-target score sets. |
| `scrape/mavedb-human/verify_scores.py` | 84 | True | Verify scores/ against the manifest and the metadata files. |
| `variant_index/build_index.py` | 97 | True | Build the High-Value Human Variant Reference Index from the acquired sources. |
| `variant_index/build_pgs_registry.py` | 32 | True | Build the PGS Catalog model registry (data/pgs_registry) from the acquired bulk metadata and the scoring-file  |
| `variant_index/check_idempotence.py` | 92 | True | Build the index twice from the same sources and code and compare the tables byte for byte. |
| `variant_index/check_requirements.py` | 97 | True | Check requirements.json against the repository and render REQUIREMENTS.md. |
| `variant_index/compare_runs.py` | 63 | True | Compare two sample runs of the matcher: the same input and index must give the same tables whatever the hash s |
| `variant_index/crosscheck_bcftools.py` | 71 | True | Cross-check a whole-genome run against bcftools, an independent normaliser. |
| `variant_index/fetch_public_testset.py` | 79 | True | Fetch the public Genome in a Bottle HG002 GRCh38 benchmark call set, a real, openly released callset used to t |
| `variant_index/make_reports.py` | 108 | True | Coverage statistics, validation, target lifecycle, change log and gap report for a built index. |
| `variant_index/make_status.py` | 53 | True | Keep the module table of STATUS.md equal to the module registry (vindex.reports.MODULES). |
| `variant_index/make_views.py` | 41 | True | Priority lists, trait groupings and pathway groupings over the matches of one sample (see vindex/views.py for  |
| `variant_index/match_wgs.py` | 80 | True | Match one sample's GRCh38 calls against the variant reference index, one chromosome at a time. |
| `variant_index/score_pgs.py` | 64 | True | Score a sample with PGS Catalog models, one chromosome at a time, with the registry's gates beside every resul |
| `variant_index/tests/bamfix.py` | 97 | True | Synthetic alignments for the read-evidence tests: reads are cut from the synthetic reference with the edits a  |
| `variant_index/tests/evidence_fixtures.py` | 33 | True | Hand-made evidence rows for matcher tests. Every row is validated against the schema vocabulary. |
| `variant_index/tests/fixtures/make_reference.py` | 49 | True | Deterministic synthetic reference ('SYN38') for the tests: small contigs with engineered repeats, so that left |
| `variant_index/tests/fixtures/make_sources.py` | 879 | True | Mini source files in the real formats, built on the synthetic reference, for the source -> index -> matcher te |
| `variant_index/tests/helpers.py` | 80 | True | Shared test helpers: synthetic reference, VCF writer, bcftools discovery. |
| `variant_index/tests/test_alleles.py` | 73 | True |  |
| `variant_index/tests/test_build_clingen.py` | 135 | True |  |
| `variant_index/tests/test_build_clinpgx.py` | 283 | True |  |
| `variant_index/tests/test_build_clinvar.py` | 225 | True |  |
| `variant_index/tests/test_build_clinvar_vcv.py` | 279 | True |  |
| `variant_index/tests/test_build_events.py` | 107 | True |  |
| `variant_index/tests/test_build_gnomad.py` | 111 | True |  |
| `variant_index/tests/test_build_gwas.py` | 330 | True |  |
| `variant_index/tests/test_build_mavedb.py` | 170 | True |  |
| `variant_index/tests/test_build_opentargets.py` | 168 | True |  |
| `variant_index/tests/test_build_opentargets_gene.py` | 130 | True |  |
| `variant_index/tests/test_build_pathways.py` | 39 | True |  |
| `variant_index/tests/test_build_pgx.py` | 279 | True |  |
| `variant_index/tests/test_build_regions.py` | 134 | True |  |
| `variant_index/tests/test_classes.py` | 66 | True |  |
| `variant_index/tests/test_compare_runs.py` | 91 | True | compare_runs: two runs of one input are identical whatever the moment; a changed table, a changed manifest fie |
| `variant_index/tests/test_crosscheck.py` | 223 | True | The cross-check recomputes calls and evidence matches with bcftools only. It must agree with a clean run and m |
| `variant_index/tests/test_effect_alleles.py` | 122 | True | An association is about one allele of a site. When that allele is the reference allele the sample is linked th |
| `variant_index/tests/test_gap_audit.py` | 200 | True |  |
| `variant_index/tests/test_genes.py` | 90 | True |  |
| `variant_index/tests/test_idempotence.py` | 78 | True |  |
| `variant_index/tests/test_intervals.py` | 78 | True |  |
| `variant_index/tests/test_known_loci.py` | 242 | True | The known-locus engine on a small in-memory index: a fact that holds passes, and every way the index can contr |
| `variant_index/tests/test_make_status.py` | 30 | True |  |
| `variant_index/tests/test_match_hap.py` | 94 | True |  |
| `variant_index/tests/test_match_repeat.py` | 63 | True |  |
| `variant_index/tests/test_match_small.py` | 268 | True |  |
| `variant_index/tests/test_match_sv.py` | 154 | True |  |
| `variant_index/tests/test_normalise.py` | 223 | True |  |
| `variant_index/tests/test_odd_records.py` | 177 | True | Records as real callsets hold them: lowercase and N bases, a REF that is not the reference, '*' alleles, gVCF  |
| `variant_index/tests/test_other_contigs.py` | 160 | True | Contigs outside the primary assembly (alternate loci, unlocalised scaffolds, decoys) are indexed, matched and  |
| `variant_index/tests/test_pgs.py` | 132 | True |  |
| `variant_index/tests/test_pgs_registry.py` | 114 | True |  |
| `variant_index/tests/test_pgs_run.py` | 126 | True |  |
| `variant_index/tests/test_provenance.py` | 284 | True |  |
| `variant_index/tests/test_reads.py` | 485 | True | Read-level evidence at a locus: allele counting from alignments independent of how the aligner wrote the event |
| `variant_index/tests/test_reads_run.py` | 178 | True | Read evidence through a whole run: the matcher's status of an unresolved locus is never changed by the reads,  |
| `variant_index/tests/test_reference_effect.py` | 252 | True | Evidence whose effect allele is the reference allele (a GWAS risk allele is often the reference allele) is rea |
| `variant_index/tests/test_reports.py` | 376 | True |  |
| `variant_index/tests/test_requirements.py` | 50 | True |  |
| `variant_index/tests/test_resolvers.py` | 136 | True | Every provenance resolver reproduces, for EVERY evidence row of its release, the digest the adapter stored: th |
| `variant_index/tests/test_rollup.py` | 113 | True | Trait roll-up through the cross-references and retired terms of the ontology table: ids the table does not hol |
| `variant_index/tests/test_runner_guards.py` | 227 | True |  |
| `variant_index/tests/test_scale_paths.py` | 128 | True | The paths that keep a whole-genome match within memory must give the same answers as the plain ones: per-conti |
| `variant_index/tests/test_summary.py` | 38 | True |  |
| `variant_index/tests/test_vcfio.py` | 162 | True |  |
| `variant_index/tests/test_wgs.py` | 405 | True | End to end on synthetic data: sources -> index -> a sample's calls -> per-chromosome matching -> result tables |
| `variant_index/validate_reads.py` | 292 | True | Check the read-evidence resolver (vindex/reads.py) against a benchmark on real alignments. |
| `variant_index/verify_index.py` | 39 | True | Check an index's provenance: release linkage, release records against the acquisition manifest, and raw-record |
| `variant_index/verify_known_loci.py` | 43 | True | Check a built index against facts that were established outside it (vindex/known_loci.py, known_loci.json). |
| `variant_index/vindex/__init__.py` | 11 | True | High-Value Human Variant Reference Index (vindex). |
| `variant_index/vindex/build/__init__.py` | 1 | True | Index build framework: adapters turn acquired source files into evidence rows and entities. |
| `variant_index/vindex/build/alleles.py` | 134 | True | Orientation of a reported effect allele against the reference: the part of harmonisation that decides whether  |
| `variant_index/vindex/build/classes.py` | 83 | True | ClinVar classification vocabulary, shared by the VCF adapter (aggregate strings such as 'Likely_pathogenic') a |
| `variant_index/vindex/build/clingen.py` | 285 | True | ClinGen -> evidence rows, each at its own level. |
| `variant_index/vindex/build/clinpgx.py` | 529 | True | ClinPGx (formerly PharmGKB) -> evidence rows. |
| `variant_index/vindex/build/clinvar_vcf.py` | 179 | True | ClinVar GRCh38 VCF -> evidence rows (variant level) with cross-references and gene relations. |
| `variant_index/vindex/build/clinvar_vcv.py` | 338 | True | ClinVar full VCV release (stream-derived Parquet tables) -> submission-level evidence, structural events and c |
| `variant_index/vindex/build/common.py` | 288 | True | Build context, row sinks and evidence-row construction shared by all adapters. |
| `variant_index/vindex/build/events.py` | 244 | True | Structural-variant and copy-number events observed in population resources (GRCh38) -> entity_event rows. |
| `variant_index/vindex/build/genes.py` | 132 | True | Gene identity: Ensembl gene id (version stripped) is the entity id; HGNC id, NCBI Gene id and symbols are cros |
| `variant_index/vindex/build/gnomad.py` | 207 | True | gnomAD v4.1 -> evidence rows. |
| `variant_index/vindex/build/gwas.py` | 310 | True | GWAS Catalog associations (release 2026-10-22) -> evidence rows, with studies and ancestry as study entities. |
| `variant_index/vindex/build/hgvs.py` | 95 | True | Genomic HGVS (g. / m.) changes -> normalised CHROM:POS:REF:ALT keys, with the reference sequence as the arbite |
| `variant_index/vindex/build/mavedb.py` | 309 | True | MaveDB multiplexed assays of variant effect -> assay entities and variant-level evidence. |
| `variant_index/vindex/build/opentargets_gene.py` | 217 | True | Open Targets 26.09 gene-level perturbation and burden evidence. Every row here is attached to the gene; none o |
| `variant_index/vindex/build/opentargets_genetics.py` | 211 | True | Open Targets Genetics (platform 26.09): fine-mapped credible sets, studies and locus-to-gene predictions. |
| `variant_index/vindex/build/pathways.py` | 48 | True | Reactome pathway membership -> gene_set_member rows. |
| `variant_index/vindex/build/pgx.py` | 330 | True | Pharmacogenomic haplotype definitions and function: CPIC (API tables) and PharmVar (complete database 6.2.29). |
| `variant_index/vindex/build/regions.py` | 158 | True | Genomic regions that change how a call or a negative result is read: GIAB difficult regions and STRchive repea |
| `variant_index/vindex/build/runner.py` | 161 | True | Run adapters and finalise an index directory: source_release, entity_variant (deduplicated from evidence), man |
| `variant_index/vindex/callability.py` | 171 | True | Reference-confidence providers: what the input can say about positions where it holds no variant record. |
| `variant_index/vindex/crosscheck.py` | 485 | True | Independent cross-check of a whole-genome run against established command-line tools. |
| `variant_index/vindex/index_reader.py` | 239 | True | Read access to a built index (a directory of Parquet tables named as in schema.TABLES, or in-memory tables). |
| `variant_index/vindex/intervals.py` | 72 | True | Overlap queries over closed integer intervals, sized for millions of intervals per chromosome and millions of  |
| `variant_index/vindex/known_loci.py` | 308 | True | Known-locus checks: facts that are established in the literature and in public curation, tested against the bu |
| `variant_index/vindex/match_hap.py` | 139 | True | Haplotype (star-allele) calling from defining variants, with explicit evidence requirements. |
| `variant_index/vindex/match_repeat.py` | 107 | True | Repeat-expansion loci: a repeat length needs a repeat genotyper. |
| `variant_index/vindex/match_small.py` | 251 | True | Matching of called SNVs and indels (and MNV components) against the index, by normalised CHROM:POS:REF:ALT onl |
| `variant_index/vindex/match_sv.py` | 148 | True | Structural variant and copy-number matching. |
| `variant_index/vindex/normalise.py` | 100 | True | Allele normalisation: left-aligned, parsimonious, validated against the reference. |
| `variant_index/vindex/pgs.py` | 242 | True | Polygenic score models, kept separate from the variant index. |
| `variant_index/vindex/pgs_registry.py` | 174 | True | PGS Catalog model registry: the metadata that gates the interpretation of a polygenic score, kept apart from t |
| `variant_index/vindex/pgs_run.py` | 118 | True | Polygenic scores for one sample, one chromosome at a time, with the registry's gates beside every result. |
| `variant_index/vindex/provenance.py` | 320 | True | Provenance checks on a built index. |
| `variant_index/vindex/reads.py` | 354 | True | Targeted read-level evidence for small variants at the loci a callset cannot answer for. |
| `variant_index/vindex/refgenome.py` | 119 | True | Reference sequence access, contig naming and assembly validation. |
| `variant_index/vindex/reports.py` | 493 | True | Reports over a built index: coverage statistics, validation, the lifecycle of every acquisition target, a chan |
| `variant_index/vindex/resolvers.py` | 211 | True | Raw-record resolvers: the exact text an adapter hashed for an evidence row, recovered again from the acquired  |
| `variant_index/vindex/schema.py` | 222 | True | Table schemas and controlled vocabulary of the High-Value Human Variant Reference Index. |
| `variant_index/vindex/summary.py` | 52 | True | Per-axis evidence summaries. A summary lists what exists on each axis; it never ranks or adds axes together. |
| `variant_index/vindex/vcfio.py` | 436 | True | Sample VCF/gVCF reading for the matcher. |
| `variant_index/vindex/views.py` | 493 | True | Rules and views over the complete catalogue of matches. |
| `variant_index/vindex/wgs.py` | 447 | True | Whole-genome matching of one sample against the index, one chromosome at a time. |
