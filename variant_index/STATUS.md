# Status of the High-Value Human Variant Reference Index

Written from test, build and run results, not from intent: every figure comes from an artefact named beside it, and `inventory/repository_audit.py` recomputes the ones it can. `BUILD_REPORT.md` holds the full tables, `REQUIREMENTS.md` maps every clause of the project goal to code, tests and a status, `inventory/PHASE_AUDIT.md` is the adversarial audit of both.

Requirements: 44, of which 30 verified, 12 partial, 1 unavailable input and 1 blocked external (each partial one states its gap in `REQUIREMENTS.md`).

## What was built and run

- **Code.** The index was built on 2026-10-28 by the build modules with digest `72efa3a5592e4ec6` (24 files, found in the history by `inventory/gap_audit.py`); the sample runs were made by commit `0e4dad0` with a clean tree, the last commit on this branch (`ccr-b13fd75a-ra9gfy`) that changed `variant_index` code; later commits change documents and audit reports only. 482 unit tests pass (`python3 -m unittest discover` in `variant_index/tests`, resource warnings are errors; the cross-checks against bcftools and plink2 run when `VINDEX_TOOLS_BIN` names them).
- **Index build** (`build_index.py`, all 22 adapters, 103 min on 4 CPUs, sampled peak 4.3 GB): 24,195,946 evidence rows from 72 source releases, 6,857,645 variant entities, 4,408,192 structural events, 7,516,808 regions, 78,797 genes, 1,410,446 studies, 3.1 GB of Parquet under `data/index` (git-ignored, rebuilt from the acquired sources by one command).
- **Validation** (`make_reports.py`): ok: 97 checks over 24,195,946 evidence rows and 50,067 reference alleles read back from the FASTA, 0 violations.
- **Provenance** (`verify_index.py`): ok: 72 releases linked and agreeing with the acquisition manifest, 32 releases with a resolver, 15,242 records re-opened in the acquired files, 15,242 digests re-derived, 0 mismatched, 0 unresolved, every release has a resolver.
- **Idempotence of builds** (`check_idempotence.py`, two builds in separate processes with different hash seeds, `--limit 30000` records per adapter): identical: 216 files re-hashed, adapter code equal True, rows, releases and problems equal, seeds 1 and 2, builds of 23.7, 18.1 min.
- **ClinVar against bcftools** (`crosscheck_bcftools.py --clinvar-fields`): ok: {"claim_value_differs": 0, "condition_ids_differ": 0, "condition_text_differs": 0, "extra_in_index": 0, "extra_json_differs": 0, "frequency_differs": 0, "gene_count_differs": 0, "gene_id_inconsistent": 0, "index_rows": 6413262, "key_differs": 0, "mapped": 4597317, "mapped_but_bcftools_rejects": 0, "missing_in_index": 0, "n_submitters_differs": 0, "origin_differs": 0, "review_stars_differ": 0, "review_status_differs": 0, "unmapped_but_bcftools_accepts": 0, "unmapped_consistent": 2009, "vcf_classifications": 4360956, "vcf_records": 4599326, "vcf_records_without_classification": 241297}
- **Known loci** (`verify_known_loci.py`, facts established outside the index): all pass: 27 loci, {"ABSENT": 0, "FAIL": 0, "PASS": 49}.
- **Dataset lifecycle** (listed, acquired, verified, mapped, integrated; `BUILD_REPORT.md`): 173 targets: 173 listed, 173 acquired, 173 verified, 32 mapped, 72 integrated (mapped as the index's lifecycle report counts it, by an evidence row mapped to a normalised entity; `inventory/repository_audit.py` also counts a passing join test and reports 122, requirement A3 states both).

## Whole-genome run on a public truth set

Input: the GIAB HG002 GRCh38 benchmark callset v4.2.1 (public, not anyone's personal genome; `fetch_public_testset.py` downloads it and checks the digests in `public_test_inputs.json`). Command: `match_wgs.py --callable-bed <benchmark high-confidence BED>`; results stay under `data/results/HG002` (git-ignored).

Run 1,374 s (commit `0e4dad0`, clean tree, code unchanged during the run). 4,095,610 small-variant alleles (3,463,511 SNVs, 323,073 deletions, 308,981 insertions, 45 complex; 47,269 multiallelic records split); 360,956 `EXACT_MATCH` (ClinVar 45,260, GWAS Catalog 125,345, Open Targets fine-mapping 272,321, ClinPGx 2,989, ClinGen 102, gnomAD 66, MaveDB 57; a call can have evidence from several sources), 17 `FILTERED_CALL`, 3,734,637 `NO_INDEX_RECORD` (absence from the index, not evidence of no effect); 1,175,013 calls lie in GIAB difficult regions (flagged, not removed).

2,615,243 evidence links: 2,445,088 to called alleles (each states the effect allele and the copies of it the sample carries) and 170,155 at 37,640 verified-reference sites, where the evidence is about the reference allele (GWAS Catalog rows; two copies of the effect allele on an autosome). Of the 1,111,374 mapped GWAS Catalog rows that state an effect allele, 400,814 (36.1%) name the reference allele. 44,718 high-value indexed alleles the file cannot answer for are listed with the evidence that would settle each (44,352 `ABSENT_FROM_CALLSET` outside the callable region, 365 `UNRESOLVED`, 1 `FILTERED_CALL`). Star alleles, 27 genes: 13 `RESOLVED`, 3 `NO_CONSISTENT_DIPLOTYPE`, 8 `REQUIRES_SPECIALISED_CALLER`, 3 `UNRESOLVED`; repeats: all 82 loci `UNRESOLVED` (no repeat genotyper output was supplied); 225 structural calls were made from long indels (184 have an equivalent event).

Checks of this run: bcftools recomputes the calls (4,095,610 on both sides, 0 differences) and the ClinVar evidence of the calls (45,260 keys on both sides, 0 differences); the verified-reference sites are recomputed from the index, the callable-region BED and bcftools (86,877 candidate sites, 82,003 inside the region, 44,363 with a record of the sample, 37,640 expected, 37,640 linked, 0 on either side only); two runs with hash seeds 1 and 777 on chromosomes 21 and 22 give byte-identical tables (8 tables) and equal manifests.

Views over the 360,956 exact matches (the sample's `views/` directory and its report, `data/results/<sample>/report.md`): association 111,804; clinical effect 402; fine-mapped 272,320; functional effect 117; pharmacogenomic 1; conflicting evidence 1,133; several axes assert an effect 94,848. 10,306 traits named by effect-asserting evidence roll up into 25 therapeutic areas (363 identifiers are not in the ontology table); 2,266 pathways.

## Structural and copy-number calls, polygenic scores, alternate contigs

- **Structural calls**: the HGSVC3 genotypes of HG002 as a VCF of 25,534 symbolic calls (9,659 DEL, 15,801 INS, 74 INV; 187 s). Every call has an equivalent event in the index. HGSVC3 is itself in the index, so this is a scale and accounting test and not a validation of the equivalence rules, which are tested on adversarial fixtures (unequal overlapping events, breakpoints, dosage).
- **Polygenic scores**: PGS000001, PGS000018 and PGS002724 on the HG002 callset (262 s) are scored with their gates stated: all three `ancestry_unassessed` and `uncalibrated`, PGS000001 also `weight_coverage_below_0.9`. No score was interpreted.
- **Alternate contigs**: a synthetic VCF of 4 records on alternate loci and an unplaced scaffold gives 3 small calls (2 exact matches, 1 without an index record), 1 rejected record and 1 unplaced record (38 s).

## Modules

<!-- modules:start (generated by make_status.py) -->

| Module | Status | Tests | Note |
|---|---|---|---|
| allele normalisation, identity keys, multiallelic splitting | **implemented** | test_normalise, test_vcfio, test_match_small |  |
| small-variant matching with verified-reference and absence semantics | **implemented** | test_match_small, test_build_clinvar, test_wgs |  |
| structural / copy-number matching (type, breakpoints, dosage; no overlap-only equivalence) | **implemented** | test_match_sv, test_wgs |  |
| star-allele / haplotype calling from defining variants (CPIC, PharmVar) | **implemented** | test_match_hap, test_build_pgx, test_wgs | structural genes (CYP2D6, HLA, KIR ...) require a dedicated caller; its result can be supplied with --external-haplotypes |
| repeat-expansion classification against STRchive thresholds | **partial** | test_match_repeat, test_build_regions, test_wgs | classifies the output of ExpansionHunter / TRGT and keeps every repeat genotype of the file (NO_INDEX_RECORD where the index has no locus); sizing a repeat from short reads is not attempted |
| contigs outside the primary assembly (alternate loci, unplaced scaffolds, decoys) | **implemented** | test_other_contigs | indexed under the shared NA partition and matched in one extra pass of a whole-file run; records on contigs the reference does not know are listed as rejected |
| evidence about the reference allele linked at verified-reference genotypes (gVCF block or callable region; two copies on an autosome) | **implemented** | test_reference_effect, test_effect_alleles | sites the sample cannot answer for, or at which another allele is called, have no link and are counted in the run manifest |
| independent cross-checks with bcftools (carried alleles, ClinVar evidence, ClinVar fields of the index, verified-reference sites) | **implemented** | test_crosscheck, test_reference_effect |  |
| idempotence (two builds, separate processes, different hash seeds) and hash-seed independence of matching | **implemented** | test_idempotence, test_wgs |  |
| polygenic scores with weight-coverage, ancestry and calibration gates | **implemented** | test_pgs, test_pgs_run, test_pgs_registry | ancestry inference and reference score distributions are not part of the repository |
| ClinVar aggregate, submission, structural and composite-record evidence | **implemented** | test_build_clinvar, test_build_clinvar_vcv |  |
| ClinGen, CPIC, PharmVar, ClinPGx, MaveDB, GWAS Catalog, Open Targets, gnomAD, DGV, HGSVC3, GIAB, STRchive, Reactome adapters | **implemented** | test_build_* |  |
| provenance verification | **implemented** | test_provenance, test_resolvers | release linkage and manifest agreement are checked for every release; the SHA-256 of the raw record is recomputed from the acquired file for a sample of 500 rows of every release that holds evidence (one rule per source in vindex/resolvers.py); entity rows carry the release of their source but no digest of their own |
| read-level evidence (BAM / CRAM pileup, targeted re-genotyping) for no-call, filtered and absent loci | **unimplemented** |  | no reads exist in this environment and no read-level caller is part of the repository; reference status comes from gVCF blocks or a declared callable-region mask, and the unresolved loci name the evidence that would settle them |
| transcript HGVS (c., n.) to genomic coordinate mapping | **unimplemented** |  | records that carry only transcript HGVS stay unmapped with the reason |
| consequence annotation of called variants (VEP / bcftools csq) | **unimplemented** |  | the index stores source-stated consequences; none is computed for a sample |
| GRCh37 -> GRCh38 lift-over of legacy ClinVar CNV records | **unimplemented** |  | those records keep their GRCh37 coordinates and are unmapped with the reason |
| ancestry inference from genotypes (reference panel such as 1000 Genomes) | **unimplemented** |  | the 1000 Genomes metadata is acquired, the phased panels are not (size); polygenic-score ancestry gates use a declared group |
| annotation resources acquired for later use (the Open Targets variant table, AlphaMissense, ProteinGym, AlphaFold DB, InterPro, SIFTS, UniProt, GTEx, HPA, ENCODE, BioGRID, 1000 Genomes metadata) | **unimplemented** |  | acquired, validated and join-tested (acquisition/compat_report.json) but no adapter writes them into the index; the lifecycle report lists them as acquired and not integrated (the Open Targets variant table is read, but only to settle the orientation of GWAS Catalog effect alleles) |
| ontology roll-up of traits (EFO / MONDO ancestors) | **partial** | test_wgs | the trait ids of a sample's matches are grouped by ancestor and therapeutic area through the Open Targets ontology table (96% of the ids of the HG002 run); MedGen, OMIM and MeSH ids are not in it and stay as the sources wrote them; the index itself keeps the sources' ids |
| GWAS Catalog full summary statistics | **not_acquired** |  | the association table and study metadata are integrated; 181,588 of 232,202 catalogue studies list full summary statistics (per-study files, mostly molecular-QTL sets), none is downloaded and no adapter for them exists; choosing the relevant studies is not done |
| Open Targets colocalisation | **not_acquired** |  | credible sets, fine-mapping and L2G are integrated; the colocalisation datasets (15.4 GiB) are not downloaded and have no adapter |
| personal 30x GRCh38 SNP/indel, CNV and SV calls | **unavailable_input** |  | no personal call file exists in this environment; matching is verified on synthetic fixtures and, when run, on a public truth set |

<!-- modules:end -->

## What this environment could not do, and what is not implemented

**Input that does not exist here**
- A personal 30x GRCh38 callset (SNPs, indels, CNVs, SVs, repeat genotypes). The matcher is verified on synthetic fixtures with adversarial cases and on a public truth set (GIAB HG002), including the HGSVC3 structural genotypes of the same sample.
- Reads (BAM / CRAM). Without a gVCF or a callable-region mask a locus that the callset does not list stays `ABSENT_FROM_CALLSET`; no read-level provider exists to settle it, and the unresolved-locus table names the evidence that would.

**Acquired in part or not at all (size, the free disk of this environment is a few GB)**
- GWAS Catalog full summary statistics: 181,588 of the 232,202 catalogue studies list them (per-study files, mostly molecular-QTL sets); none is downloaded and the relevant studies have not been chosen. Association-level evidence does not stand in for them.
- Open Targets colocalisation (15.4 GiB): credible sets, fine-mapping and locus-to-gene scores are integrated, colocalisation is not.
- gnomAD whole-chromosome VCFs: allele frequencies are indexed only at the positions covered by the acquired regional slices; gene constraint, loss-of-function curation, structural and copy-number events are integrated.
- PGS Catalog: 3,467 of 7,029 harmonised GRCh38 scoring files are acquired (a tiered subset); the others are registered as `listed_not_acquired` and cannot be scored here.

**Acquired but not written into the index** (the lifecycle report lists them as acquired, verified and mapped, not integrated): the Open Targets variant table (read only to settle the orientation of GWAS effect alleles), AlphaMissense, ProteinGym, AlphaFold DB, InterPro, SIFTS, UniProt, GTEx, HPA, ENCODE, BioGRID, 1000 Genomes metadata and phased panels.

**Not implemented**
- Transcript HGVS (c., n.) to genomic mapping: records that carry only transcript HGVS stay unmapped with the reason.
- Consequence annotation of called variants (the index stores the consequences its sources state).
- GRCh37 to GRCh38 lift-over: 112,628 ClinVar records that have only GRCh37 coordinates keep them and are unmapped with the reason.
- Ancestry inference from genotypes and reference score distributions: polygenic-score ancestry gates use a declared group, and no calibrated percentile is produced.
- Structural star-allele genes (CYP2D6, HLA, KIR) need a dedicated caller (Cyrius, Aldy, StellarPGx, PharmCAT); none is installed and its result can be supplied with `--external-haplotypes`. Repeat lengths are not estimated from reads: ExpansionHunter / TRGT output is classified against STRchive thresholds.
- Colocalisation, summary-statistic and annotation-resource adapters named above.

**Evidence the sources cannot give, kept as such** (not defects): rows kept with their reason as unmapped, ambiguous or rejected, none dropped: GWAS Catalog 776,797 (no alternate allele known 349,453; author position assembly unverified 222,114; strand unresolved palindromic 115,520; multi base risk allele 30,396; ...); ClinVar 148,120 (only grch37 location no liftover performed 112,628; no vcf key or structural event for variation 15,720; no grch38 event location of at least 50bp 7,855; no vcf key for variation 6,162; ...); ClinGen 5,803 (no clinvar or allele registry match 5,732; variant named by hgvs only 55; region location not parsed 16); ClinPGx 3,884 (hla allele not defined in index needs hla typing 1,732; location not parsed 435; named polymorphism without definition or coordinates 403; annotation spans several alternate alleles 284; ...); MaveDB 885 (reference base differs:C 430; reference base differs:G 258; reference base differs:A 197); Open Targets 322 (symbolic 316; ref mismatch 4; variant id not parsable 2); gnomAD constraint 38 (gene not in gene map 38). A record that cannot be placed on a normalised key stays in the index with its source accession (and its position where the source gives one); a rejected strand is never guessed.

**Limits of the verification**
- Digests of raw source records are recomputed for a deterministic sample of 500 rows per release (one rule per source), not for every row; entity rows (genes, studies, regions, events, haplotypes) carry the release of their source but no digest of their own.
- The idempotence of builds is checked at `--limit 30000` records per adapter, and the idempotence of the matcher on two chromosomes with two hash seeds; a second full-size build is not made.
- Shared studies: the GWAS Catalog, Open Targets, ClinPGx and MaveDB share `pmid:` lineage groups, and ClinGen expert-panel classifications share the lineage of the ClinVar variation they name. The CPIC guidelines and ClinPGx's annotations of them have no common identifier and still count as independent.
- Evidence about the reference allele: a call is linked to the rows of its alternate allele, and the link states the copies of the stated effect allele that the genotype carries (`effect_allele_dosage`). A site whose genotype is a verified reference (a gVCF block, or a callable region with its stated assumption) is linked to the rows whose effect allele is the reference allele (`call_kind` `verified_reference`; two copies on an autosome, none stated on X, Y and the mitochondrion, where the providers do not record the ploidy). A site the sample cannot answer for, or at which another allele is called, has no link; each such site is counted in the run manifest by what the sample says about it (`reference_allele_effect_sites`), and a high-value one is listed as unresolved. The views are about the site: the copies are in `evidence_links.parquet`.
- Gene perturbation evidence is gene-level (IMPC knockouts, the CRISPR and essentiality evidence that Open Targets integrates, ClinGen dosage, gnomAD constraint); the BioGRID ORCS screens are acquired, not integrated, and DepMap's own files are not acquired (access restricted). No evidence about a specific allele comes from a gene perturbation.
- The PGS trial scores three models on a public callset to exercise the gates; no score was interpreted.
- The star-allele logic is verified against the published CPIC and PharmVar definitions by tests and known-locus checks, not against another caller: a cross-check against PharmCAT was prepared and not run (the environment's permission policy did not allow the downloaded tool to be executed).
- The repository is public (see below), so the first requirement of the standing conditions cannot be met from inside this session.

## Known residuals of the data

- **Verbatim source text.** Trait, cohort, curator and criteria strings are stored as the sources wrote them: double spaces in 55,137 MaveDB, 4,408 Open Targets, 3,212 GWAS Catalog and 52 ClinGen condition texts (and in 5,186 ClinVar submitter names), no-break spaces in 253 Open Targets and 32 GWAS Catalog condition texts, 260 Open Targets cohort lists with the irregular spacing of the source (`EB ,  FinnGen`), and in the ClinVar files themselves 106 cells with mojibake and 64 with HTML entities (`&amp;`). The raw record is recoverable by the locator and its digest; a join on text should normalise whitespace first.
- **Units.** 1,532 GWAS Catalog effect rows from one curator annotation carry the unit `unit novel` (the trailing word `novel` follows the direction word in the Catalog's interval text and is kept).
- **p-values of zero.** 12,367 GWAS Catalog rows (association and effect rows of the same association) and 18 ClinPGx rows have `p_value = 0`: the source reported zero (the Catalog's value underflows a double, its `-log10 p` is kept in `extra_json.pvalue_mlog`; ClinPGx printed `0.000`). The flag `p_reported_as_zero` is in `extra_json`; the assertion counts them as below the genome-wide threshold.
- **Fine-mapping method in `review_status`.** For Open Targets credible sets the method and its LD basis are stored in `review_status` (there is no review).
- **Identifiers that submitters typed.** ClinVar condition identifiers are kept as written where merging them would guess the submitter's intent: 4,538 `Gene:` ids (ClinVar's own condition-database field), 180 `MESH:` (tree numbers and descriptor ids; the VCF writes descriptors `MeSH:`, 758,268 of them), 60 `HPO:` (`HPO:HP:0000609`, `HPO:0000407`; the 555,426 `HP:` ids are merged), 61 `OMIM phenotypic series:PS...` and a few names after a database prefix (`HP:gait ataxia`). Spellings with one possible reading are merged (`MONDO_0005393`, `MONDO:MONDO:0010178`, `Human_Phenotype_Ontology:HP:0001250`, `OMIM: 612621`, the bare MONDO numbers of the ClinGen actionability export).
- **A name among identifiers.** One Open Targets FinnGen study lists `gastric adenocarcinoma` among its EFO and MONDO ids: the 7 evidence rows of that study keep it in `extra_json.condition_ids_not_identifiers` and the study entity keeps identifiers only.
- **Ancestry without discovery samples.** 8,761 Open Targets GWAS credible-set rows (0.4%) carry the LD population code `nfe` as their ancestry, because the study reports no discovery ancestry (`extra_json.ancestry_basis` says `ld_population_structure`); the other 2,362,873 carry the Catalog's names.
- **Gene-burden rows without a project or cohort.** 225 Open Targets burden rows have `unspecified` in their lineage key for the field the source left empty; one row with an odds ratio of 0 keeps the reported value in `extra_json` and stores no ratio.

## Privacy and licences

- The repository `ziadwalieddine-create/dc` is still **public** (`inventory/check_repo_policy.py`, live check). Only its owner can change that: the session's proxy refuses repository-settings writes. No personal genome file and no sensitive third-party data is tracked or was ever added to any branch (`inventory/gap_audit.py` scans the tree and all history). Private status would not waive any source licence; every release carries its licence (`source_release.licence`) and the catalogue lists the terms.
- No personal 30x callset exists in this environment. The sample pipeline is verified on synthetic fixtures and on the public HG002 callset; a personal callset is the one input still needed (`unavailable_input` in `REQUIREMENTS.md`). Nothing in a personal run is sent anywhere: matching is a local bulk join.
