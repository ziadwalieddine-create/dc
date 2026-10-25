# Audit log

Chronological record of additions, modifications, validations, corrections and unresolved issues for the reference-library work. Newest entries last. Entries state what was observed, not what was hoped.

## 2026-10-25: library extension (catalogue, acquisition, inventory)

### Preserved
- `scrape/` (NHGRI GWAS glossary, NHGRI 1000 Genomes page, Ensembl human genome registry, MaveDB human metadata and scripts) was not restructured, renamed or deleted. A SHA-256 baseline of its 97 tracked files was taken at commit `ddafb14` (`inventory/baseline_checksums.sha256`); `python3 inventory/audit_existing.py --check` verifies it.
- The MaveDB score tables (1,980 files) remain on branch `data/mavedb-scores`; `scrape/mavedb-human/verify_scores.py` still reports 0 problems.
- Unrelated content was not touched. The dating-engine repository (`DC-extended-`) is out of scope.
- Repository `ziadwalieddine-create/dc` is **public**. Nothing restricted was committed: downloaded data live in git-ignored `data/external/`; only manifests, scripts, reports and short evidence excerpts of public web pages are committed. Licence terms still apply to any later redistribution.

### Added (new directories only)
- `catalogue/`: 34 evidence-backed source records, controlled vocabulary, evidence collector (141 probes at the time of writing), generated inventory / coverage / gap / reuse documents.
- `acquisition/`: downloader, targets, manifest, validators, compatibility tests, and `mavedb_superseded/`.
- `inventory/`: audit of existing assets and the preservation guard.
- Root `.gitignore` (`data/`, `__pycache__/`).

### One modification of an existing file (necessity demonstrated)
- `scrape/mavedb-human/README.md`, Coverage-check paragraph only. It claimed every score set the public API can return was accounted for. That was false (see "Superseded score sets"). The edit is pinned by hash in `inventory/preservation_exceptions.json`; any further change to that file fails the guard.

### Findings that changed earlier conclusions
1. **Superseded score sets were missed.** The MaveDB search (used to build the human set) omits superseded versions. The CC0 bulk dump v4 (asOf 2026-06-24) lists 52 published score sets absent from the search listing; the API serves all of them (HTTP 200, flagged superseded). 51 have a human or accession target. They are kept as an additive, labelled extension in `acquisition/mavedb_superseded/` (79 tables, all row counts equal to `numVariants`). The statistics endpoints (2,995 sets; 1,617 human targets) are consistent with 2,938 + 52 hidden + about 5 non-CC0 hidden sets, and 1,477 + 84 + 51 + about 5 targets. Consistent, not proven.
2. **The MaveDB bulk dump is much richer than assumed.** Besides scores and counts it provides per-variant annotations (MaveDB-mapped genomic/coding/protein HGVS, VRS digest, VEP consequence, gnomAD v4.1 allele frequency, ClinGen Allele Registry ID, historical ClinVar significance), mapped-variant JSON and GA4GH VA-Spec annotations. Of those, `main.json` and the 1,170 annotation CSVs (708 MB) were extracted; the 9.5 GB `mapped/` and 3.4 GB `va/` files were left in the archive.
3. **Two mapping products are complementary, not interchangeable.** Ensembl's VEP MaveDB file expands protein-level variants into every codon-compatible genomic allele; MaveDB's mapper gives genomic HGVS only to nucleotide-level assay variants. A one-to-one concordance test was the wrong test and was replaced by a recorded finding.
4. **Existing tables equal the dump.** On a sample of 31 of our sets the dump's scores are identical to the API-built tables (same rows, accessions and values).
5. **Ensembl release 116 already ships VEP plugin data** (MaveDB, LOEUF, SpliceAI, SpliceVault, ProtVar, mutfunc, Enformer and others). They are catalogued; only the MaveDB file was acquired.
6. **Open Targets 26.09 already integrates** credible sets, colocalisation, locus-to-gene, ClinGen, ClinVar (EVA), IMPC, CRISPR screens, DepMap essentiality, gene burden, pharmacogenomics and interactions (61 datasets from the FTP listing). Adopted as the primary route for those dimensions; none downloaded yet.

### Mistakes made and corrected during the work
- Guessed a wrong gnomAD file size (rounded); the downloader rejected the file and it was re-fetched with the exact size and the MD5 published on the gnomAD downloads page.
- BioGRID ORCS "download" returned an HTML landing page with HTTP 200; the downloader now checks file signatures. The correct URL was found in the page HTML. The human archive is 752.7 MB: a size on the page was first attributed to the wrong file.
- Chose `gwas-catalog-studies.tsv` (no study-accession column); removed it and its manifest entry. The associations file carries `STUDY ACCESSION`.
- An expectation of 54 tissue columns in GTEx v11 median TPM was wrong (68 columns incl. sub-sites); the check was corrected after inspecting the header.
- A ClinGen validator ignored the `#`-prefixed header line; fixed.
- A PharmVar version string had a wrong year in a draft record; fixed before generation.

### Access blocks recorded (not circumvented)
DepMap (human-verification page), DECIPHER (`Access Denied`), bio.tools API (JavaScript bot challenge; OpenEBench used instead), AlphaFold DB download page (browser rendering denied; API and FTP work), dbNSFP (institutional registration), PharmVar terms (page not found), eQTL Catalogue REST API (HTTP 410, retired), ENCODE data-use page (HTTP 403), Ensembl VEP plugin documentation (HTTP 504 at probe time).

### Unresolved
- Licences not verified live: ENCODE/SCREEN, eQTL Catalogue, InterPro, DGV, dbVar, BioGRID ORCS, PharmVar, DepMap, CPIC-specific terms, NHGRI reuse terms, GENCODE.
- Not acquired (large or gated): dbNSFP, gnomAD VCFs and SV, AlphaFold proteome, SpliceAI/SpliceVault/Enformer/ProtVar/mutfunc, Open Targets datasets, DepMap, IMPC, DGV (tier B), AlphaMissense variant file (tier B).
- The ClinVar join counts exact SNV coordinate matches for Ensembl-expanded candidate alleles; it measures coordinate compatibility, not clinical meaning.
- No whole-genome analysis engine has been built, by design.

## 2026-10-25 (second phase): tier-B files, Open Targets, per-protein and regional retrieval, licence checks, corrections

This entry supersedes the "Unresolved" list of the first entry above wherever they differ.

### Added (data are in git-ignored `data/external/`; only scripts, manifests, reports, short evidence excerpts and tables of hashes are committed)
- **AlphaMissense_hg38** (642,961,469 B; MD5 equals both the GCS object metadata and the Zenodo record): 71,697,556 rows, 19,233 transcripts, 19,117 UniProt accessions. Classes follow the documented thresholds (0.34 / 0.564) except at the two boundaries, where 3,838 rows differ because scores are printed to four decimals and classes use the unrounded value (checked: every such row has printed score 0.34 or 0.564).
- **DGV GRCh38 2025-12-01** (364,191,405 B; equals the size listed on the site; the site publishes no checksum): 1,987,094 rows.
- **gnomAD v4.1 genome SV sites VCF** (1,743,544,245 B; MD5 equals the downloads page): 2,154,486 records. The 721 `SVTYPE=CNV` records have no AF because they are multi-allelic copy-number sites (CN_COUNT / CN_FREQ instead); the first validator wrongly required an AF on every record.
- **Open Targets 26.09**: 27 Parquet datasets, 8.55 GiB, every part checked against Content-Length and PAR1 framing, our SHA-256 per part in `acquisition/manifest_files/`. All 61 datasets are accounted for: 27 acquired, 34 deferred with reasons (colocalisation 15.4 GiB, baseline_expression, association scores, literature mining, re-expressions of sources held directly, cancer and drug content, transcript).
- **UniProt humsavar** (release 2026_03 of 02-Sep-2026): 85,015 rows for 83,859 distinct variants; row and per-category counts equal the statistics printed in the file header.
- **AlphaFold DB and InterPro per protein** for the 549 UniProt accessions of the MaveDB set (`fetch_per_protein.py`): 536 of 549 accessions have AlphaFold DB prediction metadata (13 have none); 529 canonical models with per-residue confidence are stored (1,594 files, 356.7 MB; mean pLDDT 70.25); model length equals UniProt length for 520 of 520 reviewed entries checked; InterPro release 110.0 gives matches for 547 accessions (2 have none): 11,977 match records, every fragment inside its protein, protein length equal to UniProt's on all 11,775 records of reviewed entries. Only canonical `AF-<acc>-F1` models are stored; 7 very long proteins have only isoform models in AlphaFold DB and are recorded as `isoform_only`.
- **gnomAD v4.1 exome and genome records at the positions of MaveDB-mapped SNVs** (`fetch_gnomad_regional.py`, remote tabix over HTTP range requests, no whole-file download): three callsets, 76,495 exome, 22,162 genome and 82,528 joint (exomes and genomes pooled) single-nucleotide records, from 1,444 regions on chr1-22/X/Y (102 of the 1,546 regions lie on alternative contigs or MT and are not in gnomAD's VCFs); 0 read retries; in the joint slice AC and AN equal the sum of the exome and genome counts on all 82,528 rows. Each region is read until two consecutive reads agree.
- `extract_dump.py` (byte-for-byte reproduction of the earlier MaveDB dump extraction: 1,297 files, 0 differences), `requirements.txt` (pyarrow, pysam), `test_helpers.py` (11 unit tests, including the interval index against brute force), `inventory/check_repo_policy.py`.

### Licences verified from live pages in this phase (evidence in `catalogue/evidence/`)
BioGRID ORCS (MIT), GENCODE open access with the EMBL-EBI terms and MANE under the NCBI data policy, DGV (TCAG disclaimer, acknowledgement), InterPro (CC0 1.0), ENCODE (unrestricted use after release; citation page itself returns HTTP 403), eQTL Catalogue (CC BY 4.0), NHGRI (public domain, acknowledgement requested), dbVar (NCBI policy), UniProt humsavar (CC BY 4.0 per the file footer). Earlier records for these said "unverified"; the later evidence replaced them.

### Corrections of earlier results
1. **The coordinate set used by the compatibility tests was too small.** `ens_snv` was a dict keyed by accession, so only one of the codon-compatible genomic alleles Ensembl lists for each protein-level variant survived: 341,851 unique SNV alleles instead of 399,389, and 136,894 positions instead of 145,710. Joins that used it (ClinVar exact match, SCREEN, CNV, ClinGen overlaps, AlphaMissense) are replaced by the numbers in `acquisition/compat_report.json`. Found because the regional gnomAD script, which used a set, reported different counts.
2. The MANE transcript map kept only the last row per gene and dropped some MANE Select transcripts; it now uses every row.
3. The compatibility script depended on a file in a session scratch directory for its dump-versus-API sample; the sample is now derived from the extracted files.
4. A list comprehension rebuilt a set inside its loop (quadratic) and stalled a run for 18 minutes; found because the process held no data file open.
5. Validator expectations written before seeing the data were wrong in five places and were corrected after inspecting the data, not relaxed to pass: ontology identifiers (`OBA_VT…`, and 1,490 of 35,477 biosample IDs are IRIs, so that check requires 95%), the `_SUCCESS` marker (single-file exports and `l2g_prediction` have none; contiguous part numbers are required instead), AF on every SV record, AlphaMissense class boundaries, and humsavar footer and one-row-per-disease structure. AlphaFold: an isoform model had been stored as if canonical for 7 long proteins; the downloader now stores canonical models only and the validator checks it.
6. htslib printed `[E::bgzf_read_block] Failed to read BGZF header` once during the first regional run. That run was stopped, and the script now reads every region until two reads agree and raises after five attempts.
7. The first Open Targets L2G test rebuilt a 4.8M-entry hash for each of 200 parts; replaced by one set.

### Cross-source findings (all in `acquisition/compat_report.json`)
- **Open Targets quotes the gnomAD v4.1 joint (exomes + genomes pooled) frequencies.** Across nine ancestry groups, 120,892 of 120,892 comparisons equal the joint slice within 1e-3 relative, against 89,759 of 114,428 for the exome slice and 34,186 of 50,531 for the genome slice. The MaveDB dump's `gnomad.gnomad_af` is the same quantity: 4,625 of 4,625 variants equal the joint AF, whereas against the genome-only AF the median ratio is 0.95 (10th percentile 0.09) and against the exome-only AF 1.13. The joint slice is therefore the frequency of record; exome or genome AFs alone differ.
- Open Targets target constraint is the gnomAD v4.1.2 MANE-select constraint: LoF o/e is equal for 17,421 of 17,421 genes.
- **GA4GH VRS identifiers do not interoperate here.** None of the 23,667 comparable variants (17,202 exome, 6,465 genome) has MaveDB's post-mapped VRS digest among gnomAD's `VRS_Allele_IDs`. Coordinates are the reliable join key until both sides are computed with one library version.
- AlphaMissense (canonical UniProt numbering) agrees with the protein change implied by Ensembl's file (RefSeq numbering) for 237,668 of 252,714 matched alleles (94.0%); the sampled disagreements are constant offsets from isoform numbering (for example D206N against D186N).
- UniProt humsavar joined through AlphaMissense reaches 3,390 MaveDB alleles (1,492 LP/P, 331 LB/B, 1,567 US). Where ClinVar also gives a clear call (1,023 alleles) humsavar agrees in 1,006 (98.3%).
- Coordinate joins on the corrected set of 399,389 SNV alleles (145,710 positions): ClinVar 38,597 exact matches (9.7%), AlphaMissense 293,696 (73.5%), gnomAD joint 71,993 (18.0%), exomes 66,524 (16.7%), genomes 18,954 (4.7%), Open Targets variant index 32,335 (8.1%). 63,025 positions lie in a SCREEN cCRE, 137,307 inside a gnomAD genome SV interval, 128,077 inside a DGV interval and 123,127 inside a gnomAD exome CNV interval; an interval overlap is a coordinate relationship, not evidence about the variant.
- Open Targets integrity: every credible-set study ID (2,840,866) resolves in the study table, all 3,293,342 L2G loci resolve in the credible sets, all 27,396,304 interaction-evidence rows resolve in the interaction table, and 7,083,286 of 7,083,286 evidence disease IDs resolve in the disease ontology. 7,312 credible sets contain a MaveDB-mapped SNV (3,212 as lead variant); 484 MaveDB genes have an L2G score of 0.5 or more for some locus.
- The Ensembl registry (existing asset) names GRCh38.p14 as its reference assembly, with variation and regulation data and GENCODE annotation, and it resolves the only assembly MaveDB names (GRCh38, 70 sets). This is a name-level check; sequence identity needs the FASTA, which is not acquired.

### Results of the final runs
- `validate.py`: 80 of 80 targets pass (80 manifest records, 16.9 GB on disk in `data/external/`).
- `compat_tests.py`: 64 tests, 59 PASS, 0 FAIL, 5 INFO (relationship findings that carry no pass criterion).
- `python3 -m unittest acquisition.test_helpers`: 11 tests pass.
- `scrape/mavedb-human/verify_scores.py`: 0 problems. `inventory/audit_existing.py --check`: 97 baseline files, 0 missing, 0 changed, 1 documented exception. `inventory/check_repo_policy.py`: PASS (the repository is public).

### Access blocks (not circumvented)
Unchanged from the first entry, plus: ENCODE citation page (HTTP 403 again, rendered), UCSC ENCODE terms (certificate not trusted by the probe's trust store, not bypassed), Ensembl VEP plugin documentation (HTTP 504, then a timeout).

### Unresolved
- **`ziadwalieddine-create/dc` is still public.** The proxy refuses repository-settings writes. Everything tracked is world-readable, including history that contains the earlier dating-project case files; no history rewrite was done (user decision). Fix: GitHub → Settings → Danger Zone → Change visibility, or `gh repo edit ziadwalieddine-create/dc --visibility private --accept-visibility-change-consequences`. The `data/mavedb-scores` branch carries CC BY / CC BY-SA / CC BY-NC-SA tables with a notice; making the repository private does not change those licences.
- `urn:mavedb:00001251-a-1` (CC0, published 2018-06-29, modified 2026-06-05) is absent from the dump for no reason found.
- Licences still unverified: PharmVar (terms page not found), DepMap (verification wall), DECIPHER (access denied), FAIRsharing (API terms need an account). Other facts still unverified: CPIC-specific terms, the GWAS Catalog assembly field, the Ensembl VEP plugin list and per-plugin licences, dbVar currency and contents.
- Gated or blocked: dbNSFP (registration), DepMap, DECIPHER, the FAIRsharing API, the bio.tools API (bot challenge; OpenEBench is used instead), PharmVar terms.
- Large and deliberately not fetched: gnomAD whole-chromosome VCFs, Open Targets colocalisation, genome-wide SpliceAI / SpliceVault / Enformer / ProtVar / mutfunc plugin data, 1000 Genomes phased panels, the genome FASTA and GENCODE annotation.
- Slice-based frequencies cover only the positions of MaveDB-mapped SNVs; a missing allele is not a measured frequency of zero.
- No whole-genome analysis engine has been built, by design.
