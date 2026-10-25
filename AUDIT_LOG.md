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
