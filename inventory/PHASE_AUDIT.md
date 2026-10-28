# Phase gap audit

Generated 2026-10-28T02:41:12Z by `inventory/gap_audit.py` (regenerate rather than edit). 28 pass, 0 fail, 0 skipped, 1 blocked.

PASS means the check ran and held. SKIP means an input this checkout does not hold was needed, so the check proves nothing here. BLOCKED means only a person outside the session can settle it.

| Check | Result | Detail |
|---|---|---|
| requirements map is followable (paths and test modules exist, partial / missing statuses state their gap) | **PASS** | 44 requirements, 0 problems |
| REQUIREMENTS.md is the rendering of requirements.json (nobody edited it by hand) | **PASS** | identical |
| no module or command-line tool of variant_index is outside the requirements map | **PASS** | 55 files, all named |
| every test module is used by at least one requirement | **PASS** | 44 test modules, all cited |
| variant_index unit tests pass (including resource-warning errors) | **PASS** | 482 tests OK |
| acquisition helper tests pass | **PASS** | Ran 17 tests |
| every number the READMEs and the latest audit-log entry quote is recomputed and agrees | **PASS** | 35 pass, 0 skipped (data absent) |
| scripts and files that the documentation tells a reader to run or open exist | **PASS** | 25 documents |
| catalogue records validate and every acquisition target is claimed by a source record | **PASS** | 54 source records valid; 173 acquisition targets; 277 evidence probes |
| statements in the catalogue are traceable to stored evidence | **PASS** | 0 of 1254 text items contain tokens found in no stored evidence; 1 further items sit in fields the record lists as unverified |
| NHGRI, Ensembl, 1000 Genomes and MaveDB assets are unchanged apart from documented exceptions | **PASS** | baseline files: 97 / missing: 0 / changed: 0 / documented exceptions: 1 / added: 0 |
| no personal-genome file type is tracked, and none was ever added to any branch | **PASS** | 0 synthetic fixture file(s) in history, no other |
| nothing under data/ is tracked and results directories are git-ignored | **PASS** | ignored |
| repository policy: file sizes, personal-genome types, credential patterns | **PASS** | visibility: public (ziadwalieddine-create/dc) |
| the repository 'dc' is private | **BLOCKED** | visibility: public (ziadwalieddine-create/dc): only the owner can change repository visibility (the session proxy refuses settings writes) |
| the built index passed validation (vocabulary, ranges, keys, alleles, partitions, reference alleles) | **PASS** | 24,195,946 evidence rows, 50,067 REF alleles read back |
| provenance: release linkage, manifest agreement, recomputed record digests | **PASS** | 15,242 sampled records re-derived their digests across 32 releases; 0 release(s) without a resolver are listed, not counted |
| every release that holds evidence rows has a resolver, and the sampled records of all of them re-derived their digests | **PASS** | 32 releases, 15,242 records |
| lineage keys are shared between databases that describe the same study (the GWAS Catalog and Open Targets count a paper once) | **PASS** | 17,484 of 4,752,989 groups are written by more than one source; GWAS Catalog and Open Targets share 4,561 |
| two builds from the same inputs with different hash seeds are byte-identical, and they ran the code that built the index | **PASS** | 216 files, limit 30000 records per adapter, build code 72efa3a5592e4ec6 |
| the index was built by committed code: its build-code digest is that of the build modules at a commit of this history | **PASS** | commit HEAD (24 files); HEAD has the same build modules |
| ClinVar-derived rows of the index equal the ClinVar VCF as read by bcftools | **PASS** | 4,599,326 records, 4,360,956 classifications, 4,597,317 keys equal to bcftools norm |
| facts established outside the index (rsID positions, effect-allele orientation, classifications, star-allele definitions, curated assertions, retained unmapped rows) read back from it | **PASS** | 27 loci, 49 checks, all pass |
| a sample run agrees with bcftools on its calls and on the ClinVar evidence of those calls | **PASS** | 4,095,610 calls; 45,260 with ClinVar evidence |
| every unresolved high-value locus of the sample run states the evidence that would settle it | **PASS** | 44,718 unresolved loci, all with text |
| the sample run accounts for every record it read (small calls, structural, repeat, rejected, no-call and reference counts add up) | **PASS** | 4,048,342 records read; 4,095,610 small calls (0 read alleles carry no copy), 225 structural, 0 rejected |
| the sample run was made by committed code that has not changed since, against the index that is on disk now | **PASS** | commit 0e4dad0c3a, 18 modules unchanged, index manifests identical |
| the structural run on real calls (HG002 genotypes of the public HGSVC3 callset) accounts for every record and carries the summary of what it refused | **PASS** | 25,534 structural calls in 187 s; largest refusal list 9,129 characters |
| two matcher runs of one input with different hash seeds give identical result tables and manifests | **PASS** | 8 tables identical |
