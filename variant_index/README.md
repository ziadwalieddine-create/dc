# High-Value Human Variant Reference Index (`vindex`)

First product of the evidence library: a reusable index of human variant, haplotype, gene, region, assay and study evidence built from
the sources acquired under `acquisition/`, and a matcher that interprets a person's called variants against it. Status of every module
is stated in `STATUS.md` (written from test and build results, not from intent).

```
source files (data/external, pinned by acquisition/manifest.jsonl)
   -> adapters (vindex/build)  : immutable source record pointer + normalised entity + evidence row on parallel axes
   -> index (data/index)       : Parquet tables, schema in SCHEMA.md
   -> matcher (vindex/match_*) : normalised CHROM:POS:REF:ALT, structural equivalence, haplotypes, repeats, PGS (separate)
```

Rules that the code and tests enforce:
- Identity of a small variant is its normalised GRCh38 key; rsIDs and database identifiers are cross-references only.
- No composite score and no ordinal tier. Seven parallel evidence axes; a priority list is a query over the whole catalogue.
- Negative, null, uncertain and conflicting records are kept. Records that cannot be mapped are kept as `unmapped`/`rejected` rows with the reason.
- Absence from a variant-only VCF never establishes a reference genotype. Reference needs positive evidence (a gVCF block or a callable-region mask). Aligned reads (`--reads`) are asked about the small variants the calls leave unresolved and answer a different question, whether any fragment carries the indexed allele; the answer sits beside the matcher's status and never replaces it.
- Gene- and region-level evidence is never attached to a variant as variant evidence.
- Personal genomes stay local: nothing here reads from or writes to a network location for sample data, and `data/` is git-ignored.

Run the tests (from `variant_index/tests`): `python3 -m unittest discover` with `pysam`, `pyarrow` installed; set `VINDEX_TOOLS_BIN` to a directory
with `bcftools` and `plink2` to include the cross-checks against those tools (`environment.yml`).

## Workflow

```
# 1. build the index from the acquired sources (data/external) into data/index; every build starts from an empty directory
python3 variant_index/build_index.py                       # all adapters, in dependency order (--sources a,b to build a subset, --limit N for trials)

# 2. check it
python3 variant_index/verify_index.py [--sample 30000000]  # release linkage, manifest agreement, raw-record digests recomputed from the acquired files (vindex/resolvers.py holds one rule per source); the default samples 50,000 rows per release, a sample larger than the biggest release re-derives every row (27 minutes)
python3 variant_index/make_reports.py --markdown variant_index/BUILD_REPORT.md   # coverage, validation, target lifecycle, shared lineage, gap report (+ --previous OLD for a change log)
python3 variant_index/crosscheck_bcftools.py --clinvar-fields --out data/index/reports/clinvar_fields.json   # every ClinVar-VCF row of the index against the VCF read by bcftools (fields and normalised key)
python3 variant_index/verify_known_loci.py                 # facts established outside the index (an rsID's position, the direction of an association, a classification, a star allele's defining variant, curated gene and region assertions) read back from it
python3 variant_index/check_idempotence.py --limit 30000   # two builds in separate processes with different hash seeds; files, digests and the code that built them must be identical
python3 variant_index/check_requirements.py                # every clause of the project goal mapped to code, tests and a status (REQUIREMENTS.md)
python3 inventory/gap_audit.py                             # executable adversarial audit of the above (PHASE_AUDIT.md)

# 3. match one sample (local files only; results must go under data/ or outside the repository)
python3 variant_index/match_wgs.py --vcf SAMPLE.vcf.gz --out data/results/SAMPLE [--callable-bed HIGH_CONFIDENCE.bed] [--external-haplotypes calls.json] [--reads SAMPLE.cram]
python3 variant_index/make_views.py --results data/results/SAMPLE        # priority lists, trait and pathway groupings, report.md (views over the complete tables)
python3 variant_index/crosscheck_bcftools.py --vcf SAMPLE.vcf.gz --results data/results/SAMPLE   # the calls and the ClinVar evidence of those calls, recomputed with bcftools norm and compared
python3 variant_index/compare_runs.py RUN_A RUN_B           # two runs of one input (different hash seeds or days): every result table byte-identical, manifests equal apart from the clock

# 4. polygenic scores are separate models, never variant evidence
python3 variant_index/build_pgs_registry.py                # every score in the catalogue, whether its GRCh38 file was acquired, and its gates
python3 variant_index/score_pgs.py --vcf SAMPLE.vcf.gz --out data/results/SAMPLE/pgs --trait EFO:0004339 [--sample-ancestry European]

# whole-genome scale test on a public callset (GIAB HG002 benchmark; not anyone's personal genome)
python3 variant_index/fetch_public_testset.py

# the read-evidence resolver against that benchmark on public HG002 reads (read over HTTPS; see the script for the arguments and what the check does not show)
python3 variant_index/validate_reads.py --help
```

Reading the index from Python (the tables are Parquet; `Index` reads one chromosome partition at a time, so a query by key never loads the whole evidence table):

```python
import sys; sys.path.insert(0, 'variant_index')
from vindex.index_reader import Index
idx = Index.open('data/index')                                   # refuses an index written under another schema version
idx.evidence_for_variants(['1:1000:A:G'])                         # {normalised key: [evidence rows]}; each row says what it asserts, on which axis, from which source record
idx.table('entity_gene')                                          # any other table of SCHEMA.md
```

What each result table holds, and what every status means, is described in `vindex/wgs.py` and `SCHEMA.md`. A call with no index record stays in `small_calls` as
`NO_INDEX_RECORD`; an indexed high-value allele the file cannot answer for is listed in `unresolved_loci` with the evidence that would settle it, and, with `--reads`, with what the reads
show (`read_state`, the fragment counts, whether that settles the open question; `vindex/reads.py` states the policy, what a read counts as and what the result does not claim: it is not a
genotype). A CRAM is decoded against the FASTA given as `--fasta`, which has to be the reference it was aligned to (its contig names are mapped to the index's by `Reference.contig`); the
index of a BAM or CRAM is looked for beside it, or given as `--reads-index`. `reads_validation.json` holds the measured rates of the resolver on public reads of the benchmark sample, at
60x, 30x and 15x.

## Privacy and licences

- `data/` (public sources, the index, per-sample results) is git-ignored; `.gitignore` also blocks VCF / BAM / CRAM / FASTQ files anywhere in the tree. `match_wgs.py`
  and `score_pgs.py` refuse an output directory inside the tracked tree.
- A sample's calls are never sent anywhere: matching is a bulk local join against the index, never a per-variant query of an external service.
- Source licences travel with each release (`source_release.licence`, the catalogue records). Private status of this repository does not waive them.
