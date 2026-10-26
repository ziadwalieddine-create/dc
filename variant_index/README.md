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
- Absence from a variant-only VCF never establishes a reference genotype. Reference needs positive evidence (gVCF block, callable-region mask, reads).
- Gene- and region-level evidence is never attached to a variant as variant evidence.
- Personal genomes stay local: nothing here reads from or writes to a network location for sample data, and `data/` is git-ignored.

Run the tests (from `variant_index/tests`): `python3 -m unittest discover` with `pysam`, `pyarrow` installed; set `VINDEX_TOOLS_BIN` to a directory
with `bcftools` and `plink2` to include the cross-checks against those tools (`environment.yml`).
