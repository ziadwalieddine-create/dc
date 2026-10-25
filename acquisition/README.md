# Acquisition, validation and compatibility tests

Reproducible download of external resources into `data/external/` (git-ignored; never committed). Licences differ per source, so a private repository would not make redistribution safe; see `catalogue/sources/*.json`.

| File | Purpose |
|---|---|
| `targets.json` | What to download: URL, destination, tier (A = fetch now, B = large/optional), expected size, source-published checksum (MD5, base64 MD5, Unix `sum`), and `supports` (the evidence dimensions the file can support, used to credit coverage only for what was actually fetched) |
| `acquire.py` | Resumable downloader. Verifies size and checksum against the source where it publishes one, rejects HTML pages saved as data files, records everything in `manifest.jsonl`. `"kind": "http_dir"` targets (Open Targets Parquet directories) are fetched part by part with parallel connections, checked against each part's Content-Length and Parquet framing, and listed with SHA-256 values in `manifest_files/<target>.tsv` |
| `manifest.jsonl` | One record per downloaded file or directory: bytes, md5, sha256, source checksum comparison, last-modified, ETag |
| `validate.py` | Per-file checks against each source's own description (headers, record counts, documented thresholds, counts shown on the source website; Parquet footers, schemas, identifier formats and value ranges). Output: `validation_report.json` |
| `compat_tests.py`, `compat_opentargets.py` | Identifier and genomic-coordinate joins between the acquired resources and the existing MaveDB human set. Each test names the targets it exercised, so a join is credited only to the files it used. Output: `compat_report.json` |
| `extract_dump.py` | Deterministic extraction of the parts of the MaveDB bulk dump that the tests read (see below) |
| `fetch_per_protein.py` | AlphaFold DB (canonical model, per-residue confidence; isoform-only proteins are recorded, not stored) and InterPro matches for the UniProt proteins named by the MaveDB score sets. Per-file SHA-256 table in `manifest_files/` |
| `fetch_gnomad_regional.py`, `mavedb_snvs.py` | gnomAD v4.1 exome and genome records at the genomic positions of MaveDB-mapped SNVs, read by remote tabix (no whole-file download); every region is read until two reads agree. `mavedb_snvs.py` derives the SNV alleles and merged regions |
| `fetch_reference_regional.py` | GRCh38 reference sequence for the same merged regions (1,546 regions on 54 sequences), read by HTTP range requests from the Ensembl release-116 toplevel FASTA; every region is read until two reads agree. The saved `.fai` is kept in `index_cache/` for validation |
| `compat_reference.py` | Ties the mirrored FASTA and the Ensembl release-116 GTF to the other coordinate-bearing files: remote reads equal the local file, REF alleles of every ClinVar record and every AlphaMissense position, start and stop codons, MANE, Open Targets locations, the amino-acid change implied by sequence and annotation against AlphaMissense, and the gnomAD LoF curation. Test code, not an annotation engine |
| `intervals.py`, `test_helpers.py` | The exact interval-stabbing query behind every overlap count, and unit tests for it and for the other small helpers (`python3 -m unittest acquisition.test_helpers` from the repository root; no data needed) |
| `dump_equivalence_sample.json` | Fixed list of our score sets whose tables are also extracted from the dump, to check dump and API agree |
| `mavedb_superseded/` | Additive set of 52 superseded MaveDB score sets that the public search hides (see `AUDIT_LOG.md`) |
| `requirements.txt` | `pyarrow` (Open Targets Parquet datasets) and `pysam` (remote slices of gnomAD and of the reference FASTA, and the local FASTA in `compat_reference.py`) |

```
python3 acquisition/acquire.py --list
python3 acquisition/acquire.py --tier A                       # all tier-A targets
python3 acquisition/acquire.py --ids alphamissense_hg38 dgv_grch38_variants_2025_12_01      # tier-B files by id
python3 acquisition/extract_dump.py                            # after mavedb_zenodo_dump_v4 is downloaded
.venv/bin/python acquisition/fetch_per_protein.py              # AlphaFold DB and InterPro per protein (about 30 minutes, polite request rate)
.venv/bin/python acquisition/fetch_gnomad_regional.py          # gnomAD slices at MaveDB-mapped positions (about 20 minutes)
python3 acquisition/acquire.py --ids ensembl_grch38_dna_toplevel ensembl_grch38_dna_toplevel_fai ensembl_grch38_dna_toplevel_gzi ensembl_grch38_gtf_116   # reference FASTA (808 MB) with indexes, and the 141 MB GTF; BSD checksums from the Ensembl CHECKSUMS files
.venv/bin/python acquisition/fetch_reference_regional.py       # reference sequence at the MaveDB-mapped regions by remote reads (about 5 minutes)
python3 -m unittest acquisition.test_helpers                   # 17 tests, no data needed
python3 -m venv .venv && .venv/bin/pip install -r acquisition/requirements.txt
.venv/bin/python acquisition/validate.py                       # every downloaded target
.venv/bin/python acquisition/compat_tests.py                   # about 30 minutes; --skip orcs ot ref shortens it
python3 catalogue/build_catalogue.py                           # regenerate the inventory, matrices and plan
```

Order matters once: the compatibility tests read the dump extraction, so run `acquire.py` for `mavedb_zenodo_dump_v4`, then `extract_dump.py`, before `compat_tests.py`. `extract_dump.py --out DIR` writes elsewhere, which is how its output was checked byte for byte against the first extraction (1,297 files, no differences).

## Notes
- The MaveDB bulk dump (`mavedb_zenodo_dump_v4`, 1.9 GB) is kept as downloaded. `extract_dump.py` extracts only `main.json`, `README.md`, `LICENSE.txt`, the 1,170 per-set `annotations.csv` files (708 MB), and the scores/counts tables of the 52 superseded sets and of the 31-set comparison sample. The 9.5 GB of `mapped/` JSON and 3.4 GB of `va/` VA-Spec files are not extracted. The archive is untrusted input: members are selected by name and written under paths built from the base name only.
- Download sizes stated on web pages can be misread (one was); the downloader treats only source-published checksums and `Content-Length` as authoritative. Where a source publishes no checksum (Open Targets, DGV) the manifest says so and records our own SHA-256.
- Tolerances in validators are derived, not guessed. Example: AlphaMissense prints scores to four decimals but assigns classes from the unrounded value, so `validate.py` accepts a different class only at the two documented boundaries (0.34 and 0.564) and reports how many rows needed the allowance (3,838 of 71,697,556, none elsewhere).
- Data are not stored in git. They are reproducible from `targets.json` plus the manifest checksums. `manifest_files/` holds the per-part SHA-256 tables for directory targets.
- `inventory/check_repo_policy.py` (repository root level) checks that nothing tracked is oversized, personal-genome-like or credential-like, and prints the live repository visibility.
