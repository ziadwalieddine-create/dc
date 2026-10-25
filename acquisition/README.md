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
| `dump_equivalence_sample.json` | Fixed list of our score sets whose tables are also extracted from the dump, to check dump and API agree |
| `mavedb_superseded/` | Additive set of 52 superseded MaveDB score sets that the public search hides (see `AUDIT_LOG.md`) |
| `requirements.txt` | `pyarrow`, needed only for the Open Targets Parquet datasets |

```
python3 acquisition/acquire.py --list
python3 acquisition/acquire.py --tier A                       # all tier-A targets
python3 acquisition/acquire.py --ids alphamissense_hg38 dgv_grch38_variants_2025_12_01      # tier-B files by id
python3 acquisition/extract_dump.py                            # after mavedb_zenodo_dump_v4 is downloaded
python3 -m venv .venv && .venv/bin/pip install -r acquisition/requirements.txt
.venv/bin/python acquisition/validate.py                       # every downloaded target
.venv/bin/python acquisition/compat_tests.py                   # about 25 minutes; --skip orcs ot shortens it
python3 catalogue/build_catalogue.py                           # regenerate the inventory, matrices and plan
```

Order matters once: the compatibility tests read the dump extraction, so run `acquire.py` for `mavedb_zenodo_dump_v4`, then `extract_dump.py`, before `compat_tests.py`. `extract_dump.py --out DIR` writes elsewhere, which is how its output was checked byte for byte against the first extraction (1,297 files, no differences).

## Notes
- The MaveDB bulk dump (`mavedb_zenodo_dump_v4`, 1.9 GB) is kept as downloaded. `extract_dump.py` extracts only `main.json`, `README.md`, `LICENSE.txt`, the 1,170 per-set `annotations.csv` files (708 MB), and the scores/counts tables of the 52 superseded sets and of the 31-set comparison sample. The 9.5 GB of `mapped/` JSON and 3.4 GB of `va/` VA-Spec files are not extracted. The archive is untrusted input: members are selected by name and written under paths built from the base name only.
- Download sizes stated on web pages can be misread (one was); the downloader treats only source-published checksums and `Content-Length` as authoritative. Where a source publishes no checksum (Open Targets, DGV) the manifest says so and records our own SHA-256.
- Tolerances in validators are derived, not guessed. Example: AlphaMissense prints scores to four decimals but assigns classes from the unrounded value, so `validate.py` accepts a different class only at the two documented boundaries (0.34 and 0.564) and reports how many rows needed the allowance (3,838 of 71,697,556, none elsewhere).
- Data are not stored in git. They are reproducible from `targets.json` plus the manifest checksums. `manifest_files/` holds the per-part SHA-256 tables for directory targets.
- `inventory/check_repo_policy.py` (repository root level) checks that nothing tracked is oversized, personal-genome-like or credential-like, and prints the live repository visibility.
