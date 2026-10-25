# Acquisition, validation and compatibility tests

Reproducible download of external resources into `data/external/` (git-ignored; never committed). Licences differ per source, so a private repository would not make redistribution safe; see `catalogue/sources/*.json`.

| File | Purpose |
|---|---|
| `targets.json` | What to download: URL, destination, tier (A = fetch now, B = large/optional), expected size, source-published checksum (MD5, base64 MD5, Unix `sum`) |
| `acquire.py` | Resumable downloader. Verifies size and checksum against the source where it publishes one, rejects HTML pages saved as data files, records everything in `manifest.jsonl` |
| `manifest.jsonl` | One record per downloaded file: bytes, md5, sha256, source checksum comparison, last-modified, ETag |
| `validate.py` | Per-file checks against each source's own description (headers, record counts, counts shown on the source website). Output: `validation_report.json` |
| `compat_tests.py` | Identifier and genomic-coordinate joins between the new resources and the existing MaveDB human set. Output: `compat_report.json` |
| `mavedb_superseded/` | Additive set of 52 superseded MaveDB score sets that the public search hides (see its notes in `AUDIT_LOG.md`) |

```
python3 acquisition/acquire.py --list
python3 acquisition/acquire.py --tier A
python3 acquisition/validate.py
python3 acquisition/compat_tests.py
```

## Notes
- The MaveDB bulk dump (`mavedb_zenodo_dump_v4`, 1.9 GB) is kept as downloaded. Only `main.json` and the per-set `annotations.csv` files (708 MB) were extracted, plus the scores/counts of the 52 superseded sets and a 31-set comparison sample. The 9.5 GB of `mapped/` JSON and 3.4 GB of `va/` VA-Spec files were not extracted.
- Download sizes stated on web pages can be misread (one was); the downloader treats only source-published checksums and `Content-Length` as authoritative.
- Data are not stored in git. They are reproducible from `targets.json` plus the manifest checksums.
