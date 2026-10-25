# Evidence-source catalogue

A structured, evidence-backed catalogue of external human-genomics resources, built to extend (not replace) the assets under `scrape/`.

| Path | Purpose |
|---|---|
| `sources/*.json` | One curated record per resource: identity, version/date, scope, variant classes, assemblies, formats, identifiers, quantitative semantics, licence, update cadence, coverage limits, overlap with existing assets, reuse decision, acquisition targets |
| `evidence/` | Provenance for the records: excerpts of what each URL returned (`evidence_index.jsonl` keeps URL, fetch time, HTTP status, sha256, headers; append-only) |
| `probes.json` + `probe.py` | The list of URLs probed and the collector (plain HTTP, ranged reads, POST, or Playwright rendering for JavaScript-only pages) |
| `derive/` | Small scripts that turn evidence into derived facts (Open Targets dataset inventory, registry counts) |
| `vocab.json` | Controlled vocabulary: evidence dimensions, variant classes, assemblies, access states, reuse actions, lifecycle stages |
| `gaps.json` | Curated gap analysis and next actions per dimension |
| `build_catalogue.py` | Validates every record (vocabulary, evidence references, acquisition targets) and generates the documents below |
| `MASTER_INVENTORY.md`, `COVERAGE_MATRIX.md`, `GAP_MATRIX.md`, `REUSE_PLAN.md`, `catalogue.json` | Generated outputs. Edit the sources, not these |

## Rules the records follow
- A fact is marked `live` only if a probe in `evidence/` returned it; otherwise the field is listed in `unverified_fields`. Nothing is taken from memory without that flag.
- Lifecycle stages are derived, not asserted: **acquired** needs a verified entry in `acquisition/manifest.jsonl`, **validated** needs passing checks in `acquisition/validation_report.json`, **mapped** needs a demonstrated join in `acquisition/compat_report.json`.
- A catalogued or accessible source is never reported as complete. See the coverage matrix.
- Access blocks (registration walls, bot challenges, access denials) are recorded and not circumvented.

## Update procedure
```
python3 catalogue/probe.py catalogue/probes.json            # refresh evidence (append-only)
python3 catalogue/build_catalogue.py --check                 # validate records
python3 acquisition/acquire.py --tier A                      # download and checksum-verify
python3 acquisition/validate.py                              # check files against their source descriptions
python3 acquisition/compat_tests.py                          # joins to existing assets
python3 catalogue/build_catalogue.py                         # regenerate the documents
python3 inventory/audit_existing.py --check                  # confirm scrape/ is unchanged
```
