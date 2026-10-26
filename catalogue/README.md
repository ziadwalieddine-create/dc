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
| `goal_resources.json` | The resources and categories that the project goal names, each mapped to the catalogue records that cover it (state is computed in `GOAL_RESOURCES.md`, never asserted) |
| `build_catalogue.py` | Validates every record (vocabulary, evidence references, acquisition targets) and generates the documents below |
| `MASTER_INVENTORY.md`, `COVERAGE_MATRIX.md`, `GAP_MATRIX.md`, `REUSE_PLAN.md`, `ACQUIRED_FILES.md`, `ATTRIBUTION.md`, `EVIDENCE_PROFILES.md`, `MECHANISM_MATRIX.md`, `GOAL_RESOURCES.md`, `catalogue.json` | Generated outputs (`ACQUIRED_FILES.md` registers every acquired file with size, hash, source date, validation state and licence; `ATTRIBUTION.md` lists the licence and attribution terms of every source whose data we hold). Edit the sources, not these |

## Rules the records follow
- A fact is marked `live` only if a probe in `evidence/` returned it; otherwise the field is listed in `unverified_fields`. Nothing is taken from memory without that flag.
- Lifecycle stages are derived, not asserted, and evaluated per evidence dimension from the files that dimension actually uses (`supports` in `acquisition/targets.json`): **acquired** needs every supporting target verified in `acquisition/manifest.jsonl`, **validated** needs passing checks in `acquisition/validation_report.json`, **mapped** needs a passed join in `acquisition/compat_report.json` that exercised a supporting target. **Usable** additionally needs a licence verified live that permits research use; it does not mean an analysis pipeline exists. The four existing `scrape/` assets enter through `EXISTING_SOURCES` in `build_catalogue.py`.
- A catalogued or accessible source is never reported as complete. See the coverage matrix.
- Access blocks (registration walls, bot challenges, access denials) are recorded and not circumvented.
- Every record carries an `evidence_profile` with independent fields: the kind of evidence (direct measurement, perturbation, population observation, association, curated classification, case observation, prediction, reference annotation, structure), the resolution at which it refers to alleles (allele, haplotype, gene, region, protein residue), the mechanism axes it informs and on what basis, evidence strength, mechanistic confidence, phenotype implication, and what absence means. `build_catalogue.py --check` enforces that a gene- or region-level resource is never recorded as direct evidence about a specific allele or small variant, that every mechanism basis rests on a matching evidence kind, that a prediction resource carries no direct or curated mechanism unless it aggregates other sources, and that a record names only variant classes it lists.
- Absence is stated per source (`absence_means`): the absence of a study, a score or a record is never recorded as the absence of an effect. Hypotheses and competing explanations are kept apart: the mechanism axes in `vocab.json` list the competing explanations and what discriminates them (see `MECHANISM_MATRIX.md`, which also computes the variant-type by mechanism cells that have no usable empirical source).
- Overlaps between sources are declared once, by the source that depends on or restates the other, with a relation from the vocabulary: `identical_records` (count once), `derived_from`, `aggregates`, `independent_related` (compare, do not pool) and `may_conflict` (keep both). Identical records are removed when evidence is combined; independent experiments are never merged only because they concern the same gene or variant, and conflicting assertions are kept side by side.
- A statement in a record must be traceable to stored evidence, a generated report or the record's own `unverified_fields`; `derive/evidence_support.py` lists distinctive numbers and names that occur in none of these, and records are corrected until it is silent. Statements that could not be verified are worded as unverified in the text.
- The lifecycle stage `integrated` (a documented, tested query layer that returns records by variant, gene, protein or region) is defined and credited to no source: no such layer exists yet.

## Update procedure
```
python3 catalogue/probe.py catalogue/probes.json            # refresh evidence (append-only)
python3 catalogue/build_catalogue.py --check                 # validate records
python3 -m venv .venv && .venv/bin/pip install -r acquisition/requirements.txt    # pyarrow, pysam
python3 acquisition/acquire.py --tier A                      # download and checksum-verify
python3 acquisition/acquire.py --ids alphamissense_hg38 dgv_grch38_variants_2025_12_01 gnomad_sv_v4_1_sites_vcf   # tier B
python3 acquisition/extract_dump.py                          # parts of the MaveDB dump the tests read
.venv/bin/python acquisition/fetch_per_protein.py            # AlphaFold DB and InterPro for the MaveDB proteins
.venv/bin/python acquisition/fetch_gnomad_regional.py        # gnomAD v4.1 slices at MaveDB-mapped positions (remote tabix)
.venv/bin/python acquisition/validate.py                     # check files against their source descriptions
.venv/bin/python acquisition/compat_tests.py                 # joins to existing assets
python3 catalogue/build_catalogue.py                         # regenerate the documents
python3 inventory/audit_existing.py --check                  # confirm scrape/ is unchanged
python3 inventory/check_repo_policy.py                       # size, personal-data, credential and visibility check
python3 inventory/repository_audit.py                        # branches, documents, scripts, manifests, dataset lifecycle and 30+ documentation claims recomputed
python3 catalogue/derive/evidence_support.py                 # statements that no stored evidence supports
python3 inventory/gap_audit.py                               # adversarial audit of the phase requirements (writes inventory/PHASE_AUDIT.md)
```
`derive/` scripts need the evidence they read (`opentargets_*.py`, `ot_interaction_counts.py`, `ot_credible_set_counts.py`, `openebench_counts.py`, `vep_plugins.py`), pyarrow (the two `ot_*_counts.py` scripts) or pysam (`gnomad_vcf_header.py`); `local_file.py` and `pdf_text.py` store a local file or a PDF as evidence with its provenance, and `evidence_support.py` checks the records against the stored evidence.
