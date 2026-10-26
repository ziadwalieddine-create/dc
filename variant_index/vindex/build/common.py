"""Build context, row sinks and evidence-row construction shared by all adapters.

Provenance: every evidence row carries release_id (a row of source_release, pinned by the acquisition manifest), the source's own
accession, a locator inside the acquired file and the SHA-256 of the raw record text. `verify_provenance` re-opens records through
the locator and recomputes the digest, so an index row can be shown unchanged from its source.
Determinism: rows are written in source order, ids are hashes of source identity (never of time or position in the output), and no
timestamp enters a data file; `index_manifest.json` holds the timestamps and the content digests.
"""
import datetime
import hashlib
import json
import os

import pyarrow as pa
import pyarrow.parquet as pq

from .. import schema

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..'))      # variant_index/
REPO = os.path.dirname(ROOT)
DATA = os.path.join(REPO, 'data', 'external')
OUT_DEFAULT = os.path.join(REPO, 'data', 'index')

CHROMS = [str(i) for i in range(1, 23)] + ['X', 'Y', 'MT']
ALL = object()      # BuildContext.want(namespace) with no identifiers: every offered identifier is kept


def sha256_text(s):
    return hashlib.sha256(s.encode('utf-8') if isinstance(s, str) else s).hexdigest()


def evidence_id(release_id, accession, subject_id, axis, claim_type, discriminator=''):
    return hashlib.sha1(f'{release_id}|{accession}|{subject_id}|{axis}|{claim_type}|{discriminator}'.encode()).hexdigest()[:20]


class Budget:
    """Rows buffered across all sinks of one build. Evidence is partitioned by chromosome, so many sinks fill slowly and together hold far more than any one
    of them; when the total reaches the limit every sink writes what it has (a row group per flush)."""

    def __init__(self, limit=300_000):
        self.limit, self.n, self.sinks = limit, 0, []

    def flush_all(self):
        for s in self.sinks:
            s.flush()
        self.n = 0


class Sink:
    """Buffered Parquet writer for one output file."""

    def __init__(self, path, table, chunk=200_000, budget=None):
        os.makedirs(os.path.dirname(path), exist_ok=True)
        self.path, self.table, self.chunk, self.budget = path, table, chunk, budget
        self.sch = schema.TABLES[table]
        self.buf, self.n = [], 0
        self.w = None
        if budget is not None:
            budget.sinks.append(self)

    def add(self, row):
        self.buf.append(row)
        if self.budget is not None:
            self.budget.n += 1
            if self.budget.n >= self.budget.limit:
                self.budget.flush_all()
        if len(self.buf) >= self.chunk:
            self.flush()

    def flush(self):
        if not self.buf:
            return
        t = schema.rows_to_table(self.table, self.buf)
        if self.w is None:
            self.w = pq.ParquetWriter(self.path, self.sch, compression='zstd')
        self.w.write_table(t)
        self.n += len(self.buf)
        if self.budget is not None:
            self.budget.n = max(0, self.budget.n - len(self.buf))
        self.buf = []

    def close(self):
        self.flush()
        if self.w is None:           # keep an empty but valid file so readers see the table
            pq.write_table(self.sch.empty_table(), self.path, compression='zstd')
        else:
            self.w.close()


class Stats:
    def __init__(self):
        self.c = {}

    def inc(self, key, n=1):
        self.c[key] = self.c.get(key, 0) + n


class BuildContext:
    def __init__(self, out_dir, reference, manifest_path=None, validate_rows=True):
        self.out = out_dir
        self.ref = reference
        self.validate_rows = validate_rows
        self.sinks = {}
        self.budget = Budget()
        self.stats = {}
        self.releases = {}
        self.problems = []
        self.manifest = {}
        mp = manifest_path or os.path.join(REPO, 'acquisition', 'manifest.jsonl')
        if os.path.exists(mp):
            with open(mp) as f:
                for line in f:
                    if line.strip():
                        r = json.loads(line)
                        self.manifest[r['id']] = r
        self.targets = {}
        tp = os.path.join(REPO, 'acquisition', 'targets.json')
        if os.path.exists(tp):
            with open(tp) as f:
                self.targets = {t['id']: t for t in json.load(f)}
        self.adapter = None
        self.interest = {}      # namespace -> identifiers some adapter wants resolved to entities
        self.found = {}         # (namespace, identifier) -> entity id, filled by the adapter that owns the namespace

    def want(self, namespace, identifiers=None):
        """Ask the adapters that run before this one to remember the entity ids they map identifiers of `namespace` to. identifiers=None asks for all of them."""
        if identifiers is None:
            self.interest[namespace] = ALL
        elif self.interest.get(namespace) is not ALL:
            self.interest.setdefault(namespace, set()).update(str(i) for i in identifiers)

    def offer(self, namespace, identifier, entity_id):
        """Called by an adapter for every identifier it maps to an entity; kept only when another adapter asked for it."""
        w = self.interest.get(namespace)
        if w is ALL or (w and identifier in w):
            self.found.setdefault((namespace, identifier), []).append(entity_id)

    def st(self, name=None):
        name = name or (type(self.adapter).__name__ if self.adapter else '_')
        return self.stats.setdefault(name, Stats())

    # ---- releases -------------------------------------------------------------------------------------------------------------
    def register_release(self, source_id, target_id, version, file=None):
        m = self.manifest.get(target_id, {})
        t = self.targets.get(target_id, {})
        rid = f'{source_id}/{target_id}@{version}'
        if rid not in self.releases:
            self.releases[rid] = dict(release_id=rid, source_id=source_id, target_id=target_id, version=version, file=file or t.get('dest', ''),
                                      size_bytes=m.get('bytes'), sha256=m.get('sha256'), url=t.get('url', m.get('url', '')),
                                      retrieved_utc=m.get('fetched_utc'), licence=t.get('licence', ''), licence_class='', manifest_line=None)
            if not m:
                self.problems.append(f'no manifest record for target {target_id}')
        return rid

    # ---- sinks ----------------------------------------------------------------------------------------------------------------
    def sink(self, relpath, table):
        key = (relpath, table)
        if key not in self.sinks:
            self.sinks[key] = Sink(os.path.join(self.out, relpath), table, budget=self.budget)
        return self.sinks[key]

    def evidence_sink(self, level, chrom, source_id):
        c = chrom if chrom in CHROMS else 'NA'      # not '_': Arrow dataset discovery silently skips paths that start with '_' or '.'
        return self.sink(os.path.join('evidence', level, c, f'{source_id}.parquet'), 'evidence')

    def emit_evidence(self, row):
        """Validate and write one evidence row. Returns the row."""
        if self.validate_rows:
            bad = schema.validate_evidence_row(row)
            if bad:
                raise ValueError(f"invalid evidence row {row.get('source_id')}:{row.get('accession')}: {bad}")
        self.evidence_sink(row['level'], row.get('chrom') or '', row['source_id']).add(row)
        self.st().inc(f"evidence:{row['level']}:{row['mapping_status']}")
        return row

    def emit(self, table, row, name=None):
        self.sink(os.path.join(table, f'{name or (self.adapter.source_id if self.adapter else "x")}.parquet'), table).add(row)

    def close(self):
        for s in self.sinks.values():
            s.close()


def evidence_row(release_id, source_id, accession, locator, raw, level, subject_id, axis, kind, claim_type, claim_value,
                 effect_assertion, mapping_status='mapped', mapping_reason=None, discriminator='', **kw):
    row = dict(source_id=source_id, release_id=release_id, accession=accession, locator=locator, raw_sha256=sha256_text(raw), level=level,
               subject_id=subject_id, axis=axis, evidence_kind=kind, claim_type=claim_type, claim_value=claim_value, effect_assertion=effect_assertion,
               mapping_status=mapping_status, mapping_reason=mapping_reason)
    row.update(kw)
    row['evidence_id'] = evidence_id(release_id, accession, subject_id or row.get('unmapped_key') or '', axis, claim_type, discriminator)
    return row


def retrieved(ctx, target):
    """Version label for a file that carries no release number: the day it was fetched (from the acquisition manifest)."""
    return 'retrieved-' + ((ctx.manifest.get(target) or {}).get('fetched_utc') or 'undated')[:10]


def utc_now():
    return datetime.datetime.now(datetime.timezone.utc).strftime('%Y-%m-%dT%H:%M:%SZ')
