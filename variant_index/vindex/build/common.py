"""Build context, row sinks and evidence-row construction shared by all adapters.

Provenance: every evidence row carries release_id (a row of source_release, pinned by the acquisition manifest), the source's own
accession, a locator inside the acquired file and the SHA-256 of the raw record text. `verify_provenance` re-opens records through
the locator and recomputes the digest, so an index row can be shown unchanged from its source.
Determinism: rows are written in source order, ids are hashes of source identity (never of time or position in the output), and no
timestamp enters a data file; `index_manifest.json` holds the timestamps and the content digests.
"""
import datetime
import email.utils
import hashlib
import json
import os
import re

import pyarrow.parquet as pq

from .. import schema

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..'))      # variant_index/
REPO = os.path.dirname(ROOT)
DATA = os.path.join(REPO, 'data', 'external')
OUT_DEFAULT = os.path.join(REPO, 'data', 'index')

CHROMS = list(schema.PRIMARY_CONTIGS)
ALL = object()      # BuildContext.want(namespace) with no identifiers: every offered identifier is kept


def sha256_text(s):
    return hashlib.sha256(s.encode('utf-8') if isinstance(s, str) else s).hexdigest()


_DUPLICATED_PREFIX = re.compile(r'^([A-Za-z0-9.]+):\1:')
_UNDERSCORE_PREFIX = re.compile(r'^([A-Za-z][A-Za-z0-9.]*)_(?=[A-Za-z0-9])')
_SPACE_AFTER_PREFIX = re.compile(r'^([A-Za-z][A-Za-z0-9._]*):\s+(\S+)$')       # 'OMIM: 612621' (ClinVar); a name after a prefix ('HP:gait ataxia') is left as written


def curie(x):
    """One spelling for an ontology or database identifier, PREFIX:LOCAL. The sources disagree: Open Targets writes MONDO_0005393 and Orphanet_199, ClinVar's VCF repeats the prefix
    (MONDO:MONDO:0010178) and spells the Human Phenotype Ontology out (Human_Phenotype_Ontology:HP:0001250), the GWAS Catalog's URIs end in EFO_0004339. Left as written, one disease is
    several entries in every grouping and a contradiction between two sources goes unseen. Identifiers that already have the form are returned as they are."""
    if not isinstance(x, str):
        return x
    s = re.sub(r'^Human_Phenotype_Ontology:\s*', '', x.strip())
    s = _DUPLICATED_PREFIX.sub(r'\1:', s)
    s = _SPACE_AFTER_PREFIX.sub(r'\1:\2', s)
    if ':' not in s:
        s = _UNDERSCORE_PREFIX.sub(r'\1:', s)
    return s


def curies(ids):
    """curie() of every identifier, in order and without repeats; None stays None."""
    if ids is None:
        return None
    return list(dict.fromkeys(c for c in (curie(i) for i in ids) if c))


def split_identifiers(ids):
    """curies(ids) as (PREFIX:LOCAL identifiers, values that are not identifiers). A source can type a name where it has identifiers (Open Targets lists 'gastric adenocarcinoma' among
    the EFO and MONDO ids of one FinnGen study); such a value is not grouped as an identifier, and evidence_row keeps it in extra_json."""
    cs = curies(ids) or []
    return [c for c in cs if ':' in c], [c for c in cs if ':' not in c]


_ISO_PREFIX = re.compile(r'^(\d{4})-(\d{2})-(\d{2})(?:$|[T\s-])')


def iso_date(x):
    """The date of an evaluation as YYYY-MM-DD. The sources write it as 2023-05-18, 2023-05-18T16:00:00.000Z, 'Mon, 05 Jun 2023 12:00:00 -0400' (ClinGen actionability) and, in 21 ClinVar
    records, 2021-01-05-12:30. The date is the one written, in the zone written. Text that is not read as a date stays as it is, so nothing is invented."""
    if not isinstance(x, str):
        return x
    s = x.strip()
    m = _ISO_PREFIX.match(s)
    try:
        if m:
            return datetime.date(int(m.group(1)), int(m.group(2)), int(m.group(3))).isoformat()
        if re.match(r'^[A-Za-z]{3},\s', s):
            return email.utils.parsedate_to_datetime(s).date().isoformat()
        for fmt in ('%b %d, %Y', '%B %d, %Y'):
            try:
                return datetime.datetime.strptime(s, fmt).date().isoformat()
            except ValueError:
                pass
    except (ValueError, TypeError, OverflowError):
        pass
    return x


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
        self.owners = {}        # namespace -> adapters that asked for it (the offered entity ids are dropped once every one of them has run)
        self.found = {}         # (namespace, identifier) -> entity id, filled by the adapter that owns the namespace

    def want(self, namespace, identifiers=None):
        """Ask the adapters that run before this one to remember the entity ids they map identifiers of `namespace` to. identifiers=None asks for all of them."""
        self.owners.setdefault(namespace, []).append(self.adapter)
        if identifiers is None:
            self.interest[namespace] = ALL
        elif self.interest.get(namespace) is not ALL:
            self.interest.setdefault(namespace, set()).update(str(i) for i in identifiers)

    def offer(self, namespace, identifier, entity_id):
        """Called by an adapter for every identifier it maps to an entity; kept only when another adapter asked for it."""
        w = self.interest.get(namespace)
        if w is ALL or (w and identifier in w):
            self.found.setdefault((namespace, identifier), []).append(entity_id)

    def release_finished(self, finished):
        """Forget the identifier -> entity pairs of every namespace whose requesting adapters have all run: millions of ClinVar variation ids are only needed until the
        adapters that join on them are done."""
        for ns in [n for n, owners in self.owners.items() if n in self.interest and all(any(o is f for f in finished) for o in owners)]:
            for k in [k for k in self.found if k[0] == ns]:
                del self.found[k]
            del self.interest[ns]

    def st(self, name=None):
        name = name or (type(self.adapter).__name__ if self.adapter else '_')
        return self.stats.setdefault(name, Stats())

    # ---- releases -------------------------------------------------------------------------------------------------------------
    def register_release(self, source_id, target_id, version, file=None):
        m = self.manifest.get(target_id, {})
        t = self.targets.get(target_id, {})
        rid = f'{source_id}/{target_id}@{version}'
        if rid not in self.releases:
            # a target written from another acquired file carries no URL of its own: the retrieval URL is that of the file it was extracted from
            url = t.get('url') or m.get('url') or self.targets.get(t.get('derived_from'), {}).get('url') or ''
            self.releases[rid] = dict(release_id=rid, source_id=source_id, target_id=target_id, version=version, file=file or t.get('dest', ''),
                                      size_bytes=m.get('bytes'), sha256=m.get('sha256'), url=url,
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
        return self.sink(os.path.join('evidence', level, schema.partition_of(chrom), f'{source_id}.parquet'), 'evidence')

    def emit_evidence(self, row):
        """Validate and write one evidence row. Returns the row."""
        if self.validate_rows:
            bad = schema.validate_evidence_row(row)
            if bad:
                raise ValueError(f"invalid evidence row {row.get('source_id')}:{row.get('accession')}: {bad}")
        if row['level'] == 'variant' and row['mapping_status'] == 'mapped':
            c, p = row['subject_id'].split(':')[:2]
            if row.get('chrom') in (None, ''):
                row['chrom'] = c                                # the partition a matcher reads is the key's chromosome; a mapped row can never be filed anywhere else
            elif row['chrom'] != c:
                raise ValueError(f"mapped variant row {row.get('source_id')}:{row.get('accession')} has chrom {row['chrom']} but key {row['subject_id']}")
            if row.get('pos') is None:
                row['pos'] = int(p)
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
    if row.get('condition_ids') is not None:
        row['condition_ids'], words = split_identifiers(row['condition_ids'])
        if words:
            extra = json.loads(row['extra_json']) if row.get('extra_json') else {}
            if not isinstance(extra, dict):
                extra = {'extra': extra}
            extra['condition_ids_not_identifiers'] = words
            row['extra_json'] = json.dumps(extra, sort_keys=True)
    if row.get('evaluated') is not None:
        row['evaluated'] = iso_date(row['evaluated'])
    row['evidence_id'] = evidence_id(release_id, accession, subject_id or row.get('unmapped_key') or '', axis, claim_type, discriminator)
    return row


def retrieved(ctx, target):
    """Version label for a file that carries no release number: the day it was fetched (from the acquisition manifest)."""
    return 'retrieved-' + ((ctx.manifest.get(target) or {}).get('fetched_utc') or 'undated')[:10]


def utc_now():
    return datetime.datetime.now(datetime.timezone.utc).strftime('%Y-%m-%dT%H:%M:%SZ')
