"""Provenance checks on a built index.

Three questions, each answered from the index and the acquired sources and never from memory of a previous run:

  1. Linkage     every release_id used by any table is a row of source_release (no evidence from an unnamed release).
  2. Release     every release names its acquisition target; the manifest record exists, carries a digest, URL and licence; the file on disk has the recorded size
                 (and, on request, the recorded SHA-256), so the immutable source the index points to is the one that was acquired.
  3. Records     for a deterministic sample of evidence rows of each release that has a resolver, the locator is opened in the acquired file and the SHA-256 of the raw
                 record text is recomputed; it must equal the raw_sha256 stored in the row. A release without a resolver is listed as such: it is not counted as verified.
"""
import gzip
import hashlib
import io
import json
import os
import random
import re
import zipfile

import pyarrow.compute as pc
import pyarrow.dataset as pds
import pyarrow.parquet as pq

from . import schema
from .build.clinvar_vcv import canon as vcv_canon
from .build.common import DATA, REPO

LOCATOR_LINE = re.compile(r'L(\d+)$')
LOCATOR_ROW = re.compile(r'^(?P<table>[a-z_]+)/(?P<part>[^#]+)#(?P<i>\d+)$')


def sha256_text(s):
    return hashlib.sha256(s.encode('utf-8')).hexdigest()


def _open(path, newline):
    return gzip.open(path, 'rt', encoding='utf-8', newline=newline) if str(path).endswith(('.gz', '.bgz')) else open(path, encoding='utf-8', newline=newline)


def lines_iter(path, wanted, newline):
    """{line number (1 = header): text without its newline} read in one pass through a (compressed) text file."""
    out, last = {}, max(wanted)
    with _open(path, newline) as f:
        for ln, line in enumerate(f, start=1):
            if ln in wanted:
                out[ln] = line.rstrip('\n')
            if ln >= last:
                break
    return out


def lines_split(path, wanted, how):
    with _open(path, '') as f:
        text = f.read()
    lines = text.split('\n') if how == 'split' else text.splitlines()
    return {n: lines[n - 1] for n in wanted if 0 < n <= len(lines)}


def lines_gwas(path, wanted):
    """The association table inside the GWAS Catalog zip, rows padded to the header width as the adapter does."""
    out, last = {}, max(wanted)
    with zipfile.ZipFile(path) as z, z.open(z.namelist()[0]) as raw:
        f = io.TextIOWrapper(raw, encoding='utf-8', errors='replace', newline='')
        hdr = None
        for ln, line in enumerate(f, start=1):
            row = line.rstrip('\r\n').split('\t')
            if hdr is None:
                hdr = row
            if ln in wanted:
                out[ln] = '\t'.join(row + [''] * (len(hdr) - len(row)))
            if ln >= last:
                break
    return out


def resolve_lines(kind, path, locators):
    wanted = {int(LOCATOR_LINE.search(l).group(1)) for l in locators if LOCATOR_LINE.search(l)}
    if not wanted:
        return {}
    got = {'iter': lambda: lines_iter(path, wanted, None), 'iter_raw': lambda: lines_iter(path, wanted, ''), 'split': lambda: lines_split(path, wanted, 'split'),
           'splitlines': lambda: lines_split(path, wanted, 'splitlines'), 'gwas': lambda: lines_gwas(path, wanted)}[kind]()
    return {l: got.get(int(LOCATOR_LINE.search(l).group(1))) for l in locators if LOCATOR_LINE.search(l)}


def resolve_parquet_rows(root, locators):
    cache, out = {}, {}
    for l in locators:
        m = LOCATOR_ROW.match(l)
        if not m:
            continue
        p = os.path.join(root, m.group('table'), m.group('part'))
        if not os.path.exists(p):
            continue
        if p not in cache:
            cache[p] = pq.read_table(p).to_pylist()
        rows = cache[p]
        i = int(m.group('i'))
        out[l] = vcv_canon(rows[i]) if i < len(rows) else None
    return out


# target id -> how its raw record text is recovered from the acquired file (the same rule the adapter used to compute raw_sha256)
RESOLVERS = {
    'clinvar_vcf_grch38': 'iter',                      # one VCF line, newline stripped
    'gnomad_constraint_v4_1_2': 'iter_raw',
    'gnomad_v4_1_joint_regional': 'iter_raw',
    'gnomad_v4_1_exomes_regional': 'iter_raw',
    'ensembl_vep_mavedb_variants': 'iter_raw',
    'gnomad_v4_1_lof_curation_incomplete_penetrance': 'splitlines',
    'clingen_erepo_variant_classifications': 'split',
    'gwas_catalog_associations': 'gwas',
    'clinvar_vcv_xml_weekly_2026_1018': 'parquet',     # canonical JSON of one row of a stream-derived Parquet table
}


def source_releases(index_dir):
    p = os.path.join(index_dir, 'source_release.parquet')
    return {r['release_id']: r for r in pq.read_table(p).to_pylist()} if os.path.exists(p) else {}


def linked_releases(index_dir):
    """{table: set of release ids used}."""
    out = {}
    for name, sch in schema.TABLES.items():
        if 'release_id' not in sch.names or name == 'source_release':
            continue
        p = os.path.join(index_dir, name)
        if not os.path.isdir(p):
            continue
        ids = set()
        for b in pds.dataset(p, format='parquet', ignore_prefixes=[]).to_batches(columns=['release_id']):
            ids.update(x for x in pc.unique(b['release_id']).to_pylist() if x is not None)
        out[name] = ids
    return out


def check_releases(releases, manifest, targets, data_root, rehash=False, rehash_limit=3 << 30):
    issues = []
    for rid, r in sorted(releases.items()):
        m = manifest.get(r['target_id'])
        t = targets.get(r['target_id'], {})
        if m is None:
            issues.append(dict(release_id=rid, issue='no_manifest_record_for_target'))
            continue
        if not m.get('sha256'):
            issues.append(dict(release_id=rid, issue='manifest_record_without_sha256'))
        if r.get('sha256') != m.get('sha256'):
            issues.append(dict(release_id=rid, issue='release_digest_differs_from_manifest'))
        for fld in ('url', 'licence'):
            if not (r.get(fld) or m.get(fld) or t.get(fld)):
                issues.append(dict(release_id=rid, issue=f'no_{fld}'))
        path = os.path.join(data_root, m.get('dest') or t.get('dest') or '')
        if not os.path.exists(path):
            issues.append(dict(release_id=rid, issue='acquired_file_missing', path=m.get('dest')))
        elif os.path.isfile(path):
            if m.get('bytes') is not None and os.path.getsize(path) != m['bytes']:
                issues.append(dict(release_id=rid, issue='acquired_file_size_differs_from_manifest', manifest=m['bytes'], disk=os.path.getsize(path)))
            elif rehash and os.path.getsize(path) <= rehash_limit:
                h = hashlib.sha256()
                with open(path, 'rb') as f:
                    for b in iter(lambda: f.read(1 << 22), b''):
                        h.update(b)
                if h.hexdigest() != m.get('sha256'):
                    issues.append(dict(release_id=rid, issue='acquired_file_digest_differs_from_manifest'))
    return issues


def sample_rows(index_dir, release_id, k, seed):
    """Deterministic reservoir sample of (evidence_id, locator, raw_sha256, accession) rows of one release, read in index order."""
    p = os.path.join(index_dir, 'evidence')
    if not os.path.isdir(p):
        return [], 0
    rng, keep, n = random.Random(f'{seed}:{release_id}'), [], 0
    ds = pds.dataset(p, format='parquet', ignore_prefixes=[])
    for b in ds.to_batches(columns=['evidence_id', 'locator', 'raw_sha256', 'accession'], filter=pc.field('release_id') == release_id):
        for row in b.to_pylist():
            n += 1
            if len(keep) < k:
                keep.append(row)
            else:
                j = rng.randrange(n)
                if j < k:
                    keep[j] = row
    return sorted(keep, key=lambda r: r['evidence_id']), n


def verify_index(index_dir, sample=500, seed=0, data_root=None, paths=None, rehash=False, manifest_path=None, targets_path=None):
    """Run the three checks. paths: {target_id: file or directory} overrides where a resolver finds the acquired source (tests, relocated data)."""
    data_root = data_root or DATA
    paths = paths or {}
    manifest = {}
    mp = manifest_path or os.path.join(REPO, 'acquisition', 'manifest.jsonl')
    if os.path.exists(mp):
        with open(mp) as f:
            manifest = {r['id']: r for r in (json.loads(l) for l in f if l.strip())}
    tp = targets_path or os.path.join(REPO, 'acquisition', 'targets.json')
    targets = {}
    if os.path.exists(tp):
        with open(tp) as f:
            targets = {t['id']: t for t in json.load(f)}
    releases = source_releases(index_dir)
    used = linked_releases(index_dir)
    evidence_dir = os.path.join(index_dir, 'evidence')
    if os.path.isdir(evidence_dir):
        ids = set()
        for b in pds.dataset(evidence_dir, format='parquet', ignore_prefixes=[]).to_batches(columns=['release_id']):
            ids.update(pc.unique(b['release_id']).to_pylist())
        used['evidence'] = ids
    orphans = {t: sorted(i for i in ids if i not in releases) for t, ids in used.items() if any(i not in releases for i in ids)}
    report = dict(index=os.path.abspath(index_dir), releases=len(releases), tables_with_releases=sorted(used), orphan_release_ids=orphans,
                  release_issues=check_releases(releases, manifest, targets, data_root, rehash), digests={}, no_resolver=[])
    for rid, rel in sorted(releases.items()):
        tid = rel['target_id']
        kind = RESOLVERS.get(tid)
        rows, total = sample_rows(index_dir, rid, sample, seed)
        if not total:
            continue                                            # a release that only feeds entity tables has no evidence rows to re-open
        if kind is None:
            report['no_resolver'].append(dict(release_id=rid, target_id=tid, evidence_rows=total))
            continue
        base = paths.get(tid) or os.path.join(data_root, (targets.get(tid) or manifest.get(tid) or {}).get('dest', ''))
        got = resolve_parquet_rows(base, [r['locator'] for r in rows]) if kind == 'parquet' else (resolve_lines(kind, base, [r['locator'] for r in rows]) if os.path.exists(base) else {})
        verified, bad, missing = 0, [], []
        for r in rows:
            text = got.get(r['locator'])
            if text is None:
                missing.append(dict(evidence_id=r['evidence_id'], locator=r['locator']))
            elif sha256_text(text) == r['raw_sha256']:
                verified += 1
            else:
                bad.append(dict(evidence_id=r['evidence_id'], locator=r['locator']))
        report['digests'][rid] = dict(target_id=tid, resolver=kind, evidence_rows=total, sampled=len(rows), verified=verified, mismatched=bad[:20], n_mismatched=len(bad),
                                      unresolved=missing[:20], n_unresolved=len(missing))
    report['ok'] = (not orphans and not report['release_issues'] and all(d['n_mismatched'] == 0 and d['n_unresolved'] == 0 for d in report['digests'].values()))
    return report
