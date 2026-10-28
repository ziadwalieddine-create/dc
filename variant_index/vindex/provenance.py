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

from . import resolvers, schema
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
    """{locator: canonical JSON of the stream-derived row}. A part is read, the wanted rows are taken and converted, and the part is released before the next one is
    opened: the VCV tables hold millions of rows per table and the converted rows of every part must never be held at once."""
    wanted = {}
    for l in locators:
        m = LOCATOR_ROW.match(l)
        if m:
            wanted.setdefault((m.group('table'), m.group('part')), {})[int(m.group('i'))] = l
    out = {}
    for (table, part), by_row in sorted(wanted.items()):
        p = os.path.join(root, table, part)
        if not os.path.exists(p):
            continue
        t = pq.read_table(p)
        keep = sorted(i for i in by_row if i < t.num_rows)
        for i, row in zip(keep, t.take(keep).to_pylist()):
            out[by_row[i]] = vcv_canon(row)
        del t
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
    # the rules of the other adapters, each in vindex.resolvers (a string is a line rule above, a function is a resolver)
    'clingen_gene_validity': resolvers.clingen_csv_records,
    'clingen_gene_curation_list': 'iter_raw',
    'clingen_region_curation_list': 'iter_raw',
    'clingen_actionability_adult_assertions': 'iter_raw',
    'clingen_actionability_adult_scoring': 'iter_raw',
    'clingen_actionability_pediatric_assertions': 'iter_raw',
    'clingen_actionability_pediatric_scoring': 'iter_raw',
    'cpic_api_allele': resolvers.json_array('allele'),
    'cpic_api_recommendation': resolvers.json_array('rec'),
    'clinpgx_summaryAnnotations': resolvers.zip_members,
    'clinpgx_variantAnnotations': resolvers.zip_members,
    'clinpgx_drugLabels': resolvers.zip_members,
    'clinpgx_guidelineAnnotations': resolvers.zip_members,
    'opentargets_26_09_credible_set': resolvers.ot_credible_set,
    'opentargets_26_09_l2g_prediction': resolvers.ot_l2g,
    'opentargets_26_09_evidence_crispr': resolvers.ot_row_fields(('id', 'targetId', 'diseaseId', 'resourceScore')),
    'opentargets_26_09_evidence_crispr_screen': resolvers.ot_row_fields(('id', 'targetId', 'studyId', 'log2FoldChangeValue', 'resourceScore')),
    'opentargets_26_09_evidence_gene_burden': resolvers.ot_row_fields(('id', 'targetId', 'diseaseId', 'oddsRatio', 'beta', 'pValueMantissa', 'pValueExponent')),
    'opentargets_26_09_target_essentiality': resolvers.ot_target_essentiality,
    'mavedb_scores_genomic_mapped': resolvers.mavedb_native,
    'strchive_hg38_general': 'iter_raw',
    'strchive_hg38_stranger_json': resolvers.strchive_stranger,
    'opentargets_26_09_evidence_impc': resolvers.ot_impc_summary,
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
        parent = targets.get(t.get('derived_from'), {})                # a target extracted from another acquired file is retrieved from that file's URL
        for fld in ('url', 'licence'):
            if not (r.get(fld) or m.get(fld) or t.get(fld) or (fld == 'url' and parent.get('url'))):
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


CHUNK = 100_000        # records re-opened per call of a resolver: the texts of a chunk are held, compared and released before the next chunk is read
MAX_EXAMPLES = 20      # mismatched or unresolved records listed per release (all are counted)


def release_dataset(index_dir):
    p = os.path.join(index_dir, 'evidence')
    return pds.dataset(p, format='parquet', ignore_prefixes=[]) if os.path.isdir(p) else None


def count_rows(index_dir, release_id):
    ds = release_dataset(index_dir)
    return ds.count_rows(filter=pc.field('release_id') == release_id) if ds is not None else 0


def sample_rows(index_dir, release_id, k, seed):
    """Deterministic reservoir sample of (evidence_id, locator, raw_sha256, accession) rows of one release, read in index order."""
    ds = release_dataset(index_dir)
    if ds is None:
        return [], 0
    rng, keep, n = random.Random(f'{seed}:{release_id}'), [], 0
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


def row_chunks(index_dir, release_id, chunk):
    """Every (evidence_id, locator, raw_sha256, accession) row of one release in index order, in lists of at most `chunk` rows (the rows of a release of millions are never held at once)."""
    ds, buf = release_dataset(index_dir), []
    for b in ds.to_batches(columns=['evidence_id', 'locator', 'raw_sha256', 'accession'], filter=pc.field('release_id') == release_id):
        buf += b.to_pylist()
        while len(buf) >= chunk:
            yield buf[:chunk]
            buf = buf[chunk:]
    if buf:
        yield buf


def verify_release(index_dir, rid, tid, base, sample=500, seed=0, chunk=None):
    """Open the source records of one release and compare the SHA-256 of their text with the digest in the row: a deterministic sample of `sample` rows, or every row when the release has no
    more than that. The records are re-opened `chunk` at a time. None: the release has no evidence rows; no_resolver: no rule recovers its records; else the digest result."""
    chunk = chunk or CHUNK
    kind = RESOLVERS.get(tid)
    total = count_rows(index_dir, rid)
    if not total:
        return None
    if kind is None:
        return dict(target_id=tid, no_resolver=True, evidence_rows=total)
    if sample >= total:
        batches = row_chunks(index_dir, rid, chunk)
    else:
        rows, _ = sample_rows(index_dir, rid, sample, seed)
        batches = (rows[i:i + chunk] for i in range(0, len(rows), chunk))
    sampled = verified = n_bad = n_missing = 0
    bad, missing = [], []
    for rows in batches:
        locators = [r['locator'] for r in rows]
        if callable(kind):
            got = kind(base, locators) if os.path.exists(base) else {}
        else:
            got = resolve_parquet_rows(base, locators) if kind == 'parquet' else (resolve_lines(kind, base, locators) if os.path.exists(base) else {})
        sampled += len(rows)
        for r in rows:
            text = got.get(r['locator'])
            if text is None:
                n_missing += 1
                if len(missing) < MAX_EXAMPLES:
                    missing.append(dict(evidence_id=r['evidence_id'], locator=r['locator']))
            elif sha256_text(text) == r['raw_sha256']:
                verified += 1
            else:
                n_bad += 1
                if len(bad) < MAX_EXAMPLES:
                    bad.append(dict(evidence_id=r['evidence_id'], locator=r['locator']))
        del got
    return dict(target_id=tid, resolver=kind if isinstance(kind, str) else kind.__name__, evidence_rows=total, sampled=sampled, verified=verified, mismatched=bad, n_mismatched=n_bad,
                unresolved=missing, n_unresolved=n_missing)


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
    here = os.path.dirname(os.path.abspath(__file__))
    report = dict(index=os.path.abspath(index_dir), releases=len(releases), tables_with_releases=sorted(used), orphan_release_ids=orphans,
                  release_issues=check_releases(releases, manifest, targets, data_root, rehash), digests={}, no_resolver=[],
                  code={f: hashlib.sha256(open(os.path.join(here, f), 'rb').read()).hexdigest() for f in ('provenance.py', 'resolvers.py')}, sample=sample, seed=seed)
    for rid, rel in sorted(releases.items()):
        tid = rel['target_id']
        base = paths.get(tid) or os.path.join(data_root, (targets.get(tid) or manifest.get(tid) or {}).get('dest', ''))
        d = verify_release(index_dir, rid, tid, base, sample, seed)
        if d is None:
            continue                                            # a release that only feeds entity tables has no evidence rows to re-open
        if d.get('no_resolver'):
            report['no_resolver'].append(dict(release_id=rid, target_id=tid, evidence_rows=d['evidence_rows']))
            continue
        report['digests'][rid] = d
    report['ok'] = (not orphans and not report['release_issues'] and all(d['n_mismatched'] == 0 and d['n_unresolved'] == 0 for d in report['digests'].values()))
    return report
