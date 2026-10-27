"""Run adapters and finalise an index directory: source_release, entity_variant (deduplicated from evidence), manifest with digests."""
import glob
import hashlib
import inspect
import json
import os
import shutil

import pyarrow.compute as pc
import pyarrow.parquet as pq

from .. import SCHEMA_VERSION, schema
from ..normalise import parse_key, classify
from .common import BuildContext, utc_now, REPO


def file_sha256(path, chunk=1 << 22):
    h = hashlib.sha256()
    with open(path, 'rb') as f:
        for b in iter(lambda: f.read(chunk), b''):
            h.update(b)
    return h.hexdigest()


def code_digest(adapter):
    src = inspect.getsourcefile(type(adapter))
    with open(src, 'rb') as f:
        return hashlib.sha256(f.read()).hexdigest()[:16]


def build_code():
    """Digest of the code that determines the content of an index: every module of vindex.build plus the normaliser, the reference reader and the schema. The digest of an adapter's
    own source file is not enough: a helper they share (the star table, the gene map) changes the rows without changing any adapter file."""
    pkg = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    files = sorted(glob.glob(os.path.join(pkg, 'build', '*.py'))) + [os.path.join(pkg, f) for f in ('__init__.py', 'normalise.py', 'refgenome.py', 'schema.py')]
    h = hashlib.sha256()
    for p in files:
        h.update(os.path.relpath(p, pkg).encode() + b'\0')
        with open(p, 'rb') as f:
            h.update(hashlib.sha256(f.read()).digest())
    return dict(sha256_16=h.hexdigest()[:16], files=len(files))


def natural_chrom(c):
    return (0, int(c), '') if c.isdigit() else (1, 0, c)


def build_entity_variant(out):
    """Unique normalised variant keys over all mapped variant-level evidence and all haplotype definitions -> entity_variant/<partition>.parquet, one partition at a time
    (the keys of a whole index run to ten million and more; a chromosome's worth fits in memory, the whole set does not comfortably). The partition of a key is the
    one its evidence rows are filed under (schema.partition_of): its chromosome on the primary assembly, 'NA' for alternate loci and unplaced scaffolds."""
    hap = {}
    for p in glob.glob(os.path.join(out, 'haplotype_definition', '*.parquet')):
        for k in pc.unique(pq.read_table(p, columns=['variant_key'])['variant_key']).to_pylist():
            if k:                                                   # unresolved definition rows have no key
                hap.setdefault(schema.partition_of(k.split(':', 1)[0]), set()).add(k)
    ev_dir = os.path.join(out, 'evidence', 'variant')
    parts = set(hap) | (set(os.listdir(ev_dir)) if os.path.isdir(ev_dir) else set())
    target = os.path.join(out, 'entity_variant')
    os.makedirs(target, exist_ok=True)
    total = 0
    for c in sorted(parts, key=natural_chrom):
        keys = set(hap.get(c, ()))
        for p in sorted(glob.glob(os.path.join(ev_dir, c, '*.parquet'))):
            t = pq.read_table(p, columns=['subject_id', 'mapping_status'])
            keys.update(pc.unique(t.filter(pc.equal(t['mapping_status'], 'mapped'))['subject_id']).to_pylist())
        keys.discard(None)
        rows = []
        for k in sorted(keys, key=lambda k: (natural_chrom(k.split(':')[0]), int(k.split(':')[1]), k)):
            ch, pos, r, a = parse_key(k)
            rows.append(dict(variant_key=k, chrom=ch, pos=pos, ref=r, alt=a, vtype=classify(r, a)))
        if rows:
            pq.write_table(schema.rows_to_table('entity_variant', rows), os.path.join(target, f'{c}.parquet'), compression='zstd')
            total += len(rows)
    if not total:                                                   # keep an empty but valid table so readers see it
        pq.write_table(schema.empty_table('entity_variant'), os.path.join(target, 'empty.parquet'), compression='zstd')
    return total


def table_inventory(out):
    inv = {}
    for name in schema.TABLES:
        files = sorted(glob.glob(os.path.join(out, name, '**', '*.parquet'), recursive=True))
        if name == 'source_release' and os.path.exists(os.path.join(out, 'source_release.parquet')):
            files = [os.path.join(out, 'source_release.parquet')]
        if not files:
            continue
        rows = sum(pq.ParquetFile(f).metadata.num_rows for f in files)
        inv[name] = dict(files=len(files), rows=rows, sha256={os.path.relpath(f, out): file_sha256(f) for f in files})
    return inv


def catalogue_sources():
    """Source ids the catalogue knows: records in catalogue/sources and the 'source' of every acquisition target."""
    ids = {os.path.splitext(os.path.basename(p))[0] for p in glob.glob(os.path.join(REPO, 'catalogue', 'sources', '*.json'))}
    tp = os.path.join(REPO, 'acquisition', 'targets.json')
    if os.path.exists(tp):
        with open(tp) as f:
            ids |= {t['source'] for t in json.load(f)}
    return ids


def check_sources(out):
    """Every evidence, study, event and assay row must name a catalogue source; a typo would otherwise split one source into two in every summary."""
    used = {}
    for p in glob.glob(os.path.join(out, 'evidence', '*', '*', '*.parquet')):
        used.setdefault(os.path.splitext(os.path.basename(p))[0], 'evidence')
    for t in ('entity_study', 'entity_event', 'entity_assay'):
        for p in glob.glob(os.path.join(out, t, '*.parquet')):
            for sid in pc.unique(pq.read_table(p, columns=['source_id'])['source_id']).to_pylist():
                used.setdefault(sid, t)
    known = catalogue_sources()
    return [f'source id {sid!r} (first seen in {where}) is not a catalogue source' for sid, where in sorted(used.items()) if sid not in known]


def build(out, adapters, reference, limit=None, clean=True):
    if clean and os.path.isdir(out):
        shutil.rmtree(out)
    os.makedirs(out, exist_ok=True)
    code = build_code()                                              # read before the first adapter runs: the digest is of the code that built the index, not of the files as they are at the end
    ctx = BuildContext(out, reference)
    ran, done = [], set()
    for a in adapters:
        missing = [x for x in getattr(a, 'requires', ()) if x not in done]
        if missing:
            raise ValueError(f'{type(a).__name__} needs the output of {missing} first; list those adapters before it (a wrong order would silently leave its records unmapped)')
        done.add(a.source_id)
        if hasattr(a, 'prepare'):
            ctx.adapter = a
            a.prepare(ctx)
    ctx.adapter = None
    done, finished = set(), []
    for a in adapters:
        done.add(a.source_id)
        stats = a.run(ctx, limit=limit) if limit else a.run(ctx)
        ran.append(dict(adapter=type(a).__name__, source_id=a.source_id, code_sha256_16=code_digest(a), stats=dict(sorted(ctx.st(type(a).__name__).c.items()))))
        finished.append(a)
        ctx.release_finished(finished)
    ctx.adapter = None
    ctx.close()
    ctx.problems.extend(check_sources(out))
    pq.write_table(schema.rows_to_table('source_release', list(ctx.releases.values())), os.path.join(out, 'source_release.parquet'), compression='zstd')
    n_var = build_entity_variant(out)
    inv = table_inventory(out)
    man = dict(schema_version=SCHEMA_VERSION, built_utc=utc_now(), reference=dict(path=os.path.relpath(reference.path, REPO) if reference.path.startswith(REPO) else reference.path,
               name=reference.name), limit=limit, build_code=code, adapters=ran, releases=sorted(ctx.releases), problems=ctx.problems, entity_variants=n_var,
               tables={k: dict(files=v['files'], rows=v['rows']) for k, v in inv.items()})
    with open(os.path.join(out, 'index_manifest.json'), 'w') as f:
        json.dump(man, f, indent=1, sort_keys=True)
    digest = {k: v['sha256'] for k, v in inv.items()}
    with open(os.path.join(out, 'table_digests.json'), 'w') as f:
        json.dump(digest, f, indent=1, sort_keys=True)
    return man
