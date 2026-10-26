"""Run adapters and finalise an index directory: source_release, entity_variant (deduplicated from evidence), manifest with digests."""
import glob
import hashlib
import inspect
import json
import os
import shutil
import sys

import pyarrow as pa
import pyarrow.compute as pc
import pyarrow.dataset as pds
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


def build_entity_variant(out):
    """Unique normalised variant keys over all variant-level evidence and haplotype definitions -> entity_variant."""
    keys = set()
    for p in glob.glob(os.path.join(out, 'evidence', 'variant', '*', '*.parquet')):
        t = pq.read_table(p, columns=['subject_id', 'mapping_status'])
        t = t.filter(pc.equal(t['mapping_status'], 'mapped'))
        keys.update(pc.unique(t['subject_id']).to_pylist())
    for p in glob.glob(os.path.join(out, 'haplotype_definition', '*.parquet')):
        keys.update(k for k in pc.unique(pq.read_table(p, columns=['variant_key'])['variant_key']).to_pylist() if k)      # unresolved definition rows have no key
    rows = []
    for k in sorted(keys, key=lambda k: (parse_key(k)[0], parse_key(k)[1], k)):
        c, p, r, a = parse_key(k)
        rows.append(dict(variant_key=k, chrom=c, pos=p, ref=r, alt=a, vtype=classify(r, a)))
    path = os.path.join(out, 'entity_variant', 'entity_variant.parquet')
    os.makedirs(os.path.dirname(path), exist_ok=True)
    pq.write_table(schema.rows_to_table('entity_variant', rows), path, compression='zstd')
    return len(rows)


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
    ctx = BuildContext(out, reference)
    ran, done = [], set()
    for a in adapters:
        missing = [x for x in getattr(a, 'requires', ()) if x not in done]
        if missing:
            raise ValueError(f'{type(a).__name__} needs the output of {missing} first; list those adapters before it (a wrong order would silently leave its records unmapped)')
        done.add(a.source_id)
        if hasattr(a, 'prepare'):
            a.prepare(ctx)
    done = set()
    for a in adapters:
        done.add(a.source_id)
        stats = a.run(ctx, limit=limit) if limit else a.run(ctx)
        ran.append(dict(adapter=type(a).__name__, source_id=a.source_id, code_sha256_16=code_digest(a), stats=dict(sorted(ctx.st(type(a).__name__).c.items()))))
    ctx.adapter = None
    ctx.close()
    ctx.problems.extend(check_sources(out))
    pq.write_table(schema.rows_to_table('source_release', list(ctx.releases.values())), os.path.join(out, 'source_release.parquet'), compression='zstd')
    n_var = build_entity_variant(out)
    inv = table_inventory(out)
    man = dict(schema_version=SCHEMA_VERSION, built_utc=utc_now(), reference=dict(path=os.path.relpath(reference.path, REPO) if reference.path.startswith(REPO) else reference.path,
               name=reference.name), limit=limit, adapters=ran, releases=sorted(ctx.releases), problems=ctx.problems, entity_variants=n_var,
               tables={k: dict(files=v['files'], rows=v['rows']) for k, v in inv.items()})
    with open(os.path.join(out, 'index_manifest.json'), 'w') as f:
        json.dump(man, f, indent=1, sort_keys=True)
    digest = {k: v['sha256'] for k, v in inv.items()}
    with open(os.path.join(out, 'table_digests.json'), 'w') as f:
        json.dump(digest, f, indent=1, sort_keys=True)
    return man
