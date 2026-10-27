#!/usr/bin/env python3
"""Build the index twice from the same sources and code and compare the tables byte for byte.

    python3 variant_index/check_idempotence.py [--limit 30000] [--sources a,b] [--workdir data/idempotence] [--keep]

Every build starts from an empty directory and writes deterministic Parquet, so two builds must give identical per-file SHA-256 digests (recomputed from the files),
the same release list and the same adapter code digests. Each build runs in its own process with a different PYTHONHASHSEED, so an adapter that emits rows in set or
dict-of-set order cannot pass. With --limit each adapter stops after that many source records (a representative subset in minutes instead of hours).
Exit status 1 if anything differs; the report is written to <workdir>/idempotence.json.
"""
import argparse
import json
import os
import shutil
import subprocess
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

import build_index                                                 # noqa: E402
from vindex.build.common import REPO                               # noqa: E402
from vindex.build.runner import file_sha256                        # noqa: E402


def file_digests(d, recorded, rehash):
    """{(table, file): sha256}: as recorded at build time, or recomputed from the files on disk (which also catches a file changed after the build)."""
    if not rehash:
        return {(t, f): h for t, fs in recorded.items() for f, h in fs.items()}
    out = {}
    for t, fs in recorded.items():
        for f in fs:
            p = os.path.join(d, f)
            out[(t, f)] = file_sha256(p) if os.path.exists(p) else None
    return out


def compare(a_dir, b_dir, rehash=False):
    """Differences between two built indexes: per-file digests (recomputed from the files when rehash), releases, problems and adapter code digests."""
    def load(d, name):
        with open(os.path.join(d, name)) as f:
            return json.load(f)
    da, db = load(a_dir, 'table_digests.json'), load(b_dir, 'table_digests.json')
    ma, mb = load(a_dir, 'index_manifest.json'), load(b_dir, 'index_manifest.json')
    files_a, files_b = file_digests(a_dir, da, rehash), file_digests(b_dir, db, rehash)
    gone_a = {k for k, h in files_a.items() if h is None}
    gone_b = {k for k, h in files_b.items() if h is None}
    ka, kb = set(files_a) - gone_a, set(files_b) - gone_b
    return dict(files=len(files_a), rehashed=rehash, build_code=ma.get('build_code'),
                only_in_first=sorted(f for t, f in ka - kb), only_in_second=sorted(f for t, f in kb - ka),
                digest_differs=sorted(f for (t, f) in ka & kb if files_a[(t, f)] != files_b[(t, f)]),
                releases_equal=ma['releases'] == mb['releases'], problems_equal=ma['problems'] == mb['problems'],
                adapter_code_equal={x['adapter']: x['code_sha256_16'] for x in ma['adapters']} == {x['adapter']: x['code_sha256_16'] for x in mb['adapters']} and ma.get('build_code') == mb.get('build_code'),
                rows_equal={k: v['rows'] for k, v in ma['tables'].items()} == {k: v['rows'] for k, v in mb['tables'].items()})


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--limit', type=int, default=30000)
    ap.add_argument('--sources', default=','.join(build_index.ORDER))
    ap.add_argument('--workdir', default=os.path.join(REPO, 'data', 'idempotence'))
    ap.add_argument('--keep', action='store_true', help='keep both builds (default: remove them after comparing)')
    ap.add_argument('--fasta', default=build_index.FASTA)
    ap.add_argument('--assembly-report', default=build_index.REPORT)
    a = ap.parse_args(argv)
    os.makedirs(a.workdir, exist_ok=True)
    dirs, seconds = [], []
    for tag, seed in (('first', '1'), ('second', '2')):
        d = os.path.join(a.workdir, tag)
        t0 = time.time()
        # a separate process per build with a different hash seed: an adapter that emits rows in set or dict-of-set order cannot give the same files twice
        r = subprocess.run([sys.executable, os.path.join(HERE, 'build_index.py'), '--out', d, '--limit', str(a.limit), '--sources', a.sources, '--fasta', a.fasta, '--assembly-report', a.assembly_report],
                           env=dict(os.environ, PYTHONHASHSEED=seed), capture_output=True, text=True)
        if r.returncode:
            raise SystemExit(f'build {tag} failed (exit {r.returncode}):\n{r.stderr[-2000:]}')
        dirs.append(d)
        seconds.append(round(time.time() - t0, 1))
    rep = compare(*dirs, rehash=True)
    rep.update(limit=a.limit, sources=a.sources.split(','), seconds=seconds, hash_seeds=['1', '2'])
    rep['identical'] = not (rep['only_in_first'] or rep['only_in_second'] or rep['digest_differs']) and rep['releases_equal'] and rep['problems_equal'] and rep['adapter_code_equal'] and rep['rows_equal']
    with open(os.path.join(a.workdir, 'idempotence.json'), 'w') as f:
        json.dump(rep, f, indent=1, sort_keys=True)
    if not a.keep:
        for d in dirs:
            shutil.rmtree(d, ignore_errors=True)
    print(json.dumps(rep, indent=1, sort_keys=True))
    return 0 if rep['identical'] else 1


if __name__ == '__main__':
    sys.exit(main())
