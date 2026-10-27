#!/usr/bin/env python3
"""Compare two sample runs of the matcher: the same input and index must give the same tables whatever the hash seed or the moment.

    python3 variant_index/compare_runs.py RUN_A RUN_B [--out compare.json]

Every result table (*.parquet) of the two directories is compared by SHA-256 of the file, and run_manifest.json field by field except the ones that name the moment or the place of the
run (finished_utc, seconds, the output path). Exit status 0 only when the tables are byte-identical and the manifests agree.
"""
import argparse
import hashlib
import json
import os
import sys

VOLATILE = ('finished_utc', 'seconds')


def sha256_file(path, chunk=1 << 22):
    h = hashlib.sha256()
    with open(path, 'rb') as f:
        for block in iter(lambda: f.read(chunk), b''):
            h.update(block)
    return h.hexdigest()


def tables(run):
    return {f: sha256_file(os.path.join(run, f)) for f in sorted(os.listdir(run)) if f.endswith('.parquet')}


def manifest_differences(a, b):
    """Top-level fields of two run manifests that differ, leaving out the moment and duration of the run."""
    with open(os.path.join(a, 'run_manifest.json')) as f:
        ma = json.load(f)
    with open(os.path.join(b, 'run_manifest.json')) as f:
        mb = json.load(f)
    return sorted(k for k in set(ma) | set(mb) if k not in VOLATILE and ma.get(k) != mb.get(k))


def compare(a, b):
    ta, tb = tables(a), tables(b)
    differ = sorted(f for f in set(ta) & set(tb) if ta[f] != tb[f])
    only_a, only_b = sorted(set(ta) - set(tb)), sorted(set(tb) - set(ta))
    manifest = manifest_differences(a, b)
    return dict(run_a=os.path.abspath(a), run_b=os.path.abspath(b), tables=len(ta), tables_differ=differ, only_in_a=only_a, only_in_b=only_b, manifest_fields_differ=manifest,
                identical=not (differ or only_a or only_b or manifest) and bool(ta))


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('run_a')
    ap.add_argument('run_b')
    ap.add_argument('--out')
    a = ap.parse_args(argv)
    res = compare(a.run_a, a.run_b)
    if a.out:
        with open(a.out, 'w') as f:
            json.dump(res, f, indent=1, sort_keys=True)
    print(json.dumps({k: v for k, v in res.items() if k not in ('run_a', 'run_b')}, indent=1))
    return 0 if res['identical'] else 1


if __name__ == '__main__':
    sys.exit(main())
