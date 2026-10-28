#!/usr/bin/env python3
"""Check an index's provenance: release linkage, release records against the acquisition manifest, and raw-record digests recomputed from the acquired files.

    python3 variant_index/verify_index.py [--index data/index] [--sample 50000] [--rehash] [--out provenance_report.json]

Exit status 0 only when nothing is orphaned, every release record matches the manifest and every sampled record re-derives its stored digest. Releases without a
resolver are listed (they are not counted as verified). --rehash also recomputes the SHA-256 of every acquired file up to 3 GB.
"""
import argparse
import json
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

from vindex import provenance                                       # noqa: E402
from vindex.build.common import OUT_DEFAULT                          # noqa: E402


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--index', default=OUT_DEFAULT)
    ap.add_argument('--sample', type=int, default=50000)
    ap.add_argument('--seed', type=int, default=0)
    ap.add_argument('--rehash', action='store_true')
    ap.add_argument('--out', help='write the report here (default: <index>/provenance_report.json)')
    a = ap.parse_args(argv)
    rep = provenance.verify_index(a.index, sample=a.sample, seed=a.seed, rehash=a.rehash)
    with open(a.out or os.path.join(a.index, 'provenance_report.json'), 'w') as f:
        json.dump(rep, f, indent=1, sort_keys=True)
    short = dict(ok=rep['ok'], releases=rep['releases'], orphan_release_ids=rep['orphan_release_ids'], release_issues=len(rep['release_issues']),
                 verified={d['target_id']: f"{d['verified']}/{d['sampled']}" for d in rep['digests'].values()}, no_resolver=[x['target_id'] for x in rep['no_resolver']])
    print(json.dumps(short, indent=1))
    return 0 if rep['ok'] else 1


if __name__ == '__main__':
    sys.exit(main())
