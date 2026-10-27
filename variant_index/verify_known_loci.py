#!/usr/bin/env python3
"""Check a built index against facts that were established outside it (vindex/known_loci.py, known_loci.json).

    python3 variant_index/verify_known_loci.py [--index data/index] [--spec variant_index/known_loci.json] [--out <index>/reports/known_loci.json]

Each locus states a public fact (an rsID's position, a classification, the direction of an association, a star allele's defining variant, a curated gene-disease or dosage
assertion, an assay that exists) and the checks that read the index the way a user would. Exit status 0 only when every check passes. A FAIL means the index contradicts the
fact; ABSENT means the index holds nothing about the subject. Neither is settled by editing the spec until the raw source has been read.
"""
import argparse
import json
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

from vindex import known_loci                                        # noqa: E402
from vindex.build.common import OUT_DEFAULT                          # noqa: E402
from vindex.index_reader import Index                                # noqa: E402


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--index', default=OUT_DEFAULT)
    ap.add_argument('--spec', default=os.path.join(HERE, 'known_loci.json'))
    ap.add_argument('--out', help='write the report here (default: <index>/reports/known_loci.json)')
    a = ap.parse_args(argv)
    rep = known_loci.run(Index.open(a.index), known_loci.load_spec(a.spec))
    rep['spec'] = os.path.abspath(a.spec)
    out = a.out or os.path.join(a.index, 'reports', 'known_loci.json')
    os.makedirs(os.path.dirname(os.path.abspath(out)), exist_ok=True)
    with open(out, 'w') as f:
        json.dump(rep, f, indent=1, sort_keys=True)
    for locus in rep['loci']:
        for c in locus['checks']:
            print(f"{c['status']:6} {locus['id']}/{c['id']}: {c['observed']}")
    print(f"{rep['n_loci']} loci: {rep['totals']['PASS']} pass, {rep['totals']['FAIL']} fail, {rep['totals']['ABSENT']} absent")
    return 0 if rep['ok'] else 1


if __name__ == '__main__':
    sys.exit(main())
