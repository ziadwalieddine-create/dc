#!/usr/bin/env python3
"""Build the PGS Catalog model registry (data/pgs_registry) from the acquired bulk metadata and the scoring-file manifest.

    python3 variant_index/build_pgs_registry.py [--out data/pgs_registry]

The registry lists every score in the catalogue, whether its GRCh38 scoring file was acquired, and the gates that apply before a score can be interpreted.
It is a separate product from the variant index (see vindex/pgs_registry.py).
"""
import argparse
import json
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

from vindex import pgs_registry                                    # noqa: E402
from vindex.build.common import DATA, REPO                          # noqa: E402


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--out', default=os.path.join(REPO, 'data', 'pgs_registry'))
    ap.add_argument('--metadata', default=os.path.join(DATA, 'pgs_catalog', 'metadata'))
    ap.add_argument('--scoring-manifest', default=os.path.join(REPO, 'acquisition', 'manifest_files', 'pgs_scoring_files_grch38.tsv'))
    a = ap.parse_args(argv)
    print(json.dumps(pgs_registry.build_registry(a.metadata, a.scoring_manifest, a.out), indent=1))
    return 0


if __name__ == '__main__':
    sys.exit(main())
