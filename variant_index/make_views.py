#!/usr/bin/env python3
"""Priority lists, trait groupings and pathway groupings over the matches of one sample (see vindex/views.py for each view's filter and ordering).

    python3 variant_index/make_views.py --results data/results/SAMPLE [--index data/index] [--include-filtered]

Views select and order rows of small_calls.parquet; no view hides a call from the full table and no view combines evidence axes into one number.
"""
import argparse
import json
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

from vindex import views                                           # noqa: E402
from vindex.build.common import OUT_DEFAULT                        # noqa: E402
from vindex.index_reader import Index                              # noqa: E402


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--results', required=True)
    ap.add_argument('--index', default=OUT_DEFAULT)
    ap.add_argument('--include-filtered', action='store_true', help='also consider calls that failed the quality policy (their presence is not asserted)')
    a = ap.parse_args(argv)
    m = views.write_views(a.results, Index.open(a.index), ('EXACT_MATCH', 'FILTERED_CALL') if a.include_filtered else ('EXACT_MATCH',))
    print(json.dumps(m, indent=1))
    return 0


if __name__ == '__main__':
    sys.exit(main())
