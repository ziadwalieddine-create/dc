"""Extract, from the MaveDB bulk dump v4, the members that the compatibility tests and the superseded-set builder read.

  main.json, README.md, LICENSE.txt and every per-set annotations CSV   -> <out>/dump_v4_extract/<root>/
  scores/counts CSVs of the 52 superseded sets (mavedb_superseded/hidden_from_search_urns.json) and of the fixed
  comparison sample (dump_equivalence_sample.json)                        -> <out>/dump_v4_sets/<root>/csv/

The archive is untrusted input: members are selected by name, only regular files are read, and every output path is built
from the member's base name (nothing from the archive's own paths is used, and extractall is never called). The 9.5 GB of
mapped-variant JSON and 3.4 GB of VA-Spec files are not extracted. Existing files of the right size are left alone.

  python3 extract_dump.py [--out data/external/mavedb] [--dump PATH]
"""
import argparse
import json
import os
import re
import shutil
import sys
import tarfile

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, '..'))
ap = argparse.ArgumentParser()
ap.add_argument('--out', default=os.path.join(ROOT, 'data', 'external', 'mavedb'))
ap.add_argument('--dump', default=os.path.join(ROOT, 'data', 'external', 'mavedb', 'mavedb-dump.2026062418131.tar.gz'))
a = ap.parse_args()

wanted = set(json.load(open(os.path.join(HERE, 'mavedb_superseded', 'hidden_from_search_urns.json')))['urns'])
wanted |= set(json.load(open(os.path.join(HERE, 'dump_equivalence_sample.json')))['urns'])
stems = {u.replace(':', '-') for u in wanted}
top = re.compile(r'^[^/]+/(main\.json|README\.md|LICENSE\.txt)$')
ann = re.compile(r'^[^/]+/csv/(urn-mavedb-[0-9a-z-]+\.annotations\.csv)$')
tab = re.compile(r'^[^/]+/csv/((urn-mavedb-[0-9a-z-]+)\.(scores|counts)\.csv)$')
n = dict(top=0, annotations=0, set_tables=0, skipped_existing=0)
root_name = None
with tarfile.open(a.dump, 'r|gz') as tf:          # one streaming pass
    for m in tf:
        if not m.isreg() or os.path.basename(m.name).startswith('._'):
            continue
        dest = None
        if (g := top.match(m.name)):
            dest, key = os.path.join(a.out, 'dump_v4_extract', m.name.split('/')[0], g.group(1)), 'top'
        elif (g := ann.match(m.name)):
            dest, key = os.path.join(a.out, 'dump_v4_extract', m.name.split('/')[0], 'csv', g.group(1)), 'annotations'
        elif (g := tab.match(m.name)) and g.group(2) in stems:
            dest, key = os.path.join(a.out, 'dump_v4_sets', m.name.split('/')[0], 'csv', g.group(1)), 'set_tables'
        if dest is None:
            continue
        if os.path.exists(dest) and os.path.getsize(dest) == m.size:
            n['skipped_existing'] += 1
            continue
        os.makedirs(os.path.dirname(dest), exist_ok=True)
        with tf.extractfile(m) as fi, open(dest + '.part', 'wb') as fo:
            shutil.copyfileobj(fi, fo)
        os.replace(dest + '.part', dest)
        n[key] += 1
print(json.dumps(n))
