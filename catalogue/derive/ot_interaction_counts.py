"""Count the rows and distinct gene pairs per source database in the acquired Open Targets interaction dataset.

  .venv/bin/python catalogue/derive/ot_interaction_counts.py        (needs pyarrow)

Writes catalogue/evidence/opentargets/interaction_source_counts.json so that statements about how much of IntAct, STRING, Reactome and Signor reaches the
library through Open Targets are backed by a recomputable count and not by a figure typed into a record.
"""
import collections
import datetime
import glob
import json
import os

import pyarrow.parquet as pq

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, '..', '..'))
files = sorted(glob.glob(os.path.join(ROOT, 'data', 'external', 'opentargets', '26.09', 'interaction', '*.parquet')))
rows, both, pairs = collections.Counter(), collections.Counter(), collections.defaultdict(set)
for fp in files:
    tb = pq.read_table(fp, columns=['sourceDatabase', 'targetA', 'targetB'])
    for s, a, b in zip(*(tb.column(c).to_pylist() for c in ('sourceDatabase', 'targetA', 'targetB'))):
        rows[s] += 1
        if a and b:
            both[s] += 1
            pairs[s].add((a, b) if a <= b else (b, a))
out = dict(dataset='opentargets 26.09 interaction', files=len(files), generated_utc=datetime.datetime.now(datetime.timezone.utc).strftime('%Y-%m-%dT%H:%M:%SZ'),
           rows_per_source=dict(rows), rows_with_both_targets=dict(both), distinct_unordered_gene_pairs={k: len(v) for k, v in pairs.items()},
           text=json.dumps(dict(rows_per_source=dict(rows), distinct_pairs={k: len(v) for k, v in pairs.items()})))
os.makedirs(os.path.join(ROOT, 'catalogue', 'evidence', 'opentargets'), exist_ok=True)
json.dump(out, open(os.path.join(ROOT, 'catalogue', 'evidence', 'opentargets', 'interaction_source_counts.json'), 'w'), indent=1)
print(out['rows_per_source'], out['distinct_unordered_gene_pairs'])
