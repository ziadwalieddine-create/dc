"""Count the acquired Open Targets credible sets and studies by type and project.

  .venv/bin/python catalogue/derive/ot_credible_set_counts.py        (needs pyarrow)

Writes catalogue/evidence/opentargets/credible_set_counts.json: rows of the credible_set dataset by studyType, rows of the study dataset by projectId and
studyType, and the columns of the credible_set table, so that statements about how much molecular-QTL and GWAS evidence the library already holds are backed
by a recomputable count.
"""
import collections
import datetime
import glob
import json
import os

import pyarrow.parquet as pq

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, '..', '..'))
base = os.path.join(ROOT, 'data', 'external', 'opentargets', '26.09')
cs_files = sorted(glob.glob(os.path.join(base, 'credible_set', '*.parquet')))
st_files = sorted(glob.glob(os.path.join(base, 'study', '*.parquet')))
cs_type = collections.Counter()
methods = collections.Counter()
for f in cs_files:
    tb = pq.read_table(f, columns=['studyType', 'finemappingMethod'])
    cs_type.update(tb.column(0).to_pylist())
    methods.update(tb.column(1).to_pylist())
proj = collections.Counter()
proj_type = collections.defaultdict(collections.Counter)
for f in st_files:
    tb = pq.read_table(f, columns=['projectId', 'studyType'])
    for a, b in zip(tb.column(0).to_pylist(), tb.column(1).to_pylist()):
        proj[a] += 1
        proj_type[a][b] += 1
cols = pq.read_schema(cs_files[0]).names
out = dict(dataset='opentargets 26.09 credible_set and study', generated_utc=datetime.datetime.now(datetime.timezone.utc).strftime('%Y-%m-%dT%H:%M:%SZ'),
           credible_set_files=len(cs_files), credible_set_rows_by_study_type=dict(cs_type), credible_set_rows_total=sum(cs_type.values()), fine_mapping_methods={str(k): v for k, v in methods.items()}, credible_set_columns=cols,
           study_files=len(st_files), studies_by_project=dict(proj.most_common(40)), study_projects=len(proj), studies_total=sum(proj.values()), study_types_of_top_projects={k: dict(v) for k, v in list(proj_type.items())[:0]},
           text=json.dumps(dict(credible_set_rows_by_study_type=dict(cs_type), fine_mapping_methods={str(k): v for k, v in methods.items()}, top_projects=dict(proj.most_common(12)))))
json.dump(out, open(os.path.join(ROOT, 'catalogue', 'evidence', 'opentargets', 'credible_set_counts.json'), 'w'), indent=1)
print(out['credible_set_rows_by_study_type'], out['fine_mapping_methods'], out['studies_total'], out['study_projects'])
