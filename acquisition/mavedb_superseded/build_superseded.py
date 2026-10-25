"""Add the MaveDB score sets that the public search hides because they are superseded versions.

Background: POST /api/v1/score-sets/search (the call behind the MaveDB search page) omits superseded versions, so the
human set in scrape/mavedb-human (collected through that search) does not contain them. The Zenodo bulk dump v4 lists
52 published CC0 score sets absent from the search listing; the API serves all of them and flags them as superseded.
They are distinct, published experimental results, so they are kept as a separate, labelled extension. Nothing under
scrape/ is modified.

Inputs : hidden_from_search_urns.json (this folder), the extracted dump csv/ files under data/external/mavedb/dump_v4_sets/
Outputs: superseded_score_set_details.jsonl (API detail records), superseded_index.csv, superseded_manifest.jsonl,
         and gzipped tables in data/external/mavedb_superseded/ (git-ignored).

  python3 build_superseded.py
"""
import csv
import glob
import gzip
import hashlib
import json
import os
import shutil
import urllib.request

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, '..', '..'))
OUT = os.path.join(ROOT, 'data', 'external', 'mavedb_superseded')
os.makedirs(OUT, exist_ok=True)
urns = json.load(open(os.path.join(HERE, 'hidden_from_search_urns.json')))['urns']
csv_dirs = glob.glob(os.path.join(ROOT, 'data', 'external', 'mavedb', 'dump_v4_sets', '*', 'csv'))
assert csv_dirs, 'extract the dump csv files first (see acquisition/README.md)'


def api(urn):
    req = urllib.request.Request('https://api.mavedb.org/api/v1/score-sets/' + urn, headers={'Accept': 'application/json', 'User-Agent': 'dc-catalogue-acquire/1.0'})
    return json.load(urllib.request.urlopen(req, timeout=60))


def org(t):
    return ((t.get('targetSequence') or {}).get('taxonomy') or {}).get('organismName')


details, manifest, problems = [], [], []
for u in urns:
    d = api(u)
    details.append(d)
    human = any(org(t) == 'Homo sapiens' or t.get('targetAccession') for t in d['targetGenes'])
    for kind in ('scores', 'counts'):
        src = os.path.join(csv_dirs[0], u.replace(':', '-') + f'.{kind}.csv')
        if not os.path.exists(src):
            continue
        dest = os.path.join(OUT, u.replace(':', '_') + f'.{kind}.csv.gz')
        with open(src, 'rb') as fi, gzip.GzipFile(dest, 'wb', mtime=0) as fo:
            shutil.copyfileobj(fi, fo)
        with gzip.open(dest, 'rt', newline='') as f:
            rd = csv.reader(f)
            header = next(rd)
            rows = sum(1 for _ in rd)
        rec = dict(urn=u, kind=kind, file=os.path.basename(dest), rows=rows, columns=header, bytes_gz=os.path.getsize(dest),
                   sha256_gz=hashlib.sha256(open(dest, 'rb').read()).hexdigest(), expected_variants=d['numVariants'],
                   rows_match_expected=rows == d['numVariants'], human_target=human, source='MaveDB bulk dump v4 (Zenodo 20840937, CC0)')
        manifest.append(rec)
        if rows != d['numVariants']:
            problems.append((u, kind, rows, d['numVariants']))

with open(os.path.join(HERE, 'superseded_score_set_details.jsonl'), 'w') as f:
    for d in details:
        f.write(json.dumps(d) + '\n')
with open(os.path.join(HERE, 'superseded_manifest.jsonl'), 'w') as f:
    for m in manifest:
        f.write(json.dumps(m, sort_keys=True) + '\n')
with open(os.path.join(HERE, 'superseded_index.csv'), 'w', newline='') as f:
    w = csv.writer(f)
    w.writerow(['urn', 'title', 'published_date', 'num_variants', 'human_target', 'target_genes', 'target_accessions', 'license', 'superseded_by', 'replaces', 'has_scores_table', 'has_counts_table'])
    for d in details:
        human = any(org(t) == 'Homo sapiens' or t.get('targetAccession') for t in d['targetGenes'])
        sup = d.get('supersedingScoreSet') or {}
        rep = d.get('supersededScoreSet') or {}
        have = {m['kind'] for m in manifest if m['urn'] == d['urn']}
        w.writerow([d['urn'], d['title'], d['publishedDate'], d['numVariants'], human, '; '.join(t['name'] for t in d['targetGenes']),
                    '|'.join(sorted({t['targetAccession']['accession'] for t in d['targetGenes'] if t.get('targetAccession')})),
                    (d.get('license') or {}).get('shortName'), sup.get('urn') if isinstance(sup, dict) else sup, rep.get('urn') if isinstance(rep, dict) else rep,
                    'scores' in have, 'counts' in have])
print(len(details), 'detail records;', len(manifest), 'tables;', 'human/accession-target sets:', sum(1 for m in manifest if m['kind'] == 'scores' and m['human_target']))
print('row-count problems:', problems)
