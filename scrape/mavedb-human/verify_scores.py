"""Verify scores/ against the manifest and the metadata files.

Checks: metadata files agree on the set of score sets; the manifest covers exactly the expected
files; every file's sha256, header, row count and row width match the manifest and the score set's
numVariants; accessions start with the score set's own URN (known exception: tmp: accessions).
Exit status is non-zero if anything fails. Usage: python3 verify_scores.py
"""
import csv
import gzip
import hashlib
import json
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, 'scores')
problems = []

det = {}
for fn in ('score_set_details.jsonl', 'score_set_details_accession_targets.jsonl'):
    for line in open(os.path.join(HERE, fn)):
        r = json.loads(line)
        if r['urn'] in det:
            problems.append(('duplicate urn across detail files', r['urn']))
        det[r['urn']] = r

search = {s['urn'] for f in ('score_sets.json', 'score_sets_accession_targets.json')
          for s in json.load(open(os.path.join(HERE, f)))['scoreSets']}
index = {r['urn'] for r in csv.DictReader(open(os.path.join(HERE, 'score_sets_index.csv')))}
if not (search == set(det) == index):
    problems.append(('metadata files disagree', len(search), len(det), len(index)))

manifest = {}
for line in open(os.path.join(HERE, 'scores_manifest.jsonl')):
    r = json.loads(line)
    manifest[(r['urn'], r['kind'])] = r
expected = {(u, 'scores') for u in det} | {(u, 'counts') for u, d in det.items()
                                           if (d.get('datasetColumns') or {}).get('countColumns')}
if set(manifest) != expected:
    problems.append(('manifest != expected files', len(manifest), len(expected)))
if {f for f in os.listdir(OUT)} != {r['file'] for r in manifest.values()}:
    problems.append(('files on disk != manifest',))

rows_total = {'scores': 0, 'counts': 0}
tmp_acc = {}
for (urn, kind), r in sorted(manifest.items()):
    path = os.path.join(OUT, r['file'])
    if not r['ok']:
        problems.append(('not ok', urn, kind))
        continue
    if hashlib.sha256(open(path, 'rb').read()).hexdigest() != r['sha256_gz']:
        problems.append(('sha256', urn, kind))
    with gzip.open(path, 'rt', newline='') as fh:
        rd = csv.reader(fh)
        header = next(rd)
        n = 0
        widths = set()
        for row in rd:
            n += 1
            widths.add(len(row))
            if not row[0].startswith(urn + '#'):
                tmp_acc.setdefault((urn, kind), []).append(row[0])
    if header != r['columns']:
        problems.append(('header', urn, kind))
    if n != r['rows'] or n != det[urn]['numVariants']:
        problems.append(('rows', urn, kind, n, r['rows'], det[urn]['numVariants']))
    if widths and widths != {len(header)}:
        problems.append(('ragged rows', urn, kind))
    rows_total[kind] += n

foreign = {k: len(v) for k, v in tmp_acc.items()}
unexpected = {k: v for k, v in foreign.items() if k != ('urn:mavedb:00001239-a-6', 'scores')}
if unexpected:
    problems.append(('accessions outside own URN', unexpected))
bad_tmp = [a for a in tmp_acc.get(('urn:mavedb:00001239-a-6', 'scores'), []) if not a.startswith('tmp:')]
if bad_tmp:
    problems.append(('non-tmp foreign accessions in 00001239-a-6', bad_tmp[:3]))

print('score sets', len(det), '| files', len(manifest), '| score rows', rows_total['scores'], '| count rows', rows_total['counts'])
print('known tmp: accession rows:', foreign.get(('urn:mavedb:00001239-a-6', 'scores'), 0))
print('PROBLEMS:', len(problems))
for p in problems[:20]:
    print('  ', p)
sys.exit(1 if problems else 0)
