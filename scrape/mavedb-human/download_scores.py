"""Download MaveDB score (and count) tables for the human-target score sets.

Reads score_set_details.jsonl, writes gzipped CSVs under scores/ and a manifest
(scores_manifest.jsonl). Resumable: files already present and recorded in the
manifest are skipped. Usage: python3 download_scores.py [workers]
"""
import csv
import concurrent.futures as cf
import gzip
import hashlib
import json
import os
import sys
import urllib.error
import urllib.request

API = 'https://api.mavedb.org/api/v1/score-sets/'
HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, 'scores')
MANIFEST = os.path.join(HERE, 'scores_manifest.jsonl')
workers = int(sys.argv[1]) if len(sys.argv) > 1 else 6

sets = [json.loads(line) for line in open(os.path.join(HERE, 'score_set_details.jsonl'))]
jobs = []
for s in sets:
    jobs.append((s['urn'], 'scores', s['numVariants']))
    if (s.get('datasetColumns') or {}).get('countColumns'):
        jobs.append((s['urn'], 'counts', s['numVariants']))

done = {}
if os.path.exists(MANIFEST):
    for line in open(MANIFEST):
        r = json.loads(line)
        if r.get('ok') and os.path.exists(os.path.join(OUT, r['file'])):
            done[(r['urn'], r['kind'])] = r
todo = [j for j in jobs if (j[0], j[1]) not in done]
print('jobs', len(jobs), 'already done', len(done), 'todo', len(todo), flush=True)
os.makedirs(OUT, exist_ok=True)


def fetch(job, tries=4):
    urn, kind, expected = job
    name = f"{urn.replace(':', '_')}.{kind}.csv.gz"
    path = os.path.join(OUT, name)
    rec = {'urn': urn, 'kind': kind, 'file': name, 'expected_variants': expected, 'ok': False}
    err = None
    for _ in range(tries):
        try:
            req = urllib.request.Request(API + f'{urn}/{kind}', headers={'Accept': 'text/csv'})
            with urllib.request.urlopen(req, timeout=300) as resp, gzip.GzipFile(path + '.part', 'wb', mtime=0) as gz:
                while True:
                    chunk = resp.read(1 << 20)
                    if not chunk:
                        break
                    gz.write(chunk)
            with gzip.open(path + '.part', 'rt', newline='') as fh:
                rd = csv.reader(fh)
                header = next(rd)
                rows = sum(1 for _ in rd)
            os.replace(path + '.part', path)
            h = hashlib.sha256(open(path, 'rb').read()).hexdigest()
            rec.update(ok=True, rows=rows, columns=header, bytes_gz=os.path.getsize(path), sha256_gz=h,
                       rows_match_expected=(rows == expected))
            return rec
        except Exception as e:
            err = repr(e)
    rec['error'] = err
    return rec


n = bad = 0
with open(MANIFEST, 'a') as mf, cf.ThreadPoolExecutor(workers) as ex:
    for rec in ex.map(fetch, todo):
        mf.write(json.dumps(rec) + '\n')
        mf.flush()
        n += 1
        bad += (not rec['ok'])
        if n % 100 == 0:
            print('done', n, 'of', len(todo), 'failed', bad, flush=True)
print('finished', n, 'failed', bad, flush=True)
