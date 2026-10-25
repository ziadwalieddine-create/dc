"""Measure the size of every Open Targets output dataset from the EBI FTP directory indexes (read-only listing requests).

Input : evidence/opentargets/datasets_inventory.json (ftp_output_directories)
Output: evidence/opentargets/dataset_sizes.json

The Apache index prints sizes rounded to two significant figures (K/M/G), so totals are approximate (about +-5%); exact sizes
and checksums come from Content-Length at download time. Each record keeps the listing URL and the SHA-256 of the listing
HTML so the numbers can be audited.

  python3 opentargets_sizes.py
"""
import datetime
import hashlib
import html
import json
import os
import re
import urllib.request

HERE = os.path.dirname(os.path.abspath(__file__))
EVID = os.path.join(HERE, '..', 'evidence', 'opentargets')
inv = json.load(open(os.path.join(EVID, 'datasets_inventory.json')))
release = inv['release']
base = f'https://ftp.ebi.ac.uk/pub/databases/opentargets/platform/{release}/output/'
UNIT = {'K': 1024, 'M': 1024 ** 2, 'G': 1024 ** 3, 'T': 1024 ** 4}
row = re.compile(r'<a href="([^"?/][^"]*)">[^<]*</a>\s*</td><td[^>]*>\s*([0-9-]+ [0-9:]+)\s*</td><td[^>]*>\s*([0-9.]+[KMGT]?|-)\s*</td>')


def size(s):
    return int(float(s[:-1]) * UNIT[s[-1]]) if s[-1] in UNIT else int(float(s))


def listing(url):
    req = urllib.request.Request(url, headers={'User-Agent': 'dc-catalogue-probe/1.0'})
    raw = urllib.request.urlopen(req, timeout=60).read()
    return raw, html.unescape(raw.decode('utf-8', 'replace'))


out = dict(release=release, listing_base=base, fetched_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(timespec='seconds'),
           note='sizes are Apache-index approximations (two significant figures); exact sizes come from Content-Length at download time', datasets={})
for d in inv['ftp_output_directories']:
    raw, text = listing(f'{base}{d}/')
    files, subdirs = [], []
    for name, mtime, sz in row.findall(text):
        if name.endswith('/'):
            subdirs.append(name)
        elif name != '_SUCCESS' and sz != '-':
            files.append((name, size(sz)))
    out['datasets'][d] = dict(url=f'{base}{d}/', files=len(files), approx_bytes=sum(s for _, s in files), subdirectories=subdirs, listing_sha256=hashlib.sha256(raw).hexdigest(),
                              largest_file_bytes=max((s for _, s in files), default=0))
json.dump(out, open(os.path.join(EVID, 'dataset_sizes.json'), 'w'), indent=1)
tot = sum(v['approx_bytes'] for v in out['datasets'].values())
print(len(out['datasets']), 'datasets;', round(tot / 1024 ** 3, 1), 'GiB total (approx)')
for d, v in sorted(out['datasets'].items(), key=lambda kv: -kv[1]['approx_bytes']):
    print(f"{d:45} {v['files']:5} files {v['approx_bytes'] / 1024 ** 2:10.0f} MiB {v['subdirectories'] or ''}")
