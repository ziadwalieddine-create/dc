"""Extract the dataset inventory (name, category, one-line description) from the rendered Open Targets downloads page.

Input: evidence/opentargets/downloads_page.txt (produced by probe.py). Output: evidence/opentargets/datasets_inventory.json
"""
import json
import os
import re

HERE = os.path.dirname(os.path.abspath(__file__))
src = os.path.join(HERE, '..', 'evidence', 'opentargets', 'downloads_page.txt')
out = os.path.join(HERE, '..', 'evidence', 'opentargets', 'datasets_inventory.json')
lines = [l.strip() for l in open(src, encoding='utf8').read().split('\n') if l.strip() and l.strip() != '​']
release = next((lines[i + 1] for i, l in enumerate(lines) if l == 'Parquet'), None)
n_all = next((int(re.search(r'\((\d+)\)', l).group(1)) for l in lines if l.startswith('All Datasets (')), None)
start = next(i for i, l in enumerate(lines) if l.startswith('All Datasets ('))
items, i = [], start + 1
while i < len(lines):
    name = lines[i]
    j = i + 1
    desc = []
    while j < len(lines) and not lines[j].startswith('[ ... more ]'):
        desc.append(lines[j])
        j += 1
    if j >= len(lines):
        break
    cat = lines[j + 1] if j + 1 < len(lines) else None
    items.append(dict(name=name, description=' '.join(desc)[:300], category=cat))
    # skip to after 'Access Data'
    k = j + 1
    while k < len(lines) and lines[k] != 'Access Data':
        k += 1
    i = k + 1
ftp = os.path.join(HERE, '..', 'evidence', 'opentargets', 'ftp_26_09_output.txt')
dirs = re.findall(r'(\S+)/\s+\d{4}-\d{2}-\d{2}', open(ftp, encoding='utf8').read()) if os.path.exists(ftp) else []
json.dump(dict(release=release, listed_total=n_all, parsed_from_rendered_page=len(items), ftp_output_directories=dirs,
               ftp_directory_count=len(dirs), note='The rendered downloads page exposes only part of the list; the FTP listing is the complete inventory.',
               datasets_from_rendered_page=items), open(out, 'w'), indent=1)
print(release, n_all, len(items), len(dirs))
import collections
print(collections.Counter(d['category'] for d in items))
