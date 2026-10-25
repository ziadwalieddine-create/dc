"""Count registry tools relevant to human variant interpretation, from the OpenEBench monitor snapshot.

bio.tools' own API is behind a JavaScript bot challenge for non-browser clients, so this uses OpenEBench
(https://openebench.bsc.es/monitor/tool/, CC BY 4.0), which aggregates bio.tools, Bioconda, Galaxy and
OpenEBench-native entries. Counts use EDAM annotations carried in each record's semantics.{topics,operations};
EDAM labels were resolved via the EBI Ontology Lookup Service (labels recorded in the output).

  python3 openebench_counts.py [--input snapshot.json]   # downloads (about 250 MB) if --input is not given
Writes ../evidence/biotools/openebench_edam_counts.json
"""
import argparse
import collections
import datetime
import hashlib
import json
import os
import tempfile
import urllib.request

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, '..', 'evidence', 'biotools', 'openebench_edam_counts.json')
URL = 'https://openebench.bsc.es/monitor/tool/'
TOPICS = {'topic_0199': 'Genetic variation', 'topic_3574': 'Human genetics', 'topic_3517': 'GWAS study',
          'topic_3056': 'Population genetics', 'topic_3175': 'Structural variation', 'topic_0622': 'Genomics'}
OPS = {'operation_3227': 'Variant calling', 'operation_3225': 'Variant classification',
       'operation_3228': 'Structural variation detection', 'operation_0331': 'Variant effect prediction',
       'operation_3675': 'Variant filtering', 'operation_3197': 'Genetic variation analysis'}

ap = argparse.ArgumentParser()
ap.add_argument('--input')
a = ap.parse_args()
path = a.input
if not path:
    path = os.path.join(tempfile.gettempdir(), 'openebench_tools_snapshot.json')
    req = urllib.request.Request(URL, headers={'Accept': 'application/json', 'User-Agent': 'dc-catalogue-probe/1.0'})
    with urllib.request.urlopen(req, timeout=300) as r, open(path, 'wb') as f:
        while True:
            c = r.read(1 << 20)
            if not c:
                break
            f.write(c)
h = hashlib.sha256()
with open(path, 'rb') as f:
    for c in iter(lambda: f.read(1 << 20), b''):
        h.update(c)
data = json.load(open(path))
tools = [x for x in data if isinstance(x, dict) and '@id' in x and 'name' in x]
ns = collections.Counter(x.get('@nmsp') for x in tools)
dep = collections.Counter((x.get('@nmsp'), x.get('deprecated')) for x in tools)


def n(items, key, ident):
    iri = 'http://edamontology.org/' + ident
    return {'all': sum(1 for x in items if iri in ((x.get('semantics') or {}).get(key) or [])),
            'biotools_namespace': sum(1 for x in items if x.get('@nmsp') == 'biotools' and iri in ((x.get('semantics') or {}).get(key) or []))}


res = dict(
    source=URL, fetched_utc=datetime.datetime.now(datetime.timezone.utc).strftime('%Y-%m-%dT%H:%M:%SZ'),
    snapshot_bytes=os.path.getsize(path), snapshot_sha256=h.hexdigest(), licence='CC BY 4.0 (per-record @license)',
    records_total=len(tools), by_namespace=dict(ns),
    deprecated_flag_by_namespace={f'{k[0]}:{k[1]}': v for k, v in dep.items()},
    deprecated_flag_note='In this snapshot every bio.tools-namespace record has deprecated=true and no other namespace sets it, so the flag marks the namespace rather than individual tools; it is not used as a filter.',
    latest_record_timestamp=max((x.get('@timestamp') or '' for x in tools)),
    note='Counts are records carrying the EDAM term (only bio.tools-namespace records carry EDAM semantics); a tool can carry several terms, so rows overlap.',
    topics={f'{k} ({v})': n(tools, 'topics', k) for k, v in TOPICS.items()},
    operations={f'{k} ({v})': n(tools, 'operations', k) for k, v in OPS.items()})
json.dump(res, open(OUT, 'w'), indent=2)
print(json.dumps({k: res[k] for k in ('records_total', 'by_namespace', 'deprecated_flag_by_namespace', 'topics', 'operations')}, indent=1))
