"""Record the INFO and FILTER definitions of the gnomAD v4.1 SV sites VCF, read from the header of the downloaded file.

  python3 catalogue/derive/gnomad_sv_header.py

Output: catalogue/evidence/gnomad/sv_vcf_header_info.json (file name, size, SHA-256 from the acquisition manifest, the INFO and FILTER definitions, the ancestry groups).
The header is the file's own description of its fields, so it supports statements such as which INFO fields a copy-number site carries.
"""
import datetime
import gzip
import json
import os
import re

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, '..', '..'))
PATH = os.path.join(ROOT, 'data', 'external', 'gnomad', 'gnomad.v4.1.sv.sites.vcf.gz')
man = {json.loads(x)['id']: json.loads(x) for x in open(os.path.join(ROOT, 'acquisition', 'manifest.jsonl')) if x.strip()}['gnomad_sv_v4_1_sites_vcf']
info, filt, other = {}, {}, {}
with gzip.open(PATH, 'rt') as f:
    for line in f:
        if not line.startswith('##'):
            break
        m = re.match(r'##INFO=<ID=([^,]+),Number=([^,]+),Type=([^,]+),Description="(.*?)"', line)
        if m:
            info[m.group(1)] = dict(number=m.group(2), type=m.group(3), description=m.group(4))
            continue
        m = re.match(r'##FILTER=<ID=([^,]+),Description="(.*?)"', line)
        if m:
            filt[m.group(1)] = m.group(2)
            continue
        m = re.match(r'##(fileformat|fileDate|reference|source)=(.*)', line)
        if m:
            other[m.group(1)] = m.group(2)
groups = sorted({k[3:] for k in info if re.fullmatch(r'AC_[a-z]+', k)} - {'grpmax', 'raw'})
out = dict(fetched_utc=datetime.datetime.now(datetime.timezone.utc).strftime('%Y-%m-%dT%H:%M:%SZ'), file=os.path.relpath(PATH, ROOT), bytes=man['bytes'], sha256=man['sha256'], header=other,
           n_info_fields=len(info), ancestry_groups=groups, filters=filt, info=info)
json.dump(out, open(os.path.join(HERE, '..', 'evidence', 'gnomad', 'sv_vcf_header_info.json'), 'w'), indent=1)
print(len(info), 'INFO fields;', len(filt), 'filters; groups', groups, '; copy-number fields:', [k for k in info if k.startswith('CN_')])
