"""Record the INFO/FILTER definitions of the gnomAD v4.1 exome and genome sites VCFs, read from the VCF headers by remote tabix.

The headers are bgzip-compressed at the start of each chromosome file; pysam/htslib reads only the first blocks. Output:
evidence/gnomad/vcf_header_info_v4_1.json with, per callset, the source URL, size, ETag and last-modified of the file read,
the number of INFO fields, the ancestry groups, the FILTER definitions and the descriptions of the INFO fields that
acquisition/fetch_gnomad_regional.py keeps. Needs pysam (acquisition/requirements.txt).

  .venv/bin/python catalogue/derive/gnomad_vcf_header.py
"""
import datetime
import json
import os
import re
import urllib.request

import pysam

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, '..', 'evidence', 'gnomad', 'vcf_header_info_v4_1.json')
KEEP = ['AC', 'AN', 'AF', 'nhomalt', 'grpmax', 'AC_grpmax', 'AN_grpmax', 'AF_grpmax', 'fafmax_faf95_max', 'fafmax_faf95_max_gen_anc', 'cadd_phred', 'revel_max', 'spliceai_ds_max',
        'pangolin_largest_ds', 'phylop', 'sift_max', 'polyphen_max', 'variant_type', 'allele_type', 'lcr', 'segdup', 'fail_interval_qc', 'VRS_Allele_IDs']
out = dict(fetched_utc=datetime.datetime.now(datetime.timezone.utc).strftime('%Y-%m-%dT%H:%M:%SZ'), tool=f'pysam {pysam.__version__} (htslib {pysam.__samtools_version__})', callsets={})
for kind in ('exomes', 'genomes'):
    url = f'https://gnomad-public-us-east-1.s3.amazonaws.com/release/4.1/vcf/{kind}/gnomad.{kind}.v4.1.sites.chr21.vcf.bgz'
    req = urllib.request.Request(url, method='HEAD', headers={'User-Agent': 'dc-catalogue-probe/1.0'})
    h = {k.lower(): v for k, v in urllib.request.urlopen(req, timeout=60).headers.items()}
    tbx = pysam.TabixFile(url)
    info, filt = {}, {}
    for line in tbx.header:
        m = re.match(r'##INFO=<ID=([^,]+),Number=([^,]+),Type=([^,]+),Description="(.*?)"', line)
        if m:
            info[m.group(1)] = dict(number=m.group(2), type=m.group(3), description=m.group(4))
        m = re.match(r'##FILTER=<ID=([^,]+),Description="(.*?)"', line)
        if m:
            filt[m.group(1)] = m.group(2)
    groups = sorted({k[3:] for k in info if re.fullmatch(r'AC_[a-z]+', k)} - {'grpmax', 'raw'})
    out['callsets'][kind] = dict(file_read=url, size=h.get('content-length'), etag=h.get('etag'), last_modified=h.get('last-modified'), info_fields=len(info), ancestry_groups=groups,
                                 filters=filt, kept_fields={k: info[k] for k in KEEP if k in info}, kept_fields_missing=[k for k in KEEP if k not in info])
json.dump(out, open(OUT, 'w'), indent=1)
for kind, d in out['callsets'].items():
    print(kind, d['info_fields'], 'INFO fields; groups', d['ancestry_groups'], '; kept fields missing:', d['kept_fields_missing'])
