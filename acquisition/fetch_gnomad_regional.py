"""Slice the gnomAD v4.1 exome and genome sites VCFs to the positions of MaveDB-mapped SNVs, by remote tabix.

The gnomAD VCFs are 11-67 GiB per chromosome and are not mirrored. They are bgzip-compressed with tabix indexes on a public
bucket, so a genomic region can be read with HTTP range requests (pysam / htslib) without downloading a file. This script
asks only for the regions that contain a MaveDB-mapped SNV position (positions merged when closer than 150 bp), keeps the
single-nucleotide records at those exact positions, and writes one table per VCF (exomes, genomes).

  Outputs  data/external/gnomad_regional/gnomad_v4_1_<exomes|genomes>_at_mavedb_positions.tsv.gz
           data/external/gnomad_regional/mavedb_snv_regions.bed   (the regions asked for, 0-based half-open, GRCh38)
           acquisition/manifest_files/gnomad_v4_1_<exomes|genomes>_regional.tsv   (per-chromosome source file, size, ETag, last-modified, regions, records)
  Manifest records gnomad_v4_1_exomes_regional and gnomad_v4_1_genomes_regional in acquisition/manifest.jsonl.

Kept per record: filter, AC/AN/AF/nhomalt, group-maximum frequency (grpmax), per-ancestry AC/AN/AF, the in-silico scores gnomAD
carries in INFO (CADD, REVEL, SpliceAI, Pangolin, PhyloP, SIFT, PolyPhen) and a few quality flags. Multi-allelic records
are split and each allele takes its own value for fields defined per allele. Predictor scores keep the licences of their
authors; gnomAD itself is CC0.

  .venv/bin/python acquisition/fetch_gnomad_regional.py [--workers 6] [--chroms 17 X] [--kinds joint]      (requires pysam; see requirements.txt)
"""
import argparse
import collections
import concurrent.futures
import datetime
import gzip
import hashlib
import os
import re
import sys
import time
import urllib.request

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, '..'))
sys.path.insert(0, HERE)
import acquire  # noqa: E402
import mavedb_snvs  # noqa: E402

OUT = os.path.join(ROOT, 'data', 'external', 'gnomad_regional')
CACHE = os.path.join(OUT, 'tbi_cache')
BASE = 'https://gnomad-public-us-east-1.s3.amazonaws.com/release/4.1/vcf/{kind}/gnomad.{kind}.v4.1.sites.chr{chrom}.vcf.bgz'
CHROMS = [str(i) for i in range(1, 23)] + ['X', 'Y']
PER_SITE = ['grpmax', 'fafmax_faf95_max', 'fafmax_faf95_max_gen_anc', 'cadd_phred', 'revel_max', 'spliceai_ds_max', 'pangolin_largest_ds', 'phylop', 'sift_max', 'polyphen_max',
            'variant_type', 'allele_type', 'lcr', 'segdup', 'fail_interval_qc', 'VRS_Allele_IDs']
CORE = ['AC', 'AN', 'AF', 'nhomalt', 'AC_grpmax', 'AN_grpmax', 'AF_grpmax']
# The joint callset (exomes + genomes pooled) has its own field names, carries the frequencies Open Targets and MaveDB quote, and holds no predictor scores.
JOINT_CORE = ['AC_joint', 'AN_joint', 'AF_joint', 'nhomalt_joint', 'AC_grpmax_joint', 'AN_grpmax_joint', 'AF_grpmax_joint',
              'AC_exomes', 'AN_exomes', 'AF_exomes', 'AC_genomes', 'AN_genomes', 'AF_genomes']
JOINT_PER_SITE = ['grpmax_joint', 'fafmax_faf95_max_joint', 'fafmax_faf95_max_gen_anc_joint', 'exomes_filters', 'genomes_filters', 'not_called_in_exomes', 'not_called_in_genomes']
KINDS = ('exomes', 'genomes', 'joint')


def layout(kind, groups):
    """(core + per-group fields, per-site fields) kept for a callset."""
    if kind == 'joint':
        return JOINT_CORE + [f'{m}_joint_{g}' for g in groups for m in ('AC', 'AN', 'AF')], JOINT_PER_SITE
    return CORE + [f'{m}_{g}' for g in groups for m in ('AC', 'AN', 'AF')], PER_SITE


def now():
    return datetime.datetime.now(datetime.timezone.utc).strftime('%Y-%m-%dT%H:%M:%SZ')


def info_dict(info):
    d = {}
    for kv in info.split(';'):
        k, _, v = kv.partition('=')
        d[k] = v if _ else True
    return d


def pick(value, idx, n_alt):
    """Value of an INFO field for allele idx: flags give 'true', a per-allele list (one value per ALT) is indexed, anything else is shared."""
    if value is None:
        return ''
    if value is True:
        return 'true'
    if n_alt > 1:
        parts = value.split(',')
        if len(parts) == n_alt:
            return parts[idx]
    return value


def groups_of(tbx, kind):
    ids = set()
    pat = r'##INFO=<ID=AC_joint_([a-z]+),' if kind == 'joint' else r'##INFO=<ID=AC_([a-z]+),'
    for line in tbx.header:
        m = re.match(pat, line)
        if m:
            ids.add(m.group(1))
    return sorted(ids - {'grpmax', 'raw'})


def head(url):
    req = urllib.request.Request(url, method='HEAD', headers={'User-Agent': acquire.UA})
    with urllib.request.urlopen(req, timeout=60) as r:
        return {k.lower(): v for k, v in r.headers.items()}


def read_region(state, chrom, start, end, tries=5):
    """Lines of one region, accepted only when two consecutive reads agree. A transient network or BGZF error (htslib prints
    '[E::bgzf_read_block] ...') could otherwise end an iteration early and silently drop records; a read that raises or
    differs from its twin is repeated on a fresh handle, and the region fails loudly after `tries` attempts.
    `state` holds the source url and the open pysam.TabixFile under 'tbx'."""
    import pysam
    retries, err = 0, None
    for attempt in range(tries):
        try:
            a = list(state['tbx'].fetch('chr' + chrom, start - 1, end))
            b = list(state['tbx'].fetch('chr' + chrom, start - 1, end))
            if a == b:
                return a, retries
            err = f'two reads of the same region differ ({len(a)} vs {len(b)} records)'
        except (OSError, ValueError) as e:
            err = repr(e)
        retries += 1
        time.sleep(2 ** attempt)
        state['tbx'] = pysam.TabixFile(state['url'])
    raise RuntimeError(f'chr{chrom}:{start}-{end}: {err}')


def slice_chrom(args):
    kind, chrom, regs, wanted = args
    import pysam
    os.makedirs(CACHE, exist_ok=True)
    os.chdir(CACHE)             # htslib saves each remote .tbi index in the working directory; keep them out of the repository
    url = BASE.format(kind=kind, chrom=chrom)
    hdr = head(url)
    state = dict(url=url, tbx=pysam.TabixFile(url))
    groups = groups_of(state['tbx'], kind)
    fields, per_site = layout(kind, groups)
    rows, seen, retries = [], 0, 0
    for (c, s, e) in regs:
        lines, r_ = read_region(state, c, s, e)
        retries += r_
        for line in lines:
            p = line.split('\t', 8)
            pos = int(p[1])
            if pos not in wanted or len(p[3]) != 1:
                continue
            alts = p[4].split(',')
            alts_idx = [(i, a) for i, a in enumerate(alts) if len(a) == 1]
            if not alts_idx:
                continue
            info = info_dict(p[7])
            seen += 1
            for i, a in alts_idx:
                rows.append([chrom, pos, p[3], a, p[6]] + [pick(info.get(f), i, len(alts)) for f in fields] + [pick(info.get(f), i, len(alts)) for f in per_site])
    meta = dict(kind=kind, chrom=chrom, url=url, size=hdr.get('content-length'), etag=hdr.get('etag'), last_modified=hdr.get('last-modified'), regions=len(regs), source_records_at_positions=seen,
                rows=len(rows), read_retries=retries)
    return meta, ['chrom', 'pos', 'ref', 'alt', 'filter'] + fields + per_site, rows


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--workers', type=int, default=6)
    ap.add_argument('--chroms', nargs='*')
    ap.add_argument('--kinds', nargs='*', choices=KINDS, help='callsets to slice (default: all three)')
    a = ap.parse_args()
    import pysam
    os.makedirs(OUT, exist_ok=True)
    keys = mavedb_snvs.snv_alleles()
    by = collections.defaultdict(set)
    for c, p, r, al in keys:
        by[c].add(p)
    regs = mavedb_snvs.regions(by)
    in_gnomad = [r for r in regs if r[0] in CHROMS and (not a.chroms or r[0] in a.chroms)]
    bed = os.path.join(OUT, 'mavedb_snv_regions.bed')
    with open(bed, 'w') as f:
        for c, s, e in [r for r in regs if r[0] in CHROMS]:
            f.write(f'chr{c}\t{s - 1}\t{e}\n')
    region_sha = hashlib.sha256(open(bed, 'rb').read()).hexdigest()
    print(f'{len(keys)} SNV alleles, {sum(len(v) for v in by.values())} positions, {len(regs)} regions ({len(in_gnomad)} on gnomAD chromosomes)', flush=True)
    for kind in (a.kinds or KINDS):
        tid = f'gnomad_v4_1_{kind}_regional'
        per_chr = collections.defaultdict(list)
        for r in in_gnomad:
            per_chr[r[0]].append(r)
        jobs = [(kind, c, per_chr[c], by[c]) for c in sorted(per_chr, key=lambda x: -len(per_chr[x]))]
        metas, all_rows, header = [], [], None
        with concurrent.futures.ProcessPoolExecutor(max_workers=a.workers) as ex:
            for meta, hdr, rows in ex.map(slice_chrom, jobs):
                metas.append(meta)
                header = header or hdr
                if header != hdr:                       # groups differ between chromosomes only if the source changed under us
                    raise SystemExit(f'column layout differs between chromosomes for {kind}')
                all_rows += rows
                print(f"  {kind} chr{meta['chrom']}: {meta['regions']} regions, {meta['source_records_at_positions']} records at wanted positions, {meta['rows']} alleles, {meta['read_retries']} read retries", flush=True)
        all_rows.sort(key=lambda r: (CHROMS.index(r[0]), r[1], r[3]))
        outp = os.path.join(OUT, f'gnomad_v4_1_{kind}_at_mavedb_positions.tsv.gz')
        with gzip.GzipFile(outp, 'wb', mtime=0) as g:
            g.write(('\t'.join(header) + '\n').encode())
            for r in all_rows:
                g.write(('\t'.join(str(x) for x in r) + '\n').encode())
        table = os.path.join(HERE, 'manifest_files', tid + '.tsv')
        os.makedirs(os.path.dirname(table), exist_ok=True)
        with open(table, 'w') as f:
            f.write('chrom\turl\tsize\tetag\tlast_modified\tregions\trecords_at_positions\talleles\tread_retries\n')
            for m in sorted(metas, key=lambda m: CHROMS.index(m['chrom'])):
                f.write('\t'.join(str(m[k]) for k in ('chrom', 'url', 'size', 'etag', 'last_modified', 'regions', 'source_records_at_positions', 'rows', 'read_retries')) + '\n')
        rec = dict(id=tid, source='gnomad', url=BASE.format(kind=kind, chrom='N'), dest=os.path.relpath(outp, os.path.join(ROOT, 'data', 'external')), tier='A', kind='regional_slice',
                   ok=bool(all_rows) and not a.chroms, fetched_utc=now(), bytes=os.path.getsize(outp), rows=len(all_rows), n_regions=len(in_gnomad), regions_sha256=region_sha,
                   sha256=hashlib.sha256(open(outp, 'rb').read()).hexdigest(), md5=None, source_md5=None, md5_verified_against_source=False, size_verified_against_source=False,
                   files_table=os.path.relpath(table, ROOT), tool=f'pysam {pysam.__version__} (htslib {pysam.__samtools_version__}) remote tabix',
                   note='range reads of the public gnomAD v4.1 site VCFs; per-chromosome source ETag and size are in the files table; gnomAD is CC0, embedded predictor scores keep their authors licences')
        if a.chroms:
            rec['note'] += ' [partial run: --chroms given, not recorded as ok]'
        acquire.update_manifest(rec)
        print(f'{tid}: {len(all_rows)} alleles from {len(in_gnomad)} regions, {rec["bytes"] / 1e6:.2f} MB written', flush=True)


if __name__ == '__main__':
    main()
