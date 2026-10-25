"""Read the GRCh38 reference sequence for the regions that contain MaveDB-mapped SNV positions, by remote random access.

The Ensembl release-116 GRCh38 "toplevel" FASTA (808 MB) is published bgzip-compressed with .fai and .gzi indexes in
fasta/homo_sapiens/dna_index/, so a region can be read with HTTP range requests (pysam / htslib) without downloading the file. This
script asks for the same merged regions that fetch_gnomad_regional.py uses, on every sequence the FASTA holds (the toplevel file
includes the patches and alternative scaffolds that the Ensembl MaveDB file names), reads every region until two consecutive reads
agree, and writes one row per region. The .fai is kept locally (index_cache/) because htslib holds a remote index in memory and
saves none; validate.py checks the regions against it. The sequence is what the reference-allele concordance tests in
compat_tests.py compare REF alleles against, which checks that the GRCh38 coordinate files we hold (Ensembl MaveDB file, ClinVar,
AlphaMissense, gnomAD) really share one coordinate system and strand convention.

  Output   data/external/reference_sequence/ensembl_grch38_toplevel_at_mavedb_regions.tsv.gz   (chrom, start, end, sequence; 1-based inclusive)
           data/external/reference_sequence/index_cache/Homo_sapiens.GRCh38.dna.toplevel.fa.gz.fai   (copy of the published index)
           acquisition/manifest_files/ensembl_grch38_reference_regional.tsv                    (source file, size, ETag, last-modified for the FASTA, .fai and .gzi)
  Manifest record ensembl_grch38_reference_regional in acquisition/manifest.jsonl.

This script does not download the FASTA, so it cannot check the file's published checksum; what it checks is that two reads agree
and that every region has the requested length and lies inside its sequence (validate.py repeats this from the saved .fai). The
same file is mirrored as target ensembl_grch38_dna_toplevel, whose BSD checksum is verified, and compat_reference.py compares
every region read here with the mirrored copy.

  .venv/bin/python acquisition/fetch_reference_regional.py [--workers 4]      (requires pysam; see requirements.txt)
"""
import argparse
import collections
import concurrent.futures
import datetime
import gzip
import hashlib
import os
import sys
import time
import urllib.request

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, '..'))
sys.path.insert(0, HERE)
import acquire  # noqa: E402
import mavedb_snvs  # noqa: E402

OUT = os.path.join(ROOT, 'data', 'external', 'reference_sequence')
CACHE = os.path.join(OUT, 'index_cache')
URL = 'https://ftp.ensembl.org/pub/release-116/fasta/homo_sapiens/dna_index/Homo_sapiens.GRCh38.dna.toplevel.fa.gz'
FAI_COPY = os.path.join(CACHE, os.path.basename(URL) + '.fai')


def now():
    return datetime.datetime.now(datetime.timezone.utc).strftime('%Y-%m-%dT%H:%M:%SZ')


def head(url):
    req = urllib.request.Request(url, method='HEAD', headers={'User-Agent': acquire.UA})
    with urllib.request.urlopen(req, timeout=60) as r:
        return {k.lower(): v for k, v in r.headers.items()}


def keep_fai():
    """Save a copy of the published .fai, checked against the server's Content-Length; return ({sequence: length}, sha256 of the copy)."""
    req = urllib.request.Request(URL + '.fai', headers={'User-Agent': acquire.UA})
    with urllib.request.urlopen(req, timeout=60) as r:
        data = r.read()
    size = int(head(URL + '.fai')['content-length'])
    if len(data) != size:
        raise RuntimeError(f'.fai: received {len(data)} bytes, server says {size}')
    os.makedirs(CACHE, exist_ok=True)
    with open(FAI_COPY, 'wb') as f:
        f.write(data)
    fai = {}
    for line in data.decode().splitlines():
        p = line.split('\t')
        fai[p[0]] = int(p[1])
    return fai, hashlib.sha256(data).hexdigest()


def read_regions(args):
    """Worker: all regions of one sequence, each read until two reads agree."""
    chrom, regs = args
    import pysam
    os.makedirs(CACHE, exist_ok=True)
    os.chdir(CACHE)                  # htslib may save remote index files in the working directory; keep them out of the repository
    state = dict(fa=pysam.FastaFile(URL))
    rows, retries = [], 0
    for (c, s, e) in regs:
        err = None
        for attempt in range(5):
            try:
                a = state['fa'].fetch(c, s - 1, e)
                b = state['fa'].fetch(c, s - 1, e)
                if a == b and len(a) == e - s + 1:
                    rows.append((c, s, e, a))
                    break
                err = f'reads differ or wrong length ({len(a)} / {len(b)} for {e - s + 1})'
            except (OSError, ValueError, KeyError) as ex:
                err = repr(ex)
            retries += 1
            time.sleep(2 ** attempt)
            state['fa'] = pysam.FastaFile(URL)
        else:
            raise RuntimeError(f'{c}:{s}-{e}: {err}')
    return chrom, rows, retries


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--workers', type=int, default=4)
    a = ap.parse_args()
    import pysam
    os.makedirs(OUT, exist_ok=True)
    fai, fai_sha = keep_fai()
    keys = mavedb_snvs.snv_alleles()
    by = collections.defaultdict(set)
    for c, p, r, al in keys:
        by[c].add(p)
    regs = mavedb_snvs.regions(by)
    std = [r for r in regs if r[0] in fai]
    skipped = collections.Counter(r[0] for r in regs if r[0] not in fai)
    beyond = [r for r in std if r[2] > fai[r[0]]]
    if beyond:
        raise RuntimeError(f'{len(beyond)} requested regions end beyond their sequence, first: {beyond[0]}')
    per = collections.defaultdict(list)
    for r in std:
        per[r[0]].append(r)
    print(f'{len(regs)} regions, {len(std)} on {len(per)} sequences of the FASTA, {sum(skipped.values())} skipped on {len(skipped)} contigs the FASTA does not hold', flush=True)
    results, retries = [], 0
    with concurrent.futures.ProcessPoolExecutor(max_workers=a.workers) as ex:
        for chrom, rows, r_ in ex.map(read_regions, [(c, per[c]) for c in sorted(per, key=lambda x: (-len(per[x]), x))]):
            results += rows
            retries += r_
            print(f'  {chrom}: {len(rows)} regions, {sum(len(x[3]) for x in rows)} bases, {r_} read retries', flush=True)
    results.sort(key=lambda r: (mavedb_snvs.contig_rank(r[0]), r[1]))
    outp = os.path.join(OUT, 'ensembl_grch38_toplevel_at_mavedb_regions.tsv.gz')
    with gzip.GzipFile(outp, 'wb', mtime=0) as g:
        g.write(b'chrom\tstart\tend\tsequence\n')
        for c, s, e, seq in results:
            g.write(f'{c}\t{s}\t{e}\t{seq}\n'.encode())
    table = os.path.join(HERE, 'manifest_files', 'ensembl_grch38_reference_regional.tsv')
    os.makedirs(os.path.dirname(table), exist_ok=True)
    with open(table, 'w') as f:
        f.write('url\tsize\tetag\tlast_modified\tsha256_of_copy_kept\n')
        for suffix in ('', '.fai', '.gzi'):
            h = head(URL + suffix)
            f.write('\t'.join([URL + suffix, h.get('content-length', ''), h.get('etag', ''), h.get('last-modified', ''), fai_sha if suffix == '.fai' else '']) + '\n')
    regions_sha = hashlib.sha256('\n'.join(f'{c}\t{s}\t{e}' for c, s, e, _ in results).encode()).hexdigest()
    rec = dict(id='ensembl_grch38_reference_regional', source='ensembl_genomes', url=URL, dest=os.path.relpath(outp, os.path.join(ROOT, 'data', 'external')), tier='A', kind='regional_slice', ok=bool(results),
               fetched_utc=now(), bytes=os.path.getsize(outp), rows=len(results), n_regions=len(results), bases=sum(len(x[3]) for x in results), sequences_read=len(per),
               skipped_contigs=dict(skipped), read_retries=retries, regions_sha256=regions_sha, fai_sha256=fai_sha, sha256=hashlib.sha256(open(outp, 'rb').read()).hexdigest(), md5=None, source_md5=None,
               md5_verified_against_source=False, size_verified_against_source=False, files_table=os.path.relpath(table, ROOT),
               tool=f'pysam {pysam.__version__} (htslib {pysam.__samtools_version__}) remote FASTA random access',
               note='range reads of the public Ensembl release-116 bgzip FASTA; the file is not downloaded, so its published checksum cannot be verified; every region was read twice with agreement')
    acquire.update_manifest(rec)
    print(f"ensembl_grch38_reference_regional: {len(results)} regions on {len(per)} sequences, {rec['bases']} bases, {rec['bytes'] / 1e3:.0f} KB written, {retries} read retries", flush=True)


if __name__ == '__main__':
    main()
