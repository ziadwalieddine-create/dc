"""Genomic SNV positions of the MaveDB human set, from the two mapping products that exist (see compat_tests.py).

  Ensembl release-116 MaveDB_variants.tsv.gz   every codon-compatible genomic allele for protein-level variants (GRCh38)
  MaveDB bulk dump v4 annotations              the MaveDB mapper's genomic HGVS (NC_0000NN.x:g.POSREF>ALT) for nucleotide-level variants

snv_alleles() returns the set of (chrom, pos, ref, alt) over both products, chromosomes written without a chr prefix, which
is the key space compat_tests.py joins on. It needs data/external/ensembl_vep/MaveDB_variants.tsv.gz and the extracted dump
annotations (acquire.py, then extract_dump.py).
"""
import csv
import glob
import gzip
import json
import os
import re

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, '..'))
D = os.path.join(ROOT, 'data', 'external')
csv.field_size_limit(1 << 28)
NC2CHR = {f'NC_{n:06d}': (str(n) if n <= 22 else {23: 'X', 24: 'Y'}[n]) for n in range(1, 25)}
HGVS_G = re.compile(r'^(NC_\d{6})\.\d+:g\.(\d+)([ACGT])>([ACGT])$')


def our_score_sets():
    sets = set()
    for fn in ('scrape/mavedb-human/score_set_details.jsonl', 'scrape/mavedb-human/score_set_details_accession_targets.jsonl'):
        for line in open(os.path.join(ROOT, fn)):
            sets.add(json.loads(line)['urn'])
    return sets


def snv_alleles():
    sets = our_score_sets()
    keys = set()
    with gzip.open(os.path.join(D, 'ensembl_vep', 'MaveDB_variants.tsv.gz'), 'rt', errors='replace') as f:
        hdr = f.readline().lstrip('#').rstrip('\n').split('\t')
        ci = {c: i for i, c in enumerate(hdr)}
        for line in f:
            p = line.split('\t', 14)
            if p[ci['urn']] in sets and len(p[3]) == 1 and len(p[4]) == 1:
                keys.add((p[0], int(p[1]), p[3], p[4]))
    ann_dir = glob.glob(os.path.join(D, 'mavedb', 'dump_v4_extract', '*', 'csv'))[0]
    for fp in glob.glob(os.path.join(ann_dir, 'urn-mavedb-*.annotations.csv')):
        urn = os.path.basename(fp)[:-len('.annotations.csv')].replace('urn-mavedb-', 'urn:mavedb:', 1)
        if urn not in sets:
            continue
        with open(fp, newline='', encoding='utf-8', errors='replace') as f:
            for row in csv.DictReader(f):
                g = row['mavedb.post_mapped_hgvs_g']
                if g != 'NA':
                    m = HGVS_G.match(g)
                    if m and m.group(1) in NC2CHR:
                        keys.add((NC2CHR[m.group(1)], int(m.group(2)), m.group(3), m.group(4)))
    return keys


def regions(positions_by_chrom, gap=150):
    """Merge sorted positions into [start, end] (1-based, inclusive) regions whose neighbours are at most `gap` apart."""
    out = []
    for chrom in sorted(positions_by_chrom):
        ps = sorted(positions_by_chrom[chrom])
        s = e = ps[0]
        for p in ps[1:]:
            if p - e <= gap:
                e = p
            else:
                out.append((chrom, s, e))
                s = e = p
        out.append((chrom, s, e))
    return out
