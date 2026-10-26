"""Extract the score tables of the MaveDB score sets that have genomic (GRCh38) mappings from the bulk dump.

The dump (mavedb_zenodo_dump_v4, 1.8 GB) holds every score set's scores.csv; only 83 of 2,777 were unpacked in an earlier phase. MaveDB's own mapping
(annotations.csv, column mavedb.post_mapped_hgvs_g) places 372,258 variants of 200 score sets on GRCh38 with a nucleotide-level HGVS expression; their
scores are needed to turn those placements into evidence. This script streams the tar once (nothing else is unpacked), keeps the scores.csv members of
those 200 sets as gzip files under data/external/mavedb/scores_genomic_mapped/, and records size, SHA-256 and row count of each member in
acquisition/manifest_files/mavedb_scores_genomic_mapped.tsv. The member text is stored byte for byte, so a line number in the stored file is a line number
in the dump.

  python3 extract_mavedb_scores.py [--restart]
"""
import argparse
import csv
import datetime
import glob
import gzip
import hashlib
import io
import json
import os
import sys
import tarfile

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, '..'))
sys.path.insert(0, HERE)
import acquire  # noqa: E402

DATA = os.path.join(ROOT, 'data', 'external', 'mavedb')
TAR = os.path.join(DATA, 'mavedb-dump.2026062418131.tar.gz')
ANN = os.path.join(DATA, 'dump_v4_extract', 'mavedb-dump.2026062418131', 'csv')
OUT = os.path.join(DATA, 'scores_genomic_mapped')
TARGET = 'mavedb_scores_genomic_mapped'
TABLE = os.path.join(HERE, 'manifest_files', TARGET + '.tsv')


def mapped_sets():
    """urn (hyphen form) -> number of annotation rows that carry a post-mapped genomic HGVS."""
    out = {}
    for f in sorted(glob.glob(os.path.join(ANN, '*.annotations.csv'))):
        urn = os.path.basename(f).split('.')[0]
        n = 0
        with open(f, newline='', encoding='utf-8') as fh:
            for r in csv.DictReader(fh):
                if r.get('mavedb.post_mapped_hgvs_g', 'NA') not in ('NA', ''):
                    n += 1
        if n:
            out[urn] = n
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--restart', action='store_true')
    a = ap.parse_args()
    want = mapped_sets()
    print(f'{len(want)} score sets with genomic mappings, {sum(want.values()):,} mapped variants', flush=True)
    os.makedirs(OUT, exist_ok=True)
    rows, seen = [], set()
    with tarfile.open(TAR, mode='r|gz') as tf:
        for m in tf:
            base = os.path.basename(m.name)
            if base.startswith('._') or not base.endswith('.scores.csv'):
                continue
            urn = base[:-len('.scores.csv')]
            if urn not in want or not m.isfile():
                continue
            body = tf.extractfile(m).read()
            seen.add(urn)
            n_lines = body.count(b'\n') + (0 if body.endswith(b'\n') or not body else 1)
            dest = os.path.join(OUT, urn + '.scores.csv.gz')
            with gzip.GzipFile(dest, 'wb', compresslevel=6, mtime=0) as g:         # mtime 0: the stored bytes depend only on the member
                g.write(body)
            with open(dest, 'rb') as fh:
                stored = fh.read()
            rows.append((urn, m.name, len(body), n_lines - 1, hashlib.sha256(body).hexdigest(), len(stored), hashlib.sha256(stored).hexdigest(), want[urn]))
    rows.sort()
    missing = sorted(set(want) - seen)
    os.makedirs(os.path.dirname(TABLE), exist_ok=True)
    with open(TABLE, 'w') as f:
        f.write('urn\ttar_member\tmember_bytes\tscore_rows\tmember_sha256\tstored_bytes\tstored_sha256\tmapped_annotation_rows\n')
        for r in rows:
            f.write('\t'.join(str(x) for x in r) + '\n')
    rec = dict(id=TARGET, source='mavedb', url=None, dest='mavedb/scores_genomic_mapped', tier='A', kind='per_entity', ok=not missing and bool(rows),
               fetched_utc=datetime.datetime.now(datetime.timezone.utc).strftime('%Y-%m-%dT%H:%M:%SZ'), n_accessions=len(want), n_files=len(rows), bytes=sum(r[5] for r in rows),
               n_errors=len(missing), errors=missing[:5], files_table=os.path.relpath(TABLE, ROOT),
               sha256=hashlib.sha256('\n'.join('\t'.join(str(x) for x in r[:5]) for r in rows).encode()).hexdigest(), md5=None, source_md5=None,
               md5_verified_against_source=False, size_verified_against_source=False,
               note='scores.csv members of the MaveDB score sets that have post-mapped genomic HGVS in annotations.csv, streamed out of mavedb_zenodo_dump_v4 (itself MD5-verified)')
    acquire.update_manifest(rec)
    print(json.dumps({k: rec[k] for k in ('n_accessions', 'n_files', 'bytes', 'n_errors')}), flush=True)


if __name__ == '__main__':
    main()
