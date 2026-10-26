#!/usr/bin/env python3
"""Fetch the public Genome in a Bottle HG002 GRCh38 benchmark call set, a real, openly released callset used to test matching at whole-genome scale.

    python3 variant_index/fetch_public_testset.py [--dest data/external/validation/giab_hg002_grch38_v4.2.1]

HG002 (NA24385) is a reference sample whose donor consented to open release of the data, and NIST publishes the benchmark for reuse. It is not anyone's personal genome in this
project and it is not part of the evidence index; the files stay under data/ (git-ignored). The script records size and SHA-256 in variant_index/public_test_inputs.json the
first time and verifies against that record afterwards, and checks that the VCF is GRCh38, holds one sample named HG002, and is indexed. Nothing is overwritten on a mismatch.

Source: https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/AshkenazimTrio/HG002_NA24385_son/NISTv4.2.1/GRCh38/
Terms: GIAB data are released for open reuse (NIST / GIAB README: publicly released data are freely available for reuse without embargo); cite Zook et al. 2019 (Nat Biotechnol,
v4.2.1: Wagner et al. 2022, Cell Genomics) when using the benchmark.
"""
import argparse
import hashlib
import json
import os
import subprocess
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

from vindex.build.common import DATA, REPO                           # noqa: E402
from vindex.refgenome import Reference                              # noqa: E402

BASE = 'https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/AshkenazimTrio/HG002_NA24385_son/NISTv4.2.1/GRCh38/'
FILES = ['HG002_GRCh38_1_22_v4.2.1_benchmark.vcf.gz', 'HG002_GRCh38_1_22_v4.2.1_benchmark.vcf.gz.tbi', 'HG002_GRCh38_1_22_v4.2.1_benchmark_noinconsistent.bed']
RECORD = os.path.join(HERE, 'public_test_inputs.json')
FASTA = os.path.join(DATA, 'reference_sequence', 'Homo_sapiens.GRCh38.dna.toplevel.fa.gz')
REPORT = os.path.join(DATA, 'refseq', 'GCF_000001405.40_GRCh38.p14_assembly_report.txt')


def sha256_file(p):
    h = hashlib.sha256()
    with open(p, 'rb') as f:
        for b in iter(lambda: f.read(1 << 22), b''):
            h.update(b)
    return h.hexdigest()


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--dest', default=os.path.join(DATA, 'validation', 'giab_hg002_grch38_v4.2.1'))
    a = ap.parse_args(argv)
    os.makedirs(a.dest, exist_ok=True)
    rec = {}
    if os.path.exists(RECORD):
        with open(RECORD) as f:
            rec = json.load(f)
    files = rec.setdefault('files', {})
    rec.update(id='giab_hg002_grch38_v4.2.1', source=BASE, note=__doc__.split('\n\n')[0].strip())
    for name in FILES:
        p = os.path.join(a.dest, name)
        if not os.path.exists(p):
            subprocess.run(['curl', '-sSL', '--fail', '--retry', '3', '-o', p + '.part', BASE + name], check=True)
            os.replace(p + '.part', p)
        got = dict(bytes=os.path.getsize(p), sha256=sha256_file(p))
        if name in files and (files[name]['sha256'] != got['sha256'] or files[name]['bytes'] != got['bytes']):
            raise SystemExit(f'{name}: size or SHA-256 differs from the recorded values; refusing to continue (recorded {files[name]}, found {got})')
        files.setdefault(name, dict(got, retrieved_utc=time.strftime('%Y-%m-%dT%H:%M:%SZ', time.gmtime())))
    import pysam
    ref = Reference(FASTA, REPORT)
    vf = pysam.VariantFile(os.path.join(a.dest, FILES[0]))
    asm = ref.assembly_check({c.name: c.length for c in vf.header.contigs.values()})
    samples = list(vf.header.samples)
    if not asm.get('ok') or samples != ['HG002']:
        raise SystemExit(f'unexpected VCF: assembly check {asm}, samples {samples}')
    n = sum(1 for _ in vf)
    rec['validation'] = dict(assembly_ok=True, build_guess=asm.get('build_guess'), samples=samples, records=n, indexed=os.path.exists(os.path.join(a.dest, FILES[1])))
    with open(RECORD, 'w') as f:
        json.dump(rec, f, indent=1, sort_keys=True)
    print(json.dumps(rec['validation'], indent=1))
    return 0


if __name__ == '__main__':
    sys.exit(main())
