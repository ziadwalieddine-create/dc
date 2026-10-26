#!/usr/bin/env python3
"""Match one sample's GRCh38 calls against the variant reference index, one chromosome at a time.

    python3 variant_index/match_wgs.py --vcf SAMPLE.vcf.gz --out data/results/SAMPLE [--index data/index] [--sample NAME]
                                       [--callable-bed HIGH_CONFIDENCE.bed] [--contigs 1,2,X] [--min-gq 20 --min-dp 10] [--external-haplotypes calls.json]

The VCF/gVCF is read from the path given and never copied, uploaded or queried against an external service. Results go to --out, which must lie under
data/ (git-ignored) or outside the repository. An indexed VCF (tabix) is read contig by contig; an unindexed one is read through once per contig.

What the tables contain is described in vindex/wgs.py. Views over the matches (priority lists, trait and pathway groupings) are made by
variant_index/make_views.py from the tables written here.
"""
import argparse
import json
import os
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

from vindex import wgs                                            # noqa: E402
from vindex.build.common import DATA, OUT_DEFAULT                 # noqa: E402
from vindex.refgenome import Reference                            # noqa: E402
from vindex.vcfio import QCPolicy                                 # noqa: E402

FASTA = os.path.join(DATA, 'reference_sequence', 'Homo_sapiens.GRCh38.dna.toplevel.fa.gz')
REPORT = os.path.join(DATA, 'refseq', 'GCF_000001405.40_GRCh38.p14_assembly_report.txt')


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--vcf', required=True)
    ap.add_argument('--out', required=True)
    ap.add_argument('--index', default=OUT_DEFAULT)
    ap.add_argument('--sample')
    ap.add_argument('--contigs', help='comma-separated canonical contigs (default: all primary contigs)')
    ap.add_argument('--callable-bed', help='callable/high-confidence regions; supports reference only under the stated assumption that the caller reports every variant inside them')
    ap.add_argument('--external-haplotypes', help='JSON {gene: [allele1, allele2]} from a dedicated caller (Cyrius, Aldy, StellarPGx, PharmCAT) for structural genes')
    ap.add_argument('--min-gq', type=int, default=20)
    ap.add_argument('--min-dp', type=int, default=10)
    ap.add_argument('--no-require-pass', action='store_true')
    ap.add_argument('--no-hash', action='store_true', help='skip the SHA-256 of the input file')
    ap.add_argument('--allow-assembly-mismatch', action='store_true')
    ap.add_argument('--fasta', default=FASTA)
    ap.add_argument('--assembly-report', default=REPORT)
    a = ap.parse_args(argv)
    ref = Reference(a.fasta, a.assembly_report)
    qc = QCPolicy(require_pass=not a.no_require_pass, min_gq=a.min_gq, min_dp=a.min_dp)
    ext = None
    if a.external_haplotypes:
        with open(a.external_haplotypes) as f:
            ext = json.load(f)
    t0 = time.time()
    log = lambda c, n: print(f'  chr{c}: {n["small_calls"]} calls, {n["svs"]} SV, {n["records"]} records ({time.time() - t0:.0f}s)', flush=True)
    man = wgs.run_sample(a.vcf, a.index, a.out, ref, sample=a.sample, qc=qc, contigs=a.contigs.split(',') if a.contigs else None, callable_bed=a.callable_bed, hash_input=not a.no_hash,
                         external_haplotypes=ext, allow_assembly_mismatch=a.allow_assembly_mismatch, progress=log)
    c = man['counts']
    print(json.dumps(dict(out=a.out, seconds=man['seconds'], calls=c['calls'], by_status=c['by_status'], rows=man['rows'], structural=c['structural'], haplotype_status=c['haplotype_status'],
                          high_value_rules=c['high_value_rules']), indent=1))
    return 0


if __name__ == '__main__':
    sys.exit(main())
