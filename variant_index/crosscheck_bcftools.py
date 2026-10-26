#!/usr/bin/env python3
"""Cross-check a whole-genome run against bcftools, an independent normaliser.

    python3 variant_index/crosscheck_bcftools.py --vcf SAMPLE.vcf.gz --results data/results/SAMPLE --evidence-vcf data/external/clinvar/clinvar.vcf.gz [--contigs 21,22] [--sample NAME]

Recomputes, with bcftools only, the sample's carried small alleles and their intersection with a source VCF (ClinVar by default), and compares both with the tables
the matcher wrote. The report goes to stdout (or --out); it holds counts and a few example keys, never genotypes or the sample's calls. Exit status 1 if the sets differ.
"""
import argparse
import json
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

from vindex.build.common import DATA                              # noqa: E402
from vindex.crosscheck import crosscheck                          # noqa: E402
from vindex.refgenome import Reference                            # noqa: E402

FASTA = os.path.join(DATA, 'reference_sequence', 'Homo_sapiens.GRCh38.dna.toplevel.fa.gz')
REPORT = os.path.join(DATA, 'refseq', 'GCF_000001405.40_GRCh38.p14_assembly_report.txt')


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--vcf', required=True)
    ap.add_argument('--results', required=True)
    ap.add_argument('--evidence-vcf', default=os.path.join(DATA, 'clinvar', 'clinvar.vcf.gz'))
    ap.add_argument('--evidence-source', default='clinvar', help='source_id whose evidence_links the evidence VCF is compared with')
    ap.add_argument('--no-evidence', action='store_true')
    ap.add_argument('--contigs', help='comma-separated canonical contigs (default: all primary)')
    ap.add_argument('--sample')
    ap.add_argument('--fasta', default=FASTA)
    ap.add_argument('--assembly-report', default=REPORT)
    ap.add_argument('--out')
    a = ap.parse_args(argv)
    ref = Reference(a.fasta, a.assembly_report)
    rep = crosscheck(a.vcf, a.results, ref, None if a.no_evidence else a.evidence_vcf, a.evidence_source, a.contigs.split(',') if a.contigs else None, a.sample)
    text = json.dumps(rep, indent=1, sort_keys=True)
    if a.out:
        with open(a.out, 'w') as f:
            f.write(text + '\n')
    print(text)
    return 0 if rep['ok'] else 1


if __name__ == '__main__':
    sys.exit(main())
