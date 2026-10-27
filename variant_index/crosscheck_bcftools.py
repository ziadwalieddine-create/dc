#!/usr/bin/env python3
"""Cross-check a whole-genome run against bcftools, an independent normaliser.

    python3 variant_index/crosscheck_bcftools.py --vcf SAMPLE.vcf.gz --results data/results/SAMPLE --evidence-vcf data/external/clinvar/clinvar.vcf.gz [--contigs 21,22] [--sample NAME]

Recomputes, with bcftools only, the sample's carried small alleles and their intersection with a source VCF (ClinVar by default), and compares both with the tables
the matcher wrote. The report goes to stdout (or --out); it holds counts and a few example keys, never genotypes or the sample's calls. Exit status 1 if the sets differ.

    python3 variant_index/crosscheck_bcftools.py --clinvar-fields --index data/index [--evidence-vcf data/external/clinvar/clinvar.vcf.gz] [--contigs 21,22]

The second form needs no sample: it compares every ClinVar-VCF-derived evidence row of the index (classification, review status, condition text) with the VCF as read by bcftools.

    python3 variant_index/crosscheck_bcftools.py --reference-sites --vcf SAMPLE.vcf.gz --results data/results/SAMPLE --index data/index --callable-bed HIGH_CONFIDENCE.bed [--contigs 21,22]

The third form recomputes the sites at which the run verified a reference genotype by a callable region (indexed keys whose evidence is about the reference allele, inside the BED,
with no record of the sample's file and no allele it carries after bcftools norm) and compares them with the verified_reference links of the run.
"""
import argparse
import json
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

from vindex.build.common import DATA                              # noqa: E402
from vindex.crosscheck import clinvar_fields, crosscheck, reference_sites          # noqa: E402
from vindex.refgenome import Reference                            # noqa: E402

FASTA = os.path.join(DATA, 'reference_sequence', 'Homo_sapiens.GRCh38.dna.toplevel.fa.gz')
REPORT = os.path.join(DATA, 'refseq', 'GCF_000001405.40_GRCh38.p14_assembly_report.txt')


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--vcf')
    ap.add_argument('--results')
    ap.add_argument('--clinvar-fields', action='store_true', help='compare the index\'s ClinVar rows with the ClinVar VCF instead of checking a sample run')
    ap.add_argument('--index', default=os.path.join(os.path.dirname(HERE), 'data', 'index'))
    ap.add_argument('--evidence-vcf', default=os.path.join(DATA, 'clinvar', 'clinvar.vcf.gz'))
    ap.add_argument('--evidence-source', default='clinvar', help='source_id whose evidence_links the evidence VCF is compared with')
    ap.add_argument('--no-evidence', action='store_true')
    ap.add_argument('--reference-sites', action='store_true', help='recompute the sites at which the run verified a reference genotype (needs --vcf, --results, --index, --callable-bed)')
    ap.add_argument('--callable-bed', help='the callable-region BED the run was made with')
    ap.add_argument('--contigs', help='comma-separated canonical contigs (default: all primary)')
    ap.add_argument('--sample')
    ap.add_argument('--fasta', default=FASTA)
    ap.add_argument('--assembly-report', default=REPORT)
    ap.add_argument('--out')
    a = ap.parse_args(argv)
    ref = Reference(a.fasta, a.assembly_report)
    if a.clinvar_fields:
        rep = clinvar_fields(a.index, a.evidence_vcf, ref, a.contigs.split(',') if a.contigs else None)
    elif a.reference_sites:
        if not (a.vcf and a.results and a.callable_bed):
            ap.error('--reference-sites needs --vcf, --results and --callable-bed')
        rep = reference_sites(a.vcf, a.results, a.index, a.callable_bed, ref, a.contigs.split(',') if a.contigs else None, a.sample)
    else:
        if not (a.vcf and a.results):
            ap.error('--vcf and --results are required unless --clinvar-fields is given')
        rep = crosscheck(a.vcf, a.results, ref, None if a.no_evidence else a.evidence_vcf, a.evidence_source, a.contigs.split(',') if a.contigs else None, a.sample)
    text = json.dumps(rep, indent=1, sort_keys=True)
    if a.out:
        with open(a.out, 'w') as f:
            f.write(text + '\n')
    print(text)
    return 0 if rep['ok'] else 1


if __name__ == '__main__':
    sys.exit(main())
