"""Mini source files in the real formats, built on the synthetic reference, for the source -> index -> matcher tests.
Everything here is synthetic. REF alleles are read from the synthetic reference so records are valid unless a test needs them not to be."""
import json
import os

import pysam

HERE = os.path.dirname(os.path.abspath(__file__))
REF = pysam.FastaFile(os.path.join(HERE, 'syn38.fa'))


def b(c, pos):
    return REF.fetch(c, pos - 1, pos).upper()


def alt(c, pos, skip=()):
    base = b(c, pos)
    comp = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}[base]
    return next(x for x in 'ACGT' if x not in (base, comp) and x not in skip)


def mini_clinvar_vcf(path):
    hdr = ['##fileformat=VCFv4.1', '##fileDate=2026-01-01', '##source=ClinVar', '##reference=SYN38',
           '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO']
    L = []

    def rec(chrom, pos, vid, ref, a, **info):
        L.append('\t'.join([chrom, str(pos), str(vid), ref, a, '.', '.', ';'.join(f'{k}={v}' for k, v in info.items())]))
    rec('1', 402, 1001, b('1', 402), alt('1', 402), ALLELEID=2001, CLNDISDB='MedGen:C0000001,OMIM:100001', CLNDN='Test_syndrome', CLNREVSTAT='reviewed_by_expert_panel',
        CLNSIG='Pathogenic', CLNSIGSCV='SCV000000001|SCV000000002', CLNVC='single_nucleotide_variant', CLNVI='ClinGen:CA000001', GENEINFO='GENEA:1001', MC='SO:0001583|missense_variant',
        ORIGIN='1', RS='111')
    rec('1', 410, 1002, b('1', 410), alt('1', 410), ALLELEID=2002, CLNDISDB='MedGen:C0000001', CLNDN='Test_syndrome', CLNREVSTAT='criteria_provided,_multiple_submitters,_no_conflicts',
        CLNSIG='Benign/Likely_benign', CLNSIGSCV='SCV000000003|SCV000000004', CLNVC='single_nucleotide_variant', GENEINFO='GENEA:1001', ORIGIN='1')
    rec('1', 420, 1003, b('1', 420), alt('1', 420), ALLELEID=2003, CLNDISDB='MedGen:CN169374', CLNDN='not_specified', CLNREVSTAT='criteria_provided,_single_submitter',
        CLNSIG='Uncertain_significance', CLNSIGSCV='SCV000000005', CLNVC='single_nucleotide_variant', GENEINFO='GENEA:1001|GENEB:1002', ORIGIN='1')
    rec('1', 430, 1004, b('1', 430), alt('1', 430), ALLELEID=2004, CLNDISDB='MedGen:C0000002|MedGen:C0000003', CLNDN='Cond_A|Cond_B', CLNREVSTAT='criteria_provided,_conflicting_classifications',
        CLNSIG='Conflicting_classifications_of_pathogenicity', CLNSIGCONF='Pathogenic(1)|Likely_benign(1)', CLNSIGSCV='SCV000000006|SCV000000007', CLNVC='single_nucleotide_variant', ORIGIN='1')
    rec('1', 105, 1005, 'AA', 'A', ALLELEID=2005, CLNDISDB='MedGen:C0000004', CLNDN='Homopolymer_disorder', CLNREVSTAT='criteria_provided,_single_submitter', CLNSIG='Likely_pathogenic',
        CLNSIGSCV='SCV000000008', CLNVC='Deletion', ORIGIN='1')                       # same deletion as the next record after normalisation
    rec('1', 101, 1006, 'CA', 'C', ALLELEID=2006, CLNDISDB='MedGen:C0000004', CLNDN='Homopolymer_disorder', CLNREVSTAT='no_assertion_criteria_provided', CLNSIG='Pathogenic',
        CLNSIGSCV='SCV000000009', CLNVC='Deletion', ORIGIN='1')
    s = REF.fetch('1', 499, 501).upper()
    rec('1', 500, 1007, s[:2], ''.join('C' if x != 'C' else 'G' for x in s[:2]), ALLELEID=2007, CLNDISDB='MedGen:C0000005', CLNDN='Mnv_disorder', CLNREVSTAT='criteria_provided,_single_submitter',
        CLNSIG='Pathogenic', CLNSIGSCV='SCV000000010', CLNVC='Microsatellite', ORIGIN='1')
    wrong = 'A' if b('1', 440) != 'A' else 'C'
    rec('1', 440, 1008, wrong, 'T' if wrong != 'T' else 'G', ALLELEID=2008, CLNDISDB='MedGen:C0000006', CLNDN='Wrong_ref', CLNREVSTAT='criteria_provided,_single_submitter', CLNSIG='Pathogenic',
        CLNSIGSCV='SCV000000011', CLNVC='single_nucleotide_variant', ORIGIN='1')
    rec('1', 450, 1009, b('1', 450), '.', ALLELEID=2009, CLNDISDB='MedGen:C0000007', CLNDN='No_alt', CLNREVSTAT='no_assertion_criteria_provided', CLNSIG='not_provided', CLNVC='Variation', ORIGIN='1')
    rec('1', 460, 1010, b('1', 460), alt('1', 460), ALLELEID=2010, ONC='Oncogenic', ONCREVSTAT='criteria_provided,_single_submitter', ONCDN='Tumour_x', ONCDISDB='MedGen:C0000008',
        ONCSCV='SCV000000012', SCI='Tier_I_-_Strong', SCIREVSTAT='criteria_provided,_single_submitter', SCIDN='Tumour_x', SCIDISDB='MedGen:C0000008', SCISCV='SCV000000013', CLNVC='single_nucleotide_variant', ORIGIN='2')
    rec('chrUn_xx', 100, 1011, 'A', 'T', ALLELEID=2011, CLNDISDB='MedGen:C0000009', CLNDN='Unplaced', CLNREVSTAT='criteria_provided,_single_submitter', CLNSIG='Pathogenic', CLNSIGSCV='SCV000000014', CLNVC='single_nucleotide_variant', ORIGIN='1')
    rec('1', 470, 1012, b('1', 470), alt('1', 470), AF_ESP='0.00039', AF_EXAC='0.00038', AF_TGP='0.0014', ALLELEID=2012, CLNDISDB='MedGen:CN169374', CLNDN='not_specified',
        CLNREVSTAT='criteria_provided,_single_submitter', CLNSIG='Likely_benign', CLNSIGSCV='SCV000000015', CLNVC='single_nucleotide_variant', ORIGIN='1')
    rec('1', 480, 1013, b('1', 480), 'N' + alt('1', 480), ALLELEID=2013, CLNDISDB='MedGen:C0000010', CLNDN='N_allele', CLNREVSTAT='criteria_provided,_single_submitter', CLNSIG='Pathogenic', CLNSIGSCV='SCV000000016', CLNVC='Insertion', ORIGIN='1')
    with open(path, 'w') as f:
        f.write('\n'.join(hdr + L) + '\n')


def write_all(out):
    """Write every mini source file into `out` (a temporary directory: synthetic VCF-like files are never committed)."""
    os.makedirs(out, exist_ok=True)
    mini_clinvar_vcf(os.path.join(out, 'mini_clinvar.vcf'))
    return out


if __name__ == '__main__':
    import sys
    print('written', write_all(sys.argv[1] if len(sys.argv) > 1 else 'sources'))
