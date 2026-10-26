import os
import tempfile
import unittest

from helpers import syn_reference, vcf_text
from vindex.vcfio import read_vcf
from vindex.match_repeat import match_repeats, classify_allele

LOCI = [dict(locus_id='HTT', gene='HTT', chrom='1', start=1000, end=1040, normal_max=26, pathologic_min=36, inheritance='AD', disease='Huntington disease'),
        dict(locus_id='FMR1', gene='FMR1', chrom='X', start=500, end=530, normal_max=44, pathologic_min=200, inheritance='XL', disease='Fragile X')]


class TestRepeats(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.ref = syn_reference()
        cls.td = tempfile.mkdtemp()

    def read(self, recs):
        p = os.path.join(self.td, 'r.vcf')
        with open(p, 'w') as f:
            f.write(vcf_text(recs, fmt='GT:REPCN:REPCI'))
        return read_vcf(p, self.ref)

    def test_no_repeat_genotype_is_unresolved_with_exact_missing_evidence(self):
        cs = self.read([('1', 402, self.ref.fetch('1', 402, 402), 'A' if self.ref.fetch('1', 402, 402) != 'A' else 'C', 'PASS', 50, ['0/1:.:.'], '.', '.', 'GT:REPCN:REPCI')])
        res = match_repeats(cs, LOCI)
        self.assertTrue(all(r.status == 'UNRESOLVED' for r in res))
        self.assertIn('repeat-length genotype', res[0].missing)

    def test_expansionhunter_genotypes_are_classified(self):
        recs = [('1', 1000, 'N', '<STR19>,<STR41>', 'PASS', 0, ['1/2:19/41:19-19/38-44'], 'END=1040;RU=CAG;VARID=HTT', '.', 'GT:REPCN:REPCI')]
        cs = self.read(recs)
        self.assertEqual(len(cs.repeats), 1)
        self.assertEqual(cs.svs, [])
        res = {r.locus_id: r for r in match_repeats(cs, LOCI)}
        self.assertEqual(res['HTT'].status, 'GENOTYPED')
        self.assertEqual(res['HTT'].classes, ('normal_range', 'pathogenic_range'))
        self.assertEqual(res['FMR1'].status, 'UNRESOLVED')

    def test_borderline_confidence_interval(self):
        self.assertEqual(classify_allele(34, (30, 40), 26, 36), 'borderline_crossing_pathologic_minimum')
        self.assertEqual(classify_allele(30, (28, 32), 26, 36), 'intermediate')
        self.assertEqual(classify_allele(26, (24, 28), 26, 36), 'borderline_crossing_normal_maximum')
        self.assertEqual(classify_allele(None, None, 26, 36), 'unknown')


if __name__ == '__main__':
    unittest.main()
