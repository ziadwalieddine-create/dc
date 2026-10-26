import os
import tempfile
import unittest

from helpers import syn_reference, vcf_text
from vindex.vcfio import read_vcf, QCPolicy


def write(td, text, name='s.vcf'):
    p = os.path.join(td, name)
    with open(p, 'w') as f:
        f.write(text)
    return p


class TestVcfIO(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.ref = syn_reference()
        cls.td = tempfile.mkdtemp()

    def read(self, recs, **kw):
        fmt = kw.pop('fmt', 'GT:GQ:DP:AD')
        p = write(self.td, vcf_text(recs, fmt=fmt, **{k: v for k, v in kw.items() if k in ('samples', 'contigs')}))
        return read_vcf(p, self.ref, **{k: v for k, v in kw.items() if k in ('sample', 'qc')})

    def base(self, pos):
        return self.ref.fetch('1', pos, pos)

    def test_multiallelic_split_keeps_genotype_phase_quality(self):
        b = self.base(402)
        a1, a2 = [x for x in 'ACGT' if x != b][:2]
        cs = self.read([('1', 402, b, f'{a1},{a2}', 'PASS', 61.5, ['1|2:60:40:0,20,20:100'], '.', '.', 'GT:GQ:DP:AD:PS')])
        self.assertEqual(len(cs.small), 2)
        c1, c2 = cs.small
        self.assertEqual((c1.alt, c2.alt), (a1, a2))
        for c in cs.small:
            self.assertEqual((c.dosage, c.other_alt, c.ref_copies, c.missing), (1, 1, 0, 0))
            self.assertEqual((c.gt, c.phased, str(c.phase_set), c.qual, c.gq, c.dp, c.ad_ref, c.ad_alt), ('1|2', True, '100', 61.5, 60, 40, 0, 20))
            self.assertTrue(c.multiallelic)
            self.assertTrue(c.qc_pass)
        self.assertEqual((c1.hap_alleles, c2.hap_alleles), ((1, 0), (0, 1)))
        self.assertEqual(cs.stats['split_multiallelic_records'], 1)

    def test_multiallelic_snv_and_deletion_normalise_independently(self):
        # 101 CA... homopolymer: the deletion written at 105 normalises to 101
        ref = self.ref.fetch('1', 105, 106)            # AA
        cs = self.read([('1', 105, ref, f'{ref[0]},{"C" if ref[0] != "C" else "G"}A', 'PASS', 50, ['1/2:50:30:0,15,15'])])
        keys = sorted(c.key for c in cs.small)
        self.assertIn('1:101:CA:C', keys)
        self.assertEqual(len(keys), 2)
        self.assertTrue(all(c.dosage == 1 and c.other_alt == 1 for c in cs.small))
        self.assertTrue(any(c.orig[1] == 105 and c.pos == 101 for c in cs.small), 'the original coordinates stay recorded')

    def test_failed_filter_and_low_quality_calls_are_kept_and_flagged(self):
        b = self.base(402)
        a = 'A' if b != 'A' else 'C'
        cs = self.read([('1', 402, b, a, 'LowQual', 5, ['0/1:10:8:4,4']), ('1', 410, self.base(410), 'A' if self.base(410) != 'A' else 'C', 'PASS', 40, ['0/1:12:5:2,3'])])
        self.assertEqual(len(cs.small), 2)
        self.assertFalse(cs.small[0].qc_pass)
        self.assertIn('FILTER:LowQual', cs.small[0].qc_flags)
        self.assertIn('GQ<20', cs.small[0].qc_flags)
        self.assertIn('GQ<20', cs.small[1].qc_flags)
        self.assertIn('DP<10', cs.small[1].qc_flags)

    def test_ref_mismatch_is_rejected_not_matched(self):
        b = self.base(402)
        wrong = 'A' if b != 'A' else 'C'
        alt = 'T' if wrong != 'T' else 'G'
        cs = self.read([('1', 402, wrong, alt, 'PASS', 50, ['0/1:50:30:15,15'])])
        self.assertEqual(cs.small, [])
        self.assertEqual([r.reason for r in cs.rejected], ['ref_mismatch'])
        self.assertEqual(cs.stats['rejected'], 1)

    def test_unmapped_contig_goes_to_unplaced(self):
        cs = self.read([('chrUn_KI270302v1', 100, 'A', 'T', 'PASS', 50, ['0/1:50:30:15,15'])], contigs=[('chrUn_KI270302v1', 2274), ('1', 3000)])
        self.assertEqual(cs.small, [])
        self.assertEqual([r.reason for r in cs.unplaced], ['contig_unmapped'])

    def test_no_calls(self):
        b = self.base(402)
        a = 'A' if b != 'A' else 'C'
        cs = self.read([('1', 402, b, a, 'PASS', 50, ['./.:.:.:.']), ('1', 410, self.base(410), 'A' if self.base(410) != 'A' else 'C', 'PASS', 50, ['./1:30:20:0,10'])])
        self.assertEqual(len(cs.nocalls), 2)
        self.assertEqual(len(cs.small), 1)                       # the half-call keeps its called allele ...
        self.assertFalse(cs.small[0].qc_pass)                    # ... but never passes QC
        self.assertIn('partial_no_call', cs.small[0].qc_flags)
        self.assertEqual(cs.small[0].missing, 1)

    def test_gvcf_blocks_and_hom_ref_genotype_sites(self):
        b = self.base(402)
        a = 'A' if b != 'A' else 'C'
        recs = [('1', 300, 'C', '<NON_REF>', 'PASS', 0, ['0/0:50:30:30,0'], 'END=399', '.', 'GT:GQ:DP:AD'),
                ('1', 400, 'A', '<NON_REF>', 'PASS', 0, ['0/0:10:30:30,0'], 'END=401', '.', 'GT:GQ:DP:AD'),    # GQ below the reference threshold
                ('1', 402, b, f'{a},<NON_REF>', 'PASS', 50, ['0/1:50:30:15,15,0'], '.'),
                ('1', 410, self.base(410), 'A' if self.base(410) != 'A' else 'C', 'PASS', 50, ['0/0:60:35:35,0'])]
        cs = self.read(recs)
        self.assertEqual([(r.start, r.end, r.confident, r.basis) for r in cs.refblocks],
                         [(300, 399, True, 'gvcf_block'), (400, 401, False, 'gvcf_block'), (410, 410, True, 'genotype_ref')])
        self.assertEqual([c.key for c in cs.small], [f'1:402:{b}:{a}'])
        self.assertEqual(cs.svs, [])

    def test_spanning_deletion_marker(self):
        ref = self.ref.fetch('1', 105, 106)
        cs = self.read([('1', 105, ref, f'{ref[0]},*', 'PASS', 50, ['1/2:50:30:0,15,15'])])
        self.assertEqual([(d.pos, d.copies) for d in cs.spandels], [(105, 1)])
        self.assertEqual(len(cs.small), 1)

    def test_haploid_calls(self):
        b = self.ref.fetch('X', 105, 105)
        a = 'A' if b != 'A' else 'C'
        cs = self.read([('X', 105, b, a, 'PASS', 50, ['1:50:30:0,30'])])
        c = cs.small[0]
        self.assertEqual((c.ploidy, c.dosage, c.gt), (1, 1, '1'))

    def test_symbolic_and_large_sequence_svs(self):
        recs = [('1', 1000, 'N', '<DEL>', 'PASS', 90, ['0/1:50:30:15,15'], 'SVTYPE=DEL;END=1500;SVLEN=-500;CIPOS=-10,10;CIEND=-5,5', 'sv1'),
                ('1', 1700, 'N', '<DUP>', 'PASS', 90, ['0/1:50:30:15,15'], 'SVTYPE=DUP;END=1800;IMPRECISE', 'sv2')]
        cs = self.read(recs)
        self.assertEqual(cs.small, [])
        d, u = cs.svs
        self.assertEqual((d.id, d.svtype, d.start, d.end, d.cipos, d.ciend, d.precise, d.dosage), ('sv1', 'DEL', 1000, 1500, (-10, 10), (-5, 5), True, 1))
        self.assertEqual((u.svtype, u.precise, u.end), ('DUP', False, 1800))
        bad = self.read([('1', 1000, 'N', '<DEL>', 'PASS', 90, ['0/1:50:30:15,15'], 'SVTYPE=DEL;END=1500;SVLEN=-900', 'sv3')]).svs[0]
        self.assertEqual((bad.end, bad.flags), (1500, ('END_SVLEN_inconsistent',)))
        ref = self.ref.fetch('1', 1200, 1260)
        cs2 = self.read([('1', 1200, ref, ref[0], 'PASS', 50, ['0/1:50:30:15,15'])])
        self.assertEqual(len(cs2.small), 1)
        self.assertEqual(len(cs2.svs), 1)
        self.assertEqual((cs2.svs[0].svtype, cs2.svs[0].source), ('DEL', 'sequence_resolved'))

    def test_assembly_mismatch_is_reported(self):
        b = self.base(402)
        cs = self.read([('1', 402, b, 'A' if b != 'A' else 'C', 'PASS', 50, ['0/1:50:30:15,15'])], contigs=[('1', 249250621)])
        self.assertFalse(cs.assembly['ok'])

    def test_multi_sample_requires_a_name(self):
        b = self.base(402)
        a = 'A' if b != 'A' else 'C'
        recs = [('1', 402, b, a, 'PASS', 50, ['0/1:50:30:15,15', '0/0:50:30:30,0'])]
        p = write(self.td, vcf_text(recs, samples=('S1', 'S2')), 'multi.vcf')
        with self.assertRaises(ValueError):
            read_vcf(p, self.ref)
        cs = read_vcf(p, self.ref, sample='S2')
        self.assertEqual(cs.small, [])
        self.assertEqual(len(cs.refblocks), 1)

    def test_adjacent_snvs_unphased_vs_phased(self):
        s = self.ref.fetch('1', 500, 501)
        a0 = 'A' if s[0] != 'A' else 'C'
        a1 = 'A' if s[1] != 'A' else 'C'
        recs = [('1', 500, s[0], a0, 'PASS', 50, ['0/1:50:30:15,15:7'], '.', '.', 'GT:GQ:DP:AD:PS'),
                ('1', 501, s[1], a1, 'PASS', 50, ['0/1:50:30:15,15:7'], '.', '.', 'GT:GQ:DP:AD:PS')]
        cs = self.read(recs)
        self.assertTrue(all(c.phase_set is None and not c.phased and c.hap_alleles is None for c in cs.small))
        recs2 = [(r[0], r[1], r[2], r[3], r[4], r[5], [r[6][0].replace('0/1', '0|1')]) + r[7:] for r in recs]
        cs2 = self.read(recs2)
        self.assertTrue(all(c.phased and str(c.phase_set) == '7' and c.hap_alleles == (0, 1) for c in cs2.small))


if __name__ == '__main__':
    unittest.main()
