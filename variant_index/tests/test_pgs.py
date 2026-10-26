import os
import tempfile
import unittest

from helpers import syn_reference, vcf_text, find_tool, run_tool, FIXTURES
from vindex.vcfio import read_vcf
from vindex.match_small import MatchContext
from vindex.pgs import read_scoring_file, score

HEADER = '#pgs_id=PGSTEST1\n#genome_build=SYN38\n#HmPOS_build=SYN38\n'
COLS = 'rsID\tchr_name\tchr_position\teffect_allele\tother_allele\teffect_weight\thm_chr\thm_pos\n'


class TestPgs(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.ref = syn_reference()
        cls.td = tempfile.mkdtemp()
        cls.b = lambda self, pos, c='1': self.ref.fetch(c, pos, pos)
        # alternate allele that does not form a palindromic (A/T, C/G) pair with the reference base
        cls.alt = lambda self, pos, skip=(), c='1': next(x for x in 'ACGT' if x != cls.ref.fetch(c, pos, pos) and x != {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}[cls.ref.fetch(c, pos, pos)] and x not in skip)

    def scoring(self, rows):
        p = os.path.join(self.td, 's.txt')
        with open(p, 'w') as f:
            f.write(HEADER + COLS + ''.join(f'rs{i}\t{c}\t{pos}\t{e}\t{o}\t{w}\t{c}\t{pos}\n' for i, (c, pos, e, o, w) in enumerate(rows)))
        return read_scoring_file(p)

    def ctx(self, recs):
        p = os.path.join(self.td, 'g.vcf')
        with open(p, 'w') as f:
            f.write(vcf_text(recs))
        return MatchContext(read_vcf(p, self.ref))

    def test_dosage_orientation_and_missingness(self):
        # positions 402 (het alt), 410 (hom alt), 420 (reference block), 430 (absent from the VCF), 440 (filtered)
        a402, a410, a440 = self.alt(402), self.alt(410), self.alt(440)
        recs = [('1', 402, self.b(402), a402, 'PASS', 50, ['0/1:50:30:15,15']),
                ('1', 410, self.b(410), a410, 'PASS', 50, ['1/1:60:30:0,30']),
                ('1', 420, self.b(420), '<NON_REF>', 'PASS', 0, ['0/0:60:30:30,0'], 'END=420', '.', 'GT:GQ:DP:AD'),
                ('1', 440, self.b(440), a440, 'LowQual', 5, ['0/1:10:8:4,4'])]
        sf = self.scoring([('1', 402, a402, self.b(402), 0.5),       # effect = alt, het  -> 1 x 0.5
                           ('1', 410, self.b(410), a410, 0.25),      # effect = ref, hom alt -> 0 x 0.25
                           ('1', 420, self.alt(420), self.b(420), 1.0),   # effect = alt, verified reference -> 0
                           ('1', 430, self.alt(430), self.b(430), 2.0),   # absent: missing, NOT zero
                           ('1', 440, a440, self.b(440), 4.0)])      # filtered: missing
        res = score(sf, self.ctx(recs), self.ref)
        self.assertEqual((res.n_weights, res.n_scored, res.n_missing), (5, 3, 2))
        self.assertAlmostEqual(res.score_sum, 0.5)
        self.assertEqual(res.reasons, {'absent_from_callset': 1, 'call_fails_quality': 1})
        self.assertAlmostEqual(res.abs_weight_scored_fraction, (0.5 + 0.25 + 1.0) / (0.5 + 0.25 + 1.0 + 2.0 + 4.0))
        self.assertIn('weight_coverage_below_0.9', res.flags)
        self.assertIn('uncalibrated', res.flags)
        self.assertEqual(res.ancestry['status'], 'unassessed')

    def test_strand_handling(self):
        comp = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
        pos = 402
        rb, ab = self.b(pos), self.alt(pos)
        recs = [('1', pos, rb, ab, 'PASS', 50, ['0/1:50:30:15,15'])]
        flipped = self.scoring([('1', pos, comp[ab], comp[rb], 1.0)])            # alleles given on the other strand
        res = score(flipped, self.ctx(recs), self.ref)
        self.assertEqual((res.n_scored, res.score_sum), (1, 1.0))
        self.assertEqual(res.per_weight[0][-1], 'flipped')
        pal = self.scoring([('1', 600, 'A', 'T', 1.0)])                           # palindromic: refused
        res2 = score(pal, self.ctx(recs), self.ref)
        self.assertEqual(res2.reasons, {'palindromic_strand_ambiguous': 1})
        others = [x for x in 'ACGT' if x not in (rb, ab, comp[rb], comp[ab])]
        wrong = self.scoring([('1', pos, others[0], others[0] if len(others) == 1 else others[1], 1.0)]) if len(others) >= 1 else None
        if wrong is not None:
            self.assertEqual(score(wrong, self.ctx(recs), self.ref).n_scored, 0)

    def test_indel_weight(self):
        recs = [('1', 105, 'AA', 'A', 'PASS', 50, ['0/1:50:30:15,15'])]     # deletion written right-shifted
        sf = self.scoring([('1', 101, 'C', 'CA', 0.7)])                       # weight: effect allele C (deleted form) vs other CA
        res = score(sf, self.ctx(recs), self.ref)
        # effect allele 'C' is the deletion allele (REF CA -> ALT C): het => 1 copy
        self.assertEqual((res.n_scored, res.score_sum), (1, 0.7))

    def test_unsupported_weight_types_are_excluded(self):
        p = os.path.join(self.td, 'd.txt')
        with open(p, 'w') as f:
            f.write(HEADER + COLS.rstrip('\n') + '\tis_dominant\n' + f'rs1\t1\t402\t{self.alt(402)}\t{self.b(402)}\t1.0\t1\t402\tTRUE\n')
        sf = read_scoring_file(p)
        res = score(sf, self.ctx([('1', 402, self.b(402), self.alt(402), 'PASS', 50, ['0/1:50:30:15,15'])]), self.ref)
        self.assertEqual(res.reasons, {'unsupported_weight_type:dominant': 1})

    def test_ancestry_and_calibration_gates(self):
        sf = self.scoring([('1', 402, self.alt(402), self.b(402), 0.5)])
        ctx = self.ctx([('1', 402, self.b(402), self.alt(402), 'PASS', 50, ['0/1:50:30:15,15'])])
        meta = dict(development='European:100', evaluation='European:72.7|Not Reported:18.2|East Asian:9.1')
        r1 = score(sf, ctx, self.ref, sample_ancestry='African', ancestry_meta=meta)
        self.assertEqual(r1.ancestry['status'], 'not_represented')
        self.assertIn('ancestry_not_represented', r1.flags)
        r2 = score(sf, ctx, self.ref, sample_ancestry='European', ancestry_meta=meta, reference_distribution=dict(mean=0.0, sd=0.25))
        self.assertEqual(r2.ancestry['status'], 'matched')
        self.assertAlmostEqual(r2.calibration['z'], 2.0)
        self.assertNotIn('uncalibrated', r2.flags)

    @unittest.skipUnless(find_tool('plink2'), 'plink2 not available (set VINDEX_TOOLS_BIN)')
    def test_dosage_arithmetic_matches_plink2(self):
        """Established-tool check: plink2 --score on the same genotypes gives the same sum (all variants known, diploid autosome)."""
        recs, rows, rnd_gts = [], [], ['0/1', '1/1', '0/1', '1/1', '0/1', '1/1', '0/1']
        pos_list = [402, 410, 460, 470, 480, 490, 520]
        for pos, gt in zip(pos_list, rnd_gts):
            a = self.alt(pos)
            recs.append(('1', pos, self.b(pos), a, 'PASS', 50, [f'{gt}:60:30:10,20']))
        weights = [0.1, -0.2, 0.3, 0.05, -0.4, 0.25, 0.15]
        for i, (pos, w) in enumerate(zip(pos_list, weights)):
            a = self.alt(pos)
            rows.append(('1', pos, a if i % 2 == 0 else self.b(pos), self.b(pos) if i % 2 == 0 else a, w))     # alternate effect allele = alt / ref
        sf = self.scoring(rows)
        res = score(sf, self.ctx(recs), self.ref)
        self.assertEqual(res.n_missing, 0)
        td = tempfile.mkdtemp()
        vcf, sc = os.path.join(td, 'c.vcf'), os.path.join(td, 'sc.tsv')
        with open(vcf, 'w') as f:
            f.write(vcf_text(recs))
        with open(sc, 'w') as f:
            for (c, pos, e, o, w) in rows:
                f.write(f'{c}:{pos}:{self.b(pos)}:{self.alt(pos)}\t{e}\t{w}\n')
        rc, out, err = run_tool([find_tool('plink2'), '--vcf', vcf, '--set-all-var-ids', '@:#:$r:$a', '--score', sc, '1', '2', '3', 'no-mean-imputation', 'cols=+scoresums', 'ignore-dup-ids', '--out', os.path.join(td, 'o')])
        self.assertEqual(rc, 0, out + err)
        with open(os.path.join(td, 'o.sscore')) as f:
            hdr = f.readline().lstrip('#').split()
            vals = f.readline().split()
        d = dict(zip(hdr, vals))
        self.assertAlmostEqual(float(d['SCORE1_SUM']), res.score_sum, places=6)


if __name__ == '__main__':
    unittest.main()
