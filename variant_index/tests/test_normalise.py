import os
import random
import tempfile
import unittest

from helpers import syn_reference, find_tool, run_tool, vcf_text, FIXTURES  # noqa: F401
from vindex.normalise import normalise, atomise, classify
from vindex.refgenome import Reference, PRIMARY


class TestNormalise(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.ref = syn_reference()

    def n(self, pos, ref, alt, chrom='1'):
        return normalise(self.ref, chrom, pos, ref, alt)

    def test_snv_is_unchanged(self):
        base = self.ref.fetch('1', 402, 402)
        r = self.n(402, base, 'A' if base != 'A' else 'C')
        self.assertTrue(r.ok)
        self.assertEqual((r.pos, r.vtype, r.shifted), (402, 'snv', False))

    def test_homopolymer_deletion_all_representations_agree(self):
        # reference 1:100-111 = G C A A A A A A A A G T ; one A deleted, written at every possible position
        self.assertEqual(self.ref.fetch('1', 100, 111), 'GCAAAAAAAAGT')
        keys = set()
        for pos in range(101, 110):           # anchor C at 101, As at 102..109, G at 110
            ref = self.ref.fetch('1', pos, pos + 1)
            if ref[1] != 'A':
                continue
            keys.add(self.n(pos, ref, ref[0]).key)
        self.assertEqual(keys, {'1:101:CA:C'})

    def test_homopolymer_insertion_all_representations_agree(self):
        keys = set()
        for pos in range(101, 110):
            b = self.ref.fetch('1', pos, pos)
            if b == 'A' or pos == 101:
                keys.add(self.n(pos, b, b + 'A').key)
        self.assertEqual(keys, {'1:101:C:CA'})

    def test_tandem_repeat_unit_deletion_left_aligned(self):
        self.assertEqual(self.ref.fetch('1', 200, 217), 'TCAGCAGCAGCAGCAGTT')
        keys = set()
        for shift in range(0, 12):             # delete one CAG unit at any phase inside the repeat
            pos = 200 + shift
            span = self.ref.fetch('1', pos, pos + 3)
            keys.add(self.n(pos, span, span[0]).key)
        # every representation that deletes 3 bases inside the repeat normalises to one of the three phases only if the deleted
        # sequence differs; deleting a rotation of CAG leaves the same resulting sequence, so all collapse to the leftmost one
        left = {k for k in keys}
        self.assertEqual(len(left), 1, left)
        self.assertEqual(next(iter(left)), '1:200:TCAG:T')

    def test_dinucleotide_repeat_insertion(self):
        self.assertEqual(self.ref.fetch('1', 300, 313), 'CGTGTGTGTGTGAC')
        a = self.n(301, 'G', 'GTG')
        b = self.n(307, 'T', 'TGT') if self.ref.fetch('1', 307, 307) == 'T' else None
        self.assertEqual(a.key, '1:300:C:CGT')
        if b:
            self.assertEqual(b.key, a.key)

    def test_ref_mismatch_is_reported_not_normalised(self):
        base = self.ref.fetch('1', 402, 402)
        wrong = 'A' if base != 'A' else 'C'
        r = self.n(402, wrong, 'T' if wrong != 'T' else 'G')
        self.assertEqual(r.status, 'ref_mismatch')
        self.assertFalse(r.ok)

    def test_mnv_kept_and_trimmed_padding_removed(self):
        s = self.ref.fetch('1', 500, 505)              # GATTAC
        mnv = self.n(500, s[:2], 'CC' if s[:2] != 'CC' else 'GG')
        self.assertEqual(mnv.vtype, 'mnv')
        # padded SNV: shared flanks are trimmed
        padded = self.n(500, s[:3], s[0] + ('C' if s[1] != 'C' else 'G') + s[2])
        self.assertEqual((padded.pos, padded.vtype), (501, 'snv'))

    def test_atomise_mnv(self):
        s = self.ref.fetch('1', 500, 503)
        alt = ''.join('C' if b != 'C' else 'G' for b in s)
        m = self.n(500, s, alt)
        parts = atomise(m)
        self.assertEqual([p.pos for p in parts], [500, 501, 502, 503])
        self.assertTrue(all(p.vtype == 'snv' for p in parts))

    def test_lowercase_and_invalid_alleles(self):
        base = self.ref.fetch('1', 402, 402)
        alt = 'A' if base != 'A' else 'C'
        self.assertTrue(self.n(402, base.lower(), alt.lower()).ok)
        self.assertEqual(self.n(402, base, 'R').status, 'invalid_allele')
        self.assertEqual(self.n(402, base, '<DEL>').status, 'symbolic')
        self.assertEqual(self.n(402, base, '*').status, 'spanning_deletion')
        self.assertEqual(self.n(402, base, '.').status, 'no_alt_allele')
        self.assertEqual(self.n(402, base, base).status, 'not_a_variant')

    def test_out_of_contig(self):
        self.assertEqual(self.n(3000, 'AC', 'A').status, 'out_of_contig')

    def test_insertion_at_contig_start_is_right_anchored(self):
        first = self.ref.fetch('1', 1, 1)
        r = self.n(1, first, first + 'T' if first == 'T' else first + 'A')
        self.assertTrue(r.ok)
        # not a symbolic failure; if the insertion cannot be left-extended it is anchored on the right
        self.assertGreaterEqual(r.pos, 1)

    def test_idempotent(self):
        rnd = random.Random(1)
        for _ in range(500):
            pos = rnd.randint(60, 2900)
            k = rnd.randint(1, 5)
            ref = self.ref.fetch('1', pos, pos + k)
            alt = ref[0] + ''.join(rnd.choice('ACGT') for _ in range(rnd.randint(0, 4)))
            r = self.n(pos, ref, alt)
            if r.ok:
                r2 = self.n(r.pos, r.ref, r.alt)
                self.assertEqual(r2.key, r.key)
                self.assertFalse(r2.shifted)

    def test_classify(self):
        self.assertEqual(classify('A', 'G'), 'snv')
        self.assertEqual(classify('AC', 'GT'), 'mnv')
        self.assertEqual(classify('A', 'AT'), 'ins')
        self.assertEqual(classify('AT', 'A'), 'del')
        self.assertEqual(classify('AT', 'GCC'), 'complex')


class TestReference(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.ref = syn_reference()

    def test_contig_aliases(self):
        self.assertEqual(self.ref.contig('chr1'), '1')
        self.assertEqual(self.ref.contig('chrX'), 'X')
        self.assertEqual(self.ref.contig('chrM'), 'MT')
        self.assertEqual(self.ref.contig('M'), 'MT')
        self.assertIsNone(self.ref.contig('chr3'))

    def test_assembly_check(self):
        good = self.ref.assembly_check({'chr1': 3000, 'chr2': 1500, 'chrM': 400})
        self.assertTrue(good['ok'])
        self.assertEqual(good['build_guess'], 'SYN38')
        other = self.ref.assembly_check({'chr1': 249250621, 'chr2': 243199373})
        self.assertFalse(other['ok'])
        self.assertTrue(other['mismatched'])
        nolen = self.ref.assembly_check({'chr1': None})
        self.assertEqual(nolen['no_length'], ['chr1'])
        unk = self.ref.assembly_check({'chr1': 3000, 'HLA-A*01:01': 3503})
        self.assertEqual(unk['unmapped'], ['HLA-A*01:01'])

    def test_primary_constant(self):
        self.assertEqual(len(PRIMARY), 25)


@unittest.skipUnless(find_tool('bcftools'), 'bcftools not available (set VINDEX_TOOLS_BIN)')
class TestAgainstBcftools(unittest.TestCase):
    """The established tool is the arbiter: `bcftools norm -f` must give the same keys as `normalise` for random variants."""

    @classmethod
    def setUpClass(cls):
        cls.ref = syn_reference()

    def test_random_variants_match_bcftools_norm(self):
        rnd = random.Random(7)
        recs, expect = [], []
        feature_spots = [(100, 112), (200, 218), (300, 314), (500, 520), (700, 720)]
        for i in range(600):
            lo, hi = rnd.choice(feature_spots) if rnd.random() < 0.7 else (50, 2900)
            pos = rnd.randint(lo, hi)
            kind = rnd.choice(['snv', 'ins', 'del', 'mnv', 'delins', 'padded_snv'])
            if kind == 'snv':
                ref = self.ref.fetch('1', pos, pos)
                alt = rnd.choice([b for b in 'ACGT' if b != ref])
            elif kind == 'ins':
                ref = self.ref.fetch('1', pos, pos)
                alt = ref + ''.join(rnd.choice('ACGT') for _ in range(rnd.randint(1, 6)))
            elif kind == 'del':
                k = rnd.randint(1, 6)
                ref = self.ref.fetch('1', pos, pos + k)
                alt = ref[0]
            elif kind == 'mnv':
                k = rnd.randint(2, 4)
                ref = self.ref.fetch('1', pos, pos + k - 1)
                alt = ''.join(rnd.choice([b for b in 'ACGT' if b != r]) for r in ref)
            elif kind == 'delins':
                k = rnd.randint(2, 5)
                ref = self.ref.fetch('1', pos, pos + k - 1)
                alt = ''.join(rnd.choice('ACGT') for _ in range(rnd.randint(1, 4)))
                if alt == ref:
                    continue
            else:
                ref = self.ref.fetch('1', pos, pos + 2)
                mid = rnd.choice([b for b in 'ACGT' if b != ref[1]])
                alt = ref[0] + mid + ref[2]
            recs.append(('1', pos, ref, alt, 'PASS', 50, ['0/1:50:30:15,15'], '.', f'v{i}'))
        td = tempfile.mkdtemp()
        vin, vout = os.path.join(td, 'in.vcf'), os.path.join(td, 'out.vcf')
        with open(vin, 'w') as f:
            f.write(vcf_text(recs))
        rc, out, err = run_tool([find_tool('bcftools'), 'norm', '-f', os.path.join(FIXTURES, 'syn38.fa'), '-c', 'x', '-o', vout, vin])
        self.assertEqual(rc, 0, err)
        got = {}
        with open(vout) as fh:
            for line in fh:
                if line.startswith('#'):
                    continue
                f = line.split('\t')
                got[f[2]] = f'{f[0]}:{f[1]}:{f[3]}:{f[4]}'
        diffs, shifted = [], 0
        for r in recs:
            mine = normalise(self.ref, r[0], r[1], r[2], r[3])
            shifted += mine.shifted
            if got.get(r[8]) != mine.key:
                diffs.append((r[:4], got.get(r[8]), mine.key))
        self.assertGreater(shifted, 100, 'the comparison must include many representations that normalisation changes')
        self.assertGreater(len(recs), 500)
        self.assertEqual(diffs, [], f'{len(diffs)} of {len(recs)} differ from bcftools norm; first: {diffs[:3]}')


if __name__ == '__main__':
    unittest.main()
