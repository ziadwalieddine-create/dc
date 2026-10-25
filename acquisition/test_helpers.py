"""Unit tests for the small functions that every overlap and join number in this directory rests on.

  python3 -m unittest acquisition.test_helpers        (from the repository root)
  python3 acquisition/test_helpers.py

They use only synthetic inputs, so they need no downloaded data. The pysam and pyarrow code paths are covered by the
validators and compatibility tests instead.
"""
import os
import random
import sys
import unittest

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

import compat_reference  # noqa: E402
import intervals  # noqa: E402
import mavedb_snvs  # noqa: E402


class StabbingQuery(unittest.TestCase):
    def brute(self, iv, pos0):
        return any(s <= pos0 < e for s, e in iv)

    def test_matches_brute_force_on_random_intervals(self):
        rnd = random.Random(7)
        for _ in range(300):
            iv = []
            for _ in range(rnd.randint(1, 40)):
                s = rnd.randint(0, 200)
                iv.append((s, s + rnd.randint(1, 60)))     # nested, overlapping and adjacent intervals all occur
            idx = intervals.build_index(iv)
            for pos0 in range(-2, 270):
                self.assertEqual(intervals.covered(idx, pos0), self.brute(iv, pos0), (iv, pos0))

    def test_half_open_boundaries(self):
        idx = intervals.build_index([(10, 12)])
        self.assertFalse(intervals.covered(idx, 9))
        self.assertTrue(intervals.covered(idx, 10))
        self.assertTrue(intervals.covered(idx, 11))
        self.assertFalse(intervals.covered(idx, 12))

    def test_a_long_interval_hides_behind_later_short_ones(self):
        idx = intervals.build_index([(0, 100), (10, 11), (20, 21)])
        self.assertTrue(intervals.covered(idx, 50))       # needs the prefix maximum, not just the nearest start


class RegionMerging(unittest.TestCase):
    def test_gap_rule_is_inclusive(self):
        out = mavedb_snvs.regions({'1': {100, 150, 300, 1000}, 'X': {5}}, gap=150)
        self.assertEqual(out, [('1', 100, 300), ('1', 1000, 1000), ('X', 5, 5)])           # 150 -> 300 is exactly the gap
        out = mavedb_snvs.regions({'1': {100, 150, 301}}, gap=150)
        self.assertEqual(out, [('1', 100, 150), ('1', 301, 301)])                          # 151 apart: separate regions

    def test_single_position(self):
        self.assertEqual(mavedb_snvs.regions({'2': {42}}), [('2', 42, 42)])

    def test_every_position_is_inside_some_region(self):
        rnd = random.Random(3)
        pos = {'7': {rnd.randint(1, 100000) for _ in range(500)}}
        regs = mavedb_snvs.regions(pos, gap=150)
        for p in pos['7']:
            self.assertTrue(any(c == '7' and s <= p <= e for c, s, e in regs))
        starts = [s for _, s, _ in regs]
        self.assertEqual(starts, sorted(starts))


class ContigOrder(unittest.TestCase):
    def test_standard_sequences_in_karyotype_order_then_other_contigs_by_name(self):
        names = ['HSCHR6_MHC_COX_CTG1', 'MT', 'X', '10', '2', 'HG1007_PATCH', '1', 'Y', '22']
        self.assertEqual(sorted(names, key=mavedb_snvs.contig_rank), ['1', '2', '10', '22', 'X', 'Y', 'MT', 'HG1007_PATCH', 'HSCHR6_MHC_COX_CTG1'])

    def test_standard_list_has_25_sequences(self):
        self.assertEqual(len(mavedb_snvs.STANDARD), 25)
        self.assertEqual(len(set(mavedb_snvs.STANDARD)), 25)


class CodonHelpers(unittest.TestCase):
    """The genetic-code table and the genomic-to-codon helper behind the AlphaMissense consistency test, on synthetic genes."""
    class FakeFasta:
        def __init__(self, seqs):
            self.seqs = seqs

        def fetch(self, chrom, start, end):
            return self.seqs[chrom][start:end]

    def test_standard_genetic_code(self):
        c = compat_reference.CODON
        self.assertEqual(len(c), 64)
        self.assertEqual(sorted(k for k, v in c.items() if v == '*'), ['TAA', 'TAG', 'TGA'])
        self.assertEqual((c['ATG'], c['TGG'], c['GGG'], c['AAA'], c['CGG'], c['GAA']), ('M', 'W', 'G', 'K', 'R', 'E'))

    def test_reverse_complement(self):
        self.assertEqual(compat_reference.rc('ATGN'), 'NCAT')

    def test_plus_strand_gene_with_two_coding_exons(self):
        g = list('N' * 50)
        g[10:16] = 'ATGAAA'
        g[30:36] = 'TGGTAA'                                     # spliced CDS ATG AAA TGG TAA = M K W stop
        fa = self.FakeFasta({'c': ''.join(g)})
        pieces = [(11, 16), (31, 36)]
        f = lambda pos, ref, alt: compat_reference.implied_change(fa, 'c', '+', pieces, pos, ref, alt, {}, 'tx')
        self.assertEqual(f(14, 'A', 'G'), (2, 'K', 'E', True))        # AAA -> GAA
        self.assertEqual(f(31, 'T', 'C'), (3, 'W', 'R', True))        # first base of the exon-2 codon: TGG -> CGG
        self.assertEqual(f(14, 'C', 'G'), (2, 'K', 'E', False))       # a wrong REF base is reported, not hidden
        self.assertIsNone(f(20, 'N', 'A'))                            # intron

    def test_minus_strand_gene_is_read_from_the_high_coordinates(self):
        g = list('N' * 50)
        g[30:36] = 'TTTCAT'                                    # reverse complement of ATGAAA, first exon in transcription order
        g[10:16] = 'TTACCA'                                    # reverse complement of TGGTAA
        fa = self.FakeFasta({'c': ''.join(g)})
        pieces = [(11, 16), (31, 36)]
        f = lambda pos, ref, alt: compat_reference.implied_change(fa, 'c', '-', pieces, pos, ref, alt, {}, 'tx')
        self.assertEqual(f(33, 'T', 'C'), (2, 'K', 'E', True))        # + strand T>C is A>G on the transcript: AAA -> GAA
        self.assertEqual(f(16, 'A', 'G'), (3, 'W', 'R', True))        # the highest base of exon 2 is the first base of codon 3 (TGG): A>G on + is T>C on the transcript, CGG
        self.assertIsNone(f(25, 'N', 'A'))


class GenomicHgvs(unittest.TestCase):
    def test_snv_pattern(self):
        m = mavedb_snvs.HGVS_G.match('NC_000017.11:g.43044295A>G')
        self.assertEqual(m.groups(), ('NC_000017', '43044295', 'A', 'G'))
        self.assertEqual(mavedb_snvs.NC2CHR['NC_000017'], '17')
        self.assertEqual(mavedb_snvs.NC2CHR['NC_000023'], 'X')
        self.assertEqual(mavedb_snvs.NC2CHR['NC_000024'], 'Y')

    def test_non_snv_hgvs_is_rejected(self):
        for bad in ('NC_000017.11:g.43044295_43044296insA', 'NC_000017.11:g.43044295del', 'NP_003336.1:p.Pro73Gln', 'NC_000017.11:g.43044295AA>G'):
            self.assertIsNone(mavedb_snvs.HGVS_G.match(bad), bad)


class GnomadInfoParsing(unittest.TestCase):
    def setUp(self):
        # fetch_gnomad_regional imports acquire and mavedb_snvs but not pysam at module level
        import fetch_gnomad_regional as f
        self.f = f

    def test_info_dict_flags_and_values(self):
        d = self.f.info_dict('AC=3;AN=100;lcr;AF=0.03;vep=a=b')
        self.assertEqual(d['AC'], '3')
        self.assertIs(d['lcr'], True)
        self.assertEqual(d['vep'], 'a=b')                 # only the first '=' separates key and value

    def test_pick_per_allele_and_shared_values(self):
        self.assertEqual(self.f.pick('1,2', 1, 2), '2')              # one value per ALT
        self.assertEqual(self.f.pick('7', 1, 2), '7')                # site-level value shared by both alleles
        self.assertEqual(self.f.pick('0.5', 0, 1), '0.5')
        self.assertEqual(self.f.pick(True, 0, 1), 'true')            # INFO flag
        self.assertEqual(self.f.pick(None, 0, 1), '')                # field absent


class AlphaMissenseBands(unittest.TestCase):
    def test_band_function_matches_the_documented_thresholds(self):
        def band(x):
            return 'likely_benign' if x < 0.34 else ('likely_pathogenic' if x > 0.564 else 'ambiguous')
        self.assertEqual(band(0.339), 'likely_benign')
        self.assertEqual(band(0.34), 'ambiguous')
        self.assertEqual(band(0.564), 'ambiguous')
        self.assertEqual(band(0.5641), 'likely_pathogenic')


if __name__ == '__main__':
    unittest.main(verbosity=2)
