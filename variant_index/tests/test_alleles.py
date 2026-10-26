import unittest

import helpers  # noqa: F401  (puts the package on sys.path)
from vindex.build.alleles import resolve_snv, effect_is_alt_from_frequency, frequency_side


class TestOrientation(unittest.TestCase):
    def test_plus_strand_alt_allele(self):
        r = resolve_snv('C', 'T', 0.3, {'T'})
        self.assertEqual((r.status, r.alt, r.effect, r.other, r.strand), ('mapped', 'T', 'T', 'C', '+'))

    def test_effect_allele_is_reference(self):
        r = resolve_snv('C', 'C', 0.7, {'T'})
        self.assertEqual((r.status, r.alt, r.effect, r.other), ('mapped', 'T', 'C', 'T'))
        self.assertIn('effect_allele_is_reference', r.flags)

    def test_minus_strand_is_flipped(self):
        # reference C, alternate T known at the site; the study reported G (the complement of C) -> effect allele is the reference on the + strand
        r = resolve_snv('C', 'G', 0.7, {'T'})
        self.assertEqual((r.status, r.effect, r.alt, r.strand), ('mapped', 'C', 'T', '-'))
        r2 = resolve_snv('C', 'A', 0.3, {'T'})        # A is the complement of T: effect allele is the alternate T
        self.assertEqual((r2.status, r2.effect, r2.alt, r2.strand), ('mapped', 'T', 'T', '-'))

    def test_palindromic_site_needs_frequency(self):
        self.assertEqual(resolve_snv('A', 'A', None, {'T'}, {'T': [0.2, 0.25]}).reason, 'strand_unresolved_palindromic')
        # alt T is minor everywhere (0.2); the reported frequency 0.8 is that of the REF allele A -> effect allele is A on the + strand
        r = resolve_snv('A', 'A', 0.8, {'T'}, {'T': [0.2, 0.25]})
        self.assertEqual((r.status, r.effect, r.alt, r.strand), ('mapped', 'A', 'T', '+'))
        # reported A with frequency 0.2: the study counted the minor allele, so its A is the forward T (minus strand)
        r = resolve_snv('A', 'A', 0.2, {'T'}, {'T': [0.2, 0.25]})
        self.assertEqual((r.status, r.effect, r.alt, r.strand), ('mapped', 'T', 'T', '-'))
        self.assertIn('strand_resolved_by_frequency', r.flags)

    def test_frequency_near_half_or_ancestry_dependent_is_rejected(self):
        r = resolve_snv('A', 'A', 0.5, {'T'}, {'T': [0.45, 0.5]})
        self.assertEqual((r.status, r.reason), ('rejected', 'strand_unresolved_palindromic'))
        r = resolve_snv('A', 'A', 0.8, {'T'}, {'T': [0.2, 0.7]})          # alt minor in one population, major in another
        self.assertEqual(r.status, 'rejected')
        r = resolve_snv('A', 'A', 0.8, {'T'}, None)
        self.assertEqual(r.status, 'rejected')

    def test_allele_not_at_site(self):
        self.assertEqual(resolve_snv('C', 'G', 0.5, {'T'}).status, 'mapped')   # G is the complement of the reference: minus strand
        self.assertEqual(resolve_snv('T', 'G', 0.5, {'A'}).reason, 'effect_allele_not_among_site_alleles')
        r = resolve_snv('C', 'A', 0.5, {'T'})
        self.assertEqual(r.status, 'mapped')      # complement of A is T, a known alternate: minus strand
        r = resolve_snv('C', 'G', 0.5, {'A'})     # G is the complement of the reference C: effect allele = reference; alt A
        self.assertEqual((r.status, r.effect, r.strand), ('mapped', 'C', '-'))
        r = resolve_snv('C', 'T', 0.5, {'A'})     # T and its complement A... complement A is a known alt -> minus
        self.assertEqual((r.status, r.effect), ('mapped', 'A'))
        r = resolve_snv('C', 'C', 0.5, {'A', 'T'})
        self.assertEqual((r.status, r.reason), ('ambiguous', 'effect_allele_is_reference_several_alternate_alleles'))
        r = resolve_snv('C', 'C', 0.5, {'A', 'G'})      # C is the reference and its complement G is an alternate: both strands possible, several alternates
        self.assertEqual((r.status, r.reason), ('rejected', 'strand_unresolved_multiallelic'))

    def test_unknown_or_missing_information(self):
        self.assertEqual(resolve_snv('C', 'T', 0.3, set()).reason, 'no_alternate_allele_known')
        self.assertEqual(resolve_snv('C', '?', None, {'T'}).flags, ('effect_allele_unknown',))
        self.assertEqual(resolve_snv('C', '?', None, {'T', 'G'}).status, 'ambiguous')
        self.assertEqual(resolve_snv('N', 'T', None, {'T'}).status, 'unmapped')
        self.assertEqual(resolve_snv('C', 'TT', None, {'T'}).reason, 'effect_allele_not_a_single_base')

    def test_frequency_helpers(self):
        self.assertEqual(frequency_side([0.1, 0.3]), 'minor')
        self.assertEqual(frequency_side([0.7, 0.9]), 'major')
        self.assertIsNone(frequency_side([0.3, 0.7]))
        self.assertTrue(effect_is_alt_from_frequency(0.2, [0.1, 0.3]))
        self.assertFalse(effect_is_alt_from_frequency(0.8, [0.1, 0.3]))
        self.assertIsNone(effect_is_alt_from_frequency(0.5, [0.1, 0.3]))


if __name__ == '__main__':
    unittest.main()
