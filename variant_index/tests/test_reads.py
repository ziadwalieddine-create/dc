"""Read-level evidence at a locus: allele counting from alignments independent of how the aligner wrote the event, the filters, the caveats that stop absence being read as
reference, and the refusals. The alignments are synthetic (tests/bamfix.py), cut from the synthetic reference."""
import os
import shutil
import tempfile
import unittest

import pysam

import bamfix
from bamfix import make_read, span
from helpers import SYN_CONTIGS, SYN_FASTA, syn_reference
from vindex import reads as R
from vindex.reads import ReadEvidence, ReadPolicy
from vindex.refgenome import Reference

SUBST = {'A': 'G', 'C': 'T', 'G': 'A', 'T': 'C'}
OTHER_SUBST = {'A': 'C', 'C': 'A', 'G': 'T', 'T': 'G'}


class ReadsCase(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.ref = syn_reference()
        cls.fa = pysam.FastaFile(SYN_FASTA)
        cls.dir = tempfile.mkdtemp()
        cls.count = 0

    @classmethod
    def tearDownClass(cls):
        shutil.rmtree(cls.dir, ignore_errors=True)

    def base(self, pos, chrom='1'):
        return self.fa.fetch(chrom, pos - 1, pos)

    def evaluate(self, reads, key, policy=None, contigs=SYN_CONTIGS, rename=None):
        type(self).count += 1
        path = bamfix.write_bam(os.path.join(self.dir, f'r{self.count}.bam'), reads, contigs, rename)
        with ReadEvidence(path, self.ref, policy) as e:
            return e.evaluate(key)

    # a SNV at pos with reference base b and alternate SUBST[b]
    def snv_key(self, pos):
        b = self.base(pos)
        return f'1:{pos}:{b}:{SUBST[b]}'

    def snv_reads(self, pos, n_ref=0, n_alt=0, n_other=0, prefix='r', **kw):
        """Fragments of single reads of 40 bases around pos, starts spread over 20 offsets and strands alternating."""
        b, out, i = self.base(pos), [], 0
        for allele, n in (('ref', n_ref), ('alt', n_alt), ('other', n_other)):
            for _ in range(n):
                start = pos - 35 + (i % 20)
                edits = {} if allele == 'ref' else {pos: SUBST[b]} if allele == 'alt' else {pos: OTHER_SUBST[b]}
                flag = kw.get('flag', 0) | (bamfix.FLAG_REVERSE if i % 2 else 0)
                out.append(span(self.fa, '1', start, start + 39, f'{prefix}{allele}{i}', edits, **{k: v for k, v in kw.items() if k != 'flag'}, flag=flag))
                i += 1
        return out


class TestSnv(ReadsCase):
    def test_no_alternate_fragment_at_depth_is_reference_only(self):
        s = self.evaluate(self.snv_reads(405, n_ref=12), self.snv_key(405))
        self.assertEqual((s.state, s.depth, s.ref_reads, s.alt_reads, s.other_reads, s.alt_fraction), (R.REFERENCE_ONLY, 12, 12, 0, 0, 0.0))

    def test_heterozygous_and_homozygous_alternate_are_present_with_their_counts(self):
        het = self.evaluate(self.snv_reads(405, n_ref=7, n_alt=6), self.snv_key(405))
        self.assertEqual((het.state, het.depth, het.ref_reads, het.alt_reads, het.alt_fraction), (R.ALTERNATE_PRESENT, 13, 7, 6, round(6 / 13, 6)))
        hom = self.evaluate(self.snv_reads(405, n_alt=12), self.snv_key(405))
        self.assertEqual((hom.state, hom.ref_reads, hom.alt_reads, hom.alt_fraction), (R.ALTERNATE_PRESENT, 0, 12, 1.0))

    def test_the_depth_threshold_is_exact(self):
        self.assertEqual(self.evaluate(self.snv_reads(405, n_ref=10), self.snv_key(405)).state, R.REFERENCE_ONLY)
        s = self.evaluate(self.snv_reads(405, n_ref=9), self.snv_key(405))
        self.assertEqual((s.state, s.reasons), (R.INSUFFICIENT, ('fewer_informative_fragments_than_min_depth',)))

    def test_no_reads_is_insufficient_and_says_so(self):
        s = self.evaluate(self.snv_reads(900, n_ref=12), self.snv_key(405))                           # reads elsewhere on the contig
        self.assertEqual((s.state, s.reasons, s.depth), (R.INSUFFICIENT, ('no_overlapping_alignments',), 0))

    def test_alternate_fragments_below_the_thresholds_are_ambiguous_not_reference(self):
        for n_ref, n_alt in ((19, 1), (28, 2), (37, 3)):                                                # 5%, 7%, 7.5%: under both the count and the fraction, or the fraction
            s = self.evaluate(self.snv_reads(405, n_ref=n_ref, n_alt=n_alt), self.snv_key(405))
            self.assertEqual((s.state, s.reasons, s.alt_reads), (R.AMBIGUOUS, ('alternate_fragments_below_the_thresholds',), n_alt), (n_ref, n_alt))
        s = self.evaluate(self.snv_reads(405, n_ref=9, n_alt=3), self.snv_key(405))                     # 3 of 12 = 25%: present
        self.assertEqual(s.state, R.ALTERNATE_PRESENT)
        s = self.evaluate(self.snv_reads(405, n_ref=8, n_alt=2), self.snv_key(405))                     # 20% of 10, but only two fragments: not enough to say present
        self.assertEqual((s.state, s.reasons), (R.AMBIGUOUS, ('alternate_fragments_below_the_thresholds',)))
        s = self.evaluate(self.snv_reads(405, n_ref=6, n_alt=2), self.snv_key(405))                     # 25% of 8, two fragments, too few fragments to read the rest as reference
        self.assertEqual(s.state, R.INSUFFICIENT)

    def test_a_few_alternate_fragments_at_low_depth_are_present_with_a_flag(self):
        s = self.evaluate(self.snv_reads(405, n_ref=2, n_alt=4), self.snv_key(405))
        self.assertEqual((s.state, s.depth), (R.ALTERNATE_PRESENT, 6))
        self.assertIn('fewer_informative_fragments_than_min_depth', s.flags)

    def test_another_allele_at_the_site_is_not_the_indexed_one(self):
        b = self.base(405)
        s = self.evaluate(self.snv_reads(405, n_ref=8, n_other=8), self.snv_key(405))
        self.assertEqual((s.state, s.alt_reads, s.other_reads), (R.OTHER_ALLELE_PRESENT, 0, 8))
        self.assertEqual(s.detail['other_alleles'], {OTHER_SUBST[b]: 8})
        key = f'1:405:{b}:{OTHER_SUBST[b]}'                                                              # asked for the allele the reads do carry
        self.assertEqual(self.evaluate(self.snv_reads(405, n_ref=8, n_other=8), key).state, R.ALTERNATE_PRESENT)

    def test_flagged_alignments_never_count(self):
        for flag, name in ((bamfix.FLAG_DUPLICATE, 'duplicate'), (bamfix.FLAG_SECONDARY, 'secondary'), (bamfix.FLAG_SUPPLEMENTARY, 'supplementary'), (bamfix.FLAG_QCFAIL, 'qc_fail')):
            reads = self.snv_reads(405, n_ref=12) + self.snv_reads(405, n_alt=10, prefix='bad', flag=flag)
            s = self.evaluate(reads, self.snv_key(405))
            self.assertEqual((s.state, s.alt_reads, s.detail['filtered'].get(name)), (R.REFERENCE_ONLY, 0, 10), name)

    def test_unmapped_reads_placed_at_the_locus_never_count(self):
        placed = [bamfix.Read(f'u{i}', '1', 405, 'ACGT' * 5, '', flag=bamfix.FLAG_UNMAPPED, mapq=0) for i in range(5)]
        s = self.evaluate(self.snv_reads(405, n_ref=12) + placed, self.snv_key(405))
        self.assertEqual((s.state, s.detail['filtered'].get('unmapped')), (R.REFERENCE_ONLY, 5))

    def test_alternate_reads_below_the_mapping_quality_are_counted_and_stop_a_reference_call(self):
        reads = self.snv_reads(405, n_ref=12) + self.snv_reads(405, n_alt=4, prefix='low', mapq=0)
        s = self.evaluate(reads, self.snv_key(405))
        self.assertEqual((s.state, s.alt_reads, s.detail['low_mapq_by_allele']), (R.AMBIGUOUS, 0, {'alt': 4}))
        self.assertEqual(s.reasons, ('alternate_like_reads_below_the_mapping_quality_threshold',))

    def test_low_mapping_quality_that_dominates_forbids_reading_absence_as_reference(self):
        s = self.evaluate(self.snv_reads(405, n_ref=12) + self.snv_reads(405, n_ref=6, prefix='low', mapq=5), self.snv_key(405))              # 6 of 18 = 33%
        self.assertEqual((s.state, s.reasons), (R.AMBIGUOUS, ('low_mapping_quality_reads_dominate',)))
        s = self.evaluate(self.snv_reads(405, n_ref=12) + self.snv_reads(405, n_ref=2, prefix='low', mapq=5), self.snv_key(405))              # 2 of 14 = 14%
        self.assertEqual(s.state, R.REFERENCE_ONLY)

    def test_reads_that_are_all_low_mapping_quality_give_no_depth(self):
        s = self.evaluate(self.snv_reads(405, n_ref=30, mapq=0), self.snv_key(405))
        self.assertEqual((s.state, s.depth, s.reasons), (R.INSUFFICIENT, 0, ('no_informative_fragments',)))

    def test_base_quality_of_the_allele_base_decides_and_alternate_like_reads_below_it_are_counted(self):
        s = self.evaluate(self.snv_reads(405, n_ref=12) + self.snv_reads(405, n_alt=4, prefix='lq', quals=10), self.snv_key(405))
        self.assertEqual((s.state, s.alt_reads, s.detail['low_base_quality_by_allele']), (R.AMBIGUOUS, 0, {'alt': 4}))
        self.assertEqual(s.reasons, ('alternate_like_reads_below_the_base_quality_threshold',))
        s = self.evaluate(self.snv_reads(405, n_ref=12) + self.snv_reads(405, n_ref=3, prefix='lq', quals=10), self.snv_key(405))
        self.assertEqual((s.state, s.depth), (R.REFERENCE_ONLY, 12))                                      # low-quality reference-like reads are only left out
        s = self.evaluate(self.snv_reads(405, n_ref=12) + self.snv_reads(405, n_alt=4, prefix='lq', quals=10), self.snv_key(405), ReadPolicy(min_baseq=5))
        self.assertEqual((s.state, s.alt_reads), (R.ALTERNATE_PRESENT, 4))                                # the threshold is the policy's

    def test_missing_base_qualities_fail_the_base_quality_policy(self):
        reads = self.snv_reads(405, n_ref=12, quals=None)
        s = self.evaluate(reads, self.snv_key(405))
        self.assertEqual((s.state, s.depth, s.detail['filtered'].get('low_base_quality')), (R.INSUFFICIENT, 0, 12))
        self.assertEqual(s.reasons, ('no_informative_fragments', 'most_alignments_fail_the_base_quality_policy'))
        self.assertEqual(self.evaluate(reads, self.snv_key(405), ReadPolicy(min_baseq=0)).state, R.REFERENCE_ONLY)

    def test_both_mates_of_a_pair_are_one_fragment(self):
        b, reads = self.base(405), []
        for i in range(12):
            alt = {405: SUBST[b]} if i < 5 else {}
            reads.append(span(self.fa, '1', 380 + i % 5, 415, f'p{i}', alt, flag=bamfix.FLAG_PAIRED | bamfix.FLAG_PROPER | bamfix.FLAG_READ1 | bamfix.FLAG_MATE_REVERSE))
            reads.append(span(self.fa, '1', 395, 430 + i % 5, f'p{i}', alt, flag=bamfix.FLAG_PAIRED | bamfix.FLAG_PROPER | bamfix.FLAG_READ2 | bamfix.FLAG_REVERSE))
        s = self.evaluate(reads, self.snv_key(405))
        self.assertEqual((s.state, s.depth, s.ref_reads, s.alt_reads, s.detail['informative_alignments']), (R.ALTERNATE_PRESENT, 12, 7, 5, 24))
        self.assertEqual(s.detail['strand'], {'alt_forward': 5, 'ref_forward': 7})                       # the strand of the first mate of a fragment

    def test_mates_that_disagree_are_set_aside_and_forbid_a_reference_call_when_many(self):
        b = self.base(405)

        def pair(i, disagree):
            f1 = bamfix.FLAG_PAIRED | bamfix.FLAG_PROPER | bamfix.FLAG_READ1
            f2 = bamfix.FLAG_PAIRED | bamfix.FLAG_PROPER | bamfix.FLAG_READ2 | bamfix.FLAG_REVERSE
            return [span(self.fa, '1', 385, 420, f'q{i}', {}, flag=f1), span(self.fa, '1', 390, 425, f'q{i}', {405: SUBST[b]} if disagree else {}, flag=f2)]
        reads = self.snv_reads(405, n_ref=12) + [r for i in range(5) for r in pair(i, True)]
        s = self.evaluate(reads, self.snv_key(405))
        self.assertEqual((s.state, s.depth, s.detail['conflicting_fragments'], s.reasons), (R.AMBIGUOUS, 12, 5, ('overlapping_mates_disagree',)))
        reads = self.snv_reads(405, n_ref=12) + [r for i in range(1) for r in pair(i, True)]
        s = self.evaluate(reads, self.snv_key(405))
        self.assertEqual((s.state, s.detail['conflicting_fragments']), (R.REFERENCE_ONLY, 1))

    def test_alternate_fragments_on_one_strand_only_are_flagged_not_dismissed(self):
        b = self.base(405)
        reads = self.snv_reads(405, n_ref=10) + [span(self.fa, '1', 385 + i, 424 + i, f'f{i}', {405: SUBST[b]}) for i in range(6)]
        s = self.evaluate(reads, self.snv_key(405))
        self.assertEqual((s.state, s.flags), (R.ALTERNATE_PRESENT, ('alternate_fragments_on_one_strand_only',)))

    def test_a_deletion_over_the_site_in_some_reads_forbids_a_reference_call(self):
        gapped = [make_read(self.fa, '1', 380, [('M', 25), ('D', 6), ('M', 20)], f'g{i}') for i in range(6)]      # 380..404, 405..410 deleted, 411..430
        s = self.evaluate(self.snv_reads(405, n_ref=12) + gapped, self.snv_key(405))
        self.assertEqual((s.state, s.depth, s.reasons, s.detail['filtered']['window_end_in_a_gap']), (R.AMBIGUOUS, 12, ('a_deletion_covers_the_locus_in_some_reads',), 6))

    def test_alignments_that_do_not_cover_the_site_are_not_evidence_of_it(self):
        short = [span(self.fa, '1', 395 + i, 404, f's{i}') for i in range(5)]                          # end before the site
        s = self.evaluate(self.snv_reads(405, n_ref=12) + short, self.snv_key(405))
        self.assertEqual((s.state, s.depth), (R.REFERENCE_ONLY, 12))

    def test_read_order_and_names_do_not_change_the_answer(self):
        reads = self.snv_reads(405, n_ref=7, n_alt=6)
        a = self.evaluate(reads, self.snv_key(405))
        b = self.evaluate(list(reversed(reads)), self.snv_key(405))
        self.assertEqual(a.columns('ABSENT_FROM_CALLSET', ()), b.columns('ABSENT_FROM_CALLSET', ()))

    def test_more_alignments_than_the_policy_allows_is_not_evaluated(self):
        s = self.evaluate(self.snv_reads(405, n_ref=12), self.snv_key(405), ReadPolicy(max_reads=5))
        self.assertEqual((s.state, s.reasons), (R.NOT_EVALUATED, ('more_alignments_than_the_read_policy_allows',)))


def equivalent_deletions(fa, chrom, pos, ref):
    """By brute force over the sequence: every start of a deleted span that leaves the same sequence as deleting ref[1:] after the anchor at pos."""
    k, lo, hi = len(ref) - 1, pos - 40, pos + 80
    seq = fa.fetch(chrom, lo - 1, hi)
    target = seq[:pos - lo + 1] + seq[pos - lo + 1 + k:]
    return [s for s in range(pos - 5, pos + 60) if seq[:s - lo] + seq[s - lo + k:] == target]


def equivalent_insertions(fa, chrom, pos, alt):
    """By brute force: every position q such that inserting some string of the same length after q gives the sequence of inserting alt[1:] after pos; (q, inserted string)."""
    ins, lo, hi = alt[1:], pos - 40, pos + 80
    seq = fa.fetch(chrom, lo - 1, hi)
    target = seq[:pos - lo + 1] + ins + seq[pos - lo + 1:]
    out = []
    for q in range(pos - 5, pos + 60):
        x = target[q - lo + 1:q - lo + 1 + len(ins)]
        if seq[:q - lo + 1] + x + seq[q - lo + 1:] == target:
            out.append((q, x))
    return out


class TestIndels(ReadsCase):
    HOMOPOLYMER_DEL = '1:101:CA:C'            # one A of the A8 at 102..109, left-aligned
    HOMOPOLYMER_INS = '1:101:C:CA'
    DINUC_DEL = '1:300:CGT:C'                 # the GT unit of the (GT) repeat at 301..311, left-aligned

    def deletion_read(self, chrom, s, k, i, **kw):
        start = s - 25
        return make_read(self.fa, chrom, start, [('M', s - start), ('D', k), ('M', 40)], f'd{i}', **kw)

    def insertion_read(self, chrom, q, x, i, **kw):
        start = q - 24
        return make_read(self.fa, chrom, start, [('M', q - start + 1), ('I', x), ('M', 40)], f'i{i}', **kw)

    def ref_reads(self, pos, n, prefix='ref'):
        return [span(self.fa, '1', pos - 30 + i % 10, pos + 49, f'{prefix}{i}', flag=bamfix.FLAG_REVERSE if i % 2 else 0) for i in range(n)]

    def test_the_equivalent_placements_the_test_enumerates_include_the_right_end_of_the_repeat(self):
        self.assertEqual(equivalent_deletions(self.fa, '1', 101, 'CA'), list(range(102, 110)))
        self.assertEqual([q for q, _ in equivalent_insertions(self.fa, '1', 101, 'CA')], list(range(101, 110)))
        self.assertEqual(equivalent_deletions(self.fa, '1', 300, 'CGT')[0], 301)

    def test_a_deletion_is_counted_wherever_the_aligner_put_the_gap(self):
        placements = equivalent_deletions(self.fa, '1', 101, 'CA')
        self.assertGreater(len(placements), 3)
        reads = [self.deletion_read('1', s, 1, i) for i, s in enumerate(placements)] + self.ref_reads(101, 8)
        s = self.evaluate(reads, self.HOMOPOLYMER_DEL)
        self.assertEqual((s.state, s.alt_reads, s.ref_reads, s.other_reads), (R.ALTERNATE_PRESENT, len(placements), 8, 0))

    def test_an_insertion_is_counted_wherever_the_aligner_put_the_gap(self):
        placements = equivalent_insertions(self.fa, '1', 101, 'CA')
        reads = [self.insertion_read('1', q, x, i) for i, (q, x) in enumerate(placements)] + self.ref_reads(101, 8)
        s = self.evaluate(reads, self.HOMOPOLYMER_INS)
        self.assertEqual((s.state, s.alt_reads, s.ref_reads, s.other_reads), (R.ALTERNATE_PRESENT, len(placements), 8, 0))

    def test_a_dinucleotide_deletion_is_counted_at_the_right_end_of_its_repeat(self):
        placements = equivalent_deletions(self.fa, '1', 300, 'CGT')
        self.assertEqual(placements, list(range(301, 311)))                                                 # GTGTGTGTGTG at 301..311: a two-base span anywhere inside it
        reads = [self.deletion_read('1', s, 2, i) for i, s in enumerate(placements)] + self.ref_reads(300, 8)
        s = self.evaluate(reads, self.DINUC_DEL)
        self.assertEqual((s.state, s.alt_reads, s.other_reads), (R.ALTERNATE_PRESENT, len(placements), 0))

    def test_a_tandem_repeat_insertion_is_counted_at_any_unit_boundary(self):
        from vindex.normalise import normalise
        n = normalise(self.ref, '1', 215, 'G', 'GCAG')                                                      # an extra CAG written at the right end of the (CAG)5 repeat at 201..215
        self.assertEqual((n.pos, n.ref, n.alt), (200, 'T', 'TCAG'))
        placements = equivalent_insertions(self.fa, '1', n.pos, n.alt)
        self.assertGreaterEqual(len(placements), 3)
        reads = [self.insertion_read('1', q, x, i) for i, (q, x) in enumerate(placements)] + self.ref_reads(200, 8)
        s = self.evaluate(reads, n.key)
        self.assertEqual((s.state, s.alt_reads, s.other_reads), (R.ALTERNATE_PRESENT, len(placements), 0))

    def test_the_other_length_of_the_same_event_is_another_allele(self):
        reads = [self.deletion_read('1', 102, 2, i) for i in range(8)] + self.ref_reads(101, 8)           # two A deleted, not one
        s = self.evaluate(reads, self.HOMOPOLYMER_DEL)
        self.assertEqual((s.state, s.alt_reads, s.other_reads), (R.OTHER_ALLELE_PRESENT, 0, 8))
        wrong = [self.insertion_read('1', 101, 'G', i) for i in range(8)] + self.ref_reads(101, 8)        # an inserted base that is not A
        self.assertEqual(self.evaluate(wrong, self.HOMOPOLYMER_INS).state, R.OTHER_ALLELE_PRESENT)

    def test_the_indexed_insertion_with_another_variant_of_the_same_haplotype_in_the_window_is_another_sequence(self):
        # the limit of comparing whole windows, seen in real reads (the benchmark lists an insertion and a neighbouring SNV of one haplotype as two records): the reads carry the insertion,
        # but their window is neither the reference's nor the insertion's alone, so they are counted as another sequence and kept, never as the indexed allele and never as reference
        b = self.base(111)                                                                                  # a base inside the window, after the A run
        hap = [make_read(self.fa, '1', 77, [('M', 25), ('I', 'A'), ('M', 9), ('X', SUBST[b]), ('M', 30)], f'h{i}') for i in range(8)]           # 77..101, an A, 102..110, 111 changed
        s = self.evaluate(hap + self.ref_reads(101, 8), self.HOMOPOLYMER_INS)
        self.assertEqual((s.state, s.alt_reads, s.other_reads), (R.OTHER_ALLELE_PRESENT, 0, 8))
        (seq, n), = s.detail['other_alleles'].items()
        self.assertEqual(n, 8)
        self.assertEqual(len(seq), len(R.plan_locus(self.ref, self.HOMOPOLYMER_INS, ReadPolicy())[0].alt_window))     # the insertion is in it: one base longer than the reference's window
        self.assertEqual(R.interpretation('ABSENT_FROM_CALLSET', s.state), 'reads_show_another_sequence_over_the_window')

    def test_a_deletion_with_a_mismatch_inside_the_window_is_not_the_indexed_allele(self):
        b = self.base(110)
        reads = [make_read(self.fa, '1', 76, [('M', 26), ('D', 1), ('M', 7), ('X', SUBST[b]), ('M', 30)], f'm{i}') for i in range(8)] + self.ref_reads(101, 8)     # 76..101, 102 deleted, 103..109, 110 changed
        s = self.evaluate(reads, self.HOMOPOLYMER_DEL)
        self.assertEqual((s.state, s.alt_reads, s.other_reads), (R.OTHER_ALLELE_PRESENT, 0, 8))

    def test_reads_that_do_not_span_the_repeat_cannot_show_its_length(self):
        short = [span(self.fa, '1', 90 + i, 106, f's{i}') for i in range(8)]                              # end inside the A run
        s = self.evaluate(short + self.ref_reads(101, 4), self.HOMOPOLYMER_DEL)
        self.assertEqual((s.state, s.depth, s.reasons), (R.INSUFFICIENT, 4, ('fewer_informative_fragments_than_min_depth', 'most_alignments_do_not_span_the_window')))
        self.assertEqual(s.detail['filtered']['not_spanning_the_window'], 8)

    def test_soft_clipped_reads_ending_at_the_site_forbid_a_reference_call_when_many(self):
        clipped = [make_read(self.fa, '1', 70, [('M', 32), ('S', 'A' * 25)], f'c{i}') for i in range(8)]      # aligned to 101, the rest clipped: the place of the event
        s = self.evaluate(self.ref_reads(101, 12) + clipped, self.HOMOPOLYMER_INS)
        self.assertEqual((s.state, s.reasons, s.detail['soft_clipped_ending_at_the_locus']), (R.AMBIGUOUS, ('soft_clipped_reads_end_at_the_locus',), 8))
        s = self.evaluate(self.ref_reads(101, 12) + clipped[:2], self.HOMOPOLYMER_INS)
        self.assertEqual((s.state, s.detail['soft_clipped_ending_at_the_locus']), (R.REFERENCE_ONLY, 2))
        left = [make_read(self.fa, '1', 108, [('S', 'C' * 25), ('M', 40)], f'l{i}') for i in range(8)]       # clipped on the left, starting inside the window
        s = self.evaluate(self.ref_reads(101, 12) + left, self.HOMOPOLYMER_INS)
        self.assertEqual(s.reasons, ('soft_clipped_reads_end_at_the_locus',))

    def test_clipped_reads_far_from_the_site_are_not_a_caveat(self):
        far = [make_read(self.fa, '1', 20, [('M', 50), ('S', 'A' * 25)], f'c{i}') for i in range(8)]
        s = self.evaluate(self.ref_reads(101, 12) + far, self.HOMOPOLYMER_INS)
        self.assertEqual(s.state, R.REFERENCE_ONLY)

    def test_mnv_components_in_the_same_read_are_the_mnv_and_one_component_is_not(self):
        b = [self.base(p) for p in (504, 505)]
        key = f'1:504:{b[0]}{b[1]}:{SUBST[b[0]]}{SUBST[b[1]]}'
        both = [span(self.fa, '1', 480 + i, 520 + i, f'b{i}', {504: SUBST[b[0]], 505: SUBST[b[1]]}) for i in range(5)]
        first = [span(self.fa, '1', 480 + i, 520 + i, f'f{i}', {504: SUBST[b[0]]}) for i in range(5)]
        second = [span(self.fa, '1', 480 + i, 520 + i, f's{i}', {505: SUBST[b[1]]}) for i in range(5)]
        ref = [span(self.fa, '1', 480 + i, 520 + i, f'r{i}') for i in range(8)]
        s = self.evaluate(both + ref, key)
        self.assertEqual((s.state, s.alt_reads, s.ref_reads), (R.ALTERNATE_PRESENT, 5, 8))
        trans = self.evaluate(first + second + ref, key)                                                    # the two alternate bases lie on different reads: in trans
        self.assertEqual((trans.state, trans.alt_reads, trans.other_reads), (R.OTHER_ALLELE_PRESENT, 0, 10))
        self.assertEqual(sorted(trans.detail['other_alleles'].values()), [5, 5])


class TestPlanning(ReadsCase):
    def plan(self, key, **kw):
        return R.plan_locus(self.ref, key, ReadPolicy(**kw))

    def test_an_snv_is_compared_on_its_own_base_and_an_indel_over_its_whole_repeat(self):
        locus, why = self.plan('1:405:' + self.base(405) + ':' + SUBST[self.base(405)])
        self.assertEqual((why, locus.start, locus.end, len(locus.ref_window)), (None, 405, 405, 1))
        locus, _ = self.plan('1:101:CA:C')
        self.assertEqual((locus.start, locus.end), (99, 111))
        self.assertEqual((locus.ref_window, locus.alt_window), ('TGCAAAAAAAAGT', 'TGCAAAAAAAGT'))
        locus, _ = self.plan('1:101:C:CA')
        self.assertEqual((locus.start, locus.end), (99, 112))
        self.assertEqual(locus.alt_window, 'TGCAAAAAAAAAGTA')
        locus, _ = self.plan('1:300:CGT:C')
        self.assertEqual((locus.start, locus.end), (298, 313))

    def test_a_key_that_reads_cannot_be_asked_is_not_evaluated_with_the_reason(self):
        b = self.base(405)
        self.assertEqual(self.plan(f'1:405:{"T" if b != "T" else "A"}:G')[1], 'REF_disagrees_with_the_reference')
        self.assertEqual(self.plan('1:405:' + self.base(405) + 'A' * 60 + ':' + self.base(405))[1], 'allele_longer_than_the_read_policy_allows')
        self.assertEqual(self.plan('7:405:A:G')[1], 'contig_not_in_the_reference')
        self.assertIsNone(self.plan('1:2999:' + self.fa.fetch('1', 2998, 3000) + ':A')[1])                    # a REF that ends on the last base is inside the contig
        self.assertEqual(self.plan('1:3000:' + self.base(3000) + 'A:A')[1], 'locus_outside_the_contig')
        self.assertEqual(self.plan('1:101:CA:C', max_window=5)[1], 'repeat_tract_longer_than_the_window_limit')        # the A run goes on past what a window may span
        self.assertEqual(self.plan('1:101:CA:C', max_window=10)[1], 'window_longer_than_the_window_limit')           # the run ends inside the limit but the window around it does not

    def test_n_in_the_window_is_refused(self):
        d = tempfile.mkdtemp()
        try:
            fa = os.path.join(d, 'n.fa')
            with open(fa, 'w') as f:
                f.write('>1\n' + 'ACGT' * 25 + 'NACGT' + 'ACGT' * 25 + '\n')
            pysam.faidx(fa)
            ref = Reference(fa, name='N', builds={})
            self.assertEqual(R.plan_locus(ref, '1:101:N:A', ReadPolicy())[1], 'reference_has_N_in_the_window')
            self.assertEqual(R.plan_locus(ref, '1:99:GT:G', ReadPolicy())[1], 'reference_has_N_in_the_window')        # the flank of a deletion reaches the N
        finally:
            shutil.rmtree(d, ignore_errors=True)

    def test_a_key_near_the_start_of_the_contig_is_clamped_not_refused(self):
        b = self.base(2)
        locus, why = self.plan(f'1:2:{b}:{SUBST[b]}')
        self.assertEqual((why, locus.start), (None, 2))
        locus, why = self.plan('1:1:' + self.base(1) + ':' + SUBST[self.base(1)])
        self.assertEqual((why, locus.start, locus.end), (None, 1, 1))
        pair = self.fa.fetch('1', 1, 3)                                                                     # a deletion anchored at position 2: the flank is cut at the contig start
        locus, why = self.plan(f'1:2:{pair}:{pair[0]}')
        self.assertEqual((why, locus.start), (None, 1))

    def test_policy_values_are_validated(self):
        for bad in (dict(min_depth=0), dict(min_alt_reads=0), dict(min_alt_fraction=0), dict(min_alt_fraction=1.5), dict(flank=0), dict(max_reads=0), dict(min_mapq=-1)):
            with self.assertRaises(ValueError, msg=str(bad)):
                ReadPolicy(**bad)
        self.assertEqual(ReadPolicy().min_depth, 10)


class TestAlignmentAndAssembly(ReadsCase):
    def test_contig_names_in_the_header_are_mapped_to_the_reference(self):
        reads = self.snv_reads(405, n_ref=12)
        s = self.evaluate(reads, self.snv_key(405), rename={'1': 'chr1', '2': 'chr2', 'X': 'chrX', 'MT': 'chrM'})
        self.assertEqual((s.state, s.depth), (R.REFERENCE_ONLY, 12))

    def test_a_header_of_another_assembly_stops_the_run(self):
        contigs = [('1', 3001)] + SYN_CONTIGS[1:]                                                          # chr1 one base longer than the reference
        path = bamfix.write_bam(os.path.join(self.dir, 'wrong.bam'), self.snv_reads(405, n_ref=3), contigs)
        with self.assertRaises(R.ReadsAssemblyMismatch):
            ReadEvidence(path, self.ref)

    def test_a_contig_the_alignment_does_not_carry_is_not_evaluated(self):
        s = self.evaluate(self.snv_reads(405, n_ref=12), '2:105:' + self.base(105, '2') + ':' + SUBST[self.base(105, '2')], contigs=[('1', 3000)])
        self.assertEqual((s.state, s.reasons), (R.NOT_EVALUATED, ('contig_not_in_the_alignment_header',)))

    def test_a_contig_of_another_length_is_never_used_even_when_the_mismatch_is_allowed(self):
        contigs = [('1', 3000), ('2', 1501)] + SYN_CONTIGS[2:]
        path = bamfix.write_bam(os.path.join(self.dir, 'len2.bam'), self.snv_reads(405, n_ref=12), contigs)
        with self.assertRaises(R.ReadsAssemblyMismatch):
            ReadEvidence(path, self.ref)
        with ReadEvidence(path, self.ref, allow_assembly_mismatch=True) as e:
            other = e.evaluate('2:105:' + self.base(105, '2') + ':' + SUBST[self.base(105, '2')])
            same = e.evaluate(self.snv_key(405))
        self.assertEqual((other.state, other.reasons), (R.NOT_EVALUATED, ('contig_not_in_the_alignment_header',)))
        self.assertEqual(same.state, R.REFERENCE_ONLY)

    def test_a_cram_is_read_against_the_reference_and_gives_the_answer_of_the_bam(self):
        reads = self.snv_reads(405, n_ref=7, n_alt=6)
        bam = bamfix.write_bam(os.path.join(self.dir, 'c.bam'), reads, SYN_CONTIGS)
        cram = os.path.join(self.dir, 'c.cram')
        pysam.view('-C', '-T', SYN_FASTA, '-o', cram, bam, catch_stdout=False)
        pysam.index(cram)
        with ReadEvidence(bam, self.ref) as b, ReadEvidence(cram, self.ref) as c:
            self.assertEqual((b.format, c.format), ('BAM', 'CRAM'))
            from_bam, from_cram = b.evaluate(self.snv_key(405)), c.evaluate(self.snv_key(405))
        self.assertEqual((from_cram.state, from_cram.ref_reads, from_cram.alt_reads), (R.ALTERNATE_PRESENT, 7, 6))
        self.assertEqual(from_bam.columns('NO_CALL', ()), from_cram.columns('NO_CALL', ()))

    def test_an_alignment_with_its_index_elsewhere_is_opened_through_the_index_argument(self):
        path = bamfix.write_bam(os.path.join(self.dir, 'idx.bam'), self.snv_reads(405, n_ref=12), SYN_CONTIGS)
        moved = os.path.join(self.dir, 'elsewhere.bai')
        os.rename(path + '.bai', moved)
        with self.assertRaises(OSError):
            ReadEvidence(path, self.ref)
        with ReadEvidence(path, self.ref, index=moved) as e:
            self.assertEqual(e.evaluate(self.snv_key(405)).state, R.REFERENCE_ONLY)

    def test_an_unindexed_alignment_is_refused_at_open(self):
        path = bamfix.write_bam(os.path.join(self.dir, 'noindex.bam'), self.snv_reads(405, n_ref=3), SYN_CONTIGS)
        os.remove(path + '.bai')
        with self.assertRaises(OSError):
            ReadEvidence(path, self.ref)

    def test_read_window_reports_why_a_read_is_not_informative(self):
        h = pysam.AlignmentHeader.from_dict({'HD': {'VN': '1.6'}, 'SQ': [{'SN': '1', 'LN': 3000}]})

        def aln(cigar, n):                                                                                  # aligned from the 0-based reference position 100
            a = pysam.AlignedSegment(h)
            a.reference_id, a.reference_start, a.cigarstring, a.query_sequence = 0, 100, cigar, 'A' * n
            return a
        self.assertEqual(R.read_window(aln('50M', 50), 109, 111), (9, 11, None))
        self.assertEqual(R.read_window(aln('10M5I40M', 55), 109, 111), (9, 16, None))                       # the five inserted bases (query 10..14) lie inside the window
        self.assertEqual(R.read_window(aln('10M5I40M', 55), 110, 111), (15, 16, None))                      # an insertion before the window does not
        self.assertEqual(R.read_window(aln('10M5D40M', 50), 109, 120), (9, 15, None))                       # a deletion inside the window is the allele, not a failure
        self.assertEqual(R.read_window(aln('10M5D40M', 50), 111, 120)[2], 'window_end_in_a_gap')
        self.assertEqual(R.read_window(aln('10M5N40M', 50), 105, 120)[2], 'reference_skip_in_the_window')
        self.assertEqual(R.read_window(aln('50M', 50), 90, 120)[2], 'not_spanning_the_window')
        self.assertEqual(R.read_window(aln('50M', 50), 140, 160)[2], 'not_spanning_the_window')


class TestWhatTheReadsSettle(unittest.TestCase):
    def test_reads_answer_the_presence_question_for_absent_no_call_filtered_and_low_confidence_loci_only(self):
        for status, reasons in (('ABSENT_FROM_CALLSET', ('no_record_and_no_reference_evidence',)), ('NO_CALL', ('site_without_genotype',)), ('FILTERED_CALL', ('LowQual',)),
                                ('UNRESOLVED', ('reference_confidence_below_policy',)), ('UNRESOLVED', ('overlapping_call_fails_quality',)), ('UNRESOLVED', ('phase_required',))):
            for state in R.ANSWERS:
                self.assertTrue(R.settles(status, reasons, state), (status, reasons, state))
            for state in (R.AMBIGUOUS, R.INSUFFICIENT, R.NOT_EVALUATED):
                self.assertFalse(R.settles(status, reasons, state), (status, reasons, state))
        for status, reasons in (('UNRESOLVED', ('inside_called_structural_variant',)), ('UNRESOLVED', ('inside_called_deletion_on_some_haplotypes',)), ('UNRESOLVED', ('duplicate_records_disagree',)),
                                ('EXACT_MATCH', ()), ('VERIFIED_REFERENCE', ())):
            self.assertFalse(R.settles(status, reasons, R.ALTERNATE_PRESENT), (status, reasons))

    def test_interpretation_names_what_the_reads_say_about_the_callset(self):
        self.assertEqual(R.interpretation('FILTERED_CALL', R.ALTERNATE_PRESENT), 'reads_support_the_filtered_call')
        self.assertEqual(R.interpretation('FILTERED_CALL', R.REFERENCE_ONLY), 'reads_contradict_the_filtered_call')
        self.assertEqual(R.interpretation('ABSENT_FROM_CALLSET', R.ALTERNATE_PRESENT), 'reads_carry_the_allele_the_callset_does_not_list')
        self.assertEqual(R.interpretation('NO_CALL', R.REFERENCE_ONLY), 'reads_show_the_reference_allele_only')
        for state in R.READ_STATES:
            for status in ('ABSENT_FROM_CALLSET', 'FILTERED_CALL', 'UNRESOLVED'):
                self.assertTrue(R.interpretation(status, state))

    def test_columns_are_those_of_the_table_and_a_read_support_without_reads_is_empty(self):
        s = R.ReadSupport('1:5:A:G', R.NOT_EVALUATED, ('contig_not_in_the_alignment_header',))
        c = s.columns('ABSENT_FROM_CALLSET', ())
        self.assertEqual((c['read_state'], c['read_depth'], c['read_alt_fraction'], c['read_settles'], c['read_detail_json']), (R.NOT_EVALUATED, 0, None, False, None))


if __name__ == '__main__':
    unittest.main()
