"""The cross-check of the read-evidence counts against `bcftools mpileup`: that the two agree where they count the same thing, that the differences bcftools' overlap handling makes
are the ones the check names (each shown here on a synthetic pair), and that the check fails when the resolver's counts or its bookkeeping are wrong. The alignments are synthetic."""
import json
import os
import shutil
import tempfile
import unittest
from unittest import mock

import pyarrow as pa
import pyarrow.parquet as pq
import pysam

import bamfix
from bamfix import make_read, span
from helpers import SYN_CONTIGS, SYN_FASTA, find_tool, syn_reference
import crosscheck_reads_bcftools as X
from vindex import reads as R

BCFTOOLS = find_tool('bcftools')
SUBST = {'A': 'G', 'C': 'T', 'G': 'A', 'T': 'C'}
POLICY = dict(min_mapq=20, min_baseq=20)
FLAG_A = bamfix.FLAG_PAIRED | bamfix.FLAG_MATE_REVERSE | bamfix.FLAG_READ1
FLAG_B = bamfix.FLAG_PAIRED | bamfix.FLAG_REVERSE | bamfix.FLAG_READ2

# a site and, beside it, the fragments (12 reference and 3 alternate single reads) the case adds one pair or a few reads to
NAMES = ('plain', 'agree', 'disagree', 'gap', 'gap_kept', 'low_agree', 'low_disagree', 'improper', 'improper_disagree', 'excluded', 'clipped',
         'lowmapq_mate', 'dup_mate', 'qcfail_mate', 'secondary_mate', 'supplementary_mate', 'half_proper', 'boundary_disagree', 'boundary_low_agree',
         'improper_low_agree', 'improper_gap', 'gap_poor', 'improper_agree_poor_mate')
SITES = {name: 300 + 100 * i for i, name in enumerate(NAMES)}
NOT_PROPER = (False, False)
MASK = 0xFFFFFFFF


def keeps_left_mate(name):
    """Whether htslib's overlap handling (sam.c, tweak_overlap_quality) lets the left mate of a pair keep the quality where the mates agree or one has a deletion, as bcftools mpileup
    sees it: it picks the mate from a hash of the read name (X31 string hash, then Wang's integer hash), the other mate's quality is zeroed or, opposite a deletion, scaled by 0.8."""
    b = name.encode()
    h = b[0]
    for c in b[1:]:
        h = ((h << 5) - h + c) & MASK
    h = (h + (~(h << 15) & MASK)) & MASK
    h ^= h >> 10
    h = (h + (h << 3)) & MASK
    h ^= h >> 6
    h = (h + (~(h << 11) & MASK)) & MASK
    h ^= h >> 16
    return bool(h & 1)


def name_where(keeps_left, prefix):
    """A read name, from the prefix, for which keeps_left_mate gives the answer asked for."""
    return next(f'{prefix}-{i}' for i in range(1000) if keeps_left_mate(f'{prefix}-{i}') == keeps_left)


class TestAgainstBcftools(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        if not BCFTOOLS:
            raise unittest.SkipTest('bcftools not available')
        cls.dir = tempfile.mkdtemp()
        cls.fa = pysam.FastaFile(SYN_FASTA)
        cls.ref = syn_reference()
        reads = []
        for name, site in SITES.items():
            reads += cls.fragments(site, 12, 3, name)
        reads += cls.pair(SITES['agree'], 'agree-pair')
        reads += cls.pair(SITES['disagree'], 'disagree-pair', b_base=SUBST[cls.base(SITES['disagree'])], b_q=30)
        reads += cls.pair(SITES['gap'], name_where(False, 'gap-pair'), b_gap=True)               # the covering mate (a) is the one the hash zeroes
        reads += cls.pair(SITES['gap_kept'], name_where(True, 'gap-kept-pair'), b_gap=True)      # the covering mate is the one the hash scales by 0.8: Q35 becomes Q28 and is counted
        reads += cls.pair(SITES['low_agree'], 'low-agree-pair', a_q=12, b_q=12)
        reads += cls.pair(SITES['low_disagree'], 'low-disagree-pair', a_q=24, b_base=SUBST[cls.base(SITES['low_disagree'])], b_q=12)
        reads += cls.pair(SITES['improper'], 'improper-pair', proper=NOT_PROPER)
        reads += cls.pair(SITES['improper_disagree'], 'improper-disagree-pair', b_base=SUBST[cls.base(SITES['improper_disagree'])], proper=NOT_PROPER)
        # a mate that bcftools and the resolver both set aside, disagreeing with the other: it makes no pair
        alt_at = lambda name: SUBST[cls.base(SITES[name])]
        reads += cls.pair(SITES['lowmapq_mate'], 'lowmapq-mate-pair', b_base=alt_at('lowmapq_mate'), b_mapq=10)
        reads += cls.pair(SITES['dup_mate'], 'dup-mate-pair', b_base=alt_at('dup_mate'), b_flag=bamfix.FLAG_DUPLICATE)
        reads += cls.pair(SITES['qcfail_mate'], 'qcfail-mate-pair', b_base=alt_at('qcfail_mate'), b_flag=bamfix.FLAG_QCFAIL)
        reads += cls.pair(SITES['secondary_mate'], 'secondary-mate-pair', b_base=alt_at('secondary_mate'), b_flag=bamfix.FLAG_SECONDARY)
        reads += cls.pair(SITES['supplementary_mate'], 'supplementary-mate-pair', b_base=alt_at('supplementary_mate'), b_flag=bamfix.FLAG_SUPPLEMENTARY)
        # the edges of the rules: only one mate marked proper; a base at the quality minimum (Q20) disagreeing with a better one; two bases of Q10 that add up to the minimum
        reads += cls.pair(SITES['half_proper'], 'half-proper-pair', proper=(True, False))
        reads += cls.pair(SITES['boundary_disagree'], 'boundary-disagree-pair', a_q=20, b_base=alt_at('boundary_disagree'), b_q=35)
        reads += cls.pair(SITES['boundary_low_agree'], 'boundary-low-agree-pair', a_q=10, b_q=10)
        # what an improper pair or a poor covering base does not allow
        reads += cls.pair(SITES['improper_low_agree'], 'improper-low-agree-pair', a_q=12, b_q=12, proper=NOT_PROPER)
        reads += cls.pair(SITES['improper_gap'], 'improper-gap-pair', b_gap=True, proper=NOT_PROPER)
        reads += cls.pair(SITES['gap_poor'], 'gap-poor-pair', a_q=12, b_gap=True)
        reads += cls.pair(SITES['improper_agree_poor_mate'], 'improper-agree-poor-mate-pair', b_q=12, proper=NOT_PROPER)
        s = SITES['excluded']
        alt = SUBST[cls.base(s)]
        for flag, mapq, quals, label in ((0, 10, 35, 'lowmapq'), (bamfix.FLAG_DUPLICATE, 60, 35, 'dup'), (bamfix.FLAG_QCFAIL, 60, 35, 'qcfail'), (bamfix.FLAG_SECONDARY, 60, 35, 'secondary'),
                                         (bamfix.FLAG_SUPPLEMENTARY, 60, 35, 'supplementary'), (0, 60, 10, 'lowbq')):
            for i in range(2):
                reads.append(span(cls.fa, '1', s - 20 + i, s + 19 + i, f'excluded-{label}{i}', {s: alt}, flag=flag, mapq=mapq, quals=quals))
        s = SITES['clipped']
        reads.append(make_read(cls.fa, '1', s - 20, [('M', 20), ('S', 'ACGT')], 'clipped-before', quals=35))              # ends just before the site: no base there
        cls.path = bamfix.write_bam(os.path.join(cls.dir, 'pairs.bam'), reads, SYN_CONTIGS)
        cls.rows = {}
        with R.ReadEvidence(cls.path, cls.ref) as ev:
            for name, site in SITES.items():
                key = f'1:{site}:{cls.base(site)}:{SUBST[cls.base(site)]}'
                cls.rows[name] = dict(key=key, **ev.evaluate(key).columns('ABSENT_FROM_CALLSET'))
        cls.af = pysam.AlignmentFile(cls.path)
        cls.summary, lines = X.compare(list(cls.rows.values()), cls.path, None, POLICY, cls.ref, BCFTOOLS, processes=1, chunk=4)
        cls.line = {name: dict(zip(X.HEADER, next(l for l in lines if l[0] == cls.rows[name]['key']))) for name in SITES}

    @classmethod
    def tearDownClass(cls):
        cls.af.close()
        shutil.rmtree(cls.dir, ignore_errors=True)

    @classmethod
    def base(cls, pos):
        return cls.fa.fetch('1', pos - 1, pos)

    @classmethod
    def fragments(cls, site, n_ref, n_alt, label):
        """Single reads of 40 bases around the site, starts spread over 20 offsets, strands alternating."""
        out, alt = [], SUBST[cls.base(site)]
        for i in range(n_ref + n_alt):
            start = site - 35 + (i % 20)
            out.append(span(cls.fa, '1', start, start + 39, f'{label}-single{i}', {} if i < n_ref else {site: alt}, flag=bamfix.FLAG_REVERSE if i % 2 else 0))
        return out

    @classmethod
    def pair(cls, site, name, a_base=None, b_base=None, a_q=35, b_q=35, b_gap=False, proper=(True, True), b_mapq=60, b_flag=0):
        """Two mates of 40 bases that overlap the site; mate a is forward and starts 30 before it, mate b is reverse and starts 10 before it (with a 2-base deletion over the site if asked).
        `proper` says which mates carry the proper-pair flag; `b_flag` and `b_mapq` are those of mate b."""
        a_start, b_start = site - 30, site - 10
        flag_a, flag_b = FLAG_A | (bamfix.FLAG_PROPER if proper[0] else 0), FLAG_B | (bamfix.FLAG_PROPER if proper[1] else 0) | b_flag
        a = span(cls.fa, '1', a_start, a_start + 39, name, {site: a_base} if a_base else {}, flag=flag_a, quals=[35] * 30 + [a_q] + [35] * 9)
        if b_gap:
            b = make_read(cls.fa, '1', b_start, [('M', 10), ('D', 2), ('M', 30)], name, flag=flag_b, mapq=b_mapq)
        else:
            b = span(cls.fa, '1', b_start, b_start + 39, name, {site: b_base} if b_base else {}, flag=flag_b, mapq=b_mapq, quals=[35] * 10 + [b_q] + [35] * 29)
        a.mate_start1, b.mate_start1 = b_start, a_start
        return [a, b]

    def counts(self, name):
        line = self.line[name]
        return {k: int(line[k]) for k in ('resolver_ref', 'resolver_alt', 'resolver_fragments', 'informative_alignments', 'conflicting_fragments', 'bcftools_ref', 'bcftools_alt',
                                          'bcftools_fragments', 'bcftools_alignments')}

    # ---- where the two count the same thing

    def test_every_locus_agrees_alignment_by_alignment(self):
        self.assertEqual(self.summary['alignments'], {'exact': len(SITES)})
        self.assertEqual((self.summary['n_alignment_problems'], self.summary['errors'], self.summary['pileups_where_dp_differs_from_the_allele_counts']), (0, [], 0))
        for name in SITES:
            self.assertEqual(self.line[name]['alignments_verdict'], 'exact', name)
            self.assertEqual(self.counts(name)['bcftools_alignments'], self.counts(name)['informative_alignments'], name)

    def test_single_reads_and_a_pair_that_agrees_are_counted_the_same_way(self):
        self.assertEqual((self.line['plain']['fragments_verdict'], self.counts('plain')['resolver_ref'], self.counts('plain')['bcftools_ref'], self.counts('plain')['bcftools_alt']), ('exact', 12, 12, 3))
        c = self.counts('agree')                                                   # two alignments, one fragment: bcftools, with overlap detection, counts it once
        self.assertEqual((self.line['agree']['fragments_verdict'], c['resolver_ref'], c['bcftools_ref'], c['informative_alignments']), ('exact', 13, 13, 17))

    def test_what_both_set_aside_is_set_aside(self):
        c = self.counts('excluded')                                                 # alternate-carrying reads: low mapping quality, duplicate, QC fail, secondary, supplementary, low base quality
        self.assertEqual((self.line['excluded']['fragments_verdict'], c['resolver_alt'], c['bcftools_alt'], c['bcftools_ref']), ('exact', 3, 3, 12))
        self.assertEqual(self.counts('clipped')['informative_alignments'], 15)      # a clipped read that ends before the site has no base there

    # ---- where bcftools' overlap handling differs, and what the check calls it

    def test_mates_that_disagree_are_set_aside_by_the_resolver_and_resolved_by_bcftools(self):
        c, line = self.counts('disagree'), self.line['disagree']
        self.assertEqual((c['resolver_ref'], c['conflicting_fragments'], c['bcftools_ref'], c['bcftools_alt']), (12, 1, 13, 3))
        self.assertEqual((line['fragments_verdict'], line['explainable_by']), ('explained', 'mates_disagree'))

    def test_a_mate_opposite_a_deletion_can_cost_the_other_mate_its_base_in_bcftools(self):
        c, line = self.counts('gap'), self.line['gap']
        self.assertEqual((c['resolver_ref'], c['bcftools_ref'], json.loads(self.rows['gap']['read_detail_json'])['filtered']['window_end_in_a_gap']), (13, 12, 1))
        self.assertEqual((line['fragments_verdict'], line['explainable_by']), ('explained', 'mate_has_a_gap_over_the_site'))
        c, line = self.counts('gap_kept'), self.line['gap_kept']                                       # the same pair under a name whose hash scales the base instead: it is still counted
        self.assertEqual((c['resolver_ref'], c['bcftools_ref'], line['fragments_verdict']), (13, 13, 'exact'))

    def test_two_bases_below_the_quality_minimum_that_agree_add_up_in_bcftools(self):
        c, line = self.counts('low_agree'), self.line['low_agree']
        self.assertEqual((c['resolver_ref'], c['bcftools_ref'], c['informative_alignments']), (12, 13, 15))
        self.assertEqual((line['fragments_verdict'], line['explainable_by']), ('explained', 'low_quality_pair_agrees'))

    def test_a_good_base_beside_a_disagreeing_poor_one_is_penalised_by_bcftools(self):
        c, line = self.counts('low_disagree'), self.line['low_disagree']
        self.assertEqual((c['resolver_ref'], c['bcftools_ref'], c['bcftools_alt'], c['informative_alignments']), (13, 12, 3, 16))
        self.assertEqual((line['fragments_verdict'], line['explainable_by']), ('explained', 'low_quality_mate_disagrees'))

    def test_a_pair_not_marked_proper_is_not_de_overlapped_by_bcftools(self):
        c, line = self.counts('improper'), self.line['improper']
        self.assertEqual((c['resolver_ref'], c['bcftools_ref'], c['informative_alignments']), (13, 14, 17))
        self.assertEqual((line['fragments_verdict'], line['explainable_by']), ('explained', 'agreeing_pair_not_marked_proper'))
        c, line = self.counts('improper_disagree'), self.line['improper_disagree']
        self.assertEqual((c['resolver_ref'], c['resolver_alt'], c['conflicting_fragments'], c['bcftools_ref'], c['bcftools_alt']), (12, 3, 1, 13, 4))
        self.assertEqual((line['fragments_verdict'], line['explainable_by']), ('explained', 'disagreeing_pair_not_marked_proper'))        # both mates counted, each its own base

    def test_the_summary_counts_what_was_compared(self):
        s = self.summary
        self.assertEqual((s['loci'], s['snv_keys_of_the_run'], s['not_compared'], s['fragments'], s['n_unexplained'], s['ok']), (len(SITES), len(SITES), 0, {'exact': 14, 'explained': 9}, 0, True))
        self.assertEqual(s['explainable_by'], {'low_quality_mate_disagrees': 1, 'low_quality_pair_agrees': 2, 'mate_has_a_gap_over_the_site': 1, 'mates_disagree': 2,
                                               'agreeing_pair_not_marked_proper': 2, 'disagreeing_pair_not_marked_proper': 1})
        self.assertEqual({e['key']: e['pairs'] for e in s['explained']}[self.rows['gap']['key']], {'gap_mate': 1})                         # the pairs behind a difference are recorded

    # ---- the pairs, read from the alignments

    def patterns(self, name):
        return X.pair_patterns(self.af, '1', SITES[name], POLICY['min_mapq'], POLICY['min_baseq'])

    def test_the_kinds_of_pair_and_what_each_allows(self):
        # (up, down, up effects, down effects, kinds)
        self.assertEqual(self.patterns('agree'), (0, 0, {}, {}, {'agree': 1}))
        self.assertEqual(self.patterns('disagree'), (1, 0, {'mates_disagree': 1}, {}, {'disagree': 1}))
        self.assertEqual(self.patterns('gap'), (0, 1, {}, {'mate_has_a_gap_over_the_site': 1}, {'gap_mate': 1}))
        self.assertEqual(self.patterns('gap_kept'), (0, 1, {}, {'mate_has_a_gap_over_the_site': 1}, {'gap_mate': 1}))                # an allowance, whatever the hash does with it
        self.assertEqual(self.patterns('low_agree'), (1, 0, {'low_quality_pair_agrees': 1}, {}, {'agree_both_low': 1}))
        self.assertEqual(self.patterns('low_disagree'), (0, 1, {}, {'low_quality_mate_disagrees': 1}, {'disagree_low_mate': 1}))
        self.assertEqual(self.patterns('improper'), (1, 0, {'agreeing_pair_not_marked_proper': 1}, {}, {'agree/not_proper': 1}))
        self.assertEqual(self.patterns('improper_disagree'), (2, 0, {'disagreeing_pair_not_marked_proper': 2}, {}, {'disagree/not_proper': 1}))

    def test_single_reads_and_reads_the_filters_remove_make_no_pair(self):
        for name in ('plain', 'excluded', 'clipped'):
            self.assertEqual(self.patterns(name), (0, 0, {}, {}, {}), name)

    def test_a_mate_the_filters_remove_makes_no_pair(self):                         # it disagrees with the other mate, and both bcftools and the resolver set it aside
        for name in ('lowmapq_mate', 'dup_mate', 'qcfail_mate', 'secondary_mate', 'supplementary_mate'):
            self.assertEqual(self.patterns(name), (0, 0, {}, {}, {}), name)
            self.assertEqual((self.line[name]['fragments_verdict'], self.counts(name)['resolver_ref'], self.counts(name)['bcftools_ref']), ('exact', 13, 13), name)

    def test_the_edges_of_the_rules(self):
        # a pair with only one mate marked proper is not de-overlapped; a base at the minimum (Q20) is good; two bases of Q10 reach the minimum together
        self.assertEqual(self.patterns('half_proper'), (1, 0, {'agreeing_pair_not_marked_proper': 1}, {}, {'agree/not_proper': 1}))
        self.assertEqual((self.counts('half_proper')['resolver_ref'], self.counts('half_proper')['bcftools_ref'], self.line['half_proper']['fragments_verdict']), (13, 14, 'explained'))
        self.assertEqual(self.patterns('boundary_disagree'), (1, 0, {'mates_disagree': 1}, {}, {'disagree': 1}))
        self.assertEqual((self.counts('boundary_disagree')['resolver_alt'], self.counts('boundary_disagree')['bcftools_alt'], self.line['boundary_disagree']['fragments_verdict']), (3, 4, 'explained'))
        self.assertEqual(self.patterns('boundary_low_agree'), (1, 0, {'low_quality_pair_agrees': 1}, {}, {'agree_both_low': 1}))
        self.assertEqual((self.counts('boundary_low_agree')['resolver_ref'], self.counts('boundary_low_agree')['bcftools_ref'], self.line['boundary_low_agree']['fragments_verdict']), (12, 13, 'explained'))

    def test_what_an_improper_pair_or_a_poor_covering_base_does_not_allow(self):
        for name, kinds in (('improper_low_agree', {'agree_both_low/not_proper': 1}), ('improper_gap', {'gap_mate/not_proper': 1}), ('gap_poor', {'gap_mate': 1}),
                            ('improper_agree_poor_mate', {'agree/not_proper': 1})):
            self.assertEqual(self.patterns(name), (0, 0, {}, {}, kinds), name)
            self.assertEqual(self.line[name]['fragments_verdict'], 'exact', name)                     # bcftools and the resolver count the same

    # ---- what the check does when it cannot compare

    def test_a_site_whose_pileup_failed_is_not_compared_and_the_check_does_not_pass(self):
        real = X.work

        def flaky(task):                                                            # the first site of every call is lost, as when a transfer fails
            overlaps, positions, pile, errs = real(task)
            gone = positions[0]
            return overlaps, positions, {k: v for k, v in pile.items() if k != gone}, errs + [dict(site=list(gone), error='simulated')]
        with mock.patch.object(X, 'work', flaky):
            s = X.compare(list(self.rows.values()), self.path, None, POLICY, self.ref, BCFTOOLS, processes=1, chunk=4)[0]
        chunks = -(-len(SITES) // 4)
        self.assertEqual((s['loci'], s['not_compared'], len(s['errors']), s['ok']), (len(SITES) - chunks, chunks, 2 * chunks, False))
        self.assertEqual(s['alignments'], {'exact': len(SITES) - chunks})            # the sites that were compared are not affected

    def test_a_check_that_compared_nothing_does_not_pass(self):
        s = X.compare([], self.path, None, POLICY, self.ref, BCFTOOLS, processes=1)[0]
        self.assertEqual((s['loci'], s['alignments'], s['fragments'], s['ok']), (0, {}, {}, False))

    # ---- the check fails when the resolver is wrong

    def compare(self, row):
        return X.compare([row], self.path, None, POLICY, self.ref, BCFTOOLS, processes=1)[0]

    def test_a_wrong_count_is_unexplained_and_fails_both_comparisons(self):
        row = dict(self.rows['plain'], read_ref=self.rows['plain']['read_ref'] + 3, read_depth=self.rows['plain']['read_depth'] + 3)
        s = self.compare(row)
        self.assertEqual((s['alignments'], s['fragments'], s['ok']), ({'unexplained': 1}, {'unexplained': 1}, False))
        self.assertTrue(s['alignment_problems'])                                                                # listed with the reasons
        self.assertEqual(s['unexplained'][0]['differences']['ref'], -3)

    def test_a_wrong_alternate_count_is_caught_where_the_total_is_right(self):
        r = self.rows['plain']
        row = dict(r, read_ref=r['read_ref'] + 1, read_alt=r['read_alt'] - 1)                                   # same fragments, one moved from alternate to reference
        s = self.compare(row)
        self.assertEqual((s['alignments'], s['fragments']), ({'unexplained': 1}, {'unexplained': 1}))

    def test_inconsistent_bookkeeping_is_caught(self):
        r = self.rows['agree']
        detail = json.loads(r['read_detail_json'])
        for tampered in (dict(detail, informative_alignments=detail['informative_alignments'] - 1), dict(detail, conflicting_fragments=3)):
            s = self.compare(dict(r, read_detail_json=json.dumps(tampered)))
            self.assertEqual(s['alignments'], {'unexplained': 1}, tampered)
            self.assertTrue(s['alignment_problems'])

    def test_agreeing_mates_the_resolver_calls_conflicting_are_not_explained_by_anything(self):
        r = self.rows['agree']                                                      # a resolver that sets aside a pair that agrees: bcftools counts one fragment more, and the alignments hold no pair that disagrees
        detail = dict(json.loads(r['read_detail_json']), conflicting_fragments=1)
        s = self.compare(dict(r, read_ref=r['read_ref'] - 1, read_depth=r['read_depth'] - 1, read_detail_json=json.dumps(detail)))
        self.assertEqual((s['fragments'], s['alignments'], s['ok']), ({'unexplained': 1}, {'exact': 1}, False))
        self.assertEqual(s['unexplained'][0]['pairs'], {'agree': 1})

    def test_conflicting_fragments_are_compared_with_the_pairs_of_the_alignments(self):
        r = self.rows['disagree']                                                   # the counts are those of the run, the bookkeeping claims two conflicts, or none, where the alignments hold one pair
        detail = json.loads(r['read_detail_json'])
        self.assertEqual(detail['conflicting_fragments'], 1)
        for claimed in (0, 2):
            s = self.compare(dict(r, read_detail_json=json.dumps(dict(detail, conflicting_fragments=claimed))))
            self.assertEqual(s['fragments'], {'unexplained': 1}, claimed)
            self.assertEqual(s['unexplained'][0]['pairs'], {'disagree': 1})

    def test_the_total_of_alignments_is_compared(self):
        r = self.rows['plain']                                                      # one alignment more than bcftools finds, in no count of fragments
        detail = json.loads(r['read_detail_json'])
        s = self.compare(dict(r, read_detail_json=json.dumps(dict(detail, informative_alignments=detail['informative_alignments'] + 1))))
        self.assertEqual((s['alignments'], s['fragments'], s['ok']), ({'unexplained': 1}, {'exact': 1}, False))

    def test_a_fragment_of_another_base_that_bcftools_does_not_have_is_caught_with_the_total_right(self):
        r = self.rows['disagree']                                                   # the reference and total differences are inside the bounds; only the other-base difference is not
        s = self.compare(dict(r, read_other=r['read_other'] + 1, read_depth=r['read_depth'] + 1))
        self.assertEqual(s['fragments'], {'unexplained': 1})
        self.assertEqual(s['unexplained'][0]['differences'], {'ref': 1, 'alt': 0, 'total': 0, 'other': -1})

    def test_an_effect_needs_a_pair_in_the_alignments_not_a_claim_of_the_resolver(self):
        r = self.rows['plain']                                                      # the resolver counts one fragment too many; no pair at the site could make bcftools count one fewer
        s = self.compare(dict(r, read_ref=r['read_ref'] + 1, read_depth=r['read_depth'] + 1, read_detail_json=json.dumps(dict(json.loads(r['read_detail_json']), filtered=dict(window_end_in_a_gap=3)))))
        self.assertEqual(s['fragments'], {'unexplained': 1})

    def test_a_difference_beyond_the_bounds_is_unexplained(self):
        r = self.rows['disagree']                                                   # one pair that disagrees allows bcftools one more fragment than the resolver, not two or three
        s = self.compare(dict(r, read_ref=r['read_ref'] - 1, read_depth=r['read_depth'] - 1))
        self.assertEqual(s['fragments'], {'unexplained': 1})
        s = self.compare(dict(r, read_ref=r['read_ref'] - 2, read_other=r['read_other'] + 2))
        self.assertEqual(s['fragments'], {'unexplained': 1})
        s = self.compare(dict(r, read_ref=r['read_ref'] - 1, read_other=r['read_other'] + 1))
        self.assertEqual(s['fragments'], {'unexplained': 1})                       # the total is right but a base is moved: no effect moves a base from one allele to another


class TestFunctions(unittest.TestCase):
    ROW = dict(read_ref=10, read_alt=2, read_other=1, read_depth=13)

    def detail(self, **kw):
        return dict(dict(informative_alignments=13, conflicting_fragments=0), **kw)

    def test_parse_pileup(self):
        text = ('##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ts\n'
                '1\t300\t.\tN\tA,G,<*>\t0\t.\tDP=14\tPL:DP:AD\t0:14:0,9,5,0\n'
                '1\t310\t.\tN\t<*>\t0\t.\tDP=0\tPL:DP:AD\t0:0:0,0\n'
                '1\t320\t.\tN\tC,<*>\t0\t.\tDP=3\tPL:DP:AD\t0:3:0,2,0\n')
        self.assertEqual(X.parse_pileup(text, {('1', 300), ('1', 310)}), {('1', 300): ({'A': 9, 'G': 5}, 14, True), ('1', 310): ({}, 0, True)})
        self.assertEqual(X.parse_pileup(text.replace(':3:0,2,0', ':4:0,2,0'), {('1', 320)}), {('1', 320): ({'C': 2}, 2, False)})        # DP that differs from the allele counts is recorded

    def test_the_verdicts_at_their_boundaries(self):
        got = lambda a, g, other, dp=True: ({'A': a, 'G': g, 'T': other} if other else {'A': a, 'G': g}, a + g + other, dp)
        base = ('A', 'G')
        none = lambda: (0, 0, {}, {}, {})
        never = lambda: self.fail('the pairs are looked up only when the counts differ')
        self.assertEqual(X.fragments_verdict(self.ROW, base, got(10, 2, 1), never, 0)[0], 'exact')
        self.assertEqual(X.fragments_verdict(self.ROW, base, got(11, 2, 1), none, 0)[0], 'unexplained')                   # a difference no pair allows
        two_more = lambda: (2, 0, {'disagreeing_pair_not_marked_proper': 2}, {}, {'disagree/not_proper': 1})
        verdict = X.fragments_verdict(self.ROW, base, got(12, 2, 1), two_more, 1)
        self.assertEqual((verdict[0], verdict[4]), ('explained', ['disagreeing_pair_not_marked_proper']))
        self.assertEqual(X.fragments_verdict(self.ROW, base, got(13, 2, 1), two_more, 1)[0], 'unexplained')
        one_fewer = lambda: (0, 1, {}, {'mate_has_a_gap_over_the_site': 1}, {'gap_mate': 1})
        self.assertEqual(X.fragments_verdict(self.ROW, base, got(9, 2, 1), one_fewer, 0)[0], 'explained')
        self.assertEqual(X.fragments_verdict(self.ROW, base, got(8, 2, 1), one_fewer, 0)[0], 'unexplained')
        self.assertEqual(X.fragments_verdict(self.ROW, base, got(9, 3, 1), one_fewer, 0)[0], 'unexplained')               # more of one base and fewer of another: the up allowance is zero
        both = lambda: (1, 1, {'mates_disagree': 1}, {'mate_has_a_gap_over_the_site': 1}, {'disagree': 1, 'gap_mate': 1})
        self.assertEqual(X.fragments_verdict(self.ROW, base, got(9, 3, 1), both, 1)[0], 'explained')                      # one more alternate and one fewer reference

    def test_the_conflicting_fragments_are_the_pairs_whose_mates_show_different_bases(self):
        got = ({'A': 10, 'G': 2, 'T': 1}, 13, True)                                                                         # counts that agree: the pairs are not looked up
        self.assertEqual(X.fragments_verdict(self.ROW, ('A', 'G'), got, lambda: self.fail('not looked up'), 3)[0], 'exact')
        mates = lambda: (1, 0, {'mates_disagree': 1}, {}, {'disagree': 1})
        more = ({'A': 11, 'G': 2, 'T': 1}, 14, True)
        self.assertEqual(X.fragments_verdict(self.ROW, ('A', 'G'), more, mates, 1)[0], 'explained')
        self.assertEqual(X.fragments_verdict(self.ROW, ('A', 'G'), more, mates, 0)[0], 'unexplained')                      # the resolver claims no conflict, the alignments hold one
        self.assertEqual(X.fragments_verdict(self.ROW, ('A', 'G'), more, mates, 2)[0], 'unexplained')                      # the resolver claims two, the alignments hold one
        self.assertEqual(X.fragments_verdict(self.ROW, ('A', 'G'), more, lambda: (1, 0, {'mates_disagree': 1}, {}, {'agree': 1}), 1)[0], 'unexplained')    # a conflict where the mates agree

    def test_alignment_verdict_has_no_tolerance(self):
        base = ('A', 'G')
        mk = lambda a, g, other: ({'A': a, 'G': g, 'T': other} if other else {'A': a, 'G': g}, a + g + other, True)
        self.assertEqual(X.alignments_verdict(self.ROW, self.detail(), base, mk(10, 2, 1)), [])
        self.assertTrue(X.alignments_verdict(self.ROW, self.detail(), base, mk(10, 2, 2)))                               # one alignment more than the resolver saw
        self.assertTrue(X.alignments_verdict(self.ROW, self.detail(), base, mk(9, 2, 1)))                                # and one fewer
        self.assertTrue(X.alignments_verdict(self.ROW, self.detail(), base, mk(11, 1, 1)))                               # the right total, a base moved: one fragment per alignment, so equal counts
        self.assertEqual(X.alignments_verdict(self.ROW, self.detail(informative_alignments=15), base, mk(12, 2, 1)), [])  # two mates agree on a base: two alignments for one fragment
        self.assertTrue(X.alignments_verdict(self.ROW, self.detail(informative_alignments=15), base, mk(9, 4, 2)))       # fewer alignments carrying the reference base than fragments that do

    def test_a_bcftools_that_cannot_run_is_an_error_not_a_silent_miss(self):
        if not BCFTOOLS:
            self.skipTest('bcftools not available')
        X._init('/nonexistent/reads.bam', None, BCFTOOLS, 20, 20)
        with self.assertRaises(RuntimeError):
            X.pileup([('1', 300)], True)

    def test_the_rows_that_are_compared_are_the_snv_keys_the_reads_were_asked_about(self):
        d = tempfile.mkdtemp()
        try:
            keys = [('1:100:A:G', R.REFERENCE_ONLY), ('1:200:A:C', R.ALTERNATE_PRESENT), ('1:300:A:AT', R.REFERENCE_ONLY), ('1:400:AT:A', R.INSUFFICIENT), ('1:500:A:T', R.NOT_EVALUATED),
                    ('1:600:A:G', R.REFERENCE_ONLY), ('1:700:A:G', None), ('1:800:AC:GT', R.REFERENCE_ONLY)]
            schema = pa.schema([('key', pa.string()), ('read_state', pa.string())])
            pq.write_table(pa.Table.from_pylist([dict(key=k, read_state=s) for k, s in keys], schema=schema), os.path.join(d, 'u.parquet'))
            rows = X.snv_rows(os.path.join(d, 'u.parquet'))
            self.assertEqual([r['key'] for r in rows], ['1:100:A:G', '1:200:A:C', '1:600:A:G'])
            capped = X.snv_rows(os.path.join(d, 'u.parquet'), limit_per_state=1, seed=1)
            self.assertEqual(sorted(r['read_state'] for r in capped), [R.ALTERNATE_PRESENT, R.REFERENCE_ONLY])
            self.assertEqual([r['key'] for r in X.snv_rows(os.path.join(d, 'u.parquet'), limit_per_state=1, seed=1)], [r['key'] for r in capped])        # the same sample for the same seed
        finally:
            shutil.rmtree(d, ignore_errors=True)


if __name__ == '__main__':
    unittest.main()
