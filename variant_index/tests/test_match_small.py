import os
import tempfile
import unittest

from helpers import syn_reference, vcf_text
from evidence_fixtures import ev, index_from
from vindex.vcfio import read_vcf
from vindex.callability import CallableRegionCallability
from vindex.match_small import MatchContext, status_of, interpret_calls


class Base(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.ref = syn_reference()
        cls.td = tempfile.mkdtemp()

    def b(self, pos, contig='1'):
        return self.ref.fetch(contig, pos, pos)

    def alt(self, pos, contig='1', skip=()):
        base = self.b(pos, contig)
        return next(x for x in 'ACGT' if x != base and x not in skip)

    def read(self, recs, gvcf=False, **kw):
        p = os.path.join(self.td, 'x.vcf')
        with open(p, 'w') as f:
            f.write(vcf_text(recs, **kw))
        return read_vcf(p, self.ref)

    def ctx(self, cs, *extra):
        return MatchContext(cs, extra)


class TestExactAndNonMatches(Base):
    def test_equivalent_representation_matches_normalised_index_key(self):
        # a deletion of one A in the homopolymer, written at 105 (right-shifted); the index holds the left-aligned form
        cs = self.read([('1', 105, 'AA', 'A', 'PASS', 50, ['0/1:50:30:15,15'])])
        st = status_of(self.ctx(cs), '1:101:CA:C')
        self.assertEqual((st.status, st.dosage, st.zygosity), ('EXACT_MATCH', 1, 'heterozygous'))

    def test_same_position_other_allele_is_not_a_match(self):
        a1, a2 = self.alt(402), None
        a2 = self.alt(402, skip=(a1,))
        cs = self.read([('1', 402, self.b(402), a1, 'PASS', 50, ['0/1:50:30:15,15'])])
        st = status_of(self.ctx(cs), f'1:402:{self.b(402)}:{a2}')
        self.assertEqual(st.status, 'OTHER_ALLELE_AT_SITE')
        idx = index_from([ev(f'1:402:{self.b(402)}:{a2}')])
        out = interpret_calls(cs, idx, self.ctx(cs))
        self.assertEqual([o.status for o in out], ['NO_INDEX_RECORD'])

    def test_rsid_equality_is_never_identity(self):
        a1 = self.alt(402)
        a2 = self.alt(402, skip=(a1,))
        cs = self.read([('1', 402, self.b(402), a1, 'PASS', 50, ['0/1:50:30:15,15'], '.', 'rs12345')])
        idx = index_from([ev(f'1:402:{self.b(402)}:{a2}', extra_json='{"rsid": "rs12345"}')])
        out = interpret_calls(cs, idx)
        self.assertEqual(out[0].status, 'NO_INDEX_RECORD')
        self.assertEqual(out[0].evidence, [])

    def test_exact_match_carries_evidence_and_axes(self):
        a1 = self.alt(402)
        key = f'1:402:{self.b(402)}:{a1}'
        rows = [ev(key, 'clinvar', claim='Pathogenic', lineage='trial1'),
                ev(key, 'gwascat', axis='association_replication', assertion='effect', claim_type='association', kind='statistical_association', lineage='study9', p_value=1e-9),
                ev(key, 'other', axis='association_replication', assertion='no_effect', claim_type='association', kind='statistical_association', lineage='study9', accession='o1')]
        cs = self.read([('1', 402, self.b(402), a1, 'PASS', 50, ['1/1:60:30:0,30'])])
        out = interpret_calls(cs, index_from(rows))
        o = out[0]
        self.assertEqual(o.status, 'EXACT_MATCH')
        self.assertEqual(set(o.axes), {'clinical_validity_actionability', 'association_replication'})
        ar = o.axes['association_replication']
        self.assertEqual((ar['n_rows'], ar['n_independent'], ar['conflict']), (2, 1, True))     # shared lineage counts once; effect vs no_effect both kept

    def test_rejected_and_unused_alleles_are_reported_not_dropped(self):
        wrong = self.alt(402)
        cs = self.read([('1', 402, wrong, self.alt(402, skip=(wrong,)), 'PASS', 50, ['0/1:50:30:15,15']),
                        ('1', 410, self.b(410), f'{self.alt(410)},{self.alt(410, skip=(self.alt(410),))}', 'PASS', 50, ['0/1:50:30:15,15,0'])])
        self.assertEqual([r.reason for r in cs.rejected], ['ref_mismatch'])
        self.assertEqual([r.reason for r in cs.unused], ['allele_not_in_genotype'])
        self.assertEqual(len(cs.small), 1)


class TestMissingAndReference(Base):
    def test_status_taxonomy(self):
        a = self.alt(402)
        b410 = self.b(410)
        recs = [('1', 300, 'C', '<NON_REF>', 'PASS', 0, ['0/0:50:30:30,0'], 'END=399', '.', 'GT:GQ:DP:AD'),        # confident reference 300-399
                ('1', 402, self.b(402), a, 'LowQual', 5, ['0/1:10:8:4,4']),                                       # filtered
                ('1', 410, b410, self.alt(410), 'PASS', 50, ['./.:.:.:.']),                                         # no-call
                ('1', 420, 'A', '<NON_REF>', 'PASS', 0, ['0/0:10:30:30,0'], 'END=430', '.', 'GT:GQ:DP:AD')]          # low-GQ block
        cs = self.read(recs)
        ctx = self.ctx(cs)
        key = lambda pos, alt=None: f'1:{pos}:{self.b(pos)}:{alt or self.alt(pos)}'
        self.assertEqual(status_of(ctx, key(402, a)).status, 'FILTERED_CALL')
        self.assertEqual(status_of(ctx, key(410)).status, 'NO_CALL')
        self.assertEqual(status_of(ctx, key(350)).status, 'VERIFIED_REFERENCE')
        st = status_of(ctx, key(425))
        self.assertEqual(st.status, 'UNRESOLVED')
        self.assertIn('reference_confidence_below_policy', st.reasons)
        self.assertEqual(status_of(ctx, key(600)).status, 'ABSENT_FROM_CALLSET')

    def test_an_indel_across_a_confident_and_a_low_confidence_block_is_low_confidence_and_the_union_is_built_once(self):
        recs = [('1', 300, 'C', '<NON_REF>', 'PASS', 0, ['0/0:50:30:30,0'], 'END=349', '.', 'GT:GQ:DP:AD'),
                ('1', 350, 'A', '<NON_REF>', 'PASS', 0, ['0/0:10:30:30,0'], 'END=399', '.', 'GT:GQ:DP:AD')]
        ctx = self.ctx(self.read(recs))
        span = self.ref.fetch('1', 347, 352)                                    # positions 348-352: 348 and 349 are in the confident block, 350-352 in the low-confidence one
        st = status_of(ctx, f'1:348:{span}:{span[0]}')
        self.assertEqual((st.status, st.reasons), ('UNRESOLVED', ('reference_confidence_below_policy',)))
        gvcf = ctx.callability.providers[0]
        built = gvcf._union
        self.assertIsNotNone(built)
        self.assertEqual(status_of(ctx, f'1:349:{self.ref.fetch("1", 348, 352)}:{self.b(349)}').status, 'UNRESOLVED')
        self.assertIs(gvcf._union, built, 'the union of the blocks is built once, not for every site that asks')
        self.assertEqual(status_of(ctx, f'1:340:{self.ref.fetch("1", 339, 343)}:{self.b(340)}').status, 'VERIFIED_REFERENCE')

    def test_a_record_the_reader_refused_is_a_record_at_the_site_so_reference_is_not_verified_there(self):
        wrong = next(x for x in 'ACGT' if x not in (self.b(350), self.alt(350)))
        recs = [('1', 300, 'C', '<NON_REF>', 'PASS', 0, ['0/0:50:30:30,0'], 'END=399', '.', 'GT:GQ:DP:AD'),
                ('1', 350, wrong, self.alt(350), 'PASS', 50, ['0/1:60:30:15,15'])]                                  # REF disagrees with the reference, inside a confident block
        cs = self.read(recs)
        self.assertEqual([(r.reason, r.pos) for r in cs.rejected], [('ref_mismatch', 350)])
        ctx = self.ctx(cs)
        st = status_of(ctx, f'1:350:{self.b(350)}:{self.alt(350)}')
        self.assertEqual((st.status, st.reasons, st.basis), ('NO_CALL', ('site_without_genotype',), 'record_refused_by_the_reader'))
        self.assertEqual(status_of(ctx, f'1:351:{self.b(351)}:{self.alt(351)}').status, 'VERIFIED_REFERENCE', 'the neighbouring site is unaffected')
        cs2 = self.read(recs[1:])                                                                                   # the same record in a variant-only file with a callable region
        ctx2 = self.ctx(cs2, CallableRegionCallability([('1', 300, 400)], cs2))
        self.assertEqual(status_of(ctx2, f'1:350:{self.b(350)}:{self.alt(350)}').status, 'NO_CALL')
        self.assertEqual(status_of(ctx2, f'1:351:{self.b(351)}:{self.alt(351)}').status, 'VERIFIED_REFERENCE')

    def test_variant_only_vcf_never_proves_reference(self):
        cs = self.read([('1', 402, self.b(402), self.alt(402), 'PASS', 50, ['0/1:50:30:15,15'])])
        st = status_of(self.ctx(cs), f'1:350:{self.b(350)}:{self.alt(350)}')
        self.assertEqual(st.status, 'ABSENT_FROM_CALLSET')
        self.assertIn('missing_evidence', st.detail)

    def test_callable_region_gives_reference_only_with_its_assumption(self):
        cs = self.read([('1', 402, self.b(402), self.alt(402), 'PASS', 50, ['0/1:50:30:15,15'])])
        ctx = self.ctx(cs, CallableRegionCallability([('1', 300, 500)], cs))
        st = status_of(ctx, f'1:350:{self.b(350)}:{self.alt(350)}')
        self.assertEqual(st.status, 'VERIFIED_REFERENCE')
        self.assertIn('assumes_all_variants_reported', st.basis)
        self.assertEqual(status_of(ctx, f'1:600:{self.b(600)}:{self.alt(600)}').status, 'ABSENT_FROM_CALLSET')   # outside the mask

    def test_spanning_deletion_and_deletion_cover(self):
        cs = self.read([('1', 105, 'AA', 'A,*', 'PASS', 50, ['1/2:50:30:0,15,15']),
                        ('1', 410, self.ref.fetch('1', 410, 413), self.b(410), 'PASS', 50, ['1/1:50:30:0,30'])])
        ctx = self.ctx(cs)
        self.assertEqual(status_of(ctx, f'1:106:{self.b(106)}:{self.alt(106)}').status, 'DELETED_BY_UPSTREAM')
        inside = status_of(ctx, f'1:411:{self.b(411)}:{self.alt(411)}')
        self.assertEqual(inside.status, 'DELETED_BY_UPSTREAM')
        cs2 = self.read([('1', 410, self.ref.fetch('1', 410, 413), self.b(410), 'PASS', 50, ['0/1:50:30:15,15'])])
        het = status_of(self.ctx(cs2), f'1:411:{self.b(411)}:{self.alt(411)}')
        self.assertEqual(het.status, 'UNRESOLVED')

    def test_hemizygous_call(self):
        a = self.alt(105, 'X')
        cs = self.read([('X', 105, self.b(105, 'X'), a, 'PASS', 50, ['1:50:30:0,30'])])
        st = status_of(self.ctx(cs), f'X:105:{self.b(105, "X")}:{a}')
        self.assertEqual((st.status, st.dosage, st.zygosity), ('EXACT_MATCH', 1, 'hemizygous'))


class TestPhaseAndMnv(Base):
    def setUp(self):
        s = self.ref.fetch('1', 500, 501)
        self.p0, self.p1 = s[0], s[1]
        self.a0, self.a1 = self.alt(500), self.alt(501)
        self.mnv = f'1:500:{self.p0 + self.p1}:{self.a0 + self.a1}'

    def recs(self, gt0, gt1, ps=True):
        fmt = 'GT:GQ:DP:AD:PS'
        return [('1', 500, self.p0, self.a0, 'PASS', 50, [f'{gt0}:50:30:15,15:7'], '.', '.', fmt),
                ('1', 501, self.p1, self.a1, 'PASS', 50, [f'{gt1}:50:30:15,15:7'], '.', '.', fmt)]

    def test_unphased_heterozygous_components_are_unresolved_with_bounds(self):
        st = status_of(self.ctx(self.read(self.recs('0/1', '0/1'))), self.mnv)
        self.assertEqual(st.status, 'UNRESOLVED')
        self.assertIn('phase_required', st.reasons)
        self.assertEqual((st.detail['carries_all_components_min'], st.detail['carries_all_components_max']), (0, 1))

    def test_phased_cis_matches_and_trans_does_not(self):
        cis = status_of(self.ctx(self.read(self.recs('0|1', '0|1'))), self.mnv)
        self.assertEqual((cis.status, cis.dosage), ('EXACT_MATCH', 1))
        trans = status_of(self.ctx(self.read(self.recs('0|1', '1|0'))), self.mnv)
        self.assertEqual(trans.status, 'OTHER_ALLELE_AT_SITE')
        self.assertIn('mnv_components_in_trans', trans.reasons)

    def test_homozygous_components_imply_cis(self):
        st = status_of(self.ctx(self.read(self.recs('1/1', '1/1'))), self.mnv)
        self.assertEqual((st.status, st.dosage), ('EXACT_MATCH', 2))
        mixed = status_of(self.ctx(self.read(self.recs('1/1', '0/1'))), self.mnv)
        self.assertEqual((mixed.status, mixed.dosage), ('EXACT_MATCH', 1))          # the haplotype with the second variant must carry the first

    def test_partial_components_with_reference_evidence(self):
        recs = self.recs('0/1', '0/0')[:1] + [('1', 501, self.p1, 'N', 'PASS', 0, ['0/0:50:30:30,0'], 'END=501', '.', 'GT:GQ:DP:AD')]
        recs[1] = ('1', 501, self.p1, '<NON_REF>', 'PASS', 0, ['0/0:50:30:30,0'], 'END=501', '.', 'GT:GQ:DP:AD')
        cs = self.read(recs)
        st = status_of(self.ctx(cs), self.mnv)
        self.assertEqual(st.status, 'OTHER_ALLELE_AT_SITE')
        st2 = status_of(self.ctx(self.read(recs[:1])), self.mnv)
        self.assertEqual(st2.status, 'UNRESOLVED')                                      # the other component is unknown in a variant-only VCF

    def test_sample_mnv_call_supports_index_snv_component(self):
        cs = self.read([('1', 500, self.p0 + self.p1, self.a0 + self.a1, 'PASS', 50, ['0/1:50:30:15,15'])])
        st = status_of(self.ctx(cs), f'1:500:{self.p0}:{self.a0}')
        self.assertEqual((st.status, st.basis, st.dosage), ('EXACT_MATCH', 'component_of_called_mnv', 1))

    def test_call_centric_reports_mnv_context(self):
        cs = self.read(self.recs('0/1', '0/1'))
        idx = index_from([ev(self.mnv, 'clinvar')], variants=[self.mnv])
        out = interpret_calls(cs, idx, self.ctx(cs))
        self.assertTrue(all(o.status == 'NO_INDEX_RECORD' for o in out))
        self.assertEqual(out[0].mnv_context[0][0], self.mnv)
        self.assertEqual(out[0].mnv_context[0][1].status, 'UNRESOLVED')



class TestMissingEvidenceText(unittest.TestCase):
    def test_every_unresolved_status_the_matcher_can_build_states_what_would_settle_it(self):
        """Read the matcher's source: each Status(..., 'UNRESOLVED' | 'NO_CALL', reasons) either carries missing_evidence in its detail or has its reason in MISSING_BY_REASON."""
        import ast
        import vindex.match_small as ms
        with open(ms.__file__) as f:
            tree = ast.parse(f.read())
        seen = []
        for n in ast.walk(tree):
            if isinstance(n, ast.Call) and getattr(n.func, 'id', None) == 'Status' and len(n.args) >= 3 and isinstance(n.args[1], ast.Constant) and n.args[1].value in ('UNRESOLVED', 'NO_CALL'):
                reasons = [e.value for e in getattr(n.args[2], 'elts', []) if isinstance(e, ast.Constant)]
                detail = next((k.value for k in n.keywords if k.arg == 'detail'), None)
                has_own = detail is not None and any(k.arg == 'missing_evidence' for k in getattr(detail, 'keywords', []))
                seen.append((tuple(reasons), has_own))
        self.assertGreaterEqual(len(seen), 8)
        for reasons, has_own in seen:
            self.assertTrue(has_own or any(r in ms.MISSING_BY_REASON for r in reasons), reasons)

    def test_resolved_statuses_need_no_text_and_filtered_calls_name_their_flags(self):
        import vindex.match_small as ms
        self.assertEqual(ms.missing_evidence(ms.Status('k', 'EXACT_MATCH')), '')
        self.assertIn('GQ<20', ms.missing_evidence(ms.Status('k', 'FILTERED_CALL', ('GQ<20',))))
        self.assertEqual(ms.missing_evidence(ms.Status('k', 'UNRESOLVED', ('x',), detail=dict(missing_evidence='phase'))), 'phase')


class TestOverlappingCallsMatchesABruteForceScan(unittest.TestCase):
    def test_random_spans_and_queries_agree_with_a_linear_scan(self):
        import random
        from vindex.match_small import MatchContext
        from vindex.vcfio import CallSet
        rnd = random.Random(5)

        class C:                                                                                    # the three attributes the context reads
            def __init__(self, i, pos, n):
                self.key, self.chrom, self.pos, self.ref, self.dosage, self.orig = f'1:{pos}:{"A" * n}:G#{i}', '1', pos, 'A' * n, 1, ('1', pos, 'A' * n, 'G')
                self.vtype, self.i = 'del' if n > 1 else 'snv', i

        calls = [C(i, rnd.randint(1, 5000), rnd.choice([1, 1, 1, 2, 5, 40, 400])) for i in range(600)]
        cs = CallSet(sample='S', path='x', assembly={})
        cs.small = calls
        ctx = MatchContext(cs)
        for _ in range(500):
            a = rnd.randint(1, 5200)
            b = a + rnd.choice([0, 0, 1, 10, 300])
            want = [c for c in calls if c.pos <= b and c.pos + len(c.ref) - 1 >= a]
            self.assertEqual([c.i for c in ctx.overlapping_calls('1', a, b)], [c.i for c in want], (a, b))
        self.assertEqual(ctx.overlapping_calls('2', 1, 10**6), [])

if __name__ == '__main__':
    unittest.main()
