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


if __name__ == '__main__':
    unittest.main()
