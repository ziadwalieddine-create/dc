import unittest

from helpers import syn_reference, vcf_text
from evidence_fixtures import ev, index_from
from vindex import schema
from vindex.vcfio import SVCall
from vindex.index_reader import Index
from vindex.match_sv import compare, match_svs, SvPolicy


def sv(start, end, typ='DEL', cipos=(0, 0), ciend=(0, 0), precise=True, cn=None, dosage=1, svlen=None):
    return SVCall(f'{typ}:{start}-{end}', '1', start, end, typ, svlen if svlen is not None else (-(end - start) if typ == 'DEL' else end - start), cipos, ciend,
                  precise, cn, '0/1', dosage, 2, 50, 'PASS', True, 0, f'<{typ}>')


def event(eid, start, end, typ='DEL', **kw):
    d = dict(event_id=eid, source_id='gnomad_sv', release_id='r', accession=eid, event_type=typ, chrom='1', start=start, end=end, start_ci_lo=0, start_ci_hi=0,
             end_ci_lo=0, end_ci_hi=0, svlen=-(end - start) if typ == 'DEL' else end - start, copy_number=None, precise=True, population_af=0.01, n_samples=1000)
    d.update(kw)
    return d


class TestSvEquivalence(unittest.TestCase):
    def test_same_breakpoints_are_equivalent(self):
        m = compare(sv(1000, 1500), event('A', 1000, 1500))
        self.assertTrue(m.equivalent)
        self.assertEqual(m.relation, 'equivalent_breakpoints')

    def test_small_breakpoint_jitter_within_tolerance(self):
        m = compare(sv(1000, 1500), event('A', 1030, 1480))
        self.assertTrue(m.equivalent, m)

    def test_half_overlap_is_not_equivalent(self):
        m = compare(sv(1000, 1500), event('B', 1250, 1750))
        self.assertFalse(m.equivalent)
        self.assertEqual(m.relation, 'partial_overlap')
        self.assertAlmostEqual(m.metrics['reciprocal_overlap'], 0.5, places=2)

    def test_containing_event_is_not_equivalent(self):
        m = compare(sv(1000, 1500), event('C', 900, 1600))
        self.assertFalse(m.equivalent)
        self.assertEqual(m.relation, 'partial_overlap')

    def test_nested_small_event_is_overlap_only(self):
        m = compare(sv(1000, 1500), event('D', 1100, 1150))
        self.assertEqual((m.equivalent, m.relation), (False, 'overlap_only'))

    def test_type_mismatch_with_identical_coordinates(self):
        m = compare(sv(1000, 1500), event('E', 1000, 1500, 'DUP'))
        self.assertEqual((m.equivalent, m.relation), (False, 'type_mismatch'))

    def test_deletion_matches_copy_number_loss_and_not_gain(self):
        self.assertTrue(compare(sv(1000, 1500), event('F', 1000, 1500, 'CNV_LOSS')).equivalent)
        self.assertFalse(compare(sv(1000, 1500), event('G', 1000, 1500, 'CNV_GAIN')).equivalent)

    def test_imprecise_calls_use_reciprocal_overlap(self):
        wide = sv(1000, 1500, precise=False)
        self.assertEqual(compare(wide, event('H', 1100, 1560)).relation, 'partial_overlap')          # RO 400/500 = 0.80 < 0.90
        self.assertEqual(compare(wide, event('I', 1020, 1520)).relation, 'equivalent_breakpoints')   # both ends within 50 bp
        big = sv(100000, 200000, precise=False)
        self.assertEqual(compare(big, event('I2', 105000, 205000)).relation, 'equivalent_overlap')   # ends 5 kb apart, RO 0.95, no breakpoint resolution
        self.assertEqual(compare(sv(100000, 200000), event('I3', 105000, 205000)).relation, 'partial_overlap')   # a precise call needs matching breakpoints
        self.assertEqual(compare(sv(1000, 1500, precise=True), event('H2', 1100, 1600)).relation, 'partial_overlap')   # precise calls need breakpoints, not just overlap

    def test_large_confidence_interval_cannot_make_everything_equivalent(self):
        s = sv(1000, 1500, cipos=(-400, 400), ciend=(-400, 400))
        m = compare(s, event('J', 1300, 1800))
        self.assertFalse(m.equivalent)
        self.assertLessEqual(m.metrics['tol_start'], 50)            # 10% of a 500 bp event is 50 bp

    def test_dosage_rules(self):
        homdel = event('K', 1000, 1500, copy_number=0)
        self.assertEqual(compare(sv(1000, 1500, cn=1), homdel).relation, 'dosage_mismatch')
        self.assertEqual(compare(sv(1000, 1500, cn=None), homdel).relation, 'dosage_unresolved')
        self.assertTrue(compare(sv(1000, 1500, cn=0), homdel).equivalent)

    def test_generic_cnv_without_copy_number_is_unresolved(self):
        m = compare(sv(1000, 1500, 'CNV'), event('L', 1000, 1500, 'CNV_LOSS'))
        self.assertEqual(m.relation, 'dosage_unresolved')
        self.assertTrue(compare(sv(1000, 1500, 'CNV', cn=1), event('M', 1000, 1500, 'CNV_LOSS')).equivalent)

    def test_insertions(self):
        i = sv(2000, 2000, 'INS', svlen=300)
        self.assertTrue(compare(i, event('N', 2010, 2010, 'INS', svlen=290)).equivalent)
        self.assertFalse(compare(i, event('O', 2010, 2010, 'INS', svlen=100)).equivalent)
        self.assertFalse(compare(i, event('P', 2400, 2400, 'INS', svlen=300)).equivalent)

    def test_breakends_are_unsupported(self):
        b = sv(1000, 1000, 'BND')
        self.assertEqual(compare(b, event('Q', 1000, 1500)).relation, 'unsupported')


class TestSvInterpretation(unittest.TestCase):
    def test_equivalent_and_refused_events_and_context(self):
        events = [event('A', 1000, 1500), event('B', 1250, 1750), event('E', 1000, 1500, 'DUP')]
        genes = [dict(gene_id='G1', symbol='GENEA', hgnc_id='', chrom='1', start=1100, end=1400, strand='+', mane_select='', biotype='protein_coding', release_id='r'),
                 dict(gene_id='G2', symbol='GENEB', hgnc_id='', chrom='1', start=1450, end=2000, strand='+', mane_select='', biotype='protein_coding', release_id='r')]
        regions = [dict(region_id='R1', kind='clingen_dosage_region', chrom='1', start=900, end=1600, name='larger', gene_ids=['G1'], release_id='r'),
                   dict(region_id='R2', kind='clingen_dosage_region', chrom='1', start=1100, end=1400, name='inside', gene_ids=['G1'], release_id='r')]
        idx = Index.from_tables(entity_event=schema.rows_to_table('entity_event', events), entity_gene=schema.rows_to_table('entity_gene', genes),
                                entity_region=schema.rows_to_table('entity_region', regions))

        class CS:
            svs = [sv(1000, 1500)]
        out = match_svs(CS, idx)[0]
        eq = [m.event_id for m in out.matches if m.equivalent]
        self.assertEqual(eq, ['A'])
        refused = {m.event_id: m.relation for m in out.matches if not m.equivalent}
        self.assertEqual(refused, {'B': 'partial_overlap', 'E': 'type_mismatch'})
        rel = {c['id']: c['relation'] for c in out.context}
        self.assertEqual(rel, {'G1': 'encompasses_gene', 'G2': 'partially_overlaps_gene', 'R1': 'partially_overlaps_region', 'R2': 'encompasses_region'})
        self.assertTrue(all(c['applies_as'].endswith('context_only') for c in out.context))

    def test_a_homozygous_reference_structural_genotype_is_kept_but_never_matched(self):
        class CS:
            svs = [sv(1000, 1500, dosage=0), sv(1000, 1500, dosage=1)]
        out = match_svs(CS, Index.from_tables(entity_event=schema.rows_to_table('entity_event', [event('A', 1000, 1500)])))
        self.assertEqual(len(out), 2)                                                              # both records are accounted for
        self.assertEqual((out[0].matches, out[0].context), ([], []))
        self.assertIn('homozygous_reference', out[0].unresolved[0])
        self.assertEqual([m.event_id for m in out[1].matches if m.equivalent], ['A'])               # the same record, carried, does match

    def test_dosage_unknown_is_reported(self):
        class CS:
            svs = [sv(1000, 1500, dosage=None)]
        out = match_svs(CS, Index.from_tables(entity_event=schema.rows_to_table('entity_event', [event('A', 1000, 1500)])))[0]
        self.assertIn('dosage unknown', out.unresolved[0])


if __name__ == '__main__':
    unittest.main()
