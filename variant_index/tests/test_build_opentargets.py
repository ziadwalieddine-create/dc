import json
import os
import sys
import tempfile
import unittest

from helpers import syn_reference, FIXTURES
from vindex.build import runner
from vindex.build.opentargets_genetics import OpenTargetsGeneticsAdapter, p_value, split_variant_id
from vindex.index_reader import Index

sys.path.insert(0, FIXTURES)
import make_sources  # noqa: E402

SOURCES = make_sources.write_all(tempfile.mkdtemp())
KEYS = make_sources.write_all.opentargets
OT = os.path.join(SOURCES, 'opentargets')


class FakeGenes:
    by_ensg = {'ENSG_TESTG': 1, 'ENSG_OTHER': 1}


def adapters(**kw):
    return [OpenTargetsGeneticsAdapter(FakeGenes(), root=OT, **kw)]


class TestHelpers(unittest.TestCase):
    def test_p_values_and_variant_ids(self):
        self.assertAlmostEqual(p_value(2.0, -7), 2e-7)
        self.assertIsNone(p_value(1.0, -400))          # underflow stays text only
        self.assertIsNone(p_value(None, -3))
        self.assertIsNone(p_value(12.0, 1))             # not a probability
        self.assertEqual(split_variant_id('1_100_A_G'), ['1', '100', 'A', 'G'])
        self.assertIsNone(split_variant_id('OTVAR_deadbeef'))
        self.assertIsNone(split_variant_id(None))


class TestAncestry(unittest.TestCase):
    def test_the_label_holds_reported_names_only_and_the_sizes_travel_separately(self):
        from vindex.build.opentargets_genetics import ancestry_of
        names, sizes, basis = ancestry_of(dict(discoverySamples=[dict(ancestry='European', sampleSize=100), dict(ancestry='European', sampleSize=50), dict(ancestry='NR', sampleSize=7)]))
        self.assertEqual((names, sizes, basis), ('European', {'European': 150, 'NR': 7}, 'discovery_samples'))               # 'NR' is no ancestry; its size is kept
        names, sizes, basis = ancestry_of(dict(discoverySamples=[dict(ancestry='NR', sampleSize=7)], ldPopulationStructure=[dict(ldPopulation='nfe', relativeSampleSize=0.8), dict(ldPopulation='afr', relativeSampleSize=0.2)]))
        self.assertEqual((names, basis), ('afr; nfe', 'ld_population_structure'))                                              # nothing reported: the LD populations, labelled as such
        self.assertEqual(ancestry_of(dict(discoverySamples=[dict(ancestry='NR', sampleSize=7)])), (None, {'NR': 7}, 'discovery_samples'))
        self.assertEqual(ancestry_of(dict()), (None, {}, 'discovery_samples'))

    def test_a_missing_relative_size_or_sample_size_is_kept_as_missing_not_a_crash(self):
        # found by the full build: an LD population without a relative sample size (the fallback is reached for studies whose only reported ancestry is 'NR')
        from vindex.build.opentargets_genetics import ancestry_of
        names, sizes, basis = ancestry_of(dict(discoverySamples=[dict(ancestry='NR', sampleSize=None)], ldPopulationStructure=[dict(ldPopulation='nfe', relativeSampleSize=None), dict(ldPopulation='afr', relativeSampleSize=0.25)]))
        self.assertEqual((names, sizes, basis), ('afr; nfe', {'nfe': None, 'afr': 0.25}, 'ld_population_structure'))
        self.assertEqual(ancestry_of(dict(discoverySamples=[dict(ancestry='African', sampleSize=None)])), ('African', {'African': 0}, 'discovery_samples'))
        self.assertEqual(ancestry_of(dict(discoverySamples=None, ldPopulationStructure=[dict(ldPopulation=None, relativeSampleSize=None)])), (None, {}, 'discovery_samples'))


class TestOpenTargets(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.ref = syn_reference()
        cls.out = tempfile.mkdtemp()
        cls.man = runner.build(cls.out, adapters(), cls.ref)
        cls.idx = Index.open(cls.out)
        cls.ev = [r for r in cls.idx.table('evidence').to_pylist() if r['source_id'] == 'opentargets']
        cls.st = cls.man['adapters'][0]['stats']

    def rows(self, **kw):
        out = self.ev
        for k, v in kw.items():
            out = [r for r in out if r[k] == v]
        return out

    def test_thresholds_split_gwas_and_qtl_and_everything_is_accounted(self):
        s = self.st
        self.assertEqual((s['credible_sets'], s['members_total'], s['members_integrated']), (3, 9, 5))
        # SL1 (gwas, pip >= 0.1): 402, 410, rejected 440, 430 -> 4; SL2 (qtl, pip >= 0.5): 402 -> 1; SL3: none reaches 0.1
        self.assertEqual(sum(v for k, v in s.items() if k.startswith('integrated_rows:')), 5)
        self.assertEqual(self.rows(accession='SL3', level='variant'), [])

    def test_gwas_member_row(self):
        r = self.rows(accession='SL1', subject_id=KEYS['k402'])[0]
        self.assertEqual((r['level'], r['axis'], r['evidence_kind'], r['claim_type'], r['effect_assertion']), ('variant', 'fine_mapping_probability', 'statistical_association', 'gwas_fine_mapping_pip', 'not_applicable'))
        self.assertEqual((r['pip'], r['beta'], r['se'], r['effect_allele'], r['other_allele'], r['eaf']), (0.8, 0.5, 0.1, KEYS['k402'].split(':')[3], KEYS['k402'].split(':')[2], 0.25))
        self.assertAlmostEqual(r['p_value'], 1e-9)
        self.assertEqual((r['lineage_group'], r['sample_size'], r['n_cases'], r['ancestry']), ('pmid:111', 1000, 300, 'European'))                  # the label holds the name only; the size travels in extra_json
        self.assertEqual((json.loads(r['extra_json'])['ancestry_sizes'], json.loads(r['extra_json'])['ancestry_basis']), ({'European': 1000}, 'discovery_samples'))
        self.assertEqual(r['condition_ids'], ['EFO_0000001', 'MONDO_0000001'])
        extra = json.loads(r['extra_json'])
        self.assertEqual((extra['is_lead'], extra['method'], extra['beta_reference']), (True, 'SuSiE-inf', 'alternate allele of the variant id'))
        non_lead = self.rows(accession='SL1', subject_id=KEYS['k410'])[0]
        self.assertIsNone(non_lead['eaf'])
        self.assertFalse(json.loads(non_lead['extra_json'])['is_lead'])

    def test_member_below_threshold_is_not_integrated_and_a_rejected_allele_is_kept_as_rejected(self):
        self.assertEqual([r for r in self.rows(accession='SL1') if r['unmapped_key'] and r['unmapped_key'].startswith('1_420_')], [])
        rej = [r for r in self.rows(accession='SL1') if r['mapping_status'] != 'mapped' and r['level'] == 'variant']
        self.assertEqual([(r['mapping_status'], r['mapping_reason']) for r in rej], [('rejected', 'ref_mismatch')])
        self.assertTrue(rej[0]['unmapped_key'].startswith('1_440_'))

    def test_underflowing_p_value_is_text_only(self):
        r = self.rows(accession='SL1', subject_id=KEYS['k430'])[0]
        self.assertIsNone(r['p_value'])
        self.assertEqual(json.loads(r['extra_json'])['p_text'], '1.0e-400')

    def test_qtl_row_carries_tissue_gene_and_project_lineage(self):
        r = self.rows(accession='SL2', subject_id=KEYS['k402'])[0]
        self.assertEqual((r['claim_type'], r['tissue'], r['gene_ids'], r['lineage_group'], r['beta']), ('qtl_fine_mapping_pip', 'lymphoblastoid cell line', ['ENSG_TESTG'], 'qtl_project:GTEx-v10', -0.8))
        self.assertEqual(r['condition_text'], 'eqtl of ENSG_TESTG in lymphoblastoid cell line')
        self.assertEqual(self.rows(accession='SL2', subject_id=KEYS['k410']), [])                    # pip 0.4 < 0.5
        vg = [x for x in self.idx.table('variant_gene').to_pylist() if x['relation'] == 'ot_eqtl_credible_set']
        self.assertEqual([(x['variant_key'], x['gene_id']) for x in vg], [(KEYS['k402'], 'ENSG_TESTG')])

    def test_l2g_rows_gene_level_thresholded_and_unmapped_gene_kept(self):
        g = self.rows(claim_type='l2g_score')
        self.assertEqual(sorted((r['accession'], r['subject_id'] or r['unmapped_key'], r['mapping_status']) for r in g), [('SL1', 'ENSG_TESTG', 'mapped'), ('SL3', 'ENSG_UNKNOWN', 'unmapped')])
        t = [r for r in g if r['accession'] == 'SL1'][0]
        self.assertEqual((t['level'], t['evidence_kind'], t['score'], t['lineage_group']), ('gene', 'computational_prediction', 0.9, 'pmid:111'))
        self.assertEqual(json.loads(t['extra_json'])['lead_key'], KEYS['k402'])
        vg = [x for x in self.idx.table('variant_gene').to_pylist() if x['relation'] == 'ot_l2g_lead_variant']
        self.assertEqual([(x['variant_key'], x['gene_id'], x['consequence']) for x in vg], [(KEYS['k402'], 'ENSG_TESTG', '0.900')])
        self.assertEqual(self.st['l2g_rows_at_threshold'], 2)                                           # 0.3 is below the 0.5 threshold

    def test_study_entities_for_studies_with_integrated_rows(self):
        stu = {r['study_id']: r for r in self.idx.table('entity_study').to_pylist() if r['source_id'] == 'opentargets'}
        self.assertEqual(sorted(stu), ['FINNGEN_R12_TEST', 'GCST000001', 'GTEX_EQTL_X'])
        self.assertEqual((stu['GCST000001']['pubmed'], stu['GCST000001']['n_total'], stu['GCST000001']['ancestry']), ('111', 1000, 'European'))
        self.assertEqual(stu['GTEX_EQTL_X']['biosample'], 'lymphoblastoid cell line')
        self.assertEqual((stu['GCST000001']['cohorts'], stu['FINNGEN_R12_TEST']['cohorts'], stu['GTEX_EQTL_X']['cohorts']), ('TestCohort', 'FinnGen', None))     # the Catalog's 'NR' is no cohort
        self.assertEqual(stu['FINNGEN_R12_TEST']['lineage_group'], 'ot_project:FINNGEN_R12:FINNGEN_R12_TEST')

    def test_studies_described_by_another_source_are_not_duplicated(self):
        class Gwas:
            source_id = 'gwas_catalog'

            def run(self, ctx, limit=None):
                ctx.offer('study', 'GCST000001', 'GCST000001')
                return {}
        ot = OpenTargetsGeneticsAdapter(FakeGenes(), root=OT)
        out = tempfile.mkdtemp()
        man = runner.build(out, [Gwas(), ot], self.ref)
        stu = [r['study_id'] for r in Index.open(out).table('entity_study').to_pylist()]
        self.assertNotIn('GCST000001', stu)
        self.assertEqual(man['adapters'][1]['stats']['studies_described_by_another_source'], 1)

    def test_thresholds_are_parameters_and_idempotent(self):
        out = tempfile.mkdtemp()
        man = runner.build(out, adapters(gwas_pip=0.01, qtl_pip=0.3), self.ref)
        self.assertEqual(man['adapters'][0]['stats']['members_integrated'], 9)       # every member now: SL1 5 (0.05 kept), SL2 2 (0.4 kept), SL3 2
        out2 = tempfile.mkdtemp()
        runner.build(out2, adapters(), self.ref)
        with open(os.path.join(self.out, 'table_digests.json')) as f1, open(os.path.join(out2, 'table_digests.json')) as f2:
            self.assertEqual(json.load(f1), json.load(f2))


if __name__ == '__main__':
    unittest.main()
