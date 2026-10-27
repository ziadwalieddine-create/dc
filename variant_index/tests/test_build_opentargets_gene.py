import json
import os
import sys
import tempfile
import unittest

from helpers import syn_reference, FIXTURES
from vindex.build import runner
from vindex.build.opentargets_gene import OpenTargetsGeneAdapter, zygosity
from vindex.index_reader import Index

sys.path.insert(0, FIXTURES)
import make_sources  # noqa: E402

SOURCES = make_sources.write_all(tempfile.mkdtemp())
OT = os.path.join(SOURCES, 'opentargets')


class FakeGenes:
    by_ensg = {'ENSG_TESTG': 1, 'ENSG_TESTH': 1}


def adapters():
    return [OpenTargetsGeneAdapter(FakeGenes(), root=OT)]


class TestZygosity(unittest.TestCase):
    def test_allelic_compositions(self):
        self.assertEqual(zygosity('Bmp4<tm1Blh>/Bmp4<+>'), 'heterozygous')
        self.assertEqual(zygosity('Bmp5<se>/Bmp5<se>'), 'homozygous')
        self.assertEqual(zygosity('Dlx5/Dlx6<tm1Levi>/Dlx5/Dlx6<tm1Levi>'), 'homozygous')          # one allele spanning two genes, present twice
        self.assertEqual(zygosity('H4c1<em1(IMPC)Ccpcz> hom early'), 'homozygous')
        self.assertEqual(zygosity('X<em1> het'), 'heterozygous')
        self.assertEqual(zygosity('A<x>/A<y>'), 'compound heterozygous or hemizygous')
        self.assertEqual(zygosity(None), 'unspecified')


class TestGeneEvidence(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.ref = syn_reference()
        cls.out = tempfile.mkdtemp()
        cls.man = runner.build(cls.out, adapters(), cls.ref)
        cls.idx = Index.open(cls.out)
        cls.ev = [r for r in cls.idx.table('evidence').to_pylist() if r['source_id'] == 'opentargets']

    def rows(self, **kw):
        out = self.ev
        for k, v in kw.items():
            out = [r for r in out if r[k] == v]
        return out

    def test_every_row_is_gene_level(self):
        self.assertTrue(all(r['level'] == 'gene' for r in self.ev))
        self.assertEqual(self.idx.table('entity_variant').num_rows, 0)                  # nothing here ever creates or touches a variant

    def test_depmap_gene_and_tissue_rows(self):
        g = self.rows(claim_type='crispr_knockout_essentiality', subject_id='ENSG_TESTG')[0]
        self.assertEqual((g['claim_value'], g['effect_assertion'], g['perturbation'], g['sample_size']), ('essential', 'effect', 'CRISPR knockout (DepMap)', 6))
        self.assertAlmostEqual(g['score'], 5 / 6)                                                                  # 3 of 4 lung screens and both skin screens
        extra = json.loads(g['extra_json'])
        self.assertEqual(extra['dependent_tissues'], [{'tissue': 'lung', 'screens': 4, 'dependent': 3}])           # skin has only two screens
        t = self.rows(claim_type='crispr_knockout_tissue_dependency')
        self.assertEqual([(r['tissue'], r['claim_value']) for r in t], [('lung', '3/4 screens')])
        h = self.rows(claim_type='crispr_knockout_essentiality', subject_id='ENSG_TESTH')[0]
        self.assertEqual((h['claim_value'], h['effect_assertion']), ('not pan-essential', 'no_effect'))

    def test_impc_rows_are_per_model_with_zygosity_and_phenotypes(self):
        hom = [r for r in self.rows(claim_type='mouse_model_phenotypes') if r['inheritance'] == 'homozygous' and r['subject_id'] == 'ENSG_TESTG'][0]
        extra = json.loads(hom['extra_json'])
        self.assertEqual(hom['perturbation'], 'mouse allele Testg<tm1A>/Testg<tm1A> (homozygous)')
        self.assertEqual((hom['claim_value'].split('; ')[0], extra['n_phenotypes'], extra['source_rows']), ('short limbs', 3, 2))
        self.assertEqual([p['n'] for p in extra['phenotype_counts']][:1], [2])
        self.assertEqual(hom['condition_ids'], ['MONDO:1', 'MONDO:2'])
        self.assertIn('not a negative result', extra['note'])
        het = [r for r in self.rows(claim_type='mouse_model_phenotypes') if r['inheritance'] == 'heterozygous'][0]
        self.assertEqual(het['lineage_group'], 'mgi_model:MGI:2')
        self.assertEqual(len(self.rows(claim_type='mouse_model_phenotypes')), 3)

    def test_unknown_gene_is_kept_unmapped(self):
        u = [r for r in self.ev if r['mapping_status'] == 'unmapped']
        self.assertEqual([(r['unmapped_key'], r['mapping_reason'], r['subject_id']) for r in u], [('ENSG_UNKNOWN', 'gene_not_in_gene_map', None)])

    def test_crispr_rows(self):
        a = self.rows(accession='crispr:c1')[0]
        self.assertEqual((a['claim_type'], a['p_value'], a['model_system'], a['lineage_group']), ('crispr_screen_hit', 0.01, 'CELL1', 'pmid:222'))
        b = self.rows(accession='crispr_screen:s1')[0]
        self.assertEqual((b['subject_id'], b['claim_num'], b['cell_type'], b['model_system'], b['p_value']), ('ENSG_TESTH', 4.5, 'Glutamatergic Neuron', 'iPSC', 0.005))
        self.assertEqual(b['perturbation'], 'CRISPR screen Genome-wide CRISPRi-v2; neuron')

    def test_burden_rows_carry_effect_sizes_and_scope(self):
        a = self.rows(accession='burden:b1')[0]
        self.assertEqual((a['odds_ratio'], a['ci_low'], a['ci_high'], a['effect_assertion'], a['inheritance'], a['n_cases']), (3.2, 2.0, 5.1, 'effect', 'monoallelic', 1200))
        self.assertAlmostEqual(a['p_value'], 2e-9)
        self.assertIn('not an effect of one allele', json.loads(a['extra_json'])['scope'])
        b = self.rows(accession='burden:b2')[0]
        self.assertEqual((b['beta'], b['ci_low'], b['effect_assertion']), (-0.2, -0.4, 'uncertain'))

    def test_a_ratio_that_is_not_an_estimate_is_kept_as_reported_and_not_stored(self):
        # found by validating the standalone build on the real file: an AstraZeneca burden row with 0 cases carrying a qualifying variant has an odds ratio of exactly 0
        from vindex.build.opentargets_gene import positive_ratio, finite_or_text
        self.assertEqual([positive_ratio(x) for x in (2.5, 0.01, 0.0, -1.0, float('inf'), float('nan'), None, 'x')], [2.5, 0.01, None, None, None, None, None, None])
        self.assertEqual((finite_or_text(0.0), finite_or_text(float('inf')), finite_or_text(float('nan'))), (0.0, 'inf', 'nan'))
        for acc, reported in (('burden:b3', 0.0), ('burden:b4', 'inf')):
            r = self.rows(accession=acc)[0]
            extra = json.loads(r['extra_json'])
            self.assertEqual((r['odds_ratio'], r['ci_low'], r['ci_high'], r['beta']), (None, None, None, None), acc)
            self.assertEqual((extra['odds_ratio_reported'], extra['odds_ratio_not_stored']), (reported, 'not a positive finite number'))
        b3 = self.rows(accession='burden:b3')[0]
        self.assertEqual((b3['effect_assertion'], json.loads(b3['extra_json'])['odds_ratio_interval_reported'], json.loads(b3['extra_json'])['cases_with_qualifying_variants']), ('effect', [0.0, 1.5], 0))
        self.assertAlmostEqual(b3['p_value'], 4.359e-12)                                                                         # the evidence the row does carry is kept
        ok = self.rows(accession='burden:b1')[0]
        self.assertNotIn('odds_ratio_reported', json.loads(ok['extra_json']))

    def test_a_missing_project_or_cohort_is_unspecified_in_the_lineage_key_not_the_text_none(self):
        got = {a: self.rows(accession=a)[0]['lineage_group'] for a in ('burden:b1', 'burden:b3', 'burden:b4')}
        self.assertEqual(got, {'burden:b1': 'burden_cohort:Genebass:UKB', 'burden:b3': 'burden_cohort:unspecified:INTERVAL', 'burden:b4': 'burden_cohort:Genebass:unspecified'})
        self.assertFalse([r for r in self.ev if 'None' in (r['lineage_group'] or '')])

    def test_ids_unique_and_rebuild_idempotent(self):
        ids = [r['evidence_id'] for r in self.ev]
        self.assertEqual(len(ids), len(set(ids)))
        out2 = tempfile.mkdtemp()
        runner.build(out2, adapters(), self.ref)
        with open(os.path.join(self.out, 'table_digests.json')) as f1, open(os.path.join(out2, 'table_digests.json')) as f2:
            self.assertEqual(json.load(f1), json.load(f2))


if __name__ == '__main__':
    unittest.main()
