import os
import shutil
import sys
import tempfile
import unittest

import pyarrow.compute as pc

from helpers import syn_reference, FIXTURES
from vindex.build import runner
from vindex.build.clinvar_vcf import ClinVarVcfAdapter
from vindex.build.clingen import ClinGenErepoAdapter, ClinGenGeneRegionAdapter
from vindex.index_reader import Index

sys.path.insert(0, FIXTURES)
import make_sources  # noqa: E402

SOURCES = make_sources.write_all(tempfile.mkdtemp())


class FakeGenes:
    """Stands in for GeneMap: no Open Targets data are needed for the unit test."""
    rows = []
    ids = {'HGNC:20': 'ENSG_AARS1', 'HGNC:1100': 'ENSG_BRCA1'}
    syms = {'AARS': 'ENSG_AARS1', 'BRCA1': 'ENSG_BRCA1'}

    def resolve(self, hgnc=None, symbol=None, ncbi=None):
        for key, d in ((hgnc, self.ids), (symbol, self.syms)):
            if key in d:
                return d[key], 'hgnc' if d is self.ids else 'symbol'
        return None, 'unresolved'


class MiniClinVar(ClinVarVcfAdapter):
    def run(self, ctx, limit=None):
        return super().run(ctx, path=os.path.join(SOURCES, 'mini_clinvar.vcf'))


class MiniErepo(ClinGenErepoAdapter):
    def __init__(self):
        super().__init__(os.path.join(SOURCES, 'erepo_variant_classifications.tsv'))


class MiniGeneRegion(ClinGenGeneRegionAdapter):
    def __init__(self):
        super().__init__(FakeGenes(), root=SOURCES)


class TestClinGen(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.ref = syn_reference()
        cls.out = tempfile.mkdtemp()
        # ERepo registers its interest in prepare(), before any adapter runs, and must be listed after ClinVar
        cls.man = runner.build(cls.out, [MiniClinVar(), MiniErepo(), MiniGeneRegion()], cls.ref)
        cls.idx = Index.open(cls.out)
        cls.ev = [r for r in cls.idx.table('evidence').to_pylist() if r['source_id'] == 'clingen']

    @classmethod
    def tearDownClass(cls):
        shutil.rmtree(cls.out, ignore_errors=True)

    def rows(self, **kw):
        return [r for r in self.ev if all(r[k] == v for k, v in kw.items())]

    def test_erepo_variant_resolution_and_unmapped_rows(self):
        main = {r['accession']: r for r in self.rows(claim_type='vcep_classification')}
        base = self.ref.fetch('1', 402, 402)
        self.assertEqual(main['ERepo:uuid-1']['mapping_status'], 'mapped')
        self.assertTrue(main['ERepo:uuid-1']['subject_id'].startswith('1:402:'))
        self.assertEqual(main['ERepo:uuid-5']['subject_id'], main['ERepo:uuid-1']['subject_id'])       # resolved through the Allele Registry id alone
        self.assertEqual((main['ERepo:uuid-3']['mapping_status'], main['ERepo:uuid-3']['mapping_reason']), ('unmapped', 'no_clinvar_or_allele_registry_match'))
        self.assertEqual((main['ERepo:uuid-4']['mapping_status'], main['ERepo:uuid-4']['mapping_reason']), ('unmapped', 'hgvs_only'))
        self.assertEqual(main['ERepo:uuid-3']['unmapped_key'], 'ClinVar:99999')
        self.assertEqual(len(main), 5)                                                                   # nothing dropped

    def test_criteria_land_on_their_axes_with_met_and_not_met_kept(self):
        crit = {(r['accession'], r['claim_type']): r for r in self.ev if r['claim_type'].startswith('acmg_criterion:')}
        self.assertEqual(crit[('ERepo:uuid-1', 'acmg_criterion:PVS1')]['axis'], 'mechanistic_confidence')
        self.assertEqual(crit[('ERepo:uuid-1', 'acmg_criterion:PS3_Moderate')]['axis'], 'empirical_functional_support')
        self.assertEqual(crit[('ERepo:uuid-1', 'acmg_criterion:PM2_Supporting')]['axis'], 'applicability_context')
        self.assertEqual(crit[('ERepo:uuid-1', 'acmg_criterion:PVS1')]['effect_assertion'], 'effect')
        self.assertEqual((crit[('ERepo:uuid-1', 'acmg_criterion:PS4')]['claim_value'], crit[('ERepo:uuid-1', 'acmg_criterion:PS4')]['effect_assertion']), ('not_met', 'not_applicable'))
        self.assertEqual(crit[('ERepo:uuid-2', 'acmg_criterion:BP4')]['effect_assertion'], 'no_effect')
        self.assertEqual(len({r['lineage_group'] for r in self.ev if r['accession'] == 'ERepo:uuid-1'}), 1)   # components of one classification are not independent
        # an expert-panel classification that ClinVar carries is the same assertion as the ClinVar record of that variation, so it shares the lineage of the ClinVar rows
        lin = {r['accession']: r['lineage_group'] for r in self.ev if r['claim_type'] == 'vcep_classification'}
        self.assertEqual((lin['ERepo:uuid-1'], lin['ERepo:uuid-3']), ('clinvar:variation:1001', 'clinvar:variation:99999'))
        self.assertEqual(lin['ERepo:uuid-5'], 'clingen_erepo:uuid-5')                                          # no ClinVar variation id: its own lineage

    def test_gene_and_region_evidence_stay_on_genes_and_regions(self):
        v = self.rows(claim_type='gene_disease_validity')
        self.assertEqual({(r['accession'], r['mapping_status']) for r in v}, {('CGGV:assertion_1', 'mapped'), ('CGGV:assertion_2', 'unmapped'), ('CGGV:assertion_3', 'mapped')})
        self.assertEqual({r['claim_value']: r['effect_assertion'] for r in v}, {'Definitive': 'effect', 'Limited': 'uncertain', 'Disputed': 'conflicting'})
        self.assertTrue(all(r['level'] == 'gene' for r in v))
        self.assertEqual({r['accession']: r['inheritance'] for r in v}, {'CGGV:assertion_1': 'Autosomal dominant inheritance', 'CGGV:assertion_2': 'Autosomal recessive inheritance',
                                                                         'CGGV:assertion_3': 'Autosomal dominant inheritance'})        # the file's codes AD and AR, spelled as the ERepo rows spell them
        hi = self.rows(claim_type='haploinsufficiency_score', level='gene')
        self.assertEqual(sorted((r['subject_id'] or r['unmapped_key'], r['claim_num'], r['effect_assertion']) for r in hi), [('ENSG_BRCA1', 3.0, 'effect'), ('NOSUCHGENE', 30.0, 'no_effect')])
        reg = [r for r in self.ev if r['level'] == 'region']
        self.assertEqual({(r['subject_id'], r['mapping_status'], r['mapping_reason']) for r in reg},
                         {('ISCA-1', 'mapped', None), (None, 'unmapped', 'region_location_not_parsed')})
        regions = {r['region_id']: r for r in self.idx.table('entity_region').to_pylist()}
        self.assertEqual((regions['ISCA-1']['start'], regions['ISCA-1']['end']), (1000, 2000))
        self.assertTrue(any(r['kind'] == 'clingen_recurrent_cnv' and r['start'] == 1000 for r in regions.values()))          # BED start converted to 1-based

    def test_wrong_adapter_order_is_refused_not_silently_degraded(self):
        with self.assertRaises(ValueError):
            runner.build(tempfile.mkdtemp(), [MiniErepo(), MiniClinVar()], self.ref)

    def test_actionability_is_gene_condition_level_and_variant_curations_stay_unmapped(self):
        a = self.rows(claim_type='clinical_actionability_assertion_adult')
        gene = [r for r in a if r['level'] == 'gene']
        var = [r for r in a if r['level'] == 'variant']
        self.assertEqual([(r['subject_id'], r['claim_value'], r['effect_assertion']) for r in gene], [('ENSG_BRCA1', 'Strong Actionability', 'effect')])
        self.assertEqual([(r['subject_id'], r['mapping_status'], r['unmapped_key']) for r in var], [(None, 'unmapped', 'NM_000059.4(BRCA2):c.1A>G')])
        self.assertTrue(self.rows(claim_type='clinical_actionability_scoring_pediatric'))

    def test_gene_evidence_is_never_returned_as_variant_evidence(self):
        self.assertFalse([r for r in self.ev if r['level'] == 'variant' and r['claim_type'] in ('gene_disease_validity', 'haploinsufficiency_score', 'triplosensitivity_score')])
        base = self.ref.fetch('1', 402, 402)
        alt = next(x for x in 'ACGT' if x != base and x != {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}[base])
        got = self.idx.evidence_for_variants([f'1:402:{base}:{alt}'])
        self.assertTrue(got)
        self.assertTrue(all(r['level'] == 'variant' for rows in got.values() for r in rows))


if __name__ == '__main__':
    unittest.main()
