"""The known-locus engine on a small in-memory index: a fact that holds passes, and every way the index can contradict or fail to hold it is reported as FAIL or ABSENT, never PASS."""
import json
import os
import subprocess
import sys
import tempfile
import unittest

from helpers import ROOT
from vindex import known_loci, schema
from vindex.build.common import evidence_row
from vindex.index_reader import Index

KEY = '10:100:C:T'


def ev(source, acc, claim_type, claim_value, subject=KEY, level='variant', axis='association_replication', kind='statistical_association', assertion='effect', **kw):
    return evidence_row('rel', source, acc, 'L1', f'{acc}|{claim_type}|{claim_value}', level, subject, axis, kind, claim_type, claim_value, assertion, discriminator=json.dumps(kw, sort_keys=True, default=str), **kw)


def gwas(acc, allele, other, odds=None, beta=None, trait='Type 2 diabetes'):
    return ev('gwas_catalog', acc, 'gwas_effect_size', 'odds_ratio' if odds else 'beta', axis='effect_magnitude_direction', effect_allele=allele, other_allele=other,
              odds_ratio=odds, beta=beta, condition_text=trait)


def index(evidence=(), xref=(), haplotypes=(), definitions=(), genes=(), regions=(), assays=()):
    T = schema.rows_to_table
    return Index.from_tables(evidence=T('evidence', list(evidence)), xref=T('xref', list(xref)), entity_haplotype=T('entity_haplotype', list(haplotypes)),
                             haplotype_definition=T('haplotype_definition', list(definitions)), entity_gene=T('entity_gene', list(genes)), entity_region=T('entity_region', list(regions)),
                             entity_assay=T('entity_assay', list(assays)))


def xref(rsid, key=KEY, relation='clinvar_RS'):
    return dict(namespace='dbsnp', identifier=rsid, entity_level='variant', entity_id=key, release_id='rel', relation=relation)


def spec(*checks, **locus):
    return dict(loci=[dict(dict(id='L', title='t', basis='b', rsid='rs1', key=KEY, checks=list(checks)), **locus)])


def status(index_, check, **locus):
    return known_loci.run(index_, spec(check, **locus))['loci'][0]['checks'][0]['status']


class TestXrefKey(unittest.TestCase):
    check = dict(id='x', type='xref_key')

    def test_the_stated_key_passes(self):
        self.assertEqual(status(index(xref=[xref('rs1')]), self.check), 'PASS')

    def test_another_key_at_the_same_site_is_still_the_same_position_but_a_different_site_fails(self):
        self.assertEqual(status(index(xref=[xref('rs1'), xref('rs1', '10:100:C:G')]), self.check), 'PASS')                         # multiallelic site: one position
        self.assertEqual(status(index(xref=[xref('rs1'), xref('rs1', '10:200:C:T', 'gwas_catalog_snps')]), self.check), 'FAIL')    # two sources place the rsID at two sites
        self.assertEqual(status(index(xref=[xref('rs1', '10:200:C:T')]), self.check), 'FAIL')                                        # the key is not among them

    def test_an_rsid_the_index_does_not_know_is_absent(self):
        self.assertEqual(status(index(xref=[xref('rs2')]), self.check), 'ABSENT')


class TestOrientation(unittest.TestCase):
    check = dict(id='o', type='orientation', where={'source_id': 'gwas_catalog', 'claim_type': 'gwas_effect_size'}, regex={'condition_text': 'type 2 diabetes'},
                 not_regex={'condition_text': 'gestational'}, risk_allele='T', other_allele='C', min_rows=3, min_fraction=0.9)

    def rows(self, flipped=0):
        good = [gwas(f'A{i}', 'T', 'C', odds=1.3 + i / 100) for i in range(5)] + [gwas(f'B{i}', 'C', 'T', odds=0.8) for i in range(4)] + [gwas('C0', 'T', 'C', beta=0.1)]
        bad = [gwas(f'F{i}', 'T', 'C', odds=0.7) for i in range(flipped)]
        return good + bad

    def test_consistent_orientation_passes_and_a_flipped_allele_fails(self):
        self.assertEqual(status(index(self.rows()), self.check), 'PASS')
        res = known_loci.run(index(self.rows(flipped=3)), spec(self.check))['loci'][0]['checks'][0]
        self.assertEqual((res['status'], res['discordant']), ('FAIL', 3))
        self.assertEqual(len(res['discordant_examples']), 3)

    def test_excluded_traits_and_uninformative_rows_are_not_counted(self):
        rows = self.rows() + [gwas('G0', 'T', 'C', odds=0.5, trait='Gestational type 2 diabetes')] + [gwas('N0', 'T', 'C', odds=1.0), gwas('N1', 'A', 'G', odds=0.5)]
        res = known_loci.run(index(rows), spec(self.check))['loci'][0]['checks'][0]
        self.assertEqual((res['status'], res['discordant'], res['uninformative']), ('PASS', 0, 2))

    def test_too_few_informative_rows_is_absent_not_pass(self):
        self.assertEqual(status(index(self.rows()[:2]), self.check), 'ABSENT')
        self.assertEqual(status(index([]), self.check), 'ABSENT')
        self.assertEqual(status(index([gwas('Z0', 'T', 'C', odds=1.2, trait='Height')]), self.check), 'ABSENT')

    def test_beta_and_odds_ratio_are_read_the_same_way(self):
        self.assertEqual((known_loci.raises_risk(dict(odds_ratio=1.2)), known_loci.raises_risk(dict(odds_ratio=0.9)), known_loci.raises_risk(dict(beta=-0.1)),
                          known_loci.raises_risk(dict(beta=0.2)), known_loci.raises_risk(dict(odds_ratio=1.0)), known_loci.raises_risk(dict())), (1, -1, -1, 1, 0, 0))


class TestEvidence(unittest.TestCase):
    check = dict(id='e', type='evidence', where={'source_id': 'clinvar', 'axis': 'clinical_validity_actionability'}, regex={'claim_value': '^pathogenic'}, min_rows=1)

    def clinvar(self, value, subject=KEY):
        return ev('clinvar', 'VCV1', 'germline_classification', value, subject=subject, axis='clinical_validity_actionability', kind='curated_expert_classification')

    def test_a_matching_classification_passes(self):
        self.assertEqual(status(index([self.clinvar('Pathogenic')]), self.check), 'PASS')
        self.assertEqual(status(index([self.clinvar('Pathogenic/Likely pathogenic')]), self.check), 'PASS')

    def test_a_different_classification_fails_and_no_evidence_at_all_is_absent(self):
        self.assertEqual(status(index([self.clinvar('Benign')]), self.check), 'FAIL')
        self.assertEqual(status(index([self.clinvar('Pathogenic', subject='10:999:A:G')]), self.check), 'ABSENT')

    def test_numeric_bounds_text_and_list_selections(self):
        rows = [ev('x', 'a', 'k', 'v', axis='effect_magnitude_direction', score=0.3, condition_ids=['EFO:1'], condition_text='Warfarin dose')]
        idx = index(rows)
        self.assertEqual(status(idx, dict(id='b', type='evidence', bounds={'score': [None, 0.5]})), 'PASS')
        self.assertEqual(status(idx, dict(id='b', type='evidence', bounds={'score': [0.4, None]})), 'FAIL')
        self.assertEqual(status(idx, dict(id='t', type='evidence', text_regex='warfarin')), 'PASS')
        self.assertEqual(status(idx, dict(id='t', type='evidence', text_regex='clopidogrel')), 'FAIL')
        self.assertEqual(status(idx, dict(id='l', type='evidence', list_any={'condition_ids': ['EFO:1', 'EFO:2']})), 'PASS')
        self.assertEqual(status(idx, dict(id='l', type='evidence', list_any={'condition_ids': ['EFO:2']})), 'FAIL')
        self.assertEqual(status(idx, dict(id='n', type='evidence', regex={'condition_text': 'dose'}, not_regex={'condition_text': 'warfarin'})), 'FAIL')

    def test_gene_level_evidence_is_found_through_the_symbol(self):
        gene = dict(gene_id='ENSG1', symbol='GENE1', hgnc_id='', chrom='1', start=1, end=2, strand='+', mane_select='', biotype='protein_coding', release_id='rel')
        row = ev('clingen', 'g1', 'gene_disease_validity', 'Definitive', subject='ENSG1', level='gene', axis='clinical_validity_actionability', condition_text='hereditary cancer')
        check = dict(id='g', type='evidence', subject='gene', gene_symbol='GENE1', where={'claim_type': 'gene_disease_validity'}, regex={'claim_value': '^definitive$'})
        self.assertEqual(status(index([row], genes=[gene]), check, key=None), 'PASS')
        self.assertEqual(status(index([row], genes=[dict(gene, symbol='OTHER')]), check, key=None), 'ABSENT')


class TestHaplotypeRegionAndAssay(unittest.TestCase):
    def hap(self, name='*2'):
        return dict(haplotype_id='h1', system='star', gene_symbol='CYP2C19', allele_name=name, is_reference=False, structural=False, release_id='rel', accession='a', function='No function',
                    allele_level='core')

    def test_a_star_allele_is_checked_through_its_definition(self):
        check = dict(id='h', type='haplotype_allele', gene_symbol='CYP2C19', allele_regex=r'^(?:CYP2C19)?\*2$')
        definition = dict(haplotype_id='h1', variant_key=KEY, hgvs=None, release_id='rel', definition_scope='defining')
        self.assertEqual(status(index(haplotypes=[self.hap()], definitions=[definition]), check), 'PASS')
        self.assertEqual(status(index(haplotypes=[self.hap('CYP2C19*2')], definitions=[definition]), check), 'PASS')
        self.assertEqual(status(index(haplotypes=[self.hap()], definitions=[dict(definition, variant_key='10:101:A:G')]), check), 'FAIL')
        self.assertEqual(status(index(haplotypes=[self.hap('*3')], definitions=[definition]), check), 'ABSENT')

    def test_region_evidence_needs_overlap_kind_and_matching_rows(self):
        region = dict(region_id='r1', kind='clingen_dosage_region', chrom='22', start=100, end=500, name='x', gene_ids=None, release_id='rel')
        row = ev('clingen', 'r', 'haploinsufficiency_score', 'Sufficient evidence', subject='r1', level='region', axis='clinical_validity_actionability', claim_num=3.0)
        check = dict(id='r', type='region_evidence', chrom='22', start=200, end=300, kinds=['clingen_dosage_region'], where={'claim_type': 'haploinsufficiency_score'}, bounds={'claim_num': [3, 3]})
        self.assertEqual(status(index([row], regions=[region]), check), 'PASS')
        self.assertEqual(status(index([dict(row, claim_num=1.0)], regions=[region]), check), 'FAIL')
        self.assertEqual(status(index([row], regions=[dict(region, start=1000, end=2000)]), check), 'ABSENT')
        self.assertEqual(status(index([row], regions=[dict(region, kind='strchive_locus')]), check), 'ABSENT')

    def test_assays_are_found_by_the_gene_id_their_target_symbol_resolves_to(self):
        gene = dict(gene_id='ENSG1', symbol='BRCA1', hgnc_id='', chrom='17', start=1, end=2, strand='-', mane_select='', biotype='protein_coding', release_id='rel')
        assay = dict(assay_id='m1', source_id='mavedb', kind='dms', target_gene='ENSG1', model_system='HAP1', tissue_or_cell=None, perturbation=None, dose_context=None, readout='x',
                     score_direction=None, lineage_group='l', release_id='rel')
        check = dict(id='a', type='assay_exists', gene_symbol='BRCA1', min_rows=1)
        self.assertEqual(status(index(assays=[assay], genes=[gene]), check, key=None), 'PASS')
        self.assertEqual(status(index(assays=[dict(assay, target_gene='ENSG2')], genes=[gene]), check, key=None), 'ABSENT')
        self.assertEqual(status(index(assays=[assay], genes=[dict(gene, symbol='OTHER')]), check, key=None), 'ABSENT')           # the symbol does not resolve to the assay's gene
        self.assertEqual(status(index(assays=[dict(assay, target_gene='BRCA1')]), check, key=None), 'PASS')                      # an adapter that stores the symbol itself is found too


class TestSpec(unittest.TestCase):
    def test_a_spec_that_cannot_run_is_rejected_rather_than_passed(self):
        for bad, text in ((dict(loci=[dict(id='a', title='t', basis='b', key=KEY, checks=[dict(id='x', type='nonsense')])]), 'unknown check type'),
                          (dict(loci=[dict(id='a', title='t', basis='b', checks=[dict(id='x', type='xref_key')])]), 'needs a key'),
                          (dict(loci=[dict(id='a', title='t', basis='b', key=KEY, checks=[dict(id='x', type='xref_key')])]), 'needs an rsid'),
                          (dict(loci=[dict(id='a', title='t', basis='b', key=KEY, checks=[dict(id='x', type='orientation')])]), 'orientation needs'),
                          (dict(loci=[dict(id='a', title='t', basis='b', checks=[dict(id='x', type='evidence', subject='gene')])]), 'gene_symbol'),
                          (dict(loci=[dict(id='a', title='t', basis='', key=KEY, checks=[dict(id='x', type='evidence')])]), 'no basis'),
                          (dict(loci=[dict(id='a', title='t', basis='b', key=KEY, checks=[])]), 'no checks'),
                          (dict(loci=[dict(id='a', title='t', basis='b', key=KEY, checks=[dict(type='evidence')])]), 'without an id'),
                          (dict(loci=[dict(id='a', title='t', basis='b', key=KEY, checks=[dict(id='x', type='evidence')])] * 2), 'repeated')):
            with self.subTest(text):
                with self.assertRaises(ValueError) as cm:
                    known_loci.run(index(), bad)
                self.assertIn(text, str(cm.exception))

    def test_the_shipped_spec_is_valid_and_every_fact_states_its_basis(self):
        s = known_loci.load_spec(os.path.join(ROOT, 'known_loci.json'))
        self.assertEqual(known_loci.validate_spec(s), [])
        self.assertGreaterEqual(len(s['loci']), 15)
        kinds = {c['type'] for loc in s['loci'] for c in loc['checks']}
        self.assertEqual(kinds, set(known_loci.CHECKS))                                             # every kind of check is used by the shipped facts

    def test_the_command_line_writes_a_report_and_exits_nonzero_when_a_fact_fails(self):
        work = tempfile.mkdtemp()
        idx = os.path.join(work, 'idx')
        os.makedirs(idx)
        for name, rows in (('evidence', []), ('xref', [xref('rs1')])):
            import pyarrow.parquet as pq
            pq.write_table(schema.rows_to_table(name, rows), os.path.join(idx, f'{name}.parquet'))
        from vindex import SCHEMA_VERSION
        with open(os.path.join(idx, 'index_manifest.json'), 'w') as f:
            json.dump(dict(schema_version=SCHEMA_VERSION), f)
        good = os.path.join(work, 'good.json')
        with open(good, 'w') as f:
            json.dump(spec(dict(id='x', type='xref_key')), f)
        bad = os.path.join(work, 'bad.json')
        with open(bad, 'w') as f:
            json.dump(spec(dict(id='x', type='xref_key'), dict(id='e', type='evidence', where={'source_id': 'clinvar'})), f)
        run = lambda sp, out: subprocess.run([sys.executable, os.path.join(ROOT, 'verify_known_loci.py'), '--index', idx, '--spec', sp, '--out', out], capture_output=True, text=True)
        def report(name):
            with open(os.path.join(work, name)) as f:
                return json.load(f)
        r = run(good, os.path.join(work, 'good_report.json'))
        self.assertEqual(r.returncode, 0, r.stderr)
        self.assertTrue(report('good_report.json')['ok'])
        r = run(bad, os.path.join(work, 'bad_report.json'))
        self.assertEqual(r.returncode, 1)
        self.assertIn('ABSENT', r.stdout)
        self.assertEqual(report('bad_report.json')['totals'], dict(PASS=1, FAIL=0, ABSENT=1))


if __name__ == '__main__':
    unittest.main()
