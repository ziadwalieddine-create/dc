import json
import os
import sys
import tempfile
import unittest

from helpers import syn_reference, FIXTURES
from vindex.build import runner
from vindex.build.gnomad import GnomadConstraintAdapter, GnomadLofCurationAdapter, GnomadFrequencyAdapter
from vindex.index_reader import Index

sys.path.insert(0, FIXTURES)
import make_sources  # noqa: E402

SOURCES = make_sources.write_all(tempfile.mkdtemp())
GN = os.path.join(SOURCES, 'gnomad')
b, alt = make_sources.b, make_sources.alt


class FakeGenes:
    by_ensg = {'ENSG_G1': 1, 'ENSG_G2': 1, 'ENSG_G4': 1}

    def from_ncbi(self, ncbi):
        return {'1001': 'ENSG_G1', '1002': 'ENSG_G2', '1004': 'ENSG_G4'}.get(str(ncbi))


def adapters():
    return [GnomadConstraintAdapter(FakeGenes(), os.path.join(GN, 'gnomad.v4.1.2.constraint_metrics.tsv.bgz')),
            GnomadLofCurationAdapter(FakeGenes(), os.path.join(GN, 'lof_curation', 'incomplete_penetrance_curation_results.csv')),
            GnomadFrequencyAdapter(os.path.join(GN, 'regional', 'gnomad_v4_1_joint_at_mavedb_positions.tsv.gz'), os.path.join(GN, 'regional', 'gnomad_v4_1_exomes_at_mavedb_positions.tsv.gz'))]


class TestGnomad(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.ref = syn_reference()
        cls.out = tempfile.mkdtemp()
        cls.man = runner.build(cls.out, adapters(), cls.ref)
        cls.idx = Index.open(cls.out)
        cls.ev = [r for r in cls.idx.table('evidence').to_pylist() if r['source_id'] in ('gnomad', 'gnomad_constraint')]
        cls.stats = {a['adapter']: a['stats'] for a in cls.man['adapters']}
        cls.k402 = f"1:402:{b('1', 402)}:{alt('1', 402)}"
        cls.k410 = f"1:410:{b('1', 410)}:{alt('1', 410)}"

    def rows(self, **kw):
        out = self.ev
        for k, v in kw.items():
            out = [r for r in out if r[k] == v]
        return out

    def test_constraint_uses_mane_select_only_and_stays_at_gene_level(self):
        s = self.stats['GnomadConstraintAdapter']
        self.assertEqual((s['transcript_rows'], s['mane_rows']), (8, 7))
        lof = self.rows(claim_type='gnomad_lof_constraint', subject_id='ENSG_G1')[0]
        self.assertEqual((lof['level'], lof['axis'], lof['evidence_kind'], lof['effect_assertion'], lof['claim_value'], lof['claim_num'], lof['score']), ('gene', 'mechanistic_confidence', 'population_observation', 'not_applicable', 'LOEUF 0.210', 0.21, 0.098))
        x = json.loads(lof['extra_json'])
        self.assertEqual((x['observed'], x['expected'], x['pLI'], x['upper_decile'], x['transcript']), (2.0, 20.5, 1.0, 0, 'ENST1'))
        mis = self.rows(claim_type='gnomad_missense_constraint', subject_id='ENSG_G1')[0]
        self.assertEqual((mis['claim_num'], json.loads(mis['extra_json'])['z_score']), (0.52, 4.2))
        self.assertEqual(self.rows(claim_type='gnomad_lof_constraint', subject_id='ENSG_G2'), [])          # no loss-of-function metric: no row
        self.assertEqual(json.loads(self.rows(claim_type='gnomad_missense_constraint', subject_id='ENSG_G2')[0]['extra_json'])['flags'], ['no_exp_lof'])
        self.assertTrue(all(r['level'] == 'gene' for r in self.rows(source_id='gnomad_constraint')))
        un = [r for r in self.rows(source_id='gnomad_constraint') if r['mapping_status'] == 'unmapped']
        self.assertEqual(sorted((r['unmapped_key'], r['mapping_reason']) for r in un), [('ENSG_UNKNOWN', 'gene_not_in_gene_map'), ('NCBIGene:1005', 'gene_not_in_gene_map')])

    def test_each_mane_select_gene_is_one_measurement_whichever_transcript_source_lists_it(self):
        s = self.stats['GnomadConstraintAdapter']
        self.assertEqual((s['refseq_mane_rows_repeating_an_ensembl_row'], s['refseq_mane_rows_without_an_ensembl_partner']), (2, 2))
        self.assertEqual(len(self.rows(claim_type='gnomad_lof_constraint', subject_id='ENSG_G1')), 1)               # not once per transcript source
        g4 = self.rows(claim_type='gnomad_lof_constraint', subject_id='ENSG_G4')
        self.assertEqual((len(g4), g4[0]['mapping_status'], g4[0]['gene_ids']), (1, 'mapped', ['ENSG_G4']))           # reached through its NCBI gene id

    def test_lof_curation_rows(self):
        a = self.rows(claim_type='gnomad_lof_curation_verdict', subject_id=self.k402)[0]
        self.assertEqual((a['claim_value'], a['effect_assertion'], a['axis'], a['evidence_kind'], a['gene_ids']), ('LoF', 'effect', 'mechanistic_confidence', 'curated_expert_classification', ['ENSG_G1']))
        n = self.rows(claim_type='gnomad_lof_curation_verdict', subject_id=self.k410)[0]
        self.assertEqual((n['claim_value'], n['effect_assertion']), ('Not LoF', 'no_effect'))
        self.assertEqual(json.loads(n['extra_json'])['flags'], sorted(['Last Exon', 'LoF in untranslated transcript']))
        s = self.stats['GnomadLofCurationAdapter']
        self.assertEqual((s['records'], s['curated_variants:mapped'], s['curated_variants:rejected'], s['curated_variants:unmapped']), (4, 2, 1, 1))
        rej = [r for r in self.rows(claim_type='gnomad_lof_curation_verdict') if r['mapping_status'] == 'rejected'][0]
        self.assertEqual((rej['subject_id'], rej['mapping_reason']), (None, 'ref_mismatch'))

    def test_frequency_rows_keep_filters_and_ancestry_detail(self):
        f = self.rows(claim_type='gnomad_allele_frequency', subject_id=self.k402)[0]
        self.assertEqual((f['axis'], f['evidence_kind'], f['eaf'], f['ancestry'], f['sample_size'], f['review_status']), ('applicability_context', 'population_observation', 2e-05, 'afr', 100000, 'PASS'))
        x = json.loads(f['extra_json'])
        self.assertEqual((x['ac'], x['an'], x['by_ancestry']['afr']['ac']), (4, 200000, 4))
        filt = self.rows(claim_type='gnomad_allele_frequency', subject_id=self.k410)[0]
        self.assertEqual((filt['review_status'], filt['eaf'], json.loads(filt['extra_json'])['exomes_filters']), ('EXOMES_FILTERED', 0.0, 'AC0'))        # a filtered site is not 'absent'
        self.assertEqual(self.stats['GnomadFrequencyAdapter']['frequency_rows:rejected'], 1)

    def test_predictions_are_marked_as_predictions(self):
        p = self.rows(claim_type='prediction:CADD', subject_id=self.k402)[0]
        self.assertEqual((p['evidence_kind'], p['score'], p['axis'], p['effect_assertion'], p['score_unit'], p['lineage_group']), ('computational_prediction', 26.5, 'mechanistic_confidence', 'not_applicable', 'CADD phred', 'predictor:CADD'))
        names = sorted(r['claim_type'] for r in self.rows(subject_id=self.k402, evidence_kind='computational_prediction'))
        self.assertEqual(names, sorted(['prediction:CADD', 'prediction:PolyPhen-2', 'prediction:REVEL', 'prediction:SpliceAI', 'prediction:Pangolin']))
        self.assertEqual(self.rows(claim_type='prediction:SIFT'), [])                          # empty cell: no row
        self.assertEqual(self.rows(subject_id=self.k410, evidence_kind='computational_prediction'), [])

    def test_ids_unique_and_rebuild_idempotent(self):
        ids = [r['evidence_id'] for r in self.ev]
        self.assertEqual(len(ids), len(set(ids)))
        out2 = tempfile.mkdtemp()
        runner.build(out2, adapters(), self.ref)
        with open(os.path.join(self.out, 'table_digests.json')) as f1, open(os.path.join(out2, 'table_digests.json')) as f2:
            self.assertEqual(json.load(f1), json.load(f2))


if __name__ == '__main__':
    unittest.main()
