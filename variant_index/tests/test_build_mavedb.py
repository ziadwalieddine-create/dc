import json
import os
import sys
import tempfile
import unittest

from helpers import syn_reference, FIXTURES, SYN_FASTA, SYN_BUILDS
from vindex.refgenome import Reference
from vindex.build import runner
from vindex.build.mavedb import MaveDbAdapter, in_range
from vindex.index_reader import Index

sys.path.insert(0, FIXTURES)
import make_sources  # noqa: E402

SOURCES = make_sources.write_all(tempfile.mkdtemp())
KEYS = make_sources.write_all.mavedb
MV = os.path.join(SOURCES, 'mavedb')


class FakeGenes:
    def resolve(self, hgnc=None, symbol=None, ncbi=None):
        return ('ENSG_FAKE_' + symbol, 'symbol') if symbol else (None, 'unresolved')


def adapters():
    return [MaveDbAdapter(FakeGenes(), vep_path=os.path.join(MV, 'MaveDB_variants.tsv.gz'), main_json=os.path.join(MV, 'main.json'), ann_dir=os.path.join(MV, 'csv'),
                          scores_dir=os.path.join(MV, 'scores_genomic_mapped'))]


def reference():
    return Reference(SYN_FASTA, assembly_report=os.path.join(SOURCES, 'assembly_report.txt'), name='SYN38', builds=SYN_BUILDS)


class TestRanges(unittest.TestCase):
    def test_bounds_are_inclusive_only_when_stated(self):
        self.assertTrue(in_range(1.0, None, 1.0, False, True))
        self.assertFalse(in_range(1.0, None, 1.0, False, False))
        self.assertTrue(in_range(1.5, 1.5, None, True, False))
        self.assertFalse(in_range(1.5, 1.5, None, False, False))
        self.assertTrue(in_range(-1e9, None, None, False, False))


class TestMaveDb(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.ref = reference()
        cls.out = tempfile.mkdtemp()
        cls.man = runner.build(cls.out, adapters(), cls.ref)
        cls.idx = Index.open(cls.out)
        cls.ev = [r for r in cls.idx.table('evidence').to_pylist() if r['source_id'] == 'mavedb']
        cls.st = cls.man['adapters'][0]['stats']

    def rows(self, **kw):
        out = self.ev
        for k, v in kw.items():
            out = [r for r in out if r[k] == v]
        return out

    def test_row_accounting(self):
        s = self.st
        self.assertEqual(s['rows'], 11)
        self.assertEqual((s['rows_multi_variant_construct'], s['rows_superseded_score_set'], s['rows_duplicate_allele_for_accession'], s['rows_score_not_numeric'],
                          s['rows_score_set_not_in_dump_metadata']), (1, 1, 1, 1, 1))
        # integrated: 402 (kept once), 410, 420, the rejected 440, and two protein-level candidate alleles
        self.assertEqual(sum(v for k, v in s.items() if k.startswith('score_rows:') and not k.startswith('score_rows:native')), 6)
        self.assertEqual(s['rows'] - 5, 6)
        self.assertEqual(sum(v for k, v in s.items() if k.startswith('score_rows:native')), 3)        # the single-variant genomic rows of score set 5, listed in test_native_genomic_mapping

    def test_allele_level_score_row(self):
        r = self.rows(claim_type='mavedb_variant_score', subject_id=KEYS['k402'])[0]
        self.assertEqual((r['level'], r['axis'], r['evidence_kind'], r['mapping_status'], r['mapping_reason']), ('variant', 'empirical_functional_support', 'direct_experimental', 'mapped', None))
        self.assertEqual((r['score'], r['se'], r['score_unit'], r['effect_assertion']), (0.2, 0.05, 'mavedb_score:urn:mavedb:00000001-a-1', 'not_applicable'))
        self.assertEqual((r['lineage_group'], r['model_system'], r['gene_ids'], r['accession']), ('pmid:12345', 'Immortalized human cells', ['ENSG00000000001'], 'urn:mavedb:00000001-a-1#1'))
        extra = json.loads(r['extra_json'])
        self.assertEqual((extra['sd'], extra['assay_hgvs_nt'], extra['other_score_columns']), (0.1, 'c.1' + KEYS['k402'].split(':')[2] + '>' + KEYS['k402'].split(':')[3], {'score_rep1': '0.3'}))
        self.assertEqual(len(self.rows(claim_type='mavedb_variant_score', subject_id=KEYS['k402'], accession='urn:mavedb:00000001-a-1#1')), 1)         # the repeated row did not become a second row

    def test_calibrations_are_separate_interpretation_rows(self):
        cal = self.rows(claim_type='mavedb_calibrated_class', subject_id=KEYS['k402'])
        by = {r['curator']: r for r in cal}
        self.assertEqual(sorted(by), ['Research calibration', 'Test calibration'])
        t = by['Test calibration']
        self.assertEqual((t['claim_value'], t['effect_assertion'], t['axis'], t['evidence_kind']), ('abnormal|PS3|MODERATE', 'effect', 'empirical_functional_support', 'curated_expert_classification'))
        self.assertEqual((t['lineage_group'], json.loads(t['extra_json'])['points']), ('pmid:12345', 2))
        self.assertIn('not an additional study', t['lineage_note'])
        r = by['Research calibration']
        self.assertEqual((r['claim_value'], r['review_status']), ('abnormal||', 'research use only'))
        hi = self.rows(claim_type='mavedb_calibrated_class', subject_id=KEYS['k410'])
        self.assertEqual([(x['claim_value'], x['effect_assertion']) for x in hi], [('normal|BS3|SUPPORTING', 'no_effect')])
        edge = self.rows(claim_type='mavedb_calibrated_class', accession='urn:mavedb:00000001-a-1#3')              # score 1.0: inclusive upper bound of the second calibration only
        self.assertEqual([x['curator'] for x in edge], ['Research calibration'])

    def test_ref_mismatch_is_rejected_not_repaired(self):
        r = self.rows(accession='urn:mavedb:00000001-a-1#4', claim_type='mavedb_variant_score')[0]
        self.assertEqual((r['mapping_status'], r['mapping_reason'], r['subject_id']), ('rejected', 'ref_mismatch', None))
        self.assertTrue(r['unmapped_key'].startswith('1:440:'))
        self.assertTrue(all(x['mapping_status'] == 'rejected' for x in self.rows(accession='urn:mavedb:00000001-a-1#4')))

    def test_protein_level_assay_rows_say_what_they_measure(self):
        rows = self.rows(accession='urn:mavedb:00000003-a-1#7')
        self.assertEqual(sorted(r['subject_id'] for r in rows), sorted([KEYS['kmnv1'], KEYS['kmnv2']]))
        self.assertTrue(all(r['claim_type'] == 'mavedb_protein_level_score' and r['mapping_status'] == 'mapped' and r['mapping_reason'] == 'protein_level_assay_allele_inferred_from_genetic_code' for r in rows))
        extra = json.loads(rows[0]['extra_json'])
        self.assertEqual((extra['protein_hgvs'], extra['assay_hgvs_pro']), ('NP_000001.1:p.Gly170Ala', 'p.Gly170Ala'))
        self.assertIn('amino-acid change measured', extra['scope'])

    def test_superseded_and_multi_variant_sets_contribute_no_rows(self):
        self.assertEqual(self.rows(accession='urn:mavedb:00000002-a-1#1'), [])
        self.assertEqual(self.rows(accession='urn:mavedb:00000004-a-1#9'), [])

    def test_assay_entities_cover_every_score_set(self):
        a = {r['assay_id']: r for r in self.idx.table('entity_assay').to_pylist()}
        self.assertEqual(sorted(a), [f'urn:mavedb:0000000{i}-a-1' for i in (1, 2, 3, 4, 5)])
        self.assertEqual((a['urn:mavedb:00000001-a-1']['kind'], a['urn:mavedb:00000001-a-1']['target_gene'], a['urn:mavedb:00000001-a-1']['model_system']),
                         ('endogenous_locus_variant_library', 'ENSG00000000001', 'Immortalized human cells'))
        self.assertEqual(self.st['assays_superseded'], 1)
        self.assertEqual(a['urn:mavedb:00000001-a-1']['readout'], 'Loss of function')

    def test_native_genomic_mapping(self):
        n = [r for r in self.ev if r['accession'].startswith('urn:mavedb:00000005-a-1#') and r['claim_type'] == 'mavedb_variant_score']
        by = {r['accession'].split('#')[1]: r for r in n}
        self.assertEqual(sorted(by), ['1', '2', '6'])
        self.assertEqual((by['1']['subject_id'], by['1']['score'], by['1']['se'], by['1']['mapping_status']), (KEYS['k402'], 0.3, 0.02, 'mapped'))
        ex = json.loads(by['1']['extra_json'])
        self.assertEqual((ex['mapping_basis'], ex['vep_consequence'], ex['score_locator'], ex['post_mapped_hgvs_g'].split(':')[0]), ('mavedb_post_mapped_hgvs_g', 'missense_variant', 'urn-mavedb-00000005-a-1.scores.csv.gz#L2', 'NC_SYN001.1'))
        self.assertEqual((by['2']['mapping_status'], by['2']['subject_id']), ('rejected', None))                                                        # wrong REF in the HGVS: kept as rejected
        self.assertTrue(by['2']['mapping_reason'].startswith('reference_base_differs'))
        self.assertEqual(by['6']['subject_id'], '1:705:GCC:G')                                                                                         # deletion normalised
        s = self.st
        self.assertEqual((s['native_rows'], s['native_rows_multi_variant'], s['native_rows_without_numeric_score'], s['native_rows_reference_sequence']), (7, 1, 2, 1))                  # 2: one genomic row has no score row, one has a non-numeric score
        cal = [r for r in self.ev if r['accession'] == 'urn:mavedb:00000005-a-1#1' and r['claim_type'] == 'mavedb_calibrated_class']
        self.assertEqual(sorted(r['claim_value'] for r in cal), ['abnormal|PS3|MODERATE', 'abnormal||'])
        # exact record text is what the digest covers: re-read the annotation line and compare
        import hashlib
        with open(os.path.join(MV, 'csv', 'urn-mavedb-00000005-a-1.annotations.csv'), newline='') as f:
            lines = f.read().split('\r\n')                 # the fixture is written with the csv module's default line terminator
        ln = int(by['1']['locator'].rsplit('#L', 1)[1])
        self.assertEqual(hashlib.sha256(lines[ln - 1].encode()).hexdigest(), by['1']['raw_sha256'])

    def test_variant_gene_links_are_not_repeated(self):
        vg = [r for r in self.idx.table('variant_gene').to_pylist() if r['relation'] == 'mavedb_target_gene' and r['variant_key'] == KEYS['k402']]
        self.assertEqual(len(vg), 1)

    def test_evidence_ids_are_unique_and_rebuild_is_idempotent(self):
        ids = [r['evidence_id'] for r in self.ev]
        self.assertEqual(len(ids), len(set(ids)))
        out2 = tempfile.mkdtemp()
        runner.build(out2, adapters(), self.ref)
        with open(os.path.join(self.out, 'table_digests.json')) as f1, open(os.path.join(out2, 'table_digests.json')) as f2:
            self.assertEqual(json.load(f1), json.load(f2))

    def test_unexpected_header_is_an_error_not_silent(self):
        import gzip
        bad = os.path.join(tempfile.mkdtemp(), 'bad.tsv.gz')
        with gzip.open(bad, 'wt') as f:
            f.write('chr\tstart\n1\t2\n')
        with self.assertRaises(ValueError):
            runner.build(tempfile.mkdtemp(), [MaveDbAdapter(FakeGenes(), vep_path=bad, main_json=os.path.join(MV, 'main.json'))], self.ref)


if __name__ == '__main__':
    unittest.main()
