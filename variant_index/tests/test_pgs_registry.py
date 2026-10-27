import json
import os
import shutil
import sys
import tempfile
import unittest

import pyarrow.parquet as pq

from helpers import FIXTURES
from vindex import pgs, pgs_registry as reg

sys.path.insert(0, FIXTURES)
import make_sources  # noqa: E402

SOURCES = make_sources.write_all(tempfile.mkdtemp())
META = os.path.join(SOURCES, 'pgs_meta')


class TestParsing(unittest.TestCase):
    def test_metric_strings(self):
        self.assertEqual(reg.metric('1.55 [1.52,1.58]'), (None, 1.55, 1.52, 1.58))
        self.assertEqual(reg.metric('R2: 0.07 [0.05, 0.09]'), ('R2', 0.07, 0.05, 0.09))
        self.assertEqual(reg.metric('Partial-R2 = 0.07'), ('Partial-R2', 0.07, None, None))
        self.assertEqual(reg.metric('1e-3'), (None, 0.001, None, None))
        self.assertIsNone(reg.metric('AUC about seventy percent'))
        self.assertIsNone(reg.metric(''))

    def test_ancestry_distribution(self):
        self.assertEqual(reg.ancestry('European:72.7|Not Reported:18.2|East Asian:9.1'), {'European': 72.7, 'Not Reported': 18.2, 'East Asian': 9.1})
        self.assertEqual(reg.ancestry(''), {})
        self.assertEqual(reg.ancestry('garbage|European:50'), {'European': 50.0})


class TestTraitSelection(unittest.TestCase):
    models = {'PGS1': dict(efo_ids=['MONDO_0004989', 'EFO_1'], scoring_file='acquired'), 'PGS2': dict(efo_ids=['EFO_0004339'], scoring_file='acquired'),
              'PGS3': dict(efo_ids=['EFO_0004339'], scoring_file='listed_not_acquired'), 'PGS4': dict(efo_ids=None, scoring_file='acquired')}

    def test_a_trait_is_found_in_either_spelling_and_only_among_acquired_scores(self):
        # the PGS Catalog writes EFO_0004339, the index and its views EFO:0004339
        for trait in ('EFO_0004339', 'EFO:0004339', ' EFO:0004339 '):
            self.assertEqual(reg.models_for_trait(self.models, trait), ['PGS2'], trait)
        self.assertEqual(reg.models_for_trait(self.models, 'MONDO:0004989'), ['PGS1'])
        self.assertEqual(reg.models_for_trait(self.models, 'EFO:9'), [])
        self.assertEqual(reg.models_for_trait(self.models, ''), [])


class TestRegistry(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.out = tempfile.mkdtemp()
        cls.man = reg.build_registry(META, os.path.join(META, 'scoring_manifest.tsv'), cls.out)
        cls.m = {r['pgs_id']: r for r in pq.read_table(os.path.join(cls.out, 'pgs_model.parquet')).to_pylist()}

    @classmethod
    def tearDownClass(cls):
        shutil.rmtree(cls.out, ignore_errors=True)

    def test_every_catalogue_score_is_listed_with_the_state_of_its_scoring_file(self):
        self.assertEqual(sorted(self.m), ['PGS000001', 'PGS000002', 'PGS000003', 'PGS000004'])
        self.assertEqual({k: v['scoring_file'] for k, v in self.m.items()}, {'PGS000001': 'acquired', 'PGS000002': 'listed_not_acquired', 'PGS000003': 'acquired', 'PGS000004': 'listed_not_acquired'})
        self.assertEqual((self.man['models'], self.man['scoring_files_acquired']), (4, 2))
        self.assertEqual(self.man['acquired_not_in_catalogue_metadata'], ['PGS999999'])             # a stored file the catalogue metadata does not list is reported, not ignored
        self.assertEqual((self.m['PGS000001']['scoring_file_tier'], self.m['PGS000001']['scoring_file_variants']), ('A', 10))

    def test_gates_state_why_a_score_needs_a_further_step(self):
        g = lambda i: set(self.m[i]['gates'])
        self.assertEqual(g('PGS000001'), {'uncalibrated_without_reference_distribution', 'weight_type_not_reported'})
        self.assertIn('interaction_terms_not_additive', g('PGS000002'))
        self.assertIn('no_grch38_scoring_file_acquired', g('PGS000002'))
        self.assertIn('no_ancestry_information', g('PGS000003'))
        self.assertIn('scoring_file_variant_count_differs_from_header', g('PGS000003'))
        self.assertIn('weights_not_on_log_scale:OR', g('PGS000003'))
        self.assertEqual(g('PGS000004'), {'no_grch38_scoring_file_acquired', 'uncalibrated_without_reference_distribution'})         # an unfamiliar weight label is kept, not guessed to be a problem

    def test_ancestry_distributions_keep_text_and_parsed_forms_and_feed_the_scorer_gate(self):
        r = self.m['PGS000001']
        self.assertEqual(json.loads(r['ancestry_evaluation_json']), {'European': 72.7, 'Not Reported': 18.2, 'East Asian': 9.1})
        meta = reg.ancestry_meta(r)
        self.assertEqual(meta['development'], 'European:100')
        self.assertEqual(pgs._ancestry_assess(meta['development'], meta['evaluation'], 'African')['status'], 'not_represented')
        self.assertEqual(pgs._ancestry_assess(meta['development'], meta['evaluation'], 'European')['status'], 'matched')
        self.assertEqual(pgs._ancestry_assess(None, None, None)['status'], 'unassessed')
        self.assertEqual(r['efo_ids'], ['MONDO_0004989', 'EFO_1'])
        self.assertEqual(r['licence'], 'cite appropriately')

    def test_performance_metrics_keep_value_interval_and_original_text(self):
        perf = pq.read_table(os.path.join(self.out, 'pgs_performance.parquet')).to_pylist()
        by = {(r['ppm_id'], r['metric_family'], r['metric_name']): r for r in perf}
        orr = by[('PPM000001', 'OR', 'OR')]
        self.assertEqual((orr['value'], orr['ci_low'], orr['ci_high'], orr['text'], orr['covariates']), (1.55, 1.52, 1.58, '1.55 [1.52,1.58]', 'age, PCs'))
        self.assertEqual(by[('PPM000001', 'other', 'R2')]['value'], 0.07)
        weird = [r for r in perf if r['text'] == 'unparsable text'][0]
        self.assertEqual((weird['value'], weird['metric_family']), (None, 'other'))                  # kept with its text, no number invented
        self.assertEqual(by[('PPM000002', 'AUROC', 'AUROC')]['value'], 0.7)
        self.assertEqual(len(perf), 5)

    def test_samples_cover_development_and_evaluation(self):
        smp = pq.read_table(os.path.join(self.out, 'pgs_sample.parquet')).to_pylist()
        self.assertEqual(sorted((r['stage'], r['pgs_id']) for r in smp), [('Score development', 'PGS000001'), ('Variant associations (GWAS)', 'PGS000001'), ('evaluation', 'PGS000001'), ('evaluation', 'PGS000003')])
        ev = [r for r in smp if r['stage'] == 'evaluation' and r['pgs_id'] == 'PGS000001'][0]
        self.assertEqual((ev['n_individuals'], ev['n_cases'], ev['ancestry_broad'], ev['cohorts']), (67054, 33673, 'European', 'ABC|DEF'))

    def test_registry_is_deterministic(self):
        out2 = tempfile.mkdtemp()
        try:
            m2 = reg.build_registry(META, os.path.join(META, 'scoring_manifest.tsv'), out2)
            self.assertEqual(self.man['table_sha256'], m2['table_sha256'])
        finally:
            shutil.rmtree(out2, ignore_errors=True)


if __name__ == '__main__':
    unittest.main()
