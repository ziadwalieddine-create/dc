import gzip
import json
import os
import shutil
import tempfile
import unittest

import pysam
import pyarrow.parquet as pq

from helpers import syn_reference, vcf_text
from vindex import pgs, pgs_run
from vindex.match_small import MatchContext
from vindex.vcfio import read_vcf

ref = syn_reference()
b = lambda pos, c='1': ref.fetch(c, pos, pos)
comp = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
alt = lambda pos, c='1': next(x for x in 'ACGT' if x != b(pos, c) and x != comp[b(pos, c)])


def write_scoring(d, pid, build, rows):
    with gzip.open(os.path.join(d, f'{pid}_hmPOS_GRCh38.txt.gz'), 'wt') as f:
        f.write(f'#pgs_id={pid}\n#HmPOS_build={build}\nrsID\tchr_name\tchr_position\teffect_allele\tother_allele\teffect_weight\thm_chr\thm_pos\n')
        for i, (c, pos, e, o, w) in enumerate(rows):
            f.write(f'rs{i}\t{c}\t{pos}\t{e}\t{o}\t{w}\t{c}\t{pos}\n')


class TestChromosomeWiseScoring(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.work = tempfile.mkdtemp()
        recs = [('1', 402, b(402), alt(402), 'PASS', 50, ['0/1:60:30:15,15']), ('1', 410, b(410), alt(410), 'PASS', 50, ['1/1:60:30:0,30']),
                ('1', 420, b(420), '<NON_REF>', 'PASS', 0, ['0/0:60:30:30,0'], 'END=440'), ('2', 100, b(100, '2'), alt(100, '2'), 'PASS', 50, ['0/1:60:30:15,15']),
                ('X', 50, b(50, 'X'), alt(50, 'X'), 'PASS', 50, ['1:60:30:0,30'])]
        p = os.path.join(cls.work, 's.vcf')
        with open(p, 'w') as f:
            f.write(vcf_text(recs))
        pysam.tabix_index(p, preset='vcf', force=True)
        cls.vcf = p + '.gz'
        cls.sd = os.path.join(cls.work, 'scores')
        os.makedirs(cls.sd)
        cls.rows1 = [('1', 402, alt(402), b(402), 0.5), ('1', 410, alt(410), b(410), 0.25), ('1', 430, alt(430), b(430), 1.0), ('1', 600, alt(600), b(600), 2.0),
                     ('2', 100, alt(100, '2'), b(100, '2'), 0.125), ('X', 50, alt(50, 'X'), b(50, 'X'), 4.0), ('', 0, 'A', 'C', 8.0)]
        write_scoring(cls.sd, 'PGS000001', 'SYN38', cls.rows1)
        write_scoring(cls.sd, 'PGS000003', 'SYN38', [('1', 402, alt(402), b(402), 1.0)])
        write_scoring(cls.sd, 'PGS000004', 'hg19', [('1', 402, alt(402), b(402), 1.0)])
        cls.models = {'PGS000001': dict(name='one', reported_trait='t', gates=['uncalibrated_without_reference_distribution'], ancestry_development='European:100', ancestry_evaluation='European:100',
                                        scoring_file_sha256='abc'),
                      'PGS000003': dict(name='three', gates=[]), 'PGS000004': dict(name='four'), 'PGS000002': dict(name='two')}
        cls.out = os.path.join(cls.work, 'out')

    @classmethod
    def tearDownClass(cls):
        shutil.rmtree(cls.work, ignore_errors=True)

    def run_models(self, **kw):
        return pgs_run.score_models(self.vcf, ref, self.sd, self.models, ['PGS000001', 'PGS000002', 'PGS000003', 'PGS000004'], self.out, build='SYN38', repo_root=os.path.join(self.work, 'x'), **kw)

    def test_chromosome_wise_result_equals_a_single_pass_over_the_whole_callset(self):
        m = self.run_models(per_weight=True)
        got = {r['pgs_id']: r for r in pq.read_table(os.path.join(self.out, 'pgs_scores.parquet')).to_pylist()}
        whole = pgs.score(pgs.read_scoring_file(pgs_run.scoring_path(self.sd, 'PGS000001')), MatchContext(read_vcf(self.vcf, ref)), ref)
        r = got['PGS000001']
        self.assertEqual((r['n_weights'], r['n_scored'], r['n_missing']), (whole.n_weights, whole.n_scored, whole.n_missing))
        self.assertAlmostEqual(r['score_sum'], whole.score_sum)
        self.assertAlmostEqual(r['scored_weight_fraction'], whole.abs_weight_scored_fraction)
        self.assertEqual(json.loads(r['missing_reasons_json']), whole.reasons)
        # 402 het (0.5) + 410 hom (2 x 0.25) + 430 verified reference (0) + chr2:100 het (0.125); 600 absent, X excluded, no-position weight missing
        self.assertAlmostEqual(r['score_sum'], 0.5 + 0.5 + 0.0 + 0.125)
        self.assertEqual(json.loads(r['missing_reasons_json']), {'absent_from_callset': 1, 'non_autosomal_excluded': 1, 'no_position': 1})
        self.assertEqual(m['models'], 4)

    def test_scoring_in_small_groups_gives_exactly_the_result_of_one_group(self):
        self.run_models(per_weight=True)
        whole = {n: pq.read_table(os.path.join(self.out, n)) for n in ('pgs_scores.parquet', 'pgs_weights/PGS000001.parquet')}
        self.assertEqual(pgs_run.weight_groups({'a': {'scoring_file_variants': 5}, 'b': {'scoring_file_variants': 5}, 'c': {'scoring_file_variants': 100}}, ['a', 'b', 'c'], 6), [['a'], ['b'], ['c']])
        self.assertEqual(pgs_run.weight_groups({'a': {'scoring_file_variants': 5}, 'b': {'scoring_file_variants': 5}}, ['a', 'b'], 100), [['a', 'b']])
        self.run_models(per_weight=True, weights_budget=1)                                          # every score is its own group, so the callset is read once per score
        for name, t in whole.items():
            self.assertEqual(pq.read_table(os.path.join(self.out, name)), t, name)

    def test_gates_travel_with_the_score_and_nothing_is_calibrated_by_default(self):
        self.run_models()
        got = {r['pgs_id']: r for r in pq.read_table(os.path.join(self.out, 'pgs_scores.parquet')).to_pylist()}
        r = got['PGS000001']
        self.assertEqual((r['status'], r['calibration_status'], r['z_score'], r['ancestry_status']), ('scored', 'uncalibrated', None, 'unassessed'))
        self.assertIn('weight_coverage_below_0.9', r['flags'])
        self.assertIn('uncalibrated', r['flags'])
        self.assertEqual(r['registry_gates'], ['uncalibrated_without_reference_distribution'])

    def test_ancestry_and_calibration_come_only_from_what_is_supplied(self):
        self.run_models(sample_ancestry='African')
        r = {x['pgs_id']: x for x in pq.read_table(os.path.join(self.out, 'pgs_scores.parquet')).to_pylist()}['PGS000001']
        self.assertEqual(r['ancestry_status'], 'not_represented')
        self.assertIn('ancestry_not_represented', r['flags'])
        self.run_models(sample_ancestry='European', reference_distributions={'PGS000001': {'European': {'mean': 0.5, 'sd': 0.5}}})
        r = {x['pgs_id']: x for x in pq.read_table(os.path.join(self.out, 'pgs_scores.parquet')).to_pylist()}['PGS000001']
        self.assertEqual((r['ancestry_status'], r['calibration_status']), ('matched', 'calibrated_to_supplied_distribution'))
        self.assertAlmostEqual(r['z_score'], (1.125 - 0.5) / 0.5)
        self.assertNotIn('uncalibrated', r['flags'])

    def test_scores_that_cannot_be_computed_are_listed_with_the_reason(self):
        self.run_models()
        got = {r['pgs_id']: r for r in pq.read_table(os.path.join(self.out, 'pgs_scores.parquet')).to_pylist()}
        self.assertEqual((got['PGS000002']['status'], got['PGS000002']['reason']), ('not_scored', 'scoring_file_not_acquired'))
        self.assertEqual((got['PGS000004']['status'], got['PGS000004']['reason']), ('not_scored', 'scoring_file_harmonised_to_hg19_not_SYN38'))
        self.assertEqual(got['PGS000003']['status'], 'scored')

    def test_per_weight_audit_table(self):
        self.run_models(per_weight=True)
        w = pq.read_table(os.path.join(self.out, 'pgs_weights', 'PGS000001.parquet')).to_pylist()
        self.assertEqual(len(w), len(self.rows1))
        by = {(x['chrom'], x['pos']): x for x in w}
        self.assertEqual((by[('1', 402)]['dosage'], by[('1', 410)]['dosage'], by[('1', 430)]['dosage']), (1, 2, 0))
        self.assertIsNone(by[('1', 600)]['dosage'])
        self.assertEqual(by[('1', 600)]['basis_or_reason'], 'absent_from_callset')

    def test_results_stay_out_of_the_tracked_tree(self):
        repo = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..'))
        with self.assertRaises(ValueError):
            pgs_run.score_models(self.vcf, ref, self.sd, self.models, ['PGS000001'], os.path.join(repo, 'variant_index', 'out'), build='SYN38', repo_root=repo)


if __name__ == '__main__':
    unittest.main()
