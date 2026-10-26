"""The cross-check recomputes calls and evidence matches with bcftools only. It must agree with a clean run and must be able to fail."""
import os
import shutil
import tempfile
import unittest

import pyarrow as pa
import pyarrow.compute as pc
import pyarrow.parquet as pq
import pysam

from helpers import find_tool, syn_reference
from test_wgs import b, alt, build_index, combined_clinvar_vcf, sample_vcf
from vindex import crosscheck, wgs


@unittest.skipUnless(find_tool('bcftools'), 'bcftools not available (set VINDEX_TOOLS_BIN)')
class TestCrosscheck(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.work = tempfile.mkdtemp()
        cls.idx = os.path.join(cls.work, 'idx')
        build_index(cls.idx)
        other = lambda c, p: [x for x in 'ACGT' if x not in (b(c, p), alt(c, p))][0]
        extra = [('1', 160, b('1', 160), f'{alt("1", 160)},{other("1", 160)}', 'PASS', 50, ['0/2:60:30:15,0,15']),            # the sample carries only the second alt
                 ('1', 170, b('1', 170), f'{alt("1", 170)},{other("1", 170)}', 'PASS', 50, ['1/2:60:30:0,15,15'])]           # carries both
        cls.vcf = sample_vcf(cls.work, extra=extra)
        cls.res = os.path.join(cls.work, 'res')
        cls.ref = syn_reference()
        wgs.run_sample(cls.vcf, cls.idx, cls.res, cls.ref, repo_root=os.path.join(cls.work, 'not_the_repo'))
        cv = os.path.join(cls.work, 'cv.vcf')
        combined_clinvar_vcf(cv)
        with open(cv) as f:
            lines = f.read().rstrip('\n').split('\n')
        head, body = [l for l in lines if l.startswith('#')], [l.split('\t') for l in lines if not l.startswith('#')]
        body.sort(key=lambda r: (r[0], int(r[1])))                         # tabix needs coordinate order
        tags = sorted({kv.split('=')[0] for r in body for kv in r[7].split(';')})
        head = head[:-1] + [f'##INFO=<ID={t},Number=.,Type=String,Description="x">' for t in tags] + head[-1:]            # a real ClinVar header declares its INFO tags
        with open(cv, 'w') as f:
            f.write('\n'.join(head + ['\t'.join(r) for r in body]) + '\n')
        pysam.tabix_index(cv, preset='vcf', force=True)
        cls.cv = cv + '.gz'

    @classmethod
    def tearDownClass(cls):
        shutil.rmtree(cls.work, ignore_errors=True)

    def check(self, res=None):
        return crosscheck.crosscheck(self.vcf, res or self.res, self.ref, self.cv)

    def copy_results(self, name):
        d = os.path.join(self.work, name)
        shutil.copytree(self.res, d)
        return d

    def test_a_clean_run_agrees_and_a_refused_ref_mismatch_is_explained(self):
        r = self.check()
        self.assertTrue(r['ok'], r['examples'])
        c = r['calls']
        self.assertEqual(c['bcftools'], c['matcher'] + 1)
        self.assertEqual((c['only_bcftools'], c['only_bcftools_explained_by_rejected_calls'], c['only_bcftools_unexplained'], c['only_matcher']), (1, 1, 0, 0))        # 1:440, REF disagrees with the reference
        e = r['evidence']
        self.assertGreater(e['matcher'], 0)
        self.assertEqual((e['only_bcftools_unexplained'], e['only_matcher']), (0, 0))
        self.assertEqual(e['bcftools'], e['matcher'] + e['only_bcftools_explained_by_rejected_calls'])

    def test_multiallelic_records_contribute_only_the_alleles_the_sample_carries(self):
        calls = {r['key'] for r in pq.read_table(os.path.join(self.res, 'small_calls.parquet')).to_pylist()}
        a, o = alt('1', 160), [x for x in 'ACGT' if x not in (b('1', 160), alt('1', 160))][0]
        self.assertIn(f'1:160:{b("1", 160)}:{o}', calls)
        self.assertNotIn(f'1:160:{b("1", 160)}:{a}', calls)                                  # 0/2: the first alt is not carried
        self.assertEqual({k for k in calls if k.startswith('1:170:')}, {f'1:170:{b("1", 170)}:{alt("1", 170)}', f'1:170:{b("1", 170)}:{[x for x in "ACGT" if x not in (b("1", 170), alt("1", 170))][0]}'})
        self.assertTrue(self.check()['ok'])

    def test_a_dropped_call_is_reported_as_unexplained(self):
        d = self.copy_results('dropped')
        t = pq.read_table(os.path.join(d, 'small_calls.parquet'))
        victim = f'1:402:{b("1", 402)}:{alt("1", 402)}'
        pq.write_table(t.filter(pc.not_equal(t['key'], victim)), os.path.join(d, 'small_calls.parquet'))
        r = self.check(d)
        self.assertFalse(r['ok'])
        self.assertEqual(r['calls']['only_bcftools_unexplained'], 1)
        self.assertEqual(r['examples']['only_bcftools_unexplained'], [victim])

    def test_an_invented_call_and_a_missing_evidence_link_are_reported(self):
        d = self.copy_results('invented')
        t = pq.read_table(os.path.join(d, 'small_calls.parquet'))
        row = t.slice(0, 1).to_pylist()[0]
        row['key'] = f'1:1999:{b("1", 1999)}:{alt("1", 1999)}'
        pq.write_table(pa.concat_tables([t, pa.Table.from_pylist([row], schema=t.schema)]), os.path.join(d, 'small_calls.parquet'))
        links = pq.read_table(os.path.join(d, 'evidence_links.parquet'))
        gone = f'1:402:{b("1", 402)}:{alt("1", 402)}'
        pq.write_table(links.filter(pc.not_equal(links['call_id'], gone)), os.path.join(d, 'evidence_links.parquet'))
        r = self.check(d)
        self.assertFalse(r['ok'])
        self.assertEqual(r['calls']['only_matcher'], 1)
        self.assertEqual(r['evidence']['only_bcftools_unexplained'], 1)
        self.assertEqual(r['examples']['evidence_only_bcftools'], [gone])


if __name__ == '__main__':
    unittest.main()
