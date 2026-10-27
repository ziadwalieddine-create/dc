"""compare_runs: two runs of one input are identical whatever the moment; a changed table, a changed manifest field or a missing table is found."""
import contextlib
import io
import json
import os
import shutil
import tempfile
import unittest

import pyarrow as pa
import pyarrow.parquet as pq

from helpers import ROOT, syn_reference
from test_wgs import build_index, sample_vcf
from vindex import wgs

import sys
sys.path.insert(0, ROOT)
import compare_runs  # noqa: E402


class TestCompareRuns(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.work = tempfile.mkdtemp()
        idx = os.path.join(cls.work, 'idx')
        build_index(idx)
        vcf = sample_vcf(cls.work)
        cls.a, cls.b = os.path.join(cls.work, 'a'), os.path.join(cls.work, 'b')
        for out in (cls.a, cls.b):
            wgs.run_sample(vcf, idx, out, syn_reference(), repo_root=os.path.join(cls.work, 'not_the_repo'))

    @classmethod
    def tearDownClass(cls):
        shutil.rmtree(cls.work, ignore_errors=True)

    def copy(self, name):
        d = os.path.join(self.work, name)
        shutil.rmtree(d, ignore_errors=True)
        shutil.copytree(self.b, d)
        return d

    def test_two_runs_of_one_input_are_identical_although_their_clocks_differ(self):
        for run in (self.a, self.b):
            with open(os.path.join(run, 'run_manifest.json')) as f:
                self.assertTrue({'finished_utc', 'seconds'} <= set(json.load(f)))                                  # the manifests carry a moment and a duration ...
        res = compare_runs.compare(self.a, self.b)
        self.assertTrue(res['identical'], res)                                                                   # ... which are the only fields allowed to differ
        self.assertGreaterEqual(res['tables'], 8)

    def test_a_changed_table_a_changed_manifest_field_and_a_missing_table_are_found(self):
        d = self.copy('changed_table')
        t = pq.read_table(os.path.join(d, 'sv_calls.parquet'))
        i = t.schema.get_field_index('n_equivalent')
        vals = t['n_equivalent'].to_pylist()
        pq.write_table(t.set_column(i, 'n_equivalent', pa.array([(v or 0) + 1 for v in vals], t.schema.field('n_equivalent').type)), os.path.join(d, 'sv_calls.parquet'))
        res = compare_runs.compare(self.a, d)
        self.assertEqual((res['identical'], res['tables_differ']), (False, ['sv_calls.parquet']))
        d = self.copy('changed_manifest')
        with open(os.path.join(d, 'run_manifest.json')) as f:
            m = json.load(f)
        m['counts']['calls'] += 1
        with open(os.path.join(d, 'run_manifest.json'), 'w') as f:
            json.dump(m, f)
        res = compare_runs.compare(self.a, d)
        self.assertEqual((res['identical'], res['manifest_fields_differ']), (False, ['counts']))
        d = self.copy('missing_table')
        os.remove(os.path.join(d, 'repeat_calls.parquet'))
        res = compare_runs.compare(self.a, d)
        self.assertEqual((res['identical'], res['only_in_a']), (False, ['repeat_calls.parquet']))

    def test_a_directory_without_tables_is_not_identical_to_itself(self):
        d = os.path.join(self.work, 'empty')
        os.makedirs(d, exist_ok=True)
        with open(os.path.join(d, 'run_manifest.json'), 'w') as f:
            json.dump({}, f)
        self.assertFalse(compare_runs.compare(d, d)['identical'])                                                 # nothing was compared: that is not a pass

    def test_the_command_line_exits_nonzero_on_a_difference_and_writes_the_report(self):
        out = os.path.join(self.work, 'report.json')
        with contextlib.redirect_stdout(io.StringIO()):
            self.assertEqual(compare_runs.main([self.a, self.b, '--out', out]), 0)
            with open(out) as f:
                self.assertTrue(json.load(f)['identical'])
            d = self.copy('cli_missing')
            os.remove(os.path.join(d, 'sv_calls.parquet'))
            self.assertEqual(compare_runs.main([self.a, d]), 1)


if __name__ == '__main__':
    unittest.main()
