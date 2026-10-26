import glob
import os
import shutil
import tempfile
import unittest

import pyarrow.parquet as pq

from helpers import ROOT  # noqa: F401  (puts variant_index on the path)
from test_wgs import build_index
import check_idempotence as ci  # noqa: E402


class TestIdempotence(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.work = tempfile.mkdtemp()
        cls.a, cls.b = os.path.join(cls.work, 'a'), os.path.join(cls.work, 'b')
        build_index(cls.a)
        build_index(cls.b)

    @classmethod
    def tearDownClass(cls):
        shutil.rmtree(cls.work, ignore_errors=True)

    def test_two_builds_from_the_same_inputs_have_identical_files(self):
        r = ci.compare(self.a, self.b)
        self.assertGreater(r['files'], 10)
        self.assertEqual((r['only_in_first'], r['only_in_second'], r['digest_differs']), ([], [], []))
        self.assertTrue(r['releases_equal'] and r['problems_equal'] and r['adapter_code_equal'] and r['rows_equal'])

    def test_a_changed_row_and_a_missing_file_are_reported_when_the_files_are_rehashed(self):
        d = os.path.join(self.work, 'c')
        shutil.copytree(self.b, d)
        victim = sorted(glob.glob(os.path.join(d, 'evidence', 'variant', '1', '*.parquet')))[0]
        t = pq.read_table(victim)
        pq.write_table(t.slice(0, t.num_rows - 1), victim)                                             # one row fewer
        gone = sorted(glob.glob(os.path.join(d, 'xref', '*.parquet')))[0]
        os.remove(gone)
        r = ci.compare(self.a, d)
        self.assertEqual(ci.compare(self.a, d)['digest_differs'], [])                                     # the digests recorded at build time cannot see a later edit...
        r = ci.compare(self.a, d, rehash=True)                                                          # ...the files themselves can
        self.assertEqual(r['digest_differs'], [os.path.relpath(victim, d)])
        self.assertEqual(r['only_in_first'], [os.path.relpath(gone, d)])


if __name__ == '__main__':
    unittest.main()
