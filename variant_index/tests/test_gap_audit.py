import argparse
import importlib.util
import os
import unittest

from helpers import ROOT  # noqa: F401  (puts variant_index on the path)

PATH = os.path.join(ROOT, '..', 'inventory', 'gap_audit.py')


def load():
    import sys
    sys.path.insert(0, os.path.dirname(PATH))
    spec = importlib.util.spec_from_file_location('gap_audit_under_test', PATH)
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m


class TestGapAudit(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.ga = load()

    def args(self, **kw):
        return argparse.Namespace(index=os.path.join(ROOT, 'no_such_index'), results=os.path.join(ROOT, 'no_such_results'), no_tests=True, **kw)

    def test_statuses_distinguish_pass_fail_skip_and_blocked(self):
        ga = self.ga
        saved = list(ga.CHECKS)
        try:
            def ok(a): return 'fine'
            def bad(a): raise AssertionError('broken')
            def skip(a): raise ga.Skip('needs x')
            def blocked(a): raise ga.Blocked('owner only')
            def crash(a): raise RuntimeError('oops')
            ga.CHECKS[:] = [('ok', ok), ('bad', bad), ('skip', skip), ('blocked', blocked), ('crash', crash)]
            res = {r['check']: r['status'] for r in ga.run(self.args())}
        finally:
            ga.CHECKS[:] = saved
        self.assertEqual(res, {'ok': 'PASS', 'bad': 'FAIL', 'skip': 'SKIP', 'blocked': 'BLOCKED', 'crash': 'FAIL'})

    def test_checks_that_need_an_absent_artefact_are_skipped_never_passed(self):
        ga = self.ga
        saved = list(ga.CHECKS)
        try:
            ga.CHECKS[:] = [(n, f) for n, f in saved if n.startswith(('the built index passed validation', 'provenance:', 'a sample run agrees with bcftools', 'ClinVar-derived rows',
                                                                      'every unresolved high-value locus', 'the sample run accounts'))]
            self.assertEqual(len(ga.CHECKS), 6)
            res = ga.run(self.args())
        finally:
            ga.CHECKS[:] = saved
        self.assertEqual({r['status'] for r in res}, {'SKIP'})

    def test_a_broken_requirements_map_is_a_failure(self):
        import sys
        sys.path.insert(0, ROOT)
        import check_requirements as cr
        ga = self.ga
        real = cr.load
        try:
            def broken():
                d = real()
                d['requirements'][0]['code'] = ['no/such/file.py']
                return d
            cr.load = broken
            fn = next(f for n, f in ga.CHECKS if n.startswith('REQUIREMENTS.md is the rendering'))
            with self.assertRaises(AssertionError):
                fn(self.args())
        finally:
            cr.load = real

    def test_the_markdown_report_counts_every_outcome(self):
        ga = self.ga
        text = ga.markdown([dict(check='a', status='PASS', detail='x'), dict(check='b', status='FAIL', detail='y|z'), dict(check='c', status='BLOCKED', detail='w')], None)
        self.assertIn('1 pass, 1 fail, 0 skipped, 1 blocked', text)
        self.assertIn('y/z', text)


if __name__ == '__main__':
    unittest.main()
