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
                                                                      'every unresolved high-value locus', 'the sample run accounts', 'the sample run was made by committed code',
                                                                      'every release that holds evidence rows', 'lineage keys are shared', 'facts established outside the index',
                                                                      'the structural run on real calls', 'two matcher runs of one input'))]
            self.assertEqual(len(ga.CHECKS), 12)
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

    def test_a_run_made_by_code_that_has_changed_or_against_another_index_is_a_failure(self):
        import hashlib
        import json
        import shutil
        import subprocess
        import tempfile
        ga = self.ga
        fn = next(f for n, f in ga.CHECKS if n.startswith('the sample run was made by committed code'))
        head = subprocess.run(['git', 'rev-parse', 'HEAD'], cwd=ga.ROOT, capture_output=True, text=True).stdout.strip()
        if not head:
            self.skipTest('not a git checkout')
        work = tempfile.mkdtemp()
        try:
            idx, res = os.path.join(work, 'index'), os.path.join(work, 'results')
            os.makedirs(idx)
            os.makedirs(res)
            sha = lambda b: hashlib.sha256(b).hexdigest()
            for name in ('index_manifest.json', 'table_digests.json'):
                with open(os.path.join(idx, name), 'w') as f:
                    f.write('{}')
            with open(os.path.join(ga.VI, 'vindex', 'wgs.py'), 'rb') as f:
                wgs_sha = sha(f.read())
            good = dict(code=dict(commit=head, uncommitted_changes_in_variant_index=False), code_changed_during_run=False, code_files={'vindex/wgs.py': wgs_sha},
                        index={'index_manifest.json_sha256': sha(b'{}'), 'table_digests.json_sha256': sha(b'{}')})

            def run(**change):
                m = json.loads(json.dumps(good))
                m.update(change)
                with open(os.path.join(res, 'run_manifest.json'), 'w') as f:
                    json.dump(m, f)
                return fn(argparse.Namespace(index=idx, results=res, no_tests=True))
            self.assertIn('1 modules unchanged', run())
            for change, text in ((dict(code=dict(commit='0' * 40, uncommitted_changes_in_variant_index=False)), 'not in this repository'),
                                 (dict(code=dict(commit=head, uncommitted_changes_in_variant_index=True)), 'uncommitted changes'),
                                 (dict(code_changed_during_run=True), 'changed during the run'),
                                 (dict(code_changed_during_run=None), 'changed during the run'),
                                 (dict(code_files={}), 'no per-module digests'),
                                 (dict(code_files={'vindex/wgs.py': '0' * 64}), 'vindex/wgs.py'),
                                 (dict(code_files={'vindex/no_longer_exists.py': wgs_sha}), 'vindex/no_longer_exists.py'),
                                 (dict(index={'index_manifest.json_sha256': '0' * 64, 'table_digests.json_sha256': sha(b'{}')}), 'index_manifest.json of the index')):
                with self.subTest(text):
                    with self.assertRaises(AssertionError) as cm:
                        run(**change)
                    self.assertIn(text, str(cm.exception))
        finally:
            shutil.rmtree(work, ignore_errors=True)

    def test_the_markdown_report_counts_every_outcome(self):
        ga = self.ga
        text = ga.markdown([dict(check='a', status='PASS', detail='x'), dict(check='b', status='FAIL', detail='y|z'), dict(check='c', status='BLOCKED', detail='w')], None)
        self.assertIn('1 pass, 1 fail, 0 skipped, 1 blocked', text)
        self.assertIn('y/z', text)


if __name__ == '__main__':
    unittest.main()
