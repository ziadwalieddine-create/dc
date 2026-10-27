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
                                                                      'the structural run on real calls', 'two matcher runs of one input', 'the index was built by committed code'))]
            self.assertEqual(len(ga.CHECKS), 13)
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

    def test_the_digest_of_committed_build_modules_is_the_digest_the_runner_computes(self):
        import subprocess
        from vindex.build import runner
        ga = self.ga
        paths = ['variant_index/vindex/build'] + ['variant_index/vindex/' + f for f in ga.BUILD_CODE_FIXED]
        r = subprocess.run(['git', 'status', '--porcelain', '--'] + paths, cwd=ga.ROOT, capture_output=True, text=True)
        if r.returncode or r.stdout.strip():
            self.skipTest('not a git checkout, or the build modules have uncommitted changes')
        want = runner.build_code()
        self.assertEqual(ga.build_code_at('HEAD'), (want['sha256_16'], want['files']))
        self.assertEqual(ga.commit_with_build_code(want['sha256_16']), ('HEAD', want['files']))
        self.assertIsNone(ga.build_code_at('0' * 40))

    def test_an_index_built_by_code_that_is_in_no_commit_is_a_failure_and_the_idempotence_check_must_have_examined_the_same_code(self):
        import json
        import shutil
        import subprocess
        import tempfile
        from vindex.build import runner
        ga = self.ga
        if subprocess.run(['git', 'rev-parse', 'HEAD'], cwd=ga.ROOT, capture_output=True).returncode:
            self.skipTest('not a git checkout')
        here = runner.build_code()['sha256_16']
        work = tempfile.mkdtemp()
        saved_root, saved_find = ga.ROOT, ga.commit_with_build_code
        try:
            idx = os.path.join(work, 'index')
            os.makedirs(idx)
            os.makedirs(os.path.join(work, 'data', 'idempotence'))

            def write_manifest(code):
                with open(os.path.join(idx, 'index_manifest.json'), 'w') as f:
                    json.dump(dict(build_code=code) if code else {}, f)

            def write_report(**kw):
                with open(os.path.join(work, 'data', 'idempotence', 'idempotence.json'), 'w') as f:
                    json.dump(dict(identical=True, files=3, limit=10, **kw), f)
            committed = next(f for n, f in ga.CHECKS if n.startswith('the index was built by committed code'))
            idem = next(f for n, f in ga.CHECKS if n.startswith('two builds from the same inputs'))
            a = argparse.Namespace(index=idx, results=None, no_tests=True)
            ga.commit_with_build_code = lambda want, max_commits=None: saved_find(want, max_commits=3)          # a search of the whole history is not needed to see a miss
            write_manifest(None)
            with self.assertRaises(AssertionError) as cm:
                committed(a)
            self.assertIn('records no build_code', str(cm.exception))
            write_manifest(dict(sha256_16='f' * 16, files=24))
            with self.assertRaises(AssertionError) as cm:
                committed(a)
            self.assertIn('uncommitted or unknown code', str(cm.exception))
            write_manifest(dict(sha256_16=here, files=24))
            self.assertIn('commit', committed(a))
            ga.ROOT = work
            write_report()
            with self.assertRaises(AssertionError) as cm:
                idem(a)
            self.assertIn('predates that field', str(cm.exception))
            write_report(build_code=dict(sha256_16='f' * 16, files=24))
            with self.assertRaises(AssertionError) as cm:
                idem(a)
            self.assertIn('examined build code', str(cm.exception))
            write_report(build_code=dict(sha256_16=here, files=24))
            self.assertIn(here, idem(a))
        finally:
            ga.ROOT, ga.commit_with_build_code = saved_root, saved_find
            shutil.rmtree(work, ignore_errors=True)

    def test_the_markdown_report_counts_every_outcome(self):
        ga = self.ga
        text = ga.markdown([dict(check='a', status='PASS', detail='x'), dict(check='b', status='FAIL', detail='y|z'), dict(check='c', status='BLOCKED', detail='w')], None)
        self.assertIn('1 pass, 1 fail, 0 skipped, 1 blocked', text)
        self.assertIn('y/z', text)


if __name__ == '__main__':
    unittest.main()
