import copy
import unittest

from helpers import ROOT  # noqa: F401  (puts variant_index on the path)
import check_requirements as cr  # noqa: E402


class TestRequirementsMap(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.req = cr.load()

    def test_every_clause_is_mapped_to_existing_code_and_tests(self):
        probs = cr.check(self.req)
        self.assertEqual(probs, [])

    def test_the_rules_are_able_to_fail(self):
        for mutate, expect in (
            (lambda r: r['requirements'][0].update(code=['no/such/file.py']), 'does not exist'),
            (lambda r: r['requirements'][0].update(tests=['test_no_such_module']), 'does not exist'),
            (lambda r: r['requirements'][0].update(status='done'), 'vocabulary'),
            (lambda r: next(x for x in r['requirements'] if x['status'] == 'partial').update(gap=''), 'without a gap'),
            (lambda r: r['requirements'][0].update(gap='something'), 'verified but states a gap'),
            (lambda r: r['requirements'].append(copy.deepcopy(r['requirements'][0])), 'duplicate id'),
        ):
            bad = copy.deepcopy(self.req)
            mutate(bad)
            self.assertTrue(any(expect in p for p in cr.check(bad)), expect)

    def test_the_rendered_table_names_every_requirement(self):
        text = cr.render(self.req)
        for r in self.req['requirements']:
            self.assertIn(f"| {r['id']} |", text)


if __name__ == '__main__':
    unittest.main()
