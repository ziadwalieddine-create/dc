import os
import unittest

from helpers import ROOT  # noqa: F401  (puts variant_index on the path)
import make_status as ms  # noqa: E402
from vindex import reports


class TestStatusTable(unittest.TestCase):
    def test_the_block_between_the_markers_is_replaced_and_the_text_around_it_is_kept(self):
        text = f'before\n{ms.START}\nold table\n{ms.END}\nafter\n'
        new = ms.render(text)
        self.assertTrue(new.startswith('before\n'))
        self.assertTrue(new.endswith('after\n'))
        self.assertNotIn('old table', new)
        for m in reports.MODULES:
            self.assertIn(f"| {m['name']} | **{m['status']}** |", new)
        self.assertEqual(ms.render(new), new)                                                       # idempotent

    def test_the_committed_status_md_matches_the_registry(self):
        path = os.path.join(ROOT, 'STATUS.md')
        if not os.path.exists(path):
            self.skipTest('STATUS.md not written yet')
        with open(path) as f:
            text = f.read()
        self.assertEqual(ms.render(text), text, 'run: python3 variant_index/make_status.py --write')


if __name__ == '__main__':
    unittest.main()
