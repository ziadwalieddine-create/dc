import os
import shutil
import sys
import tempfile
import unittest

from helpers import syn_reference, FIXTURES
from vindex.build import runner
from vindex.build.pathways import ReactomeAdapter
from vindex.index_reader import Index

sys.path.insert(0, FIXTURES)
import make_sources  # noqa: E402

SOURCES = make_sources.write_all(tempfile.mkdtemp())


class Genes:
    by_ensg = {'ENSG_G1': 1, 'ENSG_G2': 1}


class TestReactome(unittest.TestCase):
    def test_human_gene_memberships_only_and_every_skip_counted(self):
        out = tempfile.mkdtemp()
        try:
            man = runner.build(out, [ReactomeAdapter(Genes(), root=os.path.join(SOURCES, 'reactome'))], syn_reference())
            rows = Index.open(out).table('gene_set_member').to_pylist()
        finally:
            shutil.rmtree(out, ignore_errors=True)
        self.assertEqual(sorted((r['set_id'], r['gene_id']) for r in rows), [('R-HSA-1', 'ENSG_G1'), ('R-HSA-1', 'ENSG_G2'), ('R-HSA-2', 'ENSG_G1'), ('R-HSA-4', 'ENSG_ZZZ')])      # version suffix dropped
        self.assertTrue(all(r['kind'] == 'reactome_pathway' and r['source_id'] == 'reactome' for r in rows))
        st = man['adapters'][0]['stats']
        self.assertEqual((st['rows'], st['rows_of_other_species'], st['human_rows_not_keyed_by_gene'], st['duplicate_memberships'], st['memberships_of_genes_not_in_gene_map'], st['memberships']),
                         (7, 1, 1, 1, 1, 4))
        self.assertEqual(man['problems'], [])


if __name__ == '__main__':
    unittest.main()
