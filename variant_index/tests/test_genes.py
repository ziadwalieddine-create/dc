import unittest

from helpers import ROOT  # noqa: F401  (puts the package on the path)
from vindex.build.genes import entity_id


class Map:
    """A gene map that resolves HGNC id, then NCBI id, then symbol, like GeneMap.resolve."""
    by = {('hgnc', 'HGNC:5'): 'ENSG5', ('ncbi', '7'): 'ENSG7', ('symbol', 'ABC'): 'ENSG9'}

    def resolve(self, hgnc=None, symbol=None, ncbi=None):
        for key, v in (('hgnc', hgnc), ('ncbi', ncbi), ('symbol', symbol)):
            if v and (key, v) in self.by:
                return self.by[(key, v)], key
        return None, 'unresolved'


class TestEntityId(unittest.TestCase):
    def test_resolved_genes_use_the_ensembl_id_whatever_identifier_the_source_had(self):
        self.assertEqual(entity_id(Map(), hgnc='HGNC:5'), 'ENSG5')
        self.assertEqual(entity_id(Map(), ncbi='7'), 'ENSG7')
        self.assertEqual(entity_id(Map(), symbol='ABC'), 'ENSG9')
        self.assertEqual(entity_id(Map(), symbol='OLD', ncbi='7'), 'ENSG7')                     # the stronger identifier wins over an unknown symbol

    def test_unresolved_genes_keep_their_own_identifier_with_a_namespace(self):
        self.assertEqual(entity_id(Map(), symbol='NOPE'), 'symbol:NOPE')
        self.assertEqual(entity_id(Map(), ncbi='99'), 'NCBIGene:99')
        self.assertEqual(entity_id(Map(), hgnc='HGNC:1', symbol='X'), 'HGNC:HGNC:1')
        self.assertEqual(entity_id(None, ncbi='99'), 'NCBIGene:99')
        self.assertIsNone(entity_id(Map()))


if __name__ == '__main__':
    unittest.main()
