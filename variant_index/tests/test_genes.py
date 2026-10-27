import os
import shutil
import tempfile
import unittest

import pyarrow as pa
import pyarrow.parquet as pq

from helpers import ROOT  # noqa: F401  (puts the package on the path)
from vindex.build import clinvar_vcf, clinvar_vcv
from vindex.build.genes import GeneMap, entity_id


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


class TestOneGeneUnderSeveralSymbols(unittest.TestCase):
    """ClinVar writes some NCBI Gene ids under two symbols (CFAP65-AS1:100129175 and LOC100129175:100129175). Resolved record by record, the same gene got an Ensembl id from one
    record and an NCBIGene: id from another (45 genes on the real file)."""

    def genes(self):
        g = GeneMap(None, ot_dir='/nonexistent', mane_path='/nonexistent')
        g.by_symbol = {'NEWSYM': dict(gene_id='ENSG1'), 'OTHER': dict(gene_id='ENSG2')}
        g.by_ncbi = {'5': dict(gene_id='ENSG5')}
        return g

    def both(self, g, ncbi='100'):
        return entity_id(g, symbol='OLDSYM', ncbi=ncbi), entity_id(g, symbol='NEWSYM', ncbi=ncbi)

    def test_the_resolution_does_not_depend_on_which_symbol_a_record_carries(self):
        g = self.genes()
        self.assertEqual(self.both(g), ('NCBIGene:100', 'ENSG1'))                       # before: the two records of one gene disagree
        g.learn_ncbi_symbols([('100', 'OLDSYM'), ('100', 'NEWSYM'), ('100', 'NEWSYM')])
        self.assertEqual(self.both(g), ('ENSG1', 'ENSG1'))
        self.assertEqual(g.resolve(symbol='OLDSYM', ncbi='100'), ('ENSG1', 'ncbi_symbols'))

    def test_symbols_that_name_different_genes_leave_the_gene_unresolved_for_every_record(self):
        g = self.genes()
        g.learn_ncbi_symbols([('100', 'NEWSYM'), ('100', 'OTHER')])
        self.assertEqual(self.both(g), ('NCBIGene:100', 'NCBIGene:100'))

    def test_the_ncbi_table_still_decides_first_and_a_single_symbol_resolves_as_before(self):
        g = self.genes()
        g.learn_ncbi_symbols([('5', 'NEWSYM'), ('5', 'OTHER'), ('100', 'NEWSYM')])
        self.assertEqual(entity_id(g, symbol='OTHER', ncbi='5'), 'ENSG5')
        self.assertEqual(entity_id(g, symbol='NEWSYM', ncbi='100'), 'ENSG1')
        self.assertEqual(entity_id(g, symbol='OLDSYM', ncbi='100'), 'NCBIGene:100')           # one symbol known: the record's own symbol decides, as it did

    def test_the_pairs_are_read_from_the_vcf_info_column_and_the_vcv_variation_table(self):
        d = tempfile.mkdtemp()
        try:
            vcf = os.path.join(d, 'x.vcf')
            with open(vcf, 'w') as f:
                f.write('##fileformat=VCFv4.1\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n'
                        '1\t5\t1\tA\tC\t.\t.\tALLELEID=1;GENEINFO=CFAP65-AS1:100129175|X:7;ORIGIN=1\n'
                        '1\t6\t2\tA\tC\t.\t.\tALLELEID=2;GENEINFO=LOC100129175:100129175\n'
                        '1\t7\t3\tA\tC\t.\t.\tALLELEID=3\n')
            self.assertEqual(sorted(clinvar_vcf.gene_pairs(vcf)), [('100129175', 'CFAP65-AS1'), ('100129175', 'LOC100129175'), ('7', 'X')])
            os.makedirs(os.path.join(d, 'vcv', 'variation'))
            pq.write_table(pa.table({'genes': pa.array(['GENEA:1001|GENEB:1002', None, '', 'GENEA:1001'])}), os.path.join(d, 'vcv', 'variation', 'part.parquet'))
            self.assertEqual(sorted(clinvar_vcv.gene_pairs(os.path.join(d, 'vcv'))), [('1001', 'GENEA'), ('1001', 'GENEA'), ('1002', 'GENEB')])
        finally:
            shutil.rmtree(d, ignore_errors=True)


if __name__ == '__main__':
    unittest.main()
