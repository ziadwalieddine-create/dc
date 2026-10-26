import os
import tempfile
import unittest

from helpers import syn_reference, vcf_text
from vindex import schema
from vindex.index_reader import Index
from vindex.vcfio import read_vcf
from vindex.match_small import MatchContext
from vindex.match_hap import call_gene


class TestHaplotypes(unittest.TestCase):
    """Synthetic pharmacogene TESTG with defining variants at 1:705, 1:710, 1:715.
    *1 reference; *2 = v705; *3 = v710; *4 = v705 + v710 (so v705/v710 het-het unphased fits *2/*3 and *1/*4); *5 = v715."""

    @classmethod
    def setUpClass(cls):
        cls.ref = syn_reference()
        cls.td = tempfile.mkdtemp()
        alt = lambda pos: next(x for x in 'ACGT' if x != cls.ref.fetch('1', pos, pos))
        cls.v = {p: f'1:{p}:{cls.ref.fetch("1", p, p)}:{alt(p)}' for p in (705, 710, 715)}
        H = lambda n, struct=False: dict(haplotype_id=f'cpic:TESTG:{n}', system='cpic', gene_symbol='TESTG', allele_name=n, is_reference=n == '*1', structural=struct,
                                         release_id='r', accession=n, function='')
        defs = {'*1': [], '*2': [705], '*3': [710], '*4': [705, 710], '*5': [715]}
        haps = [H(n) for n in defs]
        rows = [dict(haplotype_id=f'cpic:TESTG:{n}', variant_key=cls.v[p], hgvs='', release_id='r', definition_scope='gene') for n, ps in defs.items() for p in ps]
        cls.idx = Index.from_tables(entity_haplotype=schema.rows_to_table('entity_haplotype', haps), haplotype_definition=schema.rows_to_table('haplotype_definition', rows))

    def run_gene(self, genos, blocks=(), fmt='GT:GQ:DP:AD:PS'):
        """genos: {pos: GT string or None for unlisted}; blocks: positions given confident reference blocks."""
        recs = []
        for pos, gt in genos.items():
            if gt is not None:
                recs.append(('1', pos, self.ref.fetch('1', pos, pos), self.v[pos].split(':')[3], 'PASS', 50, [f'{gt}:50:30:15,15:9'], '.', '.', fmt))
        for pos in blocks:
            recs.append(('1', pos, self.ref.fetch('1', pos, pos), '<NON_REF>', 'PASS', 0, ['0/0:50:30:30,0'], f'END={pos}', '.', 'GT:GQ:DP:AD'))
        recs.sort(key=lambda r: r[1])
        p = os.path.join(self.td, 'h.vcf')
        with open(p, 'w') as f:
            f.write(vcf_text(recs))
        cs = read_vcf(p, self.ref)
        return call_gene(MatchContext(cs), self.idx, 'TESTG')

    def test_all_positions_known_unique(self):
        c = self.run_gene({705: '0/0', 710: '0/0', 715: '0/0'}, blocks=(705, 710, 715)) if False else self.run_gene({705: None, 710: None, 715: None}, blocks=(705, 710, 715))
        self.assertEqual((c.status, c.diplotypes), ('RESOLVED', [('*1', '*1')]))

    def test_heterozygous_single_variant(self):
        c = self.run_gene({705: '0/1', 710: None, 715: None}, blocks=(710, 715))
        self.assertEqual((c.status, c.diplotypes), ('RESOLVED', [('*1', '*2')]))

    def test_absent_positions_in_variant_only_vcf_block_the_call(self):
        c = self.run_gene({705: '0/1', 710: None, 715: None})
        self.assertEqual(c.status, 'UNRESOLVED')
        self.assertEqual(set(c.missing['positions']), {self.v[710], self.v[715]})
        self.assertTrue(all(s == 'ABSENT_FROM_CALLSET' for s in c.missing['positions'].values()))
        self.assertIn(('*1', '*2'), c.diplotypes)        # still consistent with what is known, but not a call

    def test_unphased_double_het_is_ambiguous(self):
        c = self.run_gene({705: '0/1', 710: '0/1', 715: None}, blocks=(715,))
        self.assertEqual(c.status, 'AMBIGUOUS')
        self.assertEqual(set(c.diplotypes), {('*1', '*4'), ('*2', '*3')})

    def test_phase_resolves_cis_and_trans(self):
        cis = self.run_gene({705: '0|1', 710: '0|1', 715: None}, blocks=(715,))
        self.assertEqual((cis.status, cis.diplotypes, cis.basis), ('RESOLVED', [('*1', '*4')], 'phased'))
        trans = self.run_gene({705: '0|1', 710: '1|0', 715: None}, blocks=(715,))
        self.assertEqual((trans.status, trans.diplotypes), ('RESOLVED', [('*2', '*3')]))

    def test_homozygous_combination(self):
        c = self.run_gene({705: '1/1', 710: '1/1', 715: None}, blocks=(715,))
        self.assertEqual((c.status, c.diplotypes), ('RESOLVED', [('*4', '*4')]))

    def test_undefined_phased_combination(self):
        # v705 and v715 on the same haplotype: no defined allele carries both
        c = self.run_gene({705: '0|1', 710: None, 715: '0|1'}, blocks=(710,))
        self.assertEqual(c.status, 'NO_CONSISTENT_DIPLOTYPE')

    def test_structural_gene_needs_a_dedicated_caller(self):
        H = dict(haplotype_id='cpic:CYP2D6:*5', system='cpic', gene_symbol='CYP2D6', allele_name='*5', is_reference=False, structural=True, release_id='r', accession='*5', function='')
        idx = Index.from_tables(entity_haplotype=schema.rows_to_table('entity_haplotype', [H]), haplotype_definition=schema.empty_table('haplotype_definition'))
        p = os.path.join(self.td, 'e.vcf')
        with open(p, 'w') as f:
            f.write(vcf_text([]))
        cs = read_vcf(p, self.ref)
        c = call_gene(MatchContext(cs), idx, 'CYP2D6')
        self.assertEqual(c.status, 'REQUIRES_SPECIALISED_CALLER')
        ext = call_gene(MatchContext(cs), idx, 'CYP2D6', external={'CYP2D6': ('*1', '*5')})
        self.assertEqual((ext.status, ext.flags), ('RESOLVED', ('not_derived_from_the_variant_calls',)))


if __name__ == '__main__':
    unittest.main()
