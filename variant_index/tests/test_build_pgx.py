import json
import os
import sys
import tempfile
import unittest

import pyarrow.compute as pc

from helpers import syn_reference, vcf_text, FIXTURES, SYN_FASTA, SYN_BUILDS
from vindex.build import runner
from vindex.build.pgx import CpicAdapter, PharmVarAdapter, cpic_variants
from vindex.index_reader import Index
from vindex.match_hap import call_gene
from vindex.match_small import MatchContext
from vindex.refgenome import Reference
from vindex.vcfio import read_vcf

sys.path.insert(0, FIXTURES)
import make_sources  # noqa: E402

SOURCES = make_sources.write_all(tempfile.mkdtemp())
b, alt = make_sources.b, make_sources.alt


class FakeGenes:
    """Stands in for GeneMap: no Open Targets data are needed for the unit test."""
    syms = {'TESTG': 'ENSG_TESTG', 'TESTH': 'ENSG_TESTH'}

    def resolve(self, hgnc=None, symbol=None, ncbi=None):
        return (self.syms[symbol], 'symbol') if symbol in self.syms else (None, 'unresolved')


class TestCpicParser(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.ref = syn_reference()

    def p(self, loc, allele):
        return cpic_variants(self.ref, '1', loc, allele)

    def test_substitution_variant_reference_and_multichange(self):
        r, a = b('1', 705), alt('1', 705)
        self.assertEqual(self.p(f'g.705{r}>{a}', a), ('variant', [f'1:705:{r}:{a}'], None))
        self.assertEqual(self.p(f'g.705{r}>{a}', r), ('reference', [], None))
        other = next(x for x in 'ACGT' if x not in (r, a))
        loc = f'g.705{r}>{a}; g.705{r}>{other}'
        self.assertEqual(self.p(loc, other)[1], [f'1:705:{r}:{other}'])
        self.assertEqual(self.p(f'g.705{r}>{a};g.705{r}>{other}', r)[0], 'reference')       # separator without a space
        self.assertEqual(self.p(loc, 'T' if 'T' not in (a, other, r) else 'A')[0], 'unresolved')

    def test_reference_base_conflict_is_never_guessed(self):
        wrong = 'C' if b('1', 705) != 'C' else 'G'
        st, keys, why = self.p(f'g.705{wrong}>T', 'T')
        self.assertEqual((st, keys), ('unresolved', []))
        self.assertTrue(why.startswith('reference_base_differs'))

    def test_iupac_tolerated_positions(self):
        r, a = b('1', 710), alt('1', 710)
        code = make_sources.IUPAC2[frozenset((r, a))]
        self.assertEqual(self.p(f'g.710{r}>{a}', code), ('tolerated', [f'1:710:{r}:{a}'], None))
        # a code that does not contain the reference base, or that allows a base nobody listed, is not interpreted
        excl = next(v for k, v in make_sources.IUPAC2.items() if r not in k)
        self.assertEqual(self.p(f'g.710{r}>{a}', excl)[2], 'iupac_code_excludes_reference')
        self.assertEqual(self.p(f'g.710{r}>{a}', 'N')[2], 'iupac_code_outside_listed_changes')

    def test_deletion_notations_agree(self):
        d = f"1:705:{b('1', 705)}{b('1', 706)}{b('1', 707)}:{b('1', 705)}"
        self.assertEqual(self.p('g.706_707del', 'delCC'), ('variant', [d], None))
        self.assertEqual(self.p('g.706_707delCC', 'delCC'), ('variant', [d], None))
        self.assertEqual(self.p('g.706_707del', 'CC'), ('reference', [], None))
        self.assertEqual(self.p('g.706_707delAA', 'delAA')[2], 'deleted_sequence_differs_from_reference')

    def test_duplication_and_insertion_normalise_to_one_key(self):
        """1:102-109 is an A8 run: a duplicated A and an inserted A anywhere in the run are the same variant once left-aligned."""
        self.assertEqual(self.p('g.105dup', 'AA'), ('variant', ['1:101:C:CA'], None))
        self.assertEqual(self.p('g.105dup', 'A'), ('reference', [], None))
        self.assertEqual(self.p('g.104_105insA', 'insA'), ('variant', ['1:101:C:CA'], None))
        self.assertEqual(self.p('g.104_105insA', 'del'), ('reference', [], None))
        self.assertEqual(self.p('g.104_106insA', 'insA')[2], 'insertion_flanks_not_adjacent')

    def test_delins(self):
        self.assertEqual(self.p('g.404_405delinsTT', 'insTT')[1], ['1:404:AC:TT'])
        self.assertEqual(self.p('g.404_405delinsTT', 'delAC')[0], 'reference')

    def test_unevaluable_notations_are_reported_not_guessed(self):
        self.assertEqual(self.p('g.200TCAG[5]', 'TCAG(5)'), ('unresolved', [], 'repeat_length_notation'))
        self.assertEqual(self.p('g.200TCAG[5]', 'TCAG'), ('unresolved', [], 'repeat_length_notation'))
        self.assertEqual(self.p('g.705_706inv', 'x')[2], 'unsupported_notation')
        self.assertEqual(self.p('', 'A')[2], 'empty_location')
        self.assertEqual(self.p('p.Arg1Gly', 'A')[2], 'unparsed_notation')
        self.assertEqual(self.p('g.99999A>C', 'C')[2], 'outside_contig')
        self.assertEqual(self.p('g.705G>A', 'zz')[2], 'allele_string_matches_neither_reference_nor_a_change')


class TestPgxAdapters(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.ref = Reference(SYN_FASTA, assembly_report=os.path.join(SOURCES, 'assembly_report.txt'), name='SYN38', builds=SYN_BUILDS)
        cls.out = tempfile.mkdtemp()
        cls.man = runner.build(cls.out, [CpicAdapter(FakeGenes(), root=os.path.join(SOURCES, 'cpic')), PharmVarAdapter(os.path.join(SOURCES, 'pharmvar-test.zip'))], cls.ref)
        cls.idx = Index.open(cls.out)
        cls.haps = {h['haplotype_id']: h for h in cls.idx.table('entity_haplotype').to_pylist()}
        cls.defs = {}
        for d in cls.idx.table('haplotype_definition').to_pylist():
            cls.defs.setdefault(d['haplotype_id'], []).append(d)
        cls.stats = {a['source_id']: a['stats'] for a in cls.man['adapters']}

    def scopes(self, hid):
        return sorted((d['definition_scope'], d['variant_key']) for d in self.defs.get(hid, ()))

    def test_cpic_entities_and_definitions(self):
        v705 = f"1:705:{b('1', 705)}:{alt('1', 705)}"
        v710 = f"1:710:{b('1', 710)}:{alt('1', 710)}"
        self.assertEqual(self.scopes('cpic:TESTG:*1'), [])
        self.assertEqual(self.scopes('cpic:TESTG:*2'), sorted([('defining', v705), ('tolerated', v710)]))
        self.assertEqual(self.scopes('cpic:TESTG:*3'), [('defining', v710)])
        self.assertTrue(self.haps['cpic:TESTG:*1']['is_reference'])
        self.assertEqual(self.haps['cpic:TESTG:*2']['allele_level'], 'cpic')
        self.assertEqual(self.haps['cpic:TESTG:*2']['function'], 'Decreased function')

    def test_unevaluable_cpic_alleles_are_kept_with_reasons(self):
        reasons = {n: [s for s, k in self.scopes(f'cpic:TESTH:{n}')] for n in ('*3', '*4', '*5', '*6')}
        self.assertEqual(reasons['*3'], ['unresolved:repeat_length_notation'])
        self.assertEqual(reasons['*4'], ['unresolved:no_location_rows'])
        self.assertEqual(reasons['*5'], ['unresolved:iupac_code_excludes_reference'])
        self.assertEqual(reasons['*6'], ['unresolved:reference_base_differs'])
        self.assertTrue(all(k is None for n in reasons for s, k in self.scopes(f'cpic:TESTH:{n}')))
        self.assertEqual(self.stats['cpic']['alleles_without_location_rows'], 1)

    def test_gene_without_chromosome_and_mitochondrial_gene(self):
        self.assertEqual(self.scopes('cpic:NOCHR:*2'), [('unresolved:gene_chromosome_unknown', None)])
        self.assertEqual(self.scopes('cpic:MTG:*2'), [('defining', f"MT:50:{b('MT', 50)}:{alt('MT', 50)}")])

    def test_allele_function_is_haplotype_level_evidence(self):
        ev = [r for r in self.idx.table('evidence').to_pylist() if r['source_id'] == 'cpic' and r['level'] == 'haplotype']
        by = {r['subject_id']: r for r in ev}
        self.assertEqual(by['cpic:TESTG:*1']['effect_assertion'], 'no_effect')
        self.assertEqual(by['cpic:TESTG:*2']['effect_assertion'], 'effect')
        self.assertEqual(by['cpic:TESTG:*2']['direction'], 'loss_of_function')
        self.assertEqual(by['cpic:TESTG:*2']['claim_num'], 0.5)
        self.assertEqual(by['cpic:TESTH:*3']['effect_assertion'], 'uncertain')
        self.assertIsNone(by['cpic:TESTG:*2']['chrom'])
        self.assertTrue(all(r['axis'] == 'effect_magnitude_direction' and r['mapping_status'] == 'mapped' for r in ev))
        self.assertNotIn('cpic:TESTH:*5', by)                  # no function stated: no row invented
        self.assertEqual(self.stats['cpic']['allele_function_rows'], 8)
        self.assertEqual(self.stats['cpic']['allele_rows_without_function'], 5)      # TESTH *5 and *6, NOCHR *1 and *2, MTG *1

    def test_recommendations_stay_at_gene_level_and_unmapped_genes_are_kept(self):
        rec = [r for r in self.idx.table('evidence').to_pylist() if r['claim_type'] == 'cpic_drug_recommendation']
        self.assertEqual(len(rec), 4)               # rec 1: one gene; rec 2: two genes; rec 3: one gene that is not in the gene map
        self.assertTrue(all(r['level'] == 'gene' and r['effect_assertion'] == 'not_applicable' for r in rec))
        two = [r for r in rec if r['accession'].startswith('rec:2:')]
        self.assertEqual(sorted(r['subject_id'] for r in two), ['ENSG_TESTG', 'ENSG_TESTH'])
        self.assertEqual({r['lineage_group'] for r in two}, {'cpic_guideline:1'})            # one guideline, counted once
        self.assertEqual(json.loads(two[0]['extra_json'])['combination_genes'], ['TESTG', 'TESTH'])
        gene_specific = {r['subject_id']: json.loads(r['extra_json'])['gene_specific']['phenotypes'] for r in two}
        self.assertEqual(gene_specific, {'ENSG_TESTG': 'Normal Metabolizer', 'ENSG_TESTH': 'Poor Metabolizer'})
        un = [r for r in rec if r['mapping_status'] == 'unmapped']
        self.assertEqual([(r['unmapped_key'], r['mapping_reason'], r['subject_id']) for r in un], [('UNKNOWNGENE', 'gene_symbol_not_in_gene_map', None)])

    def test_pharmvar_haplotypes_levels_and_definitions(self):
        pv = {k: v for k, v in self.haps.items() if v['system'] == 'pharmvar'}
        self.assertEqual(sorted(pv), sorted(f'pharmvar:TESTG:{n}' for n in ('*1', '*1.001', '*2', '*2.001', '*3', '*4', '*5', '*6', '*7')))
        self.assertEqual({v['allele_name']: v['allele_level'] for v in pv.values()}['*2.001'], 'suballele')
        self.assertEqual(pv['pharmvar:TESTG:*2']['allele_level'], 'core')
        self.assertTrue(pv['pharmvar:TESTG:*1']['is_reference'] and pv['pharmvar:TESTG:*1.001']['is_reference'])
        v705, v712 = f"1:705:{b('1', 705)}:{alt('1', 705)}", f"1:712:{b('1', 712)}:{alt('1', 712)}"
        self.assertEqual(self.scopes('pharmvar:TESTG:*2'), [('defining', v705)])
        self.assertEqual(self.scopes('pharmvar:TESTG:*2.001'), sorted([('defining', v705), ('defining', v712)]))
        self.assertEqual(self.scopes('pharmvar:TESTG:*4'), [('defining', '1:705:GCC:G')])
        self.assertEqual(self.scopes('pharmvar:TESTG:*5'), [('defining', '1:712:G:GTT')])
        self.assertEqual(self.scopes('pharmvar:TESTG:*6'), [('unresolved:deleted_sequence_differs_from_reference', None)])
        self.assertEqual(self.scopes('pharmvar:TESTG:*7'), [('unresolved:reference_sequence_unmapped', None)])
        self.assertEqual(self.stats['pharmvar']['haplotypes'], 9)           # the GRCh37 file is not read

    def test_star_allele_names_link_the_two_systems(self):
        x = [r for r in self.idx.table('xref').to_pylist() if r['namespace'] == 'star_allele_name' and r['identifier'] == 'TESTG*2']
        self.assertEqual(sorted(r['entity_id'] for r in x), ['cpic:TESTG:*2', 'pharmvar:TESTG:*2'])
        self.assertFalse([r for r in self.idx.table('xref').to_pylist() if r['identifier'] == 'TESTG*2.001' and r['namespace'] == 'star_allele_name'])

    def test_entity_variants_cover_resolved_definitions_only(self):
        keys = set(self.idx.table('entity_variant')['variant_key'].to_pylist())
        defined = {d['variant_key'] for ds in self.defs.values() for d in ds if d['variant_key']}
        self.assertEqual(keys, defined)
        self.assertIn('1:705:GCC:G', keys)

    def test_rebuild_is_idempotent(self):
        out2 = tempfile.mkdtemp()
        runner.build(out2, [CpicAdapter(FakeGenes(), root=os.path.join(SOURCES, 'cpic')), PharmVarAdapter(os.path.join(SOURCES, 'pharmvar-test.zip'))], self.ref)
        with open(os.path.join(self.out, 'table_digests.json')) as f1, open(os.path.join(out2, 'table_digests.json')) as f2:
            self.assertEqual(json.load(f1), json.load(f2))


def unph(gt):
    return gt, ('', '')


class TestPgxMatching(unittest.TestCase):
    """CPIC alleles built from the mini source: TESTG *1 reference, *2 defined by 705 and tolerating the alternative at 710, *3 defined by 710."""

    @classmethod
    def setUpClass(cls):
        cls.ref = Reference(SYN_FASTA, assembly_report=os.path.join(SOURCES, 'assembly_report.txt'), name='SYN38', builds=SYN_BUILDS)
        cls.out = tempfile.mkdtemp()
        runner.build(cls.out, [CpicAdapter(FakeGenes(), root=os.path.join(SOURCES, 'cpic')), PharmVarAdapter(os.path.join(SOURCES, 'pharmvar-test.zip'))], cls.ref)
        cls.idx = Index.open(cls.out)
        cls.td = tempfile.mkdtemp()

    def call(self, genos, blocks=(), gene='TESTG', **kw):
        recs = []
        for pos, (gt, extra) in genos.items():
            recs.append(('1', pos, b('1', pos), alt('1', pos), 'PASS', 50, [f'{gt}:50:30:15,15' + extra[0]], '.', '.', 'GT:GQ:DP:AD' + extra[1]))
        for pos in blocks:
            recs.append(('1', pos, b('1', pos), '<NON_REF>', 'PASS', 0, ['0/0:50:30:30,0'], f'END={pos}', '.', 'GT:GQ:DP:AD'))
        recs.sort(key=lambda r: r[1])
        path = os.path.join(self.td, 'p.vcf')
        with open(path, 'w') as f:
            f.write(vcf_text(recs))
        return call_gene(MatchContext(read_vcf(path, self.ref)), self.idx, gene, **kw)

    def test_tolerated_position_widens_the_candidates_instead_of_forcing_one(self):
        c = self.call({705: unph('0/1'), 710: unph('0/1')}, system='cpic')
        # *2 may carry the alternative at 710, so *1/*2 is as consistent as *2/*3; without the tolerance only *2/*3 would remain
        self.assertEqual((c.status, c.diplotypes), ('AMBIGUOUS', [('*1', '*2'), ('*2', '*3')]))

    def test_tolerated_position_with_reference_genotype_still_resolves(self):
        c = self.call({705: unph('0/1')}, blocks=(710,), system='cpic')
        self.assertEqual((c.status, c.diplotypes), ('RESOLVED', [('*1', '*2')]))
        c = self.call({710: unph('0/1')}, blocks=(705,), system='cpic')
        self.assertEqual((c.status, c.diplotypes), ('RESOLVED', [('*1', '*3')]))

    def test_phase_resolves_what_the_tolerance_made_ambiguous(self):
        cis = self.call({705: ('0|1', (':1', ':PS')), 710: ('0|1', (':1', ':PS'))}, system='cpic')
        self.assertEqual((cis.status, cis.diplotypes, cis.basis), ('RESOLVED', [('*1', '*2')], 'phased'))
        trans = self.call({705: ('0|1', (':1', ':PS')), 710: ('1|0', (':1', ':PS'))}, system='cpic')
        self.assertEqual((trans.status, trans.diplotypes), ('RESOLVED', [('*2', '*3')]))

    def test_missing_defining_position_blocks_the_call(self):
        c = self.call({705: unph('0/1')}, system='cpic')
        self.assertEqual(c.status, 'UNRESOLVED')
        self.assertEqual(c.missing['positions'], {f"1:710:{b('1', 710)}:{alt('1', 710)}": 'ABSENT_FROM_CALLSET'})

    def test_unevaluable_alleles_make_the_call_conditional(self):
        c = self.call({402: unph('0/1')}, gene='TESTH', system='cpic')
        self.assertEqual((c.status, c.diplotypes), ('RESOLVED_CONDITIONAL', [('*1', '*2')]))
        self.assertEqual(c.unevaluable_alleles, {'*3': ['repeat_length_notation'], '*4': ['no_location_rows'], '*5': ['iupac_code_excludes_reference'], '*6': ['reference_base_differs']})
        self.assertEqual(set(c.missing['not_excluded_alleles']), {'*3', '*4', '*5', '*6'})
        # ambiguity keeps its status but still reports what could not be evaluated
        self.assertIn('unevaluable_allele_definitions', self.call({}, gene='TESTH', system='cpic').missing)

    def test_system_selection_and_levels(self):
        c = self.call({705: unph('0/1')}, blocks=(710,), system=None)
        self.assertEqual((c.system, c.status, c.diplotypes), ('cpic', 'RESOLVED', [('*1', '*2')]))
        self.assertIn('system_defaulted:cpic_of:cpic,pharmvar', c.flags)
        # PharmVar is called at the core level by default (a sub-allele repeats its core allele's variants); *6 and *7 cannot be evaluated
        pv = self.call({710: unph('0/1')}, blocks=(705, 706, 707, 712, 713), system='pharmvar')
        self.assertEqual((pv.system, pv.status, pv.diplotypes), ('pharmvar', 'RESOLVED_CONDITIONAL', [('*1', '*3')]))
        self.assertEqual(set(pv.unevaluable_alleles), {'*6', '*7'})
        sub = self.call({}, blocks=(705, 712), system='pharmvar', levels=('suballele',))
        self.assertEqual((sub.status, sub.diplotypes), ('RESOLVED', [('*1.001', '*1.001')]))

    def test_unknown_gene_and_level_without_alleles(self):
        self.assertEqual(self.call({}, gene='NOSUCH').basis, 'no_haplotype_definitions_in_index')
        self.assertEqual(self.call({}, system='cpic', levels=('core',)).basis, 'no_haplotype_definitions_at_requested_level')


if __name__ == '__main__':
    unittest.main()
