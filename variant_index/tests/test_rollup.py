"""Trait roll-up through the cross-references and retired terms of the ontology table: ids the table does not hold are placed through the terms that list or replace them, the mapping
is kept and its many-to-one cases are flagged, a prefix is part of the match, a direct term wins, and what maps to nothing is listed, not dropped."""
import os
import shutil
import tempfile
import unittest

import pyarrow as pa
import pyarrow.parquet as pq

import helpers  # noqa: F401  (puts variant_index on the path)
from vindex import views

SCHEMA = pa.schema([('id', pa.string()), ('name', pa.string()), ('ancestors', pa.list_(pa.string())), ('therapeuticAreas', pa.list_(pa.string())), ('dbXRefs', pa.list_(pa.string())),
                    ('obsoleteTerms', pa.list_(pa.string()))])
ROWS = [
    dict(id='EFO_0000200', name='area A', ancestors=['EFO_0000100'], therapeuticAreas=[], dbXRefs=[], obsoleteTerms=[]),
    dict(id='EFO_0000300', name='area B', ancestors=[], therapeuticAreas=[], dbXRefs=[], obsoleteTerms=[]),
    dict(id='EFO_0000001', name='trait one', ancestors=['EFO_0000200'], therapeuticAreas=['EFO_0000200'], dbXRefs=['UMLS:C0001', 'OMIM:100100', 'MESH:D000001', 'Orphanet:555'], obsoleteTerms=[]),
    dict(id='MONDO_0000002', name='trait two', ancestors=['EFO_0000200', 'EFO_0000300'], therapeuticAreas=['EFO_0000200', 'EFO_0000300'], dbXRefs=['UMLS:C0002', 'OMIM:100200', 'OMIMPS:300300'],
         obsoleteTerms=['EFO_0000777', 'Orphanet_999']),
    dict(id='MONDO_0000003', name='trait three', ancestors=['EFO_0000300'], therapeuticAreas=['EFO_0000300'], dbXRefs=['UMLS:C0002'], obsoleteTerms=[]),
    dict(id='MONDO_0000004', name='trait with no area', ancestors=[], therapeuticAreas=[], dbXRefs=['MESH:D000004', 'UMLS:C0004'], obsoleteTerms=[]),
    dict(id='Orphanet_100', name='an Orphanet term', ancestors=[], therapeuticAreas=['EFO_0000300'], dbXRefs=['OMIM:100100'], obsoleteTerms=[]),
]


def traits(*ids):
    return {i: dict(texts=[i], axes=['association_replication'], n_calls=1, calls=[f'call-{n}']) for n, i in enumerate(ids)}


class TestCrossReferences(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.dir = tempfile.mkdtemp()
        pq.write_table(pa.Table.from_pylist(ROWS, schema=SCHEMA), os.path.join(cls.dir, 'part-0.parquet'))
        cls.onto = views.ontology_terms(cls.dir)

    @classmethod
    def tearDownClass(cls):
        shutil.rmtree(cls.dir, ignore_errors=True)

    def roll(self, *ids):
        return views.trait_rollup(traits(*ids), self.onto)

    def test_the_table_is_read_with_its_cross_references_and_retired_terms(self):
        t = self.onto['MONDO:0000002']
        self.assertEqual((t['xrefs'], t['replaces']), (['OMIM:100200', 'OMIMPS:300300', 'UMLS:C0002'], ['EFO:0000777', 'Orphanet:999']))
        self.assertEqual(self.onto['EFO:0000200']['xrefs'], [])

    def test_a_medgen_concept_is_placed_through_the_umls_cross_reference(self):
        r = self.roll('MedGen:C0001')
        self.assertEqual(r['cross_referenced']['MedGen:C0001'], dict(terms=['EFO:0000001'], how=['dbXRef'], many_to_one=False))
        self.assertEqual({a: (g['n_traits'], g['n_traits_cross_referenced']) for a, g in r['by_therapeutic_area'].items()}, {'EFO:0000200': (1, 1)})
        self.assertEqual((r['not_rolled_up'], r['ancestors']['MedGen:C0001']), ([], ['EFO:0000200']))

    def test_omim_mesh_and_orphanet_ids_are_placed_by_their_own_prefix(self):
        r = self.roll('OMIM:100200', 'MeSH:D000001', 'OMIM:300300', 'Orphanet:555')
        self.assertEqual({k: v['terms'] for k, v in r['cross_referenced'].items()}, {'OMIM:100200': ['MONDO:0000002'], 'MeSH:D000001': ['EFO:0000001'], 'OMIM:300300': ['MONDO:0000002'],
                                                                                    'Orphanet:555': ['EFO:0000001']})            # an OMIM phenotypic series is listed as OMIMPS and still reaches its term
        self.assertEqual(r['not_rolled_up'], [])

    def test_a_retired_term_is_placed_through_the_term_that_replaces_it(self):
        r = self.roll('EFO:0000777', 'Orphanet:999')
        self.assertEqual({k: (v['terms'], v['how']) for k, v in r['cross_referenced'].items()}, {'EFO:0000777': (['MONDO:0000002'], ['replaces_a_retired_term']),
                                                                                                 'Orphanet:999': (['MONDO:0000002'], ['replaces_a_retired_term'])})

    def test_a_cross_reference_shared_by_two_terms_is_flagged_and_counted_in_the_areas_of_both(self):
        r = self.roll('MedGen:C0002')
        self.assertEqual(r['cross_referenced']['MedGen:C0002'], dict(terms=['MONDO:0000002', 'MONDO:0000003'], how=['dbXRef'], many_to_one=True))
        self.assertEqual(set(r['by_therapeutic_area']), {'EFO:0000200', 'EFO:0000300'})
        self.assertEqual(r['ancestors']['MedGen:C0002'], ['EFO:0000200', 'EFO:0000300'])

    def test_the_prefix_is_part_of_the_match(self):
        r = self.roll('OMIM:555', 'Orphanet:100100', 'MeSH:C0001', 'MedGen:D000001', 'OMIM:C0001', 'HP:100100')
        self.assertEqual(r['cross_referenced'], {})                                                          # 555 is an Orphanet cross-reference, 100100 an OMIM one, C0001 a UMLS one
        self.assertEqual(len(r['not_rolled_up']), 6)

    def test_a_trait_the_table_holds_is_placed_through_its_own_term_even_when_another_lists_it(self):
        r = self.roll('Orphanet:100')                                                                        # a term of the table; another term lists OMIM:100100, not this
        self.assertEqual(r['cross_referenced'], {})
        self.assertEqual(set(r['by_therapeutic_area']), {'EFO:0000300'})
        r = self.roll('EFO:0000001')
        self.assertEqual((r['cross_referenced'], r['by_therapeutic_area']['EFO:0000200']['n_traits_cross_referenced']), ({}, 0))

    def test_a_term_without_an_area_is_its_own_group_whether_reached_directly_or_by_cross_reference(self):
        direct = self.roll('MONDO:0000004')['by_therapeutic_area']
        crossed = self.roll('MeSH:D000004')['by_therapeutic_area']
        self.assertEqual(set(direct), {'MONDO:0000004'})
        self.assertEqual(set(crossed), {'MONDO:0000004'})

    def test_what_maps_to_nothing_is_listed_not_dropped_and_calls_are_unchanged(self):
        t = traits('MedGen:C9999', 'GeneReviews:NBK1', 'MedGen:C0001')
        r = views.trait_rollup(t, self.onto)
        self.assertEqual(r['not_rolled_up'], ['GeneReviews:NBK1', 'MedGen:C9999'])
        self.assertEqual(r['by_therapeutic_area']['EFO:0000200']['calls'], ['call-2'])                      # the calls of the cross-referenced trait, as trait_view listed them

    def test_a_table_without_the_columns_still_rolls_up_what_it_holds(self):
        d = tempfile.mkdtemp()
        try:
            pq.write_table(pa.Table.from_pylist([dict(id='EFO_0000001', name='one', ancestors=[], therapeuticAreas=['EFO_0000001'])],
                                                schema=pa.schema([('id', pa.string()), ('name', pa.string()), ('ancestors', pa.list_(pa.string())), ('therapeuticAreas', pa.list_(pa.string()))])),
                           os.path.join(d, 'p.parquet'))
            onto = views.ontology_terms(d)
            self.assertEqual((onto['EFO:0000001']['xrefs'], onto['EFO:0000001']['replaces']), ([], []))
            r = views.trait_rollup(traits('EFO:0000001', 'MedGen:C0001'), onto)
            self.assertEqual((r['cross_referenced'], r['not_rolled_up']), ({}, ['MedGen:C0001']))
        finally:
            shutil.rmtree(d, ignore_errors=True)


if __name__ == '__main__':
    unittest.main()
