import unittest

from helpers import ROOT  # noqa: F401  (puts the package on the path)
from vindex.summary import summarise


def row(**kw):
    base = dict(axis='clinical_validity_actionability', effect_assertion='effect', evidence_kind='curated_expert_classification', lineage_group=None, source_id='clinvar', accession='A',
                evidence_id='e', review_stars=None, p_value=None, pip=None, condition_ids=None)
    base.update(kw)
    return base


class TestSummary(unittest.TestCase):
    def test_replication_is_counted_from_independent_lineages_asserting_the_effect(self):
        rows = [row(evidence_id='1', lineage_group='L1', review_stars=1), row(evidence_id='2', lineage_group='L1', review_stars=3),          # one lineage, two rows
                row(evidence_id='3', lineage_group='L2', effect_assertion='no_effect'), row(evidence_id='4', lineage_group='L3', evidence_kind='aggregated_integration')]
        s = summarise(rows)['clinical_validity_actionability']
        self.assertEqual((s['n_rows'], s['n_independent']), (4, 3))
        self.assertEqual(s['n_independent_by_assertion'], {'effect': 1, 'no_effect': 1})          # the aggregate restates other rows and is not an independent source
        self.assertEqual((s['max_review_stars'], s['conflict']), (3, True))

    def test_descriptors_stay_per_axis_and_missing_values_stay_missing(self):
        rows = [row(axis='association_replication', source_id='gwas_catalog', p_value=3e-9, lineage_group='g1', condition_ids=['EFO_1']),
                row(axis='association_replication', source_id='gwas_catalog', p_value=2e-6, lineage_group='g2', effect_assertion='uncertain', condition_ids=['EFO_2']),
                row(axis='fine_mapping_probability', source_id='opentargets', pip=0.8, effect_assertion='not_applicable')]
        s = summarise(rows)
        self.assertEqual((s['association_replication']['min_p_value'], s['association_replication']['max_pip'], s['association_replication']['max_review_stars']), (3e-9, None, None))
        self.assertEqual(s['association_replication']['condition_ids_effect'], ['EFO_1'])                      # only rows that assert an effect name a condition here
        self.assertEqual((s['fine_mapping_probability']['max_pip'], s['fine_mapping_probability']['min_p_value']), (0.8, None))
        self.assertNotIn('clinical_validity_actionability', s)                                              # an axis without rows is absent, not zero

    def test_empty(self):
        self.assertEqual(summarise([]), {})


if __name__ == '__main__':
    unittest.main()
