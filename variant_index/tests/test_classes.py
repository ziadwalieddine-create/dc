import unittest

from helpers import ROOT  # noqa: F401  (puts the package on the path)
from vindex.build.classes import effect_of, stars, unknown_tokens


class TestClassificationVocabulary(unittest.TestCase):
    CASES = {
        # strings as the ClinVar VCF writes them
        'Pathogenic': 'effect', 'Likely_pathogenic': 'effect', 'Pathogenic/Likely_pathogenic': 'effect', 'Benign/Likely_benign': 'no_effect', 'Uncertain_significance': 'uncertain',
        'Conflicting_classifications_of_pathogenicity': 'conflicting', 'drug_response': 'effect', 'risk_factor': 'effect', 'association': 'effect', 'Affects': 'effect', 'protective': 'effect',
        'Pathogenic,_low_penetrance': 'effect', 'Likely_pathogenic,_low_penetrance': 'effect', 'Pathogenic/Likely_pathogenic/Pathogenic,_low_penetrance': 'effect',
        'Uncertain_significance/Uncertain_risk_allele': 'uncertain', 'Pathogenic/Likely_risk_allele': 'effect', 'Uncertain_significance/VUS-high': 'uncertain', 'VUS-low': 'uncertain',
        'Pathogenic|other': 'effect', 'Conflicting_classifications_of_pathogenicity|other': 'conflicting', 'Benign|risk_factor': 'conflicting', 'Uncertain_significance|risk_factor': 'effect',
        'association_not_found': 'no_effect', 'Tier_IV_-_Benign/Likely_benign': 'no_effect', 'Tier_I_-_Strong': 'effect', 'Tier_III_-_Unknown': 'uncertain', 'Likely_oncogenic': 'effect',
        # assertions of no assertion
        'not_provided': 'not_applicable', 'other': 'not_applicable', 'no_classification_for_the_single_variant': 'not_applicable', 'no_classifications_from_unflagged_records': 'not_applicable', '': 'not_applicable',
        # labels as submitters spell them in the VCV release
        'Likely Pathogenic': 'effect', 'likely pathogenic': 'effect', 'PATHOGENIC': 'effect', 'Uncertain Significance': 'uncertain', 'Likely Benign': 'no_effect', 'Pathogenic, low penetrance': 'effect',
        'conflicting data from submitters': 'conflicting', 'Tier I - Strong': 'effect', 'Established risk allele': 'effect', 'no known pathogenicity': 'no_effect', 'Uncertain risk allele': 'uncertain',
    }

    def test_known_labels(self):
        for label, want in self.CASES.items():
            with self.subTest(label=label):
                self.assertEqual(effect_of(label), want)
                self.assertEqual(unknown_tokens(label), [])

    def test_an_unknown_label_is_uncertain_and_reported_not_hidden(self):
        self.assertEqual(effect_of('Totally new label'), 'uncertain')
        self.assertEqual(unknown_tokens('Totally new label'), ['totally_new_label'])
        self.assertEqual(effect_of('Pathogenic|Totally new label'), 'effect')
        self.assertEqual(unknown_tokens('Pathogenic|Totally new label'), ['totally_new_label'])

    def test_review_status_stars_read_both_spellings(self):
        self.assertEqual(stars('criteria_provided,_single_submitter'), 1)
        self.assertEqual(stars('criteria provided, single submitter'), 1)
        self.assertEqual(stars('reviewed by expert panel'), 3)
        self.assertEqual(stars('practice guideline'), 4)
        self.assertEqual(stars('flagged submission'), 0)
        self.assertIsNone(stars('a status that does not exist'))
        self.assertIsNone(stars(None))


if __name__ == '__main__':
    unittest.main()
