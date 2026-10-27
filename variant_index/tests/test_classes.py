import unittest

from helpers import ROOT  # noqa: F401  (puts the package on the path)
from vindex.build.classes import effect_of, inheritance_label, review_label, stars, unknown_tokens


class TestClassificationVocabulary(unittest.TestCase):
    CASES = {
        # strings as the ClinVar VCF writes them
        'Pathogenic': 'effect', 'Likely_pathogenic': 'effect', 'Pathogenic/Likely_pathogenic': 'effect', 'Benign/Likely_benign': 'no_effect', 'Uncertain_significance': 'uncertain',
        'Conflicting_classifications_of_pathogenicity': 'conflicting', 'drug_response': 'effect', 'risk_factor': 'effect', 'association': 'effect', 'Affects': 'effect', 'protective': 'effect',
        'Pathogenic,_low_penetrance': 'effect', 'Likely_pathogenic,_low_penetrance': 'effect', 'Pathogenic/Likely_pathogenic/Pathogenic,_low_penetrance': 'effect',
        'Uncertain_significance/Uncertain_risk_allele': 'uncertain', 'Pathogenic/Likely_risk_allele': 'effect', 'Uncertain_significance/VUS-high': 'uncertain', 'VUS-low': 'uncertain',
        'Pathogenic|other': 'effect', 'Conflicting_classifications_of_pathogenicity|other': 'conflicting', 'Benign|risk_factor': 'conflicting', 'Uncertain_significance|risk_factor': 'effect',
        'association_not_found': 'no_effect', 'Tier_IV_-_Benign/Likely_benign': 'no_effect', 'Tier_I_-_Strong': 'effect', 'Tier_III_-_Unknown': 'uncertain', 'Likely_oncogenic': 'effect',
        # semicolon-joined lists as the 2026 VCV release writes them
        'Pathogenic;_drug_response': 'effect', 'Pathogenic;_risk_factor': 'effect', 'Pathogenic;_protective;_risk_factor': 'effect', 'Pathogenic;_affects': 'effect', 'Pathogenic;_other': 'effect',
        'Uncertain_significance;_pathogenic': 'effect', 'Benign;_risk_factor': 'conflicting', 'Pathogenic; drug response': 'effect',
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

    def test_one_review_status_is_one_value_whichever_file_wrote_it(self):
        self.assertEqual(review_label('criteria_provided,_single_submitter'), 'criteria provided, single submitter')                  # the VCF's spelling becomes the VCV tables'
        self.assertEqual(review_label('criteria provided, single submitter'), 'criteria provided, single submitter')
        self.assertEqual(review_label(None), '')

    def test_inheritance_codes_become_the_labels_the_other_rows_use_and_unknown_text_is_kept(self):
        self.assertEqual(inheritance_label('AD'), 'Autosomal dominant inheritance')
        self.assertEqual(inheritance_label('AD,AR'), 'Autosomal dominant inheritance; Autosomal recessive inheritance')
        self.assertEqual(inheritance_label('XR'), 'X-linked inheritance (recessive (HP:0001419))')
        self.assertEqual(inheritance_label('Autosomal recessive inheritance'), 'Autosomal recessive inheritance')            # already a label
        self.assertEqual(inheritance_label('AD,something new'), 'AD,something new')                                          # one unknown token: nothing is guessed
        self.assertEqual(inheritance_label('dominant'), 'dominant')                                                           # Open Targets' allelic requirements are not ClinGen codes
        self.assertIsNone(inheritance_label(''))
        self.assertIsNone(inheritance_label(None))

    def test_review_status_stars_read_both_spellings(self):
        self.assertEqual(stars('criteria_provided,_single_submitter'), 1)
        self.assertEqual(stars('criteria provided, single submitter'), 1)
        self.assertEqual(stars('reviewed by expert panel'), 3)
        self.assertEqual(stars('practice guideline'), 4)
        self.assertEqual(stars('flagged submission'), 0)
        self.assertEqual(stars('no_classifications_from_unflagged_records'), 0)                 # in the real file for 140 records
        self.assertEqual(stars('criteria_provided,_multiple_submitters'), 2)                    # in the real file for somatic clinical impact (9 records)
        self.assertIsNone(stars('a status that does not exist'))
        self.assertIsNone(stars(None))


if __name__ == '__main__':
    unittest.main()
