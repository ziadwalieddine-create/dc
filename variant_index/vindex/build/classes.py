"""ClinVar classification vocabulary, shared by the VCF adapter (aggregate strings such as 'Likely_pathogenic') and the VCV adapter (submitted
labels such as 'Likely Pathogenic', 'Pathogenic, low penetrance', 'Tier I - Strong').

effect_of() reduces a classification string to the five-valued effect_assertion of the index. The string itself is always kept in claim_value:
the reduction exists for filtering and contradiction checks, never as a replacement for the label.
  * a qualifier after a comma ('Pathogenic, low penetrance') does not change the class it qualifies
  * compound labels ('Benign/Likely benign', 'Uncertain significance/Uncertain risk allele', 'Pathogenic; drug response') and '|'-joined lists are read token by token
  * 'not provided', 'other' and the 'no classification' labels make no assertion about effect (not_applicable)
  * a label the table does not know is read as 'uncertain' and reported by unknown_tokens(), so vocabulary drift in a new release is counted, not hidden
"""
import re

STARS = {'practice_guideline': 4, 'reviewed_by_expert_panel': 3, 'criteria_provided,_multiple_submitters,_no_conflicts': 2, 'criteria_provided,_single_submitter': 1,
         'criteria_provided,_conflicting_classifications': 1, 'criteria_provided,_conflicting_interpretations': 1, 'no_assertion_criteria_provided': 0,
         'no_classification_provided': 0, 'no_classification_for_the_single_variant': 0, 'no_assertion_provided': 0, 'flagged_submission': 0,
         'no_classifications_from_unflagged_records': 0,                       # every submission flagged: no star
         'criteria_provided,_multiple_submitters': 2}                        # somatic clinical impact: conflicts do not arise, so the label has no ',_no_conflicts'

EFFECT = {
    'pathogenic': 'effect', 'likely_pathogenic': 'effect', 'benign': 'no_effect', 'likely_benign': 'no_effect', 'uncertain_significance': 'uncertain',
    'risk_factor': 'effect', 'association': 'effect', 'protective': 'effect', 'drug_response': 'effect', 'affects': 'effect', 'likely_risk_allele': 'effect',
    'established_risk_allele': 'effect', 'confers_sensitivity': 'effect', 'uncertain_risk_allele': 'uncertain', 'association_not_found': 'no_effect',
    'oncogenic': 'effect', 'likely_oncogenic': 'effect', 'tier_i_-_strong': 'effect', 'tier_i_-_potential': 'effect', 'tier_ii_-_strong': 'effect',
    'tier_ii_-_potential': 'effect', 'tier_iii_-_unknown': 'uncertain', 'tier_iv_-_benign': 'no_effect', 'tier_iv_-_likely_benign': 'no_effect',
    'vus-high': 'uncertain', 'vus-mid': 'uncertain', 'vus-low': 'uncertain', 'no_known_pathogenicity': 'no_effect',
    'not_provided': 'not_applicable', 'other': 'not_applicable', 'no_classification_provided': 'not_applicable', 'no_classification_for_the_single_variant': 'not_applicable',
    'no_classifications_from_unflagged_records': 'not_applicable',
    'conflicting_classifications_of_pathogenicity': 'conflicting', 'conflicting_classifications_of_oncogenicity': 'conflicting', 'conflicting_classification_of_pathogenicity': 'conflicting',
    'conflicting_interpretations_of_pathogenicity': 'conflicting', 'conflicting_data_from_submitters': 'conflicting',
}


def review_label(rev):
    """ClinVar review status in the spelling of the VCV tables (spaces). The VCF writes underscores; one status is one value in the index, whichever file it came from."""
    return (rev or '').replace('_', ' ')


# Mode-of-inheritance codes (ClinGen gene-disease validity, STRchive) as the HPO-style labels ClinGen's other tables use. Only a code that means one thing is expanded.
INHERITANCE_CODES = {'AD': 'Autosomal dominant inheritance', 'AR': 'Autosomal recessive inheritance', 'XL': 'X-linked inheritance', 'SD': 'Semidominant inheritance',
                     'MT': 'Mitochondrial inheritance', 'UD': 'Undetermined mode of inheritance', 'XD': 'X-linked inheritance (dominant (HP:0001423))',
                     'XR': 'X-linked inheritance (recessive (HP:0001419))'}


def inheritance_label(value):
    """The label for a code or comma-separated codes; a value that contains anything else is returned as the source wrote it."""
    if not value:
        return None
    codes = [c.strip() for c in value.split(',')]
    return '; '.join(INHERITANCE_CODES[c] for c in codes) if all(c in INHERITANCE_CODES for c in codes) else value


def norm_key(s):
    """Lower case with runs of whitespace as '_': the spelling the VCF uses, so VCF strings and VCV labels meet in one table."""
    return re.sub(r'\s+', '_', (s or '').strip().lower())


def stars(review_status):
    return STARS.get(norm_key(review_status))


def tokens(sig):
    s = re.sub(r',_?low_penetrance', '', norm_key(sig))
    parts = (re.sub(r'\(.*', '', p).strip('_') for p in re.split(r'[|/,;]', s))
    return [p for p in parts if p]


def unknown_tokens(sig):
    return [t for t in tokens(sig) if t not in EFFECT]


def effect_of(sig):
    toks = tokens(sig)
    if not toks:
        return 'not_applicable'
    vals = {EFFECT.get(t, 'uncertain') for t in toks}
    if vals == {'not_applicable'}:
        return 'not_applicable'
    vals.discard('not_applicable')
    if 'conflicting' in vals or {'effect', 'no_effect'} <= vals:
        return 'conflicting'
    if len(vals) == 1:
        return next(iter(vals))
    return 'effect' if 'effect' in vals else 'uncertain'
