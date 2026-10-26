import json
import os
import shutil
import sys
import tempfile
import unittest

from helpers import syn_reference, FIXTURES
from vindex.build import runner
from vindex.build.gwas import GwasCatalogAdapter, parse_effect, efo_ids
from vindex.index_reader import Index

sys.path.insert(0, FIXTURES)
import make_sources  # noqa: E402

SOURCES = make_sources.write_all(tempfile.mkdtemp())
CASES = make_sources.write_all.gwas_cases
COMP = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}


class StubVariants:
    """Alleles known at each site (what the Open Targets variant index would supply)."""

    def __init__(self, known, freqs=None):
        self.known, self.freqs = known, freqs or {}

    def load(self, positions):
        return len(positions)

    def snv_alts(self, chrom, pos, ref=None):
        return [(a, rs, f'{chrom}_{pos}_{ref}_{a}') for a, rs in self.known.get((chrom, pos), ())]

    def load_frequencies(self, ids):
        return {i: self.freqs[i] for i in ids if i in self.freqs}


def variants():
    c = CASES
    known = {('1', c['plus_alt']['pos']): [(c['plus_alt']['alt'], ('rs101',))],
             ('1', c['ref_effect']['pos']): [(c['ref_effect']['alt'], ('rs102',))],
             ('1', c['minus_strand']['pos']): [(c['minus_strand']['alt'], ('rs103',))],
             ('1', c['palindromic_resolved']['pos']): [(c['palindromic_resolved']['alt'], ('rs104',))],
             ('1', c['palindromic_unresolved']['pos']): [(c['palindromic_unresolved']['alt'], ('rs105',))],
             ('1', c['unknown_allele']['pos']): [(c['unknown_allele']['alt'], ('rs106',))],
             ('1', c['chr_verified']['pos']): [(c['chr_verified']['alt'], ())],
             ('1', c['merged_rsid_p_zero']['pos']): [(c['merged_rsid_p_zero']['alt'], ('rs110',))]}
    pr = c['palindromic_resolved']
    freqs = {f"1_{pr['pos']}_{pr['ref']}_{pr['alt']}": [0.2, 0.25]}
    return StubVariants(known, freqs)


class TestGwas(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.ref = syn_reference()
        cls.out = tempfile.mkdtemp()

        class Genes:
            rows = []
        cls.man = runner.build(cls.out, [GwasCatalogAdapter(Genes(), variants(), root=SOURCES)], cls.ref)
        cls.idx = Index.open(cls.out)
        cls.ev = cls.idx.table('evidence').to_pylist()

    @classmethod
    def tearDownClass(cls):
        shutil.rmtree(cls.out, ignore_errors=True)

    def assoc(self, **kw):
        return [r for r in self.ev if r['claim_type'] == 'gwas_association' and all(r.get(k) == v for k, v in kw.items())]

    def by_pos(self, pos):
        return [r for r in self.ev if r['pos'] == pos]

    def test_every_record_is_kept_once_as_an_association_row(self):
        st = self.man['adapters'][0]['stats']
        self.assertEqual(st['records'], 17)
        self.assertEqual(len(self.assoc()), 17)
        reasons = {r['mapping_reason'] for r in self.assoc() if r['mapping_status'] != 'mapped'}
        self.assertEqual(reasons, {'strand_unresolved_palindromic', 'interaction_between_variants', 'multi_variant_record', 'non_rsid_identifier_or_hla_or_amino_acid', 'author_position_assembly_unverified',
                                   'indel_notation_in_author_identifier', 'multi_base_risk_allele', 'no_position_in_catalog', 'no_alternate_allele_known', 'position_beyond_contig_length'})

    def test_an_impossible_position_is_kept_as_unmapped_not_coerced(self):
        r = self.assoc(unmapped_key='chr1:99999999999')
        self.assertEqual(len(r), 1)
        self.assertEqual((r[0]['mapping_status'], r[0]['mapping_reason'], r[0]['subject_id'], r[0]['chrom'], r[0]['pos']), ('unmapped', 'position_beyond_contig_length', None, None, None))
        self.assertEqual(self.man['adapters'][0]['stats']['positions_beyond_contig_length'], 1)
        self.assertTrue(any(x['claim_type'] == 'gwas_effect_size' and x['unmapped_key'] == 'chr1:99999999999' for x in self.ev))        # its effect size row is kept as well

    def test_plus_minus_and_reference_effect_alleles(self):
        c = CASES
        a = self.assoc(subject_id=f"1:{c['plus_alt']['pos']}:{c['plus_alt']['ref']}:{c['plus_alt']['alt']}")[0]
        self.assertEqual((a['effect_allele'], a['other_allele'], a['p_value'], a['effect_assertion'], a['sample_size'], a['n_cases'], a['ancestry']), (c['plus_alt']['alt'], c['plus_alt']['ref'], 3e-9, 'effect', 1500, 400, 'East Asian; European'))
        m = self.assoc(subject_id=f"1:{c['minus_strand']['pos']}:{c['minus_strand']['ref']}:{c['minus_strand']['alt']}")[0]
        self.assertEqual(m['effect_allele'], c['minus_strand']['alt'])
        self.assertIn('allele_reported_on_minus_strand', m['extra_json'])
        r = self.assoc(subject_id=f"1:{c['ref_effect']['pos']}:{c['ref_effect']['ref']}:{c['ref_effect']['alt']}")[0]
        self.assertEqual((r['effect_allele'], r['other_allele'], r['effect_assertion']), (c['ref_effect']['ref'], c['ref_effect']['alt'], 'uncertain'))    # p 2e-6 is suggestive only

    def test_effect_size_rows_keep_units_signs_and_direction(self):
        c = CASES
        eff = {r['subject_id']: r for r in self.ev if r['claim_type'] == 'gwas_effect_size' and r['subject_id']}
        plus = eff[f"1:{c['plus_alt']['pos']}:{c['plus_alt']['ref']}:{c['plus_alt']['alt']}"]
        self.assertEqual((plus['odds_ratio'], plus['ci_low'], plus['ci_high'], plus['direction'], plus['axis']), (1.25, 1.1, 1.4, 'increases', 'effect_magnitude_direction'))
        ref = eff[f"1:{c['ref_effect']['pos']}:{c['ref_effect']['ref']}:{c['ref_effect']['alt']}"]
        self.assertEqual((ref['beta'], ref['direction'], ref['score_unit']), (-0.04, 'decreases', 'unit decrease'))      # a decrease is stored with a negative sign
        self.assertEqual(plus['lineage_group'], 'gwas_study:GCST000001')
        pair = [r for r in self.by_pos(c['plus_alt']['pos']) if r['claim_type'] in ('gwas_association', 'gwas_effect_size')]
        self.assertEqual(len({r['lineage_group'] for r in pair}), 1)                                                     # one measurement, not two independent ones

    def test_unknown_effect_allele_maps_to_the_site_without_direction(self):
        c = CASES['unknown_allele']
        rows = self.by_pos(c['pos'])
        self.assertTrue(rows and all(r['subject_id'] == f"1:{c['pos']}:{c['ref']}:{c['alt']}" and r['effect_allele'] is None for r in rows))
        e = [r for r in rows if r['claim_type'] == 'gwas_effect_size'][0]
        self.assertIsNone(e['direction'])
        self.assertIn('effect_allele_unknown', rows[0]['extra_json'])

    def test_palindromic_sites_need_frequency(self):
        ok = CASES['palindromic_resolved']
        r = self.assoc(subject_id=f"1:{ok['pos']}:{ok['ref']}:{ok['alt']}")[0]
        self.assertEqual(r['effect_allele'], ok['ref'])
        bad = self.by_pos(CASES['palindromic_unresolved']['pos'])
        self.assertTrue(bad and all(x['mapping_status'] == 'rejected' and x['subject_id'] is None for x in bad))

    def test_author_positions_are_verified_and_merged_rsids_and_zero_p_values_handled(self):
        c = CASES['chr_verified']
        self.assertEqual(self.assoc(subject_id=f"1:{c['pos']}:{c['ref']}:{c['alt']}")[0]['mapping_status'], 'mapped')
        m = CASES['merged_rsid_p_zero']
        row = self.assoc(subject_id=f"1:{m['pos']}:{m['ref']}:{m['alt']}")[0]
        self.assertEqual(row['p_value'], 0.0)
        self.assertIn('p_reported_as_zero', row['extra_json'])
        xr = [x for x in self.idx.table('xref').to_pylist() if x['identifier'] == 'rs110']
        self.assertEqual(len(xr), 1)                                    # the current rsID of the merged record, not the retired one

    def test_studies_become_entities(self):
        st = self.idx.table('entity_study').to_pylist()
        self.assertEqual([(s['study_id'], s['n_total'], s['n_cases'], s['has_sumstats'], s['lineage_group']) for s in st], [('GCST000001', 1500, 400, True, 'gwas_study:GCST000001')])

    def test_effect_parsing(self):
        self.assertEqual(parse_effect('1.25', '[1.1-1.4]')[:5], ('odds_ratio', 1.25, 1.1, 1.4, None))
        self.assertEqual(parse_effect('0.5', '[0.4-0.6] mmol/L decrease')[:2], ('beta', -0.5))
        self.assertEqual(parse_effect('13.6', 'NR unit decrease')[:2], ('beta', -13.6))
        self.assertEqual(parse_effect('', '')[0], None)
        self.assertEqual(efo_ids('http://www.ebi.ac.uk/efo/EFO_0007660, http://purl.obolibrary.org/obo/MONDO_0004989'), ['EFO:0007660', 'MONDO:0004989'])


if __name__ == '__main__':
    unittest.main()
