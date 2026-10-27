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
             ('1', c['merged_rsid_p_zero']['pos']): [(c['merged_rsid_p_zero']['alt'], ('rs110',))],
             ('1', c['rsid_of_another_variant']['pos']): [(c['rsid_of_another_variant']['alt'], ('rs999',))]}
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
        self.assertEqual(st['records'], 18)
        self.assertEqual(len(self.assoc()), 18)
        reasons = {r['mapping_reason'] for r in self.assoc() if r['mapping_status'] != 'mapped'}
        self.assertEqual(reasons, {'strand_unresolved_palindromic', 'interaction_between_variants', 'multi_variant_record', 'non_rsid_identifier_or_hla_or_amino_acid', 'author_position_assembly_unverified',
                                   'indel_notation_in_author_identifier', 'multi_base_risk_allele', 'no_position_in_catalog', 'no_alternate_allele_known', 'position_beyond_contig_length'})

    def test_an_rsid_is_never_given_the_alternate_allele_of_another_variant_at_its_position(self):
        # found by the known-locus check of CFTR F508del: the Catalog's rs113993960 (a three-base deletion, given at one position) was mapped to a SNV T>G that Open Targets lists under
        # another rsID at that position, and the xref table then said rs113993960 was that SNV; 650 rows of the real file were mapped this way
        c = CASES['rsid_of_another_variant']
        r = self.assoc(unmapped_key='rs112')
        self.assertEqual(len(r), 1)
        self.assertEqual((r[0]['mapping_status'], r[0]['mapping_reason'], r[0]['subject_id']), ('unmapped', 'no_alternate_allele_known', None))
        self.assertEqual((r[0]['chrom'], r[0]['pos']), ('1', c['pos']))                                                                # still findable by position
        self.assertEqual([x for x in self.idx.table('xref').to_pylist() if x['identifier'] == 'rs112'], [])
        self.assertEqual([x for x in self.ev if x['subject_id'] == f"1:{c['pos']}:{c['ref']}:{c['alt']}"], [])                         # nothing is attached to the other variant
        self.assertFalse([x for x in self.ev if 'rsid_not_listed_at_position' in (x['extra_json'] or '')])

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
        self.assertEqual((ref['beta'], ref['direction'], ref['score_unit']), (-0.04, 'decreases', 'unit'))               # a decrease is stored with a negative sign; the unit is not 'unit decrease'
        self.assertEqual((ref['ci_low'], ref['ci_high']), (-0.06, -0.02))                                              # ... and so is its interval ([0.02-0.06] in the file)
        self.assertEqual(plus['lineage_group'], 'pmid:1')                                                              # PMID 1 of the fixture's association rows
        pair = [r for r in self.by_pos(c['plus_alt']['pos']) if r['claim_type'] in ('gwas_association', 'gwas_effect_size')]
        self.assertEqual(len({r['lineage_group'] for r in pair}), 1)                                                     # one measurement, not two independent ones

    def test_a_paper_has_one_lineage_whichever_database_describes_it(self):
        # real data: 4,328 of the 4,385 publications behind Open Targets' GWAS credible sets are also in the Catalog; under the Catalog's study accession their rows counted as independent
        from vindex.build.gwas import lineage_of
        from vindex.build.opentargets_genetics import OpenTargetsGeneticsAdapter
        ot = OpenTargetsGeneticsAdapter(None).lineage(dict(studyType='gwas', pubmedId='33139745', projectId='GCST', studyId='GCST90012345'))
        self.assertEqual(lineage_of('33139745', 'GCST90012345'), ot)
        self.assertEqual(lineage_of('33139745', 'GCST90012346'), ot)                              # a second Catalog study of the same paper: the same lineage (never more independent)
        self.assertEqual(lineage_of('', 'GCST9'), 'gwas_study:GCST9')                              # no PubMed id: the accession is all there is
        self.assertEqual(lineage_of('NR', 'GCST9'), 'gwas_study:GCST9')

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
        self.assertEqual([(s['study_id'], s['n_total'], s['n_cases'], s['has_sumstats'], s['lineage_group']) for s in st], [('GCST000001', 1500, 400, True, 'pmid:1')])

    def test_effect_parsing(self):
        self.assertEqual(parse_effect('1.25', '[1.1-1.4]')[:5], ('odds_ratio', 1.25, 1.1, 1.4, None))
        self.assertEqual(parse_effect('0.5', '[0.4-0.6] mmol/L decrease')[:2], ('beta', -0.5))
        self.assertEqual(parse_effect('13.6', 'NR unit decrease')[:2], ('beta', -13.6))
        self.assertEqual(parse_effect('', '')[0], None)
        self.assertIsNone(parse_effect('nan', '[1.1-1.4]')[0])                                                    # float() reads 'nan' and 'inf'; the real file has such cells
        self.assertIsNone(parse_effect('inf', '')[0])

    def test_the_interval_of_a_beta_contains_the_signed_beta(self):
        # real file: for 393,000 decreases the interval was left as the bounds of the magnitude, so a beta of -0.032 sat beside [0.018, 0.046]
        v = parse_effect('0.032', '[0.018-0.046] unit decrease')
        self.assertEqual((v[0], v[1], v[2], v[3]), ('beta', -0.032, -0.046, -0.018))
        v = parse_effect('0.0687', '[-0.05--0.087] unit decrease')                                                     # bounds already signed, in the order they were typed
        self.assertEqual((v[1], v[2], v[3]), (-0.0687, -0.087, -0.05))
        v = parse_effect('0.02', '[0.01-0.03] unit increase')
        self.assertEqual((v[1], v[2], v[3]), (0.02, 0.01, 0.03))
        v = parse_effect('0.02', '[-0.01-0.05] unit increase')                                                         # an interval that includes zero stays as it is
        self.assertEqual((v[1], v[2], v[3]), (0.02, -0.01, 0.05))
        v = parse_effect('1.5', '[1.8-1.2]')                                                                          # an odds ratio interval typed in the wrong order is sorted
        self.assertEqual((v[0], v[2], v[3]), ('odds_ratio', 1.2, 1.8))
        for value, text in (('0.032', '[0.018-0.046] unit decrease'), ('0.0687', '[-0.05--0.087] unit decrease'), ('0.02', '[0.01-0.03] unit increase')):
            kind, b, lo, hi = parse_effect(value, text)[:4]
            self.assertTrue(lo <= b <= hi, (value, text, lo, b, hi))

    def test_the_row_says_why_it_is_a_beta_or_an_odds_ratio(self):
        from vindex.build.gwas import effect_type_basis
        self.assertEqual(effect_type_basis('beta', '[0.1-0.2] unit increase'), 'unit_in_interval_text')
        self.assertEqual(effect_type_basis('odds_ratio', '[1.1-1.4]'), 'interval_text_without_unit')
        self.assertEqual(effect_type_basis('odds_ratio', ''), 'assumed_odds_ratio_no_interval_text')
        eff = [r for r in self.ev if r['claim_type'] == 'gwas_effect_size' and r['subject_id']]
        self.assertTrue(all(json.loads(r['extra_json'])['effect_type_basis'] for r in eff))

    def test_a_frequency_outside_zero_to_one_is_not_a_frequency(self):
        from vindex.build.gwas import frequency
        self.assertEqual((frequency('0.25'), frequency('0'), frequency('1')), (0.25, 0.0, 1.0))
        self.assertEqual((frequency('36'), frequency('-0.1'), frequency('NR'), frequency(''), frequency('nan')), (None,) * 5)

    def test_not_reported_is_the_absence_of_a_unit_or_a_category_not_a_value(self):
        # real file: 2,097 odds ratios carried the unit 'NR' and 34,838 rows the ancestry 'NR'
        self.assertEqual(parse_effect('1.2', 'NR'), ('odds_ratio', 1.2, None, None, None, 'increases'))
        self.assertIsNone(parse_effect('1.2', '[NR]')[4])
        self.assertEqual(parse_effect('13.6', 'NR unit decrease')[4], 'unit')                                # the unit after the missing interval is kept (without its direction word)
        self.assertEqual(parse_effect('13.6', '[NR] mmol/L decrease')[4], 'mmol/L')
        self.assertEqual(parse_effect('1.1', 'NRF2 activity increase')[4], 'NRF2 activity')                  # only the token NR is the marker
        self.assertIsNone(parse_effect('1.2', '[1.1-1.4]')[4])

    def test_the_direction_word_is_the_direction_column_not_part_of_the_unit(self):
        # real file: 1,625 distinct units before the words were taken out ('unit increase' on 426,732 rows and 'unit decrease' on 411,273)
        for ci, unit, direction in (('[0.02-0.043] kg/m2 increase', 'kg/m2', 'increases'), ('[0.34 - 0.46] kg/m2 increase (adults)', 'kg/m2 (adults)', 'increases'),
                                    ('[0.1-0.2] mmol/L decrease', 'mmol/L', 'decreases'), ('unit increase', 'unit', 'increases'), ('SD unit decrease', 'SD unit', 'decreases'),
                                    ('[0.09-0.11] standard deviation unit per year increase', 'standard deviation unit per year', 'increases'), ('increase', None, 'increases')):
            v = parse_effect('0.1', ci)
            self.assertEqual((v[0], v[4], v[5]), ('beta', unit, direction), ci)

    def test_intervals_in_the_other_notations_of_the_catalogue_are_read_as_intervals(self):
        # real file: 1,474 rows had an interval typed without brackets, with a comma, an en dash, a unicode minus or a closing brace, and kept it as their 'unit'
        for ci, lo, hi in (('[0.980, 0.990]', 0.98, 0.99), ('1.04-1.08', 1.04, 1.08), ('[1.48\u20131.75]', 1.48, 1.75), ('[1.17023-1.33871}', 1.17023, 1.33871), ('(1.2-1.9)', 1.2, 1.9),
                           ('[0.18\u22120.22]', 0.18, 0.22), ('[12.5-\u221e]', 12.5, None)):
            v = parse_effect('1.5', ci)
            self.assertEqual((v[0], v[2], v[3], v[4]), ('odds_ratio', lo, hi, None), ci)
        v = parse_effect('0.3', '0.18-0.45 unit increase')
        self.assertEqual((v[0], v[2], v[3], v[4]), ('beta', 0.18, 0.45, 'unit'))
        v = parse_effect('0.3', '[-0.5--0.1] unit increase')                                                       # the strict notation still reads a negative range
        self.assertEqual((v[2], v[3]), (-0.5, -0.1))
        self.assertIsNone(parse_effect('1.5', 'kg/m2')[2])                                                              # text that merely starts with a letter is a unit, not an interval
        for ci, lo, hi in (('[1.18\u20101.34]', 1.18, 1.34), ('[1\u00b700\u20131\u00b723]', 1.0, 1.23), ('[2\u00b731-4\u00b755]', 2.31, 4.55)):          # a hyphen of another kind, a middle dot as the decimal point
            v = parse_effect('1.5', ci)
            self.assertEqual((v[0], v[2], v[3], v[4]), ('odds_ratio', lo, hi, None), ci)
        self.assertIsNone(parse_effect('13.6', 'NR NR decrease')[4])                                                   # the unit that is itself 'NR' is no unit
        self.assertEqual(parse_effect('13.6', '[NR] 1,000/\u03bcL decrease')[4], '1,000/\u03bcL')

    def test_a_blank_trait_and_a_cohort_of_nr_are_missing_values_not_text(self):
        # real file: 2,395 studies had the cohort 'NR' and a few associations a blank DISEASE/TRAIT, kept as '' in condition_text
        from vindex.build.gwas import text_or_none
        self.assertEqual([text_or_none(x) for x in ('', ' ', 'NR', 'nr', 'NA', '-', None, 'UK Biobank', ' Type 2 diabetes ')], [None] * 7 + ['UK Biobank', 'Type 2 diabetes'])
        d = tempfile.mkdtemp()
        try:
            with open(os.path.join(d, 'gwas-catalog-download-studies-v1.0.3.1.tsv'), 'w') as f:
                f.write('STUDY ACCESSION\tPUBMED ID\tDISEASE/TRAIT\tMAPPED_TRAIT\tMAPPED_TRAIT_URI\tCOHORT\tFULL SUMMARY STATISTICS\tINITIAL SAMPLE SIZE\n'
                        'GCST9\t1\t\tt\t\tNR\tno\tx\nGCST8\t2\tType 2 diabetes\tt\t\tUKB\tno\ty\n')
            with open(os.path.join(d, 'gwas-catalog-download-ancestries-v1.0.3.1.tsv'), 'w') as f:
                f.write('STUDY ACCESSION\tBROAD ANCESTRAL CATEGORY\tSTAGE\tNUMBER OF INDIVIDUALS\tNUMBER OF CASES\n')
            rows = []

            class Ctx:
                def emit(self, table, row, name=None):
                    rows.append(row)

            class St:
                def inc(self, k, n=1):
                    pass
            GwasCatalogAdapter(None, None, root=d).studies(Ctx(), 'rel_s', 'rel_a', St())
        finally:
            shutil.rmtree(d, ignore_errors=True)
        by = {r['study_id']: r for r in rows}
        self.assertEqual((by['GCST9']['cohorts'], by['GCST9']['trait_text']), (None, None))
        self.assertEqual((by['GCST8']['cohorts'], by['GCST8']['trait_text']), ('UKB', 'Type 2 diabetes'))
        blank = self.assoc(unmapped_key='chr1:99999999999')[0]                                                         # the fixture record whose DISEASE/TRAIT and MAPPED_TRAIT are blank
        self.assertEqual((blank['condition_text'], blank['claim_value']), (None, None))
        self.assertEqual({r['condition_text'] for r in self.ev if r['unmapped_key'] == 'chr1:99999999999'}, {None})
        self.assertFalse([r for r in self.ev if r['condition_text'] == ''])

    def test_impossible_cells_are_not_taken_as_measurements(self):
        import zipfile
        c = CASES
        d = tempfile.mkdtemp()
        try:
            p1, p2, p3 = c['plus_alt'], c['ref_effect'], c['minus_strand']
            rows = [make_sources.gwas_row(**{'CHR_ID': '1', 'CHR_POS': str(p1['pos']), 'STRONGEST SNP-RISK ALLELE': f"rs101-{p1['alt']}", 'SNPS': 'rs101', 'RISK ALLELE FREQUENCY': '36', 'OR or BETA': '1.25', '95% CI (TEXT)': '[1.1-1.4]'}),
                    make_sources.gwas_row(**{'CHR_ID': '1', 'CHR_POS': str(p2['pos']), 'STRONGEST SNP-RISK ALLELE': f"rs102-{p2['ref']}", 'SNPS': 'rs102', 'OR or BETA': 'nan', '95% CI (TEXT)': '[0.02-0.06] unit decrease'}),
                    make_sources.gwas_row(**{'CHR_ID': '1', 'CHR_POS': str(p3['pos']), 'STRONGEST SNP-RISK ALLELE': f"rs103-{COMP[p3['alt']]}", 'SNPS': 'rs103', 'RISK ALLELE FREQUENCY': '0.3', 'OR or BETA': '0.04', '95% CI (TEXT)': '[0.02-0.06] unit decrease'})]
            zp = os.path.join(d, 'assoc.zip')
            with zipfile.ZipFile(zp, 'w') as z:
                z.writestr('gwas-catalog-associations.tsv', '\t'.join(make_sources.GWAS_COLS) + '\n' + '\n'.join(rows) + '\n')

            class Genes:
                rows = []
            out = os.path.join(d, 'idx')
            man = runner.build(out, [GwasCatalogAdapter(Genes(), variants(), root=SOURCES, assoc_zip=zp)], self.ref)
            ev = Index.open(out).table('evidence').to_pylist()
        finally:
            shutil.rmtree(d, ignore_errors=True)
        a1 = [r for r in ev if r['claim_type'] == 'gwas_association' and r['pos'] == c['plus_alt']['pos']][0]
        self.assertEqual((a1['mapping_status'], a1['eaf']), ('mapped', None))                        # 36 is not a frequency: kept out of eaf, and not used to settle a strand
        self.assertEqual(json.loads(a1['extra_json'])['risk_allele_frequency_reported'], '36')
        self.assertEqual(man['adapters'][0]['stats']['risk_allele_frequency_outside_0_1'], 1)
        self.assertEqual([r for r in ev if r['claim_type'] == 'gwas_effect_size' and r['pos'] == c['ref_effect']['pos']], [])      # 'nan' is no effect size; the association row stays
        self.assertEqual(len([r for r in ev if r['claim_type'] == 'gwas_association' and r['pos'] == c['ref_effect']['pos']]), 1)
        e3 = [r for r in ev if r['claim_type'] == 'gwas_effect_size' and r['pos'] == c['minus_strand']['pos']][0]
        self.assertEqual((e3['beta'], e3['ci_low'], e3['ci_high']), (-0.04, -0.06, -0.02))

    def test_a_study_whose_ancestry_is_not_reported_has_no_ancestry_but_keeps_its_sample_size(self):
        d = tempfile.mkdtemp()
        try:
            with open(os.path.join(d, 'gwas-catalog-download-studies-v1.0.3.1.tsv'), 'w') as f:
                f.write('STUDY ACCESSION\tPUBMED ID\tDISEASE/TRAIT\tMAPPED_TRAIT\tMAPPED_TRAIT_URI\tCOHORT\tFULL SUMMARY STATISTICS\tINITIAL SAMPLE SIZE\n'
                        'GCST9\t1\tt\tt\t\t\tno\tx\nGCST8\t2\tt\tt\t\t\tno\ty\nGCST7\t3\tt\tt\t\t\tno\tz\n')
            with open(os.path.join(d, 'gwas-catalog-download-ancestries-v1.0.3.1.tsv'), 'w') as f:
                f.write('STUDY ACCESSION\tBROAD ANCESTRAL CATEGORY\tSTAGE\tNUMBER OF INDIVIDUALS\tNUMBER OF CASES\n'
                        'GCST9\tNR\tinitial\t1000\t\nGCST8\tNR\tinitial\t200\t\nGCST8\tEuropean\tinitial\t300\t\nGCST7\tAfrican\treplication\t50\t\n')

            class Ctx:
                def __init__(self):
                    self.rows = []

                def emit(self, table, row, name=None):
                    self.rows.append(row)

            class St:
                def __init__(self):
                    self.c = {}

                def inc(self, k, n=1):
                    self.c[k] = self.c.get(k, 0) + n
            ctx, st = Ctx(), St()
            GwasCatalogAdapter(None, None, root=d).studies(ctx, 'rel_s', 'rel_a', st)
        finally:
            shutil.rmtree(d, ignore_errors=True)
        by = {r['study_id']: r for r in ctx.rows}
        self.assertEqual((by['GCST9']['ancestry'], by['GCST9']['n_total']), (None, 1000))
        self.assertEqual((by['GCST8']['ancestry'], by['GCST8']['n_total']), ('European', 500))               # a partly reported study keeps the categories it reported
        self.assertEqual((by['GCST7']['ancestry'], by['GCST7']['n_total']), (None, None))                    # replication-stage rows do not describe the initial sample
        self.assertEqual(st.c['ancestry_rows_with_category_not_reported'], 2)
        self.assertEqual(efo_ids('http://www.ebi.ac.uk/efo/EFO_0007660, http://purl.obolibrary.org/obo/MONDO_0004989'), ['EFO:0007660', 'MONDO:0004989'])


if __name__ == '__main__':
    unittest.main()
