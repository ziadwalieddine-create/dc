import json
import os
import shutil
import sys
import tempfile
import unittest

from helpers import syn_reference, FIXTURES
from vindex.build import runner
from vindex.build.genes import GeneMap
from vindex.build.regions import GiabRegionAdapter, StrchiveAdapter, giab_kind
from vindex.index_reader import Index
from vindex.match_repeat import loci_from_index, match_repeats
from vindex.match_small import is_difficult

sys.path.insert(0, FIXTURES)
import make_sources  # noqa: E402

SOURCES = make_sources.write_all(tempfile.mkdtemp())
ROOT = os.path.join(SOURCES, 'regions')


def build(out, adapters):
    ref = syn_reference()
    return runner.build(out, [a(GeneMap(ref, ot_dir='/nonexistent', mane_path='/nonexistent')) for a in adapters], ref)


class TestGiabRegions(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.out = tempfile.mkdtemp()
        cls.man = build(cls.out, [lambda g: GiabRegionAdapter(root=ROOT)])
        cls.idx = Index.open(cls.out)
        cls.rows = cls.idx.table('entity_region').to_pylist()
        cls.stats = cls.man['adapters'][0]['stats']

    @classmethod
    def tearDownClass(cls):
        shutil.rmtree(cls.out, ignore_errors=True)

    def test_kinds_separate_call_difficulty_from_sex_chromosome_context(self):
        self.assertEqual(giab_kind('segdups'), 'giab_difficult_segdups')
        self.assertEqual(giab_kind('chrx_par'), 'giab_xy_chrx_par')
        self.assertEqual({r['kind'] for r in self.rows}, {'giab_difficult_segdups', 'giab_difficult_alldifficultregions', 'giab_xy_chrx_par'})
        self.assertFalse(is_difficult([r for r in self.rows if r['kind'] == 'giab_xy_chrx_par']))
        self.assertTrue(is_difficult([r for r in self.rows if r['kind'] == 'giab_difficult_segdups']))

    def test_bed_intervals_become_one_based_inclusive_and_unusable_ones_are_counted(self):
        seg = sorted((r['chrom'], r['start'], r['end']) for r in self.rows if r['kind'] == 'giab_difficult_segdups')
        self.assertEqual(seg, [('1', 101, 300), ('1', 1001, 1200), ('2', 1, 50)])
        self.assertEqual(self.stats['intervals_on_other_contigs_or_unusable'], 2)            # an unplaced contig and an empty interval
        self.assertEqual(self.stats['intervals:giab_xy_chrx_par'], 1)
        self.assertEqual(next(r['region_id'] for r in self.rows if r['chrom'] == '2'), 'giab:segdups:2:1-50')

    def test_overlap_lookup_returns_components_next_to_the_union(self):
        kinds = lambda c, s, e: sorted(r['kind'] for r in self.idx.regions_overlapping(c, s, e))
        self.assertEqual(kinds('1', 150, 150), ['giab_difficult_alldifficultregions', 'giab_difficult_segdups'])
        self.assertEqual(kinds('1', 600, 600), [])
        self.assertEqual(kinds('X', 15, 15), ['giab_difficult_alldifficultregions', 'giab_xy_chrx_par'])

    def test_files_that_are_not_acquired_are_reported_not_ignored(self):
        self.assertGreater(self.stats['stratification_files_not_acquired'], 10)
        self.assertEqual(self.man['problems'], [])


class TestStrchive(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.out = tempfile.mkdtemp()
        cls.man = build(cls.out, [lambda g: StrchiveAdapter(g, root=os.path.join(ROOT, 'strchive'))])
        cls.idx = Index.open(cls.out)
        cls.stats = cls.man['adapters'][0]['stats']
        cls.ev = cls.idx.table('evidence').to_pylist()

    @classmethod
    def tearDownClass(cls):
        shutil.rmtree(cls.out, ignore_errors=True)

    def claims(self, locus):
        return {r['claim_type']: r for r in self.ev if r['accession'] == locus}

    def test_regions_and_accounting(self):
        regs = {r['region_id']: r for r in self.idx.table('entity_region').to_pylist()}
        self.assertEqual(set(regs), {'strchive:TEST1_GENEA', 'strchive:TEST2_GENEB', 'strchive:TEST3_MULTI', 'strchive:TEST5_NOJSON'})
        r = regs['strchive:TEST1_GENEA']
        self.assertEqual((r['kind'], r['chrom'], r['start'], r['end'], r['gene_ids']), ('strchive_locus', '1', 1001, 1030, ['symbol:GENEA']))
        self.assertEqual((self.stats['loci'], self.stats['loci_on_other_contigs'], self.stats['loci_without_stranger_record'], self.stats['stranger_records_without_general_row']), (5, 1, 1, 1))

    def test_thresholds_are_separate_claims_with_units(self):
        c = self.claims('TEST1_GENEA')
        self.assertEqual(set(c), {'strchive_disease_locus', 'repeat_normal_max', 'repeat_pathologic_min'})
        self.assertEqual((c['repeat_normal_max']['claim_num'], c['repeat_normal_max']['effect_assertion'], c['repeat_normal_max']['score_unit']), (35.0, 'no_effect', 'repeat_units_of_CAG'))
        self.assertEqual((c['repeat_pathologic_min']['claim_num'], c['repeat_pathologic_min']['effect_assertion']), (40.0, 'effect'))
        self.assertEqual((c['strchive_disease_locus']['condition_text'], c['strchive_disease_locus']['inheritance'], c['strchive_disease_locus']['level']), ('Test disease one', 'AD', 'region'))
        self.assertTrue(all(r['lineage_group'] == 'strchive:TEST1_GENEA' for r in c.values()))

    def test_a_missing_threshold_is_missing_not_zero(self):
        c = self.claims('TEST2_GENEB')
        self.assertEqual(set(c), {'strchive_disease_locus', 'repeat_normal_max'})
        self.assertEqual(self.stats['repeat_pathologic_min_missing'], 1)
        self.assertEqual(set(self.claims('TEST5_NOJSON')), {'strchive_disease_locus'})            # the BED's own minimum is kept as context only
        self.assertEqual(json.loads(self.claims('TEST5_NOJSON')['strchive_disease_locus']['extra_json'])['bed_pathogenic_min'], '35')

    def test_complex_locus_keeps_structure_and_the_region_the_threshold_applies_to(self):
        ex = json.loads(self.claims('TEST3_MULTI')['repeat_pathologic_min']['extra_json'])
        self.assertEqual((ex['pathologic_region'], ex['locus_structure']), ('chr2:500-530', '(AAAAT)*(TGAAA)*'))
        loc = json.loads(self.claims('TEST3_MULTI')['strchive_disease_locus']['extra_json'])
        self.assertEqual((loc['trgt_motifs'], loc['reference_region']), ('AAAAT,TGAAA', ['chr2:500-530', 'chr2:530-560']))

    def test_loci_round_trip_into_the_repeat_matcher(self):
        loci = {L['locus_id']: L for L in loci_from_index(self.idx)}
        self.assertEqual(set(loci), {'TEST1_GENEA', 'TEST2_GENEB', 'TEST3_MULTI', 'TEST5_NOJSON'})
        self.assertEqual((loci['TEST1_GENEA']['normal_max'], loci['TEST1_GENEA']['pathologic_min'], loci['TEST1_GENEA']['gene'], loci['TEST1_GENEA']['inheritance']), (35, 40, 'GENEA', 'AD'))
        self.assertEqual((loci['TEST2_GENEB']['normal_max'], loci['TEST2_GENEB']['pathologic_min']), (20, None))
        self.assertEqual((loci['TEST5_NOJSON']['normal_max'], loci['TEST5_NOJSON']['pathologic_min']), (None, None))

        class NoRepeats:
            repeats = []
        res = {r.locus_id: r for r in match_repeats(NoRepeats(), list(loci.values()))}
        self.assertTrue(all(r.status == 'UNRESOLVED' for r in res.values()))                # no repeat genotyper output: nothing can be said about length
        self.assertEqual(res['TEST1_GENEA'].thresholds, {'normal_max': 35, 'pathologic_min': 40})

    def test_build_is_idempotent(self):
        out2 = tempfile.mkdtemp()
        try:
            build(out2, [lambda g: StrchiveAdapter(g, root=os.path.join(ROOT, 'strchive'))])
            with open(os.path.join(self.out, 'table_digests.json')) as f1, open(os.path.join(out2, 'table_digests.json')) as f2:
                self.assertEqual(json.load(f1), json.load(f2))
        finally:
            shutil.rmtree(out2, ignore_errors=True)


if __name__ == '__main__':
    unittest.main()
