import hashlib
import json
import os
import shutil
import sys
import tempfile
import unittest

import pyarrow.parquet as pq

from helpers import syn_reference, FIXTURES
from vindex.build import runner
from vindex.build.clinvar_vcf import ClinVarVcfAdapter
from vindex.build.clinvar_vcv import ClinVarVcvAdapter, canon
from vindex.build.genes import GeneMap
from vindex.index_reader import Index
from vindex.match_hap import call_gene

sys.path.insert(0, FIXTURES)
import make_sources  # noqa: E402

SOURCES = make_sources.write_all(tempfile.mkdtemp())
VCV = os.path.join(SOURCES, 'vcv')


class MiniVcf(ClinVarVcfAdapter):
    def run(self, ctx, limit=None):
        return super().run(ctx, path=os.path.join(SOURCES, 'mini_clinvar.vcf'), limit=limit)


class Seed:
    """Stands in for the VCF adapter and offers exactly the identifiers a test needs resolved."""
    source_id = 'clinvar'

    def __init__(self, offers):
        self.offers = offers

    def run(self, ctx, limit=None):
        ctx.adapter = self
        for ns, ident, entity in self.offers:
            ctx.offer(ns, ident, entity)
        return {}


def k(pos):
    ref = make_sources.REF
    base = ref.fetch('1', pos - 1, pos).upper()
    return f'1:{pos}:{base}:{make_sources.alt("1", pos)}'


def build(out, first=None):
    ref = syn_reference()
    genes = GeneMap(ref, ot_dir='/nonexistent', mane_path='/nonexistent')
    return runner.build(out, [first or MiniVcf(genes), ClinVarVcvAdapter(genes, root=VCV)], ref)


class TestVcvSlice(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.out = tempfile.mkdtemp()
        cls.man = build(cls.out)
        cls.idx = Index.open(cls.out)
        cls.stats = next(a['stats'] for a in cls.man['adapters'] if a['adapter'] == 'ClinVarVcvAdapter')
        cls.rows = [r for r in cls.idx.table('evidence').to_pylist() if r['release_id'].startswith('clinvar/clinvar_vcv')]
        cls.scv = {r['accession']: r for r in cls.rows if r['accession'].startswith('SCV')}
        cls.agg = {(r['accession'], r['claim_type']): r for r in cls.rows if r['accession'].startswith('VCV')}
        cls.events = {r['event_id']: r for r in cls.idx.table('entity_event').to_pylist()}

    @classmethod
    def tearDownClass(cls):
        shutil.rmtree(cls.out, ignore_errors=True)

    def s(self, n):
        return self.scv[f'SCV{n:09d}']

    # ---- accounting ---------------------------------------------------------------------------------------------------------------
    def test_every_submission_is_kept_with_a_mapping_status_and_reason(self):
        self.assertEqual(self.stats['scv_rows'], make_sources.write_all.vcv['scv_count'])
        self.assertEqual(len(self.scv), self.stats['scv_rows'])
        by = {}
        for r in self.scv.values():
            by[r['mapping_status']] = by.get(r['mapping_status'], 0) + 1
            if r['mapping_status'] != 'mapped':
                self.assertTrue(r['mapping_reason'] and r['unmapped_key'].startswith('VariationID:') and r['subject_id'] is None)
        self.assertEqual(by, {'mapped': 16, 'unmapped': 6})
        self.assertEqual(self.stats['scv_evidence:curated_expert_classification:mapped'] + self.stats['scv_evidence:experimental_annotation:mapped']
                         + self.stats['scv_evidence:case_observation:mapped'], 16)

    def test_release_is_the_vcv_release_not_the_vcf(self):
        self.assertEqual({r['release_id'] for r in self.rows}, {'clinvar/clinvar_vcv_xml_weekly_2026_1018@2026-01-02'})
        self.assertIn('clinvar/clinvar_vcf_grch38@2026-01-01', self.man['releases'])
        self.assertEqual(self.man['problems'], [])

    def test_variants_the_vcf_adapter_handled_get_no_second_aggregate_row(self):
        accessions = {r['accession'] for r in self.agg.values()}
        self.assertFalse({f'VCV{n:09d}' for n in range(1001, 1016)} & accessions)        # mapped, rejected and unmapped VCF records alike
        self.assertEqual(self.stats['aggregate_rows:event:mapped'], 6)

    # ---- submissions --------------------------------------------------------------------------------------------------------------
    def test_submission_attaches_to_the_vcf_key_and_keeps_its_own_fields(self):
        r = self.s(1)
        self.assertEqual((r['level'], r['subject_id'], r['chrom'], r['pos']), ('variant', k(402), '1', 402))
        self.assertEqual((r['claim_type'], r['claim_value'], r['effect_assertion'], r['evidence_kind'], r['axis']),
                         ('scv_germline_classification', 'Pathogenic', 'effect', 'curated_expert_classification', 'clinical_validity_actionability'))
        self.assertEqual((r['review_status'], r['review_stars'], r['curator'], r['assertion_criteria']), ('reviewed by expert panel', 3, 'Lab 1', 'ACMG 2015'))
        self.assertEqual((r['condition_ids'], r['condition_text'], r['sample_size'], r['origin']), (['MedGen:C0000001', 'OMIM:100001'], 'Test syndrome', 52, 'germline'))
        self.assertEqual(json.loads(r['extra_json'])['citations'], 'PubMed:111|PubMed:222')

    def test_submissions_and_their_aggregate_share_one_lineage_group(self):
        vcf_agg = [r for r in self.idx.table('evidence').to_pylist() if r['accession'].startswith('VariationID:1001|') and r['claim_type'] == 'germline_classification']
        groups = {r['lineage_group'] for r in vcf_agg} | {self.s(n)['lineage_group'] for n in (1, 2)}
        self.assertEqual(groups, {'clinvar:variation:1001'})

    def test_disagreeing_submitters_stay_visible_as_separate_rows(self):
        a, b = self.s(6), self.s(7)
        self.assertEqual((a['effect_assertion'], b['effect_assertion']), ('effect', 'no_effect'))
        self.assertEqual((a['subject_id'], b['subject_id']), (k(430), k(430)))
        self.assertEqual(b['condition_ids'], ['HP:0002282', 'MedGen:CN002072'])       # the source's repeated 'HP:HP:' prefix is read as one
        self.assertEqual(b['claim_value'], 'Likely Benign')                             # the label stays as submitted

    def test_submission_flags_and_unparsed_counts_are_kept_not_guessed(self):
        r = self.s(2)
        self.assertEqual((r['claim_value'], r['effect_assertion'], r['origin'], r['sample_size']), ('likely pathogenic', 'effect', 'de_novo', None))
        ex = json.loads(r['extra_json'])
        self.assertEqual((ex['not_in_aggregate'], ex['number_tested']), (True, '118|331'))
        self.assertEqual(self.s(4)['origin'], 'germline')                               # 'germline|not applicable'

    def test_label_vocabulary_is_read_not_assumed(self):
        self.assertEqual(self.s(20)['effect_assertion'], 'effect')                      # 'Pathogenic, low penetrance'
        self.assertEqual(self.s(22)['effect_assertion'], 'not_applicable')              # 'not provided'
        self.assertEqual(self.s(21)['effect_assertion'], 'uncertain')                   # unknown label
        self.assertEqual(self.stats['unrecognised_classification:mystery_class'], 1)

    def test_somatic_submissions_are_separate_claims_with_somatic_origin(self):
        onc, sci = self.s(12), self.s(13)
        self.assertEqual((onc['claim_type'], onc['origin'], onc['effect_assertion']), ('scv_oncogenicity_classification', 'somatic', 'effect'))
        self.assertEqual((sci['claim_type'], sci['origin'], sci['effect_assertion']), ('scv_somatic_clinical_impact', 'somatic', 'effect'))

    def test_submissions_without_classification_are_observations_not_classifications(self):
        f, p, two = self.s(30), self.s(31), self.s(32)
        self.assertEqual((f['axis'], f['evidence_kind'], f['claim_type'], f['claim_value'], f['effect_assertion']),
                         ('empirical_functional_support', 'experimental_annotation', 'scv_observation_without_classification', 'in vitro', 'not_applicable'))
        self.assertEqual(f['lineage_group'], 'pmid:31275557')                           # a published study is the unit of independence, whatever database repeats it
        self.assertIsNone(f['origin'])
        self.assertEqual((p['axis'], p['evidence_kind'], p['lineage_group']), ('applicability_context', 'case_observation', 'clinvar:scv:SCV000000031'))
        self.assertEqual(two['lineage_group'], 'clinvar:scv:SCV000000032')               # two citations: no single study to share
        self.assertTrue(all(r['effect_assertion'] == 'not_applicable' for r in (f, p, two)))

    def test_unmapped_submissions_say_why(self):
        self.assertEqual((self.s(42)['level'], self.s(42)['mapping_reason']), ('event', 'only_grch37_location_no_liftover_performed'))
        self.assertEqual(self.s(43)['mapping_reason'], 'no_vcf_key_or_structural_event_for_variation')
        self.assertEqual(self.s(44)['mapping_reason'], 'vcf_record_unmapped:ref_n')
        self.assertEqual(self.s(45)['mapping_reason'], 'vcf_record_unmapped:ref_mismatch')
        self.assertEqual(self.s(46)['unmapped_key'], 'VariationID:9999')
        self.assertEqual(self.s(47)['mapping_reason'], 'no_vcf_key_or_structural_event_for_variation')

    # ---- structural events -------------------------------------------------------------------------------------------------------
    def test_precise_deletion_event(self):
        e = self.events['clinvar:VCV000002001']
        self.assertEqual((e['event_type'], e['chrom'], e['start'], e['end'], e['svlen'], e['precise'], e['assembly'], e['origin']), ('DEL', '1', 1000, 1099, 100, True, 'GRCh38', 'clinical_submission'))
        self.assertEqual((e['start_ci_lo'], e['start_ci_hi'], e['end_ci_lo'], e['end_ci_hi']), (None, None, None, None))

    def test_uncertain_extent_uses_inner_as_nominal_and_outer_as_interval(self):
        e = self.events['clinvar:VCV000002002']
        self.assertEqual((e['event_type'], e['start'], e['end'], e['precise'], e['copy_number'], e['svlen']), ('CNV_GAIN', 1500, 1799, False, 3, 300))
        self.assertEqual((e['start_ci_lo'], e['start_ci_hi'], e['end_ci_lo'], e['end_ci_hi']), (-10, 0, 0, 11))

    def test_unbounded_limits_leave_the_interval_open(self):
        e = self.events['clinvar:VCV000002003']
        self.assertEqual((e['event_type'], e['start'], e['end'], e['precise'], e['copy_number']), ('CNV_LOSS', 2000, 2300, False, 1))
        self.assertEqual((e['start_ci_lo'], e['end_ci_hi']), (None, None))               # '?' is not a number

    def test_copy_number_is_not_guessed_from_an_ambiguous_name(self):
        e = self.events['clinvar:VCV000002006']
        self.assertEqual((e['chrom'], e['copy_number']), ('2', None))
        self.assertEqual(json.loads(e['extra_json'])['copy_state_text'], ')x2,3')

    def test_inversions_translocations_and_exclusions(self):
        self.assertEqual(self.events['clinvar:VCV000002007']['event_type'], 'INV')
        t = self.events['clinvar:VCV000002008']
        self.assertEqual((t['event_type'], t['start'], t['end'], t['svlen']), ('CTX', 2500, 2500, None))
        self.assertEqual(set(self.events), {f'clinvar:VCV00000200{n}' for n in (1, 2, 3, 6, 7, 8)})
        self.assertEqual(self.stats['locations_below_event_length'], 1)                   # the 30 bp deletion is a small variant, not an SV

    def test_event_level_aggregates_and_the_records_that_could_not_become_events(self):
        a = self.agg[('VCV000002003', 'germline_classification')]
        self.assertEqual((a['level'], a['subject_id'], a['mapping_status'], a['review_stars'], a['chrom'], a['pos'], a['end']), ('event', 'clinvar:VCV000002003', 'mapped', 3, '1', 2000, 2300))
        self.assertEqual(a['gene_ids'], ['NCBIGene:1001'])
        g37 = self.agg[('VCV000002004', 'germline_classification')]
        self.assertEqual((g37['level'], g37['mapping_status'], g37['mapping_reason'], g37['unmapped_key'], g37['subject_id']),
                         ('event', 'unmapped', 'only_grch37_location_no_liftover_performed', 'VariationID:2004', None))
        self.assertEqual(json.loads(g37['extra_json'])['legacy_location'], 'GRCh37:1:900-1000')
        for acc in ('VCV000002005', 'VCV000002009'):
            self.assertEqual(self.agg[(acc, 'germline_classification')]['mapping_reason'], 'no_grch38_event_location_of_at_least_50bp')

    def test_submission_on_a_keyless_event_attaches_to_the_event(self):
        r = self.s(41)
        self.assertEqual((r['level'], r['subject_id'], r['mapping_status']), ('event', 'clinvar:VCV000002001', 'mapped'))
        self.assertEqual(r['lineage_group'], 'clinvar:variation:2001')

    # ---- composite records -------------------------------------------------------------------------------------------------------
    def test_composite_records_are_haplotype_entities_with_resolved_and_unresolved_members(self):
        haps = {h['haplotype_id']: h for h in self.idx.haplotypes(system='clinvar')}
        self.assertEqual(set(haps), {f'clinvar:VCV00000300{n}' for n in (1, 2, 3, 4)})
        self.assertEqual((haps['clinvar:VCV000003001']['gene_symbol'], haps['clinvar:VCV000003001']['allele_level']), ('GENEA', 'composite:Haplotype'))
        self.assertIsNone(haps['clinvar:VCV000003002']['gene_symbol'])                  # two genes: not one
        defs = self.idx.haplotype_definitions(list(haps))
        d1 = {(d['variant_key'], d['definition_scope']) for d in defs['clinvar:VCV000003001']}
        self.assertEqual(d1, {(k(402), 'defining'), (k(410), 'defining')})
        d2 = {(d['variant_key'], d['hgvs'], d['definition_scope']) for d in defs['clinvar:VCV000003002']}
        self.assertEqual(d2, {(k(420), None, 'member:genotype'), (None, 'AlleleID:9999', 'unresolved:member_has_no_vcf_key')})      # a repeated allele is one definition row
        self.assertEqual({d['definition_scope'] for d in defs['clinvar:VCV000003004']}, {'member:phase_unknown'})
        self.assertEqual((self.stats['composite_members_resolved'], self.stats['composite_members_unresolved']), (7, 1))

    def test_composite_classification_is_about_the_combination(self):
        a = self.agg[('VCV000003002', 'germline_classification')]
        self.assertEqual((a['level'], a['subject_id'], a['claim_value'], a['review_stars']), ('haplotype', 'clinvar:VCV000003002', 'drug response', 4))
        self.assertIn('combination', json.loads(a['extra_json'])['note'])
        self.assertEqual(json.loads(a['extra_json'])['unresolved_members'], ['9999'])
        self.assertNotIn(('VCV000003004', 'germline_classification'), self.agg)         # nothing classified, nothing invented
        self.assertEqual(self.s(40)['level'], 'haplotype')

    def test_clinvar_composites_are_never_called_as_star_alleles(self):
        r = call_gene(None, self.idx, 'GENEA')
        self.assertEqual((r.status, r.basis), ('UNRESOLVED', 'no_haplotype_definitions_in_index'))
        r = call_gene(None, self.idx, 'GENEA', system='clinvar')
        self.assertEqual((r.status, r.basis), ('UNRESOLVED', 'system_is_not_a_star_allele_system'))

    # ---- provenance and determinism ---------------------------------------------------------------------------------------------
    def test_locators_reopen_the_source_row_and_the_digest_matches(self):
        cache = {}
        for r in self.rows:
            table, _, rest = r['locator'].partition('/')
            part, _, i = rest.partition('#')
            path = os.path.join(VCV, table, part)
            if path not in cache:
                cache[path] = pq.read_table(path).to_pylist()
            row = cache[path][int(i)]
            self.assertEqual(hashlib.sha256(canon(row).encode()).hexdigest(), r['raw_sha256'], r['locator'])
            self.assertIn(r['accession'], (row.get('scv'), row.get('accession')))

    def test_build_is_idempotent(self):
        out2 = tempfile.mkdtemp()
        try:
            build(out2)
            with open(os.path.join(self.out, 'table_digests.json')) as f1, open(os.path.join(out2, 'table_digests.json')) as f2:
                self.assertEqual(json.load(f1), json.load(f2))
        finally:
            shutil.rmtree(out2, ignore_errors=True)


class TestVcvResolution(unittest.TestCase):
    def test_ambiguous_variation_id_is_not_attached_to_either_key(self):
        out = tempfile.mkdtemp()
        try:
            build(out, Seed([('clinvar_variation', '1001', k(402)), ('clinvar_variation', '1001', k(410))]))
            r = next(r for r in Index.open(out).table('evidence').to_pylist() if r['accession'] == 'SCV000000001')
        finally:
            shutil.rmtree(out, ignore_errors=True)
        self.assertEqual((r['mapping_status'], r['mapping_reason'], r['subject_id'], r['unmapped_key']), ('ambiguous', 'variation_id_maps_to_several_keys', None, 'VariationID:1001'))
        self.assertEqual(json.loads(r['extra_json'])['candidate_keys'], sorted([k(402), k(410)]))

    def test_event_that_also_has_a_vcf_key_keeps_submissions_on_the_variant_and_links_the_key(self):
        out = tempfile.mkdtemp()
        try:
            build(out, Seed([('clinvar_variation', '2001', k(500))]))
            idx = Index.open(out)
            ev = idx.table('evidence').to_pylist()
            events = {e['event_id']: e for e in idx.table('entity_event').to_pylist()}
        finally:
            shutil.rmtree(out, ignore_errors=True)
        r = next(r for r in ev if r['accession'] == 'SCV000000041')
        self.assertEqual((r['level'], r['subject_id']), ('variant', k(500)))              # the exact allele is matchable by key
        self.assertEqual(json.loads(events['clinvar:VCV000002001']['extra_json'])['variant_key'], k(500))
        self.assertTrue(any(a['level'] == 'event' and a['subject_id'] == 'clinvar:VCV000002001' and a['claim_type'] == 'germline_classification' for a in ev))


if __name__ == '__main__':
    unittest.main()
