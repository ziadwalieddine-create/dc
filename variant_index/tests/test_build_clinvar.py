import json
import os
import shutil
import sys
import tempfile
import unittest

import pyarrow.compute as pc

from helpers import syn_reference, vcf_text, FIXTURES
from vindex.build import runner
from vindex.build.clinvar_vcf import ClinVarVcfAdapter
from vindex.build.genes import GeneMap
from vindex.index_reader import Index
from vindex.vcfio import read_vcf
from vindex.match_small import MatchContext, status_of, interpret_calls


sys.path.insert(0, FIXTURES)
import make_sources  # noqa: E402

SOURCES = make_sources.write_all(tempfile.mkdtemp())


class MiniClinVarAdapter(ClinVarVcfAdapter):
    def run(self, ctx, limit=None):
        return super().run(ctx, path=os.path.join(SOURCES, 'mini_clinvar.vcf'), limit=limit)


class TestClinVarSlice(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.ref = syn_reference()
        cls.out = tempfile.mkdtemp()
        cls.man = runner.build(cls.out, [MiniClinVarAdapter(GeneMap(cls.ref, ot_dir='/nonexistent', mane_path='/nonexistent'))], cls.ref)
        cls.idx = Index.open(cls.out)

    @classmethod
    def tearDownClass(cls):
        shutil.rmtree(cls.out, ignore_errors=True)

    def evidence(self):
        return self.idx.table('evidence')

    def test_row_accounting_preserves_every_record(self):
        # regression: rows without a chromosome used to sit in a '_' folder that Arrow dataset discovery silently skips
        st = self.man['adapters'][0]['stats']
        self.assertEqual(st['records'], 15)
        e = self.evidence()
        self.assertEqual(e.num_rows, sum(v for k, v in st.items() if k.startswith('evidence:')))
        status = {r['values']: r['counts'] for r in pc.value_counts(e['mapping_status']).to_pylist()}
        reasons = {r['values']: r['counts'] for r in pc.value_counts(e['mapping_reason']).to_pylist() if r['values']}
        self.assertEqual(status.get('rejected'), 1)            # REF mismatch
        self.assertEqual(reasons, {'ref_mismatch': 1, 'no_alt_allele': 1, 'contig_unmapped': 1, 'ref_n': 1})
        self.assertEqual(status.get('unmapped'), 3)
        self.assertEqual(status.get('mapped'), e.num_rows - 4)

    def test_records_without_any_classification_are_kept_as_present_in_clinvar(self):
        rows = [r for r in self.evidence().to_pylist() if r['claim_type'] == 'clinvar_record_without_classification']
        self.assertEqual(sorted(r['accession'].split('|')[0] for r in rows), ['VariationID:1014', 'VariationID:1015'])
        self.assertTrue(all(r['effect_assertion'] == 'not_applicable' and r['evidence_kind'] == 'aggregated_integration' and r['mapping_status'] == 'mapped' for r in rows))
        self.assertEqual(self.man['adapters'][0]['stats']['records_without_classification'], 2)
        keys = {v['variant_key'] for v in self.idx.table('entity_variant').to_pylist()}
        self.assertTrue({r['subject_id'] for r in rows} <= keys)                      # such a variant is an entity, not absent from the index

    def test_dbsnp_ids_are_exported_per_record_never_carried_over_from_another_record(self):
        class Probe:
            source_id = 'probe'
            requires = ('clinvar',)

            def prepare(self, ctx):
                ctx.want('dbsnp', ['rs111', 'rs222', 'rs333', 'rs999'])

            def run(self, ctx, limit=None):
                self.found = {k: list(v) for k, v in ctx.found.items() if k[0] == 'dbsnp'}
                return {}

        probe = Probe()
        out = tempfile.mkdtemp()
        try:
            runner.build(out, [MiniClinVarAdapter(GeneMap(self.ref, ot_dir='/nonexistent', mane_path='/nonexistent')), probe], self.ref)
        finally:
            shutil.rmtree(out, ignore_errors=True)
        b = lambda p: self.ref.fetch('1', p, p)
        a = lambda p: next(x for x in 'ACGT' if x != b(p) and x != {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}[b(p)])
        self.assertEqual(probe.found, {('dbsnp', 'rs111'): [f'1:402:{b(402)}:{a(402)}'], ('dbsnp', 'rs222'): [f'1:490:{b(490)}:{a(490)}'], ('dbsnp', 'rs333'): [f'1:490:{b(490)}:{a(490)}']})
        xr = [r for r in self.idx.table('xref').to_pylist() if r['namespace'] == 'dbsnp']
        self.assertEqual(sorted(r['identifier'] for r in xr), ['rs111', 'rs222', 'rs333'])

    def test_normalisation_deduplicates_identity_and_keeps_both_records(self):
        rows = [r for r in self.evidence().to_pylist() if r['subject_id'] == '1:101:CA:C' and r['claim_type'] == 'germline_classification']
        self.assertEqual(sorted(r['accession'].split('|')[0] for r in rows), ['VariationID:1005', 'VariationID:1006'])
        self.assertEqual(sorted(r['claim_value'] for r in rows), ['Likely_pathogenic', 'Pathogenic'])
        ev = self.idx.table('entity_variant').to_pylist()
        self.assertEqual(sum(1 for v in ev if v['variant_key'] == '1:101:CA:C'), 1)          # one entity, two evidence rows
        self.assertEqual(len(ev), len({v['variant_key'] for v in ev}))

    def test_axes_assertions_and_review_stars(self):
        rows = {r['accession'].split('|')[0]: r for r in self.evidence().to_pylist() if r['claim_type'] == 'germline_classification'}
        self.assertEqual((rows['VariationID:1001']['effect_assertion'], rows['VariationID:1001']['review_stars']), ('effect', 3))
        self.assertEqual(rows['VariationID:1002']['effect_assertion'], 'no_effect')
        self.assertEqual(rows['VariationID:1003']['effect_assertion'], 'uncertain')
        self.assertEqual(rows['VariationID:1004']['effect_assertion'], 'conflicting')
        self.assertIn('Pathogenic(1)', rows['VariationID:1004']['extra_json'])
        self.assertEqual(rows['VariationID:1006']['review_stars'], 0)
        self.assertEqual(rows['VariationID:1001']['condition_ids'], ['MedGen:C0000001', 'OMIM:100001'])
        self.assertTrue(all(r['axis'] == 'clinical_validity_actionability' for r in rows.values()))

    def test_somatic_and_population_rows_are_separate_claims(self):
        e = self.evidence().to_pylist()
        onc = [r for r in e if r['claim_type'] == 'oncogenicity_classification']
        sci = [r for r in e if r['claim_type'] == 'somatic_clinical_impact']
        self.assertEqual((len(onc), len(sci)), (1, 1))
        self.assertEqual((onc[0]['origin'], onc[0]['effect_assertion']), ('somatic', 'effect'))
        af = [r for r in e if r['claim_type'] == 'population_allele_frequency']
        self.assertEqual(sorted(r['claim_value'] for r in af), ['1000G', 'ExAC', 'GO-ESP'])
        self.assertTrue(all(r['axis'] == 'applicability_context' and r['effect_assertion'] == 'not_applicable' for r in af))

    def test_provenance_rows_resolve_and_digest_matches_source(self):
        path = os.path.join(SOURCES, 'mini_clinvar.vcf')
        with open(path) as fh:
            lines = fh.read().split('\n')
        for r in self.evidence().to_pylist():
            ln = int(r['locator'][1:])
            import hashlib
            self.assertEqual(hashlib.sha256(lines[ln - 1].encode()).hexdigest(), r['raw_sha256'])
            self.assertTrue(r['release_id'].startswith('clinvar/clinvar_vcf_grch38@2026-01-01'))
        self.assertEqual(self.man['releases'], ['clinvar/clinvar_vcf_grch38@2026-01-01'])

    def test_build_is_idempotent(self):
        out2 = tempfile.mkdtemp()
        try:
            runner.build(out2, [MiniClinVarAdapter(GeneMap(self.ref, ot_dir='/nonexistent', mane_path='/nonexistent'))], self.ref)
            with open(os.path.join(self.out, 'table_digests.json')) as f1, open(os.path.join(out2, 'table_digests.json')) as f2:
                self.assertEqual(json.load(f1), json.load(f2))
        finally:
            shutil.rmtree(out2, ignore_errors=True)

    def test_matcher_on_built_index(self):
        base = lambda p: self.ref.fetch('1', p, p)
        alt = lambda p: next(x for x in 'ACGT' if x != base(p) and x != {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}[base(p)])
        recs = [('1', 402, base(402), alt(402), 'PASS', 50, ['0/1:60:30:15,15']),
                ('1', 105, 'AA', 'A', 'PASS', 50, ['1/1:60:30:0,30']),            # right-shifted deletion: matches both ClinVar records after normalisation
                ('1', 430, base(430), alt(430), 'PASS', 50, ['0/1:60:30:15,15']),
                ('1', 600, base(600), alt(600), 'PASS', 50, ['0/1:60:30:15,15'])]       # not in the index
        p = os.path.join(tempfile.mkdtemp(), 's.vcf')
        with open(p, 'w') as f:
            f.write(vcf_text(sorted(recs, key=lambda r: r[1])))
        cs = read_vcf(p, self.ref)
        out = {o.call.key: o for o in interpret_calls(cs, self.idx, MatchContext(cs))}
        self.assertEqual(out[f'1:402:{base(402)}:{alt(402)}'].status, 'EXACT_MATCH')
        o = out['1:101:CA:C']
        self.assertEqual((o.status, o.axes['clinical_validity_actionability']['n_rows'], o.axes['clinical_validity_actionability']['n_independent']), ('EXACT_MATCH', 2, 2))
        self.assertTrue(o.axes['clinical_validity_actionability']['conflict'] is False)      # P and LP are both 'effect'
        c = out[f'1:430:{base(430)}:{alt(430)}']
        self.assertEqual(c.axes['clinical_validity_actionability']['assertions'], {'conflicting': 1})
        self.assertEqual(out[f'1:600:{base(600)}:{alt(600)}'].status, 'NO_INDEX_RECORD')
        self.assertEqual(out[f'1:600:{base(600)}:{alt(600)}'].evidence, [])


class TestClinVarPlaceholders(unittest.TestCase):
    """CLNDN=. is how ClinVar writes a record that names no condition: the lone '.' is the file's marker for no value, not a condition text. A bit of ORIGIN the
    header does not define must stay visible. (Found by comparing the index with the file as read by bcftools: 80 of 4,356,442 real rows carried the '.'.)"""

    @classmethod
    def setUpClass(cls):
        cls.ref = syn_reference()
        cls.out = tempfile.mkdtemp()
        vcf = os.path.join(cls.out, 'placeholders.vcf')
        rows = [(402, 1, dict(CLNDN='.', CLNDISDB='.', ORIGIN='1025')),
                (410, 2, dict(CLNDN='Cond_A|.', CLNDISDB='MedGen:C0000001|.', ORIGIN='1073741825')),
                (420, 3, dict(CLNDN='Test_syndrome', CLNDISDB='MedGen:C0000001,OMIM:100001', ORIGIN='0')),
                (430, 4, dict(CLNDN='Somatic_x', CLNDISDB='MedGen:C0000002', ORIGIN='3')),
                (440, 5, dict(CLNDN='Dual', CLNDISDB='MedGen:C0000003', ORIGIN='1', ONC='Oncogenic', ONCREVSTAT='criteria_provided,_single_submitter', CLNSIGSCV='SCV000000051|SCV000000052')),
                (450, 6, dict(CLNDN='Odd', CLNDISDB='MedGen:C0000004', ORIGIN='1', CLNREVSTAT='a_status_clinvar_has_not_used'))]
        with open(vcf, 'w') as f:
            f.write('##fileformat=VCFv4.1\n##fileDate=2026-02-02\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n')
            for pos, vid, info in rows:
                info = {**dict(ALLELEID=5000 + vid, CLNREVSTAT='criteria_provided,_single_submitter', CLNSIG='Pathogenic', CLNSIGSCV=f'SCV{vid:09d}'), **info}
                f.write('\t'.join(['1', str(pos), str(vid), make_sources.b('1', pos), make_sources.alt('1', pos), '.', '.', ';'.join(f'{k}={v}' for k, v in info.items())]) + '\n')

        class Adapter(ClinVarVcfAdapter):
            def run(self, ctx, limit=None):
                return super().run(ctx, path=vcf, limit=limit)
        cls.man = runner.build(os.path.join(cls.out, 'idx'), [Adapter(GeneMap(cls.ref, ot_dir='/nonexistent', mane_path='/nonexistent'))], cls.ref)
        cls.all = Index.open(os.path.join(cls.out, 'idx')).table('evidence').to_pylist()
        cls.rows = {r['accession'].split('|')[0].split(':')[1]: r for r in cls.all if r['claim_type'] == 'germline_classification'}

    @classmethod
    def tearDownClass(cls):
        shutil.rmtree(cls.out, ignore_errors=True)

    def test_a_lone_dot_is_no_condition_but_a_dot_inside_a_list_keeps_its_place(self):
        self.assertIsNone(self.rows['1']['condition_text'])
        self.assertEqual(self.rows['1']['condition_ids'], [])
        self.assertEqual(self.rows['2']['condition_text'], 'Cond A|.')                    # positions in CLNDN and CLNDISDB stay aligned
        self.assertEqual(self.rows['2']['condition_ids'], ['MedGen:C0000001'])
        self.assertEqual(self.rows['3']['condition_text'], 'Test syndrome')
        self.assertEqual(self.rows['3']['condition_ids'], ['MedGen:C0000001', 'OMIM:100001'])

    def test_the_review_status_is_spelled_as_the_vcv_tables_spell_it(self):
        # the VCF writes criteria_provided,_single_submitter; the VCV tables write criteria provided, single submitter: one status, one value in the index
        self.assertEqual({r['review_status'] for r in self.rows.values() if r['accession'].split(':')[1].startswith(('1|', '2|', '3|'))}, {'criteria provided, single submitter'})
        self.assertEqual(self.rows['3']['review_stars'], 1)

    def test_a_row_counts_the_submissions_of_its_own_classification_only(self):
        onc = next(r for r in self.all if r['claim_type'] == 'oncogenicity_classification' and r['accession'].startswith('VariationID:5|'))
        self.assertIsNone(onc['n_submitters'])                                                 # no ONCSCV: the two germline SCVs are not its submissions
        self.assertEqual(self.rows['5']['n_submitters'], 2)

    def test_an_unknown_review_status_is_counted_in_the_build_manifest(self):
        st = self.man['adapters'][0]['stats']
        self.assertEqual(st.get('unrecognised_review_status:a_status_clinvar_has_not_used'), 1)
        self.assertIsNone(self.rows['6']['review_stars'])
        self.assertEqual([k for k in st if k.startswith('unrecognised_review_status:')], ['unrecognised_review_status:a_status_clinvar_has_not_used'])

    def test_origin_bits_are_named_and_an_undefined_bit_stays_visible(self):
        self.assertEqual(self.rows['1']['origin'], 'germline+unrecognised_bit_1024')
        self.assertEqual(self.rows['2']['origin'], 'germline+other')
        self.assertEqual(self.rows['3']['origin'], 'unknown')
        self.assertEqual(self.rows['4']['origin'], 'germline+somatic')


if __name__ == '__main__':
    unittest.main()
