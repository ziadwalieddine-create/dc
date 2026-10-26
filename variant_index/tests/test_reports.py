import glob
import json
import os
import shutil
import tempfile
import unittest

import pyarrow as pa
import pyarrow.parquet as pq

from helpers import syn_reference
from test_wgs import build_index
from vindex import reports


class TestReports(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.work = tempfile.mkdtemp()
        cls.idx = os.path.join(cls.work, 'idx')
        build_index(cls.idx)
        cls.cov = reports.coverage(cls.idx)

    @classmethod
    def tearDownClass(cls):
        shutil.rmtree(cls.work, ignore_errors=True)

    def test_coverage_counts_agree_with_the_tables(self):
        ev = self.cov['evidence']
        self.assertEqual(ev['total'], self.cov['tables']['evidence']['rows'])
        self.assertEqual(sum(ev['source_id'].values()), ev['total'])
        self.assertEqual(sum(ev['mapping_status'].values()), ev['total'])
        self.assertEqual(sum(ev['source_by_axis'].values()), ev['total'])
        self.assertGreater(ev['mapping_status']['mapped'], 0)
        self.assertIn('clinvar', ev['source_id'])
        self.assertIn('clinical_validity_actionability', ev['axis'])
        self.assertGreater(ev['variants_with_evidence_by_axis']['clinical_validity_actionability'], 3)

    def test_rows_kept_but_not_mapped_carry_their_reasons(self):
        nm = self.cov['evidence']['not_mapped']
        self.assertTrue(any(k.startswith('clinvar|unmapped|only_grch37_location_no_liftover_performed') for k in nm))
        self.assertTrue(any(k.startswith('clinvar|rejected|ref_mismatch') for k in nm))
        self.assertTrue(all(k.count('|') == 2 and not k.endswith('|None') for k in nm))
        gap = reports.gap_report(self.idx, life={'targets': [], 'acquired_not_integrated': []}, cov=self.cov)
        self.assertGreater(gap['rows_kept_but_not_mapped']['clinvar']['total'], 0)

    def test_a_clean_index_validates_including_reference_alleles(self):
        v = reports.validate(self.idx, syn_reference())
        self.assertTrue(v['ok'], {k: c for k, c in v['checks'].items() if c['violations']})
        self.assertGreater(v['evidence_rows_checked'], 0)
        self.assertGreater(v['entity_variant_ref_checked'], 0)

    def damaged(self, mutate):
        d = os.path.join(self.work, 'damaged')
        shutil.rmtree(d, ignore_errors=True)
        shutil.copytree(self.idx, d)
        mutate(d)
        return reports.validate(d, syn_reference())

    def rewrite(self, d, pattern, col, fn):
        f = sorted(glob.glob(os.path.join(d, pattern)))[0]
        t = pq.read_table(f)
        vals = t[col].to_pylist()
        vals = fn(vals)
        pq.write_table(t.set_column(t.schema.get_field_index(col), col, pa.array(vals, t.schema.field(col).type)), f)

    def test_vocabulary_range_and_format_violations_are_caught(self):
        v = self.damaged(lambda d: self.rewrite(d, 'evidence/variant/1/clinvar.parquet', 'effect_assertion', lambda xs: ['bogus'] + xs[1:]))
        self.assertEqual(v['checks']['effect_assertion_outside_vocabulary']['violations'], 1)
        self.assertFalse(v['ok'])
        v = self.damaged(lambda d: self.rewrite(d, 'evidence/variant/1/clinvar.parquet', 'p_value', lambda xs: [1.5] + xs[1:]))
        self.assertEqual(v['checks']['p_value_outside_0_1']['violations'], 1)
        v = self.damaged(lambda d: self.rewrite(d, 'evidence/variant/1/clinvar.parquet', 'raw_sha256', lambda xs: ['xyz'] + xs[1:]))
        self.assertEqual(v['checks']['raw_sha256_not_64_hex']['violations'], 1)

    def test_mapped_rows_need_a_subject_and_a_wellformed_key(self):
        def drop_subject(d):
            def fn(xs):
                i = next(i for i, x in enumerate(xs) if x)
                return xs[:i] + [None] + xs[i + 1:]
            self.rewrite(d, 'evidence/variant/1/clinvar.parquet', 'subject_id', fn)
        self.assertGreaterEqual(self.damaged(drop_subject)['checks']['mapped_without_subject']['violations'], 1)

        def bad_key(d):
            def fn(xs):
                i = next(i for i, x in enumerate(xs) if x and x.count(':') == 3)
                return xs[:i] + ['1:5:A'] + xs[i + 1:]
            self.rewrite(d, 'evidence/variant/1/clinvar.parquet', 'subject_id', fn)
        self.assertGreaterEqual(self.damaged(bad_key)['checks']['variant_key_not_chrom_pos_ref_alt']['violations'], 1)

    def test_a_mapped_row_filed_under_the_wrong_chromosome_is_caught(self):
        d = os.path.join(self.work, 'misfiled')
        shutil.rmtree(d, ignore_errors=True)
        shutil.copytree(self.idx, d)
        src = os.path.join(d, 'evidence', 'variant', '1', 'clinvar.parquet')
        os.makedirs(os.path.join(d, 'evidence', 'variant', '2'), exist_ok=True)
        shutil.move(src, os.path.join(d, 'evidence', 'variant', '2', 'clinvar.parquet'))
        v = reports.validate(d, syn_reference())
        self.assertGreater(v['checks']['mapped_variant_row_filed_under_another_chromosome']['violations'], 0)
        self.assertFalse(v['ok'])

    def test_duplicate_ids_and_orphan_subjects_are_caught(self):
        v = self.damaged(lambda d: self.rewrite(d, 'evidence/variant/1/clinvar.parquet', 'evidence_id', lambda xs: [xs[0]] + [xs[0]] + xs[2:]))
        self.assertEqual(v['checks']['evidence_id_not_unique']['violations'], 1)
        v = self.damaged(lambda d: self.rewrite(d, 'entity_variant/1.parquet', 'variant_key', lambda xs: ['9:9:A:C'] + xs[1:]))
        self.assertGreaterEqual(v['checks']['mapped_variant_subject_not_an_entity']['violations'], 1)

    def test_a_reference_allele_that_differs_from_the_fasta_is_caught(self):
        def corrupt(d):
            self.rewrite(d, 'entity_variant/1.parquet', 'ref', lambda xs: ['A' if x != 'A' else 'C' for x in xs])
        v = self.damaged(corrupt)
        self.assertGreater(v['checks']['entity_variant_ref_differs_from_reference']['violations'], 0)

    def test_lifecycle_separates_listed_acquired_verified_mapped_and_integrated(self):
        tp, mp, vp = (os.path.join(self.work, n) for n in ('targets.json', 'manifest.jsonl', 'validation.json'))
        rels = {r['target_id'] for r in pq.read_table(os.path.join(self.idx, 'source_release.parquet')).to_pylist()}
        some = sorted(rels)[0]
        ids = [some, 'listed_only', 'acquired_unverified', 'acquired_verified_not_integrated']
        with open(tp, 'w') as f:
            json.dump([dict(id=i, source='clinvar', dest='x') for i in ids], f)
        with open(mp, 'w') as f:
            f.write('\n'.join(json.dumps(dict(id=i, ok=True)) for i in (some, 'acquired_unverified', 'acquired_verified_not_integrated')) + '\n')
        with open(vp, 'w') as f:
            json.dump({some: dict(passed=True), 'acquired_unverified': dict(passed=False), 'acquired_verified_not_integrated': dict(passed=True)}, f)
        life = reports.lifecycle(self.idx, tp, mp, vp, catalogue_dir=os.path.join(self.work, 'none'))
        by = {r['target_id']: r for r in life['targets']}
        flags = lambda r: (r['listed'], r['acquired'], r['verified'], r['integrated'], r['mapped'])
        self.assertEqual(flags(by[some]), (True, True, True, True, True))
        self.assertEqual(flags(by['listed_only']), (True, False, False, False, False))
        self.assertEqual(flags(by['acquired_unverified']), (True, True, False, False, False))
        self.assertEqual(flags(by['acquired_verified_not_integrated']), (True, True, True, False, False))
        self.assertEqual(life['acquired_not_integrated'], ['acquired_unverified', 'acquired_verified_not_integrated'])
        self.assertEqual(life['summary']['total'], 4)

    def test_change_log_names_what_changed_between_builds(self):
        other = os.path.join(self.work, 'idx2')
        shutil.copytree(self.idx, other)
        with open(os.path.join(other, 'index_manifest.json')) as f:
            m = json.load(f)
        m['releases'] = [r for r in m['releases'] if not r.startswith('reactome/')] + ['clinvar/clinvar_vcf_grch38@2099-01-01']
        m['releases'] = [r for r in m['releases'] if r != 'clinvar/clinvar_vcf_grch38@2026-01-01']
        m['adapters'][0]['code_sha256_16'] = 'changed'
        m['tables']['evidence']['rows'] += 5
        with open(os.path.join(other, 'index_manifest.json'), 'w') as f:
            json.dump(m, f)
        c = reports.change_log(self.idx, other)
        self.assertEqual(c['releases']['removed'], ['reactome/reactome_ensembl2reactome'])
        self.assertEqual(c['releases']['reversioned']['clinvar/clinvar_vcf_grch38'], dict(previous='clinvar/clinvar_vcf_grch38@2026-01-01', current='clinvar/clinvar_vcf_grch38@2099-01-01'))
        self.assertEqual(c['adapters']['code_changed'], [m['adapters'][0]['adapter']])
        self.assertEqual(c['tables']['evidence']['delta'], 5)
        self.assertEqual(reports.change_log(self.idx, self.idx)['tables'], {})

    def test_the_module_registry_names_what_is_not_implemented(self):
        status = {m['name']: m['status'] for m in reports.MODULES}
        self.assertEqual(status['personal 30x GRCh38 SNP/indel, CNV and SV calls'], 'unavailable_input')
        self.assertEqual(status['transcript HGVS (c., n.) to genomic coordinate mapping'], 'unimplemented')
        self.assertTrue(all(m['status'] in ('implemented', 'partial', 'unimplemented', 'external', 'not_acquired', 'unavailable_input') for m in reports.MODULES))
        self.assertTrue(all(m.get('note') for m in reports.MODULES if m['status'] != 'implemented'))


if __name__ == '__main__':
    unittest.main()
