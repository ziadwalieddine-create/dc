"""The cross-check recomputes calls and evidence matches with bcftools only. It must agree with a clean run and must be able to fail."""
import glob
import gzip
import os
import shutil
import tempfile
import unittest

import pyarrow as pa
import pyarrow.compute as pc
import pyarrow.parquet as pq
import pysam

from helpers import find_tool, syn_reference
from test_wgs import b, alt, build_index, combined_clinvar_vcf, sample_vcf
from vindex import crosscheck, wgs


@unittest.skipUnless(find_tool('bcftools'), 'bcftools not available (set VINDEX_TOOLS_BIN)')
class TestCrosscheck(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.work = tempfile.mkdtemp()
        cls.idx = os.path.join(cls.work, 'idx')
        build_index(cls.idx)
        other = lambda c, p: [x for x in 'ACGT' if x not in (b(c, p), alt(c, p))][0]
        extra = [('1', 160, b('1', 160), f'{alt("1", 160)},{other("1", 160)}', 'PASS', 50, ['0/2:60:30:15,0,15']),            # the sample carries only the second alt
                 ('1', 170, b('1', 170), f'{alt("1", 170)},{other("1", 170)}', 'PASS', 50, ['1/2:60:30:0,15,15'])]           # carries both
        cls.vcf = sample_vcf(cls.work, extra=extra)
        cls.res = os.path.join(cls.work, 'res')
        cls.ref = syn_reference()
        wgs.run_sample(cls.vcf, cls.idx, cls.res, cls.ref, repo_root=os.path.join(cls.work, 'not_the_repo'))
        cv = os.path.join(cls.work, 'cv.vcf')
        combined_clinvar_vcf(cv)
        with open(cv) as f:
            lines = f.read().rstrip('\n').split('\n')
        head, body = [l for l in lines if l.startswith('#')], [l.split('\t') for l in lines if not l.startswith('#')]
        body.sort(key=lambda r: (r[0], int(r[1])))                         # tabix needs coordinate order
        tags = sorted({kv.split('=')[0] for r in body for kv in r[7].split(';')})
        head = head[:-1] + [f'##INFO=<ID={t},Number=.,Type=String,Description="x">' for t in tags] + head[-1:]            # a real ClinVar header declares its INFO tags
        with open(cv, 'w') as f:
            f.write('\n'.join(head + ['\t'.join(r) for r in body]) + '\n')
        pysam.tabix_index(cv, preset='vcf', force=True)
        cls.cv = cv + '.gz'
        with gzip.open(cls.cv, 'rt') as f:
            lines = f.read().rstrip('\n').split('\n')
        old = []
        for l in lines:
            if l.startswith('##INFO=<ID=ONC') or l.startswith('##INFO=<ID=SCI'):
                continue
            if not l.startswith('#'):
                c = l.split('\t')
                c[7] = ';'.join(kv for kv in c[7].split(';') if not kv.startswith(('ONC', 'SCI')))
                l = '\t'.join(c)
            old.append(l)
        cls.cv_chr = os.path.join(cls.work, 'cv_chr.vcf')                                      # the same records under UCSC-style contig names, in a header that defines no contigs
        with open(cls.cv_chr, 'w') as f:
            f.write('\n'.join(l if l.startswith('#') else 'chr' + l for l in lines) + '\n')
        pysam.tabix_index(cls.cv_chr, preset='vcf', force=True)
        cls.cv_chr += '.gz'
        cls.cv_old = os.path.join(cls.work, 'cv_old.vcf')
        with open(cls.cv_old, 'w') as f:
            f.write('\n'.join(old) + '\n')
        pysam.tabix_index(cls.cv_old, preset='vcf', force=True)
        cls.cv_old += '.gz'

    @classmethod
    def tearDownClass(cls):
        shutil.rmtree(cls.work, ignore_errors=True)

    def check(self, res=None):
        return crosscheck.crosscheck(self.vcf, res or self.res, self.ref, self.cv)

    def copy_results(self, name):
        d = os.path.join(self.work, name)
        shutil.copytree(self.res, d)
        return d

    def test_a_clean_run_agrees_and_a_refused_ref_mismatch_is_explained(self):
        r = self.check()
        self.assertTrue(r['ok'], r['examples'])
        c = r['calls']
        self.assertEqual(c['bcftools'], c['matcher'] + 1)
        self.assertEqual((c['only_bcftools'], c['only_bcftools_explained_by_rejected_calls'], c['only_bcftools_unexplained'], c['only_matcher']), (1, 1, 0, 0))        # 1:440, REF disagrees with the reference
        e = r['evidence']
        self.assertGreater(e['matcher'], 0)
        self.assertEqual((e['only_bcftools_unexplained'], e['only_matcher']), (0, 0))
        self.assertEqual(e['bcftools'], e['matcher'] + e['only_bcftools_explained_by_rejected_calls'])

    def test_multiallelic_records_contribute_only_the_alleles_the_sample_carries(self):
        calls = {r['key'] for r in pq.read_table(os.path.join(self.res, 'small_calls.parquet')).to_pylist()}
        a, o = alt('1', 160), [x for x in 'ACGT' if x not in (b('1', 160), alt('1', 160))][0]
        self.assertIn(f'1:160:{b("1", 160)}:{o}', calls)
        self.assertNotIn(f'1:160:{b("1", 160)}:{a}', calls)                                  # 0/2: the first alt is not carried
        self.assertEqual({k for k in calls if k.startswith('1:170:')}, {f'1:170:{b("1", 170)}:{alt("1", 170)}', f'1:170:{b("1", 170)}:{[x for x in "ACGT" if x not in (b("1", 170), alt("1", 170))][0]}'})
        self.assertTrue(self.check()['ok'])

    def test_a_dropped_call_is_reported_as_unexplained(self):
        d = self.copy_results('dropped')
        t = pq.read_table(os.path.join(d, 'small_calls.parquet'))
        victim = f'1:402:{b("1", 402)}:{alt("1", 402)}'
        pq.write_table(t.filter(pc.not_equal(t['key'], victim)), os.path.join(d, 'small_calls.parquet'))
        r = self.check(d)
        self.assertFalse(r['ok'])
        self.assertEqual(r['calls']['only_bcftools_unexplained'], 1)
        self.assertEqual(r['examples']['only_bcftools_unexplained'], [victim])

    def test_an_invented_call_and_a_missing_evidence_link_are_reported(self):
        d = self.copy_results('invented')
        t = pq.read_table(os.path.join(d, 'small_calls.parquet'))
        row = t.slice(0, 1).to_pylist()[0]
        row['key'] = f'1:1999:{b("1", 1999)}:{alt("1", 1999)}'
        pq.write_table(pa.concat_tables([t, pa.Table.from_pylist([row], schema=t.schema)]), os.path.join(d, 'small_calls.parquet'))
        links = pq.read_table(os.path.join(d, 'evidence_links.parquet'))
        gone = f'1:402:{b("1", 402)}:{alt("1", 402)}'
        pq.write_table(links.filter(pc.not_equal(links['call_id'], gone)), os.path.join(d, 'evidence_links.parquet'))
        r = self.check(d)
        self.assertFalse(r['ok'])
        self.assertEqual(r['calls']['only_matcher'], 1)
        self.assertEqual(r['evidence']['only_bcftools_unexplained'], 1)
        self.assertEqual(r['examples']['evidence_only_bcftools'], [gone])


    # ---- the ClinVar rows of the index against the ClinVar VCF, read by bcftools ----------------------------------------------------------
    def fields(self, idx=None):
        return crosscheck.clinvar_fields(idx or self.idx, self.cv, self.ref)

    def tampered_index(self, name, change, partition='1'):
        d = os.path.join(self.work, name)
        shutil.copytree(self.idx, d)
        path = os.path.join(d, 'evidence', 'variant', partition, 'clinvar.parquet')
        t = pq.read_table(path)
        rows = change(t.to_pylist())
        pq.write_table(pa.Table.from_pylist(rows, schema=t.schema), path)
        return d

    @staticmethod
    def pick(rows, variation, claim_type='germline_classification'):
        return next(i for i, r in enumerate(rows) if r['claim_type'] == claim_type and r['locator'].startswith('L') and r['accession'].startswith(f'VariationID:{variation}|'))

    def test_every_clinvar_vcf_row_of_the_index_equals_the_file(self):
        r = self.fields()
        self.assertTrue(r['ok'], r['examples'])
        t = r['totals']
        self.assertGreaterEqual(t['vcf_classifications'], 10)
        self.assertGreater(t['index_rows'], t['vcf_classifications'])                       # the frequencies and the record without a classification are rows too
        self.assertEqual((t['vcf_records_without_classification'] > 0, t['mapped'] > 0, t['unmapped_consistent'] > 0), (True, True, True))
        self.assertEqual(r['unmapped_by_reason'].get('ref_mismatch'), 1)                    # the record bcftools drops under --check-ref x is the one the index refused

    def test_each_field_is_compared_and_a_change_to_one_is_reported_as_that_one(self):
        def germline(i, edit):
            def change(rows):
                edit(rows[self.pick(rows, i)])
                return rows
            return change

        def frequency(edit):
            def change(rows):
                edit(next(r for r in rows if r['claim_type'] == 'population_allele_frequency' and r['claim_value'] == 'ExAC'))
                return rows
            return change
        ghost = f'1:403:{b("1", 403)}:{alt("1", 403)}'
        cases = {
            'claim_value_differs': germline(1001, lambda r: r.update(claim_value='Benign')),
            'review_status_differs': germline(1001, lambda r: r.update(review_status='no_assertion_criteria_provided')),
            'condition_text_differs': germline(1001, lambda r: r.update(condition_text='.')),                  # the real defect: the VCF's marker for no condition kept as a condition
            'condition_ids_differ': germline(1001, lambda r: r.update(condition_ids=[])),
            'n_submitters_differs': germline(1001, lambda r: r.update(n_submitters=5)),
            'review_stars_differ': germline(1001, lambda r: r.update(review_stars=None)),                       # the real defect: statuses the star table lacked
            'origin_differs': germline(1001, lambda r: r.update(origin='somatic')),
            'gene_count_differs': germline(1001, lambda r: r.update(gene_ids=[])),
            'extra_json_differs': germline(1001, lambda r: r.update(extra_json='{}')),
            'key_differs': germline(1001, lambda r: r.update(subject_id=ghost)),
            'unmapped_but_bcftools_accepts': germline(1001, lambda r: r.update(mapping_status='unmapped', mapping_reason='contig_unmapped')),
            'mapped_but_bcftools_rejects': germline(1008, lambda r: r.update(mapping_status='mapped', mapping_reason=None, subject_id=ghost)),
            'frequency_differs': frequency(lambda r: r.update(score=0.5)),
            'missing_in_index': lambda rows: [r for i, r in enumerate(rows) if i != self.pick(rows, 1001)],
            'extra_in_index': lambda rows: rows + [dict(rows[self.pick(rows, 1001)], accession='VariationID:777|AlleleID:888')],
        }
        for counter, change in cases.items():
            with self.subTest(counter):
                r = self.fields(self.tampered_index('t_' + counter, change))
                self.assertFalse(r['ok'])
                self.assertEqual({k: v for k, v in r['totals'].items() if k in crosscheck.DIFFERENCES and v}, {counter: 1}, r['examples'])

    def test_the_example_names_the_record_and_the_claim_type(self):
        r = self.fields(self.tampered_index('t_examples', lambda rows: [dict(x, claim_value='Benign') if i == self.pick(rows, 1001) else x for i, x in enumerate(rows)]))
        self.assertEqual(r['examples']['claim_value_differs'], ['VariationID:1001|AlleleID:2001'])
        r = self.fields(self.tampered_index('t_examples2', lambda rows: [x for i, x in enumerate(rows) if i != self.pick(rows, 1010, 'oncogenicity_classification')]))
        self.assertEqual(r['examples']['missing_in_index'], ['oncogenicity_classification:VariationID:1010|AlleleID:2010'])

    def test_one_ncbi_gene_id_must_map_to_one_gene_entity(self):
        def split(rows):
            rows[self.pick(rows, 1001)]['gene_ids'] = ['ENSG_ELSEWHERE']                              # GENEA:1001 is also on VariationID 1002 and 1003
            return rows
        r = self.fields(self.tampered_index('t_genes', split))
        self.assertEqual(r['totals']['gene_id_inconsistent'], 1)
        self.assertTrue(r['examples']['gene_id_inconsistent'][0].startswith('NCBI Gene 1001'))

    def test_records_on_a_contig_the_file_names_differently_are_renamed_before_normalising(self):
        # found on the real file: its records on RefSeq-named contigs failed because the header defines no contig and bcftools annotate --rename-chrs needs one
        r = crosscheck.clinvar_fields(self.idx, self.cv_chr, self.ref)
        base = self.fields()
        self.assertTrue(r['ok'], r['examples'])
        self.assertEqual({k: v for k, v in r['totals'].items() if k != 'vcf_records'}, {k: v for k, v in base['totals'].items() if k != 'vcf_records'})
        self.assertGreater(r['totals']['mapped'], 0)

    def test_a_file_whose_header_lacks_the_somatic_tags_reads_them_as_absent(self):
        r = crosscheck.clinvar_fields(self.idx, self.cv_old, self.ref)                       # a ClinVar release older than the oncogenicity and somatic-impact tags
        self.assertEqual({k: v for k, v in r['totals'].items() if k in crosscheck.DIFFERENCES and v}, {'extra_in_index': 2, 'missing_in_index': 1})       # ONC and SCI rows; the record that now has no classification
        self.assertEqual(sorted(r['examples']['extra_in_index']), ['oncogenicity_classification:VariationID:1010|AlleleID:2010', 'somatic_clinical_impact:VariationID:1010|AlleleID:2010'])


if __name__ == '__main__':
    unittest.main()
