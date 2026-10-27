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
from vindex import reports, schema


class TestReports(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.work = tempfile.mkdtemp()
        cls.idx = os.path.join(cls.work, 'idx')
        build_index(cls.idx)
        cls.cov = reports.coverage(cls.idx)
        cls.clean = reports.validate(cls.idx, syn_reference())

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

    def test_lineage_groups_shared_between_sources_are_counted(self):
        # the fixture index has no group two sources write: the figure is 0 here and the pairs are empty ...
        lin = self.cov['lineage']
        self.assertGreater(lin['groups'], 5)
        self.assertEqual((lin['shared_between_sources'], lin['source_pairs']), (0, {}))
        # ... and a group written by two sources is found, with the pair that shares it
        d = os.path.join(self.work, 'shared')
        shutil.copytree(self.idx, d)
        for src, group in (('cpic', 'clinvar:variation:1001'), ('strchive', 'clinvar:variation:1001')):
            for p in glob.glob(os.path.join(d, 'evidence', '*', '*', f'{src}.parquet')):
                t = pq.read_table(p)
                rows = t.to_pylist()
                rows[0]['lineage_group'] = group
                pq.write_table(pa.Table.from_pylist(rows, schema=t.schema), p)
                break
        for p in glob.glob(os.path.join(d, 'evidence', '*', '*', 'clinvar.parquet')):
            t = pq.read_table(p)
            rows = t.to_pylist()
            rows[0]['lineage_group'] = 'clinvar:variation:1001'
            pq.write_table(pa.Table.from_pylist(rows, schema=t.schema), p)
            break
        lin = reports.lineage_sharing(d)
        self.assertEqual(lin['shared_between_sources'], 1)
        self.assertEqual(lin['source_pairs'], {'clinvar|cpic': 1, 'clinvar|strchive': 1, 'cpic|strchive': 1})

    def test_dbsnp_identifiers_placed_at_two_sites_are_counted_with_the_sources_that_disagree(self):
        x = self.cov['xref_coherence']
        self.assertGreater(x['identifiers'], 0)
        self.assertEqual((x['multi_position'], x['multi_position_across_sources'], x['source_pairs']), (0, 0, {}))              # the fixture sources agree on every rsID
        d = os.path.join(self.work, 'xrefs')
        shutil.copytree(self.idx, d)
        def row(rs, key, relation):
            return dict(namespace='dbsnp', identifier=rs, entity_level='variant', entity_id=key, release_id='r', relation=relation)
        rows = [row('rs900', '1:100:A:G', 'clinvar_RS'), row('rs900', '1:100:A:T', 'gwas_catalog_snps'),                       # one site, two alleles: not a disagreement
                row('rs901', '1:100:A:G', 'clinvar_RS'), row('rs901', '2:555:C:T', 'gwas_catalog_snps'),                       # two sources, two sites
                row('rs902', '3:7:A:G', 'clinvar_RS'), row('rs902', '3:9:A:G', 'clinvar_RS'),                                   # one source, two sites
                row('rs903', '4:1:A:G', 'pharmvar_rsID'), row('rs903', '5:1:A:G', 'clinvar_RS'), row('rs903', '4:1:A:C', 'gwas_catalog_snps')]
        pq.write_table(pa.Table.from_pylist(rows, schema=schema.XREF), os.path.join(d, 'xref', 'zz_disagree.parquet'))
        x2 = reports.xref_coherence(d)
        self.assertEqual((x2['identifiers'] - x['identifiers'], x2['multi_position'], x2['multi_position_across_sources'], x2['multi_position_within_one_source']), (4, 3, 2, 1))
        self.assertEqual(x2['source_pairs'], {'clinvar_RS|gwas_catalog_snps': 2, 'clinvar_RS|pharmvar_rsID': 1})                  # PharmVar and the GWAS Catalog agree on rs903
        self.assertEqual([e['identifier'] for e in x2['examples']], ['rs901', 'rs903'])
        self.assertEqual(x2['examples'][0]['placed'], {'clinvar_RS': '1:100', 'gwas_catalog_snps': '2:555'})

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

    def add_rows(self, d, table, rows):
        os.makedirs(os.path.join(d, table), exist_ok=True)
        pq.write_table(schema.rows_to_table(table, rows), os.path.join(d, table, 'zz_added.parquet'))

    def test_dangling_references_between_tables_are_counted(self):
        def mutate(d):
            self.add_rows(d, 'xref', [dict(namespace='dbsnp', identifier='rs1', entity_level=level, entity_id=f'{level}:nowhere', release_id='r', relation='x')
                                      for level in ('variant', 'event', 'haplotype', 'gene')])
            self.add_rows(d, 'haplotype_definition', [dict(haplotype_id='nowhere', variant_key='1:1:A:C', hgvs=None, release_id='r', definition_scope='defining')])
            self.add_rows(d, 'variant_gene', [dict(variant_key='1:99999:A:C', gene_id='ENSG_NOT_HELD', relation='overlaps', consequence=None, release_id='r')])
            self.add_rows(d, 'gene_set_member', [dict(set_id='s', source_id='reactome', kind='pathway', name='p', gene_id='ENSG_NOT_HELD', release_id='r')])
        v = self.damaged(mutate)
        for name in ('xref_variant_not_an_entity', 'xref_event_not_an_entity', 'xref_haplotype_not_an_entity', 'xref_gene_not_an_entity', 'haplotype_definition_haplotype_not_an_entity',
                     'variant_gene_variant_not_an_entity'):
            self.assertEqual(v['checks'][name]['violations'], 1, name)
        self.assertFalse(v['ok'])
        self.assertEqual(v['informational']['variant_gene_gene_not_in_the_gene_table'] - self.clean['informational']['variant_gene_gene_not_in_the_gene_table'], 1)
        self.assertEqual(v['informational']['gene_set_member_gene_not_in_the_gene_table'] - self.clean['informational']['gene_set_member_gene_not_in_the_gene_table'], 1)
        for name in ('variant_gene_gene_not_in_the_gene_table', 'gene_set_member_gene_not_in_the_gene_table'):
            self.assertNotIn(name, v['checks'])                                                            # reported, not failing

    def test_entity_invariants_and_marker_text_in_entity_attributes_are_caught(self):
        def mutate(d):
            self.add_rows(d, 'entity_study', [dict(study_id='GCST1', source_id='gwas_catalog', lineage_group='l', study_type='gwas', pubmed='1', trait_ids=None, trait_text='NR', n_total=100, n_cases=200,
                                                   ancestry='NR', cohorts=None, biosample=None, has_sumstats=False, release_id='r'),
                                              dict(study_id='GCST2', source_id='gwas_catalog', lineage_group='l', study_type='gwas', pubmed='1', trait_ids=None, trait_text='x', n_total=0, n_cases=None,
                                                   ancestry='', cohorts=None, biosample=None, has_sumstats=False, release_id='r')])
            self.add_rows(d, 'entity_event', [dict(event_id='e_bad', source_id='x', release_id='r', accession='a', event_type='DEL', chrom='1', start=500, end=100, population_af=1.5, copy_number=-1)])
            self.add_rows(d, 'entity_gene', [dict(gene_id='ENSG_BAD', symbol='B', hgnc_id='', chrom='1', start=900, end=100, strand='?', mane_select='', biotype='protein_coding', release_id='r')])
            self.add_rows(d, 'entity_variant', [dict(variant_key='1:0:A:C', chrom='1', pos=0, ref='A', alt='C', vtype='snv'), dict(variant_key='1:5:A:G', chrom='1', pos=5, ref='A', alt='C', vtype='snv')])
        v = self.damaged(mutate)
        expected = {'entity_study_cases_exceed_total': 1, 'entity_study_total_not_positive': 1, 'entity_study_ancestry_is_a_missing_value_marker': 1, 'entity_study_trait_text_is_a_missing_value_marker': 1,
                    'entity_event_start_after_end': 1, 'entity_event_frequency_outside_0_1': 1, 'entity_event_negative_copy_number': 1, 'entity_gene_start_after_end': 1,
                    'entity_gene_strand_not_plus_or_minus': 1, 'entity_variant_position_not_positive': 1, 'entity_variant_key_differs_from_its_parts': 1}
        got = {k: v['checks'][k]['violations'] for k in expected}
        self.assertEqual(got, expected)
        self.assertEqual(v['checks']['entity_gene_hgnc_id_is_a_missing_value_marker']['violations'], 0)          # an empty string is an empty attribute, not a marker

    def test_vocabulary_range_and_format_violations_are_caught(self):
        v = self.damaged(lambda d: self.rewrite(d, 'evidence/variant/1/clinvar.parquet', 'effect_assertion', lambda xs: ['bogus'] + xs[1:]))
        self.assertEqual(v['checks']['effect_assertion_outside_vocabulary']['violations'], 1)
        self.assertFalse(v['ok'])
        v = self.damaged(lambda d: self.rewrite(d, 'evidence/variant/1/clinvar.parquet', 'p_value', lambda xs: [1.5] + xs[1:]))
        self.assertEqual(v['checks']['p_value_outside_0_1']['violations'], 1)
        v = self.damaged(lambda d: self.rewrite(d, 'evidence/variant/1/clinvar.parquet', 'raw_sha256', lambda xs: ['xyz'] + xs[1:]))
        self.assertEqual(v['checks']['raw_sha256_not_64_hex']['violations'], 1)

    def test_a_missing_value_marker_kept_as_a_value_and_impossible_numbers_are_caught(self):
        # each of these was found on the real index by comparing it with its sources
        p = 'evidence/variant/1/clinvar.parquet'
        v = self.damaged(lambda d: self.rewrite(d, p, 'condition_text', lambda xs: ['.'] + xs[1:]))
        self.assertEqual(v['checks']['condition_text_is_a_missing_value_marker']['violations'], 1)             # CLNDN=. kept as a condition
        v = self.damaged(lambda d: self.rewrite(d, p, 'ancestry', lambda xs: ['NR'] + xs[1:]))
        self.assertEqual(v['checks']['ancestry_is_a_missing_value_marker']['violations'], 1)                   # 'not reported' kept as an ancestry
        v = self.damaged(lambda d: self.rewrite(d, p, 'review_status', lambda xs: ['n/a'] + xs[1:]))
        self.assertEqual(v['checks']['review_status_is_a_missing_value_marker']['violations'], 1)
        v = self.damaged(lambda d: self.rewrite(d, p, 'eaf', lambda xs: [1.5] + xs[1:]))
        self.assertEqual(v['checks']['eaf_outside_0_1']['violations'], 1)                                      # 74 Catalog cells above 1
        v = self.damaged(lambda d: self.rewrite(d, p, 'se', lambda xs: [0.0] + xs[1:]))
        self.assertEqual(v['checks']['standard_error_not_positive']['violations'], 1)                          # MaveDB se = 0
        v = self.damaged(lambda d: self.rewrite(d, p, 'beta', lambda xs: [float('nan')] + xs[1:]))
        self.assertEqual(v['checks']['beta_not_finite']['violations'], 1)                                      # float('nan') and float('inf') parse
        v = self.damaged(lambda d: (self.rewrite(d, p, 'ci_low', lambda xs: [3.0] + xs[1:]), self.rewrite(d, p, 'ci_high', lambda xs: [2.0] + xs[1:])))
        self.assertEqual(v['checks']['ci_low_above_ci_high']['violations'], 1)                                 # 1,746 inverted betas

    def test_an_estimate_beside_an_interval_that_does_not_contain_it_is_caught_but_a_rounding_slip_is_not(self):
        p = 'evidence/variant/1/clinvar.parquet'

        def set_rows(d, beta, lo, hi, n):
            for col, val in (('beta', beta), ('ci_low', lo), ('ci_high', hi)):
                self.rewrite(d, p, col, lambda xs, val=val: [val] * n + xs[n:])
        v = self.damaged(lambda d: set_rows(d, -0.032, 0.018, 0.046, 3))                                       # the signed beta beside the bounds of its magnitude
        c = v['checks']['estimate_outside_its_interval']
        self.assertEqual(c['violations'], 3)
        self.assertEqual(c['detail']['outside'], 3)
        v = self.damaged(lambda d: set_rows(d, 0.0321, 0.018, 0.032, 3))                                       # within the 1% slack of the interval's width
        self.assertEqual(v['checks']['estimate_outside_its_interval']['violations'], 0)

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

    def test_an_effect_allele_that_is_not_an_allele_of_the_snv_key_is_caught(self):
        clean = reports.validate(self.idx, syn_reference())
        self.assertEqual((clean['checks']['effect_allele_not_an_allele_of_the_snv_key']['violations'], clean['checks']['other_allele_not_an_allele_of_the_snv_key']['violations']), (0, 0))
        v = self.damaged(lambda d: self.rewrite(d, 'evidence/variant/1/clinvar.parquet', 'effect_allele', lambda xs: ['Z'] * len(xs)))
        self.assertGreater(v['checks']['effect_allele_not_an_allele_of_the_snv_key']['violations'], 0)
        self.assertFalse(v['ok'])
        v = self.damaged(lambda d: self.rewrite(d, 'evidence/variant/1/clinvar.parquet', 'other_allele', lambda xs: ['Z'] * len(xs)))
        self.assertGreater(v['checks']['other_allele_not_an_allele_of_the_snv_key']['violations'], 0)

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
        self.assertFalse(c['build_code']['changed'])
        m['build_code'] = dict(m['build_code'], sha256_16='0' * 16)                         # a shared helper changed, no adapter file did
        with open(os.path.join(other, 'index_manifest.json'), 'w') as f:
            json.dump(m, f)
        self.assertTrue(reports.change_log(self.idx, other)['build_code']['changed'])
        self.assertEqual(c['tables']['evidence']['delta'], 5)
        self.assertEqual(reports.change_log(self.idx, self.idx)['tables'], {})

    def test_the_markdown_build_report_names_every_section_and_is_free_of_sample_data(self):
        import make_reports
        tp, mp, vp = (os.path.join(self.work, n) for n in ('t2.json', 'm2.jsonl', 'v2.json'))
        rels = sorted({r['target_id'] for r in pq.read_table(os.path.join(self.idx, 'source_release.parquet')).to_pylist()})
        with open(tp, 'w') as f:
            json.dump([dict(id=i, source='clinvar', dest='x') for i in rels + ['not_downloaded']], f)
        with open(mp, 'w') as f:
            f.write('\n'.join(json.dumps(dict(id=i, ok=True)) for i in rels) + '\n')
        with open(vp, 'w') as f:
            json.dump({i: dict(passed=True) for i in rels}, f)
        life = reports.lifecycle(self.idx, tp, mp, vp, catalogue_dir=os.path.join(self.work, 'none'))
        cov = reports.coverage(self.idx)
        gaps = reports.gap_report(self.idx, life, cov)
        with open(os.path.join(self.idx, 'index_manifest.json')) as f:
            manifest = json.load(f)
        text = make_reports.md(cov, reports.validate(self.idx, syn_reference()), life, gaps, None, manifest)
        for section in ('## Table sizes', '## Evidence rows', '## dbSNP identifiers placed at more than one site', '## Validation', '## Target lifecycle', 'Listed but not acquired', '`not_downloaded`', '## Records kept but not mapped, by source', '## Modules and gaps'):
            self.assertIn(section, text)
        self.assertIn('passed', text)
        self.assertIn('| `' + rels[0] + '` |', text)

    def test_an_index_built_under_another_schema_version_is_refused_unless_knowingly_opened(self):
        from vindex.index_reader import Index, SchemaMismatch
        other = os.path.join(self.work, 'idx_other_schema')
        shutil.copytree(self.idx, other)
        with open(os.path.join(other, 'index_manifest.json')) as f:
            m = json.load(f)
        m['schema_version'] = '0.0.1'
        with open(os.path.join(other, 'index_manifest.json'), 'w') as f:
            json.dump(m, f)
        with self.assertRaises(SchemaMismatch):
            Index.open(other)
        self.assertGreater(Index.open(other, allow_schema_mismatch=True).table('entity_gene').num_rows, 0)
        os.remove(os.path.join(other, 'index_manifest.json'))
        with self.assertRaises(SchemaMismatch):
            Index.open(other)
        Index.open(self.idx)                                                                         # the matching index opens

    def test_the_module_registry_names_what_is_not_implemented(self):
        status = {m['name']: m['status'] for m in reports.MODULES}
        self.assertEqual(status['personal 30x GRCh38 SNP/indel, CNV and SV calls'], 'unavailable_input')
        self.assertEqual(status['transcript HGVS (c., n.) to genomic coordinate mapping'], 'unimplemented')
        self.assertTrue(all(m['status'] in ('implemented', 'partial', 'unimplemented', 'external', 'not_acquired', 'unavailable_input') for m in reports.MODULES))
        self.assertTrue(all(m.get('note') for m in reports.MODULES if m['status'] != 'implemented'))


if __name__ == '__main__':
    unittest.main()
