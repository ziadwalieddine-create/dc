import glob
import json
import os
import shutil
import sys
import tempfile
import unittest

import pyarrow as pa
import pyarrow.parquet as pq

from helpers import syn_reference, FIXTURES
from test_build_gwas import variants as gwas_variants
from vindex import provenance
from vindex.build import runner
from vindex.build.clinvar_vcf import ClinVarVcfAdapter
from vindex.build.clinvar_vcv import ClinVarVcvAdapter
from vindex.build.genes import GeneMap
from vindex.build.gnomad import GnomadConstraintAdapter, GnomadFrequencyAdapter, GnomadLofCurationAdapter
from vindex.build.gwas import GwasCatalogAdapter

sys.path.insert(0, FIXTURES)
import make_sources  # noqa: E402

SOURCES = make_sources.write_all(tempfile.mkdtemp())


class FakeGenes:
    by_ensg = {'ENSG_G1': 1, 'ENSG_G2': 1}
    rows = []

    def from_ncbi(self, x):
        return None


def make_world():
    """A private copy of the mini sources, an index built from them, and a manifest/targets pair describing exactly those files."""
    work = tempfile.mkdtemp()
    src = os.path.join(work, 'src')
    shutil.copytree(SOURCES, src)
    ref = syn_reference()
    genes = GeneMap(ref, ot_dir='/nonexistent', mane_path='/nonexistent')
    vcf = os.path.join(src, 'mini_clinvar.vcf')
    gn = os.path.join(src, 'gnomad')

    class Vcf(ClinVarVcfAdapter):
        def run(self, ctx, limit=None):
            return super().run(ctx, path=vcf, limit=limit)
    idx = os.path.join(work, 'idx')
    runner.build(idx, [Vcf(genes), ClinVarVcvAdapter(genes, root=os.path.join(src, 'vcv')), GwasCatalogAdapter(FakeGenes(), gwas_variants(), root=src),
                       GnomadConstraintAdapter(FakeGenes(), os.path.join(gn, 'gnomad.v4.1.2.constraint_metrics.tsv.bgz')),
                       GnomadLofCurationAdapter(FakeGenes(), os.path.join(gn, 'lof_curation', 'incomplete_penetrance_curation_results.csv')),
                       GnomadFrequencyAdapter(os.path.join(gn, 'regional', 'gnomad_v4_1_joint_at_mavedb_positions.tsv.gz'), os.path.join(gn, 'regional', 'gnomad_v4_1_exomes_at_mavedb_positions.tsv.gz'))], ref)
    paths = {'clinvar_vcf_grch38': vcf, 'clinvar_vcv_xml_weekly_2026_1018': os.path.join(src, 'vcv'), 'gwas_catalog_associations': os.path.join(src, 'gwas-catalog-associations_ontology-annotated-full.zip'),
             'gnomad_constraint_v4_1_2': os.path.join(gn, 'gnomad.v4.1.2.constraint_metrics.tsv.bgz'),
             'gnomad_v4_1_lof_curation_incomplete_penetrance': os.path.join(gn, 'lof_curation', 'incomplete_penetrance_curation_results.csv'),
             'gnomad_v4_1_joint_regional': os.path.join(gn, 'regional', 'gnomad_v4_1_joint_at_mavedb_positions.tsv.gz'),
             'gnomad_v4_1_exomes_regional': os.path.join(gn, 'regional', 'gnomad_v4_1_exomes_at_mavedb_positions.tsv.gz')}
    rels = provenance.source_releases(idx)
    man, tgt = [], []
    for r in rels.values():
        p = paths.get(r['target_id'])
        man.append(dict(id=r['target_id'], sha256=r['sha256'], bytes=os.path.getsize(p) if p and os.path.isfile(p) else None, dest=p or os.path.join(src, 'unused', r['target_id']), url='https://x', licence='open'))
        tgt.append(dict(id=r['target_id'], dest=p or os.path.join(src, 'unused', r['target_id'])))
        if not (p and os.path.exists(p)):
            os.makedirs(os.path.join(src, 'unused', r['target_id']), exist_ok=True)           # releases with no evidence rows to re-open still need their file to exist
    mp, tp = os.path.join(work, 'manifest.jsonl'), os.path.join(work, 'targets.json')
    with open(mp, 'w') as f:
        f.write('\n'.join(json.dumps(m) for m in man) + '\n')
    with open(tp, 'w') as f:
        json.dump(tgt, f)
    return dict(work=work, idx=idx, paths=paths, manifest=mp, targets=tp)


def verify(w, **kw):
    return provenance.verify_index(w['idx'], sample=1000, paths=w['paths'], manifest_path=w['manifest'], targets_path=w['targets'], data_root=w['work'], **kw)


class TestProvenance(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.w = make_world()
        cls.report = verify(cls.w)

    @classmethod
    def tearDownClass(cls):
        shutil.rmtree(cls.w['work'], ignore_errors=True)

    def test_a_clean_index_verifies_and_every_resolver_family_is_exercised(self):
        r = self.report
        self.assertTrue(r['ok'], json.dumps({k: v for k, v in r.items() if k in ('orphan_release_ids', 'release_issues')}) + json.dumps({k: (v['n_mismatched'], v['n_unresolved']) for k, v in r['digests'].items()}))
        by_target = {d['target_id']: d for d in r['digests'].values()}
        self.assertEqual(set(by_target), {'clinvar_vcf_grch38', 'clinvar_vcv_xml_weekly_2026_1018', 'gwas_catalog_associations', 'gnomad_constraint_v4_1_2',
                                          'gnomad_v4_1_lof_curation_incomplete_penetrance', 'gnomad_v4_1_joint_regional', 'gnomad_v4_1_exomes_regional'})
        for t, d in by_target.items():
            self.assertGreater(d['sampled'], 0, t)
            self.assertEqual((d['verified'], d['n_mismatched'], d['n_unresolved']), (d['sampled'], 0, 0), t)
        self.assertEqual({d['resolver'] for d in by_target.values()}, {'iter', 'iter_raw', 'splitlines', 'gwas', 'parquet'})
        self.assertEqual(r['no_resolver'], [])

    def test_the_report_names_the_code_that_made_it_and_the_sample_it_asked_for(self):
        import hashlib
        here = os.path.dirname(os.path.abspath(provenance.__file__))
        for f, digest in self.report['code'].items():
            with open(os.path.join(here, f), 'rb') as fh:
                self.assertEqual(digest, hashlib.sha256(fh.read()).hexdigest(), f)
        self.assertEqual((sorted(self.report['code']), self.report['sample'], self.report['seed']), (['provenance.py', 'resolvers.py'], 1000, 0))

    def test_parquet_rows_are_resolved_part_by_part_and_out_of_range_rows_stay_unresolved(self):
        root = tempfile.mkdtemp()
        try:
            os.makedirs(os.path.join(root, 'scv'))
            for k in range(3):
                pq.write_table(pa.table({'scv_id': [f'SCV{k}{i}' for i in range(4)], 'n': list(range(4))}), os.path.join(root, 'scv', f'part-000{k}.parquet'))
            got = provenance.resolve_parquet_rows(root, ['scv/part-0000.parquet#1', 'scv/part-0002.parquet#3', 'scv/part-0001.parquet#0', 'scv/part-0001.parquet#9', 'scv/part-0007.parquet#0', 'bad'])
            self.assertEqual(sorted(got), ['scv/part-0000.parquet#1', 'scv/part-0001.parquet#0', 'scv/part-0002.parquet#3'])
            self.assertEqual(got['scv/part-0002.parquet#3'], provenance.vcv_canon({'scv_id': 'SCV23', 'n': 3}))
        finally:
            shutil.rmtree(root, ignore_errors=True)

    def test_sampling_is_deterministic_and_bounded(self):
        small = provenance.verify_index(self.w['idx'], sample=3, paths=self.w['paths'], manifest_path=self.w['manifest'], targets_path=self.w['targets'], data_root=self.w['work'])
        again = provenance.verify_index(self.w['idx'], sample=3, paths=self.w['paths'], manifest_path=self.w['manifest'], targets_path=self.w['targets'], data_root=self.w['work'])
        self.assertEqual(small['digests'], again['digests'])
        self.assertTrue(all(d['sampled'] <= 3 for d in small['digests'].values()))
        self.assertTrue(small['ok'])

    def test_the_records_are_re_opened_in_chunks_and_the_chunk_size_changes_nothing(self):
        saved = provenance.CHUNK
        try:
            provenance.CHUNK = 2
            small = verify(self.w)
        finally:
            provenance.CHUNK = saved
        self.assertEqual(small['digests'], self.report['digests'])
        self.assertTrue(any(d['sampled'] > 2 for d in small['digests'].values()))                        # more than one chunk was needed

    def test_no_resolver_is_asked_for_more_records_than_one_chunk(self):
        sizes, saved = [], (provenance.resolve_lines, provenance.resolve_parquet_rows, provenance.CHUNK)

        def spy(fn):
            def wrapper(*args, **kw):
                sizes.append(len(args[-1]))
                return fn(*args, **kw)
            return wrapper
        try:
            provenance.resolve_lines, provenance.resolve_parquet_rows, provenance.CHUNK = spy(saved[0]), spy(saved[1]), 2
            verify(self.w)                                                                                  # every row of each release
            every = list(sizes)
            del sizes[:]
            provenance.verify_index(self.w['idx'], sample=3, paths=self.w['paths'], manifest_path=self.w['manifest'], targets_path=self.w['targets'], data_root=self.w['work'])    # a sample
        finally:
            provenance.resolve_lines, provenance.resolve_parquet_rows, provenance.CHUNK = saved
        self.assertTrue(every and max(every) <= 2 and 2 in every, every)
        self.assertTrue(sizes and max(sizes) <= 2 and 2 in sizes and 1 in sizes, sizes)                      # three sampled rows are asked for as two and one

    def test_every_row_is_re_opened_when_the_sample_is_at_least_the_release_and_a_smaller_sample_is_exactly_that_size(self):
        full = provenance.verify_index(self.w['idx'], sample=10 ** 7, paths=self.w['paths'], manifest_path=self.w['manifest'], targets_path=self.w['targets'], data_root=self.w['work'])
        part = provenance.verify_index(self.w['idx'], sample=3, paths=self.w['paths'], manifest_path=self.w['manifest'], targets_path=self.w['targets'], data_root=self.w['work'])
        for rid, d in full['digests'].items():
            self.assertEqual((d['sampled'], d['verified']), (d['evidence_rows'], d['evidence_rows']), rid)
            self.assertEqual(part['digests'][rid]['sampled'], min(3, d['evidence_rows']), rid)

    def test_every_damaged_record_is_counted_when_only_the_first_few_are_listed(self):
        def edit(w):
            p = w['paths']['clinvar_vcf_grch38']
            with open(p) as f:
                text = f.read()
            with open(p, 'w') as f:
                f.write(text.replace('Test_syndrome', 'Test_syndrome_EDITED'))
        saved = provenance.MAX_EXAMPLES, provenance.CHUNK
        try:
            provenance.MAX_EXAMPLES, provenance.CHUNK = 1, 2
            r = self.damaged(edit)
        finally:
            provenance.MAX_EXAMPLES, provenance.CHUNK = saved
        d = next(d for d in r['digests'].values() if d['target_id'] == 'clinvar_vcf_grch38')
        self.assertGreater(d['n_mismatched'], 1)                                                         # several records carry the edited text, in more than one chunk
        self.assertEqual(len(d['mismatched']), 1)

    # ---- adversarial: each kind of damage is caught ---------------------------------------------------------------------------------------------
    def damaged(self, mutate):
        w = make_world()
        try:
            mutate(w)
            return verify(w)
        finally:
            shutil.rmtree(w['work'], ignore_errors=True)

    def test_an_edited_source_record_is_caught(self):
        def edit(w):
            p = w['paths']['clinvar_vcf_grch38']
            with open(p) as f:
                text = f.read()
            with open(p, 'w') as f:
                f.write(text.replace('Test_syndrome', 'Test_syndrome_EDITED', 1))
        r = self.damaged(edit)
        d = next(d for d in r['digests'].values() if d['target_id'] == 'clinvar_vcf_grch38')
        self.assertGreater(d['n_mismatched'], 0)
        self.assertFalse(r['ok'])

    def test_a_missing_source_file_leaves_rows_unresolved_not_verified(self):
        def drop(w):
            os.remove(w['paths']['gnomad_constraint_v4_1_2'])
        r = self.damaged(drop)
        d = next(d for d in r['digests'].values() if d['target_id'] == 'gnomad_constraint_v4_1_2')
        self.assertEqual((d['verified'], d['n_unresolved'] > 0), (0, True))
        self.assertTrue(any(i['issue'] == 'acquired_file_missing' for i in r['release_issues']))
        self.assertFalse(r['ok'])

    def test_a_tampered_digest_in_the_index_is_caught(self):
        def tamper(w):
            f = sorted(glob.glob(os.path.join(w['idx'], 'evidence', 'gene', '*', 'gnomad_constraint.parquet')))[0]
            t = pq.read_table(f)
            col = t['raw_sha256'].to_pylist()
            col[0] = '0' * 64
            pq.write_table(t.set_column(t.schema.get_field_index('raw_sha256'), 'raw_sha256', pa.array(col, pa.string())), f)
        r = self.damaged(tamper)
        d = next(d for d in r['digests'].values() if d['target_id'] == 'gnomad_constraint_v4_1_2')
        self.assertEqual(d['n_mismatched'], 1)

    def test_evidence_from_an_unnamed_release_is_caught(self):
        def drop_release(w):
            p = os.path.join(w['idx'], 'source_release.parquet')
            t = pq.read_table(p)
            keep = [i for i, r in enumerate(t['target_id'].to_pylist()) if r != 'gnomad_constraint_v4_1_2']
            pq.write_table(t.take(keep), p)
        r = self.damaged(drop_release)
        self.assertIn('evidence', r['orphan_release_ids'])
        self.assertTrue(any('gnomad_constraint_v4_1_2' in x for x in r['orphan_release_ids']['evidence']))
        self.assertFalse(r['ok'])

    def test_manifest_disagreement_and_size_change_are_caught(self):
        def edit(w):
            with open(w['manifest']) as f:
                recs = [json.loads(l) for l in f if l.strip()]
            for m in recs:
                if m['id'] == 'clinvar_vcf_grch38':
                    m['sha256'] = 'f' * 64
            with open(w['manifest'], 'w') as f:
                f.write('\n'.join(json.dumps(m) for m in recs) + '\n')
            with open(w['paths']['gnomad_constraint_v4_1_2'], 'ab') as f:
                f.write(b'x')
        r = self.damaged(edit)
        issues = {(i['release_id'].split('@')[0].split('/')[1], i['issue']) for i in r['release_issues']}
        self.assertIn(('clinvar_vcf_grch38', 'release_digest_differs_from_manifest'), issues)
        self.assertIn(('gnomad_constraint_v4_1_2', 'acquired_file_size_differs_from_manifest'), issues)

    def test_releases_without_a_resolver_are_listed_not_counted_as_verified(self):
        w = make_world()
        try:
            saved = dict(provenance.RESOLVERS)
            del provenance.RESOLVERS['gnomad_v4_1_joint_regional']
            try:
                r = verify(w)
            finally:
                provenance.RESOLVERS.update(saved)
        finally:
            shutil.rmtree(w['work'], ignore_errors=True)
        self.assertEqual([x['target_id'] for x in r['no_resolver']], ['gnomad_v4_1_joint_regional'])
        self.assertNotIn('gnomad_v4_1_joint_regional', {d['target_id'] for d in r['digests'].values()})


class TestDerivedTargets(unittest.TestCase):
    """A target written from another acquired file has no URL of its own; its retrieval URL is that of the file it was extracted from."""

    def issues(self, targets):
        with tempfile.TemporaryDirectory() as d:
            with open(os.path.join(d, 'derived'), 'w') as f:
                f.write('x')
            rel = {'mavedb/derived@v1': dict(release_id='mavedb/derived@v1', target_id='derived', sha256='a' * 64, url='', licence='CC0')}
            manifest = {'derived': dict(sha256='a' * 64, bytes=1, dest='derived')}
            return [i['issue'] for i in provenance.check_releases(rel, manifest, targets, d)]

    def test_the_url_comes_from_the_parent_and_only_from_a_parent_that_has_one(self):
        parent = {'parent': dict(url='https://example.org/dump.tar.gz')}
        self.assertEqual(self.issues({'derived': dict(licence='CC0', derived_from='parent'), **parent}), [])
        self.assertEqual(self.issues({'derived': dict(licence='CC0'), **parent}), ['no_url'])                      # nothing says where it came from
        self.assertEqual(self.issues({'derived': dict(licence='CC0', derived_from='missing')}), ['no_url'])        # the parent is not a target
        self.assertEqual(self.issues({'derived': dict(licence='CC0', derived_from='parent'), 'parent': dict(url=None)}), ['no_url'])


if __name__ == '__main__':
    unittest.main()
