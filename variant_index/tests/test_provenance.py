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


if __name__ == '__main__':
    unittest.main()
