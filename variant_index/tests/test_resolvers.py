"""Every provenance resolver reproduces, for EVERY evidence row of its release, the digest the adapter stored: the releases are built by the real adapters from the mini sources and each
row is re-opened by locator in the file it came from. A resolver that disagrees with its adapter fails here, and a row it cannot find is reported as unresolved, never as verified."""
import functools
import os
import shutil
import tempfile
import unittest

from helpers import SYN_BUILDS, SYN_FASTA, syn_reference
import test_build_clingen as tcg
import test_build_clinpgx as tcp
import test_build_mavedb as tmv
import test_build_opentargets as tot
import test_build_opentargets_gene as tog
import test_build_pgx as tpg
import test_build_regions as trg
from vindex import provenance, resolvers
from vindex.build import runner
from vindex.build.genes import GeneMap
from vindex.build.regions import StrchiveAdapter
from vindex.refgenome import Reference


class ResolverCase(unittest.TestCase):
    """build() an index with `adapters`, then verify every row of each release whose target is in `paths`."""

    def verify(self, adapters, ref, paths, override=None):
        out = tempfile.mkdtemp()
        try:
            runner.build(out, adapters, ref)
            saved = dict(provenance.RESOLVERS)
            provenance.RESOLVERS.update(override or {})
            try:
                results = {}
                for rid, rel in provenance.source_releases(out).items():
                    tid = rel['target_id']
                    if tid not in paths:
                        continue
                    d = provenance.verify_release(out, rid, tid, paths[tid], sample=10 ** 7)
                    if d is not None:
                        results[tid] = d
                return results
            finally:
                provenance.RESOLVERS.clear()
                provenance.RESOLVERS.update(saved)
        finally:
            shutil.rmtree(out, ignore_errors=True)

    def assertAllVerified(self, results, targets):
        self.assertEqual(sorted(results), sorted(targets))
        for tid, d in results.items():
            self.assertGreater(d['evidence_rows'], 0, tid)
            self.assertEqual((d['verified'], d['n_mismatched'], d['n_unresolved']), (d['evidence_rows'], 0, 0), (tid, d))


class TestClinGen(ResolverCase):
    def test_validity_dosage_and_actionability_records_are_recovered(self):
        S = tcg.SOURCES
        ac = lambda name: os.path.join(S, 'actionability', name)
        paths = {'clingen_gene_validity': os.path.join(S, 'clingen_gene_disease_validity.csv'), 'clingen_gene_curation_list': os.path.join(S, 'ClinGen_gene_curation_list_GRCh38.tsv'),
                 'clingen_region_curation_list': os.path.join(S, 'ClinGen_region_curation_list_GRCh38.tsv'), 'clingen_erepo_variant_classifications': os.path.join(S, 'erepo_variant_classifications.tsv'),
                 'clingen_actionability_adult_assertions': ac('adult_assertions_all.tsv'), 'clingen_actionability_adult_scoring': ac('adult_scoring.tsv'),
                 'clingen_actionability_pediatric_assertions': ac('pediatric_assertions_all.tsv'), 'clingen_actionability_pediatric_scoring': ac('pediatric_scoring.tsv')}
        r = self.verify([tcg.MiniClinVar(), tcg.MiniErepo(), tcg.MiniGeneRegion()], syn_reference(), paths)
        self.assertIn('clingen_gene_validity', r)
        self.assertIn('clingen_gene_curation_list', r)
        self.assertAllVerified(r, list(r))
        self.assertGreaterEqual(len(r), 6)

    def test_a_changed_record_is_a_mismatch_not_a_pass(self):
        S = tcg.SOURCES
        d = tempfile.mkdtemp()
        try:
            edited = os.path.join(d, 'validity.csv')
            with open(os.path.join(S, 'clingen_gene_disease_validity.csv'), encoding='utf-8') as f:
                text = f.read()
            with open(edited, 'w', encoding='utf-8') as f:
                f.write(text.replace('Definitive', 'Disputed', 1))
            r = self.verify([tcg.MiniClinVar(), tcg.MiniErepo(), tcg.MiniGeneRegion()], syn_reference(), {'clingen_gene_validity': edited})
        finally:
            shutil.rmtree(d, ignore_errors=True)
        self.assertEqual(r['clingen_gene_validity']['n_mismatched'], 1)


class TestCpic(ResolverCase):
    def test_allele_function_and_recommendation_records_are_recovered(self):
        S = os.path.join(tpg.SOURCES, 'cpic', 'api')
        ref = Reference(SYN_FASTA, assembly_report=os.path.join(tpg.SOURCES, 'assembly_report.txt'), name='SYN38', builds=SYN_BUILDS)
        r = self.verify([tpg.CpicAdapter(tpg.FakeGenes(), root=os.path.join(tpg.SOURCES, 'cpic'))], ref,
                        {'cpic_api_allele': os.path.join(S, 'allele.json'), 'cpic_api_recommendation': os.path.join(S, 'recommendation.json')})
        self.assertAllVerified(r, ['cpic_api_allele', 'cpic_api_recommendation'])


class TestClinPgx(ResolverCase):
    def test_tsv_lines_and_json_members_of_the_zip_exports_are_recovered(self):
        P = os.path.join(tcp.SOURCES, 'pharmgkb')
        ref = Reference(SYN_FASTA, assembly_report=os.path.join(tcp.SOURCES, 'assembly_report.txt'), name='SYN38', builds=SYN_BUILDS)
        paths = {'clinpgx_summaryAnnotations': os.path.join(P, 'summaryAnnotations.zip'), 'clinpgx_variantAnnotations': os.path.join(P, 'variantAnnotations.zip'),
                 'clinpgx_drugLabels': os.path.join(P, 'drugLabels.zip'), 'clinpgx_guidelineAnnotations': os.path.join(P, 'guidelineAnnotations.json.zip')}
        r = self.verify(tcp.make_adapters(), ref, paths)
        self.assertAllVerified(r, list(paths))


class TestOpenTargets(ResolverCase):
    def test_credible_set_members_and_l2g_predictions_are_recovered(self):
        paths = {'opentargets_26_09_credible_set': os.path.join(tot.OT, 'credible_set'), 'opentargets_26_09_l2g_prediction': os.path.join(tot.OT, 'l2g_prediction')}
        r = self.verify(tot.adapters(), syn_reference(), paths)
        self.assertAllVerified(r, list(paths))

    def test_gene_evidence_rows_and_the_impc_summaries_are_recovered(self):
        paths = {'opentargets_26_09_evidence_crispr': os.path.join(tog.OT, 'evidence_crispr'), 'opentargets_26_09_evidence_crispr_screen': os.path.join(tog.OT, 'evidence_crispr_screen'),
                 'opentargets_26_09_evidence_gene_burden': os.path.join(tog.OT, 'evidence_gene_burden'), 'opentargets_26_09_target_essentiality': os.path.join(tog.OT, 'target_essentiality'),
                 'opentargets_26_09_evidence_impc': os.path.join(tog.OT, 'evidence_impc')}
        r = self.verify(tog.adapters(), syn_reference(), paths)
        self.assertAllVerified(r, list(paths))                                                # including the IMPC summaries, recomputed from every row their key selects


class TestMaveDbNative(ResolverCase):
    def test_annotation_records_of_the_native_genomic_mapping_are_recovered(self):
        mapped = functools.partial(resolvers.mavedb_native, annotations_dir=os.path.join(tmv.MV, 'csv'))
        mapped.__name__ = 'mavedb_native'
        r = self.verify(tmv.adapters(), tmv.reference(), {'mavedb_scores_genomic_mapped': os.path.join(tmv.MV, 'scores_genomic_mapped')}, override={'mavedb_scores_genomic_mapped': mapped})
        self.assertAllVerified(r, ['mavedb_scores_genomic_mapped'])


class TestStrchive(ResolverCase):
    def test_general_bed_lines_and_stranger_json_objects_are_recovered(self):
        root = os.path.join(trg.ROOT, 'strchive')
        ref = syn_reference()
        r = self.verify([StrchiveAdapter(GeneMap(ref, ot_dir='/nonexistent', mane_path='/nonexistent'), root=root)], ref,
                        {'strchive_hg38_general': os.path.join(root, 'STRchive-disease-loci.hg38.general.bed'), 'strchive_hg38_stranger_json': os.path.join(root, 'STRchive-disease-loci.hg38.stranger.json')})
        self.assertAllVerified(r, ['strchive_hg38_general', 'strchive_hg38_stranger_json'])


if __name__ == '__main__':
    unittest.main()
