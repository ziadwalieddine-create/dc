import json
import os
import sys
import tempfile
import unittest

from helpers import FIXTURES, SYN_FASTA, SYN_BUILDS
from vindex.build import runner
from vindex.build.clinpgx import ClinPgxAdapter, choose_alt, letters_of, parse_p
from vindex.build.pgx import CpicAdapter, PharmVarAdapter
from vindex.index_reader import Index
from vindex.refgenome import Reference

sys.path.insert(0, FIXTURES)
import make_sources  # noqa: E402

SOURCES = make_sources.write_all(tempfile.mkdtemp())
KEYS = make_sources.write_all.clinpgx


class FakeGenes:
    syms = {'TESTG': 'ENSG_TESTG', 'TESTH': 'ENSG_TESTH', 'TESTK': 'ENSG_TESTK'}

    def resolve(self, hgnc=None, symbol=None, ncbi=None):
        return (self.syms[symbol], 'symbol') if symbol in self.syms else (None, 'unresolved')


def make_adapters():
    return [CpicAdapter(FakeGenes(), root=os.path.join(SOURCES, 'cpic')), PharmVarAdapter(os.path.join(SOURCES, 'pharmvar-test.zip')),
            ClinPgxAdapter(FakeGenes(), root=os.path.join(SOURCES, 'pharmgkb'))]


class TestChooseAlt(unittest.TestCase):
    """Which alternate allele of a multiallelic site an annotation is about, from the bases it names."""
    sites = ['1:710:A:C', '1:710:A:T']

    def test_single_candidate_records_orientation(self):
        self.assertEqual(choose_alt(['1:705:G:A'], {'G', 'A'}), ('mapped', '1:705:G:A', None, ('annotation_alleles_plus_strand',)))
        self.assertEqual(choose_alt(['1:705:G:A'], {'C', 'T'}), ('mapped', '1:705:G:A', None, ('annotation_alleles_minus_strand',)))
        self.assertEqual(choose_alt(['1:705:G:A'], {'G', 'C'})[3], ('annotation_alleles_not_among_site_alleles',))
        self.assertEqual(choose_alt(['1:705:G:A'], {'T'})[3], ('annotation_alleles_minus_strand',))
        self.assertEqual(choose_alt(['1:705:G:A'], set())[3], ())
        self.assertEqual(choose_alt([], {'A'})[:3], ('unmapped', None, 'no_alternate_allele_known'))

    def test_multiallelic_follows_the_documented_positive_strand(self):
        st, key, _r, flags = choose_alt(self.sites, {'C', 'A'})
        self.assertEqual((st, key), ('mapped', '1:710:A:C'))
        self.assertIn('strand_by_clinpgx_positive_strand_convention', flags)
        self.assertNotIn('annotation_alleles_also_fit_the_opposite_strand', flags)
        # T is an allele of the site on the positive strand, and its complement A is the reference: both readings fit, the convention decides
        st, key, _r, flags = choose_alt(self.sites, {'T'})
        self.assertEqual((st, key), ('mapped', '1:710:A:T'))
        self.assertIn('annotation_alleles_also_fit_the_opposite_strand', flags)

    def test_multiallelic_ambiguous_rejected_and_contradicting_the_convention(self):
        self.assertEqual(choose_alt(self.sites, {'C', 'T'})[:3], ('ambiguous', None, 'annotation_spans_several_alternate_alleles'))
        self.assertEqual(choose_alt(self.sites, set())[:3], ('ambiguous', None, 'several_alternate_alleles_annotation_names_no_bases'))
        self.assertEqual(choose_alt(self.sites, {'A'})[:3], ('ambiguous', None, 'annotation_names_only_the_reference_allele'))
        self.assertEqual(choose_alt(self.sites, {'G', 'C'})[:3], ('rejected', None, 'annotation_alleles_not_among_site_alleles'))    # neither strand reads as alleles of the site
        self.assertEqual(choose_alt(self.sites, {'G'})[:3], ('ambiguous', None, 'annotation_alleles_fit_only_the_opposite_strand'))
        self.assertEqual(choose_alt(['1:710:A:C', '1:710:A:T'], {'G', 'T'})[:3], ('ambiguous', None, 'annotation_alleles_fit_only_the_opposite_strand'))

    def test_mixed_snv_and_indel_site(self):
        st, key, _r, flags = choose_alt(['1:710:A:C', '1:710:AT:A'], {'C'})
        self.assertEqual((st, key), ('mapped', '1:710:A:C'))
        self.assertIn('site_also_has_non_snv_alleles', flags)
        self.assertEqual(choose_alt(['1:710:A:C', '1:710:AT:A'], set())[:3], ('ambiguous', None, 'several_alternate_alleles_annotation_names_no_bases'))
        self.assertEqual(choose_alt(['1:710:AT:A', '1:710:AT:AG'], {'C'})[2], 'several_alternate_alleles_no_single_base_reference')

    def test_letters_and_p_values(self):
        self.assertEqual(letters_of('CT + TT', '*1/*2', 'del/del', 'A'), {'C', 'T', 'A'})
        self.assertEqual(parse_p('= 0.012'), (0.012, '= 0.012'))
        self.assertEqual(parse_p('< 0.001'), (None, '< 0.001'))
        self.assertEqual(parse_p('= 1E-5')[0], 1e-5)
        self.assertEqual(parse_p('= 5'), (None, '= 5'))          # not a probability: kept as text only
        self.assertEqual(parse_p(''), (None, None))


class TestClinPgx(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.ref = Reference(SYN_FASTA, assembly_report=os.path.join(SOURCES, 'assembly_report.txt'), name='SYN38', builds=SYN_BUILDS)
        cls.out = tempfile.mkdtemp()
        cls.man = runner.build(cls.out, make_adapters(), cls.ref)
        cls.idx = Index.open(cls.out)
        cls.ev = [r for r in cls.idx.table('evidence').to_pylist() if r['source_id'] == 'pharmgkb']
        cls.stats = {a['adapter']: a['stats'] for a in cls.man['adapters']}['ClinPgxAdapter']

    def rows(self, acc, **kw):
        out = [r for r in self.ev if r['accession'] == acc]
        for k, v in kw.items():
            out = [r for r in out if r[k] == v]
        return out

    def one(self, acc, **kw):
        out = self.rows(acc, **kw)
        self.assertEqual(len(out), 1, (acc, kw, [(r['subject_id'], r['mapping_status']) for r in out]))
        return out[0]

    # ---- summary annotations ------------------------------------------------------------------------------------------------------
    def test_summary_rsid_biallelic(self):
        r = self.one('summary_annotation:1001')
        self.assertEqual((r['level'], r['subject_id'], r['mapping_status']), ('variant', KEYS['k705'], 'mapped'))
        self.assertEqual((r['axis'], r['evidence_kind'], r['claim_type'], r['claim_value']), ('association_replication', 'aggregated_integration', 'clinpgx_level_of_evidence', '1A'))
        self.assertEqual((r['effect_assertion'], r['score'], r['score_unit'], r['pos'], r['chrom']), ('effect', 12.5, 'clinpgx_annotation_score', 705, '1'))
        self.assertEqual(r['gene_ids'], ['ENSG_TESTG'])
        self.assertEqual(r['lineage_group'], 'clinpgx_summary:1001')
        extra = json.loads(r['extra_json'])
        self.assertEqual([a['Genotype/Allele'] for a in extra['genotype_annotations']], ['GG', 'AG', 'AA'])
        self.assertEqual(extra['evidence'][0]['PMID'], '111')
        self.assertIn('annotation_alleles_plus_strand', extra['flags'])

    def test_multiallelic_site_resolved_by_the_named_bases_or_kept_ambiguous(self):
        r = self.one('summary_annotation:1002')
        self.assertEqual((r['mapping_status'], r['subject_id']), ('mapped', KEYS['k710c']))             # genotypes CC and AC name only C as the alternate
        t = self.one('summary_annotation:1003')                                                           # genotype T fits both strands; the documented convention decides
        self.assertEqual((t['mapping_status'], t['subject_id']), ('mapped', KEYS['k710t']))
        self.assertIn('annotation_alleles_also_fit_the_opposite_strand', json.loads(t['extra_json'])['flags'])
        amb = self.one('summary_annotation:1017')                                                         # genotypes CC and TT: two alternate alleles
        self.assertEqual((amb['mapping_status'], amb['subject_id'], amb['mapping_reason'], amb['unmapped_key']), ('ambiguous', None, 'annotation_spans_several_alternate_alleles', 'rs710'))
        self.assertEqual(json.loads(amb['extra_json'])['candidate_keys'], [KEYS['k710c'], KEYS['k710t']])
        self.assertEqual((amb['chrom'], amb['pos']), ('1', 710))

    def test_unlocated_and_unresolvable_rsids_stay_unmapped(self):
        self.assertEqual(self.one('summary_annotation:1004')['mapping_reason'], 'no_grch38_location_in_clinpgx')
        r = self.one('summary_annotation:1005')
        self.assertEqual((r['mapping_status'], r['mapping_reason'], r['pos']), ('ambiguous', 'alternate_allele_not_given', 402))
        self.assertEqual(self.one('summary_annotation:1013')['mapping_reason'], 'location_on_a_contig_of_another_assembly:NC_999999.1')
        rep = self.one('summary_annotation:1015')
        self.assertEqual((rep['mapping_status'], rep['mapping_reason'], rep['pos']), ('unmapped', 'repeat_length_polymorphism', 200))      # position kept for the unresolved-loci report
        self.assertEqual(self.one('summary_annotation:1014')['subject_id'], KEYS['kdel'])                 # a deletion synonym resolves to its normalised key

    def test_haplotypes_resolve_through_the_definitions_and_reference_alleles_are_not_the_subject(self):
        r = self.one('summary_annotation:1006')
        self.assertEqual((r['level'], r['subject_id'], r['mapping_status']), ('haplotype', 'cpic:TESTG:*2', 'mapped'))
        self.assertEqual(json.loads(r['extra_json'])['also_defined_as'], ['pharmvar:TESTG:*2'])
        self.assertEqual(self.rows('summary_annotation:1006', subject_id='cpic:TESTG:*1'), [])
        un = self.one('summary_annotation:1009')
        self.assertEqual((un['level'], un['mapping_status'], un['mapping_reason'], un['unmapped_key']), ('haplotype', 'unmapped', 'allele_not_defined_in_cpic_or_pharmvar', 'TESTG*9'))
        self.assertEqual(self.one('summary_annotation:1018')['mapping_reason'], 'hla_allele_not_defined_in_index_needs_hla_typing')

    def test_named_alleles_with_commas_are_rejoined(self):
        r = self.one('summary_annotation:1008')
        self.assertEqual((r['subject_id'], r['mapping_status']), ('cpic:TESTK:Odd, name', 'mapped'))

    def test_gene_level_descriptor_stays_on_the_gene(self):
        r = self.one('summary_annotation:1007')
        self.assertEqual((r['level'], r['subject_id'], r['mapping_reason']), ('gene', 'ENSG_TESTG', 'gene_level_descriptor'))
        self.assertEqual(self.one('summary_annotation:1016')['mapping_reason'], 'named_polymorphism_without_definition_or_coordinates')
        self.assertEqual(self.one('summary_annotation:1012')['mapping_reason'], 'transcript_or_protein_hgvs_not_resolved')

    def test_several_variants_in_one_annotation_each_get_a_row_sharing_one_lineage_group(self):
        rows = self.rows('summary_annotation:1010')
        self.assertEqual(sorted((r['subject_id'] or r['unmapped_key']) for r in rows), sorted([KEYS['k705'], 'rs710']))
        self.assertEqual(len({r['lineage_group'] for r in rows}), 1)
        self.assertEqual(rows[0]['gene_ids'], ['ENSG_TESTG', 'ENSG_TESTH'])

    def test_level_four_is_uncertain(self):
        self.assertEqual(self.one('summary_annotation:1011')['effect_assertion'], 'uncertain')

    # ---- variant annotations ------------------------------------------------------------------------------------------------------
    def test_variant_annotation_association_row(self):
        r = self.one('variant_annotation:2001')
        self.assertEqual((r['level'], r['subject_id'], r['axis'], r['evidence_kind'], r['claim_type']), ('variant', KEYS['k705'], 'association_replication', 'statistical_association', 'pharmacogenomic_association'))
        self.assertEqual((r['effect_assertion'], r['lineage_group'], r['claim_value']), ('effect', 'pmid:111', 'sentence 2001'))
        self.assertIsNone(r['effect_allele'])                  # no effect allele or direction is derived from genotype text
        self.assertIsNone(r['direction'])
        extra = json.loads(r['extra_json'])
        self.assertEqual((extra['alleles'], extra['comparison'], extra['direction_of_effect']), ('AA', 'GG', 'increased'))
        self.assertEqual(len(extra['studies']), 2)
        self.assertEqual(r['ancestry'], 'East Asian; European')

    def test_ratio_rows_are_separate_effect_size_evidence(self):
        orr = self.one('variant_annotation:2001:params:3001')
        self.assertEqual((orr['axis'], orr['claim_type'], orr['odds_ratio'], orr['ci_low'], orr['ci_high'], orr['p_value']), ('effect_magnitude_direction', 'ratio_statistic:OR', 2.5, 1.1, 5.7, 0.012))
        self.assertEqual((orr['sample_size'], orr['n_cases'], orr['ancestry'], orr['subject_id']), (300, 100, 'European', KEYS['k705']))
        self.assertEqual(json.loads(orr['extra_json'])['reference_unit'], 'genotype_comparison_not_alt_allele')
        hr = self.one('variant_annotation:2001:params:3002')
        self.assertEqual((hr['odds_ratio'], hr['score'], hr['score_unit'], hr['p_value']), (None, 1.8, 'HR', None))      # '< 0.001' is not a p-value
        self.assertEqual(json.loads(hr['extra_json'])['p_value_text'], '< 0.001')
        self.assertEqual(self.rows('variant_annotation:2002:params:3003'), [])                                    # no usable ratio: no effect-size row
        self.assertEqual(self.one('variant_annotation:2101:params:3004')['odds_ratio'], 0.4)

    def test_negative_finding_is_kept_and_shares_the_paper_lineage(self):
        r = self.one('variant_annotation:2002')
        self.assertEqual((r['effect_assertion'], r['lineage_group']), ('no_effect', 'pmid:111'))

    def test_not_stated_significance_is_uncertain_and_the_allele_of_interest_is_the_subject(self):
        r = self.one('variant_annotation:2003')
        self.assertEqual((r['effect_assertion'], r['subject_id']), ('uncertain', 'cpic:TESTG:*2'))
        only_ref = self.one('variant_annotation:2005')
        self.assertEqual(only_ref['subject_id'], 'cpic:TESTG:*1')                                                  # nothing else named: the finding is about the reference allele

    def test_annotation_that_contradicts_the_strand_convention_is_not_used_to_pick_an_allele(self):
        r = self.one('variant_annotation:2004')                       # genotype GG: only the opposite strand fits the two alternate alleles C and T
        self.assertEqual((r['mapping_status'], r['subject_id'], r['mapping_reason']), ('ambiguous', None, 'annotation_alleles_fit_only_the_opposite_strand'))

    def test_functional_assay_goes_on_the_functional_axis(self):
        r = self.one('variant_annotation:2201')
        self.assertEqual((r['axis'], r['evidence_kind'], r['claim_type'], r['model_system'], r['cell_type'], r['subject_id']),
                         ('empirical_functional_support', 'direct_experimental', 'pharmacogenomic_functional_assay', 'in human liver microsomes', 'hepatocyte', 'cpic:TESTG:*2'))

    def test_studies_are_entities_with_the_paper_lineage(self):
        st = {r['study_id']: r for r in self.idx.table('entity_study').to_pylist() if r['source_id'] == 'pharmgkb'}
        self.assertEqual(sorted(st), ['clinpgx_va:2001', 'clinpgx_va:2002', 'clinpgx_va:2003', 'clinpgx_va:2004', 'clinpgx_va:2005', 'clinpgx_va:2101', 'clinpgx_va:2201'])
        self.assertEqual((st['clinpgx_va:2001']['pubmed'], st['clinpgx_va:2001']['lineage_group'], st['clinpgx_va:2001']['study_type']), ('111', 'pmid:111', 'case/control'))
        self.assertEqual(st['clinpgx_va:2101']['n_total'], 22)
        self.assertEqual(st['clinpgx_va:2201']['study_type'], 'functional_assay')

    # ---- guidelines, labels ------------------------------------------------------------------------------------------------------
    def test_guideline_rows_gene_and_named_alleles(self):
        gene = self.one('guideline_annotation:PA1000', level='gene')
        self.assertEqual((gene['subject_id'], gene['claim_value'], gene['effect_assertion'], gene['curator']), ('ENSG_TESTG', 'CPIC: recommendation provided', 'not_applicable', 'CPIC'))
        rows = self.rows('guideline_annotation:PA1000')
        self.assertEqual(sorted((r['level'], r['subject_id']) for r in rows), [('gene', 'ENSG_TESTG'), ('haplotype', 'cpic:TESTG:*2'), ('variant', KEYS['k705'])])
        no = self.rows('guideline_annotation:PA1001')
        self.assertEqual({r['claim_value'] for r in no}, {'DPWG: no recommendation'})
        self.assertEqual(sorted(r['mapping_status'] for r in no), ['mapped', 'unmapped'])           # TESTH is in the gene map, UNKNOWNX is not

    def test_drug_label_rows_are_gene_level(self):
        rows = self.rows('drug_label:PA2000')
        self.assertEqual(sorted(r['subject_id'] for r in rows), ['ENSG_TESTG', 'ENSG_TESTH'])
        self.assertEqual({(r['claim_value'], r['curator'], r['axis']) for r in rows}, {('Actionable PGx', 'FDA', 'clinical_validity_actionability')})

    # ---- cross-checks and release versions ------------------------------------------------------------------------------------------
    def test_release_versions_come_from_the_archives(self):
        rels = set(self.man['releases'])
        self.assertIn('pharmgkb/clinpgx_summaryAnnotations@2026-10-05', rels)
        self.assertIn('pharmgkb/clinpgx_variantAnnotations@2026-10-05', rels)

    def test_legacy_release_is_compared_not_integrated(self):
        s = self.stats
        self.assertEqual((s['legacy_ids_not_in_summary'], s['summary_ids_not_in_legacy'], s['legacy_vs_summary_changed:Level of Evidence']), (1, 1, 1))
        self.assertEqual(s['clinical_variants_rows'], 2)
        self.assertEqual(s['clinical_variants_without_matching_summary'], 1)
        self.assertFalse([r for r in self.ev if r['accession'] == 'summary_annotation:1099'])
        self.assertTrue(any('only the newer export is integrated' in p for p in self.man['problems']))

    def test_sites_link_dbsnp_names_to_keys_and_genes(self):
        x = [r for r in self.idx.table('xref').to_pylist() if r['namespace'] == 'dbsnp' and r['entity_id'] == KEYS['k705'] and r['release_id'].startswith('pharmgkb/')]
        self.assertEqual([(r['identifier'], r['relation']) for r in x], [('rs705', 'clinpgx_variant_name')])
        every = [(r['identifier'], r['release_id'], r['relation']) for r in self.idx.table('xref').to_pylist() if r['namespace'] == 'dbsnp' and r['entity_id'] == KEYS['k705']]
        self.assertEqual(len(every), len(set(every)))            # one row per source statement, never one per definition row
        vg = [r for r in self.idx.table('variant_gene').to_pylist() if r['variant_key'] == KEYS['k705']]
        self.assertEqual([(r['gene_id'], r['relation']) for r in vg], [('ENSG_TESTG', 'clinpgx_listed_gene')])

    def test_every_row_validates_and_ids_are_unique(self):
        ids = [r['evidence_id'] for r in self.ev]
        self.assertEqual(len(ids), len(set(ids)))
        self.assertTrue(all(r['raw_sha256'] and r['locator'] for r in self.ev))

    def test_rebuild_is_idempotent(self):
        out2 = tempfile.mkdtemp()
        runner.build(out2, make_adapters(), self.ref)
        with open(os.path.join(self.out, 'table_digests.json')) as f1, open(os.path.join(out2, 'table_digests.json')) as f2:
            self.assertEqual(json.load(f1), json.load(f2))

    def test_adapter_order_is_enforced(self):
        with self.assertRaises(ValueError):
            runner.build(tempfile.mkdtemp(), [ClinPgxAdapter(FakeGenes(), root=os.path.join(SOURCES, 'pharmgkb'))], self.ref)


if __name__ == '__main__':
    unittest.main()
