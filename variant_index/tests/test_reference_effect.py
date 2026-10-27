"""Evidence whose effect allele is the reference allele (about half of the GWAS Catalog's associations) is reached by a sample through a verified reference genotype as well as through an
alternate call. Without that the effect of such an allele would be reported for the people who carry the alternate allele and never for the people who carry two copies of the effect allele."""
import os
import shutil
import sys
import tempfile
import unittest

import pyarrow as pa
import pyarrow.compute as pc
import pyarrow.parquet as pq
import pysam

from helpers import syn_reference, vcf_text, FIXTURES, SYN_CONTIGS
from vindex import crosscheck, schema, wgs
from vindex.build import runner
from vindex.build.common import evidence_row
from vindex.index_reader import Index

sys.path.insert(0, FIXTURES)
import make_sources  # noqa: E402

b, alt = make_sources.b, make_sources.alt


class Associations:
    """An adapter that writes the associations of the test: (contig, position, [effect allele is the 'ref' allele or the 'alt' allele of the key, ...]) and, for an indel key, one row whose
    effect allele is the whole reference allele."""
    source_id = 'gwas_catalog'

    def __init__(self, sites, indel=None):
        self.sites, self.indel = sites, indel

    def run(self, ctx, limit=None):
        ctx.adapter = self
        rel = ctx.register_release('gwas_catalog', 'gwas_catalog_associations', 'test')
        for chrom, pos, effects in self.sites:
            ref, alt_ = b(chrom, pos), alt(chrom, pos)
            key = f'{chrom}:{pos}:{ref}:{alt_}'
            for i, which in enumerate(effects):
                effect, other = (ref, alt_) if which == 'ref' else (alt_, ref)
                ctx.emit_evidence(evidence_row(rel, 'gwas_catalog', f'GCST{chrom}{pos}:{i}', f'L{pos}.{i}', f'raw{chrom}{pos}.{i}', 'variant', key, 'effect_magnitude_direction', 'statistical_association',
                                               'gwas_effect_size', 'odds_ratio', 'effect', discriminator=f'L{pos}.{i}', chrom=chrom, pos=pos, effect_allele=effect, other_allele=other, odds_ratio=1.4,
                                               condition_text='test trait'))
        if self.indel:
            chrom, pos = self.indel
            ref2, alt2 = b(chrom, pos) + b(chrom, pos + 1), b(chrom, pos)
            ctx.emit_evidence(evidence_row(rel, 'gwas_catalog', f'GCST{chrom}{pos}:indel', f'L{pos}', f'raw{chrom}{pos}', 'variant', f'{chrom}:{pos}:{ref2}:{alt2}', 'effect_magnitude_direction',
                                           'statistical_association', 'gwas_effect_size', 'odds_ratio', 'effect', discriminator=f'L{pos}', chrom=chrom, pos=pos, effect_allele=ref2, other_allele=alt2,
                                           odds_ratio=1.4, condition_text='test trait'))
        return {}


class StubGenes:
    source_id = 'opentargets'

    def run(self, ctx, limit=None):
        ctx.adapter = self
        ctx.register_release('opentargets', 'opentargets_26_09_target', '26.09')
        return {}


SITES = [('1', 200, ['ref']),            # homozygous reference, inside the callable region: verified reference
         ('1', 210, ['ref']),            # outside the callable region and without a record: absent from the callset
         ('1', 220, ['ref']),            # a no-call record
         ('1', 230, ['alt']),            # about the alternate allele: a verified reference says nothing about it
         ('1', 240, ['ref']),            # homozygous alternate: linked through the call, with no copy of the effect allele
         ('1', 250, ['ref']),            # a filtered call
         ('1', 260, ['ref', 'alt']),     # one row about each allele: only the one about the reference allele is reached at a reference genotype
         ('X', 100, ['ref'])]            # a verified reference on X: the ploidy is not recorded


def key(chrom, pos):
    return f'{chrom}:{pos}:{b(chrom, pos)}:{alt(chrom, pos)}'


class TestVerifiedReferenceGenotypes(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.work = tempfile.mkdtemp()
        cls.idx = os.path.join(cls.work, 'idx')
        cls.ref = syn_reference()
        runner.build(cls.idx, [StubGenes(), Associations(SITES, indel=('1', 300))], cls.ref)
        gt = lambda g: f'{g}:60:30:10,10'
        recs = [('1', 220, b('1', 220), alt('1', 220), 'PASS', 50, [gt('./.')]),
                ('1', 240, b('1', 240), alt('1', 240), 'PASS', 50, [gt('1/1')]),
                ('1', 250, b('1', 250), alt('1', 250), 'LowQual', 5, [gt('0/1')])]
        cls.vcf = os.path.join(cls.work, 'sample.vcf')
        with open(cls.vcf, 'w') as f:
            f.write(vcf_text(recs, contigs=SYN_CONTIGS))
        pysam.tabix_index(cls.vcf, preset='vcf', force=True)
        cls.bed = os.path.join(cls.work, 'callable.bed')
        with open(cls.bed, 'w') as f:
            f.write('1\t194\t205\n1\t214\t265\nX\t89\t110\n')
        cls.out = os.path.join(cls.work, 'res')
        cls.manifest = wgs.run_sample(cls.vcf + '.gz', cls.idx, cls.out, cls.ref, callable_bed=cls.bed, repo_root=os.path.join(cls.work, 'not_the_repo'))
        cls.links = pq.read_table(os.path.join(cls.out, 'evidence_links.parquet')).to_pylist()

    @classmethod
    def tearDownClass(cls):
        shutil.rmtree(cls.work, ignore_errors=True)

    def reference_links(self, links=None):
        return [r for r in (links or self.links) if r['call_kind'] == 'verified_reference']

    def test_a_verified_reference_is_linked_to_the_rows_about_the_reference_allele_and_to_nothing_else(self):
        got = sorted({r['call_id'] for r in self.reference_links()})
        self.assertEqual(got, sorted([key('1', 200), key('1', 260), key('X', 100)]))
        for r in self.reference_links():
            self.assertEqual(r['effect_allele'], r['call_id'].split(':')[2], 'every link is about the reference allele')
            self.assertEqual((r['level'], r['source_id'], r['relation']), ('variant', 'gwas_catalog', 'verified_reference_genotype:callable_region_assumes_all_variants_reported'))

    def test_the_row_about_the_alternate_allele_of_the_same_site_is_not_reached(self):
        both = [r for r in self.links if r['call_id'] == key('1', 260)]
        self.assertEqual(len(both), 1)
        self.assertEqual(both[0]['effect_allele'], b('1', 260))
        evidence = Index.open(self.idx).evidence_for_variants([key('1', 260)])[key('1', 260)]
        self.assertEqual(sorted(e['effect_allele'] for e in evidence), sorted([b('1', 260), alt('1', 260)]), 'the index holds both rows')

    def test_two_copies_on_an_autosome_and_none_stated_where_the_ploidy_is_not_recorded(self):
        got = {r['call_id'].split(':')[0]: (r['effect_allele_dosage'], r['effect_allele_basis']) for r in self.reference_links() if r['call_id'].split(':')[1] in ('200', '100')}
        self.assertEqual(got['1'], (2, 'reference_genotype_diploid_autosome'))
        self.assertEqual(got['X'], (None, 'ploidy_of_a_verified_reference_not_recorded_outside_the_autosomes'))

    def test_a_site_the_sample_cannot_answer_for_has_no_link_and_is_counted_by_what_the_sample_says(self):
        linked = {r['call_id'] for r in self.reference_links()}
        for pos in (210, 220, 230, 240, 250):
            self.assertNotIn(key('1', pos), linked, pos)
        self.assertEqual(self.manifest['counts']['reference_allele_effect_sites'],
                         {'VERIFIED_REFERENCE': 3, 'ABSENT_FROM_CALLSET': 1, 'NO_CALL': 1, 'EXACT_MATCH': 1, 'FILTERED_CALL': 1})

    def test_a_site_that_is_called_is_linked_once_through_its_call(self):
        at_240 = [r for r in self.links if r['call_id'] == key('1', 240)]
        self.assertEqual([(r['call_kind'], r['effect_allele_dosage'], r['effect_allele_basis']) for r in at_240], [('small', 0, 'reference_allele_copies')])

    def test_an_indel_key_is_never_linked_this_way(self):
        self.assertFalse([r for r in self.links if r['call_id'].startswith('1:300:')])
        self.assertNotIn('1:300', ' '.join(self.manifest['counts']['reference_allele_effect_sites']))
        self.assertEqual(sum(self.manifest['counts']['reference_allele_effect_sites'].values()), 7)

    def test_without_a_callable_region_or_reference_blocks_nothing_is_verified(self):
        out = os.path.join(self.work, 'res_variant_only')
        m = wgs.run_sample(self.vcf + '.gz', self.idx, out, self.ref, repo_root=os.path.join(self.work, 'not_the_repo'))
        links = pq.read_table(os.path.join(out, 'evidence_links.parquet')).to_pylist()
        self.assertEqual(self.reference_links(links), [])
        self.assertEqual(m['counts']['reference_allele_effect_sites'], {'ABSENT_FROM_CALLSET': 4, 'NO_CALL': 1, 'EXACT_MATCH': 1, 'FILTERED_CALL': 1})

    def test_a_reference_block_verifies_by_its_own_basis_and_a_low_confidence_block_does_not(self):
        gt = lambda g, gq: f'{g}:{gq}:30:30,0'
        recs = [('1', 195, b('1', 195), '<NON_REF>', 'PASS', 0, [gt('0/0', 50)], 'END=205', '.', 'GT:GQ:DP:AD'),           # confident, covers 200
                ('1', 255, b('1', 255), '<NON_REF>', 'PASS', 0, [gt('0/0', 5)], 'END=265', '.', 'GT:GQ:DP:AD')]           # GQ far below the policy, covers 260
        vcf = os.path.join(self.work, 'gvcf.vcf')
        with open(vcf, 'w') as f:
            f.write(vcf_text(recs, contigs=SYN_CONTIGS))
        pysam.tabix_index(vcf, preset='vcf', force=True)
        out = os.path.join(self.work, 'res_gvcf')
        m = wgs.run_sample(vcf + '.gz', self.idx, out, self.ref, repo_root=os.path.join(self.work, 'not_the_repo'))
        links = self.reference_links(pq.read_table(os.path.join(out, 'evidence_links.parquet')).to_pylist())
        self.assertEqual([(r['call_id'], r['relation']) for r in links], [(key('1', 200), 'verified_reference_genotype:gvcf_reference_block')])
        self.assertEqual(m['counts']['reference_allele_effect_sites']['VERIFIED_REFERENCE'], 1)
        self.assertEqual(m['counts']['reference_allele_effect_sites']['UNRESOLVED'], 1)

    def test_the_run_is_reproducible_and_the_links_are_in_a_fixed_order(self):
        out = os.path.join(self.work, 'res_again')
        wgs.run_sample(self.vcf + '.gz', self.idx, out, self.ref, callable_bed=self.bed, repo_root=os.path.join(self.work, 'not_the_repo'))
        again = pq.read_table(os.path.join(out, 'evidence_links.parquet')).to_pylist()
        self.assertEqual(again, self.links)


@unittest.skipUnless(crosscheck.find_bcftools(), 'bcftools is not available (set VINDEX_TOOLS_BIN)')
class TestTheVerifiedSitesAreRecomputedWithBcftools(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.t = TestVerifiedReferenceGenotypes
        cls.t.setUpClass()

    @classmethod
    def tearDownClass(cls):
        cls.t.tearDownClass()

    def check(self, bed=None, results=None):
        return crosscheck.reference_sites(self.t.vcf + '.gz', results or self.t.out, self.t.idx, bed or self.t.bed, self.t.ref)

    def test_the_sites_the_run_verified_are_the_sites_bcftools_leaves_unblocked_inside_the_region(self):
        rep = self.check()
        self.assertTrue(rep['ok'], rep)
        # 8 indexed sites carry a row about the reference allele on a single-nucleotide key (200, 210, 220, 240, 250, 260, X:100 and 230's row is about the alternate allele: not one)
        self.assertEqual((rep['candidates'], rep['in_callable_region'], rep['blocked_by_a_record'], rep['expected'], rep['matcher']), (7, 6, 3, 3, 3))
        self.assertEqual(rep['per_contig']['1'], dict(candidates=6, callable=5, blocked=3, expected=2, matcher=2))

    def test_a_disagreement_is_reported_in_both_directions(self):
        narrower = os.path.join(self.t.work, 'narrower.bed')
        with open(narrower, 'w') as f:
            f.write('1\t214\t265\nX\t89\t110\n')                      # site 200 is no longer inside the region the check reads
        rep = self.check(bed=narrower)
        self.assertFalse(rep['ok'])
        self.assertEqual((rep['only_matcher'], rep['only_expected']), (1, 0))
        self.assertEqual(rep['examples']['only_matcher'], [key('1', 200)])
        wider = os.path.join(self.t.work, 'wider.bed')
        with open(wider, 'w') as f:
            f.write('1\t194\t300\nX\t89\t110\n')                      # the check reads 210 as inside, the run did not
        rep = self.check(bed=wider)
        self.assertEqual((rep['only_matcher'], rep['only_expected']), (0, 1))
        self.assertEqual(rep['examples']['only_expected'], [key('1', 210)])

    def test_a_link_the_matcher_should_not_have_made_is_found(self):
        forged = os.path.join(self.t.work, 'forged')
        shutil.copytree(self.t.out, forged)
        t = pq.read_table(os.path.join(forged, 'evidence_links.parquet'))
        row = [r for r in t.to_pylist() if r['call_kind'] == 'small'][0]                     # the link of a called site (240), presented as a verified reference
        row.update(call_kind='verified_reference', relation='verified_reference_genotype:forged')
        pq.write_table(pa.Table.from_pylist(t.to_pylist() + [row], schema=t.schema), os.path.join(forged, 'evidence_links.parquet'))
        rep = self.check(results=forged)
        self.assertFalse(rep['ok'])
        self.assertEqual(rep['examples']['only_matcher'], [row['call_id']])


class TestTheRowsThatAreAboutTheReferenceAllele(unittest.TestCase):
    """Index.reference_allele_effect_rows on a table in memory: the selection is a property of the rows, whatever source wrote them."""

    @staticmethod
    def row(subject, effect, other, status='mapped', level='variant', evidence_id=None):
        r = {n: None for n in schema.EVIDENCE.names}
        r.update(subject_id=subject, evidence_id=evidence_id or f'e-{subject}-{effect}', level=level, axis='effect_magnitude_direction', source_id='gwas_catalog', claim_type='c', effect_assertion='effect',
                lineage_group='g', effect_allele=effect, other_allele=other, mapping_status=status, chrom=subject.split(':')[0])
        return r

    def test_only_mapped_variant_rows_of_a_single_nucleotide_key_whose_effect_allele_is_the_reference(self):
        rows = [self.row('1:100:A:G', 'A', 'G'),                           # reference allele: selected
                self.row('1:100:A:G', 'G', 'A', evidence_id='e-alt'),      # alternate allele
                self.row('1:101:C:T', None, None, evidence_id='e-none'),   # no effect allele
                self.row('1:102:CA:C', 'CA', 'C'),                         # an indel key
                self.row('1:103:A:G', 'A', 'G', status='ambiguous', evidence_id='e-amb'),
                self.row('1:104:A:G', 'A', 'G', level='gene', evidence_id='e-gene'),
                self.row('2:100:A:G', 'A', 'G'),                           # another partition
                self.row('chr1_KI270706v1_random:50:A:G', 'A', 'G', evidence_id='e-alt-contig')]
        t = pa.Table.from_pylist(rows, schema=schema.EVIDENCE)
        idx = Index.from_tables(evidence=t)
        got = idx.reference_allele_effect_rows('1')
        self.assertEqual({k: [r['evidence_id'] for r in v] for k, v in got.items()}, {'1:100:A:G': ['e-1:100:A:G-A']})
        self.assertEqual(sorted(idx.reference_allele_effect_rows('2')), ['2:100:A:G'])
        self.assertEqual(sorted(idx.reference_allele_effect_rows(schema.OTHER_PARTITION)), ['chr1_KI270706v1_random:50:A:G'])
        self.assertEqual(idx.reference_allele_effect_rows('3'), {})

    def test_the_copies_by_chromosome(self):
        self.assertEqual(wgs.reference_effect_copies('22')[0], 2)
        for c in ('X', 'Y', 'MT', 'chr1_KI270706v1_random'):
            self.assertIsNone(wgs.reference_effect_copies(c)[0], c)


if __name__ == '__main__':
    unittest.main()
