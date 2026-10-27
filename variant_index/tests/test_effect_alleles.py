"""An association is about one allele of a site. When that allele is the reference allele the sample is linked through its alternate call, so the link has to say how many copies of the
effect allele the sample carries: the effect of the reference allele points the other way from the effect of the alternate allele."""
import os
import shutil
import sys
import tempfile
import unittest
from types import SimpleNamespace

import pyarrow.parquet as pq
import pysam

from helpers import syn_reference, vcf_text, FIXTURES, SYN_CONTIGS
from vindex import wgs
from vindex.build import runner
from vindex.build.common import evidence_row

sys.path.insert(0, FIXTURES)
import make_sources  # noqa: E402

b, alt = make_sources.b, make_sources.alt


def call(ref='C', alt_='T', dosage=1, ref_copies=1, missing=0, vtype='snv'):
    return SimpleNamespace(ref=ref, alt=alt_, dosage=dosage, ref_copies=ref_copies, missing=missing, vtype=vtype)


class TestEffectAlleleState(unittest.TestCase):
    def test_the_alternate_allele_is_counted_from_the_genotype(self):
        self.assertEqual(wgs.effect_allele_state(call(dosage=1), 'T'), (1, 'alternate_allele_dosage'))
        self.assertEqual(wgs.effect_allele_state(call(dosage=2, ref_copies=0), 'T'), (2, 'alternate_allele_dosage'))

    def test_the_reference_allele_is_counted_as_the_copies_without_the_alternate(self):
        self.assertEqual(wgs.effect_allele_state(call(dosage=1, ref_copies=1), 'C'), (1, 'reference_allele_copies'))
        self.assertEqual(wgs.effect_allele_state(call(dosage=2, ref_copies=0), 'C'), (0, 'reference_allele_copies'))
        self.assertEqual(wgs.effect_allele_state(call(dosage=1, ref_copies=0), 'C'), (0, 'reference_allele_copies'), 'a genotype with another alternate allele has no reference copy')

    def test_nothing_is_guessed_when_the_effect_allele_or_the_genotype_does_not_allow_it(self):
        self.assertEqual(wgs.effect_allele_state(call(), None), (None, 'effect_allele_not_stated'))
        self.assertEqual(wgs.effect_allele_state(call(), ''), (None, 'effect_allele_not_stated'))
        self.assertEqual(wgs.effect_allele_state(call(), 'G'), (None, 'effect_allele_is_neither_allele_of_the_key'))
        self.assertEqual(wgs.effect_allele_state(call(missing=1), 'T'), (None, 'genotype_incomplete'))
        self.assertEqual(wgs.effect_allele_state(call(missing=1), 'C'), (None, 'genotype_incomplete'))
        for vtype, ref, alt_ in (('del', 'CA', 'C'), ('ins', 'C', 'CA'), ('mnv', 'CA', 'TG')):
            self.assertEqual(wgs.effect_allele_state(call(ref=ref, alt_=alt_, vtype=vtype), 'C'), (None, 'effect_allele_defined_for_single_nucleotide_keys_only'))
        self.assertEqual(wgs.effect_allele_state(None, 'C'), (None, 'effect_allele_defined_for_single_nucleotide_keys_only'))
        self.assertEqual(wgs.effect_allele_state(call(), 'CT'), (None, 'effect_allele_defined_for_single_nucleotide_keys_only'))


class Association:
    """An adapter that writes the associations of the test: (position, effect allele is the reference allele or the alternate allele)."""
    source_id = 'gwas_catalog'

    def __init__(self, sites):
        self.sites = sites

    def run(self, ctx, limit=None):
        ctx.adapter = self
        rel = ctx.register_release('gwas_catalog', 'gwas_catalog_associations', 'test')
        for pos, which in self.sites:
            ref, alt_ = b('1', pos), alt('1', pos)
            key = f'1:{pos}:{ref}:{alt_}'
            effect, other = (ref, alt_) if which == 'ref' else (alt_, ref)
            ctx.emit_evidence(evidence_row(rel, 'gwas_catalog', f'GCST{pos}:{pos}:trait', f'L{pos}', f'raw{pos}', 'variant', key, 'effect_magnitude_direction', 'statistical_association', 'gwas_effect_size',
                                           'odds_ratio', 'effect', discriminator=f'L{pos}', chrom='1', pos=pos, effect_allele=effect, other_allele=other, odds_ratio=1.4, condition_text='test trait'))
        return {}


class StubGenes:
    source_id = 'opentargets'

    def run(self, ctx, limit=None):
        ctx.adapter = self
        ctx.register_release('opentargets', 'opentargets_26_09_target', '26.09')
        return {}


class TestLinksCarryTheEffectAlleleAndTheSamplesCopies(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.work = tempfile.mkdtemp()
        cls.idx = os.path.join(cls.work, 'idx')
        cls.ref = syn_reference()
        cls.sites = [(150, 'alt'), (160, 'ref'), (170, 'ref'), (180, 'alt'), (190, 'ref')]
        runner.build(cls.idx, [StubGenes(), Association(cls.sites)], cls.ref)
        gt = lambda g, n=2: f'{g}:60:30:' + ','.join(['10'] * n)
        other = lambda pos, skip: alt('1', pos, skip=skip)
        recs = [('1', 150, b('1', 150), alt('1', 150), 'PASS', 50, [gt('1/1')]),                                   # effect = alternate, homozygous: 2 copies
                ('1', 160, b('1', 160), alt('1', 160), 'PASS', 50, [gt('0/1')]),                                   # effect = reference, heterozygous: 1 copy
                ('1', 170, b('1', 170), f"{alt('1', 170)},{other(170, (alt('1', 170),))}", 'PASS', 50, [gt('1/2', 3)]),   # effect = reference, but the genotype has two alternate alleles: 0 copies
                ('1', 180, b('1', 180), alt('1', 180), 'PASS', 50, [gt('./1')]),                                   # effect = alternate, half the genotype missing: not known
                ('1', 190, b('1', 190), alt('1', 190), 'PASS', 50, [gt('1/1')]),                                   # effect = reference, homozygous alternate: 0 copies
                ]
        recs.sort(key=lambda r: (r[0], r[1]))
        vcf = os.path.join(cls.work, 'sample.vcf')
        with open(vcf, 'w') as f:
            f.write(vcf_text(recs, contigs=SYN_CONTIGS))
        pysam.tabix_index(vcf, preset='vcf', force=True)
        cls.out = os.path.join(cls.work, 'res')
        wgs.run_sample(vcf + '.gz', cls.idx, cls.out, cls.ref, repo_root=os.path.join(cls.work, 'not_the_repo'))
        cls.links = {r['call_id'].split(':')[1]: r for r in pq.read_table(os.path.join(cls.out, 'evidence_links.parquet')).to_pylist()}

    @classmethod
    def tearDownClass(cls):
        shutil.rmtree(cls.work, ignore_errors=True)

    def test_each_link_states_the_effect_allele_and_the_copies_the_sample_carries(self):
        got = {pos: (r['effect_allele'], r['other_allele'], r['effect_allele_dosage'], r['effect_allele_basis']) for pos, r in self.links.items()}
        r150, a150 = b('1', 150), alt('1', 150)
        self.assertEqual(got['150'], (a150, r150, 2, 'alternate_allele_dosage'))
        self.assertEqual(got['160'], (b('1', 160), alt('1', 160), 1, 'reference_allele_copies'))
        self.assertEqual(got['170'], (b('1', 170), alt('1', 170), 0, 'reference_allele_copies'), 'with 1/2 no chromosome carries the reference allele (ploidy minus dosage would say 1)')
        self.assertEqual(got['180'][2:], (None, 'genotype_incomplete'))
        self.assertEqual(got['190'], (b('1', 190), alt('1', 190), 0, 'reference_allele_copies'))

    def test_a_link_to_an_event_has_no_effect_allele(self):
        rows = wgs.link_rows('ev', 'sv', [dict(evidence_id='e', level='event', axis='a', source_id='s', claim_type='c', effect_assertion='effect', lineage_group='g')], 'equivalent_event')
        self.assertEqual((rows[0]['effect_allele'], rows[0]['effect_allele_dosage'], rows[0]['effect_allele_basis']), (None, None, 'not_applicable'))


if __name__ == '__main__':
    unittest.main()
