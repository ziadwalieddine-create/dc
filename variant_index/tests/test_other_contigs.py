"""Contigs outside the primary assembly (alternate loci, unlocalised scaffolds, decoys) are indexed, matched and reported, never dropped:
  * a mapped ClinVar record on an alternate locus is filed under the shared 'NA' partition, becomes an entity_variant and is found by the reader;
  * a whole-file run reads the sample's calls on such contigs in one extra pass; a contig the reference does not know is listed as rejected, not skipped;
  * the validator notices a mapped row filed under the wrong partition, in either direction."""
import glob
import os
import random
import shutil
import tempfile
import unittest

import pysam
import pyarrow.parquet as pq

from helpers import SYN_FASTA, SYN_BUILDS, vcf_text
from vindex import reports, schema, views, wgs
from vindex.build import runner
from vindex.build.clinvar_vcf import ClinVarVcfAdapter
from vindex.build.genes import GeneMap
from vindex.index_reader import Index
from vindex.refgenome import Reference

ALT = 'ALT1_CTG1'
UNKNOWN = 'UNK_CTG'
MULTI = 'criteria_provided,_multiple_submitters,_no_conflicts'


def make_reference(work):
    """The synthetic primary contigs plus one alternate locus."""
    rng = random.Random(7)
    path = os.path.join(work, 'syn_with_alt.fa')
    with open(SYN_FASTA) as f, open(path, 'w') as g:
        g.write(f.read())
        seq = ''.join(rng.choice('ACGT') for _ in range(600))
        g.write(f'>{ALT}\n' + '\n'.join(seq[i:i + 60] for i in range(0, 600, 60)) + '\n')
    pysam.faidx(path)
    return Reference(path, name='SYN38', builds=SYN_BUILDS)


def clinvar_vcf(path, ref):
    def rec(c, pos, vid, ref_allele, alt, sig='Pathogenic', rev=MULTI):
        return '\t'.join([c, str(pos), str(vid), ref_allele, alt, '.', '.', f'ALLELEID={vid + 1000};CLNDISDB=MedGen:C9;CLNDN=Alt_locus_disorder;CLNREVSTAT={rev};CLNSIG={sig};CLNSIGSCV=SCV9;CLNVC=single_nucleotide_variant;ORIGIN=1'])
    other = lambda c, pos: 'T' if ref.fetch(c, pos, pos) != 'T' else 'G'
    base = lambda c, pos: ref.fetch(c, pos, pos)
    wrong = 'A' if base(ALT, 300) != 'A' else 'C'
    lines = ['##fileformat=VCFv4.1', '##fileDate=2026-01-01', '##source=ClinVar', '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO',
             rec('1', 402, 1, base('1', 402), other('1', 402)),
             rec(ALT, 100, 2, base(ALT, 100), other(ALT, 100)),                      # carried by the sample
             rec(ALT, 200, 3, base(ALT, 200), other(ALT, 200)),                      # not in the sample's file at all
             rec(ALT, 300, 4, wrong, other(ALT, 300))]                               # REF disagrees with the reference: rejected, never an entity
    with open(path, 'w') as f:
        f.write('\n'.join(lines) + '\n')
    return {'carried': f'{ALT}:100:{base(ALT, 100)}:{other(ALT, 100)}', 'absent': f'{ALT}:200:{base(ALT, 200)}:{other(ALT, 200)}', 'primary': f'1:402:{base("1", 402)}:{other("1", 402)}',
            'carried_pos': 100}


def sample_vcf(work, ref, keys):
    gt = lambda g: [f'{g}:60:30:15,15']
    contigs = [('1', 3000), (ALT, 600), (UNKNOWN, 500)]
    c, p, r, a = keys['carried'].split(':')
    recs = [('1', 402, ref.fetch('1', 402, 402), 'T' if ref.fetch('1', 402, 402) != 'T' else 'G', 'PASS', 50, gt('0/1')),
            (ALT, int(p), r, a, 'PASS', 50, gt('0/1')),
            (UNKNOWN, 10, 'A', 'C', 'PASS', 50, gt('0/1'))]
    path = os.path.join(work, 'sample.vcf')
    with open(path, 'w') as f:
        f.write(vcf_text(recs, contigs=contigs))
    pysam.tabix_index(path, preset='vcf', force=True)
    return path + '.gz'


class TestOtherContigs(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.work = tempfile.mkdtemp()
        cls.ref = make_reference(cls.work)
        cv = os.path.join(cls.work, 'cv.vcf')
        cls.keys = clinvar_vcf(cv, cls.ref)
        genes = GeneMap(cls.ref, ot_dir='/nonexistent', mane_path='/nonexistent')

        class Vcf(ClinVarVcfAdapter):
            def run(self, ctx, limit=None):
                return super().run(ctx, path=cv, limit=limit)
        cls.idx = os.path.join(cls.work, 'idx')
        cls.man = runner.build(cls.idx, [Vcf(genes)], cls.ref)

    @classmethod
    def tearDownClass(cls):
        shutil.rmtree(cls.work, ignore_errors=True)

    def entities(self):
        return {os.path.basename(f): pq.read_table(f)['variant_key'].to_pylist() for f in sorted(glob.glob(os.path.join(self.idx, 'entity_variant', '*.parquet')))}

    def test_alternate_locus_records_are_filed_in_the_shared_partition_and_become_entities(self):
        self.assertEqual(sorted(os.listdir(os.path.join(self.idx, 'evidence', 'variant'))), ['1', schema.OTHER_PARTITION])
        self.assertEqual(self.entities(), {'1.parquet': [self.keys['primary']], 'NA.parquet': [self.keys['carried'], self.keys['absent']]})
        rows = pq.read_table(glob.glob(os.path.join(self.idx, 'evidence', 'variant', 'NA', '*.parquet'))[0]).to_pylist()
        self.assertEqual({r['chrom'] for r in rows if r['mapping_status'] == 'mapped'}, {ALT})            # the row keeps its own contig name inside the shared partition

    def test_a_ref_mismatch_on_an_alternate_locus_is_not_a_mapped_variant(self):
        rej = [r for f in glob.glob(os.path.join(self.idx, 'evidence', 'variant', '*', '*.parquet')) for r in pq.read_table(f).to_pylist() if r['mapping_status'] != 'mapped']
        self.assertEqual([(r['mapping_status'], r['mapping_reason']) for r in rej], [('rejected', 'ref_mismatch')])
        self.assertTrue(all(k.split(':')[1] != '300' for ks in self.entities().values() for k in ks))

    def test_the_reader_finds_evidence_for_an_alternate_locus_key(self):
        index = Index.open(self.idx)
        ev = index.evidence_for_variants([self.keys['carried'], self.keys['primary'], f'{ALT}:5:A:C'])
        self.assertEqual(sorted(ev), sorted([self.keys['carried'], self.keys['primary']]))
        self.assertEqual(set(views.high_value_keys(index, schema.OTHER_PARTITION)), {self.keys['carried'], self.keys['absent']})
        self.assertEqual(set(views.high_value_keys(index, ALT)), {self.keys['carried'], self.keys['absent']})
        self.assertEqual(set(views.high_value_keys(index, '1')), {self.keys['primary']})

    def test_a_clean_index_validates(self):
        v = reports.validate(self.idx, self.ref)
        self.assertTrue(v['ok'], v['checks'])
        self.assertEqual(v['checks']['mapped_variant_row_filed_under_another_chromosome']['violations'], 0)

    def test_misfiled_rows_are_caught_in_both_directions(self):
        for moved_from, moved_to in (('NA', '1'), ('1', 'NA')):
            d = os.path.join(self.work, f'misfiled_{moved_from}')
            shutil.copytree(self.idx, d)
            src = glob.glob(os.path.join(d, 'evidence', 'variant', moved_from, '*.parquet'))[0]
            os.makedirs(os.path.join(d, 'evidence', 'variant', moved_to), exist_ok=True)
            shutil.move(src, os.path.join(d, 'evidence', 'variant', moved_to, os.path.basename(src)))
            v = reports.validate(d, self.ref)
            self.assertGreater(v['checks']['mapped_variant_row_filed_under_another_chromosome']['violations'], 0, moved_from)
            self.assertFalse(v['ok'])

    # ---- the sample ------------------------------------------------------------------------------------------------------------------
    def run_sample(self, **kw):
        out = os.path.join(self.work, 'res_' + str(len(os.listdir(self.work))))
        vcf = sample_vcf(self.work, self.ref, self.keys)
        man = wgs.run_sample(vcf, self.idx, out, self.ref, repo_root=os.path.join(self.work, 'not_the_repo'), **kw)
        return man, lambda name: pq.read_table(os.path.join(out, f'{name}.parquet')).to_pylist()

    def test_a_whole_file_run_reads_matches_and_reports_calls_off_the_primary_assembly(self):
        man, t = self.run_sample()
        self.assertEqual(man['input']['contigs'][-1], wgs.OTHER_PASS)
        pc = man['counts']['per_chromosome'][wgs.OTHER_PASS]
        self.assertEqual((pc['records'], pc['small_calls']), (2, 1))                               # the alternate-locus call and the record on a contig the reference lacks
        calls = {r['key']: r for r in t('small_calls')}
        self.assertEqual(set(calls), {self.keys['primary'], self.keys['carried']})
        self.assertEqual(calls[self.keys['carried']]['status'], 'EXACT_MATCH')
        self.assertGreaterEqual(calls[self.keys['carried']]['n_evidence'], 1)
        self.assertEqual([(r['chrom'], r['reason']) for r in t('rejected_calls')], [(UNKNOWN, 'contig_unmapped')])
        un = {r['key']: r for r in t('unresolved_loci')}
        self.assertEqual(set(un), {self.keys['absent']})
        self.assertEqual((un[self.keys['absent']]['chrom'], un[self.keys['absent']]['status']), (ALT, 'ABSENT_FROM_CALLSET'))
        rules = man['counts']['high_value_rules']['clinvar_pathogenic_reviewed']
        self.assertEqual(rules, {'EXACT_MATCH': 2, 'ABSENT_FROM_CALLSET': 1})                      # the primary allele and the alternate-locus one both carried; the third absent

    def test_a_restricted_run_does_not_read_off_primary_contigs(self):
        man, t = self.run_sample(contigs=['1'])
        self.assertEqual(man['input']['contigs'], ['1'])
        self.assertEqual({r['key'] for r in t('small_calls')}, {self.keys['primary']})
        self.assertEqual(t('rejected_calls'), [])
        self.assertEqual(t('unresolved_loci'), [])


if __name__ == '__main__':
    unittest.main()
