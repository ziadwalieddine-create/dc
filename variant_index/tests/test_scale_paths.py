"""The paths that keep a whole-genome match within memory must give the same answers as the plain ones: per-contig VCF reads, on-disk filtered index
reads against in-memory tables, and chunked interpretation against one pass."""
import os
import shutil
import sys
import tempfile
import unittest

import pysam

from helpers import syn_reference, vcf_text, FIXTURES
from vindex.build import runner
from vindex.build.clinvar_vcf import ClinVarVcfAdapter
from vindex.build.clinvar_vcv import ClinVarVcvAdapter
from vindex.build.genes import GeneMap
from vindex.build.regions import GiabRegionAdapter
from vindex.index_reader import Index
from vindex.match_small import MatchContext, iter_interpretations, interpret_calls
from vindex.vcfio import read_vcf

sys.path.insert(0, FIXTURES)
import make_sources  # noqa: E402

SOURCES = make_sources.write_all(tempfile.mkdtemp())
b, alt = make_sources.b, make_sources.alt


class MiniVcf(ClinVarVcfAdapter):
    def run(self, ctx, limit=None):
        return super().run(ctx, path=os.path.join(SOURCES, 'mini_clinvar.vcf'), limit=limit)


def write_vcf(tmp, indexed):
    recs = [('1', 402, b('1', 402), alt('1', 402), 'PASS', 50, ['0/1:60:30:15,15']), ('1', 430, b('1', 430), alt('1', 430), 'PASS', 50, ['1/1:60:30:0,30']),
            ('1', 600, b('1', 600), alt('1', 600), 'PASS', 50, ['0/1:60:30:15,15']), ('2', 100, b('2', 100), alt('2', 100), 'PASS', 50, ['0/1:60:30:15,15']),
            ('2', 200, b('2', 200), alt('2', 200), 'PASS', 50, ['1/1:60:30:0,30']), ('X', 50, b('X', 50), alt('X', 50), 'PASS', 50, ['1/1:60:30:0,30'])]
    p = os.path.join(tmp, 's.vcf')
    with open(p, 'w') as f:
        f.write(vcf_text(recs))
    if indexed:
        pysam.tabix_index(p, preset='vcf', force=True)
        return p + '.gz'
    return p


class TestPerContigReading(unittest.TestCase):
    def test_per_contig_reads_add_up_to_the_whole_file_indexed_or_not(self):
        ref = syn_reference()
        for indexed in (True, False):
            tmp = tempfile.mkdtemp()
            try:
                p = write_vcf(tmp, indexed)
                whole = read_vcf(p, ref)
                parts = {c: read_vcf(p, ref, contigs=[c]) for c in ('1', '2', 'X', 'MT')}
                self.assertEqual(sorted(c.key for c in whole.small), sorted(c.key for cs in parts.values() for c in cs.small), indexed)
                self.assertEqual({c: len(cs.small) for c, cs in parts.items()}, {'1': 3, '2': 2, 'X': 1, 'MT': 0})
                self.assertEqual(sum(cs.stats['records'] for cs in parts.values()), whole.stats['records'])
                two = read_vcf(p, ref, contigs=['2', 'X'])
                self.assertEqual(sorted(c.chrom for c in two.small), ['2', '2', 'X'])
            finally:
                shutil.rmtree(tmp, ignore_errors=True)


class TestIndexReadsAgree(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.ref = syn_reference()
        cls.out = tempfile.mkdtemp()
        g = GeneMap(cls.ref, ot_dir='/nonexistent', mane_path='/nonexistent')
        runner.build(cls.out, [MiniVcf(g), ClinVarVcvAdapter(g, root=os.path.join(SOURCES, 'vcv')), GiabRegionAdapter(root=os.path.join(SOURCES, 'regions'))], cls.ref)
        cls.disk = Index.open(cls.out)
        cls.mem = Index.from_tables(**{n: Index.open(cls.out).table(n) for n in ('evidence', 'entity_event', 'entity_region', 'entity_gene', 'entity_variant')})

    @classmethod
    def tearDownClass(cls):
        shutil.rmtree(cls.out, ignore_errors=True)

    def test_overlap_queries_agree(self):
        for chrom, hi in (('1', 3000), ('2', 1500), ('X', 1500), ('MT', 400)):
            for s in range(1, hi, 37):
                for width in (0, 25, 400):
                    q = (chrom, s, s + width)
                    key = lambda r: (r['kind'] if 'kind' in r else '', r.get('region_id') or r.get('event_id') or r.get('gene_id'))
                    self.assertEqual(sorted(self.disk.regions_overlapping(*q), key=key), sorted(self.mem.regions_overlapping(*q), key=key), q)
                    self.assertEqual(sorted(self.disk.events_overlapping(*q), key=key), sorted(self.mem.events_overlapping(*q), key=key), q)
        self.assertTrue(self.disk.events_overlapping('1', 1050, 1060))              # the VCV deletion: the comparison above was not between two empty answers
        self.assertTrue(self.disk.regions_overlapping('1', 150, 150))

    def test_region_kind_filter_and_cache_eviction(self):
        kinds = lambda i, c, s, e, k=None: sorted(r['kind'] for r in i.regions_overlapping(c, s, e, k))
        self.assertEqual(kinds(self.disk, '1', 150, 150, ['giab_difficult_segdups']), ['giab_difficult_segdups'])
        for c in ('1', '2', 'X', 'MT', '1'):                                     # more chromosomes than are kept, then back again
            self.disk.regions_overlapping(c, 1, 10)
        self.assertEqual(kinds(self.disk, '1', 150, 150), kinds(self.mem, '1', 150, 150))
        self.assertLessEqual(len([k for k in self.disk._c if k[0] == 'entity_region']), 3)

    def test_evidence_reads_agree_across_chromosomes(self):
        keys = {f'1:402:{b("1", 402)}:{alt("1", 402)}', f'1:410:{b("1", 410)}:{alt("1", 410)}', f'X:1:A:C', f'2:5:A:C'}
        d, m = self.disk.evidence_for_variants(keys), self.mem.evidence_for_variants(keys)
        self.assertEqual({k: [r['evidence_id'] for r in v] for k, v in d.items()}, {k: [r['evidence_id'] for r in v] for k, v in m.items()})
        self.assertTrue(d)
        self.assertEqual(self.disk.evidence_for_variants([]), {})
        self.assertEqual(self.disk.evidence_for_variants(['MT:1:A:C', 'GL000:5:A:C']), {})                    # partitions that do not exist hold nothing
        ids = ['clinvar:VCV000002001', 'clinvar:VCV000003001']
        self.assertEqual({k: len(v) for k, v in self.disk.evidence_for('event', ids).items()}, {k: len(v) for k, v in self.mem.evidence_for('event', ids).items()})


class TestChunkedInterpretation(unittest.TestCase):
    def test_chunk_size_does_not_change_the_answer(self):
        ref = syn_reference()
        tmp, out = tempfile.mkdtemp(), tempfile.mkdtemp()
        try:
            runner.build(out, [MiniVcf(GeneMap(ref, ot_dir='/nonexistent', mane_path='/nonexistent'))], ref)
            idx = Index.open(out)
            cs = read_vcf(write_vcf(tmp, True), ref)
            whole = [(o.call.key, o.status, o.axes and sorted(o.axes), o.flags) for o in interpret_calls(cs, idx, MatchContext(cs))]
            for chunk in (1, 2, 1000):
                got = [(o.call.key, o.status, o.axes and sorted(o.axes), o.flags) for o in iter_interpretations(cs, idx, MatchContext(cs), chunk=chunk)]
                self.assertEqual(got, whole, chunk)
            self.assertIn('EXACT_MATCH', {w[1] for w in whole})
            self.assertIn('NO_INDEX_RECORD', {w[1] for w in whole})
        finally:
            shutil.rmtree(tmp, ignore_errors=True)
            shutil.rmtree(out, ignore_errors=True)


if __name__ == '__main__':
    unittest.main()
