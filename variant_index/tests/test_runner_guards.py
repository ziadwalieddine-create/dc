import glob
import os
import tempfile
import unittest

import pyarrow.parquet as pq

from helpers import syn_reference
from vindex.build import runner
from vindex.build.common import evidence_row


class Emitter:
    requires = ()

    def __init__(self, source_id, sid):
        self.source_id, self.sid = source_id, sid

    def run(self, ctx, limit=None):
        ctx.adapter = self
        rel = ctx.register_release(self.sid, 'clinvar_vcf_grch38', 'x')
        ctx.emit_evidence(evidence_row(rel, self.sid, 'a1', 'L1', 'raw', 'variant', '1:10:A:C', 'clinical_validity_actionability', 'curated_expert_classification', 'c', 'v', 'effect'))
        return {}


class TestGuards(unittest.TestCase):
    def test_unknown_source_id_is_reported(self):
        man = runner.build(tempfile.mkdtemp(), [Emitter('typo', 'clinvr')], syn_reference())
        self.assertTrue(any("'clinvr'" in p and 'not a catalogue source' in p for p in man['problems']), man['problems'])

    def test_catalogue_source_id_passes(self):
        man = runner.build(tempfile.mkdtemp(), [Emitter('clinvar', 'clinvar')], syn_reference())
        self.assertEqual([p for p in man['problems'] if 'catalogue source' in p], [])

    def test_the_source_ids_used_by_the_real_adapters_are_catalogue_sources(self):
        known = runner.catalogue_sources()
        for sid in ('clinvar', 'clingen', 'gwas_catalog', 'cpic', 'pharmvar', 'pharmgkb', 'mavedb', 'opentargets', 'gnomad', 'gnomad_constraint', 'dgv', 'hgsvc3'):
            self.assertIn(sid, known)


class Offerer:
    source_id, requires = 'clinvar', ()

    def run(self, ctx, limit=None):
        ctx.adapter = self
        ctx.offer('ns_a', '1', 'E1')
        ctx.offer('ns_b', '1', 'E2')
        ctx.offer('ns_unwanted', '1', 'E3')
        return {}


class Wants:
    requires = ('clinvar',)

    def __init__(self, source_id, names):
        self.source_id, self.names, self.seen = source_id, names, None

    def prepare(self, ctx):
        for n in self.names:
            ctx.want(n)

    def run(self, ctx, limit=None):
        ctx.adapter = self
        self.seen = {k: list(v) for k, v in ctx.found.items()}
        return {}


class TestReleaseOfOfferedIds(unittest.TestCase):
    def test_ids_are_kept_until_every_adapter_that_asked_for_them_has_run(self):
        a, ab, last = Wants('probe_a', ['ns_a']), Wants('probe_ab', ['ns_a', 'ns_b']), Wants('probe_c', [])
        runner.build(tempfile.mkdtemp(), [Offerer(), a, ab, last], syn_reference())
        both = {('ns_a', '1'): ['E1'], ('ns_b', '1'): ['E2']}
        self.assertEqual(a.seen, both)                                                          # ns_a is asked for by two adapters, so it is still there for the second
        self.assertEqual(ab.seen, both)
        self.assertEqual(last.seen, {})                                                         # nobody left who asked for them; unrequested ids were never kept


class Rows:
    """Emits mapped variant-level rows with a chosen chrom column, to exercise the filing guard and the streaming entity_variant build."""
    requires = ()

    def __init__(self, source_id, specs, haplotype_keys=()):
        self.source_id, self.specs, self.haplotype_keys = source_id, specs, haplotype_keys

    def run(self, ctx, limit=None):
        ctx.adapter = self
        rel = ctx.register_release(self.source_id, 'clinvar_vcf_grch38', 'x')
        for i, (key, chrom) in enumerate(self.specs):
            row = evidence_row(rel, self.source_id, f'a{i}', f'L{i}', f'raw{i}', 'variant', key, 'clinical_validity_actionability', 'curated_expert_classification', 'c', 'v', 'effect')
            if chrom != '__unset__':
                row['chrom'] = chrom
            ctx.emit_evidence(row)
        for k in self.haplotype_keys:
            ctx.emit('haplotype_definition', dict(haplotype_id='H1', variant_key=k, release_id=rel, definition_scope='defining'))
        return {}


class TestEntityVariantStreaming(unittest.TestCase):
    def build(self, specs, **kw):
        out = tempfile.mkdtemp()
        man = runner.build(out, [Rows('clinvar', specs, **kw)], syn_reference())
        return out, man

    def keys_by_file(self, out):
        return {os.path.basename(f): pq.read_table(f)['variant_key'].to_pylist() for f in sorted(glob.glob(os.path.join(out, 'entity_variant', '*.parquet')))}

    def test_every_mapped_key_is_filed_under_its_own_chromosome_and_files_are_position_sorted(self):
        out, man = self.build([('2:5:A:C', '2'), ('1:300:A:C', '1'), ('1:20:G:T', '1'), ('X:7:A:T', 'X'), ('10:3:C:G', '10'), ('2:5:A:C', '2')])
        self.assertEqual(self.keys_by_file(out), {'1.parquet': ['1:20:G:T', '1:300:A:C'], '10.parquet': ['10:3:C:G'], '2.parquet': ['2:5:A:C'], 'X.parquet': ['X:7:A:T']})
        self.assertEqual(man['problems'], [])

    def test_a_row_without_chrom_is_filed_from_its_key_and_a_contradicting_chrom_is_refused(self):
        out, _ = self.build([('3:11:A:C', '__unset__')])
        self.assertEqual(self.keys_by_file(out), {'3.parquet': ['3:11:A:C']})
        files = glob.glob(os.path.join(out, 'evidence', 'variant', '3', '*.parquet'))
        self.assertEqual(len(files), 1)
        self.assertEqual(pq.read_table(files[0]).to_pylist()[0]['pos'], 11)
        with self.assertRaises(ValueError):
            runner.build(tempfile.mkdtemp(), [Rows('clinvar', [('3:11:A:C', '4')])], syn_reference())

    def test_haplotype_definition_keys_without_variant_evidence_still_become_entities(self):
        out, _ = self.build([('1:20:G:T', '1')], haplotype_keys=('7:9:A:G',))
        self.assertEqual(self.keys_by_file(out), {'1.parquet': ['1:20:G:T'], '7.parquet': ['7:9:A:G']})

    def test_an_index_without_variants_keeps_an_empty_valid_table(self):
        out, _ = self.build([])
        self.assertEqual(self.keys_by_file(out), {'empty.parquet': []})


if __name__ == '__main__':
    unittest.main()
