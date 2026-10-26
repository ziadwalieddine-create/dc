import os
import tempfile
import unittest

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


if __name__ == '__main__':
    unittest.main()
