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


if __name__ == '__main__':
    unittest.main()
