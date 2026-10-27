import glob
import json
import os
import tempfile
import unittest

import pyarrow.parquet as pq

from helpers import syn_reference
from vindex.build import runner
from vindex.build.common import REPO, curie, curies, evidence_row, iso_date, split_identifiers


class Emitter:
    requires = ()

    def __init__(self, source_id, sid):
        self.source_id, self.sid = source_id, sid

    def run(self, ctx, limit=None):
        ctx.adapter = self
        rel = ctx.register_release(self.sid, 'clinvar_vcf_grch38', 'x')
        ctx.emit_evidence(evidence_row(rel, self.sid, 'a1', 'L1', 'raw', 'variant', '1:10:A:C', 'clinical_validity_actionability', 'curated_expert_classification', 'c', 'v', 'effect'))
        return {}


class DerivedTargetAdapter:
    requires = ()
    source_id = 'mavedb'

    def run(self, ctx, limit=None):
        ctx.adapter = self
        ctx.register_release('mavedb', 'mavedb_scores_genomic_mapped', 'x')
        return {}


class TestOntologyIdentifiersHaveOneSpelling(unittest.TestCase):
    """Found by reading the identifiers of the final adapters side by side: Open Targets wrote MONDO_0005393, the GWAS Catalog's URIs gave EFO:0000384 and ClinVar's VCF MONDO:MONDO:0010178 and
    Human_Phenotype_Ontology:HP:0001250, so one disease was several entries in every grouping and two sources could not contradict each other."""

    def test_the_spellings_the_sources_use_become_prefix_colon_local(self):
        for written, spelled in (('MONDO_0005393', 'MONDO:0005393'), ('Orphanet_199', 'Orphanet:199'), ('EFO_0004339', 'EFO:0004339'), ('MONDO:MONDO:0010178', 'MONDO:0010178'),
                                 ('HP:HP:0002282', 'HP:0002282'), ('Human_Phenotype_Ontology:HP:0001250', 'HP:0001250'), ('GO_0000050', 'GO:0000050'), (' EFO_1 ', 'EFO:1'), ('OMIM: 612621', 'OMIM:612621'), (' OMIM:   612621 ', 'OMIM:612621'), ('Human_Phenotype_Ontology: HP:0001250', 'HP:0001250')):
            self.assertEqual(curie(written), spelled, written)

    def test_an_identifier_that_has_the_form_is_not_touched(self):
        for same in ('EFO:0000384', 'MedGen:C0403814', 'OMIM:277180', 'Orphanet:48', 'MeSH:D012345', 'RxNorm:1048', 'ClinGen:CA000001', 'HP:0000118', 'pharmgkb:PA123', 'ENSG00000141510',
                     'OMIM phenotypic series:PS156200', 'HP:gait ataxia', 'HP:0002836 Bladder exstrophy', 'HP: gait ataxia'):         # a name after a prefix is the submitter's, not a spelling
            self.assertEqual(curie(same), same)
        self.assertEqual((curie(None), curie(''), curie(7)), (None, '', 7))

    def test_a_list_keeps_its_order_and_loses_its_repeats(self):
        self.assertEqual(curies(['MONDO:MONDO:1', 'OMIM:2', 'MONDO_1', 'MONDO:1']), ['MONDO:1', 'OMIM:2'])
        self.assertEqual((curies(None), curies([]), curies(['', None])), (None, [], []))

    def test_text_typed_where_the_source_has_identifiers_is_split_off(self):
        self.assertEqual(split_identifiers(['EFO_1', 'gastric adenocarcinoma', 'MONDO:MONDO:2', 'EFO:1']), (['EFO:1', 'MONDO:2'], ['gastric adenocarcinoma']))
        self.assertEqual((split_identifiers(None), split_identifiers([])), (([], []), ([], [])))

    def test_such_text_stays_in_extra_json_beside_what_the_adapter_wrote(self):
        mk = lambda **kw: evidence_row('r', 'opentargets', 'a', 'l', 'raw', 'gene', None, 'ax', 'k', 'c', 'v', 'effect', **kw)
        row = mk(condition_ids=['EFO_1', 'gastric adenocarcinoma'], extra_json=json.dumps({'lead_key': '1:2:A:C'}))
        self.assertEqual(row['condition_ids'], ['EFO:1'])
        self.assertEqual(json.loads(row['extra_json']), {'lead_key': '1:2:A:C', 'condition_ids_not_identifiers': ['gastric adenocarcinoma']})
        self.assertEqual(json.loads(mk(condition_ids=['stomach cancer'])['extra_json']), {'condition_ids_not_identifiers': ['stomach cancer']})
        self.assertEqual(mk(condition_ids=['stomach cancer'])['condition_ids'], [])
        clean = mk(condition_ids=['EFO_1'], extra_json=json.dumps({'k': 1}))
        self.assertEqual(clean['extra_json'], json.dumps({'k': 1}))                                                  # a row with only identifiers keeps its extra_json byte for byte
        self.assertNotIn('extra_json', mk(condition_ids=['EFO_1']))

    def test_every_evidence_row_is_written_with_the_one_spelling(self):
        row = evidence_row('r', 'clinvar', 'a', 'l', 'raw', 'variant', None, 'ax', 'k', 'c', 'v', 'effect', condition_ids=['MONDO:MONDO:0010178', 'MedGen:C0403814', 'MONDO_0005393'])
        self.assertEqual(row['condition_ids'], ['MONDO:0010178', 'MedGen:C0403814', 'MONDO:0005393'])
        self.assertNotIn('condition_ids', evidence_row('r', 's', 'a', 'l', 'raw', 'variant', None, 'ax', 'k', 'c', 'v', 'effect'))                       # nothing is invented for a row without conditions
        self.assertEqual(evidence_row('r', 's', 'a', 'l', 'raw', 'variant', None, 'ax', 'k', 'c', 'v', 'effect', condition_ids=[])['condition_ids'], [])


class TestDatesHaveOneSpelling(unittest.TestCase):
    def test_the_formats_the_sources_use_become_iso_dates(self):
        for written, spelled in (('2023-05-18', '2023-05-18'), ('2023-05-18T16:00:00.000Z', '2023-05-18'), ('2023-05-18T16:00:00', '2023-05-18'), ('2021-01-05-12:30', '2021-01-05'),
                                 ('Mon, 05 Jun 2023 12:00:00 -0400', '2023-06-05'), ('Tue, 31 Oct 2023 23:30:00 -0500', '2023-10-31'), ('Oct 29, 2025', '2025-10-29'), (' 2024-02-29 ', '2024-02-29')):
            self.assertEqual(iso_date(written), spelled, written)

    def test_text_that_is_not_a_date_is_left_alone(self):
        for same in ('2023-13-45', '2023-02-30', 'not evaluated', 'daily export 2026-10-25', '', 'Mon, 99 Foo 2023 12:00:00 -0400'):
            self.assertEqual(iso_date(same), same)
        self.assertEqual((iso_date(None), iso_date(5)), (None, 5))

    def test_every_evidence_row_is_written_with_an_iso_date(self):
        row = evidence_row('r', 'clingen', 'a', 'l', 'raw', 'gene', None, 'ax', 'k', 'c', 'v', 'effect', evaluated='Mon, 05 Jun 2023 12:00:00 -0400')
        self.assertEqual(row['evaluated'], '2023-06-05')
        self.assertNotIn('evaluated', evidence_row('r', 's', 'a', 'l', 'raw', 'variant', None, 'ax', 'k', 'c', 'v', 'effect'))


class TestGuards(unittest.TestCase):
    def test_a_release_of_a_derived_target_carries_the_url_of_the_file_it_was_extracted_from(self):
        # the provenance verifier found this release without any URL: its target is written from another acquired file and has none of its own
        out = tempfile.mkdtemp()
        runner.build(out, [DerivedTargetAdapter()], syn_reference())
        (rel,) = pq.read_table(os.path.join(out, 'source_release.parquet')).to_pylist()
        with open(os.path.join(REPO, 'acquisition', 'targets.json')) as f:
            targets = {t['id']: t for t in json.load(f)}
        t = targets['mavedb_scores_genomic_mapped']
        self.assertIsNone(t['url'])
        self.assertEqual(rel['url'], targets[t['derived_from']]['url'])
        self.assertTrue(rel['url'].startswith('https://'))

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

    def __init__(self, source_id, specs, haplotype_keys=(), link_keys=(), xref_keys=()):
        self.source_id, self.specs, self.haplotype_keys, self.link_keys, self.xref_keys = source_id, specs, haplotype_keys, link_keys, xref_keys

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
        for k in self.link_keys:
            ctx.emit('variant_gene', dict(variant_key=k, gene_id='ENSG1', relation='ot_l2g_lead_variant', consequence=None, release_id=rel))
        for i, k in enumerate(self.xref_keys):
            ctx.emit('xref', dict(namespace='dbsnp', identifier=f'rs{i + 1}', entity_level='variant', entity_id=k, release_id=rel, relation='clinpgx_variant_name'))
        ctx.emit('xref', dict(namespace='hgnc', identifier='HGNC:1', entity_level='gene', entity_id='ENSG1', release_id=rel, relation='x'))             # a gene identifier names no variant
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

    def test_a_variant_that_only_a_gene_link_names_is_an_entity_too(self):
        # real data: 41,504 Open Targets lead variants were linked to genes (L2G) without being entities, so the link pointed at nothing
        out, _ = self.build([('1:20:G:T', '1')], link_keys=('1:20:G:T', '12:44:C:A', 'X:9:A:G'))
        self.assertEqual(self.keys_by_file(out), {'1.parquet': ['1:20:G:T'], '12.parquet': ['12:44:C:A'], 'X.parquet': ['X:9:A:G']})

    def test_a_variant_that_only_an_identifier_names_is_an_entity_too(self):
        # real data: 398 ClinPGx variants had a dbSNP name resolved to a key but no annotation, so 796 xref rows pointed at nothing in entity_variant
        out, _ = self.build([('1:20:G:T', '1')], xref_keys=('1:20:G:T', '10:103209821:A:T', 'X:9373808:C:G'))
        self.assertEqual(self.keys_by_file(out), {'1.parquet': ['1:20:G:T'], '10.parquet': ['10:103209821:A:T'], 'X.parquet': ['X:9373808:C:G']})

    def test_an_index_without_variants_keeps_an_empty_valid_table(self):
        out, _ = self.build([])
        self.assertEqual(self.keys_by_file(out), {'empty.parquet': []})


if __name__ == '__main__':
    unittest.main()
