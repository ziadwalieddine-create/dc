import json
import os
import sys
import tempfile
import unittest

from helpers import syn_reference, FIXTURES
from vindex.build import runner
from vindex.build.events import GnomadSvAdapter, GnomadCnvAdapter, DgvAdapter, Hgsvc3Adapter, info_dict
from vindex.index_reader import Index

sys.path.insert(0, FIXTURES)
import make_sources  # noqa: E402

SOURCES = make_sources.write_all(tempfile.mkdtemp())
EV = os.path.join(SOURCES, 'events')


def adapters():
    return [GnomadSvAdapter(os.path.join(EV, 'gnomad.v4.1.sv.sites.vcf.gz')), GnomadCnvAdapter(os.path.join(EV, 'gnomad.v4.1.cnv.all.vcf.gz')),
            DgvAdapter(os.path.join(EV, 'GRCh38_hg38_variants_2025-12-01.txt')),
            Hgsvc3Adapter([('hgsvc3_grch38_sv_insdel_sym_vcf', os.path.join(EV, 'hgsvc3_insdel.vcf.gz')), ('hgsvc3_grch38_sv_inv_sym_vcf', os.path.join(EV, 'hgsvc3_inv.vcf.gz'))])]


class TestEvents(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.ref = syn_reference()
        cls.out = tempfile.mkdtemp()
        cls.man = runner.build(cls.out, adapters(), cls.ref)
        cls.idx = Index.open(cls.out)
        cls.ev = {r['event_id']: r for r in cls.idx.table('entity_event').to_pylist()}
        cls.stats = {a['adapter']: a['stats'] for a in cls.man['adapters']}

    def test_info_parsing(self):
        self.assertEqual(info_dict('A=1;FLAG;B=x=y'), {'A': '1', 'FLAG': True, 'B': 'x=y'})

    def test_gnomad_sv_events(self):
        d = self.ev['gnomad_sv:gnomAD-SV_v3_DEL_chr1_aaa']
        self.assertEqual((d['event_type'], d['chrom'], d['start'], d['end'], d['svlen'], d['precise'], d['population_af'], d['n_samples'], d['assembly']), ('DEL', '1', 1000, 2000, 1000, True, 0.02, 1000, 'GRCh38'))
        x = json.loads(d['extra_json'])
        self.assertEqual((x['af_by_ancestry'], x['predicted_lof'], x['algorithms'], x['n_homalt']), ({'afr': 0.05, 'nfe': 0.01}, 'GENEA', ['depth', 'manta'], 2))
        dup = self.ev['gnomad_sv:gnomAD-SV_v3_DUP_chr1_bbb']
        self.assertFalse(dup['precise'])                                                  # depth evidence only
        self.assertTrue(json.loads(dup['extra_json'])['low_confidence_repetitive_large_dup'])
        self.assertEqual(json.loads(dup['extra_json'])['filter'], 'HIGH_NCR')

    def test_breakend_and_multiallelic_cnv(self):
        b = self.ev['gnomad_sv:gnomAD-SV_v3_BND_chr1_ccc']
        self.assertEqual((b['event_type'], b['start'], b['end'], b['svlen']), ('BND', 100, 100, None))
        self.assertEqual(json.loads(b['extra_json'])['chr2'], 'chr2')
        m = self.ev['gnomad_sv:gnomAD-SV_v3_CNV_chr1_ddd']
        self.assertEqual((m['event_type'], m['population_af'], m['n_samples'], m['precise'], m['copy_number']), ('CNV', 0.3, 900, False, None))
        self.assertTrue(json.loads(m['extra_json'])['multiallelic'])

    def test_unplaced_contigs_are_counted_not_dropped_silently(self):
        self.assertNotIn('gnomad_sv:gnomAD-SV_v3_INS_chrUn_eee', self.ev)
        self.assertEqual(self.stats['GnomadSvAdapter']['events_on_unplaced_or_unknown_contigs'], 1)
        ins = self.ev['gnomad_sv:gnomAD-SV_v3_INS_chr2_fff']
        self.assertEqual((ins['chrom'], ins['start'], ins['end'], ins['svlen'], ins['precise']), ('2', 300, 300, 300, True))

    def test_gnomad_cnv_confidence_intervals_come_from_observed_ranges(self):
        c = self.ev['gnomad_cnv:variant_is_1_1__DEL']
        self.assertEqual((c['event_type'], c['start'], c['end'], c['start_ci_lo'], c['start_ci_hi'], c['end_ci_lo'], c['end_ci_hi'], c['precise'], c['population_af'], c['n_samples']),
                         ('DEL', 700, 900, -10, 0, 0, 10, False, 0.003, 1000))
        self.assertEqual(json.loads(c['extra_json'])['sf_by_ancestry'], {'nfe': 0.004})

    def test_dgv_merged_variants_only_and_types_from_subtypes(self):
        self.assertEqual(sorted(k for k in self.ev if k.startswith('dgv:')), ['dgv:nsv1', 'dgv:nsv2', 'dgv:nsv3'])         # nssv9 is a sample record, nsv4 has end < start
        self.assertEqual((self.ev['dgv:nsv1']['event_type'], self.ev['dgv:nsv2']['event_type'], self.ev['dgv:nsv3']['event_type']), ('DEL', 'CNV_GAIN', 'UNUSUAL THING'))
        self.assertIsNone(self.ev['dgv:nsv2']['population_af'])                            # DGV merged records state no allele frequency
        self.assertEqual(json.loads(self.ev['dgv:nsv2']['extra_json'])['carrier_fraction'], 0.02)
        x = json.loads(self.ev['dgv:nsv1']['extra_json'])
        self.assertEqual((x['observed_losses'], x['study'], x['pubmed']), (5, 'Study_A', '111'))
        s = self.stats['DgvAdapter']
        self.assertEqual((s['records_not_merged_variants'], s['records_with_unmapped_subtype'], s['events_with_unusable_interval']), (1, 1, 1))

    def test_hgsvc3_allele_frequency_over_called_haplotypes(self):
        d = self.ev['hgsvc3:chr1-1200-DEL-300']
        self.assertEqual((d['event_type'], d['start'], d['end'], d['svlen'], d['precise']), ('DEL', 1200, 1500, 300, True))
        self.assertAlmostEqual(d['population_af'], 2 / 5)                                   # alleles 0|1, 0|0, 1|. -> 2 alt of 5 called
        self.assertEqual(d['n_samples'], 2)
        ins = self.ev['hgsvc3:chr1-1300-INS-150']
        self.assertEqual((ins['event_type'], ins['start'], ins['end'], ins['population_af']), ('INS', 1300, 1300, None))      # no called allele: no frequency, never 0
        inv = self.ev['hgsvc3:chr1-600-INV-200']
        self.assertEqual((inv['event_type'], inv['end'], inv['population_af']), ('INV', 800, 0.5))

    def test_no_evidence_rows_and_idempotent(self):
        self.assertEqual(self.idx.table('evidence').num_rows, 0)
        out2 = tempfile.mkdtemp()
        runner.build(out2, adapters(), self.ref)
        with open(os.path.join(self.out, 'table_digests.json')) as f1, open(os.path.join(out2, 'table_digests.json')) as f2:
            self.assertEqual(json.load(f1), json.load(f2))

    def test_matcher_uses_the_built_events(self):
        from vindex.match_sv import compare, SvPolicy

        class Sv:
            svtype, chrom, start, end, cn, svlen, precise = 'DEL', '1', 1010, 2005, None, None, True
            cipos = ciend = None
        d = self.ev['gnomad_sv:gnomAD-SV_v3_DEL_chr1_aaa']
        m = compare(Sv(), d, SvPolicy())
        self.assertEqual((m.relation, m.equivalent), ('equivalent_breakpoints', True))


if __name__ == '__main__':
    unittest.main()
