"""End to end on synthetic data: sources -> index -> a sample's calls -> per-chromosome matching -> result tables and views."""
import json
import os
import shutil
import subprocess
import sys
import tempfile
import unittest

import pysam
import pyarrow.parquet as pq

from helpers import syn_reference, vcf_text, FIXTURES, SYN_CONTIGS
from vindex import views, wgs
from vindex.build import runner
from vindex.build.clinvar_vcf import ClinVarVcfAdapter
from vindex.build.clinvar_vcv import ClinVarVcvAdapter
from vindex.build.genes import GeneMap
from vindex.build.gnomad import GnomadConstraintAdapter
from vindex.build.pathways import ReactomeAdapter
from vindex.build.pgx import CpicAdapter
from vindex.build.regions import GiabRegionAdapter, StrchiveAdapter
from vindex.index_reader import Index

sys.path.insert(0, FIXTURES)
import make_sources  # noqa: E402

SOURCES = make_sources.write_all(tempfile.mkdtemp())
b, alt = make_sources.b, make_sources.alt
ROOT = os.path.join(SOURCES, 'regions')


def combined_clinvar_vcf(path):
    """The mini ClinVar VCF plus records that need reference evidence to settle: reviewed pathogenic alleles the sample does and does not cover."""
    src = os.path.join(SOURCES, 'mini_clinvar.vcf')
    with open(src) as f:
        lines = f.read().rstrip('\n').split('\n')

    def rec(pos, vid, sig, rev):
        return '\t'.join(['1', str(pos), str(vid), b('1', pos), alt('1', pos), '.', '.', f'ALLELEID={vid + 1000};CLNDISDB=MedGen:C9;CLNDN=Extra;CLNREVSTAT={rev};CLNSIG={sig};CLNSIGSCV=SCV9;CLNVC=single_nucleotide_variant;ORIGIN=1'])
    lines += [rec(710, 2101, 'Pathogenic', 'criteria_provided,_multiple_submitters,_no_conflicts'), rec(730, 2102, 'Likely_pathogenic', 'criteria_provided,_multiple_submitters,_no_conflicts'),
              rec(750, 2103, 'Pathogenic', 'reviewed_by_expert_panel'), rec(1050, 2104, 'Pathogenic', 'criteria_provided,_multiple_submitters,_no_conflicts')]
    with open(path, 'w') as f:
        f.write('\n'.join(lines) + '\n')


class StubGenes:
    """Gene entities the synthetic sources do not provide (the real ones come from Open Targets)."""
    source_id = 'opentargets'

    def run(self, ctx, limit=None):
        ctx.adapter = self
        rel = ctx.register_release('opentargets', 'opentargets_26_09_target', '26.09')
        ctx.emit('entity_gene', dict(gene_id='ENSG_G1', symbol='GENE1', hgnc_id='', chrom='1', start=395, end=445, strand='+', mane_select='', biotype='protein_coding', release_id=rel), name='stub_genes')
        return {}


class Fake(GeneMap):
    pass


class FakeGenes:
    by_ensg = {'ENSG_G1': 1, 'ENSG_G2': 1}
    syms = {'TESTG': 'ENSG_TESTG', 'TESTH': 'ENSG_TESTH'}

    def resolve(self, hgnc=None, symbol=None, ncbi=None):
        return (self.syms[symbol], 'symbol') if symbol in self.syms else (None, 'unresolved')

    def from_ncbi(self, x):
        return None

    def from_symbol(self, s):
        return None


def sample_vcf(tmp, extra=()):
    wrong = 'A' if b('1', 440) != 'A' else 'C'
    gt = lambda g: [f'{g}:60:30:15,15']
    recs = [('1', 105, 'AA', 'A', 'PASS', 50, gt('1/1')), ('1', 150, b('1', 150), alt('1', 150), 'PASS', 50, gt('0/1')), ('1', 402, b('1', 402), alt('1', 402), 'PASS', 50, gt('0/1')),
            ('1', 410, b('1', 410), alt('1', 410), 'LowQual', 50, gt('0/1')), ('1', 430, b('1', 430), alt('1', 430), 'PASS', 50, gt('0/1')),
            ('1', 440, wrong, 'T' if wrong != 'T' else 'G', 'PASS', 50, gt('0/1')), ('1', 600, b('1', 600), alt('1', 600), 'PASS', 50, gt('0/1')),
            ('1', 700, b('1', 700), '<NON_REF>', 'PASS', 50, ['0/0:60:30:30,0'], 'END=720'), ('1', 750, b('1', 750), alt('1', 750), 'PASS', 50, ['./.:.:.:.']),
            ('1', 1000, b('1', 1000), '<DEL>', 'PASS', 50, gt('0/1'), 'SVTYPE=DEL;END=1099;SVLEN=-100'), ('1', 1500, b('1', 1500), '<DUP>', 'PASS', 50, gt('0/1'), 'SVTYPE=DUP;END=1799;SVLEN=300'),
            ('2', 100, b('2', 100), alt('2', 100), 'PASS', 50, gt('0/1')), ('X', 50, b('X', 50), alt('X', 50), 'PASS', 50, ['1:60:30:0,30'])] + list(extra)
    order = {c: i for i, (c, _) in enumerate(SYN_CONTIGS)}
    recs.sort(key=lambda r: (order[r[0]], r[1]))                       # tabix needs coordinate order; extras may fall on any contig
    p = os.path.join(tmp, 'sample.vcf')
    with open(p, 'w') as f:
        f.write(vcf_text(recs))
    pysam.tabix_index(p, preset='vcf', force=True)
    return p + '.gz'


def build_index(out):
    ref = syn_reference()
    genes = GeneMap(ref, ot_dir='/nonexistent', mane_path='/nonexistent')
    cv = os.path.join(os.path.dirname(out), 'cv.vcf')
    combined_clinvar_vcf(cv)

    class Vcf(ClinVarVcfAdapter):
        def run(self, ctx, limit=None):
            return super().run(ctx, path=cv, limit=limit)
    runner.build(out, [StubGenes(), Vcf(genes), ClinVarVcvAdapter(genes, root=os.path.join(SOURCES, 'vcv')), GiabRegionAdapter(root=ROOT), StrchiveAdapter(genes, root=os.path.join(ROOT, 'strchive')),
                       ReactomeAdapter(FakeGenes(), root=os.path.join(SOURCES, 'reactome')),
                       GnomadConstraintAdapter(FakeGenes(), os.path.join(SOURCES, 'gnomad', 'gnomad.v4.1.2.constraint_metrics.tsv.bgz')),
                       CpicAdapter(FakeGenes(), root=os.path.join(SOURCES, 'cpic'))], ref)


class TestWholeGenomeRun(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.work = tempfile.mkdtemp()
        cls.idx_dir = os.path.join(cls.work, 'idx')
        build_index(cls.idx_dir)
        cls.vcf = sample_vcf(cls.work)
        cls.out = os.path.join(cls.work, 'res')
        cls.man = wgs.run_sample(cls.vcf, cls.idx_dir, cls.out, syn_reference(), repo_root=os.path.join(cls.work, 'not_the_repo'))
        cls.t = lambda self, name: pq.read_table(os.path.join(self.out, f'{name}.parquet')).to_pylist()
        cls.calls = {r['key']: r for r in cls.t(cls, 'small_calls')}

    @classmethod
    def tearDownClass(cls):
        shutil.rmtree(cls.work, ignore_errors=True)

    def k(self, c, pos, ref=None, a=None):
        return f"{c}:{pos}:{ref or b(c, pos)}:{a or alt(c, pos)}"

    def test_every_call_is_accounted_for_with_its_status(self):
        c = self.man['counts']
        self.assertEqual((c['calls'], c['by_status']), (8, {'EXACT_MATCH': 3, 'FILTERED_CALL': 1, 'NO_INDEX_RECORD': 4}))
        self.assertEqual(len(self.calls), 8)
        self.assertEqual(self.calls[self.k('1', 402)]['status'], 'EXACT_MATCH')
        self.assertEqual(self.calls[self.k('1', 410)]['status'], 'FILTERED_CALL')               # a filtered call is not a confirmed presence, and it is still listed
        self.assertEqual(self.calls[self.k('1', 600)]['status'], 'NO_INDEX_RECORD')
        self.assertEqual(self.calls['1:101:CA:C']['status'], 'EXACT_MATCH')                       # the right-shifted deletion matches by normalised identity
        x = self.calls[self.k('X', 50)]
        self.assertEqual((x['zygosity'], x['ploidy'], x['status']), ('hemizygous', 1, 'NO_INDEX_RECORD'))
        self.assertEqual(self.calls[self.k('1', 402)]['zygosity'], 'heterozygous')

    def test_rejected_records_are_listed_with_the_reason(self):
        rej = self.t('rejected_calls')
        self.assertEqual([(r['chrom'], r['pos'], r['reason']) for r in rej], [('1', 440, 'ref_mismatch')])

    def test_matched_calls_carry_per_axis_summaries_and_links(self):
        c = self.calls[self.k('1', 402)]
        axes = json.loads(c['axes_json'])
        self.assertEqual(axes['clinical_validity_actionability']['assertions'], {'effect': 3})                          # aggregate + two submissions
        self.assertEqual(axes['clinical_validity_actionability']['n_independent'], 1)                                       # one ClinVar lineage, counted once
        self.assertEqual(axes['clinical_validity_actionability']['max_review_stars'], 3)
        self.assertEqual(c['n_evidence'], len(c['evidence_ids']))
        links = [r for r in self.t('evidence_links') if r['call_id'] == self.k('1', 402)]
        self.assertEqual(sorted(r['evidence_id'] for r in links), c['evidence_ids'])
        self.assertTrue(all(r['relation'] == 'exact_variant' and r['level'] == 'variant' for r in links))
        conflict = json.loads(self.calls[self.k('1', 430)]['axes_json'])['clinical_validity_actionability']
        self.assertTrue(conflict['conflict'] or 'conflicting' in conflict['assertions'])

    def test_difficult_region_is_flagged_not_filtered(self):
        c = self.calls[self.k('1', 150)]
        self.assertIn('difficult_region', c['flags'])
        self.assertEqual(c['status'], 'NO_INDEX_RECORD')
        self.assertIn('giab_difficult_segdups', c['region_kinds'])
        self.assertNotIn('difficult_region', self.calls[self.k('1', 600)]['flags'])

    def test_gene_level_evidence_is_context_not_variant_evidence(self):
        ctx = self.t('gene_context')
        self.assertTrue(ctx)
        self.assertEqual({r['applies_as'] for r in ctx}, {'gene_level_context_only'})
        g1 = [r for r in ctx if r['gene_id'] == 'ENSG_G1']
        self.assertEqual({r['source_id'] for r in g1}, {'gnomad_constraint'})
        self.assertEqual({r['n_matched_calls'] for r in g1}, {2})                                  # 402 and 430 sit in the gene with index evidence and a confirmed call; 410 is filtered
        self.assertEqual(g1[0]['symbol'], 'GENE1')
        gene_ids = {r['evidence_id'] for r in g1}
        self.assertFalse(gene_ids & {e for c in self.calls.values() for e in (c['evidence_ids'] or ())})       # never attached to a variant call

    def test_structural_calls_equivalent_only_when_type_breakpoints_and_dosage_agree(self):
        sv = {r['svtype']: r for r in self.t('sv_calls')}
        d = sv['DEL']
        self.assertEqual((d['n_equivalent'], d['equivalent_event_ids']), (1, ['clinvar:VCV000002001']))
        self.assertEqual(len(d['event_evidence_ids']), 2)                                          # the aggregate and one submission, both at event level
        links = [r for r in self.t('evidence_links') if r['call_kind'] == 'sv']
        self.assertEqual({r['relation'] for r in links}, {'equivalent_event'})
        u = sv['DUP']
        self.assertEqual(u['n_equivalent'], 0)                                                    # same type and breakpoints, but the index states a copy number the call does not
        self.assertIn('dosage_unresolved', u['refused_json'])

    def test_unresolved_loci_state_the_exact_missing_evidence(self):
        un = {r['key']: r for r in self.t('unresolved_loci')}
        self.assertEqual({k.split(':')[1]: r['status'] for k, r in un.items()}, {'730': 'ABSENT_FROM_CALLSET', '750': 'NO_CALL', '1050': 'UNRESOLVED'})
        self.assertTrue(all(r['rules'] == ['clinvar_pathogenic_reviewed'] for r in un.values()))
        self.assertTrue(all((r['missing_evidence'] or '').strip() for r in un.values()))           # every unresolved locus says what would settle it
        self.assertIn('genotype at this site', un[self.k('1', 750)]['missing_evidence'])
        absent = next(r for r in un.values() if r['pos'] == 730)
        self.assertIn('reference block', absent['missing_evidence'])
        inside_sv = next(r for r in un.values() if r['pos'] == 1050)
        self.assertIn('inside_called_structural_variant', inside_sv['reasons'])
        rules = self.man['counts']['high_value_rules']['clinvar_pathogenic_reviewed']
        self.assertEqual(rules, {'EXACT_MATCH': 1, 'VERIFIED_REFERENCE': 1, 'ABSENT_FROM_CALLSET': 1, 'NO_CALL': 1, 'UNRESOLVED': 1})          # 402 carried; 710 covered by the reference block

    def test_haplotype_and_repeat_calls_say_what_is_missing(self):
        hap = {r['gene']: r for r in self.t('haplotype_calls')}
        self.assertEqual(set(hap), {'TESTG', 'TESTH', 'TESTK', 'NOCHR', 'MTG'})
        self.assertEqual(hap['TESTG']['status'], 'RESOLVED')                                       # both defining positions lie in the reference block
        self.assertEqual(json.loads(hap['TESTG']['diplotypes_json']), [['*1', '*1']])
        self.assertEqual(hap['MTG']['status'], 'UNRESOLVED')                                       # nothing known on the mitochondrial contig
        rep = self.t('repeat_calls')
        self.assertEqual(sorted(r['locus_id'] for r in rep), ['TEST1_GENEA', 'TEST2_GENEB', 'TEST3_MULTI', 'TEST5_NOJSON'])
        self.assertTrue(all(r['status'] == 'UNRESOLVED' and 'repeat-length genotype' in r['missing'] for r in rep))

    def test_repeat_genotypes_at_loci_the_index_does_not_know_are_kept(self):
        rec = ('2', 700, 'N', '<STR10>,<STR12>', 'PASS', 0, ['1/2:10/12:10-10/11-13'], 'END=720;RU=CAG;VARID=UNLISTED1', '.', 'GT:REPCN:REPCI')
        p = os.path.join(self.work, 'rep.vcf')
        with open(p, 'w') as f:
            f.write(vcf_text([rec]))
        pysam.tabix_index(p, preset='vcf', force=True)
        out = os.path.join(self.work, 'res_rep')
        man = wgs.run_sample(p + '.gz', self.idx_dir, out, syn_reference(), repo_root=os.path.join(self.work, 'not_the_repo'))
        rep = {r['locus_id']: r for r in pq.read_table(os.path.join(out, 'repeat_calls.parquet')).to_pylist()}
        self.assertEqual({k for k, r in rep.items() if r['status'] == 'UNRESOLVED'}, {'TEST1_GENEA', 'TEST2_GENEB', 'TEST3_MULTI', 'TEST5_NOJSON'})
        u = rep['UNLISTED1']
        self.assertEqual((u['status'], u['counts'], u['classes'], u['chrom'], u['start'], u['end'], u['motif'], u['caller']), ('NO_INDEX_RECORD', [10, 12], [], '2', 700, 720, 'CAG', 'expansionhunter'))
        self.assertEqual(man['counts']['repeat_status'], {'UNRESOLVED': 4, 'NO_INDEX_RECORD': 1})

    def test_manifest_records_input_index_and_code(self):
        m = self.man
        self.assertEqual(m['input']['sha256'], wgs.sha256_file(self.vcf))
        self.assertTrue(m['input']['assembly']['ok'])
        self.assertEqual(m['index']['schema_version'], m['schema_version'])
        self.assertTrue(m['index']['table_digests.json_sha256'] and m['index']['index_manifest.json_sha256'])
        self.assertEqual(m['rows']['small_calls'], 8)
        self.assertEqual(set(m['rules']), {r.name for r in views.HIGH_VALUE_RULES})
        with open(os.path.join(self.out, 'run_manifest.json')) as f:
            self.assertEqual(json.load(f)['rows'], m['rows'])

    def test_views_filter_and_order_without_dropping_calls(self):
        rows = list(views.read_calls(self.out))
        v = views.build_views(rows)
        clin = [r['key'] for r in v['clinical_effect']['rows']]
        self.assertEqual(clin[0], self.k('1', 402))                                                # three review stars first
        self.assertIn('1:101:CA:C', clin)
        self.assertIn(self.k('1', 430), [r['key'] for r in v['conflicting_evidence']['rows']])
        self.assertTrue(v['clinical_effect']['ordering'])
        self.assertEqual(len(list(views.read_calls(self.out, statuses=None))), 8)                      # the full table is untouched by any view
        traits = views.trait_view(rows)
        self.assertIn('MedGen:C0000001', traits)
        self.assertIn(self.k('1', 402), traits['MedGen:C0000001']['calls'])
        members = Index.open(self.idx_dir).table('gene_set_member').to_pylist()
        paths = views.pathway_view(rows, members)
        self.assertEqual({p: e['n_calls'] for p, e in paths.items()}, {'R-HSA-1': 2, 'R-HSA-2': 2})
        self.assertEqual(paths['R-HSA-1']['genes'], ['ENSG_G1'])

    def test_written_views_list_rows_of_the_full_table_and_say_how_they_are_ordered(self):
        out = os.path.join(self.work, 'viewres')
        shutil.copytree(self.out, out)
        m = views.write_views(out, Index.open(self.idx_dir))
        self.assertEqual(m['n_calls_considered'], 3)
        self.assertEqual(m['views']['clinical_effect']['n_calls'], 3)
        self.assertTrue(all(v['ordering'] for v in m['views'].values()))
        t = pq.read_table(os.path.join(out, 'views', 'clinical_effect.parquet')).to_pylist()
        self.assertEqual(t[0]['key'], self.k('1', 402))
        full = {r['key'] for r in pq.read_table(os.path.join(out, 'small_calls.parquet')).to_pylist()}
        self.assertTrue({r['key'] for r in t} <= full and len(full) == 8)
        with open(os.path.join(out, 'views', 'pathway_view.json')) as f:
            self.assertEqual(sorted(json.load(f)), ['R-HSA-1', 'R-HSA-2'])
        with open(os.path.join(out, 'report.md')) as f:
            rep = f.read()
        self.assertIn(self.k('1', 402), rep)                                                        # the three-star pathogenic call heads the clinical list
        self.assertIn('Ordering: review stars (highest first)', rep)
        self.assertIn('absence from the index, not evidence of no effect', rep)
        self.assertIn('| `1:750:', rep)                                                              # an unresolved locus with its missing evidence
        self.assertIn('NO_CALL', rep)
        self.assertIn('lineage', rep)

    def test_the_catalogue_view_lists_the_high_value_variants_without_a_sample(self):
        out = os.path.join(self.work, 'catview')
        m = views.write_catalogue_views(Index.open(self.idx_dir), out)
        self.assertEqual(m['variants'], 5)                                                          # 402 (3 stars) and the four reviewed alleles of the fixture
        self.assertEqual(m['per_rule'], {'clinvar_pathogenic_reviewed': 5})
        t = pq.read_table(os.path.join(out, 'high_value_variants', '1.parquet')).to_pylist()
        self.assertEqual([r['pos'] for r in t], sorted(r['pos'] for r in t))
        r = next(x for x in t if x['key'] == self.k('1', 402))
        axes = json.loads(r['axes_json'])
        self.assertEqual(axes['clinical_validity_actionability']['max_review_stars'], 3)
        self.assertEqual(r['n_evidence'], len(r['evidence_ids']))
        self.assertTrue(r['gene_ids'] and all(g.startswith('ENSG') or ':' in g for g in r['gene_ids']))
        self.assertIn('nothing here is a score or a tier', m['note'])

    def test_run_is_deterministic(self):
        out2 = os.path.join(self.work, 'res2')
        wgs.run_sample(self.vcf, self.idx_dir, out2, syn_reference(), repo_root=os.path.join(self.work, 'not_the_repo'))
        for name in wgs.SCHEMAS:
            self.assertEqual(pq.read_table(os.path.join(self.out, f'{name}.parquet')), pq.read_table(os.path.join(out2, f'{name}.parquet')), name)

    def test_results_do_not_depend_on_the_hash_seed(self):
        """Separate processes with different PYTHONHASHSEED values: an output that depended on set order would differ."""
        cli = os.path.join(os.path.dirname(__file__), '..', 'match_wgs.py')
        outs = []
        for seed in ('1', '2'):
            out = os.path.join(self.work, f'seed{seed}')
            r = subprocess.run([sys.executable, cli, '--vcf', self.vcf, '--index', self.idx_dir, '--out', out, '--fasta', os.path.join(FIXTURES, 'syn38.fa'),
                                '--assembly-report', os.path.join(SOURCES, 'assembly_report.txt'), '--no-hash', '--allow-assembly-mismatch'],
                               env=dict(os.environ, PYTHONHASHSEED=seed), capture_output=True, text=True)
            self.assertEqual(r.returncode, 0, r.stderr[-1500:])
            outs.append(out)
        for name in wgs.SCHEMAS:
            self.assertEqual(pq.read_table(os.path.join(outs[0], f'{name}.parquet')), pq.read_table(os.path.join(outs[1], f'{name}.parquet')), name)

    def test_personal_results_cannot_be_written_into_the_tracked_tree(self):
        repo = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..'))
        with self.assertRaises(ValueError):
            wgs.check_output_location(os.path.join(repo, 'variant_index', 'results'), repo)
        wgs.check_output_location(os.path.join(repo, 'data', 'results', 'x'), repo)
        wgs.check_output_location(self.out, repo)

    def test_a_different_assembly_is_refused(self):
        tmp = tempfile.mkdtemp()
        try:
            p = os.path.join(tmp, 'wrong.vcf')
            with open(p, 'w') as f:
                f.write(vcf_text([('1', 402, b('1', 402), alt('1', 402), 'PASS', 50, ['0/1:60:30:15,15'])], contigs=[('1', 248956422), ('2', 242193529), ('X', 156040895), ('MT', 16569)]))
            with self.assertRaises(wgs.AssemblyMismatch):
                wgs.run_sample(p, self.idx_dir, os.path.join(tmp, 'o'), syn_reference(), repo_root=os.path.join(tmp, 'x'))
        finally:
            shutil.rmtree(tmp, ignore_errors=True)


if __name__ == '__main__':
    unittest.main()
