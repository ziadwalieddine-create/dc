"""Read evidence through a whole run: the matcher's status of an unresolved locus is never changed by the reads, the read_* columns sit beside it, the manifest accounts for them, and
a run with reads is reproducible. The index and the calls are those of test_wgs; the alignments are synthetic (tests/bamfix.py)."""
import json
import os
import shutil
import tempfile
import unittest

import pyarrow.parquet as pq
import pysam

import bamfix
from bamfix import span
from helpers import SYN_CONTIGS, SYN_FASTA, syn_reference
from test_wgs import b, alt, build_index, sample_vcf
from vindex import reads as R
from vindex import views, wgs
from vindex.index_reader import Index


class TestReadsThroughARun(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.work = tempfile.mkdtemp()
        cls.idx = os.path.join(cls.work, 'idx')
        build_index(cls.idx)
        cls.fa = pysam.FastaFile(SYN_FASTA)
        cls.vcf = sample_vcf(cls.work)
        cls.ref = syn_reference()
        cls.repo = os.path.join(cls.work, 'not_the_repo')
        cls.plain = cls.match(cls, cls.vcf, 'plain')
        cls.bam = cls.make_bam(cls, 'sample.bam', {730: (14, 0), 750: (8, 8), 1050: (6, 6)})
        cls.with_reads = cls.match(cls, cls.vcf, 'reads', reads=cls.bam)

    @classmethod
    def tearDownClass(cls):
        shutil.rmtree(cls.work, ignore_errors=True)

    def make_bam(self, name, plan, mapq=60):
        """plan {position on contig 1: (reference-allele fragments, alternate-allele fragments)}: reads of 24 bases that cover their own position and no other of the plan (the positions
        of the test are 20 bases apart)."""
        reads = []
        for pos, (n_ref, n_alt) in plan.items():
            for i in range(n_ref + n_alt):
                start = pos - 12 + i % 8
                reads.append(span(self.fa, '1', start, start + 23, f'p{pos}_{i}', {pos: alt('1', pos)} if i >= n_ref else {}, mapq=mapq, flag=bamfix.FLAG_REVERSE if i % 2 else 0))
        return bamfix.write_bam(os.path.join(self.work, name), reads, SYN_CONTIGS)

    def match(self, vcf, name, **kw):
        out = os.path.join(self.work, name)
        man = wgs.run_sample(vcf, self.idx, out, self.ref, repo_root=self.repo, **kw)
        return dict(out=out, man=man, un={r['key']: r for r in pq.read_table(os.path.join(out, 'unresolved_loci.parquet')).to_pylist()})

    def key(self, pos):
        return f'1:{pos}:{b("1", pos)}:{alt("1", pos)}'

    def test_the_matcher_status_is_untouched_and_the_read_columns_sit_beside_it(self):
        self.assertEqual(set(self.plain['un']), set(self.with_reads['un']))
        callset_columns = [n for n in wgs.SCHEMAS['unresolved_loci'].names if not n.startswith('read_')]
        for k, r in self.plain['un'].items():
            self.assertEqual({c: r[c] for c in callset_columns}, {c: self.with_reads['un'][k][c] for c in callset_columns}, k)
            self.assertTrue(all(r[c] is None for c in wgs.SCHEMAS['unresolved_loci'].names if c.startswith('read_')), k)           # a run without reads leaves them empty
        self.assertEqual({k.split(':')[1]: r['status'] for k, r in self.with_reads['un'].items()}, {'730': 'ABSENT_FROM_CALLSET', '750': 'NO_CALL', '1050': 'UNRESOLVED'})

    def test_each_locus_states_what_the_reads_say_and_whether_that_settles_the_open_question(self):
        un = self.with_reads['un']
        absent, no_call, inside = un[self.key(730)], un[self.key(750)], un[self.key(1050)]
        self.assertEqual((absent['read_state'], absent['read_depth'], absent['read_ref'], absent['read_alt'], absent['read_settles'], absent['read_interpretation']),
                         (R.REFERENCE_ONLY, 14, 14, 0, True, 'reads_show_the_reference_allele_only'))
        self.assertEqual((no_call['read_state'], no_call['read_ref'], no_call['read_alt'], no_call['read_alt_fraction'], no_call['read_settles'], no_call['read_interpretation']),
                         (R.ALTERNATE_PRESENT, 8, 8, 0.5, True, 'reads_carry_the_allele_the_callset_does_not_list'))
        self.assertEqual((inside['read_state'], inside['read_settles']), (R.ALTERNATE_PRESENT, False))                       # inside a called deletion: haplotype and dosage are open whatever the reads say
        self.assertIn('inside_called_structural_variant', inside['reasons'])
        self.assertTrue(inside['missing_evidence'])                                                                           # the matcher's own statement of what is missing stays
        detail = json.loads(no_call['read_detail_json'])
        self.assertEqual((detail['window'], detail['informative_alignments']), (['1', 750, 750], 16))

    def test_the_manifest_accounts_for_the_reads(self):
        m = self.with_reads['man']
        c = m['counts']['reads']
        self.assertEqual((c['loci_evaluated'], c['by_state']), (3, {R.ALTERNATE_PRESENT: 2, R.REFERENCE_ONLY: 1}))
        self.assertEqual(c['by_status_and_state'], {f'ABSENT_FROM_CALLSET|{R.REFERENCE_ONLY}': 1, f'NO_CALL|{R.ALTERNATE_PRESENT}': 1, f'UNRESOLVED|{R.ALTERNATE_PRESENT}': 1})
        self.assertEqual(c['settles_the_open_question'], {'yes': 2, 'no': 1})
        r = m['input']['reads']
        self.assertEqual((r['sha256'], r['size_bytes'], r['assembly']['ok'], r['policy']), (wgs.sha256_file(self.bam), os.path.getsize(self.bam), True, vars(R.ReadPolicy())))
        self.assertIsNone(self.plain['man']['input']['reads'])
        self.assertIsNone(self.plain['man']['counts']['reads'])
        with open(os.path.join(self.with_reads['out'], 'run_manifest.json')) as f:
            self.assertEqual(json.load(f)['counts']['reads'], c)

    def test_the_report_says_what_the_reads_say_and_a_report_without_reads_has_no_such_section(self):
        texts = {}
        for name, run in (('reads', self.with_reads), ('plain', self.plain)):
            out = os.path.join(self.work, f'report_{name}')
            shutil.copytree(run['out'], out)
            views.write_views(out, Index.open(self.idx))
            with open(os.path.join(out, 'report.md')) as f:
                texts[name] = f.read()
        t = texts['reads']
        self.assertIn('### Read evidence for the unresolved loci', t)
        self.assertIn('Aligned reads were asked about 3 of them', t)
        self.assertIn('the question it leaves open for 2 loci', t)
        self.assertIn('| Callset status | ALTERNATE_PRESENT | REFERENCE_ONLY |', t)
        self.assertIn('| ABSENT_FROM_CALLSET | 0 | 1 |', t)
        self.assertIn('| NO_CALL | 1 | 0 |', t)
        self.assertIn('reads show the reference allele only 1', t)
        self.assertIn('| Locus | Status | Rules | Missing evidence | Reads |', t)
        self.assertIn('ALTERNATE_PRESENT (8/16 alternate)', t)
        self.assertNotIn('Read evidence', texts['plain'])
        self.assertNotIn('| Reads |', texts['plain'])

    def test_a_filtered_call_is_supported_or_contradicted_by_the_reads(self):
        extra = [('1', 730, b('1', 730), alt('1', 730), 'LowQual', 50, ['0/1:60:30:15,15'])]
        vcf = sample_vcf(tempfile.mkdtemp(dir=self.work), extra=extra)                                                       # its own directory: sample_vcf writes sample.vcf.gz
        for plan, state, meaning in (({730: (8, 8)}, R.ALTERNATE_PRESENT, 'reads_support_the_filtered_call'), ({730: (14, 0)}, R.REFERENCE_ONLY, 'reads_contradict_the_filtered_call')):
            run = self.match(vcf, f'filtered_{state}', reads=self.make_bam(f'f_{state}.bam', plan))
            row = run['un'][self.key(730)]
            self.assertEqual((row['status'], row['read_state'], row['read_interpretation'], row['read_settles']), ('FILTERED_CALL', state, meaning, True))

    def test_two_runs_with_reads_give_the_same_tables(self):
        again = self.match(self.vcf, 'reads_again', reads=self.bam)
        a = pq.read_table(os.path.join(self.with_reads['out'], 'unresolved_loci.parquet'))
        b_ = pq.read_table(os.path.join(again['out'], 'unresolved_loci.parquet'))
        self.assertTrue(a.equals(b_))
        self.assertEqual(self.with_reads['man']['counts']['reads'], again['man']['counts']['reads'])

    def test_reads_below_the_policy_leave_the_locus_unsettled_and_say_why(self):
        bam = self.make_bam('low.bam', {730: (14, 0), 750: (14, 0), 1050: (14, 0)}, mapq=3)
        run = self.match(self.vcf, 'lowmapq', reads=bam)
        row = run['un'][self.key(730)]
        self.assertEqual((row['read_state'], row['read_depth'], row['read_reasons'], row['read_settles']), (R.INSUFFICIENT, 0, ['no_informative_fragments'], False))
        self.assertEqual(json.loads(row['read_detail_json'])['filtered'], {'low_mapq': 14})

    def test_a_policy_is_applied_and_recorded(self):
        run = self.match(self.vcf, 'strict', reads=self.bam, read_policy=R.ReadPolicy(min_depth=15))
        self.assertEqual(run['un'][self.key(730)]['read_state'], R.INSUFFICIENT)                                              # 14 informative fragments are fewer than 15
        self.assertEqual(run['man']['input']['reads']['policy']['min_depth'], 15)

    def test_the_command_line_passes_the_alignment_and_the_policy_through(self):
        import contextlib
        import io
        import match_wgs
        report = os.path.join(self.work, 'empty_assembly_report.txt')
        open(report, 'w').close()
        out = os.path.join(self.work, 'cli')
        with contextlib.redirect_stdout(io.StringIO()) as printed:
            rc = match_wgs.main(['--vcf', self.vcf, '--out', out, '--index', self.idx, '--fasta', SYN_FASTA, '--assembly-report', report, '--allow-assembly-mismatch',
                                 '--reads', self.bam, '--read-min-depth', '15', '--read-min-alt-reads', '4'])
        self.assertEqual(rc, 0)
        with open(os.path.join(out, 'run_manifest.json')) as f:
            m = json.load(f)
        self.assertEqual((m['input']['reads']['policy']['min_depth'], m['input']['reads']['policy']['min_alt_reads']), (15, 4))
        un = {r['key']: r for r in pq.read_table(os.path.join(out, 'unresolved_loci.parquet')).to_pylist()}
        self.assertEqual(un[self.key(730)]['read_state'], R.INSUFFICIENT)                                                     # 14 fragments, 15 asked for
        self.assertIn('"reads"', printed.getvalue())

    def test_an_index_kept_elsewhere_is_passed_through_and_recorded(self):
        bam = self.make_bam('moved.bam', {730: (14, 0), 750: (8, 8), 1050: (6, 6)})
        moved = os.path.join(self.work, 'elsewhere.bai')
        os.rename(bam + '.bai', moved)
        with self.assertRaises(OSError):
            self.match(self.vcf, 'noindex', reads=bam)
        run = self.match(self.vcf, 'withindex', reads=bam, reads_index=moved)
        self.assertEqual(run['man']['input']['reads']['index'], moved)
        self.assertEqual({k: r['read_state'] for k, r in run['un'].items()}, {k: r['read_state'] for k, r in self.with_reads['un'].items()})

    def test_an_alignment_of_another_assembly_stops_the_run_before_anything_is_written(self):
        reads = [span(self.fa, '1', 700, 739, f'x{i}') for i in range(3)]
        bam = bamfix.write_bam(os.path.join(self.work, 'wrong.bam'), reads, [('1', 3001)] + SYN_CONTIGS[1:])
        out = os.path.join(self.work, 'refused')
        with self.assertRaises(R.ReadsAssemblyMismatch):
            wgs.run_sample(self.vcf, self.idx, out, self.ref, repo_root=self.repo, reads=bam)
        self.assertFalse(os.path.exists(os.path.join(out, 'unresolved_loci.parquet')))
        self.assertFalse(os.path.exists(os.path.join(out, 'run_manifest.json')))


if __name__ == '__main__':
    unittest.main()
