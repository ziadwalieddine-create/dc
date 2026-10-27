"""Records as real callsets hold them: lowercase and N bases, a REF that is not the reference, '*' alleles, gVCF blocks, an END before the start, breakends, genotypes that are missing
or half missing, contigs the reference does not know and contigs written with a 'chr' prefix. A whole-genome run must not stop on any of them, must account for every record it
read, and must give the same tables when it is run twice."""
import os
import random
import shutil
import sys
import tempfile
import unittest

import pyarrow.parquet as pq
import pysam

from helpers import syn_reference, vcf_text, SYN_CONTIGS, FIXTURES
from vindex import wgs
from test_wgs import build_index

sys.path.insert(0, FIXTURES)
import make_sources  # noqa: E402

REF = make_sources.REF
LENGTHS = dict(SYN_CONTIGS)
EXTRA_CONTIGS = [('chrUn_decoy1', 800), ('HLA-A*01:01:01:01', 600), ('chr1', 3000)]
ALL_CONTIGS = SYN_CONTIGS + EXTRA_CONTIGS
GENOTYPES = ['0/1', '1/1', '0|1', '1|0', '1|1', './.', '.', '0/0', '1', '0', './1', '1/.', '0/.', '0|0']
MULTI_GENOTYPES = ['1/2', '0/2', '2/2', '1|2', '2|1', '0/1', './2', '2/.']
KINDS = ['snv', 'snv', 'snv', 'wrong_ref', 'lower', 'n_ref', 'n_alt', 'ins', 'del', 'long_del', 'mnv', 'multi', 'star', 'block', 'dot_alt', 'non_ref_alt', 'sv_del', 'sv_dup', 'sv_inv',
         'sv_ins', 'sv_cnv', 'sv_bad_end', 'bnd', 'unknown_contig', 'prefixed', 'no_format', 'missing_sample', 'edge_start', 'edge_end', 'duplicate']


def span(c, pos, n):
    return REF.fetch(c, pos - 1, pos - 1 + n).upper()


def other(base, rng):
    return rng.choice([x for x in 'ACGT' if x != base])


def record(rng, kind):
    """(chrom, pos, ref, alt, filter, qual, sample fields, info, id, format) for one odd record of the stated kind."""
    c = rng.choices(['1', '2', 'X', 'MT'], weights=[6, 2, 1, 1])[0]
    n = LENGTHS[c]
    pos = rng.randint(1, n - 70)
    gt = rng.choice(GENOTYPES)
    ref, alt, info, fmt = span(c, pos, 1), None, '.', 'GT:GQ:DP:AD'
    filt, qual = rng.choice(['PASS', 'PASS', '.', 'LowQual']), rng.choice(['50', '.', '12.5'])
    if kind in ('snv', 'duplicate'):
        alt = other(ref, rng)
    elif kind == 'wrong_ref':
        ref = other(ref, rng)
        alt = other(ref, rng)
    elif kind == 'lower':
        ref, alt = ref.lower(), other(ref, rng).lower()
    elif kind == 'n_ref':
        ref, alt = 'N', rng.choice('ACGT')
    elif kind == 'n_alt':
        alt = 'N'
    elif kind == 'ins':
        alt = ref + ''.join(rng.choice('ACGT') for _ in range(rng.randint(1, 8)))
    elif kind == 'del':
        k = rng.randint(2, 6)
        ref, alt = span(c, pos, k), span(c, pos, 1)
    elif kind == 'long_del':
        ref, alt = span(c, pos, 60), span(c, pos, 1)
    elif kind == 'mnv':
        k = rng.randint(2, 4)
        ref = span(c, pos, k)
        alt = ''.join(other(x, rng) for x in ref)
    elif kind == 'multi':
        alts = list({other(ref, rng), other(ref, rng), ref + rng.choice('ACGT')})
        alt, gt = ','.join(alts), rng.choice(MULTI_GENOTYPES[:6])
        if len(alts) < 2:
            gt = '0/1'
    elif kind == 'star':
        alt, gt = other(ref, rng) + ',*', rng.choice(['1/2', '2/2', '0/2', '2|1'])
    elif kind == 'block':
        alt, gt, info = '<NON_REF>', '0/0', f'END={min(pos + rng.randint(0, 40), n)}'
    elif kind == 'dot_alt':
        alt, gt, info = '.', '0/0', f'END={min(pos + rng.randint(0, 40), n)}'
    elif kind == 'non_ref_alt':
        alt = other(ref, rng) + ',<NON_REF>'
    elif kind in ('sv_del', 'sv_dup', 'sv_inv', 'sv_cnv'):
        t = {'sv_del': 'DEL', 'sv_dup': 'DUP', 'sv_inv': 'INV', 'sv_cnv': 'CNV'}[kind]
        end = pos + rng.randint(5, 60)
        alt, info = f'<{t}>', f'SVTYPE={t};END={end};SVLEN={-(end - pos) if t == "DEL" else end - pos}'
        if kind == 'sv_cnv':
            fmt, gt = 'GT:CN', rng.choice(['0/1', '1/1', './.'])
    elif kind == 'sv_ins':
        alt, info = '<INS>', f'SVTYPE=INS;SVLEN={rng.randint(5, 300)}'
    elif kind == 'sv_bad_end':
        alt, info = '<DEL>', f'SVTYPE=DEL;END={max(pos - rng.randint(1, 5), 1)}'
    elif kind == 'bnd':
        alt, info = f'{ref}]2:{rng.randint(1, 1000)}]', f'SVTYPE=BND;MATEID=m{rng.randint(1, 99)}'
    elif kind == 'unknown_contig':
        c = rng.choice(['chrUn_decoy1', 'HLA-A*01:01:01:01'])
        pos = rng.randint(1, 500)
        ref, alt = 'A', 'G'                                                      # the reference does not hold this contig: the run cannot know the base
    elif kind == 'prefixed':
        c, ref, alt = 'chr1', span('1', pos, 1), other(span('1', pos, 1), rng)
    elif kind == 'no_format':
        alt, fmt = other(ref, rng), '.'
    elif kind == 'missing_sample':
        alt = other(ref, rng)
    elif kind == 'edge_start':
        pos, alt = 1, other(span(c, 1, 1), rng)
        ref = span(c, 1, 1)
    elif kind == 'edge_end':
        pos = n
        ref = span(c, n, 1)
        alt = other(ref, rng)
    sample = f'{gt}:60:30:15,15' if fmt == 'GT:GQ:DP:AD' else (f'{gt}:2' if fmt == 'GT:CN' else '.')
    if kind == 'missing_sample':
        sample = '.'
    elif kind == 'multi' or kind == 'star':
        sample = f'{gt}:60:30:10,10,10'
    return (c, pos, ref, alt, filt, qual, [sample], info, '.', fmt)


def sample_vcf_from(rng, n, tmp):
    recs = [record(rng, rng.choice(KINDS)) for _ in range(n)]
    recs += [recs[0]]                                                      # a record repeated exactly
    order = {c: i for i, (c, _) in enumerate(ALL_CONTIGS)}
    recs.sort(key=lambda r: (order[r[0]], r[1]))
    p = os.path.join(tmp, 'odd.vcf')
    with open(p, 'w') as f:
        f.write(vcf_text(recs, contigs=ALL_CONTIGS))
    pysam.tabix_index(p, preset='vcf', force=True)
    return p + '.gz', len(recs)


class TestOddRecords(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.work = tempfile.mkdtemp()
        cls.idx = os.path.join(cls.work, 'idx')
        build_index(cls.idx)

    @classmethod
    def tearDownClass(cls):
        shutil.rmtree(cls.work, ignore_errors=True)

    def run_once(self, vcf, tag):
        out = os.path.join(self.work, tag)
        man = wgs.run_sample(vcf, self.idx, out, syn_reference(), repo_root=os.path.join(self.work, 'not_the_repo'))
        return out, man

    def test_a_run_over_odd_records_reads_every_record_and_accounts_for_every_call(self):
        for seed in range(6):
            with self.subTest(seed=seed):
                rng = random.Random(seed)
                tmp = os.path.join(self.work, f'in{seed}')
                os.makedirs(tmp)
                vcf, n = sample_vcf_from(rng, 160, tmp)
                out, man = self.run_once(vcf, f'out{seed}')
                v, rows = man['counts']['vcf'], man['rows']
                self.assertEqual(v['records'], n, 'a record went unread')
                self.assertEqual(rows['small_calls'], man['counts']['calls'])
                self.assertEqual(rows['rejected_calls'], v.get('rejected', 0) + v.get('unplaced_records', 0))
                self.assertGreaterEqual(v.get('small', 0), man['counts']['calls'])
                for name in ('small_calls', 'rejected_calls', 'sv_calls', 'repeat_calls', 'haplotype_calls', 'unresolved_loci', 'evidence_links', 'gene_context'):
                    t = pq.read_table(os.path.join(out, f'{name}.parquet'))
                    self.assertEqual(t.num_rows, rows[name], name)

    def test_the_same_odd_records_give_byte_identical_tables_twice(self):
        import compare_runs
        rng = random.Random(99)
        tmp = os.path.join(self.work, 'in_twice')
        os.makedirs(tmp)
        vcf, _ = sample_vcf_from(rng, 200, tmp)
        a, _ = self.run_once(vcf, 'twice_a')
        b, _ = self.run_once(vcf, 'twice_b')
        res = compare_runs.compare(a, b)
        self.assertTrue(res['identical'], res)


if __name__ == '__main__':
    unittest.main()
