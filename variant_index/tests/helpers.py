"""Shared test helpers: synthetic reference, VCF writer, bcftools discovery."""
import atexit
import os
import shutil
import subprocess
import sys
import tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, '..'))
sys.path.insert(0, ROOT)

from vindex.refgenome import Reference  # noqa: E402

# Test modules build indexes and source trees in directories made with tempfile.mkdtemp() at import time and keep them for the whole run, so no single test can remove them. Every
# directory made while the tests run is removed when the interpreter exits (a full run used to leave about 20 MB in /tmp, which adds up on a small disk).
_MADE_DIRS = []
_mkdtemp = tempfile.mkdtemp


def _tracked_mkdtemp(*args, **kwargs):
    d = _mkdtemp(*args, **kwargs)
    _MADE_DIRS.append(d)
    return d


tempfile.mkdtemp = _tracked_mkdtemp
atexit.register(lambda: [shutil.rmtree(d, ignore_errors=True) for d in _MADE_DIRS])

FIXTURES = os.path.join(HERE, 'fixtures')
SYN_FASTA = os.path.join(FIXTURES, 'syn38.fa')
SYN_BUILDS = {3000: 'SYN38'}


def syn_reference():
    return Reference(SYN_FASTA, name='SYN38', builds=SYN_BUILDS)


def find_tool(name):
    """Path of an established command-line tool (bcftools, plink2, tabix): $VINDEX_TOOLS_BIN, then PATH, else None."""
    d = os.environ.get('VINDEX_TOOLS_BIN')
    if d and os.path.exists(os.path.join(d, name)):
        return os.path.join(d, name)
    return shutil.which(name)


def run_tool(args):
    p = subprocess.run(args, capture_output=True, text=True)
    return p.returncode, p.stdout, p.stderr


SYN_CONTIGS = [('1', 3000), ('2', 1500), ('X', 1500), ('MT', 400)]


def vcf_text(records, samples=('S1',), contigs=SYN_CONTIGS, extra_header=(), fmt='GT:GQ:DP:AD'):
    """records: list of (chrom, pos, ref, alt_csv, filter, qual, [sample strings...]) or with an id as 8th element."""
    h = ['##fileformat=VCFv4.2'] + [f'##contig=<ID={c},length={n}>' for c, n in contigs]
    h += ['##FILTER=<ID=LowQual,Description="low">',
          '##INFO=<ID=END,Number=1,Type=Integer,Description="end">', '##INFO=<ID=SVTYPE,Number=1,Type=String,Description="t">',
          '##INFO=<ID=SVLEN,Number=.,Type=Integer,Description="l">', '##INFO=<ID=CIPOS,Number=2,Type=Integer,Description="c">',
          '##INFO=<ID=CIEND,Number=2,Type=Integer,Description="c">', '##INFO=<ID=IMPRECISE,Number=0,Type=Flag,Description="i">',
          '##INFO=<ID=MATEID,Number=.,Type=String,Description="m">',
          '##FORMAT=<ID=GT,Number=1,Type=String,Description="gt">', '##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="gq">',
          '##FORMAT=<ID=DP,Number=1,Type=Integer,Description="dp">', '##FORMAT=<ID=MIN_DP,Number=1,Type=Integer,Description="mdp">',
          '##FORMAT=<ID=AD,Number=R,Type=Integer,Description="ad">', '##FORMAT=<ID=PS,Number=1,Type=Integer,Description="ps">',
          '##FORMAT=<ID=PGT,Number=1,Type=String,Description="pgt">', '##FORMAT=<ID=PID,Number=1,Type=String,Description="pid">',
          '##FORMAT=<ID=CN,Number=1,Type=Integer,Description="cn">', '##FORMAT=<ID=REPCN,Number=1,Type=String,Description="r">', '##FORMAT=<ID=REPCI,Number=1,Type=String,Description="r">',
          '##INFO=<ID=VARID,Number=1,Type=String,Description="v">', '##INFO=<ID=RU,Number=1,Type=String,Description="r">',
          '##ALT=<ID=DEL,Description="d">', '##ALT=<ID=DUP,Description="d">', '##ALT=<ID=INV,Description="d">', '##ALT=<ID=INS,Description="d">',
          '##ALT=<ID=CNV,Description="d">', '##ALT=<ID=NON_REF,Description="n">']
    h += list(extra_header)
    h.append('#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t' + '\t'.join(samples))
    lines = []
    for r in records:
        chrom, pos, ref, alt, filt, qual, rest = r[:7]
        info = r[7] if len(r) > 7 else '.'
        vid = r[8] if len(r) > 8 else '.'
        f = r[9] if len(r) > 9 else fmt
        lines.append('\t'.join([chrom, str(pos), vid, ref, alt, str(qual), filt, info, f] + list(rest)))
    return '\n'.join(h + lines) + '\n'
