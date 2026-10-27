"""Sample VCF/gVCF reading for the matcher.

Every ALT allele of every record becomes one `SmallCall` keyed by its normalised CHROM:POS:REF:ALT, with the sample's genotype
re-expressed relative to that allele: `dosage` (copies of this ALT among the called alleles), `other_alt` (copies of another ALT of
the same site), `ref_copies`, `missing`. The original GT string, phase flag, phase set, QUAL, FILTER, GQ, DP and AD are kept, so
splitting a multiallelic record loses nothing. Nothing is dropped silently: rejected alleles (REF mismatch, invalid, unmapped contig),
no-calls, reference blocks, symbolic/structural records and spanning-deletion markers go to their own lists with reasons.
Personal VCFs are read from local paths only; nothing here writes or transmits them.
"""
import re
from dataclasses import dataclass, field

import pysam

from .normalise import normalise

SV_SYMBOLIC = re.compile(r'^<([A-Z]+)(?::([A-Za-z0-9:_]+))?>$')


@dataclass
class QCPolicy:
    require_pass: bool = True      # FILTER must be PASS or '.'
    min_gq: int = 20
    min_dp: int = 10
    min_qual: float = 0.0
    min_ref_gq: int = 20           # a reference call/block counts as confident only with GQ and DP at least this
    min_ref_dp: int = 10
    sv_min_len: int = 50           # sequence-resolved indels at least this long are also offered to the SV matcher


@dataclass
class SmallCall:
    key: str
    chrom: str
    pos: int
    ref: str
    alt: str
    vtype: str
    gt: str
    alleles: tuple
    ploidy: int
    phased: bool
    phase_set: object
    dosage: int
    other_alt: int
    ref_copies: int
    missing: int
    qual: object
    filter: str
    gq: object
    dp: object
    ad_ref: object
    ad_alt: object
    qc_pass: bool
    qc_flags: tuple
    rec_index: int
    orig: tuple
    multiallelic: bool
    detail: str = ''
    hap_alleles: object = None     # allele indices per haplotype in phased order when known (GT with '|', or GATK PGT)


@dataclass
class Rejected:
    reason: str
    contig: str
    pos: int
    ref: str
    alt: str
    rec_index: int
    detail: str = ''
    canon: str = ''        # the contig as the reference names it, where the record was placed on one


@dataclass
class NoCall:
    chrom: str
    pos: int
    end: int
    rec_index: int
    reason: str


@dataclass
class RefBlock:
    chrom: str
    start: int
    end: int
    gq: object
    dp: object
    confident: bool
    basis: str           # gvcf_block | genotype_ref
    rec_index: int


@dataclass
class SVCall:
    id: str
    chrom: str
    start: int
    end: int
    svtype: str
    svlen: object
    cipos: tuple
    ciend: tuple
    precise: bool
    cn: object
    gt: str
    dosage: object
    ploidy: int
    qual: object
    filter: str
    qc_pass: bool
    rec_index: int
    alt: str
    mate: object = None
    source: str = 'symbolic'
    flags: tuple = ()


@dataclass
class RepeatCall:
    locus_id: str
    chrom: str
    start: int
    end: int
    motif: str
    counts: tuple          # repeat units per allele (None when unknown)
    ci: tuple              # (lo, hi) per allele
    ploidy: int
    source: str            # expansionhunter | trgt
    qc_pass: bool
    rec_index: int


@dataclass
class SpanDel:
    chrom: str
    pos: int
    end: int
    copies: int
    rec_index: int


@dataclass
class CallSet:
    sample: str
    path: str
    assembly: dict
    small: list = field(default_factory=list)
    rejected: list = field(default_factory=list)
    nocalls: list = field(default_factory=list)
    refblocks: list = field(default_factory=list)
    svs: list = field(default_factory=list)
    repeats: list = field(default_factory=list)
    spandels: list = field(default_factory=list)
    unplaced: list = field(default_factory=list)
    unused: list = field(default_factory=list)      # ALT alleles listed in a record but absent from the sample's genotype
    stats: dict = field(default_factory=dict)


NONREF = ('<NON_REF>', '<*>')


def _fmt(s, key):
    """One FORMAT value of a pysam sample (None when absent); length-1 tuples are unwrapped."""
    try:
        v = s[key]
    except (KeyError, ValueError):
        return None
    if isinstance(v, tuple) and len(v) == 1:
        v = v[0]
    if isinstance(v, float) and v != v:
        return None
    return v


def _raw_info(rec):
    """INFO key -> raw string, read from the record text (pysam hides END behind its own interval arithmetic)."""
    f = str(rec).rstrip('\n').split('\t')
    out = {}
    for kv in (f[7].split(';') if len(f) > 7 and f[7] != '.' else []):
        k, _, v = kv.partition('=')
        out[k] = v
    return out


def _is_symbolic(a):
    return a.startswith('<') or '[' in a or ']' in a


def _gt_string(gt, phased):
    return ('|' if phased else '/').join('.' if a is None else str(a) for a in gt)


def _filter_name(rec):
    return ';'.join(rec.filter.keys()) if len(rec.filter) else 'PASS'


def _qc_flags(qc, filt, qual, gq, dp, gt):
    flags = []
    if qc.require_pass and filt not in ('PASS', '.'):
        flags.append('FILTER:' + filt)
    if gq is not None and gq < qc.min_gq:
        flags.append(f'GQ<{qc.min_gq}')
    if dp is not None and dp < qc.min_dp:
        flags.append(f'DP<{qc.min_dp}')
    if qual is not None and qual < qc.min_qual:
        flags.append(f'QUAL<{qc.min_qual}')
    if gt is None:
        flags.append('no_GT')
    return flags


def _ref_confident(qc, filt, gq, mdp):
    return (gq is None or gq >= qc.min_ref_gq) and (mdp is None or mdp >= qc.min_ref_dp) and not (qc.require_pass and filt not in ('PASS', '.'))


OTHER_CONTIGS = '*other*'      # read_vcf(contigs=[OTHER_CONTIGS]): every record whose contig is not on the primary assembly, including contigs the reference does not know


def _wanted(reference, contigs):
    wanted, other = set(contigs), OTHER_CONTIGS in contigs

    def test(name):
        canon = reference.contig(name)
        return (canon is not None and canon in wanted) or (other and not (canon is not None and reference.is_primary(canon)))
    return test


def _records(vf, reference, contigs):
    """(running index, record) for the whole file, or only for the records on `contigs` (canonical names, or OTHER_CONTIGS for everything off the primary
    assembly). An indexed file is fetched per contig; an unindexed one is read through and filtered, so the result is the same either way."""
    if contigs is None:
        yield from enumerate(vf)
        return
    test = _wanted(reference, contigs)
    names = [c for c in vf.header.contigs if test(c)]
    idx = 0
    try:
        for n in names:
            for rec in vf.fetch(n):
                yield idx, rec
                idx += 1
        return
    except (ValueError, OSError):           # no index: fall through to a filtered pass
        if idx:
            raise
    vf.reset()
    for rec in vf:
        if test(rec.chrom):
            yield idx, rec
            idx += 1


def read_vcf(path, reference, sample=None, qc=None, contigs=None):
    """Read one sample of a VCF/gVCF into a CallSet (see module docstring). Raises only for unusable input (no/ambiguous sample).
    contigs: canonical contig names to read (default all); a whole-genome file is matched one chromosome at a time to bound memory."""
    qc = qc or QCPolicy()
    vf = pysam.VariantFile(path)
    samples = list(vf.header.samples)
    if sample is None:
        if len(samples) != 1:
            raise ValueError(f'{path}: {len(samples)} samples {samples}; name the sample to read')
        sample = samples[0]
    elif sample not in samples:
        raise ValueError(f'{path}: sample {sample!r} not in {samples}')
    header_contigs = {c.name: c.length for c in vf.header.contigs.values()}
    has_trid = 'TRID' in vf.header.info
    cs = CallSet(sample=sample, path=path, assembly=reference.assembly_check(header_contigs))
    st = dict(records=0, alleles=0, small=0, split_multiallelic_records=0, rejected=0, nocall_records=0, ref_blocks=0, sv_records=0,
              spanning_markers=0, unplaced_records=0, mnv_calls=0, repeat_records=0)
    for idx, rec in _records(vf, reference, contigs):
        st['records'] += 1
        canon = reference.contig(rec.chrom)
        alts = tuple(rec.alts or ())
        if canon is None:
            st['unplaced_records'] += 1
            cs.unplaced.append(Rejected('contig_unmapped', rec.chrom, rec.pos, rec.ref, ','.join(alts), idx))
            continue
        s = rec.samples[sample]
        gt = tuple(s['GT']) if 'GT' in rec.format and s['GT'] is not None else None
        phased = bool(s.phased) if gt else False
        gq, dp = _fmt(s, 'GQ'), _fmt(s, 'DP')
        filt, qual = _filter_name(rec), rec.qual
        end = rec.stop
        called = bool(gt) and all(a is not None for a in gt)
        rep = _repeat_call(rec, s, canon, idx, gt, filt, qc, has_trid)
        if rep is not None:
            cs.repeats.append(rep)
            st['repeat_records'] += 1
            continue
        real = [(i, a) for i, a in enumerate(alts, start=1) if a not in NONREF]
        # gVCF reference blocks: no ALT, or only <NON_REF>/<*>
        if not real:
            if called and all(a == 0 for a in gt):
                mdp = _fmt(s, 'MIN_DP') if 'MIN_DP' in rec.format else dp
                cs.refblocks.append(RefBlock(canon, rec.pos, end, gq, mdp, bool(_ref_confident(qc, filt, gq, mdp)), 'gvcf_block', idx))
                st['ref_blocks'] += 1
            else:
                cs.nocalls.append(NoCall(canon, rec.pos, end, idx, 'reference block without a reference genotype'))
                st['nocall_records'] += 1
            continue
        flags = _qc_flags(qc, filt, qual, gq, dp, gt)
        if any(_is_symbolic(a) for _, a in real):
            st['sv_records'] += 1
            cs.svs.extend(_sv_calls(rec, idx, canon, [(i, a) for i, a in real if _is_symbolic(a)], gt, filt, qual, called, qc, sample))
        real = [(i, a) for i, a in real if not _is_symbolic(a)]
        ploidy = len(gt) if gt else 0
        for i, a in real:
            if a == '*' and gt and any(x == i for x in gt):
                cs.spandels.append(SpanDel(canon, rec.pos, rec.pos + len(rec.ref) - 1, sum(1 for x in gt if x == i), idx))
                st['spanning_markers'] += 1
        if not called:
            cs.nocalls.append(NoCall(canon, rec.pos, rec.pos + len(rec.ref) - 1, idx, 'GT absent' if not gt else 'GT has a missing allele'))
            st['nocall_records'] += 1
            if not gt or all(a is None for a in gt):
                continue
        if called and all(a == 0 for a in gt):   # sample is homozygous reference at a site that is variant in others
            cs.refblocks.append(RefBlock(canon, rec.pos, rec.pos + len(rec.ref) - 1, gq, dp, bool(_ref_confident(qc, filt, gq, dp)), 'genotype_ref', idx))
            st['ref_blocks'] += 1
            continue
        alleles = [(i, a) for i, a in real if a != '*']
        multi = len(alleles) > 1
        st['split_multiallelic_records'] += multi
        ad = tuple(s['AD']) if 'AD' in rec.format and s['AD'] is not None else None
        ps = _fmt(s, 'PS') if 'PS' in rec.format else None
        pid = _fmt(s, 'PID') if 'PID' in rec.format else None
        pgt = _fmt(s, 'PGT') if 'PGT' in rec.format else None
        # PS is meaningful only for phased genotypes; an unphased GT with GATK PGT/PID carries read-backed phase relative to its PID block
        phase_set = (ps if ps is not None else (pid if pid is not None else 'implicit')) if phased else None
        hap = tuple(gt) if phased else None
        f_hap = []
        if hap is None and isinstance(pgt, str) and re.match(r'^\d\|\d$', pgt) and len(gt) == 2 and pid is not None:
            hap = tuple(int(x) for x in pgt.split('|'))
            phase_set = pid
            f_hap = ['phase_from_PGT']
        for ai, alt in alleles:
            st['alleles'] += 1
            if called and not any(x == ai for x in gt):
                cs.unused.append(Rejected('allele_not_in_genotype', rec.chrom, rec.pos, rec.ref, alt, idx))
                st['alleles_not_in_genotype'] = st.get('alleles_not_in_genotype', 0) + 1
                continue
            n = normalise(reference, canon, rec.pos, rec.ref, alt)
            if not n.ok:
                st['rejected'] += 1
                cs.rejected.append(Rejected(n.status, rec.chrom, rec.pos, rec.ref, alt, idx, n.detail, canon))
                continue
            f2 = list(flags) + ([] if called else ['partial_no_call'])
            call = SmallCall(n.key, n.chrom, n.pos, n.ref, n.alt, n.vtype, _gt_string(gt, phased), tuple(gt), ploidy, phased, phase_set,
                             sum(1 for x in gt if x == ai), sum(1 for x in gt if x not in (None, 0, ai)), sum(1 for x in gt if x == 0),
                             sum(1 for x in gt if x is None), qual, filt, gq, dp, ad[0] if ad else None, ad[ai] if ad and len(ad) > ai else None,
                             not f2, tuple(f2), idx, (rec.chrom, rec.pos, rec.ref, alt), multi, n.detail,
                             tuple(None if x is None else int(x == ai) for x in hap) if hap else None)
            call.qc_flags = tuple(f2 + f_hap)
            call.qc_pass = not f2
            cs.small.append(call)
            st['small'] += 1
            st['mnv_calls'] += n.vtype == 'mnv'
            if abs(len(n.ref) - len(n.alt)) >= qc.sv_min_len:
                cs.svs.append(_seq_sv(call, idx, qc))
    cs.stats = st
    vf.close()
    return cs


def _repeat_call(rec, s, canon, idx, gt, filt, qc, has_trid):
    """ExpansionHunter (FORMAT REPCN) and TRGT (INFO TRID with FORMAT AL/MC) repeat genotypes; None for any other record."""
    raw = None
    if 'REPCN' in rec.format and s['REPCN'] is not None:
        raw = _raw_info(rec)
        counts = tuple(None if x in ('.', '') else int(x) for x in str(s['REPCN']).split('/'))
        ci = []
        for x in (str(s['REPCI']).split('/') if 'REPCI' in rec.format and s['REPCI'] is not None else []):
            lo, _, hi = x.partition('-')
            ci.append((int(lo), int(hi)) if lo.isdigit() and hi.isdigit() else (None, None))
        end = int(raw['END']) if raw.get('END', '').isdigit() else rec.stop
        return RepeatCall(raw.get('VARID') or raw.get('REPID') or f'{canon}:{rec.pos}', canon, rec.pos, end, raw.get('RU', ''), counts, tuple(ci), len(counts),
                          'expansionhunter', filt in ('PASS', '.') or not qc.require_pass, idx)
    raw = _raw_info(rec) if has_trid and rec.info.get('TRID') is not None else None
    if raw is not None and 'AL' in rec.format and s['AL'] is not None:
        motifs = (raw.get('MOTIFS') or '').split(',')
        unit = len(motifs[0]) if motifs and motifs[0] else 0
        als = [x for x in (s['AL'] if isinstance(s['AL'], tuple) else (s['AL'],))]
        counts = tuple(round(a / unit) if unit and a is not None else None for a in als)
        end = int(raw['END']) if raw.get('END', '').isdigit() else rec.stop
        return RepeatCall(raw.get('TRID', ''), canon, rec.pos, end, motifs[0] if motifs else '', counts, tuple(), len(counts), 'trgt',
                          filt in ('PASS', '.') or not qc.require_pass, idx)
    return None


def _seq_sv(c, idx, qc):
    typ = 'DEL' if len(c.ref) > len(c.alt) else 'INS'
    L = abs(len(c.ref) - len(c.alt))
    end = c.pos + len(c.ref) - 1
    return SVCall(c.key, c.chrom, c.pos, end, typ, -L if typ == 'DEL' else L, (0, 0), (0, 0), True, None, c.gt, c.dosage, c.ploidy,
                  c.qual, c.filter, c.qc_pass, idx, c.alt[:30], None, 'sequence_resolved')


def _sv_calls(rec, idx, canon, symbolic, gt, filt, qual, called, qc, sample):
    out = []
    info = rec.info
    ploidy = len(gt) if gt else 0
    for ai, alt in symbolic:
        m = SV_SYMBOLIC.match(alt)
        svtype = info.get('SVTYPE') or (m.group(1) if m else 'BND')
        svtype = svtype[0] if isinstance(svtype, tuple) else svtype
        svlen = info.get('SVLEN')
        svlen = (svlen[0] if svlen else None) if isinstance(svlen, tuple) else svlen
        raw = _raw_info(rec)
        flags = []
        end_info = int(raw['END']) if raw.get('END', '').isdigit() else None
        if svlen is not None and svtype in ('DEL', 'DUP', 'INV', 'CNV') and end_info is not None and abs(abs(svlen) - (end_info - rec.pos)) > 1:
            flags.append('END_SVLEN_inconsistent')
        if end_info is not None:
            end = end_info
        elif svlen not in (None, 0) and svtype in ('DEL', 'DUP', 'INV', 'CNV'):
            end = rec.pos + abs(svlen)
        else:
            end = rec.pos
        cipos = tuple(info['CIPOS']) if 'CIPOS' in info else (0, 0)
        ciend = tuple(info['CIEND']) if 'CIEND' in info else (0, 0)
        cn = None
        if 'CN' in rec.format:
            v = rec.samples[sample]['CN']
            cn = v[0] if isinstance(v, tuple) else v
        dosage = sum(1 for a in gt if a == ai) if called else None
        mate = None
        if svtype == 'BND':
            mm = re.search(r'[\[\]]([^\[\]:]+):(\d+)[\[\]]', alt)
            mate = (mm.group(1), int(mm.group(2))) if mm else None
        ok_filter = filt in ('PASS', '.') or not qc.require_pass
        out.append(SVCall(rec.id or f'{canon}:{rec.pos}:{svtype}:{end}', canon, rec.pos, int(end), svtype, svlen, cipos, ciend,
                          'IMPRECISE' not in info, cn, _gt_string(gt, False) if gt else '.', dosage, ploidy, qual, filt, bool(ok_filter and called),
                          idx, alt, mate, 'symbolic', tuple(flags)))
    return out
