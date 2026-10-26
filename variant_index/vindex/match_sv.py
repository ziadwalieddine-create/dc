"""Structural variant and copy-number matching.

Interval overlap alone is never equivalence. Two events are equivalent only when ALL of the following hold:
  1. compatible event type (DEL~copy-number loss, DUP~copy-number gain, INV~INV, INS~INS; anything else is not equivalent);
  2. breakpoints agree: both ends within a tolerance that grows with the stated confidence intervals of either event but never
     beyond `max_tol_frac` of the event length, or, for events without breakpoint resolution, reciprocal overlap >= `min_ro`;
  3. dosage is consistent where both sides state it (copy number or zygosity); if the index states a copy number the sample
     does not, the pair is flagged 'dosage_unresolved' and is not equivalent.
Everything else that overlaps is returned with its metrics under a non-equivalent relation (partial_overlap, overlap_only,
type_mismatch, dosage_mismatch) so a reviewer can see why it was refused. Gene and region evidence reaches a structural call only
through containment/overlap context rows, which are labelled as gene- or region-level and never counted as variant-level matches.
"""
from dataclasses import dataclass

LOSS_TYPES = {'DEL', 'CNV_LOSS'}
GAIN_TYPES = {'DUP', 'CNV_GAIN'}


@dataclass
class SvPolicy:
    base_tol: int = 50           # bp: allowed breakpoint difference for precise calls
    max_tol_frac: float = 0.10   # confidence intervals may widen the tolerance up to this fraction of the event length
    min_ro: float = 0.90         # reciprocal overlap required when breakpoints are not both precise
    ins_tol: int = 50
    ins_len_ratio: float = 0.80


@dataclass
class SvMatch:
    event_id: str
    relation: str       # equivalent_breakpoints | equivalent_overlap | partial_overlap | overlap_only | type_mismatch | dosage_mismatch | dosage_unresolved | unsupported
    metrics: dict
    reasons: tuple = ()
    equivalent: bool = False


@dataclass
class SvInterpretation:
    sv: object
    matches: list
    context: list
    unresolved: tuple = ()


def norm_type(t):
    t = (t or '').upper()
    return {'DELETION': 'DEL', 'COPY NUMBER LOSS': 'CNV_LOSS', 'LOSS': 'CNV_LOSS', 'DUPLICATION': 'DUP', 'TANDEM DUPLICATION': 'DUP',
            'COPY NUMBER GAIN': 'CNV_GAIN', 'GAIN': 'CNV_GAIN', 'INVERSION': 'INV', 'INSERTION': 'INS', 'DUP:TANDEM': 'DUP'}.get(t, t)


def sample_type(sv):
    """Event class of a sample call. A generic CNV is resolved from its copy number; without one it stays 'CNV' (unresolved)."""
    t = norm_type(sv.svtype)
    if t == 'CNV':
        if sv.cn is not None:
            return 'CNV_LOSS' if sv.cn < 2 else ('CNV_GAIN' if sv.cn > 2 else 'CNV')
        return 'CNV'
    return t


def types_compatible(a, b):
    if a in LOSS_TYPES and b in LOSS_TYPES:
        return True
    if a in GAIN_TYPES and b in GAIN_TYPES:
        return True
    return a == b and a in ('INV', 'INS')


def _ci(c):
    return max(abs(c[0]), abs(c[1])) if c else 0


def compare(sv, ev, policy=None):
    """Compare one sample call (SVCall) with one index event row (dict of entity_event). Returns SvMatch."""
    p = policy or SvPolicy()
    st, et = sample_type(sv), norm_type(ev['event_type'])
    L1, L2 = max(1, sv.end - sv.start), max(1, ev['end'] - ev['start'])
    ov = max(0, min(sv.end, ev['end']) - max(sv.start, ev['start']))
    ro = ov / max(L1, L2)
    ds, de = abs(sv.start - ev['start']), abs(sv.end - ev['end'])
    metrics = dict(sample_type=st, index_type=et, sample_len=L1, index_len=L2, overlap_bp=ov, reciprocal_overlap=round(ro, 4), start_diff=ds, end_diff=de,
                   size_ratio=round(min(L1, L2) / max(L1, L2), 4))
    if st == 'BND' or et == 'BND':
        return SvMatch(ev['event_id'], 'unsupported', metrics, ('breakend_records_are_not_matched',))
    if st == 'CNV' and et in LOSS_TYPES | GAIN_TYPES:
        return SvMatch(ev['event_id'], 'dosage_unresolved', metrics, ('sample_cnv_without_copy_number',))
    if not types_compatible(st, et):
        return SvMatch(ev['event_id'], 'type_mismatch', metrics, (f'{st}_vs_{et}',))
    if st == 'INS':
        dpos = abs(sv.start - ev['start'])
        ls, le = abs(sv.svlen or 0), abs(ev.get('svlen') or 0)
        ratio = (min(ls, le) / max(ls, le)) if ls and le else None
        metrics.update(pos_diff=dpos, ins_len_ratio=ratio)
        if dpos <= p.ins_tol and ratio is not None and ratio >= p.ins_len_ratio:
            return SvMatch(ev['event_id'], 'equivalent_breakpoints', metrics, ('inserted_sequence_not_compared',), True)
        return SvMatch(ev['event_id'], 'overlap_only' if dpos <= 10 * p.ins_tol else 'type_mismatch', metrics, ('insertion_position_or_length_differs',))
    tol_s = min(max(p.base_tol, _ci(sv.cipos), ev.get('start_ci_hi') or 0, abs(ev.get('start_ci_lo') or 0)), max(p.base_tol, int(p.max_tol_frac * min(L1, L2))))
    tol_e = min(max(p.base_tol, _ci(sv.ciend), ev.get('end_ci_hi') or 0, abs(ev.get('end_ci_lo') or 0)), max(p.base_tol, int(p.max_tol_frac * min(L1, L2))))
    metrics.update(tol_start=tol_s, tol_end=tol_e)
    breakpoints_ok = ds <= tol_s and de <= tol_e
    precise = bool(sv.precise) and bool(ev.get('precise', True))
    if breakpoints_ok:
        rel = 'equivalent_breakpoints'
    elif not precise and ro >= p.min_ro:
        rel = 'equivalent_overlap'
    elif ro >= 0.5:
        return SvMatch(ev['event_id'], 'partial_overlap', metrics, ('breakpoints_differ_beyond_tolerance' if precise else 'reciprocal_overlap_below_threshold',))
    else:
        return SvMatch(ev['event_id'], 'overlap_only', metrics, ('reciprocal_overlap_below_0.5',))
    # dosage
    ecn = ev.get('copy_number')
    if ecn is not None:
        if sv.cn is not None and sv.cn != ecn:
            return SvMatch(ev['event_id'], 'dosage_mismatch', metrics, (f'copy_number_{sv.cn}_vs_{ecn}',))
        if sv.cn is None:
            return SvMatch(ev['event_id'], 'dosage_unresolved', metrics, ('index_states_copy_number_sample_does_not',))
    return SvMatch(ev['event_id'], rel, metrics, (), True)


def match_svs(callset, index, policy=None, flank=10000):
    """Interpret every structural call: equivalent events (with their evidence), refused near-misses (with metrics) and gene/region context."""
    policy = policy or SvPolicy()
    out = []
    for sv in callset.svs:
        if sv.dosage == 0:                      # a homozygous-reference genotype at a structural site (joint-called files): the record is kept, but the sample does not carry the allele
            out.append(SvInterpretation(sv, [], [], ('genotype_is_homozygous_reference: the sample does not carry this structural allele',)))
            continue
        if sv.svtype == 'BND':
            out.append(SvInterpretation(sv, [], [], ('breakend_records_are_not_matched; need paired breakpoint evidence',)))
            continue
        cands = index.events_overlapping(sv.chrom, sv.start - flank, sv.end + flank) if sv.svtype != 'INS' else index.events_overlapping(sv.chrom, sv.start - 500, sv.start + 500)
        matches = [compare(sv, e, policy) for e in cands]
        ctx = []
        for g in index.genes_overlapping(sv.chrom, sv.start, sv.end):
            full = sv.start <= g['start'] and sv.end >= g['end']
            ctx.append(dict(level='gene', id=g['gene_id'], symbol=g['symbol'], relation='encompasses_gene' if full else 'partially_overlaps_gene',
                            applies_as='gene_level_context_only'))
        for r in index.regions_overlapping(sv.chrom, sv.start, sv.end):
            full = sv.start <= r['start'] and sv.end >= r['end']
            ctx.append(dict(level='region', id=r['region_id'], kind=r['kind'], relation='encompasses_region' if full else 'partially_overlaps_region',
                            applies_as='region_level_context_only'))
        unresolved = []
        if sv.svtype in ('DEL', 'DUP', 'CNV') and sv.cn is None and sv.dosage is None:
            unresolved.append('dosage unknown: no copy number and no usable genotype')
        if sv.flags if hasattr(sv, 'flags') else False:
            unresolved.extend(sv.flags)
        out.append(SvInterpretation(sv, matches, ctx, tuple(unresolved)))
    return out
