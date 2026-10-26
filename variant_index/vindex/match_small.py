"""Matching of called SNVs and indels (and MNV components) against the index, by normalised CHROM:POS:REF:ALT only.

Two views over the same logic:
  * index-centric  `status_of(key)`: what does this sample say about this indexed allele? (EXACT_MATCH, FILTERED_CALL, NO_CALL,
    VERIFIED_REFERENCE, OTHER_ALLELE_AT_SITE, DELETED_BY_UPSTREAM, ABSENT_FROM_CALLSET, UNRESOLVED; see schema.MATCH_STATUS)
  * call-centric   `interpret_calls`: what does the index hold for each called allele? Calls with no index record are kept and say so.

Rules enforced here: identity is the key (never an rsID); a record that is absent from a variant-only VCF is ABSENT_FROM_CALLSET and
never reference; reference needs positive evidence from a callability provider; a multi-base substitution is matched from its
components only when the cis/trans relation is known, and otherwise returns UNRESOLVED with the bounds of what the dosages allow.
"""
from collections import defaultdict
from dataclasses import dataclass, field

from .callability import GvcfCallability, CompositeCallability
from .normalise import parse_key, classify, atomise, Norm
from . import summary as _summary


@dataclass
class Status:
    key: str
    status: str
    reasons: tuple = ()
    dosage: object = None
    zygosity: object = None
    basis: str = ''
    call: object = None
    detail: dict = field(default_factory=dict)


@dataclass
class CallInterpretation:
    call: object
    status: str                  # EXACT_MATCH | FILTERED_CALL | NO_INDEX_RECORD | NO_CALL
    evidence: list
    axes: dict
    genes: list
    regions: list
    flags: tuple
    mnv_context: list = field(default_factory=list)


def zygosity(dosage, ploidy):
    if ploidy == 1:
        return 'hemizygous'
    if dosage >= ploidy:
        return 'homozygous'
    return 'heterozygous'


class MatchContext:
    def __init__(self, callset, extra_providers=()):
        """The callset's own records (variants, no-calls, spanning deletions, structural calls, gVCF blocks) always count; extra providers
        (callable-region mask, reads) can only add reference evidence."""
        self.cs = callset
        self.callability = CompositeCallability([GvcfCallability(callset)] + list(extra_providers))
        self.by_key = defaultdict(list)
        self.by_span = defaultdict(list)       # (chrom) -> list of (start, end, call) for calls carrying an allele
        self.mnv_parts = defaultdict(list)     # component SNV key -> [sample MNV calls containing it]
        for c in callset.small:
            self.by_key[c.key].append(c)
            if c.dosage > 0:
                self.by_span[c.chrom].append((c.pos, c.pos + len(c.ref) - 1, c))
                if c.vtype == 'mnv':
                    for p in atomise(Norm(c.chrom, c.pos, c.ref, c.alt, 'mnv', 'ok', False)):
                        self.mnv_parts[p.key].append(c)

    def overlapping_calls(self, chrom, start, end):
        return [c for s, e, c in self.by_span.get(chrom, ()) if s <= end and e >= start]


def _best(calls):
    return sorted(calls, key=lambda c: (not c.qc_pass, c.missing, -(c.gq or 0), c.rec_index))[0]


def _evaluate_call(key, calls):
    c = _best(calls)
    if len({(x.dosage, x.missing) for x in calls if x.qc_pass}) > 1:
        return Status(key, 'UNRESOLVED', ('duplicate_records_disagree',), basis='call', call=c)
    z = zygosity(c.dosage, c.ploidy)
    if c.missing:
        return Status(key, 'NO_CALL', ('genotype_incomplete',), c.dosage, z, 'call', c)
    if not c.qc_pass:
        return Status(key, 'FILTERED_CALL', c.qc_flags, c.dosage, z, 'call', c)
    return Status(key, 'EXACT_MATCH', (), c.dosage, z, 'call', c, dict(gt=c.gt, phased=c.phased, phase_set=c.phase_set, gq=c.gq, dp=c.dp, multiallelic=c.multiallelic))


def _carry_vectors(parts):
    return [p.hap_alleles for p in parts]


def _mnv_from_components(ctx, key, c, pos, ref, alt):
    comps = atomise(Norm(c, pos, ref, alt, 'mnv', 'ok', False))
    found, missing = [], []
    for p in comps:
        calls = [x for x in ctx.by_key.get(p.key, []) if x.dosage > 0]
        (found if calls else missing).append((p, _best(calls) if calls else None))
    if missing:
        sts = [ctx.callability.status(p.chrom, p.pos, p.pos)[0] for p, _ in missing]
        if found and all(s == 'reference_confident' for s in sts):
            return Status(key, 'OTHER_ALLELE_AT_SITE', ('only_some_mnv_components_present',), basis='components',
                          detail=dict(present=[p.key for p, _ in found], missing=[p.key for p, _ in missing]))
        if found:
            return Status(key, 'UNRESOLVED', ('some_mnv_components_present_others_unknown',), basis='components',
                          detail=dict(present=[p.key for p, _ in found], missing=[p.key for p, _ in missing], callability=sts))
        return None            # no component is called: fall through to site-level reference logic
    parts = [x for _, x in found]
    if any((not x.qc_pass) or x.missing for x in parts):
        return Status(key, 'FILTERED_CALL', ('mnv_component_fails_quality',), basis='components', detail=dict(components=[p.key for p, _ in found]))
    P = {x.ploidy for x in parts}
    if len(P) != 1:
        return Status(key, 'UNRESOLVED', ('components_have_different_ploidy',), basis='components')
    ploidy = P.pop()
    d = [x.dosage for x in parts]
    lo, hi = max(0, sum(d) - (len(d) - 1) * ploidy), min(d)
    same_set = len({str(x.phase_set) for x in parts}) == 1 and parts[0].phase_set is not None
    if same_set and all(x.hap_alleles is not None for x in parts):
        both = sum(1 for h in range(ploidy) if all(x.hap_alleles[h] == 1 for x in parts))
        if both:
            return Status(key, 'EXACT_MATCH', (), both, zygosity(both, ploidy), 'components_in_phase_set', parts[0], dict(phase_set=str(parts[0].phase_set)))
        return Status(key, 'OTHER_ALLELE_AT_SITE', ('mnv_components_in_trans',), basis='components_in_phase_set', detail=dict(phase_set=str(parts[0].phase_set)))
    if lo == hi:
        if hi == 0:
            return Status(key, 'OTHER_ALLELE_AT_SITE', ('mnv_components_in_trans',), basis='dosage_bounds')
        return Status(key, 'EXACT_MATCH', (), hi, zygosity(hi, ploidy), 'components_dosage_implies_cis', parts[0], dict(dosages=d))
    return Status(key, 'UNRESOLVED', ('phase_required',), basis='components',
                  detail=dict(component_dosages=d, ploidy=ploidy, carries_all_components_min=lo, carries_all_components_max=hi,
                              missing_evidence='phase information (PS or read-backed phasing) linking the component calls'))


def status_of(ctx, key):
    c, pos, ref, alt = parse_key(key)
    calls = ctx.by_key.get(key)
    if calls and any(x.dosage > 0 for x in calls):
        return _evaluate_call(key, [x for x in calls if x.dosage > 0])
    if classify(ref, alt) == 'mnv':
        st = _mnv_from_components(ctx, key, c, pos, ref, alt)
        if st is not None:
            return st
    if key in ctx.mnv_parts:
        m = _best(ctx.mnv_parts[key])
        st = _evaluate_call(key, [m])
        st.basis = 'component_of_called_mnv'
        st.detail = dict(mnv=m.key)
        return st
    span = (pos, pos + len(ref) - 1)
    state, basis = ctx.callability.status(c, *span)
    if state == 'reference_confident':
        return Status(key, 'VERIFIED_REFERENCE', (), 0, 'homozygous_reference', basis)
    if state == 'reference_low_confidence':
        return Status(key, 'UNRESOLVED', ('reference_confidence_below_policy',), basis=basis, detail=dict(missing_evidence='reference block with GQ/depth above the policy, or reads'))
    if state == 'no_call':
        return Status(key, 'NO_CALL', ('site_without_genotype',), basis=basis)
    if state == 'deleted_by_upstream':
        return Status(key, 'DELETED_BY_UPSTREAM', ('spanning_deletion_allele',), basis=basis)
    if state == 'structural_overlap':
        return Status(key, 'UNRESOLVED', ('inside_called_structural_variant',), basis=basis,
                      detail=dict(missing_evidence='breakpoint-resolved structural call and phase to say which haplotype carries the indexed allele'))
    if state == 'variant_site':
        ov = ctx.overlapping_calls(c, *span)
        same_site = [x for x in ov if x.pos == pos or (x.pos <= span[1] and x.pos + len(x.ref) - 1 >= span[0] and x.vtype in ('snv', 'mnv') and len(ref) == 1)]
        indel_cover = [x for x in ov if x.vtype in ('del', 'complex') and x.pos < pos <= x.pos + len(x.ref) - 1]
        bad = [x for x in ov if (not x.qc_pass) or x.missing]
        if indel_cover:
            if all(x.dosage >= x.ploidy for x in indel_cover) and all(x.qc_pass for x in indel_cover):
                return Status(key, 'DELETED_BY_UPSTREAM', ('inside_homozygous_deletion_call',), basis='call', detail=dict(deletion=[x.key for x in indel_cover]))
            return Status(key, 'UNRESOLVED', ('inside_called_deletion_on_some_haplotypes',), basis='call', detail=dict(deletion=[x.key for x in indel_cover],
                          missing_evidence='phase to say whether the indexed allele would sit on the deleted haplotype'))
        if bad:
            return Status(key, 'UNRESOLVED', ('overlapping_call_fails_quality',), basis='call', detail=dict(calls=[x.key for x in bad], flags=[x.qc_flags for x in bad]))
        return Status(key, 'OTHER_ALLELE_AT_SITE', ('different_allele_called_at_site',), basis='call', detail=dict(called=[x.key for x in (same_site or ov)]))
    return Status(key, 'ABSENT_FROM_CALLSET', ('no_record_and_no_reference_evidence',), basis=basis,
                  detail=dict(missing_evidence='gVCF reference block, callable-region mask with all-variants assumption, or read-level evidence at this position'))


def match_index_variants(ctx, keys):
    """{key: Status} for a list of indexed variant keys."""
    return {k: status_of(ctx, k) for k in keys}


def is_difficult(regions):
    """True when a region of a difficult-region stratification overlaps (the XY context regions are about ploidy, not call quality)."""
    return any(r['kind'].startswith('giab_') and not r['kind'].startswith('giab_xy_') for r in regions)


def iter_interpretations(callset, index, ctx=None, chunk=200_000):
    """Call-centric interpretation of every small call, as a generator so a whole-genome callset never needs all its evidence rows in memory at once:
    evidence is fetched for `chunk` calls at a time. Calls with dosage 0 are not calls and do not appear."""
    ctx = ctx or MatchContext(callset)
    calls = [c for c in callset.small if c.dosage > 0]
    mnv_idx = index.mnv_components()
    for off in range(0, len(calls), chunk):
        part = calls[off:off + chunk]
        ev = index.evidence_for_variants({c.key for c in part})
        for c in part:
            rows = ev.get(c.key, [])
            mnv_ctx = []
            for mk in mnv_idx.get(c.key, []):
                mnv_ctx.append((mk, status_of(ctx, mk)))
            genes = index.genes_overlapping(c.chrom, c.pos, c.pos + len(c.ref) - 1)
            regions = index.regions_overlapping(c.chrom, c.pos, c.pos + len(c.ref) - 1)
            if c.missing:
                status = 'NO_CALL'
            elif not c.qc_pass:
                status = 'FILTERED_CALL'
            else:
                status = 'EXACT_MATCH' if rows else 'NO_INDEX_RECORD'
            flags = tuple(c.qc_flags) + (('multiallelic_split',) if c.multiallelic else ()) + (('difficult_region',) if is_difficult(regions) else ())
            yield CallInterpretation(c, status, rows, _summary.summarise(rows), genes, regions, flags, mnv_ctx)


def interpret_calls(callset, index, ctx=None):
    """List form of iter_interpretations (tests and small callsets)."""
    return list(iter_interpretations(callset, index, ctx))
