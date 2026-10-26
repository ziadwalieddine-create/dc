"""Polygenic score models, kept separate from the variant index.

A model is applied only to dosages that the input can support. For every weight the result states why it was or was not used;
a weight whose genotype is unknown is MISSING, never zero: a site absent from a variant-only VCF carries no reference evidence.
Checks reported with every score (they gate interpretation, they do not change the arithmetic):
  * allele dosage   effect-allele copies from the normalised call or from reference-confidence evidence; third alleles, palindromic
                    (A/T, C/G) SNVs without a resolvable strand, REF mismatches and weight types this module cannot score are excluded
  * missingness     counts and the share of total absolute weight that was scored
  * ancestry        development and evaluation ancestry from the PGS Catalog metadata against the sample's genetically inferred group
                    (unassessed unless supplied; inferring ancestry needs a reference panel that is not part of this repository)
  * calibration     percentile or z-score only when a reference distribution for the same score and ancestry is supplied; otherwise raw only
Scores are research outputs. They are not risk predictions.
"""
import gzip
import re
from dataclasses import dataclass, field

from .match_small import status_of
from .normalise import normalise, parse_key

COMP = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
PALINDROMIC = {frozenset('AT'), frozenset('CG')}
AUTOSOMES = {str(i) for i in range(1, 23)}


@dataclass
class Weight:
    chrom: str
    pos: int
    effect: str
    other: object
    weight: float
    row: int
    rsid: str = ''
    kind: str = 'additive'      # additive | dominant | recessive | interaction | haplotype


@dataclass
class ScoringFile:
    meta: dict
    weights: list
    problems: list = field(default_factory=list)


@dataclass
class PgsResult:
    pgs_id: str
    score_sum: float
    n_weights: int
    n_scored: int
    n_missing: int
    abs_weight_scored_fraction: float
    reasons: dict
    ancestry: dict
    calibration: dict
    flags: tuple
    per_weight: list = field(default_factory=list)


def _open(p):
    return gzip.open(p, 'rt', encoding='utf-8') if str(p).endswith('.gz') else open(p, encoding='utf-8')


def read_scoring_file(path):
    meta, weights, problems = {}, [], []
    cols = None
    with _open(path) as f:
        for i, line in enumerate(f):
            line = line.rstrip('\r\n')
            if line.startswith('#'):
                if '=' in line:
                    k, _, v = line[1:].partition('=')
                    meta[k.strip()] = v.strip()
                continue
            if not line:
                continue
            if cols is None:
                cols = line.split('\t')
                continue
            r = dict(zip(cols, line.split('\t')))
            use_hm = r.get('hm_chr') not in (None, '', 'NA') and r.get('hm_pos') not in (None, '', 'NA')
            chrom = r.get('hm_chr') if use_hm else r.get('chr_name')
            pos = r.get('hm_pos') if use_hm else r.get('chr_position')
            ea = (r.get('effect_allele') or '').upper()
            oa = (r.get('other_allele') or r.get('hm_inferOtherAllele') or '').upper() or None
            try:
                w = float(r.get('effect_weight'))
            except (TypeError, ValueError):
                problems.append((i, 'effect_weight not numeric'))
                continue
            kind = 'additive'
            if (r.get('is_interaction') or '').upper() in ('TRUE', '1'):
                kind = 'interaction'
            elif (r.get('is_dominant') or '').upper() in ('TRUE', '1'):
                kind = 'dominant'
            elif (r.get('is_recessive') or '').upper() in ('TRUE', '1'):
                kind = 'recessive'
            elif 'haplotype' in cols and r.get('haplotype'):
                kind = 'haplotype'
            if not chrom or not pos or not str(pos).isdigit():
                weights.append(Weight(chrom or '', 0, ea, oa, w, i, r.get('rsID', ''), kind))
                continue
            weights.append(Weight(chrom, int(pos), ea, oa, w, i, r.get('hm_rsID') or r.get('rsID', ''), kind))
    return ScoringFile(meta, weights, problems)


def _ancestry_assess(meta_dev, meta_eval, sample_group):
    if sample_group is None:
        return dict(status='unassessed', note='no inferred ancestry supplied; inference needs a reference panel (for example 1000 Genomes) projected onto the sample',
                    development=meta_dev, evaluation=meta_eval)
    def share(txt):
        out = {}
        for part in (txt or '').split('|'):
            k, _, v = part.partition(':')
            try:
                out[k.strip()] = float(v)
            except ValueError:
                pass
        return out
    dev, ev = share(meta_dev), share(meta_eval)
    dev_share = dev.get(sample_group, 0.0)
    ev_share = ev.get(sample_group, 0.0)
    status = 'matched' if dev_share >= 20 or ev_share >= 20 else 'not_represented'
    return dict(status=status, sample_group=sample_group, development_share=dev_share, evaluation_share=ev_share, development=meta_dev, evaluation=meta_eval)


def score(sf, ctx, reference, sample_ancestry=None, ancestry_meta=None, reference_distribution=None, autosomes_only=True):
    reasons, per, tot, scored_abs, total_abs, s, scored = {}, [], 0, 0.0, 0.0, 0.0, 0
    for w in sf.weights:
        tot += 1
        total_abs += abs(w.weight)
        d, why = _dosage(w, ctx, reference, autosomes_only)
        per.append((w.row, w.chrom, w.pos, w.effect, w.other, w.weight, d, why))
        if d is None:
            reasons[why] = reasons.get(why, 0) + 1
            continue
        scored += 1
        scored_abs += abs(w.weight)
        s += d * w.weight
    flags = []
    cov = scored_abs / total_abs if total_abs else 0.0
    if cov < 0.9:
        flags.append('weight_coverage_below_0.9')
    am = ancestry_meta or {}
    anc = _ancestry_assess(am.get('development'), am.get('evaluation'), sample_ancestry)
    if anc['status'] in ('unassessed', 'not_represented'):
        flags.append('ancestry_' + anc['status'])
    cal = dict(status='uncalibrated', note='raw weighted sum only; a percentile or z-score needs a reference distribution for this score in the same ancestry')
    if reference_distribution:
        mu, sd = reference_distribution['mean'], reference_distribution['sd']
        cal = dict(status='calibrated_to_supplied_distribution', z=(s - mu) / sd if sd else None, reference=reference_distribution)
    else:
        flags.append('uncalibrated')
    return PgsResult(sf.meta.get('pgs_id', ''), s, tot, scored, tot - scored, cov, dict(sorted(reasons.items())), anc, cal, tuple(flags), per)


def _dosage(w, ctx, ref, autosomes_only):
    if w.kind != 'additive':
        return None, 'unsupported_weight_type:' + w.kind
    if not w.pos:
        return None, 'no_position'
    canon = ref.contig(w.chrom) or ref.contig('chr' + w.chrom)
    if canon is None:
        return None, 'contig_unmapped'
    if autosomes_only and canon not in AUTOSOMES:
        return None, 'non_autosomal_excluded'
    ea, oa = w.effect, w.other
    if not re.match(r'^[ACGT]+$', ea or '') or (oa and not re.match(r'^[ACGT]+$', oa)):
        return None, 'allele_not_ACGT'
    if len(ea) == 1 and (oa is None or len(oa) == 1):
        return _snv_dosage(w, canon, ea, oa, ctx, ref)
    return _indel_dosage(w, canon, ea, oa, ctx, ref)


def _snv_dosage(w, canon, ea, oa, ctx, ref):
    rb = ref.fetch(canon, w.pos, w.pos)
    if oa and frozenset((ea, oa)) in PALINDROMIC:
        return None, 'palindromic_strand_ambiguous'
    flip = False
    if oa:
        if rb in (ea, oa):
            flip = False
        elif rb in (COMP[ea], COMP[oa]):
            flip = True
        else:
            return None, 'alleles_match_neither_strand'
    elif ea != rb and COMP[ea] == rb:
        return None, 'strand_ambiguous_single_allele'
    e = COMP[ea] if flip else ea
    return _site_effect_dosage(ctx, canon, w.pos, e, rb, 'flipped' if flip else ('ok' if oa else 'strand_assumed_forward'))


def _site_effect_dosage(ctx, canon, pos, e, rb, basis):
    calls = [c for c in ctx.overlapping_calls(canon, pos, pos) if c.dosage > 0]
    if calls:
        if any(c.vtype != 'snv' or c.pos != pos for c in calls):
            return None, 'overlapping_non_snv_call'
        if len({c.rec_index for c in calls}) > 1:
            return None, 'several_records_at_site'
        if any((not c.qc_pass) or c.missing for c in calls):
            return None, 'call_fails_quality'
        if e == rb:
            return calls[0].ref_copies, basis
        return sum(c.dosage for c in calls if c.alt == e), basis
    state, _ = ctx.callability.status(canon, pos, pos)
    if state == 'reference_confident':
        return (2 if e == rb else 0), basis
    return None, {'unknown': 'absent_from_callset', 'reference_low_confidence': 'reference_confidence_below_policy'}.get(state, state)


def _indel_dosage(w, canon, ea, oa, ctx, ref):
    if oa is None:
        return None, 'indel_without_other_allele'
    n = normalise(ref, canon, w.pos, oa, ea)
    eff_is_alt = True
    if not n.ok:
        n = normalise(ref, canon, w.pos, ea, oa)   # the weight may be written with the effect allele as REF
        eff_is_alt = False
        if not n.ok:
            return None, 'indel_alleles_do_not_match_reference'
    st = status_of(ctx, n.key)
    if st.status == 'EXACT_MATCH':
        c = st.call
        return (c.dosage if eff_is_alt else c.ref_copies), 'ok'
    if st.status == 'VERIFIED_REFERENCE':
        return (0 if eff_is_alt else 2), 'ok'
    return None, st.status.lower()
