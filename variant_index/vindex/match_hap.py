"""Haplotype (star-allele) calling from defining variants, with explicit evidence requirements.

A diplotype is only called when EVERY defining position of the gene has positive evidence: a called genotype, or a verified reference
(an absent record in a variant-only VCF is not reference). Missing positions return UNRESOLVED with the positions and the states
found, plus the diplotypes still consistent with the known positions. Unphased heterozygous combinations that fit several diplotypes
return AMBIGUOUS with all of them; phase (PS / read-backed) resolves them only when every heterozygous defining call is in one phase
set. Genes with structural alleles (deletions, duplications, hybrids) need a dedicated caller; their result is REQUIRES_SPECIALISED_CALLER
unless the caller's diplotype is supplied through `external`.
Reuse note: PharmCAT (CPIC-endorsed) is the production engine; this module is the verified reference implementation used on fixtures
and for explaining why a call could not be made.
"""
from dataclasses import dataclass, field
from itertools import combinations_with_replacement

from .match_small import status_of

KNOWN_STRUCTURAL_GENES = {'CYP2D6', 'CYP2A6', 'GSTM1', 'GSTT1', 'UGT2B17', 'SULT1A1', 'HLA-A', 'HLA-B', 'HLA-C', 'HLA-DRB1', 'HLA-DQB1', 'HLA-DPB1', 'CYP2B6'}


@dataclass
class HaplotypeCall:
    gene: str
    system: str
    status: str                      # RESOLVED | AMBIGUOUS | UNRESOLVED | NO_CONSISTENT_DIPLOTYPE | REQUIRES_SPECIALISED_CALLER
    diplotypes: list = field(default_factory=list)
    basis: str = ''
    missing: dict = field(default_factory=dict)
    flags: tuple = ()
    position_states: dict = field(default_factory=dict)


def _alt_dosage(st):
    if st.status == 'EXACT_MATCH':
        return st.dosage
    if st.status in ('VERIFIED_REFERENCE', 'OTHER_ALLELE_AT_SITE'):
        return 0
    return None


def call_gene(ctx, index, gene, system=None, external=None, specialised_genes=KNOWN_STRUCTURAL_GENES):
    haps = index.haplotypes(gene, system)
    if not haps:
        return HaplotypeCall(gene, system or '', 'UNRESOLVED', basis='no_haplotype_definitions_in_index', missing={'definitions': gene})
    sysname = haps[0]['system']
    defs = index.haplotype_definitions([h['haplotype_id'] for h in haps])
    structural = any(h.get('structural') for h in haps) or gene in specialised_genes
    if structural:
        if external and gene in external:
            return HaplotypeCall(gene, sysname, 'RESOLVED', [tuple(external[gene])], 'external_specialised_caller', flags=('not_derived_from_the_variant_calls',))
        return HaplotypeCall(gene, sysname, 'REQUIRES_SPECIALISED_CALLER', basis='structural_alleles',
                             missing={'evidence': 'copy number / hybrid-allele calls from a dedicated caller (for example Cyrius, Aldy, StellarPGx) or targeted long reads'})
    sets = {h['haplotype_id']: frozenset(d['variant_key'] for d in defs.get(h['haplotype_id'], ())) for h in haps}
    universe = sorted(set().union(*sets.values())) if sets else []
    states = {u: status_of(ctx, u) for u in universe}
    dos = {u: _alt_dosage(s) for u, s in states.items()}
    missing = {u: states[u].status for u in universe if dos[u] is None}
    other = [u for u, s in states.items() if s.status == 'OTHER_ALLELE_AT_SITE']
    flags = tuple(['other_allele_at_defining_site:' + ','.join(other)] if other else [])
    ploidy = 2
    for s in states.values():
        if s.call is not None:
            ploidy = s.call.ploidy
            break
    pos_states = {u: (s.status, s.dosage) for u, s in states.items()}
    names = {h['haplotype_id']: h['allele_name'] for h in haps}

    def pairs_consistent(known):
        out = []
        for combo in combinations_with_replacement(sorted(sets), ploidy):
            if all(sum(1 for h in combo if u in sets[h]) == d for u, d in known.items()):
                out.append(tuple(sorted(names[h] for h in combo)))
        return sorted(set(out))

    known = {u: d for u, d in dos.items() if d is not None}
    if missing:
        return HaplotypeCall(gene, sysname, 'UNRESOLVED', pairs_consistent(known), 'missing_defining_positions',
                             {'positions': missing, 'evidence': 'a genotype call or reference-confidence evidence at each listed position'}, flags, pos_states)
    cands = pairs_consistent(known)
    # phase: every heterozygous defining call must carry haplotype assignments from one phase set
    hets = [(k, s) for k, s in states.items() if s.status == 'EXACT_MATCH' and 0 < s.dosage < ploidy]
    homs = [k for k, s in states.items() if s.status == 'EXACT_MATCH' and s.dosage == ploidy]
    phased_ok = bool(hets) and all(s.call is not None and s.call.hap_alleles is not None for _, s in hets) and len({str(s.call.phase_set) for _, s in hets}) == 1
    if phased_ok:
        obs = [frozenset([k for k, s in hets if s.call.hap_alleles[i] == 1] + homs) for i in range(ploidy)]
        by_name = {names[h]: sets[h] for h in sets}
        want = _multiset(obs)
        keep = [d for d in cands if _multiset(by_name[n] for n in d) == want]
        if not keep:
            return HaplotypeCall(gene, sysname, 'NO_CONSISTENT_DIPLOTYPE', [], 'phased_haplotypes_match_no_defined_pair',
                                 {'evidence': 'a defined haplotype carrying the observed cis combination (novel or undefined allele)', 'observed_haplotypes': [sorted(o) for o in obs]}, flags, pos_states)
        return HaplotypeCall(gene, sysname, 'RESOLVED' if len(keep) == 1 else 'AMBIGUOUS', keep, 'phased', {}, flags, pos_states)
    if not cands:
        return HaplotypeCall(gene, sysname, 'NO_CONSISTENT_DIPLOTYPE', [], 'defining_genotypes_match_no_pair_of_defined_haplotypes',
                             {'evidence': 'a defined haplotype that carries the observed variant combination (novel or undefined allele)'}, flags, pos_states)
    if len(cands) == 1:
        return HaplotypeCall(gene, sysname, 'RESOLVED', cands, 'unique_consistent_pair', {}, flags, pos_states)
    return HaplotypeCall(gene, sysname, 'AMBIGUOUS', cands, 'unphased_heterozygous_combination',
                         {'evidence': 'phase (PS or read-backed) linking the heterozygous defining calls'}, flags, pos_states)


def _multiset(sets_iter):
    return sorted(tuple(sorted(s)) for s in sets_iter)
