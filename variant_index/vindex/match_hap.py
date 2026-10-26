"""Haplotype (star-allele) calling from defining variants, with explicit evidence requirements.

A diplotype is only called when EVERY defining position of the gene has positive evidence: a called genotype, or a verified reference
(an absent record in a variant-only VCF is not reference). Missing positions return UNRESOLVED with the positions and the states
found, plus the diplotypes still consistent with the known positions. Unphased heterozygous combinations that fit several diplotypes
return AMBIGUOUS with all of them; phase (PS / read-backed) resolves them only when every heterozygous defining call is in one phase
set. Genes with structural alleles (deletions, duplications, hybrids) need a dedicated caller; their result is REQUIRES_SPECIALISED_CALLER
unless the caller's diplotype is supplied through `external`.
A definition row with scope 'tolerated' (an IUPAC code in the source: the allele may carry the alternative or the reference base) never constrains
the dosage at its position; it widens the allowed range for that allele instead.
Alleles whose definition the index could not turn into normalised variants (repeat-count notation, IUPAC ambiguity codes, no location rows,
reference-base conflicts) are never silently dropped: they are listed in `unevaluable_alleles` with the reasons. A unique diplotype among
the evaluable alleles is then RESOLVED_CONDITIONAL, not RESOLVED, because those alleles could not be excluded.
Reuse note: PharmCAT (CPIC-endorsed) is the production engine; this module is the verified reference implementation used on fixtures
and for explaining why a call could not be made.
"""
from dataclasses import dataclass, field
from itertools import combinations_with_replacement, permutations

from .match_small import status_of

KNOWN_STRUCTURAL_GENES = {'CYP2D6', 'CYP2A6', 'GSTM1', 'GSTT1', 'UGT2B17', 'SULT1A1', 'HLA-A', 'HLA-B', 'HLA-C', 'HLA-DRB1', 'HLA-DQB1', 'HLA-DPB1', 'CYP2B6'}
NON_STAR_SYSTEMS = {'clinvar'}      # haplotype systems whose rows are classified combinations of alleles, never candidate star alleles for a diplotype call


@dataclass
class HaplotypeCall:
    gene: str
    system: str
    status: str                      # RESOLVED | RESOLVED_CONDITIONAL | AMBIGUOUS | UNRESOLVED | NO_CONSISTENT_DIPLOTYPE | REQUIRES_SPECIALISED_CALLER
    diplotypes: list = field(default_factory=list)
    basis: str = ''
    missing: dict = field(default_factory=dict)
    flags: tuple = ()
    position_states: dict = field(default_factory=dict)
    unevaluable_alleles: dict = field(default_factory=dict)      # allele name -> reasons its definition could not be evaluated


def _alt_dosage(st):
    if st.status == 'EXACT_MATCH':
        return st.dosage
    if st.status in ('VERIFIED_REFERENCE', 'OTHER_ALLELE_AT_SITE'):
        return 0
    return None


def call_gene(ctx, index, gene, system=None, external=None, specialised_genes=KNOWN_STRUCTURAL_GENES, levels=None):
    """levels: restrict to allele levels (for example ('core',) for PharmVar so that sub-alleles do not multiply the candidates)."""
    if system in NON_STAR_SYSTEMS:
        return HaplotypeCall(gene, system, 'UNRESOLVED', basis='system_is_not_a_star_allele_system',
                             missing={'evidence': f'{system} composite records are classified combinations of alleles (see match_composite), not star-allele definitions to call a diplotype from'})
    all_haps = [h for h in index.haplotypes(gene, system) if h['system'] not in NON_STAR_SYSTEMS]
    if not all_haps:
        return HaplotypeCall(gene, system or '', 'UNRESOLVED', basis='no_haplotype_definitions_in_index', missing={'definitions': gene})
    present = {h['system'] for h in all_haps}
    sysname = system or next(iter([x for x in ('cpic', 'pharmvar') if x in present] + sorted(present)))
    if levels is None and sysname == 'pharmvar':
        levels = ('core',)           # a sub-allele repeats its core allele's variants and adds its own: both together are not distinguishable candidates
    haps = [h for h in all_haps if h['system'] == sysname and (not levels or h.get('allele_level') in levels)]
    if not haps:
        return HaplotypeCall(gene, sysname, 'UNRESOLVED', basis='no_haplotype_definitions_at_requested_level', missing={'definitions': gene, 'levels': levels})
    defs = index.haplotype_definitions([h['haplotype_id'] for h in haps])
    structural = any(h.get('structural') for h in haps) or gene in specialised_genes
    if structural:
        if external and gene in external:
            return HaplotypeCall(gene, sysname, 'RESOLVED', [tuple(external[gene])], 'external_specialised_caller', flags=('not_derived_from_the_variant_calls',))
        return HaplotypeCall(gene, sysname, 'REQUIRES_SPECIALISED_CALLER', basis='structural_alleles',
                             missing={'evidence': 'copy number / hybrid-allele calls from a dedicated caller (for example Cyrius, Aldy, StellarPGx) or targeted long reads'})
    unevaluable = {}
    for h in haps:
        why = sorted({(d.get('definition_scope') or '')[len('unresolved:'):] for d in defs.get(h['haplotype_id'], ()) if (d.get('definition_scope') or '').startswith('unresolved')})
        if why:
            unevaluable[h['allele_name']] = why
    haps = [h for h in haps if h['allele_name'] not in unevaluable]
    scope = lambda d: d.get('definition_scope') or ''
    sets = {h['haplotype_id']: frozenset(d['variant_key'] for d in defs.get(h['haplotype_id'], ()) if d['variant_key'] and scope(d) != 'tolerated') for h in haps}
    tol = {h['haplotype_id']: frozenset(d['variant_key'] for d in defs.get(h['haplotype_id'], ()) if d['variant_key'] and scope(d) == 'tolerated') - sets[h['haplotype_id']] for h in haps}
    if not sets:
        return HaplotypeCall(gene, sysname, 'UNRESOLVED', basis='no_allele_definition_could_be_evaluated', missing={'allele_definitions': unevaluable}, unevaluable_alleles=unevaluable)
    universe = sorted(set().union(*sets.values())) if sets else []
    states = {u: status_of(ctx, u) for u in universe}
    dos = {u: _alt_dosage(s) for u, s in states.items()}
    missing = {u: states[u].status for u in universe if dos[u] is None}
    other = [u for u, s in states.items() if s.status == 'OTHER_ALLELE_AT_SITE']
    flags = tuple((['other_allele_at_defining_site:' + ','.join(other)] if other else []) + ([f'system_defaulted:{sysname}_of:' + ','.join(sorted(present))] if system is None and len(present) > 1 else []))
    ploidy = 2
    for s in states.values():
        if s.call is not None:
            ploidy = s.call.ploidy
            break
    pos_states = {u: (s.status, s.dosage) for u, s in states.items()}
    names = {h['haplotype_id']: h['allele_name'] for h in haps}

    def done(status, diplotypes, basis, miss=None):
        if status == 'RESOLVED' and unevaluable:
            status = 'RESOLVED_CONDITIONAL'
            miss = dict(miss or {}, not_excluded_alleles=unevaluable)
        elif unevaluable and status in ('UNRESOLVED', 'AMBIGUOUS', 'NO_CONSISTENT_DIPLOTYPE'):
            miss = dict(miss or {}, unevaluable_allele_definitions=unevaluable)
        return HaplotypeCall(gene, sysname, status, diplotypes, basis, miss or {}, flags, pos_states, unevaluable)

    def pairs_consistent(known):
        out = []
        for combo in combinations_with_replacement(sorted(sets), ploidy):
            for u, d in known.items():
                lo = sum(1 for h in combo if u in sets[h])
                if not lo <= d <= lo + sum(1 for h in combo if u in tol[h]):
                    break
            else:
                out.append(tuple(sorted(names[h] for h in combo)))
        return sorted(set(out))

    known = {u: d for u, d in dos.items() if d is not None}
    if missing:
        return done('UNRESOLVED', pairs_consistent(known), 'missing_defining_positions',
                    {'positions': missing, 'evidence': 'a genotype call or reference-confidence evidence at each listed position'})
    cands = pairs_consistent(known)
    # phase: every heterozygous defining call must carry haplotype assignments from one phase set
    hets = [(k, s) for k, s in states.items() if s.status == 'EXACT_MATCH' and 0 < s.dosage < ploidy]
    homs = [k for k, s in states.items() if s.status == 'EXACT_MATCH' and s.dosage == ploidy]
    phased_ok = bool(hets) and all(s.call is not None and s.call.hap_alleles is not None for _, s in hets) and len({str(s.call.phase_set) for _, s in hets}) == 1
    if phased_ok:
        obs = [frozenset([k for k, s in hets if s.call.hap_alleles[i] == 1] + homs) for i in range(ploidy)]
        by_name = {names[h]: (sets[h], tol[h]) for h in sets}
        keep = [d for d in cands if any(all(by_name[n][0] <= o <= by_name[n][0] | by_name[n][1] for n, o in zip(perm, obs)) for perm in set(permutations(d)))]
        if not keep:
            return done('NO_CONSISTENT_DIPLOTYPE', [], 'phased_haplotypes_match_no_defined_pair',
                        {'evidence': 'a defined haplotype carrying the observed cis combination (novel or undefined allele)', 'observed_haplotypes': [sorted(o) for o in obs]})
        return done('RESOLVED' if len(keep) == 1 else 'AMBIGUOUS', keep, 'phased')
    if not cands:
        return done('NO_CONSISTENT_DIPLOTYPE', [], 'defining_genotypes_match_no_pair_of_defined_haplotypes',
                    {'evidence': 'a defined haplotype that carries the observed variant combination (novel or undefined allele)'})
    if len(cands) == 1:
        return done('RESOLVED', cands, 'unique_consistent_pair')
    return done('AMBIGUOUS', cands, 'unphased_heterozygous_combination', {'evidence': 'phase (PS or read-backed) linking the heterozygous defining calls'})


def _multiset(sets_iter):
    return sorted(tuple(sorted(s)) for s in sets_iter)
