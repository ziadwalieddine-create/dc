"""Orientation of a reported effect allele against the reference: the part of harmonisation that decides whether a record can be placed on
a normalised CHROM:POS:REF:ALT key. A record whose strand cannot be settled is rejected, never guessed.

Inputs are a reference base, the reported single-base effect allele (or unknown), the reported effect-allele frequency, and the alternate
alleles known at the site from an external allele source (Open Targets variant index) with their population frequencies.

  * plus strand is possible when the effect allele is the reference base or a known alternate allele
  * minus strand is possible when the complement of the effect allele is the reference base or a known alternate allele
  * when both are possible (a palindromic A/T or C/G site, or a multiallelic site holding an allele and its complement) the frequency decides,
    and only if the alternate allele is clearly minor or clearly major in every population (so ancestry cannot flip the answer)
"""
import glob
import os
from dataclasses import dataclass, field

COMP = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
MINOR, MAJOR = 0.4, 0.6


@dataclass
class Resolution:
    status: str                  # mapped | ambiguous | rejected | unmapped
    alt: object = None           # ALT allele of the key (forward strand)
    effect: object = None        # forward-strand effect allele, None when the source gives none
    other: object = None
    strand: str = 'unknown'
    reason: object = None
    flags: tuple = ()


def frequency_side(afs):
    """'minor' when the alternate allele is below 0.4 in every population, 'major' when above 0.6 in every population, else None."""
    afs = [a for a in (afs or []) if a is not None]
    if not afs:
        return None
    if max(afs) < MINOR:
        return 'minor'
    if min(afs) > MAJOR:
        return 'major'
    return None


def effect_is_alt_from_frequency(raf, afs):
    side = frequency_side(afs)
    if side is None or raf is None:
        return None
    if side == 'minor':
        return True if raf < MINOR else (False if raf > MAJOR else None)
    return True if raf > MAJOR else (False if raf < MINOR else None)


def resolve_snv(ref, risk, raf, alts, alt_afs=None):
    """ref: reference base; risk: reported effect allele ('A', 'C', 'G', 'T', '?' or None); raf: reported frequency of that allele or None;
    alts: set of single-base alternate alleles known at the site; alt_afs: {alt: [frequency per population]}."""
    if ref not in COMP:
        return Resolution('unmapped', reason='reference_base_not_ACGT')
    alts = {a for a in (alts or ()) if a in COMP and a != ref}
    if not alts:
        return Resolution('unmapped', reason='no_alternate_allele_known')
    if risk in (None, '', '?'):
        if len(alts) == 1:
            return Resolution('mapped', alt=next(iter(alts)), flags=('effect_allele_unknown',))
        return Resolution('ambiguous', reason='several_alternate_alleles_known', flags=('effect_allele_unknown',))
    x = risk.upper()
    if x not in COMP:
        return Resolution('unmapped', reason='effect_allele_not_a_single_base')
    comp = COMP[x]
    plus = x == ref or x in alts
    minus = comp == ref or comp in alts
    strand_note = ()
    if plus and not minus:
        strand, e = '+', x
    elif minus and not plus:
        strand, e = '-', comp
        strand_note = ('allele_reported_on_minus_strand',)
    elif not plus and not minus:
        return Resolution('rejected', reason='effect_allele_not_among_site_alleles')
    else:
        if len(alts) != 1:
            return Resolution('rejected', reason='strand_unresolved_multiallelic')
        a = next(iter(alts))
        eia = effect_is_alt_from_frequency(raf, (alt_afs or {}).get(a))
        if eia is None:
            return Resolution('rejected', reason='strand_unresolved_palindromic' if COMP[ref] == a else 'strand_unresolved_both_orientations_possible')
        # plus: effect allele is x; it is the alternate allele exactly when x != ref
        if (x != ref) == eia:
            strand, e = '+', x
        else:
            strand, e = '-', comp
        strand_note = ('strand_resolved_by_frequency',)
    if e == ref:
        if len(alts) == 1:
            a = next(iter(alts))
            return Resolution('mapped', alt=a, effect=e, other=a, strand=strand, flags=strand_note + ('effect_allele_is_reference',))
        return Resolution('ambiguous', reason='effect_allele_is_reference_several_alternate_alleles', effect=e, strand=strand, flags=strand_note)
    return Resolution('mapped', alt=e, effect=e, other=ref, strand=strand, flags=strand_note)


class OtVariantIndex:
    """Alternate alleles, rsIDs and population frequencies of the Open Targets variant index at requested positions."""

    def __init__(self, ot_dir):
        self.files = sorted(glob.glob(os.path.join(ot_dir, '*.parquet')))
        self.by_pos = {}

    def load(self, positions):
        """positions: iterable of (chrom, pos). Reads only the matching rows."""
        import pyarrow as pa
        import pyarrow.compute as pc
        import pyarrow.parquet as pq
        want = pa.array(sorted({f'{c}:{p}' for c, p in positions}), pa.string())
        for f in self.files:
            t = pq.read_table(f, columns=['variantId', 'chromosome', 'position', 'referenceAllele', 'alternateAllele', 'rsIds'])
            key = pc.binary_join_element_wise(t['chromosome'], pc.cast(t['position'], pa.string()), ':')
            t = t.filter(pc.is_in(key, value_set=want))
            for r in t.to_pylist():
                self.by_pos.setdefault((r['chromosome'], r['position']), []).append((r['referenceAllele'], r['alternateAllele'], tuple(r['rsIds'] or ()), r['variantId']))
        return len(self.by_pos)

    def snv_alts(self, chrom, pos, ref=None):
        return [(a, rs, vid) for r, a, rs, vid in self.by_pos.get((chrom, pos), ()) if len(r) == 1 and len(a) == 1 and (ref is None or r == ref)]

    def load_frequencies(self, variant_ids):
        import pyarrow as pa
        import pyarrow.compute as pc
        import pyarrow.parquet as pq
        want = pa.array(sorted(set(variant_ids)), pa.string())
        out = {}
        for f in self.files:
            t = pq.read_table(f, columns=['variantId', 'alleleFrequencies'])
            t = t.filter(pc.is_in(t['variantId'], value_set=want))
            for r in t.to_pylist():
                out[r['variantId']] = [x['alleleFrequency'] for x in (r['alleleFrequencies'] or [])]
        return out
