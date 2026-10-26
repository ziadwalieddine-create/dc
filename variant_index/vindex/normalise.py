"""Allele normalisation: left-aligned, parsimonious, validated against the reference.

Identity of a small variant is the normalised key 'CHROM:POS:REF:ALT' (canonical contig name, 1-based POS). rsIDs and other
database identifiers are cross-references, never identity. Algorithm (Tan et al. 2015): trim the common suffix; when an
allele becomes empty, extend both alleles one reference base to the left; repeat; then trim the common prefix keeping one
anchor base. Equal-length multi-base substitutions (MNVs) keep their identity; `atomise` offers their per-base view for lookups.
"""
import re
from dataclasses import dataclass

_ACGT = re.compile(r'^[ACGT]+$')


@dataclass(frozen=True)
class Norm:
    chrom: str
    pos: int
    ref: str
    alt: str
    vtype: str      # snv | mnv | ins | del | complex | '' when not a normalised small variant
    status: str     # ok | ref_mismatch | ref_n | invalid_allele | not_a_variant | symbolic | no_alt_allele | spanning_deletion | out_of_contig
    shifted: bool   # the normalised representation differs from the input
    detail: str = ''

    @property
    def key(self):
        return f'{self.chrom}:{self.pos}:{self.ref}:{self.alt}'

    @property
    def ok(self):
        return self.status == 'ok'


def classify(ref, alt):
    if len(ref) == 1 and len(alt) == 1:
        return 'snv'
    if len(ref) == len(alt):
        return 'mnv'
    if len(ref) == 1 and alt[0] == ref:
        return 'ins'
    if len(alt) == 1 and ref[0] == alt:
        return 'del'
    return 'complex'


def parse_key(key):
    c, p, r, a = key.split(':')
    return c, int(p), r, a


def normalise(reference, chrom, pos, ref, alt):
    """Normalise one ALT allele. `chrom` must already be the canonical contig name. Never raises for bad input: returns a status."""
    pos0, r0, a0 = pos, ref.upper(), alt.upper()
    ref, alt = ref.upper(), alt.upper()
    if alt == '*':
        return Norm(chrom, pos, ref, alt, '', 'spanning_deletion', False, 'ALT * marks a haplotype deleted by an upstream variant')
    if alt in ('.', ''):
        return Norm(chrom, pos, ref, alt, '', 'no_alt_allele', False, 'ALT is missing')
    if alt.startswith('<') or '[' in alt or ']' in alt:
        return Norm(chrom, pos, ref, alt, '', 'symbolic', False)
    if not ref or not _ACGT.match(ref) or not _ACGT.match(alt):
        if set(ref) <= set('ACGTN') and set(alt) <= set('ACGTN') and ref and alt and 'N' in ref + alt:
            return Norm(chrom, pos, ref, alt, '', 'ref_n', False, 'N in allele')
        return Norm(chrom, pos, ref, alt, '', 'invalid_allele', False, 'allele outside A/C/G/T')
    state, seq = reference.check_ref(chrom, pos, ref)
    if state != 'match':
        detail = {'mismatch': f'REF {ref} but the reference has {seq}', 'n': 'N in the reference or allele: not verifiable',
                  'out_of_contig': 'REF span leaves the contig'}[state]
        return Norm(chrom, pos, ref, alt, '', {'mismatch': 'ref_mismatch', 'n': 'ref_n', 'out_of_contig': 'out_of_contig'}[state], False, detail)
    if ref == alt:
        return Norm(chrom, pos, ref, alt, '', 'not_a_variant', False)
    if len(ref) == 1 and len(alt) == 1:
        return Norm(chrom, pos, ref, alt, 'snv', 'ok', False)
    detail = ''
    while True:
        changed = False
        if ref and alt and ref[-1] == alt[-1]:
            ref, alt, changed = ref[:-1], alt[:-1], True
        if not ref or not alt:
            if pos == 1:  # cannot extend left of the contig: anchor on the right instead
                nxt = reference.fetch(chrom, pos + len(ref), pos + len(ref))
                ref, alt = ref + nxt, alt + nxt
                detail = 'right_anchored_at_contig_start'
                break
            pos -= 1
            b = reference.fetch(chrom, pos, pos)
            ref, alt, changed = b + ref, b + alt, True
        if not changed:
            break
    while len(ref) > 1 and len(alt) > 1 and ref[0] == alt[0]:
        ref, alt, pos = ref[1:], alt[1:], pos + 1
    return Norm(chrom, pos, ref, alt, classify(ref, alt), 'ok', (pos, ref, alt) != (pos0, r0, a0), detail)


def atomise(n):
    """Per-base SNV view of an MNV (other types are returned unchanged). The parts are weaker evidence than the MNV itself: a lookup
    that matches a part must say so."""
    if n.vtype != 'mnv':
        return [n]
    return [Norm(n.chrom, n.pos + i, r, a, 'snv', 'ok', True, 'atomised_from_mnv') for i, (r, a) in enumerate(zip(n.ref, n.alt)) if r != a]
