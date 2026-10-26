"""Genomic HGVS (g. / m.) changes -> normalised CHROM:POS:REF:ALT keys, with the reference sequence as the arbiter.

Only what can be decided from the reference sequence is turned into a key: substitutions, deletions, duplications, insertions and
deletion-insertions. The stated reference bases are checked against the FASTA (a conflict is reported, never trusted), '=' is the reference,
and anything that needs a repeat count or is not a plain sequence change comes back unresolved with a reason.
"""
import re

from ..normalise import normalise

LOC = re.compile(r'^([gm])\.(\d+)(?:_(\d+))?(.*)$')


def split_changes(text):
    """'g.18138969C>T; g.18138969C>G' -> [(start, end, rest), ...]; None when any part is not g./m. notation."""
    out = []
    for part in (text or '').split(';'):
        part = part.strip()
        if not part:
            continue
        m = LOC.match(part)
        if not m:
            return None
        out.append((int(m.group(2)), int(m.group(3) or m.group(2)), m.group(4)))
    return out


def _norm(ref, chrom, pos, r, a):
    n = normalise(ref, chrom, pos, r, a)
    return ('variant', [n.key], None) if n.ok else ('unresolved', [], n.status)


def change_to_key(ref, chrom, start, end, rest):
    """One change at start..end (1-based, inclusive) with the notation after the position in `rest`.
    Returns (status, keys, reason): 'variant' | 'reference' | 'unresolved'."""
    try:
        span = ref.fetch(chrom, start, end)
    except ValueError:
        return 'unresolved', [], 'outside_contig'
    if rest == '=':
        return 'reference', [], None
    eq = re.match(r'^([ACGT]+)=$', rest)
    if eq:
        return ('reference', [], None) if eq.group(1) == span else ('unresolved', [], f'reference_base_differs:{span}')
    sub = re.match(r'^([ACGT])>([ACGT])$', rest)
    if sub:
        if span != sub.group(1):
            return 'unresolved', [], f'reference_base_differs:{span}'
        return _norm(ref, chrom, start, sub.group(1), sub.group(2))
    if re.match(r'^[ACGT]>[A-Z]$', rest):
        return 'unresolved', [], 'alternate_not_a_defined_base'
    dl = re.match(r'^del([ACGT]*|\d+)$', rest)
    if dl:
        stated = dl.group(1)
        if stated.isdigit() and int(stated) != end - start + 1:
            return 'unresolved', [], 'deleted_length_differs_from_range'
        if stated and not stated.isdigit() and stated != span:
            return 'unresolved', [], 'deleted_sequence_differs_from_reference'
        if start == 1:
            return 'unresolved', [], 'deletion_at_contig_start'
        anchor = ref.fetch(chrom, start - 1, start - 1)
        return _norm(ref, chrom, start - 1, anchor + span, anchor)
    dp = re.match(r'^dup([ACGT]*|\d+)$', rest)
    if dp:
        stated = dp.group(1)
        if stated.isdigit() and int(stated) != end - start + 1:
            return 'unresolved', [], 'duplicated_length_differs_from_range'
        if stated and not stated.isdigit() and stated != span:
            return 'unresolved', [], 'duplicated_sequence_differs_from_reference'
        return _norm(ref, chrom, end, span[-1], span[-1] + span)
    ins = re.match(r'^ins([ACGT]+)$', rest)
    if ins:
        if end != start + 1:
            return 'unresolved', [], 'insertion_flanks_not_adjacent'
        b = ref.fetch(chrom, start, start)
        return _norm(ref, chrom, start, b, b + ins.group(1))
    di = re.match(r'^del([ACGT]*)ins([ACGT]+)$', rest)
    if di:
        if di.group(1) and di.group(1) != span:
            return 'unresolved', [], 'deleted_sequence_differs_from_reference'
        return _norm(ref, chrom, start, span, di.group(2))
    if '[' in rest or '(' in rest:
        return 'unresolved', [], 'repeat_length_notation'
    return 'unresolved', [], 'unsupported_notation'


def hgvs_to_key(ref, chrom, text):
    """A single-change genomic HGVS expression (the accession prefix already removed: 'g.123A>G') -> (status, keys, reason)."""
    changes = split_changes(text)
    if not changes:
        return 'unresolved', [], 'unparsed_notation'
    if len(changes) > 1:
        return 'unresolved', [], 'several_changes_in_one_expression'
    s, e, rest = changes[0]
    return change_to_key(ref, chrom, s, e, rest)
