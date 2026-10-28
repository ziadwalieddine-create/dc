"""Targeted read-level evidence for small variants at the loci a callset cannot answer for.

A variant-only VCF says nothing about a position it does not list, and a no-call or a filtered call says little more. When the sample's aligned reads (BAM or CRAM) are at hand they
show whether the indexed allele is carried at a locus whatever the caller did. For one indexed key (CHROM:POS:REF:ALT) this module counts the fragments that carry the reference
allele, the indexed alternate allele or some other allele, and states what the counts allow:

  READS_ALTERNATE_PRESENT     at least min_alt_reads fragments carry the indexed allele and they are at least min_alt_fraction of the informative fragments
  READS_OTHER_ALLELE_PRESENT  none does, but at least min_alt_reads fragments have a sequence over the window that is neither the reference's nor the indexed alternate's, at that
                              fraction: another allele, or the indexed allele with another variant of the same haplotype inside the window (an indel beside an SNV, in a repeat), which
                              the comparison of whole windows cannot tell apart; the sequences seen are kept in the detail
  READS_REFERENCE_ONLY        at least min_depth informative fragments, none carries the indexed allele, and no caveat applies (below)
  READS_AMBIGUOUS             the reads neither show nor exclude the allele: a few alternate fragments below the thresholds, or a caveat that forbids reading their absence as
                              reference (soft-clipped reads ending at the locus, a deletion over it, overlapping mates that disagree, low-mapping-quality reads that dominate or
                              look like the alternate allele, alternate-like reads below the base-quality threshold)
  READS_INSUFFICIENT          fewer than min_depth informative fragments (none overlap the locus, or the filters removed them) and no allele shown present
  READS_NOT_EVALUATED         the key cannot be asked of reads here: an allele longer than max_allele_length, a repeat tract longer than max_window, a REF that disagrees with the
                              reference, an N in the window, a contig the alignment does not carry, more alignments than max_reads

The presence of an allele is shown by a few reliable fragments; its absence needs depth and the absence of every caveat. The asymmetry is deliberate.

What the allele of a read is. The read's own sequence over a window of the reference is compared with that window as it reads under the reference allele and under the indexed
alternate allele. For an insertion or a deletion the window reaches `flank` bases beyond the last position an equivalent placement of the event can occupy, so a read in which the
aligner put the gap at the right end of a repeat is counted as the allele it is, not as another allele because its CIGAR differs from the left-aligned key. A read counts only if
its alignment covers the window from end to end with aligned bases at both ends. A read whose sequence over the window matches neither is 'other' (the sequences seen are kept).

What counts as a read. Unmapped, secondary, supplementary, QC-fail and duplicate-flagged alignments are excluded; so are alignments below min_mapq and windows with a base below
min_baseq (a missing quality fails). Both mates of a pair that cover the window are one fragment; mates that disagree are set aside and counted. Every exclusion is counted by reason
in the detail, and the alignments set aside for mapping or base quality are classified by allele and counted too, so that a filter cannot turn evidence into silence.

What this is not. Reads are not a genotyper: no zygosity, copy number, phase (beyond the comparison of one window) or contamination model is attempted, and the default policy is
tuned for short reads. A fragment count is not a dosage. With min_depth informative fragments a heterozygous site shows no alternate fragment with probability 2^-min_depth if
fragments are sampled without allele bias (0.1% at 10); alignment bias, which the caveats only partly catch, is not modelled. Structural, copy-number and repeat-length evidence from
reads is not implemented. Duplicates are excluded only where the aligner or a marking tool flagged them.
"""
import json
from collections import Counter
from dataclasses import dataclass, field

import pysam

from .normalise import classify, parse_key

NOT_EVALUATED = 'READS_NOT_EVALUATED'
INSUFFICIENT = 'READS_INSUFFICIENT'
REFERENCE_ONLY = 'READS_REFERENCE_ONLY'
ALTERNATE_PRESENT = 'READS_ALTERNATE_PRESENT'
OTHER_ALLELE_PRESENT = 'READS_OTHER_ALLELE_PRESENT'
AMBIGUOUS = 'READS_AMBIGUOUS'
READ_STATES = (NOT_EVALUATED, INSUFFICIENT, REFERENCE_ONLY, ALTERNATE_PRESENT, OTHER_ALLELE_PRESENT, AMBIGUOUS)
ANSWERS = (REFERENCE_ONLY, ALTERNATE_PRESENT, OTHER_ALLELE_PRESENT)           # the states that answer whether the indexed allele is carried

# What the reads can settle: the presence of the indexed allele, which is what a variant-only file, a no-call, a filtered call and an overlapping call that fails quality leave open.
# Not its copy number or its phase, and not the structural questions (a site inside a called deletion or duplication), which reads at this resolution do not answer.
PRESENCE_UNRESOLVED_REASONS = frozenset(('reference_confidence_below_policy', 'overlapping_call_fails_quality', 'some_mnv_components_present_others_unknown', 'phase_required'))
PRESENCE_STATUSES = frozenset(('ABSENT_FROM_CALLSET', 'NO_CALL', 'FILTERED_CALL'))


class ReadsAssemblyMismatch(ValueError):
    pass


@dataclass(frozen=True)
class ReadPolicy:
    min_mapq: int = 20
    min_baseq: int = 20
    min_depth: int = 10                    # informative fragments needed before the absence of the allele is read as reference
    min_alt_reads: int = 3                 # fragments needed before an allele counts as present
    min_alt_fraction: float = 0.2          # of the informative fragments
    max_allele_length: int = 50            # longer REF or ALT alleles are structural-variant questions
    max_window: int = 140                  # reference bases a window may span: a read has to cover all of them
    flank: int = 2                         # bases beyond an insertion, a deletion or a complex allele (an SNV or an MNV is compared on its own bases)
    clip_slack: int = 2                    # a clipped read ends 'at the locus' when its clip lies within this many bases of the window
    max_clipped_fraction: float = 0.2      # clipped reads ending at the locus, of informative plus clipped, that forbid reading absence as reference
    max_low_mapq_fraction: float = 0.25    # alignments below min_mapq, of those covering the locus that are informative or below min_mapq, above which absence is not read as reference
    max_reads: int = 100_000               # alignments inspected per locus: above this the locus is a collapsed or pile-up region and is not evaluated

    def __post_init__(self):
        if not (self.min_depth >= 1 and self.min_alt_reads >= 1 and 0 < self.min_alt_fraction <= 1 and self.flank >= 1 and self.max_window >= 3 and self.max_reads >= 1
                and self.min_mapq >= 0 and self.min_baseq >= 0 and self.max_allele_length >= 1):
            raise ValueError(f'invalid read policy {self}')


@dataclass
class ReadSupport:
    key: str
    state: str
    reasons: tuple = ()
    flags: tuple = ()
    depth: int = 0                         # informative fragments (reference + alternate + other)
    ref_reads: int = 0
    alt_reads: int = 0
    other_reads: int = 0
    detail: dict = field(default_factory=dict)

    @property
    def alt_fraction(self):
        return round(self.alt_reads / self.depth, 6) if self.depth else None

    def columns(self, status=None, reasons=()):
        """The unresolved_loci columns for this evidence; `status` and `reasons` are the matcher's own for the locus (they decide what the reads settle)."""
        return dict(read_state=self.state, read_reasons=list(self.reasons) or None, read_flags=list(self.flags) or None, read_depth=self.depth, read_ref=self.ref_reads,
                    read_alt=self.alt_reads, read_other=self.other_reads, read_alt_fraction=self.alt_fraction, read_settles=settles(status, reasons, self.state),
                    read_interpretation=interpretation(status, self.state), read_detail_json=json.dumps(self.detail, sort_keys=True) if self.detail else None)


def asks_presence(status, reasons=()):
    """Whether the callset's open question at a locus is whether the indexed allele is carried (which reads can answer)."""
    return status in PRESENCE_STATUSES or (status == 'UNRESOLVED' and bool(PRESENCE_UNRESOLVED_REASONS & set(reasons or ())))


def settles(status, reasons, state):
    return bool(status) and asks_presence(status, reasons) and state in ANSWERS


def interpretation(status, state):
    """A short code for what the reads say about the callset's status at the locus."""
    if state == NOT_EVALUATED:
        return 'not_evaluated'
    if state == INSUFFICIENT:
        return 'reads_do_not_settle_it'
    if state == AMBIGUOUS:
        return 'reads_ambiguous'
    if status == 'FILTERED_CALL':
        return {ALTERNATE_PRESENT: 'reads_support_the_filtered_call', REFERENCE_ONLY: 'reads_contradict_the_filtered_call', OTHER_ALLELE_PRESENT: 'reads_show_another_sequence_over_the_window'}[state]
    return {ALTERNATE_PRESENT: 'reads_carry_the_allele_the_callset_does_not_list', REFERENCE_ONLY: 'reads_show_the_reference_allele_only', OTHER_ALLELE_PRESENT: 'reads_show_another_sequence_over_the_window'}[state]


@dataclass(frozen=True)
class Locus:
    key: str
    chrom: str
    pos: int
    ref: str
    alt: str
    vtype: str
    start: int                             # 1-based window, inclusive
    end: int
    ref_window: str
    alt_window: str


def plan_locus(reference, key, policy):
    """(Locus, None) or (None, reason): the window of the reference to compare reads over and the sequences it has under the reference and the alternate allele."""
    chrom, pos, ref, alt = parse_key(key)
    if chrom not in reference.lengths:
        return None, 'contig_not_in_the_reference'
    if max(len(ref), len(alt)) > policy.max_allele_length:
        return None, 'allele_longer_than_the_read_policy_allows'
    vtype = classify(ref, alt)
    length = reference.length(chrom)
    if pos < 1 or pos + len(ref) - 1 > length:
        return None, 'locus_outside_the_contig'
    last = pos + len(ref) - 1
    flank = 0 if vtype in ('snv', 'mnv') else policy.flank
    if vtype in ('ins', 'del'):
        tail = reference.fetch(chrom, pos, min(length, pos + policy.max_window))      # tail[0] is the base at pos, the anchor of the allele
        capped = pos + policy.max_window < length
        if vtype == 'del':
            s, t = 1, len(ref) - 1                                   # the deleted span as offsets into tail; it moves right while the base after it repeats its first base
            while t + 1 < len(tail) and tail[t + 1] == tail[s]:
                s, t = s + 1, t + 1
            ran_out = t + 1 >= len(tail)
            last = pos + t
        else:
            ins, p = alt[1:], 0                                      # the insertion sits after tail[p]; it moves right while the next base repeats its first base
            while p + 1 < len(tail) and tail[p + 1] == ins[0]:
                ins, p = ins[1:] + ins[0], p + 1
            ran_out = p + 1 >= len(tail)
            last = pos + p + 1                                       # the base after the last possible placement
        if ran_out and capped:
            return None, 'repeat_tract_longer_than_the_window_limit'
    start, end = max(1, pos - flank), min(length, last + flank)
    if end - start + 1 > policy.max_window:
        return None, 'window_longer_than_the_window_limit'
    window = reference.fetch(chrom, start, end)
    if 'N' in window:
        return None, 'reference_has_N_in_the_window'
    off = pos - start
    if window[off:off + len(ref)] != ref:
        return None, 'REF_disagrees_with_the_reference'
    return Locus(key, chrom, pos, ref, alt, vtype, start, end, window, window[:off] + alt + window[off + len(ref):]), None


def read_window(a, w0, w1):
    """(query_start, query_end inclusive, None) for the query bases aligned to the 0-based reference positions w0 and w1 (inserted bases between them lie inside), or
    (None, None, reason): the alignment does not cover the window, an end of the window lies in a deletion or a reference skip, or a reference skip lies inside it."""
    cig = a.cigartuples
    if not cig:
        return None, None, 'no_alignment'
    r, q, qs, qe, skip = a.reference_start, 0, None, None, False
    for op, n in cig:
        if op in (0, 7, 8):                                          # M = X: aligned bases
            if r <= w0 < r + n:
                qs = q + (w0 - r)
            if r <= w1 < r + n:
                qe = q + (w1 - r)
            r, q = r + n, q + n
        elif op in (1, 4):                                           # I, S: query only
            q += n
        elif op in (2, 3):                                           # D, N: reference only
            if op == 3 and r <= w1 and r + n - 1 >= w0:
                skip = True
            r += n
        if r > w1 and qe is not None:
            break
    if skip:
        return None, None, 'reference_skip_in_the_window'
    if qs is None or qe is None:
        return None, None, 'not_spanning_the_window' if (w0 < a.reference_start or w1 >= a.reference_end) else 'window_end_in_a_gap'
    return qs, qe, None


def clipped_at(a, w0, w1, slack):
    """Whether the alignment is clipped at an end that lies within slack bases of the 0-based window w0..w1: the clipped part may hold the allele the aligner could not place."""
    cig = a.cigartuples
    if not cig:
        return False
    if cig[0][0] in (4, 5) and w0 - slack <= a.reference_start <= w1 + slack:
        return True
    return cig[-1][0] in (4, 5) and w0 - slack <= a.reference_end - 1 <= w1 + slack


class ReadEvidence:
    """An aligned BAM or CRAM (indexed), read through pysam at the loci asked of it. A CRAM is decoded against the reference FASTA the matcher uses. `index` names the index file
    when it is not beside the alignment (a public alignment read over HTTPS with a local copy of its index, for example)."""

    def __init__(self, path, reference, policy=None, allow_assembly_mismatch=False, index=None):
        self.path, self.reference, self.policy = path, reference, policy or ReadPolicy()
        self.af = pysam.AlignmentFile(path, 'r', reference_filename=reference.path, index_filename=index, require_index=True)
        self.assembly = reference.assembly_check(dict(zip(self.af.references, self.af.lengths)))
        if not self.assembly.get('ok', True) and not allow_assembly_mismatch:
            self.af.close()
            raise ReadsAssemblyMismatch(f'the alignment header does not match the reference ({self.assembly}); counting reads on a different assembly would give wrong answers')
        self._contig = {}
        for name, n in zip(self.af.references, self.af.lengths):
            canon = reference.contig(name)
            if canon is not None and reference.length(canon) == n:
                self._contig.setdefault(canon, name)
        self.format = self.af.format.decode() if isinstance(self.af.format, bytes) else str(self.af.format)

    def close(self):
        self.af.close()

    def __enter__(self):
        return self

    def __exit__(self, *exc):
        self.close()

    def evaluate(self, key):
        pol = self.policy
        locus, why = plan_locus(self.reference, key, pol)
        if locus is None:
            return ReadSupport(key, NOT_EVALUATED, (why,))
        name = self._contig.get(locus.chrom)
        if name is None:
            return ReadSupport(key, NOT_EVALUATED, ('contig_not_in_the_alignment_header',))
        w0, w1 = locus.start - 1, locus.end - 1                      # 0-based, inclusive
        tally = dict(filtered=Counter(), low_baseq=Counter(), low_mapq=Counter(), frags={}, clipped=0, informative_alignments=0)
        n = 0
        for a in self.af.fetch(name, w0, w1 + 1):
            n += 1
            if n > pol.max_reads:
                return ReadSupport(key, NOT_EVALUATED, ('more_alignments_than_the_read_policy_allows',), detail=dict(window=[locus.chrom, locus.start, locus.end]))
            skip = ('unmapped' if a.is_unmapped else 'secondary' if a.is_secondary else 'supplementary' if a.is_supplementary else 'qc_fail' if a.is_qcfail
                    else 'duplicate' if a.is_duplicate else 'no_sequence' if not a.query_sequence else None)
            if skip:
                tally['filtered'][skip] += 1
            else:
                self._count(a, locus, w0, w1, tally)
        return self._conclude(locus, n, tally)

    def _count(self, a, locus, w0, w1, tally):
        pol, filtered = self.policy, tally['filtered']
        qs, qe, why = read_window(a, w0, w1)
        if why:
            filtered[why] += 1
            if why in ('not_spanning_the_window', 'window_end_in_a_gap') and a.mapping_quality >= pol.min_mapq and clipped_at(a, w0, w1, pol.clip_slack):
                tally['clipped'] += 1
            return
        seq = a.query_sequence[qs:qe + 1].upper()
        allele = 'ref' if seq == locus.ref_window else 'alt' if seq == locus.alt_window else 'other'
        if a.mapping_quality < pol.min_mapq:
            filtered['low_mapq'] += 1
            tally['low_mapq'][allele] += 1
            return
        quals = a.query_qualities
        if pol.min_baseq and (quals is None or min(quals[qs:qe + 1]) < pol.min_baseq):
            filtered['low_base_quality'] += 1
            tally['low_baseq'][allele] += 1
            return
        tally['informative_alignments'] += 1
        frags, got = tally['frags'], tally['frags'].get(a.query_name)
        if got is None:
            frags[a.query_name] = (allele, a.is_reverse, seq)
        elif got[0] != allele:
            frags[a.query_name] = ('conflict', a.is_reverse, '')
        # the same allele from the second mate of a fragment is the same evidence: counted once

    def _conclude(self, locus, overlapping, tally):
        pol, filtered, low_baseq, low_mapq, frags = self.policy, tally['filtered'], tally['low_baseq'], tally['low_mapq'], tally['frags']
        clipped = tally['clipped']
        n = Counter(v[0] for v in frags.values())
        other_seen = Counter(v[2] for v in frags.values() if v[0] == 'other')
        ref, alt, other, conflicts = n['ref'], n['alt'], n['other'], n['conflict']
        depth = ref + alt + other
        strand = Counter((v[0], 'reverse' if v[1] else 'forward') for v in frags.values() if v[0] in ('ref', 'alt'))
        gap = filtered['window_end_in_a_gap']
        mapq_low = sum(low_mapq.values())
        detail = dict(window=[locus.chrom, locus.start, locus.end], vtype=locus.vtype, overlapping_alignments=overlapping,
                      filtered={k: v for k, v in sorted(filtered.items()) if v}, informative_alignments=tally['informative_alignments'],
                      low_mapq_by_allele={k: v for k, v in sorted(low_mapq.items()) if v}, low_base_quality_by_allele={k: v for k, v in sorted(low_baseq.items()) if v},
                      conflicting_fragments=conflicts, soft_clipped_ending_at_the_locus=clipped, strand={f'{a}_{s}': v for (a, s), v in sorted(strand.items())},
                      other_alleles=dict(sorted(other_seen.items(), key=lambda kv: (-kv[1], kv[0]))[:5]))
        flags = []
        if alt >= pol.min_alt_reads and (strand[('alt', 'forward')] == 0 or strand[('alt', 'reverse')] == 0):
            flags.append('alternate_fragments_on_one_strand_only')
        out = lambda state, reasons=(): ReadSupport(locus.key, state, tuple(reasons), tuple(flags), depth, ref, alt, other, detail)
        shows = lambda k: k >= pol.min_alt_reads and k / depth >= pol.min_alt_fraction
        if shows(alt):
            if depth < pol.min_depth:
                flags.append('fewer_informative_fragments_than_min_depth')
            return out(ALTERNATE_PRESENT)
        if shows(other):
            if alt:
                flags.append('alternate_fragments_below_the_thresholds')
            if depth < pol.min_depth:
                flags.append('fewer_informative_fragments_than_min_depth')
            return out(OTHER_ALLELE_PRESENT)
        if depth < pol.min_depth:
            reasons = ['no_overlapping_alignments' if not overlapping else 'no_informative_fragments' if not depth else 'fewer_informative_fragments_than_min_depth']
            if overlapping and filtered['not_spanning_the_window'] * 2 > overlapping:
                reasons.append('most_alignments_do_not_span_the_window')
            if overlapping and sum(low_baseq.values()) * 2 > overlapping:
                reasons.append('most_alignments_fail_the_base_quality_policy')
            return out(INSUFFICIENT, reasons)
        if alt:
            return out(AMBIGUOUS, ['alternate_fragments_below_the_thresholds'])
        caveats = []
        if clipped >= pol.min_alt_reads and clipped / (depth + clipped) >= pol.max_clipped_fraction:
            caveats.append('soft_clipped_reads_end_at_the_locus')
        if gap >= pol.min_alt_reads and gap / (depth + gap) >= pol.min_alt_fraction:
            caveats.append('a_deletion_covers_the_locus_in_some_reads')
        if conflicts >= pol.min_alt_reads and conflicts / (depth + conflicts) >= pol.min_alt_fraction:
            caveats.append('overlapping_mates_disagree')
        if mapq_low and mapq_low / (mapq_low + tally['informative_alignments']) > pol.max_low_mapq_fraction:
            caveats.append('low_mapping_quality_reads_dominate')
        if low_mapq['alt'] >= pol.min_alt_reads:
            caveats.append('alternate_like_reads_below_the_mapping_quality_threshold')
        if low_baseq['alt'] >= pol.min_alt_reads:
            caveats.append('alternate_like_reads_below_the_base_quality_threshold')
        if caveats:
            return out(AMBIGUOUS, caveats)
        return out(REFERENCE_ONLY)
