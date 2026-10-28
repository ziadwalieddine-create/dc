"""Reference-confidence providers: what the input can say about positions where it holds no variant record.

A variant-only VCF says nothing about a position it does not list. A genotype of 'reference' needs positive evidence: a gVCF
reference block with adequate GQ and depth, or a callable-region mask for a caller that emits every variant inside it. Reads can also answer for a locus, but not through
this interface: vindex/reads.py counts the fragments that carry an indexed allele and the matcher writes that beside the status of the locus (it never changes it), so a position that
only reads could settle stays unresolved here, and the match says what the reads showed.
Each provider answers `status(chrom, start, end)` for a span with one of:

  reference_confident      the whole span has positive reference evidence (basis says which)
  reference_low_confidence reference evidence exists but fails the GQ/depth policy somewhere in the span
  variant_site             a variant record overlaps the span (the caller then looks at which alleles)
  no_call                  a record or block overlaps the span without a usable genotype, or the file holds a record there that the reader refused (its REF disagrees with the reference)
  deleted_by_upstream      a spanning-deletion allele covers the span
  structural_overlap       a called structural variant (deletion, duplication) overlaps the span
  unknown                  no evidence either way
"""
import bisect


class IntervalSet:
    """Closed 1-based integer intervals per contig, merged when overlapping or adjacent."""

    def __init__(self):
        self._raw = {}
        self._m = None

    def add(self, chrom, start, end):
        self._raw.setdefault(chrom, []).append((start, end))
        self._m = None

    def _merged(self):
        if self._m is None:
            m = {}
            for c, iv in self._raw.items():
                iv = sorted(iv)
                out = [list(iv[0])]
                for s, e in iv[1:]:
                    if s <= out[-1][1] + 1:
                        out[-1][1] = max(out[-1][1], e)
                    else:
                        out.append([s, e])
                m[c] = ([a for a, _ in out], [b for _, b in out])
            self._m = m
        return self._m

    def covers(self, chrom, start, end):
        starts, ends = self._merged().get(chrom, ([], []))
        i = bisect.bisect_right(starts, start) - 1
        return i >= 0 and ends[i] >= end

    def overlaps(self, chrom, start, end):
        starts, ends = self._merged().get(chrom, ([], []))
        i = bisect.bisect_right(starts, end) - 1
        return i >= 0 and ends[i] >= start

    def __len__(self):
        return sum(len(v) for v in self._raw.values())


class NoCallability:
    name = 'none'

    def status(self, chrom, start, end):
        return 'unknown', 'none'


class GvcfCallability:
    """Reference blocks and records of a read CallSet (a gVCF, or a joint VCF carrying 0/0 genotypes)."""
    name = 'gvcf'

    def __init__(self, callset):
        self.conf, self.low, self.nocall, self.site, self.spandel, self.sv = (IntervalSet() for _ in range(6))
        self._union = None            # confident and low-confidence blocks as one set, built on first use
        self.refused = IntervalSet()   # records the reader refused (the REF disagrees with the reference, an N in REF): something was called there that cannot be read
        for b in callset.refblocks:
            (self.conf if b.confident else self.low).add(b.chrom, b.start, b.end)
        for n in callset.nocalls:
            self.nocall.add(n.chrom, n.pos, n.end)
        for c in callset.small:
            if c.dosage > 0 or c.other_alt > 0:
                o = c.orig
                self.site.add(c.chrom, min(c.pos, o[1]), max(c.pos + len(c.ref) - 1, o[1] + len(o[2]) - 1))
        for rj in callset.rejected:
            if rj.canon:
                self.refused.add(rj.canon, rj.pos, rj.pos + len(rj.ref) - 1)
        for d in callset.spandels:
            self.spandel.add(d.chrom, d.pos, d.end)
        for v in callset.svs:
            if v.dosage and v.svtype in ('DEL', 'DUP', 'CNV', 'INV') and v.source == 'symbolic':
                self.sv.add(v.chrom, v.start, v.end)

    def status(self, chrom, start, end):
        if self.spandel.overlaps(chrom, start, end):
            return 'deleted_by_upstream', 'spanning_deletion_allele'
        if self.sv.overlaps(chrom, start, end):
            return 'structural_overlap', 'called_structural_variant'
        if self.nocall.overlaps(chrom, start, end):
            return 'no_call', 'no_call_record'
        if self.refused.overlaps(chrom, start, end):
            return 'no_call', 'record_refused_by_the_reader'      # reference is not verified where the file holds a record that could not be read
        if self.site.overlaps(chrom, start, end):
            return 'variant_site', 'variant_record'
        if self.conf.covers(chrom, start, end):
            return 'reference_confident', 'gvcf_reference_block'
        if self.low.overlaps(chrom, start, end) and (self.low.covers(chrom, start, end) or self._union_covers(chrom, start, end)):
            return 'reference_low_confidence', 'gvcf_reference_block_below_threshold'
        return 'unknown', 'no_block'

    def _union_covers(self, chrom, start, end):
        if self._union is None:        # once: a gVCF holds millions of blocks, and an indel next to a low-confidence block asks this for every such site
            both = IntervalSet()
            for src in (self.conf, self.low):
                for c, iv in src._raw.items():
                    for s, e in iv:
                        both.add(c, s, e)
            self._union = both
        return self._union.covers(chrom, start, end)


class CallableRegionCallability:
    """A callable-region mask (for example a benchmark's high-confidence BED). It supports 'reference_confident' only for spans inside
    the mask and only with the stated assumption that the caller reports every variant inside the mask; the basis records it.
    Records of the callset still take precedence (a variant, no-call or spanning deletion inside the mask is never reference)."""
    name = 'callable_region'

    def __init__(self, intervals, callset=None, label='callable_region'):
        self.mask = IntervalSet()
        for c, s, e in intervals:
            self.mask.add(c, s, e)
        self.inner = GvcfCallability(callset) if callset is not None else None
        self.label = label

    @classmethod
    def from_bed(cls, path, reference, callset=None, label='callable_region'):
        iv = []
        with open(path) as f:
            for line in f:
                if not line.strip() or line.startswith(('#', 'track', 'browser')):
                    continue
                c, s, e = line.split('\t')[:3]
                canon = reference.contig(c)
                if canon:
                    iv.append((canon, int(s) + 1, int(e)))   # BED is 0-based half-open
        return cls(iv, callset, label)

    def status(self, chrom, start, end):
        if self.inner is not None:
            st = self.inner.status(chrom, start, end)
            if st[0] in ('deleted_by_upstream', 'structural_overlap', 'no_call', 'variant_site'):
                return st
        if self.mask.covers(chrom, start, end):
            return 'reference_confident', f'{self.label}_assumes_all_variants_reported'
        return 'unknown', 'outside_callable_region'


class CompositeCallability:
    """First provider that gives a non-unknown answer, except that any 'variant_site/no_call/deleted/structural' answer wins."""
    name = 'composite'

    def __init__(self, providers):
        self.providers = providers

    def status(self, chrom, start, end):
        results = [p.status(chrom, start, end) for p in self.providers]
        for r in results:
            if r[0] in ('deleted_by_upstream', 'structural_overlap', 'no_call', 'variant_site'):
                return r
        for r in results:
            if r[0] in ('reference_confident', 'reference_low_confidence'):
                return r
        return results[0] if results else ('unknown', 'no_provider')
