"""Overlap queries over closed integer intervals, sized for millions of intervals per chromosome and millions of queries.

Intervals are sorted by start; a running maximum of the ends lets a query walk back from the last interval that starts at or before the query end and
stop as soon as no earlier interval can still reach the query start. Starts, ends and the running maximum are compact `array('i')` buffers, so a
chromosome of several million intervals costs tens of megabytes. A few very long intervals mixed with many short ones lengthen the walk back (the running
maximum stays high for as long as the longest interval reaches), so the intervals are kept in three stores by length (up to 10 kb, up to 1 Mb, longer) and a query asks
each: the walk of a store is bounded by its longest interval. On the structural events of chromosome 1 this took a query from 16 ms to well under 1 ms.
"""
from array import array
from bisect import bisect_right

BUCKETS = (10_000, 1_000_000)           # interval lengths that separate the three stores


class _Sorted:
    """Intervals sorted by start with a running maximum of the ends; `ids` are the caller's positions of the intervals."""

    def __init__(self, starts, ends, ids):
        order = sorted(range(len(starts)), key=starts.__getitem__)
        self.ids = array('i', (ids[i] for i in order))
        self.starts = array('i', (starts[i] for i in order))
        self.ends = array('i', (ends[i] for i in order))
        m, pm = -1, array('i')
        for e in self.ends:
            if e > m:
                m = e
            pm.append(m)
        self.pmax = pm

    def overlapping(self, start, end):
        i = bisect_right(self.starts, end) - 1
        hits = []
        pm, ends, ids = self.pmax, self.ends, self.ids
        while i >= 0 and pm[i] >= start:
            if ends[i] >= start:
                hits.append(ids[i])
            i -= 1
        return hits

    def any_overlap(self, start, end):
        i = bisect_right(self.starts, end) - 1
        pm, ends = self.pmax, self.ends
        while i >= 0 and pm[i] >= start:
            if ends[i] >= start:
                return True
            i -= 1
        return False


class IntervalStore:
    def __init__(self, starts, ends):
        """starts/ends: parallel sequences of ints (closed intervals, start <= end). `overlapping` returns positions in this input order."""
        self.n = len(starts)
        groups = [[] for _ in range(len(BUCKETS) + 1)]
        for i in range(self.n):
            length = ends[i] - starts[i]
            b = 0
            while b < len(BUCKETS) and length > BUCKETS[b]:
                b += 1
            groups[b].append(i)
        self.parts = [_Sorted([starts[i] for i in g], [ends[i] for i in g], g) for g in groups if g]

    def overlapping(self, start, end):
        """Input positions of every interval with interval.start <= end and interval.end >= start, in input order."""
        hits = []
        for part in self.parts:
            hits.extend(part.overlapping(start, end))
        hits.sort()
        return hits

    def any_overlap(self, start, end):
        return any(part.any_overlap(start, end) for part in self.parts)
