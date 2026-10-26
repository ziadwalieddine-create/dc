"""Overlap queries over closed integer intervals, sized for millions of intervals per chromosome and millions of queries.

Intervals are sorted by start; a running maximum of the ends lets a query walk back from the last interval that starts at or before the query end and
stop as soon as no earlier interval can still reach the query start. Starts, ends and the running maximum are compact `array('i')` buffers, so a
chromosome of several million intervals costs tens of megabytes. Keep intervals of one kind (one stratification, one event set) in one store: a few
very long intervals mixed with many short ones lengthen the walk back.
"""
from array import array
from bisect import bisect_right


class IntervalStore:
    def __init__(self, starts, ends):
        """starts/ends: parallel sequences of ints (closed intervals, start <= end). `overlapping` returns positions in this input order."""
        n = len(starts)
        order = sorted(range(n), key=starts.__getitem__)
        self.n = n
        self.order = array('i', order)
        self.starts = array('i', (starts[i] for i in order))
        self.ends = array('i', (ends[i] for i in order))
        m, pm = -1, array('i')
        for e in self.ends:
            if e > m:
                m = e
            pm.append(m)
        self.pmax = pm

    def overlapping(self, start, end):
        """Input positions of every interval with interval.start <= end and interval.end >= start, in input order."""
        i = bisect_right(self.starts, end) - 1
        hits = []
        pm, ends, order = self.pmax, self.ends, self.order
        while i >= 0 and pm[i] >= start:
            if ends[i] >= start:
                hits.append(order[i])
            i -= 1
        hits.sort()
        return hits

    def any_overlap(self, start, end):
        i = bisect_right(self.starts, end) - 1
        pm, ends = self.pmax, self.ends
        while i >= 0 and pm[i] >= start:
            if ends[i] >= start:
                return True
            i -= 1
        return False
