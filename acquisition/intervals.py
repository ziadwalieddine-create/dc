"""Exact stabbing queries over half-open intervals, used by every genomic-overlap test in compat_tests.py."""
import bisect


def build_index(iv):
    """Index for half-open intervals [s, e): sorted starts plus the prefix maximum of ends."""
    iv = sorted(iv)
    starts = [s for s, e in iv]
    pm, m = [], -1
    for s_, e_ in iv:
        m = max(m, e_)
        pm.append(m)
    return starts, pm


def covered(idx, pos0):
    """True if the 0-based position pos0 lies inside at least one indexed interval."""
    starts, pm = idx
    i = bisect.bisect_right(starts, pos0) - 1
    return i >= 0 and pm[i] > pos0
