"""dc_extended.variety — controlled randomness. Seeds logged for reproducibility.

Randomness lives ONLY in selection among critic-surviving candidates.
Gates, reads, and critic kills stay deterministic. A logged seed + fixed input
reproduces any output exactly.
"""
import hashlib
import math
import random


def derive_seed(*parts):
    """Deterministic seed from case + turn, so each turn varies but replays."""
    h = hashlib.sha256("|".join(str(p) for p in parts).encode())
    return int(h.hexdigest()[:8], 16)


def weighted_sample(items, weights, rng, temp=1.0):
    """Sample an index proportional to exp(w/temp).
    Falls back to argmax when rng is None (backward compatible)."""
    if not items:
        return None
    if len(items) == 1 or rng is None:
        return max(range(len(items)), key=lambda i: weights[i])
    mx = max(weights)
    exps = [math.exp((w - mx) / max(temp, 1e-6)) for w in weights]
    total = sum(exps)
    r = rng.random() * total
    acc = 0.0
    for i, e in enumerate(exps):
        acc += e
        if r <= acc:
            return i
    return len(items) - 1


def is_synonym(a, b):
    """Differ only by a trailing/leading adjective — the synonym-swap failure shape.
    Guards the two-draft rule: legitimate variants differ by phrasing of a
    PASSING frame; synonyms re-roll a killed one."""
    wa, wb = a.split(), b.split()
    if abs(len(wa) - len(wb)) != 1 or not wa or not wb:
        return False
    short, long = (wa, wb) if len(wa) < len(wb) else (wb, wa)
    return long[:len(short)] == short or long[1:] == short
