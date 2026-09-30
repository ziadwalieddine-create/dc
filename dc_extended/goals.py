"""dc_extended.goals — endpoint lock (kernel O*) and receptivity triage.

Endpoint: the operator's goal this turn wins (kernel O*, dc rule 14:
substituted_end). Never ladder-inferred. Stored with source, like claim.
Triage: derived from reciprocity + recent outcomes. Computed, never stored.
"""
BASE_SIM = {"extends": 4, "answers": 3, "mirrors_length": 2, "counters": 2,
            "dry": 1, "literal": 1, "topic_change": 0, "silence": -1, "refuses": -3}

ENDPOINT_DELTA = {
    "off_app":   {"answers": 1},                            # her answer IS the move
    "date_lock": {"counters": 2, "answers": 1, "dry": -1},  # a counter is a kept plan
    "banter":    {"extends": 1, "answers": -1},             # charge over logistics
    "close":     {},
}

_ENDPOINT_PHRASES = [
    (("get to know", "just banter", "play the bit"), "banter"),
    (("off app", "to ig", "her number", "instagram"), "off_app"),
    (("lock the date", "book it", "set the date", "the date"), "date_lock"),
    (("close it", "wrap it up"), "close"),
]


def lock_endpoint(operator_text, stored):
    """His words this turn set the endpoint; stored survives a silent turn."""
    if operator_text:
        t = operator_text.lower()
        for phrases, ep in _ENDPOINT_PHRASES:
            if any(p in t for p in phrases):
                return {"type": ep, "source": "operator_this_turn"}
        if stored and stored.get("type"):
            return {"type": stored["type"], "source": "stored"}
        return {"type": None, "source": None}
    if stored and stored.get("type"):
        return {"type": stored["type"], "source": "stored"}
    return {"type": None, "source": None}


def sim_rank(endpoint_type):
    """Endpoint-locked ranking. The same reply is worth different amounts
    against different endpoints (kernel: trajectory, not local cleverness)."""
    rank = dict(BASE_SIM)
    for cls, d in ENDPOINT_DELTA.get(endpoint_type or "", {}).items():
        rank[cls] = rank.get(cls, 0) + d
    return rank


def triage(case):
    """Receptive / neutral / cooling / unreceptive — the thread's disposition.

    Computed per turn from reciprocity + recent outcomes. Never stored.
    cooling = two flat outcomes on a neutral thread: his one-move budget is spent.
    """
    if case.get("closed"):
        return "closed"
    rec = case.get("reciprocity", 0)
    recent = [o.get("response") for o in case.get("outcomes", [])][-3:]
    if any(r == "refuses" for r in recent):
        return "unreceptive"
    if rec >= 2:
        return "receptive"
    if rec <= -2 or (len(recent) >= 2 and all(r in ("dry", "silence") for r in recent[-2:])):
        return "cooling"
    return "neutral"


def queue_weight(case):
    """Triage as queue priority: receptive first, cooling last."""
    return {"receptive": 0, "neutral": 1, "cooling": 2,
            "unreceptive": 3, "closed": 4}.get(triage(case), 1)
