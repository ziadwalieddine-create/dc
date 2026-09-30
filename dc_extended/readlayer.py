"""dc_extended.readlayer — the read layer: base read, window reads, split
confidence, competing reads, grading.

Phase 1 upgrade: the bubble still gets one dominant read, but
- confidence is REAL: split when one bubble supports two kinds that lead to
  different actions (ledger A4 competing-reading procedure)
- the last 3 bubbles form a window; known sequences annotate the read
- a split read is PENDING: her confirmed reply grades which reading was right,
  making the read layer a scored classifier from data that arrives anyway
"""
import json
import os
import re
from pathlib import Path

from .lexicon import DRY_TOKENS, REFUSAL_TOKENS, WITHDRAW, CHARGE


def _has_word(b, words):
    """Word-boundary match: \bnow\b does not match 'know'; \bbet\b not 'between'."""
    return any(re.search(r"\b" + re.escape(w) + r"\b", b) for w in words)


def read_bubble(bubble, case):
    """One dominant read from dc_extended/read.md. Returns a read dict."""
    if not bubble:
        return {"kind": "none", "heat": None, "she_extended": False,
                "confidence": "high", "reversal": "any bubble from her"}

    b = bubble.strip().lower()
    words = set(b.split())

    # refusal first — a tease in the same bubble loses
    if _has_word(b, REFUSAL_TOKENS) and not any(
            k in b for k in ("jk", "just kidding")):
        return {"kind": "refuse", "heat": 0, "she_extended": False,
                "confidence": "high", "reversal": "a clear counter-offer"}

    # dry token
    if b in DRY_TOKENS or (len(words) <= 2 and not any(c in b for c in "?!")):
        return {"kind": "dry", "heat": 0, "she_extended": False,
                "confidence": "high", "reversal": "she adds words to the bit"}

    # logistics — asked about a day (before generic question)
    if any(w in b for w in ("when", "what day", "what time", "friday", "saturday", "sunday")):
        return {"kind": "logistics", "heat": 1, "she_extended": False,
                "confidence": "high", "reversal": "the window is already dead"}

    # question
    if "?" in b:
        live = any(r.get("status") == "live" for r in case.get("referents", []))
        return {"kind": "question", "heat": 1 if live else 0,
                "she_extended": live,
                "confidence": "high",
                "reversal": "the question is rhetorical / she refuses in next bubble"}

    # withdraw — busy/later/no counter
    if _has_word(b, ("busy", "later", "not sure", "maybe")) and not any(
            w in b for w in ("but", "how about", "what about", "should", "let's")):
        return {"kind": "withdraw", "heat": 0, "she_extended": False,
                "confidence": "high", "reversal": "she proposes a specific day"}

    # escalate — charge or challenge
    if any(w in b for w in ("come over", "tonight", "prove it", "show me")) or _has_word(b, ("now", "bet")):
        return {"kind": "escalate", "heat": 2, "she_extended": True,
                "confidence": "high", "reversal": "she says stop / not actually"}

    # tease — touched the live bit and added something
    live_refs = [r for r in case.get("referents", []) if r.get("status") == "live"]
    for r in live_refs:
        tokens = [w for w in r["text"].lower().split() if len(w) > 3]
        if tokens and any(t in b for t in tokens):
            return {"kind": "tease", "heat": 2, "she_extended": True,
                    "confidence": "high", "reversal": "the touch is accidental / she withdraws"}

    # offer — she named the move
    if any(w in b for w in ("we should", "let's", "drinks", "dinner", "hang out")):
        return {"kind": "offer", "heat": 1, "she_extended": True,
                "confidence": "high", "reversal": "she hedges when he plays it"}

    # topic — live bit exists but this is not it
    if live_refs:
        return {"kind": "topic", "heat": 0, "she_extended": False,
                "confidence": "high", "reversal": "she returns to the bit herself"}

    # invest — new specific, no bit
    if len(words) >= 4:
        return {"kind": "invest", "heat": 1, "she_extended": False,
                "confidence": "high", "reversal": "the specific is a complaint"}

    return {"kind": "invest", "heat": 1, "she_extended": False,
            "confidence": "high", "reversal": "no reply"}

WINDOW_SIZE = 3

# A split = candidate kinds from ONE bubble that select DIFFERENT actions.
SPLIT_PAIRS = [
    ({"escalate", "tease", "offer"}, {"withdraw"}),
    ({"question"}, {"withdraw"}),
]

WINDOW_PATTERNS = [
    (("dry", "dry", "question"), "testing_investment"),
    (("tease", "dry", "withdraw"), "push_pull"),
    (("question", "dry", "dry"), "fading"),
]

# Her confirmed reply grades the split (affordance-as-probe).
ENGAGED = {"extends", "counters", "mirrors_length"}
DEFLECTED = {"dry", "silence", "topic_change"}
NO_COUNT = {"literal", "refuses", "answers"}

HEAT_RANK = {"escalate": 2, "tease": 2, "offer": 1, "question": 1,
             "withdraw": 0, "dry": 0, "topic": 0, "invest": 1}

_STATS_PATH = Path(__file__).resolve().parents[1] / "read_stats.json"


def push_window(case, bubble):
    """Maintain the last 3 of her bubbles. Consecutive duplicates collapse."""
    if not bubble:
        return case
    win = case.setdefault("her_recent_bubbles", [])
    if win and win[-1] == bubble:
        return case
    win.append(bubble)
    del win[:-WINDOW_SIZE]
    return case


def _has(b, words):
    return any(w in b for w in words)


def _split(primary, secondary, bubble):
    return {"kind": primary, "heat": HEAT_RANK.get(primary, 1),
            "she_extended": primary in ("escalate", "tease"),
            "confidence": "split", "alternates": [primary, secondary],
            "signal": bubble,
            "reversal": "her next reply confirms %s over %s" % (primary, secondary)}


def read_bubble_v2(bubble, case):
    """One dominant read + real confidence + window annotation."""
    override = case.get("last_read_override") if case else None
    if override and override.get("bubble") == bubble:
        if case is not None:
            case.pop("last_read_override", None)
        r = dict(override["read"])
        r["confidence"] = "high"
        return r

    read = read_bubble(bubble, case)
    b = (bubble or "").strip().lower()
    kind = read["kind"]

    # competing readings: engagement markers and withdrawal in one breath
    withdraw = _has_word(b, WITHDRAW)
    charge = _has_word(b, CHARGE)
    for hot, cold in SPLIT_PAIRS:
        if kind in hot and withdraw:
            return _split(kind, "withdraw", bubble)
        if kind == "withdraw" and charge:
            return _split("escalate", "withdraw", bubble)
        if kind in cold and charge:
            return _split("escalate", kind, bubble)

    # touched the live bit but the base read went cold: tease competes with cold
    if kind in ("withdraw", "dry"):
        for r in (case.get("referents", []) if case else []):
            if r.get("status") != "live":
                continue
            toks = [w for w in r["text"].lower().split() if len(w) > 3]
            if toks and any(tk in b for tk in toks):
                return _split("tease", kind, bubble)

    # window patterns annotate low-energy reads (matched on read kinds)
    if case is not None:
        win = (case.get("her_recent_bubbles") or [])[-WINDOW_SIZE:]
        seq = tuple(read_bubble(w, case)["kind"] for w in win)
        for pattern, signal in WINDOW_PATTERNS:
            if seq == pattern:
                read = dict(read)
                read["window_signal"] = signal
                if signal == "push_pull":
                    read["heat"] = max(read.get("heat") or 0, 1)
                return read
    return read


def needs_deep(read):
    """Compute contract: deep mode fires only on split reads (and integrity stops)."""
    return read.get("confidence") == "split"


def resolve_split(case, chosen):
    """The model's judgment, written as state. chosen must be one of the alternates."""
    pending = case.get("pending_read")
    if not pending:
        raise ValueError("no pending split read")
    if chosen not in pending["alternates"]:
        raise ValueError("chosen kind %r is not an alternate %r" % (chosen, pending["alternates"]))
    read = {"kind": chosen, "heat": HEAT_RANK.get(chosen, 1),
            "she_extended": chosen in ("escalate", "tease"),
            "confidence": "high", "signal": pending["signal"],
            "reversal": "resolved by operator/model from split"}
    case["last_read_override"] = {"bubble": pending["signal"], "read": read}
    case["last_split"] = pending          # archived for grading by her reply
    case.pop("pending_read", None)
    case["last_read"] = read
    return read


def grade_pending_read(case, response_class, stats_path=None):
    """Her confirmed reply grades the pending split. Returns the winning kind or None.

    extends/answers -> the engaged reading was right; dry/silence -> the
    deflection reading was right; literal/refuses -> no count.
    stats_path is REQUIRED: production passes the real path, tests pass a temp.
    """
    if stats_path is None:
        raise ValueError("stats_path required: pass a production path or a temp path for tests")
    pending = case.pop("last_split", None) or case.pop("pending_read", None)
    if not pending:
        return None
    primary = pending["primary"]
    alternates = pending["alternates"]
    if response_class in ENGAGED:
        winner = primary
    elif response_class in DEFLECTED:
        winner = [a for a in alternates if a != primary][0]
    else:
        return None  # NO_COUNT
    stats_file = Path(stats_path) if stats_path else _STATS_PATH
    try:
        stats = json.loads(stats_file.read_text())
    except Exception:
        stats = {}
    row = stats.setdefault(primary, {"n": 0, "correct": 0})
    row["n"] += 1
    if winner == primary:
        row["correct"] += 1
    stats_file.write_text(json.dumps(stats, indent=2) + "\n")
    return winner
