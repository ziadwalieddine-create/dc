"""dc_extended.intake — raw trigger-text parsing and case binding.

The triggering message IS the context. Do not wait for clean fields.
Ports: open_turn, bind_case (ckososher engine.py), claim precedence by scan order.
"""
import re
from .switchboard import classify

_HER_MARKS = ("her last:", "she said:", "she said ", "she:", "her:")
_HIS_MARKS = ("me:", "him:", "i said:", "i:")

# Last-match wins, so "don't book thursday" locks bit, not plan.
_PHRASES = (
    "don't book",
    "do not book",
    "get to know her",
    "play the bit",
    "just the bit",
    "leave it",
    "don't contact",
    "do not contact",
    "stop",
    "instagram",
    "to ig",
    "her number",
    "move to sms",
    "lock thursday",
    "lock it",
    "set the date",
    "the plan",
    "book",
)


def trigger_has_her(text):
    low = (text or "").lower()
    return any(mark in low for mark in _HER_MARKS)


def _instruction_tail(text):
    """Split a her-line into (bubble, instruction-note) by last phrase match.

    Containment rule: a match whose span sits inside a longer match's span
    is dropped, so "don't book thursday" locks bit, not plan.
    """
    low = text.lower()
    matches = []
    for phrase in _PHRASES:
        start = 0
        while True:
            i = low.find(phrase, start)
            if i < 0:
                break
            matches.append((i, phrase))
            start = i + 1
    if not matches:
        return text, None
    kept = []
    for i, p in matches:
        contained = any(
            (j, q) != (i, p) and j <= i and i + len(p) <= j + len(q)
            for j, q in matches
        )
        if not contained:
            kept.append((i, p))
    if not kept:
        kept = matches
    at, hit = max(kept, key=lambda m: m[0])
    return text[:at].strip(" .,-"), text[at:].strip()


def lock_claim_intake(operator_text, stored):
    """His words this turn replace the stored claim. Ordered precedence."""
    if operator_text:
        t = operator_text.lower()
        base = classify(operator_text) or {}
        fam = base.get("family")
        if any(re.search(r"\b" + re.escape(w) + r"\b", t) for w in ("stop", "block")) or "leave it" in t or "do not contact" in t or "don't contact" in t:
            return {"type": "close", "source": "operator_this_turn", "text": operator_text,
                    **({"family": fam} if fam else {})}
        if any(w in t for w in ("don't book", "do not book", "get to know", "play the bit", "just the bit")):
            return {"type": "bit", "source": "operator_this_turn", "text": operator_text,
                    **({"family": fam} if fam else {})}
        if any(w in t for w in ("instagram", "to ig", "her number", "move to sms")):
            return {"type": "channel", "source": "operator_this_turn", "text": operator_text,
                    **({"family": fam, "channel": "ig"} if fam else {"channel": "ig"})}
        if any(w in t for w in ("book", "lock", "thursday", "the plan", "set the date")):
            return {"type": "plan", "source": "operator_this_turn", "text": operator_text,
                    **({"family": fam} if fam else {})}
        if stored and stored.get("type"):
            return {"type": stored["type"], "source": "stored", "text": stored.get("text")}
        claim = classify(operator_text)
        if claim and claim.get("type") != "audit":
            claim["source"] = "operator_this_turn"
            return claim
        return {"type": None, "source": None, "text": None}
    if stored and stored.get("type"):
        return {"type": stored["type"], "source": "stored", "text": stored.get("text")}
    return {"type": None, "source": None, "text": None}


def open_turn(text, case=None):
    """Parse the raw triggering message into her bubble + his claim.

    A `she:`/`her last:` line is her bubble; his instruction in the same
    message sets the claim. Never invent a bubble that was not marked.
    """
    if case is None:
        case = {"her_last_bubble": None, "claim": {"type": None, "source": None, "text": None},
                "newer_bubble": False}
    raw = (text or "").strip()
    if not raw:
        return case
    her_lines, notes = [], []
    for line in raw.splitlines():
        line = line.strip()
        if not line:
            continue
        low = line.lower()
        found = None
        for mark in _HER_MARKS:
            if low.startswith(mark):
                found = line[len(mark):].strip()
                break
        if found is not None:
            bubble, note_txt = _instruction_tail(found)
            if bubble:
                her_lines.append(bubble)
            if note_txt:
                notes.append(note_txt)
            continue
        if any(low.startswith(mark) for mark in _HIS_MARKS):
            notes.append(line.split(":", 1)[-1].strip())
            continue
        notes.append(line)
    if her_lines:
        case["her_last_bubble"] = her_lines[-1]
        case["newer_bubble"] = True
    for note_txt in notes:
        if note_txt:
            case["claim"] = lock_claim_intake(note_txt, case.get("claim"))
    return case


def bind_case(case, open_ids):
    """Two open files for the same woman are a mistake, not two cases."""
    case_id = case.get("id")
    if not case_id or list(open_ids).count(case_id) != 1:
        raise ValueError("case binding: id missing or duplicated: %r" % case_id)
    return case
