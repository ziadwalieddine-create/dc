"""dc_extended.recovery — one deviation function, forced resolution, reciprocity.

Merges: operator rejection transaction, transitions table, incident intake.
Mechanism stays UNKNOWN/null unless named — never invent a tag.
"""
from datetime import datetime

RESPONSE_DELTA = {
    "extends": 1, "answers": 1, "counters": 1, "mirrors_length": 1,
    "literal": 0, "topic_change": 0, "dry": -1, "silence": -1, "refuses": -2,
}

# Event -> state mutation. One table, one place that writes recovery state.
EVENT_MUTATIONS = {
    "joke_misunderstood": {"close_referent": True, "explain": False,
                            "next": "Trade if nothing else live, else Play on a different referent"},
    "flirt_ignored_once": {"amp_cap": 0, "next": "Play, no tease"},
    "flirt_ignored_twice": {"close_referent": True},
    "literal_answer":     {"store_fact": True, "next": "Do not repeat the tease she stepped over"},
    "dry_token":          {"heat": 0},
    "topic_change":       {"keep_old_live": True},
    "logistics_unanswered": {"field_stays_null": True},
    "too_forward":        {"amp_cap_delta": -1},
    "too_boring":         {"next_must_pass": "missed_heat"},
    "operator_rejected":  {"store_last_candidate": True, "ban_family": True},
    "misread":            {"delete_inferences": True, "next": "Rerun from facts; no apology unless sent"},
    "gate_stop":          {},  # mechanism = the gate that fired (established, not UNKNOWN)
}

KNOWN_EVENTS = set(EVENT_MUTATIONS)


def record_deviation(case, event, mechanism=None, family=None, gate=None,
                     referent_id=None, ts=None, context=None):
    """Something went wrong -> update state so it cannot recur.

    mechanism: canonical gate/tag name, or None (UNKNOWN — legal, never invented).
    family:    object family to ban (operator_rejected only).
    gate:      named defect when the operator states one.
    Returns the deviation row.
    """
    if event not in KNOWN_EVENTS:
        raise ValueError("unknown deviation event: %r" % event)
    row = {
        "event": event,
        "mechanism": mechanism,   # UNKNOWN stays null
        "family": family,
        "gate": gate,
        "context": context or {},
        "ts": ts or datetime.now().isoformat(timespec="seconds"),
    }
    case.setdefault("deviations", []).append(row)

    mut = EVENT_MUTATIONS[event]
    if mut.get("close_referent") and referent_id:
        for r in case.get("referents", []):
            if r.get("id") == referent_id:
                r["status"] = "closed"
    if "amp_cap" in mut:
        case["amp_cap"] = mut["amp_cap"]
    if "amp_cap_delta" in mut:
        cur = case.get("amp_cap")
        base_amp = case.get("amplitude", 1) if cur is None else cur
        case["amp_cap"] = base_amp + mut["amp_cap_delta"]
    if mut.get("ban_family") and family:
        case.setdefault("banned", []).append(
            {"string": (case.get("last_send") or {}).get("text", "") if isinstance(case.get("last_send"), dict) else "",
             "family": family, "gate": gate})
    row["next"] = mut.get("next")
    return row


# Stops that must never degrade into a line. Everything else may.
INTEGRITY_STOPS = {"name_clash", "refusal", "closed", "null_claim",
                   "stated_boundary", "wrong_person"}


def forced_resolution(case, why):
    """After a degraded hold, the nearest legal line. None if illegal.

    Integrity stops (wrong woman, refusal, closed claim, null claim) never
    resurrect. Degraded holds (unanswered, dead clock) downgrade: amplitude 0
    Play on a live referent, else Trade on a FACT, else None.
    """
    if why in INTEGRITY_STOPS:
        return None
    lives = [r for r in case.get("referents", []) if r.get("status") == "live"]
    refused = bool(case.get("refused_object"))
    usable = [] if (refused and len(lives) <= 1) else lives
    if usable:
        return {"policy": "Play", "amplitude": 0, "via": "forced:%s" % why}
    facts = [f for f in case.get("facts", []) if f.get("class") == "FACT"]         if case.get("facts") else [t for t in case.get("unique_facts", [])]
    if facts and case.get("claim", {}).get("type") in ("bit", None) and not refused:
        return {"policy": "Trade", "amplitude": 0, "via": "forced:%s" % why}
    return None


def update_reciprocity(case, response_class):
    """Per-case leading indicator. Thread-level signal vs the bubble-level read.

    Receipt discipline: her reply counts only after a confirmed send.
    Two dry rows in a row cap amplitude at 0. Integer, clamped, no decimals.
    """
    if response_class not in RESPONSE_DELTA:
        raise ValueError("unknown response class: %r" % response_class)
    last_send = case.get("last_send")
    if not (last_send and last_send.get("confirmed")):
        raise ValueError("no receipt: outcomes require a confirmed send")
    case["reciprocity"] = max(-3, min(3, case.get("reciprocity", 0) + RESPONSE_DELTA[response_class]))
    outcomes = case.setdefault("outcomes", [])
    outcomes.append({"response": response_class,
                     "sent": last_send.get("text", "")[:80],
                     "ts": datetime.now().isoformat(timespec="seconds")})
    if [o["response"] for o in outcomes[-2:]] == ["dry", "dry"]:
        case["amp_cap"] = 0
    return case["reciprocity"]
