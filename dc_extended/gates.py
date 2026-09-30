"""dc_extended.gates — the twelve. Code only; no model judgment."""
from datetime import datetime, timedelta

NAMES = ["wrong_person","wrong_platform","cross_thread","dead_clock","restack",
         "dual_plan","stale_goal","unconfirmed_send","no_write","goal_smuggle",
         "stated_boundary","paste_without_send"]

WEEKDAYS = ["monday","tuesday","wednesday","thursday","friday","saturday","sunday"]

def _day_index(word):
    return WEEKDAYS.index(word)

def day_dead(word, now, cutoff_hour=18):
    """A named weekday is dead if it is past, or it is today past the cutoff."""
    t = _day_index(word)
    today = now.weekday()
    if t < today:
        return True
    if t == today and now.hour >= cutoff_hour:
        return True
    return False

GATE_LAYER_GATES = ["wrong_person", "wrong_platform", "cross_thread", "dead_clock",
                    "restack", "dual_plan", "stale_goal", "unconfirmed_send",
                    "goal_smuggle", "stated_boundary"]
SWITCHBOARD_STOPS = ["paste_without_send", "refusal_close", "staleness"]
ALL_GATES = GATE_LAYER_GATES + SWITCHBOARD_STOPS


def run_gates(case, draft, now, claim=None, all_cases=None, log=None):
    """First matching row wins. Later rows do not generate.

    log: optional list; receives {"gate", "result", "detail"} for EVERY gate,
    observed not fabricated. Integrity stops never degrade.
    """
    hits = []
    fired = set()

    def hit(gate, detail):
        hits.append((gate, detail))
        fired.add(gate)
        record(gate, "FAIL", detail)

    def record(gate, result, detail=""):
        if log is not None:
            log.append({"gate": gate, "result": result, "detail": detail})

    text = (draft or "").lower()

    # stale_goal — closed threads generate nothing
    if case.get("closed"):
        hit("stale_goal", "thread closed: %s" % case.get("close_cause", "unknown"))
    else:
        record("stale_goal", "PASS")

    # dead_clock — weekday tokens in the draft
    dc_detail = ""
    for w in WEEKDAYS:
        if w in text and day_dead(w, now):
            dc_detail = "%s is dead at %s" % (w, now.isoformat(timespec="seconds"))
            break
    # claim-level: the plan's own day word is dead (window derived this turn)
    if not dc_detail:
        op = case.get("open_plan") or {}
        if (claim or {}).get("type") == "plan" and op.get("when") and day_dead(op["when"], now):
            dc_detail = "%s is dead at %s" % (op["when"], now.isoformat(timespec="seconds"))
    # never propose the adjacent day as a synonym for a missed window
    if not dc_detail:
        for mw in case.get("missed_windows", []):
            adj = WEEKDAYS[(_day_index(mw["day"]) + 1) % 7]
            if adj in text:
                dc_detail = "never propose %s as synonym for missed %s" % (adj, mw["day"])
    if dc_detail:
        hit("dead_clock", dc_detail)
    else:
        record("dead_clock", "PASS")

    fam = (claim or {}).get("family")
    # restack — refused families stay dead
    if fam and fam in case.get("banned_families", []):
        hit("restack", "family refused earlier: %s" % fam)
    else:
        record("restack", "PASS")

    # dual_plan
    if case.get("open_plan") and (claim or {}).get("type") == "plan":
        hit("dual_plan", "open plan exists: %s" % case["open_plan"].get("what"))
    else:
        record("dual_plan", "PASS")

    # goal_smuggle — turn-locked endpoint first, stored goal fallback
    ep = (case.get("endpoint") or {}).get("type")
    want_type = ep or ((case.get("goal") or {}).get("type"))
    if (claim or {}).get("type") == "plan" and want_type in ("get_to_know", "no_dates"):
        hit("goal_smuggle", "plan claim vs operator want: %s" % want_type)
    else:
        record("goal_smuggle", "PASS")

    # unconfirmed_send — any second bid on an unconfirmed send
    ls = case.get("last_send")
    if ls and not ls.get("confirmed") and not case.get("newer_bubble"):
        hit("unconfirmed_send", "last send unconfirmed: %s" % str(ls.get("text", ""))[:40])
    elif ls and not ls.get("confirmed") and (claim or {}).get("builds_on_last"):
        hit("unconfirmed_send", "builds on unconfirmed send: %s" % str(ls.get("text", ""))[:40])
    else:
        record("unconfirmed_send", "PASS")

    # wrong_platform
    ch = (claim or {}).get("channel")
    if ch and ch != case.get("channel"):
        hit("wrong_platform", "thread lives on %s, draft targets %s" % (case.get("channel"), ch))
    else:
        record("wrong_platform", "PASS")

    # stated_boundary
    sb_detail = ""
    for b in case.get("boundaries", []):
        if b.lower() in text:
            sb_detail = "draft touches boundary topic: %s" % b
            break
    if sb_detail:
        hit("stated_boundary", sb_detail)
    else:
        record("stated_boundary", "PASS")

    # wrong_person + cross_thread via other cases' names and facts
    wp_detail = ""
    xt_detail = ""
    if all_cases:
        for other in all_cases:
            if other["id"] == case["id"]:
                continue
            oname = (other.get("name") or "").lower()
            if oname and oname in text and not wp_detail:
                wp_detail = "draft addresses %s" % other.get("name")
            for fact in other.get("unique_facts", []):
                if fact.lower() in text and not xt_detail:
                    xt_detail = "uses %s's fact: %s" % (other.get("name"), fact)
    if wp_detail:
        hit("wrong_person", wp_detail)
    else:
        record("wrong_person", "PASS")
    if xt_detail:
        hit("cross_thread", xt_detail)
    else:
        record("cross_thread", "PASS")

    return hits
