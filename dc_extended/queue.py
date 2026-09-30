"""dc_extended.queue — windows, decay, the nightly order of operations."""
from datetime import timedelta
from .state import parse_ts
from .goals import queue_weight

STAGE_DECAY_HOURS = {"matched":72,"responded":96,"off_app":144,
                     "plan_pending":48,"met":48,"followup":72}
PLAN_CONFIRM_HOURS = 48
STALE_HOURS = 48


def derive_stage(case):
    """Derived from live fields, never stored (the stage field was abolished).
    matched      — no confirmed send yet
    plan_pending — an open plan exists
    off_app      — conversation lives on ig/sms
    responded    — has a confirmed send, no open plan, app channel
    """
    if case.get("closed"):
        return "closed"
    if (case.get("open_plan") or {}).get("what"):
        return "plan_pending"
    if case.get("channel") in ("ig", "sms"):
        return "off_app"
    if (case.get("last_send") or {}).get("confirmed"):
        return "responded"
    return "matched"

def _confirmed_dm(case):
    ls = case.get("last_send")
    return bool(ls and ls.get("confirmed") and ls.get("family") in ("dm","callback","first_dm"))

def build_queue(cases, now):
    weights = {c["id"]: queue_weight(c) for c in cases}
    items = []
    for c in cases:
        if derive_stage(c) == "closed":
            continue
        updated = parse_ts(c["updated"])
        cid = c["id"]


        p = c.get("open_plan")
        if p and p.get("proposed_ts") and not p.get("confirmed"):
            dl = parse_ts(p["proposed_ts"]) + timedelta(hours=PLAN_CONFIRM_HOURS)
            items.append({"id":cid,"action":"confirm_or_repair_plan",
                          "deadline":dl.isoformat(timespec="seconds"),
                          "note":"never propose the adjacent day as a synonym"})

        dl2 = updated + timedelta(hours=STAGE_DECAY_HOURS.get(derive_stage(c),96))
        if now > dl2:
            items.append({"id":cid,"action":"archive",
                          "deadline":dl2.isoformat(timespec="seconds"),
                          "note":"stage decay exceeded; close with a cause"})

        if (now - updated).total_seconds() > STALE_HOURS*3600:
            items.append({"id":cid,"action":"reverify",
                          "deadline":(updated+timedelta(hours=STALE_HOURS)).isoformat(timespec="seconds"),
                          "note":"state stale >48h; re-verify before any line"})
    items.sort(key=lambda i: (i["deadline"], weights.get(i["id"], 1)))
    return items
