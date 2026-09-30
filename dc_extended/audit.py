"""dc_extended.audit — the system audits itself from its own emissions.

- honesty: UNKNOWN-mechanism ratio over the deviation ledger (is it learning
  faster than it fails?)
- charter query: families with >= 2 logged incidents propose new gates
  (dc's admission rule, executable)
- gate hit rates: which gates fire (never-firing gates are deletion candidates)
- deviations with full context convert to replay fixtures mechanically
"""
import json
import os
import glob
from collections import Counter
from .queue import build_queue
from .state import load_all


def collect_deviations(cases_dir):
    rows = []
    for p in sorted(glob.glob(os.path.join(cases_dir, "*.json"))):
        base_name = os.path.basename(p)
        if base_name.startswith("_"):
            continue
        try:
            case = json.load(open(p))
        except Exception:
            continue
        for row in case.get("deviations", []):
            r = dict(row)
            r["case_id"] = case.get("id") or os.path.basename(p)[:-5]
            rows.append(r)
    return rows


def honesty(rows):
    """The epistemic health number. Climbing ratio = failing faster than learning."""
    if not rows:
        return {"unknown": 0, "total": 0, "ratio": 0.0}
    unknown = sum(1 for r in rows if not r.get("mechanism"))
    return {"unknown": unknown, "total": len(rows),
            "ratio": round(unknown / len(rows), 3)}


def charter_query(rows, threshold=2, banned_families=frozenset()):
    """dc's admission rule as a query: two logged incidents of a new family.
    Families already banned (gated) do not re-propose."""
    fams = Counter(r["family"] for r in rows if r.get("family"))
    return sorted(f for f, n in fams.items() if n >= threshold and f not in banned_families)


def all_banned_families(cases_dir):
    fams = set()
    for p in glob.glob(os.path.join(cases_dir, "*.json")):
        bn = os.path.basename(p)
        if bn.startswith("_"):
            continue
        try:
            case = json.load(open(p))
        except Exception:
            continue
        for b in case.get("banned", []) or []:
            if b.get("family"):
                fams.add(b["family"])
    return fams


def gate_hit_rates(rows):
    hits = Counter(r["mechanism"] for r in rows
                   if r.get("event") == "gate_stop" and r.get("mechanism"))
    return dict(hits)


def nightly_checkin(cases_dir, now=None):
    """The nightly queue plus one line of truth: honesty + proposed gates."""
    from datetime import datetime
    now = now or datetime.now()
    rows = collect_deviations(cases_dir)
    return {"queue": build_queue(load_all(cases_dir), now),
            "honesty": honesty(rows),
            "proposed_gates": charter_query(rows, banned_families=all_banned_families(cases_dir)),
            "gate_hit_rates": gate_hit_rates(rows)}


def deviation_to_fixture(row):
    """Mechanical fixture from a gate_stop deviation with full context.
    Rows without context stay in the ledger but cannot replay — said so, not guessed."""
    if row.get("event") != "gate_stop" or not row.get("mechanism"):
        return None
    ctx = row.get("context") or {}
    if "bubble" not in ctx:
        return None
    return {"name": "dev:%s:%s" % (row.get("case_id"), row["mechanism"]),
            "mechanism": row["mechanism"],
            "context": ctx}


def fixtures_from(rows):
    return [f for f in (deviation_to_fixture(r) for r in rows) if f]


def replay_fixture(fx, cases_dir, corpus_path, now):
    """Replay one deviation fixture: the named gate must fire again."""
    from .switchboard import run_turn
    ctx = fx["context"]
    cid = "_devcase"
    case = {"id": cid, "name": "Dev", "stage": "chat", "channel": ctx.get("channel", "hinge"),
            "goal": {"type": "date"}, "updated": "2026-09-19T10:00:00",
            "banned_families": ctx.get("banned_families", []), "boundaries": [],
            "unique_facts": ctx.get("unique_facts", []), "missed_windows": [],
            "referents": ctx.get("referents", []), "frame_limits": [], "banned": [],
            "her_recent_bubbles": [], "closed": False}
    case["banned_families"] = ctx.get("banned_families", case["banned_families"])
    if ctx.get("open_plan") is not None:
        case["open_plan"] = ctx["open_plan"]
    with open(os.path.join(cases_dir, cid + ".json"), "w") as f:
        json.dump(case, f)
    op = ctx.get("operator_text") or {"plan": "thursday drinks", "bit": "mean",
                                      "answer": "still waiting"}.get(ctx.get("claim_type"), "mean")
    res = run_turn(cid, cases_dir, corpus_path, op, now,
                   her_bubble=ctx.get("bubble"), force_fast=True)
    return res.get("gate") == fx["mechanism"] or res.get("why") == fx["mechanism"]
