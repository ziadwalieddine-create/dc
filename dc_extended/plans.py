"""dc_extended.plans — store the day-word fact, derive the clock every turn.

window_future is a DERIVED value (when + now), never stored. A stored bool
would go stale by Friday — that is dc's Jisu failure in miniature.
"""
from .gates import WEEKDAYS, day_dead


def materialize_plan(case, claim):
    """claim.type == 'plan' and a day word in the claim text -> store when only."""
    if (claim or {}).get("type") != "plan":
        return case
    text = ((claim or {}).get("text") or "").lower()
    op = case.get("open_plan") or {}
    case["open_plan"] = op
    if op.get("when"):
        return case
    for w in WEEKDAYS:
        if w in text:
            op["when"] = w
            break
    return case


def derive_plan_window(case, now):
    """window_future computed fresh per turn. None when no day named -> Hold.
    Never manufactures an open_plan where there was none (queue safety)."""
    op = case.get("open_plan")
    if not op:
        return case
    when = op.get("when")
    op["window_future"] = (not day_dead(when, now)) if when else None
    case["open_plan"] = op
    return case
