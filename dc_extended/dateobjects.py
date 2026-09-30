"""dc_extended.dateobjects — canonical dating objects.

Plan: absolute datetime + venue + meeting point + full lifecycle. Weekday WORDS
resolve to absolute dates (Sunday->Wednesday is legal; Monday-on-Tuesday means
NEXT Monday; missed windows expire instead of poisoning the weekday forever).

ChannelTransition: from/to/status — source != destination IS the transition,
never a wrong_platform gate.

StateDelta: propose -> validate -> commit. Nothing durable mutates before
legality. An unsent candidate leaves conversation state unchanged.

final_audit: runs on the EXACT rendered payload, never a precursor.
"""
from datetime import datetime, timedelta

DAYS = {"monday":0,"tuesday":1,"wednesday":2,"thursday":3,"friday":4,"saturday":5,"sunday":6}
PLAN_LIFECYCLE = ("proposed","countered","accepted","locked","confirmation_due",
                  "confirmed","completed","cancelled","missed","failed","repaired")


def resolve_day(word, now):
    """Weekday word -> absolute date. Always the NEXT occurrence strictly after
    today; today itself only if the hour is still ahead (handled by caller window).
    Sunday->Wednesday is legal (finding 27); Monday on Tuesday = next Monday (28)."""
    target = DAYS[word]
    delta = (target - now.weekday()) % 7
    if delta == 0:
        delta = 7          # same weekday means NEXT week unless time still ahead
    day = (now + timedelta(days=delta)).date()
    return datetime.combine(day, datetime.min.time())


class Plan:
    def __init__(self, what=None, absolute_dt=None, venue=None, meeting_point=None):
        if absolute_dt is not None and not isinstance(absolute_dt, datetime):
            raise ValueError("Plan requires an absolute datetime")
        self.what = what
        self.absolute_dt = absolute_dt
        self.venue = venue
        self.meeting_point = meeting_point
        self.status = "proposed" if absolute_dt else "proposed"
        self.confirmation_due = None
        if absolute_dt:
            self.confirmation_due = absolute_dt.replace(hour=9, minute=0)  # morning-of

    def legal(self):
        return self.absolute_dt is not None and self.status in PLAN_LIFECYCLE

    def advance(self, new_status):
        if new_status not in PLAN_LIFECYCLE:
            raise ValueError("unknown plan status: %r" % new_status)
        i, j = PLAN_LIFECYCLE.index(self.status), PLAN_LIFECYCLE.index(new_status)
        if j < i and new_status != "repaired":
            raise ValueError("illegal transition %s -> %s" % (self.status, new_status))
        self.status = new_status
        return self

    @property
    def morning_confirmation_due(self):
        """The morning-of rule is now representable: an absolute datetime."""
        return self.confirmation_due if self.absolute_dt else None


class ChannelTransition:
    """Hinge -> IG is progression within one case (a Person->Thread model), never
    wrong_platform. source != destination is the POINT of a transition."""
    def __init__(self, from_channel, to):
        if from_channel == to:
            raise ValueError("a transition requires different endpoints")
        self.from_channel = from_channel
        self.to = to
        self.status = "requested"

    def accept(self):
        self.status = "accepted"
        return self.to


class StateDelta:
    """propose -> validate -> commit. validate returns (ok, reasons)."""
    def __init__(self, mutations=None):
        self.mutations = mutations or {}   # {case_field: new_value}

    def validate(self, case, gates=None):
        reasons = []
        for g in (gates or []):
            ok, why = g(case, self.mutations)
            if not ok:
                reasons.append(why)
        return (not reasons, reasons)

    def commit(self, case, gates=None):
        ok, reasons = self.validate(case, gates)
        if not ok:
            return False, reasons          # durable state UNTOUCHED
        case.update(self.mutations)
        return True, []


def no_duplicate_plan(case, mutations):
    if "open_plan" in mutations and case.get("open_plan"):
        return False, "dual_plan: an open plan exists; coordinate it, don't stack"
    return True, ""


def final_audit(rendered_payload, case, other_cases=None):
    """Audit the EXACT rendered text, not a precursor. Returns (ok, failures)."""
    fails = []
    text = rendered_payload if isinstance(rendered_payload, str) else " ".join(rendered_payload)
    low = text.lower()
    if case.get("closed"):
        fails.append("closed case cannot release")
    for b in case.get("boundaries", []):
        bl = b.lower()
        # phrase-level semantics: multi-word boundaries match as phrases;
        # single-word boundaries match whole tokens; "ex" never fires on "sexy"
        if " " in bl:
            if bl in low:
                fails.append("boundary touched: %s" % b)
        else:
            import re as _re
            if _re.search(r"\b%s\b" % _re.escape(bl), low):
                fails.append("boundary touched: %s" % b)
    # "home" boundary family: any come-over/come-home phrasing trips it
    for hb in case.get("home_boundaries", []):
        if hb in low:
            fails.append("home boundary: %s" % hb)
    for oc in (other_cases or []):
        if oc.get("id") == case.get("id"):
            continue
        for f in (oc.get("identity_fingerprint") or []):
            if f and f.lower() in low:
                fails.append("cross-thread fingerprint (%s): %s" % (oc.get("id"), f))
        if (oc.get("name") or "").lower() and (oc.get("name") or "").lower() in low:
            fails.append("wrong woman: %s" % oc.get("name"))
    return (not fails, fails)
