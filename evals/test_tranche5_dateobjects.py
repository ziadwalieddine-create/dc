"""Tranche 5: canonical dating objects. Every audit finding in this cluster is a test."""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from datetime import datetime
from dc_extended.dateobjects import (Plan, ChannelTransition, StateDelta,
                                     resolve_day, no_duplicate_plan, final_audit)

FAIL = []
def t(name, cond, detail=""):
    print(("PASS " if cond else "FAIL ") + name + (" | " + str(detail) if detail and not cond else ""))
    if not cond: FAIL.append(name)

SUN = datetime(2026, 9, 20, 13, 0)   # a Sunday
TUE = datetime(2026, 9, 22, 13, 0)

# --- finding 27: the gold Sunday->Wednesday progression is LEGAL ---
wed = resolve_day("wednesday", SUN)
t("gold: Sunday->Wednesday resolves to the coming Wednesday", wed.weekday() == 2 and wed.date() > SUN.date(), wed.date())

# --- finding 28: Monday on Tuesday means NEXT Monday ---
mon = resolve_day("monday", TUE)
t("Monday-on-Tuesday = NEXT Monday", mon.date() > TUE.date() and mon.weekday() == 0, mon.date())

# --- absolute datetime before any confirmation scheduling ---
p = Plan("drinks", absolute_dt=datetime(2026, 9, 23, 17, 50), venue="Wynyard", meeting_point="York St side")
t("plan: morning confirmation due is representable", p.morning_confirmation_due is not None
  and p.morning_confirmation_due.date() == p.absolute_dt.date())
t("plan: lifecycle transitions legal", p.advance("accepted").advance("locked").status == "locked")
try:
    p.advance("proposed"); t("NEG: backwards lifecycle rejected", False)
except ValueError:
    t("NEG: backwards lifecycle rejected", True)

# --- finding 11: channel transition is legal, never wrong_platform ---
ct = ChannelTransition("hinge", "ig")
t("transition: Hinge->IG is progression", ct.from_channel == "hinge" and ct.accept() == "ig")
try:
    ChannelTransition("ig", "ig"); t("NEG: same-endpoint 'transition' rejected", False)
except ValueError:
    t("NEG: same-endpoint 'transition' rejected", True)

# --- finding 23: propose -> validate -> commit ---
case = {"open_plan": None, "closed": False, "boundaries": ["ex"]}
delta = StateDelta({"open_plan": "drinks wednesday"})
ok, reasons = delta.commit(case, gates=[no_duplicate_plan])
t("commit: legal delta applies", ok and case["open_plan"] == "drinks wednesday")
case2 = {"open_plan": {"what": "dinner"}, "closed": False, "boundaries": []}
ok2, reasons2 = StateDelta({"open_plan": "drinks"}).commit(case2, gates=[no_duplicate_plan])
t("dual_plan: illegal delta leaves state UNTOUCHED", (not ok2) and case2["open_plan"]["what"] == "dinner", reasons2)

# --- final audit on the EXACT rendered payload ---
case_g = {"id": "gabby", "name": "Gabrielle-Rose", "closed": False, "boundaries": ["ex"]}
case_k = {"id": "kathe", "name": "Katherine", "identity_fingerprint": ["Damien Hirst"]}
ok, fails = final_audit("Right outside Wynyard station. You'll find me \U0001f602", case_g, [case_k])
t("final audit: clean payload passes", ok, fails)
ok, fails = final_audit("meet me at Damien Hirst's gallery", case_g, [case_k])
t("NEG: cross-thread fingerprint in RENDERED text is caught", not ok and any("fingerprint" in f for f in fails), fails)
ok, fails = final_audit("hey Katherine \U0001f602", case_g, [case_k])
t("NEG: wrong woman in RENDERED text is caught", not ok and any("wrong woman" in f for f in fails), fails)
case_g["closed"] = True
ok, fails = final_audit("any line at all", case_g, [])
t("NEG: closed case cannot release", not ok, fails)

print()
print("ALL TRANCHE-5 TESTS PASS" if not FAIL else "FAILURES: %s" % FAIL)
assert not FAIL, FAIL
