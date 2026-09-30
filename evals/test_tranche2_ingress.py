"""Tranche 2: structured ingress. Every previously-demonstrated collision is a test;
each fix carries its negative control (the playful/nearby case that must NOT trip it)."""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from dc_extended.ingress import parse_operator, split_trigger, parse_trigger

FAIL = []
def t(name, cond, detail=""):
    print(("PASS " if cond else "FAIL ") + name + (" | " + str(detail) if detail and not cond else ""))
    if not cond: FAIL.append(name)

# ---- the six demonstrated substring collisions ----
t("big tease is NOT instagram", parse_operator("big tease").channel is None)
t("make it bolder tonight is NOT instagram", parse_operator("make it bolder tonight").channel is None)
t("handle this better is NOT instagram", parse_operator("handle this better").channel is None)
t("preview this is NOT audit (token, not substring)",
  parse_operator("preview this").audit is False)
t("go overboard is NOT audit (token, not substring)",
  parse_operator("go overboard").audit is False)
t("audit the thread IS audit", parse_operator("audit the thread").audit is True)

# ---- negation scope ----
d = parse_operator("don't set the date")
t("don't set the date -> no plan claim", d.claim_type is None and d.endpoint is None,
  repr(d))
d = parse_operator("not off app")
t("not off app -> no off_app endpoint", d.endpoint is None, repr(d))
d = parse_operator("don't stop texting her")
t("don't stop is NOT a stop", d.stop is False)
d = parse_operator("stop")
t("bare stop IS a stop + close", d.stop and d.endpoint == "close")

# ---- channels distinct ----
t("move to sms targets sms, not ig", parse_operator("move to sms").channel == "sms")
t("get her number targets sms, not ig", parse_operator("get her number").channel == "sms")
t("move to ig targets ig", parse_operator("move to ig").channel == "ig")

# ---- actor separation: HER text is never control ----
her, op = split_trigger("she: stop being cute \U0001f602")
t("her bubble extracted", her == "stop being cute \U0001f602")
t("operator side empty", op is None)
her, d = parse_trigger("she: stop being cute \U0001f602")
t("her 'stop being cute' is NOT a directive", d is None or d.stop is False, repr(d))
her, d = parse_trigger("she: instagram? \nmove to ig")
t("her 'instagram?' is evidence; operator 'move to ig' is the directive",
  her == "instagram?" and d.channel == "ig")
her, d = parse_trigger("she: book me then \U0001f602")
t("her 'book me then' is NOT a plan command", d is None or d.claim_type != "plan", repr(d))

# ---- STOP dominates: locked first, nothing else parsed ----
d = parse_operator("stop, and also move to ig and book wednesday")
t("stop locks: channel/plan never parsed", d.stop and d.channel is None and d.claim_type == "close")

# ---- negative controls ----
t("NEG: 'book wednesday drinks' still parses as plan",
  parse_operator("book wednesday drinks").claim_type == "plan")
t("NEG: bare 'review the case' still audit", parse_operator("review the case").audit is True)
t("NEG: empty text -> empty directive, no guesses",
  parse_operator("").claim_type is None and parse_operator("").stop is False)

print()
print("ALL TRANCHE-2 TESTS PASS" if not FAIL else "FAILURES: %s" % FAIL)
assert not FAIL, FAIL
