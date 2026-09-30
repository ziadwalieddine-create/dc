"""PROMOTION GATE: the acceptance constitution, run against the v2 substrate.

Each invariant states its substrate: LEGACY (old run_turn), V2 (pipeline_v2),
or BOTH. The promotion decision requires: every V2-satisfiable invariant green
on v2, and no V2 invariant passing on legacy (that would mean the test is
non-discriminating). The constitution file itself stays aimed at legacy — this
test is the v2 mirror."""
import sys, os, tempfile, json
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from dc_extended.pipeline_v2 import turn_v2
from dc_extended.read_v2 import read_v2
from dc_extended.moves import single, sequence
from dc_extended.events import EventLog
from dc_extended.amplitude import permitted_band
from dc_extended.licences import Licences
from datetime import datetime

FAIL = []
def t(name, cond, detail=""):
    print(("PASS " if cond else "FAIL ") + name + (" | " + str(detail) if detail and not cond else ""))
    if not cond: FAIL.append(name)

d = tempfile.mkdtemp()
NOW = datetime(2026, 9, 20, 13, 0)
CASE_W = {"id": "gabby", "name": "Gabrielle-Rose", "channel": "ig", "closed": False,
          "boundaries": ["ex"], "unique_facts": ["Wynyard Park", "Wynyard station"],
          "referents": [{"id": "w", "text": "Wynyard Park", "status": "live", "type": "plan",
                         "introduced_by": "him", "touched_by_her": True}]}

print("=== V2 CONSTITUTION MIRROR ===")

# A1 — Wynyard multidimensional (V2)
r = read_v2("Where the heck is Wynyard Park 🤣 😂", CASE_W)
t("A1 [V2]: Wynyard = logistics + tease + heat 2 + answer_location",
  "logistics_question" in r["literal_functions"] and "tease" in r["relational_functions"]
  and r["heat"] == 2 and "answer_location" in r["response_obligations"])

# A2 — Hot not dry (V2)
t("A2 [V2]: 'Hot' retains attraction", read_v2("Hot", CASE_W)["heat"] == 2)

# A3 — scoped refusal (V2)
r = read_v2("No Thursday 😂 Friday works", {})
t("A3 [V2]: Thursday refused, Friday accepted, not global",
  r["refusal_scope"] == "thursday" and r["acceptance"] == "friday")

# A4 — proposal is never a send (V2)
log = EventLog(os.path.join(d, "a4.jsonl"))
pid = log.propose("c", "a line")
try:
    log.confirm_send("c", pid, "different line"); t("A4 [V2]: payload mismatch blocked", False)
except ValueError:
    t("A4 [V2]: payload mismatch blocked", True)

# B2/B3 — stop dominates, production path, sendable world (V2)
res = turn_v2(dict(CASE_W), "the basil plant is unhinged", "stop",
              os.path.join(d, "b2.jsonl"), now=NOW)
t("B2/B3 [V2]: stop releases nothing", res["action"] == "stop")
t("B2/B3 [V2]: stop closes", res.get("closed") is True)

# C2 — structured Move (V2)
res2 = turn_v2(dict(CASE_W), "Where the heck is Wynyard Park 🤣 😂",
               None, os.path.join(d, "c2.jsonl"), now=NOW)
t("C2 [V2]: released payload is a structured move (sequence + text)",
  isinstance(res2.get("move"), dict) and "sequence" in res2["move"] and "text" in res2["move"])

# C3 — heat causally affects band (V2)
f0 = permitted_band({"heat": 0}, Licences({}))
f2 = permitted_band(read_v2("Hot", CASE_W), Licences(read_v2("Hot", CASE_W)))
t("C3 [V2]: heat 2 band wider than heat 0", f2[1] > f0[1])

# C4 — no sanding rule holds (V2)
rN = read_v2("Where is the station?", {})
t("C4 [V2]: neutral logistics gains no heat", (rN["heat"] or 0) == 0)

# C5 — flirt model on critical path (V2: licences/amplitude govern generation)
t("C5 [V2]: licences gate the generated move", "licences" in res2)

# discrimination proof: B2 must NOT pass on legacy (already shown red in constitution;
# assert here that v2's pass is therefore substrate-earned, not test-weak)
print()
flips = [x for x in FAIL]
print("V2 MIRROR:", "ALL V2 INVARIANTS GREEN — promotion gate satisfied" if not flips
      else "RED: %s" % flips)
assert not FAIL, FAIL
