"""PROMOTION: the controller serves the primary path. run_primary releases a
Move on the canonical turn; explicit STOP is hard-gated; legacy remains callable."""
import sys, os, tempfile
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from dc_extended.controller import Controller
from datetime import datetime

FAIL = []
def t(name, cond, detail=""):
    print(("PASS " if cond else "FAIL ") + name + (" | " + str(detail) if detail and not cond else ""))
    if not cond: FAIL.append(name)

d = tempfile.mkdtemp()
CASE = {"id": "gabby", "name": "Gabrielle-Rose", "channel": "ig", "closed": False,
        "boundaries": ["ex"], "unique_facts": ["Wynyard Park", "Wynyard station"],
        "referents": [{"id": "w", "text": "Wynyard Park", "status": "live", "type": "plan",
                       "introduced_by": "him", "touched_by_her": True}]}
NOW = datetime(2026, 9, 20, 13, 0)
ctl = Controller(os.path.join(d, "gabby.jsonl"))

res = ctl.run_primary(dict(CASE), "Where the heck is Wynyard Park 🤣 😂", None, now=NOW)
t("primary: Wynyard turn releases a two-bubble Move", res["action"] == "send"
  and len(res["move"]["sequence"]) == 2, res.get("why"))
t("primary: truth stream recorded the proposal", res["shadow_proposal_id"] is not None)
t("primary: event chain verifies", ctl.log.verify())

res2 = ctl.run_primary(dict(CASE), "the basil plant is unhinged", "stop", now=NOW)
t("primary: STOP hard-gated, closes", res2["action"] == "stop" and res2.get("closed") is True)

t("legacy still importable as reference", hasattr(ctl, "run_shadow"))

print()
print("PROMOTION PASS — v2 is the primary path" if not FAIL else "FAILURES: %s" % FAIL)
assert not FAIL, FAIL
