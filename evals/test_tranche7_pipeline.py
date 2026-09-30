"""Tranche 7: end-to-end new substrate in shadow. The Wynyard turn, composed:
answer the literal obligation AND preserve the heat, in one Move, on the truth
stream, with the final audit on the exact payload."""
import sys, os, tempfile, json
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from dc_extended.pipeline_v2 import turn_v2
from datetime import datetime

FAIL = []
def t(name, cond, detail=""):
    print(("PASS " if cond else "FAIL ") + name + (" | " + str(detail) if detail and not cond else ""))
    if not cond: FAIL.append(name)

d = tempfile.mkdtemp()
CASE = {"id": "gabby", "name": "Gabrielle-Rose", "channel": "ig", "closed": False,
        "boundaries": ["ex"], "unique_facts": ["Wynyard Park", "Wynyard station"],
        "referents": [{"id": "w", "text": "Wynyard Park", "status": "live", "type": "plan",
                       "introduced_by": "him", "touched_by_her": True}],
        "location_fact": "Right outside Wynyard station"}
NOW = datetime(2026, 9, 20, 13, 0)

# --- the canonical composed turn ---
res = turn_v2(dict(CASE), "Where the heck is Wynyard Park 🤣 😂",
              None, os.path.join(d, "gabby.jsonl"), now=NOW)
t("send: composed turn releases a Move", res["action"] == "send", res.get("why"))
t("obligation satisfied: location answer is bubble 1",
  "Wynyard" in res["move"]["sequence"][0], res["move"]["sequence"])
t("heat preserved: playful continuation is bubble 2", "find me" in res["move"]["sequence"][1])
t("read preserved all Wynyard dimensions",
  "logistics_question" in res["read"]["literal_functions"]
  and "tease" in res["read"]["relational_functions"] and res["read"]["heat"] == 2)
t("band: heat 2 opened the ceiling", res["band"][1] == 2, res["band"])
t("truth stream: proposal recorded", res["shadow_proposal_id"] is not None)
t("delta committed", res["delta_committed"] is True)

# --- STOP dominates on the new path too ---
res2 = turn_v2(dict(CASE), "the basil plant is unhinged", "stop",
               os.path.join(d, "g2.jsonl"), now=NOW)
t("NEG: explicit STOP releases nothing, closes", res2["action"] == "stop"
  and res2.get("closed") is True, res2.get("why"))

# --- no heat, no licence: no move, honest stop ---
res3 = turn_v2({"id": "x", "closed": False, "boundaries": [], "unique_facts": [],
                "referents": []}, "haha", None, os.path.join(d, "x.jsonl"), now=NOW)
t("NEG: dry token with nothing live -> honest stop, no invented line",
  res3["action"] == "stop", res3.get("why"))

# --- final audit blocks a boundary in the rendered payload ---
CASE_B = dict(CASE); CASE_B["boundaries"] = ["Wynyard"]
res4 = turn_v2(dict(CASE_B), "Where even is Wynyard Park 🤣", None,
               os.path.join(d, "g4.jsonl"), now=NOW,
               corpus_path=os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "corpus", "sends.jsonl"))
t("NEG: rendered payload touching a boundary is stopped at final audit",
  res4["action"] == "stop" and res4.get("why") == "final_audit", res4.get("failures"))

# --- event chain integrity across all turns ---
from dc_extended.events import EventLog
allok = all(EventLog(os.path.join(d, f)).verify() for f in os.listdir(d))
t("all event chains verify", allok)

print()
print("ALL TRANCHE-7 TESTS PASS" if not FAIL else "FAILURES: %s" % FAIL)
assert not FAIL, FAIL
