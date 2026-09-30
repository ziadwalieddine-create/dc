"""Phase 2 acceptance: endpoint lock, endpoint-weighted ranking, triage."""
import sys, os, json, tempfile, shutil
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from dc_extended.goals import lock_endpoint, sim_rank, triage, queue_weight, BASE_SIM
from dc_extended.switchboard import run_turn
from dc_extended.queue import build_queue
from datetime import datetime

SRC = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "cases")
NOW = datetime(2026, 9, 18, 18, 11)


def t(name, cond, detail=""):
    print(("PASS " if cond else "FAIL ") + name + (" | " + str(detail) if detail and not cond else ""))
    if not cond:
        sys.exit(1)


# ---------- lock_endpoint ----------
t("get to know -> banter", lock_endpoint("get to know her", None)["type"] == "banter")
t("off app phrase", lock_endpoint("move off app", None)["type"] == "off_app")
t("lock the date", lock_endpoint("lock the date", None)["type"] == "date_lock")
t("stored survives silent turn", lock_endpoint(None, {"type": "date_lock", "source": "x"})["type"] == "date_lock")
t("operator wins over stored", lock_endpoint("get to know her", {"type": "date_lock"})["type"] == "banter")
t("unknown text -> None (UNRESOLVED legal)", lock_endpoint("mean", None)["type"] is None)
t("source tracked", lock_endpoint("off app", None)["source"] == "operator_this_turn")

# ---------- sim_rank ----------
t("None -> base", sim_rank(None) == BASE_SIM)
t("banter boosts extends", sim_rank("banter")["extends"] == BASE_SIM["extends"] + 1)
t("banter demotes answers", sim_rank("banter")["answers"] == BASE_SIM["answers"] - 1)
t("date_lock boosts counters", sim_rank("date_lock")["counters"] == BASE_SIM["counters"] + 2)
t("date_lock demotes dry", sim_rank("date_lock")["dry"] == BASE_SIM["dry"] - 1)
t("does not mutate base", BASE_SIM["extends"] == 4)

# ---------- triage ----------
def tc(rec=0, outcomes=None, closed=False):
    return {"reciprocity": rec, "outcomes": outcomes or [], "closed": closed}

t("closed", triage(tc(closed=True)) == "closed")
t("refuses -> unreceptive", triage(tc(0, [{"response": "refuses"}])) == "unreceptive")
t("rec >= 2 -> receptive", triage(tc(2, [{"response": "extends"}])) == "receptive")
t("rec <= -2 -> cooling", triage(tc(-2, [{"response": "dry"}, {"response": "dry"}])) == "cooling")
t("two flat -> cooling", triage(tc(0, [{"response": "dry"}, {"response": "silence"}])) == "cooling")
t("default neutral", triage(tc(0, [{"response": "extends"}])) == "neutral")
t("single dry not cooling", triage(tc(0, [{"response": "dry"}])) == "neutral")
t("weights ordered", queue_weight(tc(2)) < queue_weight(tc(0)) < queue_weight(tc(-2)))

# ---------- run_turn: cooling gate + endpoint lock ----------
d = tempfile.mkdtemp()
shutil.copytree(SRC, os.path.join(d, "cases"))
cp = os.path.join(d, "sends.jsonl")
shutil.copy(os.path.join(os.path.dirname(SRC), "corpus", "sends.jsonl"), cp)
rp = os.path.join(d, "cases", "rae.json")
rae = json.load(open(rp))

# cooling: two flat outcomes, no extension
rae["outcomes"] = [{"response": "dry"}, {"response": "silence"}]
rae["reciprocity"] = -2
rae["last_send"] = {"text": "x", "confirmed": True}
json.dump(rae, open(rp, "w"))
res = run_turn("rae", os.path.join(d, "cases"), cp, "mean", NOW, her_bubble="haha")
t("cooling -> stop no line", res["action"] == "stop" and res["why"] == "cooling" and res.get("text") is None, res)
t("skeleton shows triage", res["skeleton"].get("triage") == "cooling")

# cooling lifted by her extension
res = run_turn("rae", os.path.join(d, "cases"), cp, "mean", NOW,
               her_bubble="the mean prompt is unhinged")
t("extension lifts cooling", res["action"] == "send", res)

# force_fast bypasses cooling (operator override) — confirm prior send first
from dc_extended.live import confirm_send as _cs
_cs(os.path.join(d, "cases"), "rae")
res = run_turn("rae", os.path.join(d, "cases"), cp, "mean", NOW, her_bubble="haha", force_fast=True)
t("force_fast overrides cooling", res["action"] == "send", res.get("gate") or res.get("why"))

# endpoint lock: operator's goal this turn wins (substituted_end)
rae = json.load(open(rp))
rae["endpoint"] = {"type": "date_lock", "source": "stored"}
rae["outcomes"] = [{"response": "extends"}]
rae["reciprocity"] = 1
json.dump(rae, open(rp, "w"))
res = run_turn("rae", os.path.join(d, "cases"), cp, "get to know her, mean", NOW,
               her_bubble="the mean prompt is unhinged")
t("operator endpoint wins", res["skeleton"].get("endpoint") == "banter", res["skeleton"].get("endpoint"))
rae2 = json.load(open(rp))
t("endpoint stored with source", rae2["endpoint"]["type"] == "banter" and rae2["endpoint"]["source"] == "operator_this_turn", rae2["endpoint"])

# silent turn keeps stored endpoint
res = run_turn("rae", os.path.join(d, "cases"), cp, "mean", NOW, her_bubble="the mean prompt again")
t("stored endpoint survives", res["skeleton"].get("endpoint") == "banter", res["skeleton"].get("endpoint"))

# explicit endpoint kwarg
res = run_turn("rae", os.path.join(d, "cases"), cp, "mean", NOW,
               her_bubble="the basil plant again", endpoint="off app")
t("endpoint kwarg locks", res["skeleton"].get("endpoint") == "off_app", res["skeleton"].get("endpoint"))

# ---------- queue: triage tiebreak ----------
def mk(cid, rec):
    return {"id": cid, "stage": "plan_pending", "updated": "2026-09-17T10:00:00",
            "reciprocity": rec, "outcomes": [], "closed": False,
            "open_plan": {"confirmed": False, "proposed_ts": "2026-09-17T10:00:00"}}
qq = build_queue([mk("cool", -2), mk("hot", 2)], NOW)
ids = [i["id"] for i in qq]
t("receptive before cooling at same deadline", ids.index("hot") < ids.index("cool"), ids)

print("ALL PHASE 2 TESTS PASS")
