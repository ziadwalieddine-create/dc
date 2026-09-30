"""Fix plan verification: one regression per item from the two audits."""
import sys, os, json, tempfile, shutil
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from datetime import datetime
from dc_extended.readlayer import read_bubble_v2
from dc_extended.switchboard import run_turn, classify
from dc_extended.critic import critic, KILL_IDS
from dc_extended.gates import run_gates
from dc_extended.plans import materialize_plan, derive_plan_window
from dc_extended.recovery import record_deviation, forced_resolution, update_reciprocity
from dc_extended.live import confirm_send, her_reply
from dc_extended.audit import charter_query, nightly_checkin
from dc_extended.variety import is_synonym
from dc_extended.frames import certify_corpus
from dc_extended.state import load_case, save_case

BASE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SRC = os.path.join(BASE, "cases")
CORPUS = os.path.join(BASE, "corpus", "sends.jsonl")
NOW = datetime(2026, 9, 18, 18, 11)
WED = datetime(2026, 9, 16, 12, 0)


def t(name, cond, detail=""):
    print(("PASS " if cond else "FAIL ") + name + (" | " + str(detail) if detail and not cond else ""))
    if not cond:
        sys.exit(1)


def fixture():
    d = tempfile.mkdtemp()
    shutil.copytree(SRC, os.path.join(d, "cases"))
    cp = os.path.join(d, "sends.jsonl")
    shutil.copy(CORPUS, cp)
    return d, cp


CASE = {"referents": [{"id": "b", "text": "mean prompt", "status": "live", "type": "bit",
                       "introduced_by": "him", "touched_by_her": True}]}

# ===== FIX 1 — read layer =====
t("F1 dont-be-silly not refuse", read_bubble_v2("don't be silly", CASE)["kind"] not in ("refuse",))
t("F1 I know right not escalate", read_bubble_v2("I know right", CASE)["kind"] != "escalate")
t("F1 between not escalate", read_bubble_v2("between us", CASE)["kind"] != "escalate")
t("F1 no-stop-haha still refuse", read_bubble_v2("no stop it haha", CASE)["kind"] == "refuse")
t("F1 charge still fires", read_bubble_v2("come over tonight", CASE)["kind"] == "escalate")
t("F1 classify habit not tease (boundary)", (classify("that habit of hers") or {}).get("type") is None)
t("F1 classify arbitrary not tease (boundary)", (classify("arbitrary plan") or {}).get("type") is None)

# ===== FIX 2 — clock layer =====
c = {"open_plan": {}}
materialize_plan(c, {"type": "plan", "text": "thursday drinks"})
t("F2 when stored", c["open_plan"].get("when") == "thursday")
derive_plan_window(c, WED)
t("F2 window future Wed", c["open_plan"]["window_future"] is True)
derive_plan_window(c, NOW)  # Friday evening: thursday is dead
t("F2 window dead Fri", c["open_plan"]["window_future"] is False)
d, cp2 = fixture()
res = run_turn("rae", os.path.join(d, "cases"), cp2, "thursday drinks", WED)
t("F2 Wednesday plan reaches writer", res["action"] in ("send", "stop"), res.get("gate") or res.get("why"))
hits = run_gates({"open_plan": {"when": "thursday"}, "closed": False}, "", NOW,
                 claim={"type": "plan", "family": "plan"}, all_cases=[])
t("F2 dead-Friday plan gated", any(g == "dead_clock" for g, _ in hits), hits)
t("F2 critic clock live", critic("friday still works", {"_now": NOW, "unique_facts": []}) == "clock")

# ===== FIX 7 — critic kills =====
t("F7 27 kills", len(KILL_IDS) == 27)
c7 = {"_claim": {"family": "basil"}, "banned": [{"family": "basil"}]}
t("F7 negotiates_no from banned", critic("the same bit again", c7) == "negotiates_no")
c7b = {"_policy": "Play", "_amplitude": 0}
t("F7 missed_heat wired", critic("flat", c7b, read={"she_extended": True}) == "missed_heat")
t("F7 dead_end wired", critic("closed joke", {"_policy": "Play"},
                              sim_classes={"silence", "topic_change"}, other_wider=True) == "dead_end")
t("F7 corpus certifies", all(v == [] for v in certify_corpus(CORPUS).values()))

# ===== FIX 3 — seams =====
d, cp3 = fixture()
res = run_turn("rae", os.path.join(d, "cases"), cp3, None, NOW,
               trigger_text="she: the mean prompt is unhinged\nplay the bit", rng_seed=5)
t("F3 trigger_text parses", res["action"] == "send", res.get("gate") or res.get("why"))
t("F3 claim from trigger", res["claim"] == "bit", res.get("claim"))
c3 = {"referents": [{"id": "b", "text": "mean prompt", "status": "live"}],
      "claim": {"type": "bit"}, "unique_facts": ["mean prompt"], "amp_cap": 0}
deg = forced_resolution(c3, "dead_clock")
t("F3 dead_clock degrades", deg is not None and deg["amplitude"] == 0, deg)
t("F3 integrity never degrades", forced_resolution(c3, "stated_boundary") is None)
t("F3 forced no side-effect cap", "amp_cap" in c3 and c3["amp_cap"] == 0)  # untouched by the check itself
t("F3 is_synonym guard", is_synonym("the mean prompt is a lot", "the mean prompt is a lot tall"))

# ===== FIX 4 — gate honesty =====
log = []
hits = run_gates({"closed": True, "open_plan": None}, "", NOW, claim=None, all_cases=[], log=log)
t("F4 FAIL recorded", any(g["gate"] == "stale_goal" and g["result"] == "FAIL" for g in log), log[:2])
t("F4 PASS recorded too", any(g["result"] == "PASS" for g in log))
d, cp4 = fixture()
res = run_turn("camille", os.path.join(d, "cases"), cp4, "wyd on friday", NOW)
skel_gates = res["skeleton"]["gates"]
t("F4 skeleton has FAIL", any(g["result"] == "FAIL" for g in skel_gates),
  [(g["gate"], g["result"]) for g in skel_gates])

# ===== FIX 5 — live integrity =====
d, cp5 = fixture()
res = run_turn("rae", os.path.join(d, "cases"), cp5, "mean", NOW,
               her_bubble="the mean prompt is unhinged", rng_seed=7)
confirm_send(os.path.join(d, "cases"), "rae")
o1 = her_reply(os.path.join(d, "cases"), "rae", "haha yeah you wish", "extends", stats_path=os.path.join(d, "frame_stats.json"), read_stats_path=os.path.join(d, "read_stats.json"))
o2 = her_reply(os.path.join(d, "cases"), "rae", "haha yeah you wish", "extends", stats_path=os.path.join(d, "frame_stats.json"), read_stats_path=os.path.join(d, "read_stats.json"))
t("F5 duplicate detected", o2.get("duplicate") is True, o2)
case5 = load_case(os.path.join(d, "cases", "rae.json"))
t("F5 no double count", case5["reciprocity"] == 1, case5["reciprocity"])
t("F5 charter suppression", charter_query([{"family": "wyd_plan"}, {"family": "wyd_plan"}],
                                          banned_families={"wyd_plan"}) == [])

# ===== FIX 6 — hygiene =====
d, cp6 = fixture()
res = run_turn("rae", os.path.join(d, "cases"), cp6, "mean", NOW,
               her_bubble="the mean prompt is unhinged", rng_seed=7)
case6 = load_case(os.path.join(d, "cases", "rae.json"))
unders = [k for k in case6 if k.startswith("_")]
t("F6 no transient keys persisted", not unders, unders)

print("ALL FIX TESTS PASS")
