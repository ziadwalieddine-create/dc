"""Phase 4 acceptance: self-audit over the deviation ledger."""
import sys, os, json, tempfile, shutil
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from dc_extended.recovery import record_deviation
from dc_extended.audit import (collect_deviations, honesty, charter_query,
                               gate_hit_rates, nightly_checkin, deviation_to_fixture,
                               fixtures_from, replay_fixture)
from dc_extended.frames import justification_card, card_for
from dc_extended.switchboard import run_turn
from datetime import datetime

BASE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SRC = os.path.join(BASE, "cases")
CORPUS = os.path.join(BASE, "corpus", "sends.jsonl")
NOW = datetime(2026, 9, 18, 18, 11)


def t(name, cond, detail=""):
    print(("PASS " if cond else "FAIL ") + name + (" | " + str(detail) if detail and not cond else ""))
    if not cond:
        sys.exit(1)


# ---------- gate_stop event + context ----------
c = {"referents": [], "reciprocity": 0, "outcomes": [], "claim": {"type": "plan"}}
row = record_deviation(c, "gate_stop", mechanism="dead_clock", family="plan",
                       context={"bubble": "friday?", "claim_type": "plan", "channel": "hinge"})
t("gate_stop row records mechanism", row["mechanism"] == "dead_clock")
t("context stored", row["context"]["bubble"] == "friday?")
t("row in case ledger", c["deviations"][-1] is row or c["deviations"][-1]["mechanism"] == "dead_clock")

# ---------- honesty ----------
rows = [{"mechanism": None}, {"mechanism": "x"}, {"mechanism": None}, {"mechanism": "y"}]
h = honesty(rows)
t("honesty ratio", h == {"unknown": 2, "total": 4, "ratio": 0.5}, h)
t("honesty empty", honesty([])["ratio"] == 0.0)

# ---------- charter query ----------
t("two incidents propose", charter_query(rows[:0] + [{"family": "wyd_plan"}, {"family": "wyd_plan"}]) == ["wyd_plan"])
t("one incident does not", charter_query([{"family": "wyd_plan"}]) == [])
t("mixed families", charter_query([{"family": "a"}, {"family": "a"}, {"family": "b"}]) == ["a"])

# ---------- gate hit rates ----------
gr = gate_hit_rates([{"event": "gate_stop", "mechanism": "restack"},
                     {"event": "gate_stop", "mechanism": "restack"},
                     {"event": "gate_stop", "mechanism": "dead_clock"},
                     {"event": "operator_rejected", "mechanism": "x"}])
t("hit rates counted", gr == {"restack": 2, "dead_clock": 1}, gr)

# ---------- nightly checkin ----------
d = tempfile.mkdtemp()
shutil.copytree(SRC, os.path.join(d, "cases"))
cp = os.path.join(d, "sends.jsonl")
shutil.copy(CORPUS, cp)
# cause two restack deviations on the same family across two cases
for cid in ["rae", "jisu"]:
    p = os.path.join(d, "cases", cid + ".json")
    case = json.load(open(p))
    case.setdefault("banned_families", [])
    if "wyd_plan" not in case["banned_families"]:
        case["banned_families"].append("wyd_plan")
    record_deviation(case, "gate_stop", mechanism="restack", family="wyd_plan",
                     context={"bubble": "wyd on friday", "claim_type": "plan", "channel": "hinge",
                              "operator_text": "wyd on friday", "banned_families": ["wyd_plan"]})
    json.dump(case, open(p, "w"))

rep = nightly_checkin(os.path.join(d, "cases"), NOW)
for k in ("queue", "honesty", "proposed_gates", "gate_hit_rates"):
    t("nightly has " + k, k in rep)
t("nightly proposes the family", "wyd_plan" in rep["proposed_gates"], rep["proposed_gates"])
t("nightly honesty computed", rep["honesty"]["total"] >= 2, rep["honesty"])
t("nightly gate rates", rep["gate_hit_rates"].get("restack") == 2)

# collect_deviations walks the dir
collected = collect_deviations(os.path.join(d, "cases"))
t("collect walks cases", len([r for r in collected if r.get("event") == "gate_stop"]) == 2, len(collected))

# ---------- deviations -> fixtures -> replay ----------
fx = fixtures_from(collected)
t("fixtures extracted", len(fx) == 2 and all(f["mechanism"] == "restack" for f in fx), fx)
t("insufficient context stays out",
  deviation_to_fixture({"event": "gate_stop", "mechanism": "restack", "context": {}}) is None)
t("non-gate events stay out", deviation_to_fixture({"event": "operator_rejected", "mechanism": "x", "context": {"bubble": "b"}}) is None)

ok = all(replay_fixture(f, os.path.join(d, "cases"), cp, NOW) for f in fx)
t("recorded real failure replays as regression guard", ok)

# ---------- justification cards ----------
card = justification_card({"id": "tease-mean-1", "slots": ["referent", "challenge"],
                           "forecast": {"extends": 3}, "protected": ["mean"]})
t("card line per slot", len(card["lines"]) == 2 and "referent" in card["lines"][0], card)
t("unknown slot honest", "declared slot, reason at authoring time" in
  justification_card({"id": "x", "slots": ["mystery"]})["lines"][0])

c4 = card_for("tease-mean-1", CORPUS)
t("card_for real frame", c4 and c4["frame"] == "tease-mean-1" and len(c4["lines"]) == 2)
t("card_for unknown -> None", card_for("nope", CORPUS) is None)

# ---------- run_turn wires it: gate stop writes deviation; send carries card ----------
shutil.copy(os.path.join(SRC, "rae.json"), os.path.join(d, "cases", "rae.json"))
res = run_turn("rae", os.path.join(d, "cases"), cp, "wyd on friday", NOW)
t("gate stop on wyd_plan", res["action"] == "stop" and res.get("gate") in ("restack", "dead_clock"), res.get("gate"))
rae = json.load(open(os.path.join(d, "cases", "rae.json")))
devs = [r for r in rae.get("deviations", []) if r.get("event") == "gate_stop"]
t("gate hit written to ledger", len(devs) == 1 and devs[0]["mechanism"] == res.get("gate"), devs)

shutil.copy(os.path.join(SRC, "rae.json"), os.path.join(d, "cases", "rae.json"))
res2 = run_turn("rae", os.path.join(d, "cases"), cp, "mean", NOW,
                her_bubble="the mean prompt is unhinged", rng_seed=7)
t("send carries justification card", res2["skeleton"].get("frame", {}).get("card", {}).get("frame") ==
  res2["frame"]["id"], res2["skeleton"].get("frame"))

print("ALL PHASE 4 TESTS PASS")
