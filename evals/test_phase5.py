"""Phase 5 acceptance: the live-thread harness, shrunk stats, full lifecycle."""
import sys, os, json, tempfile, shutil
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from dc_extended.live import (confirm_send, her_reply, shrunk_rate,
                              record_frame_outcome, load_stats)
from dc_extended.switchboard import run_turn
from dc_extended.audit import nightly_checkin
from dc_extended.recovery import record_deviation
from dc_extended.ledger import load_corpus, family_stats
from datetime import datetime

BASE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SRC = os.path.join(BASE, "cases")
CORPUS = os.path.join(BASE, "corpus", "sends.jsonl")
NOW = datetime(2026, 9, 18, 18, 11)


def t(name, cond, detail=""):
    print(("PASS " if cond else "FAIL ") + name + (" | " + str(detail) if detail and not cond else ""))
    if not cond:
        sys.exit(1)


d = tempfile.mkdtemp()
shutil.copytree(SRC, os.path.join(d, "cases"))
cp = os.path.join(d, "sends.jsonl")
shutil.copy(CORPUS, cp)
fstats = os.path.join(d, "frame_stats.json")
rstats = os.path.join(d, "read_stats.json")

# ---------- shrunk_rate ----------
t("no stats -> family prior", shrunk_rate("x", None, 0.6) == 0.6)
t("empty row -> family prior", shrunk_rate("x", {"n": 0}, 0.6) == 0.6)
row = {"n": 10, "extends": 3, "answers": 2, "dry": 5}
r = shrunk_rate("x", row, 0.6, k=5)
t("shrunk math", abs(r - (10 * 0.5 + 5 * 0.6) / 15) < 1e-9, r)
big = {"n": 100, "extends": 50, "answers": 10, "dry": 40}
t("large n approaches own rate", abs(shrunk_rate("x", big, 0.9, k=5) - 0.6) < 0.02,
  shrunk_rate("x", big, 0.9, k=5))

# ---------- record_frame_outcome ----------
record_frame_outcome("tease-mean-1", "extends", stats_path=fstats)
record_frame_outcome("tease-mean-1", "dry", stats_path=fstats)
s = load_stats(fstats)
t("stats accumulate", s["tease-mean-1"] == {"n": 2, "extends": 1, "dry": 1}, s)
t("null frame skipped", record_frame_outcome(None, "extends", stats_path=fstats) is None)

# ---------- full live lifecycle ----------
shutil.copy(os.path.join(SRC, "rae.json"), os.path.join(d, "cases", "rae.json"))
res = run_turn("rae", os.path.join(d, "cases"), cp, "mean", NOW,
               her_bubble="the mean prompt is unhinged", rng_seed=7)
t("turn emits line", res["action"] == "send" and res["frame"]["id"], res.get("frame"))
fid = res["frame"]["id"]

try:
    confirm_send(os.path.join(d, "cases"), "rae", text="a different line")
    t("confirm rejects wrong text", False)
except ValueError:
    t("confirm rejects wrong text", True)

before = load_stats(fstats).get(fid, {"n": 0}).get("n", 0)
confirm_send(os.path.join(d, "cases"), "rae")
out = her_reply(os.path.join(d, "cases"), "rae", "haha yeah you wish", "extends",
                stats_path=fstats, read_stats_path=rstats)
t("reply updates reciprocity", out["reciprocity"] == 1, out)
t("reply attributes frame", out["frame"] == fid, out)
s = load_stats(fstats)
t("frame outcome recorded", s[fid]["n"] == before + 1 and s[fid]["extends"] == s[fid].get("extends", 0), s.get(fid))

# her_reply requires a confirmed send (receipt discipline)
res2 = run_turn("rae", os.path.join(d, "cases"), cp, "mean", NOW,
                her_bubble="the mean prompt again", rng_seed=7)
try:
    her_reply(os.path.join(d, "cases"), "rae", "x", "extends", stats_path=fstats)
    t("reply requires receipt", False)
except ValueError:
    t("reply requires receipt", True)

# pending split read graded by her reply
shutil.copy(os.path.join(SRC, "rae.json"), os.path.join(d, "cases", "rae.json"))
res3 = run_turn("rae", os.path.join(d, "cases"), cp, "mean", NOW,
                her_bubble="busy this week but prove it")
t("split turn went deep", res3["action"] == "deep", res3["action"])
from dc_extended.readlayer import resolve_split
case = json.load(open(os.path.join(d, "cases", "rae.json")))
resolve_split(case, "escalate")
json.dump(case, open(os.path.join(d, "cases", "rae.json"), "w"))
run_turn("rae", os.path.join(d, "cases"), cp, "mean", NOW,
         her_bubble="busy this week but prove it", force_fast=True)
confirm_send(os.path.join(d, "cases"), "rae")
out = her_reply(os.path.join(d, "cases"), "rae", "come over then", "extends",
                stats_path=fstats, read_stats_path=rstats)
rs = json.load(open(rstats))
t("split graded by reply", rs.get("escalate", {}).get("n") == 1 and rs["escalate"]["correct"] == 1, rs)

# ---------- writer consumes shrunk rates ----------
shutil.copy(os.path.join(SRC, "rae.json"), os.path.join(d, "cases", "rae.json"))
record_frame_outcome("tease-mean-1", "dry", stats_path=fstats)
record_frame_outcome("tease-mean-1", "dry", stats_path=fstats)
record_frame_outcome("tease-mean-1", "dry", stats_path=fstats)
record_frame_outcome("tease-mean-1", "dry", stats_path=fstats)
s = load_stats(fstats)
n, good = s["tease-mean-1"]["n"], s["tease-mean-1"].get("extends", 0) + s["tease-mean-1"].get("answers", 0)
own = good / n
fam = family_stats(load_corpus(cp)).get("tease", {}).get("rate", 0.0)
expected = shrunk_rate("tease-mean-1", s["tease-mean-1"], fam, k=5)
res4 = run_turn("rae", os.path.join(d, "cases"), cp, "mean", NOW,
                her_bubble="the mean prompt is unhinged", rng_seed=7,
                frame_stats_path=fstats)
fr = [c["family_rate"] for c in res4["skeleton"]["candidates"] if c.get("frame_id") == "tease-mean-1"]
t("writer uses shrunk rate", fr and abs(fr[0] - expected) < 1e-9, (fr, expected))

# ---------- nightly: honesty over a live ledger ----------
case = json.load(open(os.path.join(d, "cases", "rae.json")))
record_deviation(case, "operator_rejected", mechanism=None, family="wyd_plan")
record_deviation(case, "operator_rejected", mechanism="too_passive", family="plan")
json.dump(case, open(os.path.join(d, "cases", "rae.json"), "w"))
rep = nightly_checkin(os.path.join(d, "cases"), NOW)
t("honesty over live ledger", rep["honesty"]["total"] >= 2 and rep["honesty"]["unknown"] >= 1, rep["honesty"])

print("ALL PHASE 5 TESTS PASS")
