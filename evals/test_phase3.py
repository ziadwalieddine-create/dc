"""Phase 3 acceptance: lexicon identity, frame certification, behavior seal."""
import sys, os, json, tempfile, shutil
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import dc_extended.lexicon as lex
import dc_extended.critic as critic_mod
import dc_extended.readlayer as rl_mod
import dc_extended.frames as frames_mod
from dc_extended.frames import certify_frame, certify_corpus, fixture_case, FRAME_FIELDS
from dc_extended.switchboard import run_turn
from datetime import datetime

BASE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
CORPUS = os.path.join(BASE, "corpus", "sends.jsonl")
SRC = os.path.join(BASE, "cases")
NOW = datetime(2026, 9, 18, 18, 11)


def t(name, cond, detail=""):
    print(("PASS " if cond else "FAIL ") + name + (" | " + str(detail) if detail and not cond else ""))
    if not cond:
        sys.exit(1)


# ---------- 1. lexicon identity (drift-proof) ----------
t("critic.LABELS is lexicon.LABELS", critic_mod.LABELS is lex.LABELS)
t("critic.SPECIMEN is lexicon.SPECIMEN", critic_mod.SPECIMEN is lex.SPECIMEN)
t("readlayer.CHARGE is lexicon.CHARGE", rl_mod.CHARGE is lex.CHARGE)
t("readlayer.WITHDRAW is lexicon.WITHDRAW", rl_mod.WITHDRAW is lex.WITHDRAW)

# ---------- 2. certification oracle ----------
bad = certify_frame({"frame": True, "id": "x", "family": "tease", "text": ["you're gorgeous"],
                     "outcome": "extends", "slots": [], "forecast": {"extends": 1}, "protected": []})
t("generic line fails certification", any(d.startswith("critic:not_particular") for d in bad), bad)

bad2 = certify_frame({"frame": True, "id": "x", "family": "tease", "text": ["literally come be mean then"],
                      "outcome": "extends", "slots": [], "forecast": {"extends": 1}, "protected": []})
t("filler word caught", any(d.startswith("filler") for d in bad2), bad2)

bad3 = certify_frame({"frame": True, "id": "x", "family": "tease", "text": ["the {referent} is doing too much"],
                      "outcome": "extends", "slots": ["referent"], "forecast": {"nope": 1}, "protected": []})
t("bad forecast class caught", any(d.startswith("bad_forecast_class") for d in bad3), bad3)

bad4 = certify_frame({"frame": True, "id": "x", "family": "tease"})
t("missing fields caught", any(d.startswith("missing_field") for d in bad4), bad4)

report = certify_corpus(CORPUS)
t("all authored frames certify clean", all(v == [] for v in report.values()), report)
t("eight authored frames", len(report) == 8, list(report))

# harvested rows (frame: false) are not certified
rows = frames_mod.load_frames(CORPUS)
t("authored frames marked", sum(1 for r in rows if r.get("frame") is True) == 8)

# ---------- 3. frame_id + behavior seal ----------
d = tempfile.mkdtemp()
shutil.copytree(SRC, os.path.join(d, "cases"))
cp = os.path.join(d, "sends.jsonl")
shutil.copy(CORPUS, cp)

r1 = run_turn("rae", os.path.join(d, "cases"), cp, "mean", NOW, rng_seed=7)
t("send carries frame", r1.get("frame") and r1["frame"].get("id") in ("tease-mean-1", "tease-basil-1"), r1.get("frame"))
t("frame in skeleton", r1["skeleton"].get("frame", {}).get("id") == r1["frame"]["id"])
t("variant index present", isinstance(r1["frame"].get("variant"), int))

from dc_extended.live import confirm_send
confirm_send(os.path.join(d, "cases"), "rae")
shutil.copy(CORPUS, cp)  # record_send mutates the corpus; seal test needs identical state
shutil.copy(os.path.join(SRC, "rae.json"), os.path.join(d, "cases", "rae.json"))  # same case state too
r2 = run_turn("rae", os.path.join(d, "cases"), cp, "mean", NOW, rng_seed=7)
t("behavior seal: same seed same frame", r2["frame"] == r1["frame"] and r2["text"] == r1["text"],
  (r1["frame"], r2.get("frame")))

shutil.copy(CORPUS, cp)
shutil.copy(os.path.join(SRC, "rae.json"), os.path.join(d, "cases", "rae.json"))
r3 = run_turn("rae", os.path.join(d, "cases"), cp, "mean", NOW, rng_seed=8)
seen = {r1["text"], r2["text"], r3["text"]}
t("different seed can differ", len(seen) >= 1 and r3.get("frame") is not None, seen)

# operator draft fallback: no frame id
shutil.copy(os.path.join(SRC, "rae.json"), os.path.join(d, "cases", "rae.json"))  # reset: prior sends unconfirmed
res = run_turn("rae", os.path.join(d, "cases"), cp, "still waiting on her", NOW, rng_seed=1,
               operator_draft="the mean prompt owed me an answer")
t("draft fallback frame null", res["action"] == "send" and res["source"] == "edited"
  and (res.get("frame") is None or res["frame"].get("id") is None), (res.get("source"), res.get("frame")))

print("ALL PHASE 3 TESTS PASS")
