"""Skeleton tests — every field populated, candidates logged, select reasoned."""
import sys, os, shutil, tempfile, json
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from datetime import datetime
from dc_extended.switchboard import run_turn

SRC = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "cases")
FRI_EVEN = datetime(2026, 9, 18, 18, 11)


def fixture():
    d = tempfile.mkdtemp()
    shutil.copytree(SRC, os.path.join(d, "cases"))
    corpus = os.path.join(d, "sends.jsonl")
    shutil.copy(os.path.join(os.path.dirname(SRC), "corpus", "sends.jsonl"), corpus)
    return d, corpus


def t(name, cond, detail=""):
    print(("PASS " if cond else "FAIL ") + name + (" | " + str(detail) if detail and not cond else ""))
    if not cond:
        sys.exit(1)


# 1. Legal tease returns a skeleton with all fields
d, cp = fixture()
r = run_turn("rae", os.path.join(d, "cases"), cp, "mean", FRI_EVEN)
t("send returns skeleton", "skeleton" in r, r)
skel = r["skeleton"]
t("skeleton has case", "case" in skel and skel["case"]["name"] == "Rae")
t("skeleton has read", "read" in skel and skel["read"]["kind"] == "none")
t("skeleton has claim", "claim" in skel and skel["claim"]["type"] == "bit")
t("skeleton has gates", "gates" in skel and len(skel["gates"]) >= 10)
t("all gates PASS", all(g["result"] == "PASS" for g in skel["gates"]))
t("skeleton has candidates", "candidates" in skel and len(skel["candidates"]) >= 1)
t("skeleton has select", "select" in skel and skel["select"]["winner"] is not None)
t("skeleton has surface", "surface" in skel and skel["surface"]["protected_tokens"] == "PASS")
t("skeleton has release", "release" in skel and skel["release"] == r["text"])

# 2. Candidate fields are complete
c0 = skel["candidates"][0]
t("candidate has text", bool(c0.get("text")))
t("candidate has source", c0.get("source") in ("retrieved", "edited"))
t("candidate has family_rate", isinstance(c0.get("family_rate"), (int, float)))
t("candidate has heat_fit", isinstance(c0.get("heat_fit"), (int, float)))
t("candidate has amp_fit", isinstance(c0.get("amp_fit"), (int, float)))
t("candidate has critic", "critic" in c0)
t("candidate has simulate", isinstance(c0.get("simulate"), list))

# 3. Gate failure returns skeleton with FAIL gates
d, cp = fixture()
r = run_turn("camille", os.path.join(d, "cases"), cp, "wyd on friday", FRI_EVEN)
t("stop returns skeleton", "skeleton" in r, r)
skel = r["skeleton"]
t("gate FAIL logged", any(g["result"] == "FAIL" for g in skel["gates"]))
t("no candidates on gate stop", len(skel["candidates"]) == 0)
t("no select on gate stop", skel["select"] is None)
t("no release on gate stop", skel["release"] is None)

# 4. Refusal returns skeleton with read logged
d, cp = fixture()
r = run_turn("rae", os.path.join(d, "cases"), cp, "mean", FRI_EVEN,
             her_bubble="no, I'm not interested")
t("refusal returns skeleton", "skeleton" in r, r)
t("refusal read logged", r["skeleton"]["read"]["kind"] == "refuse")

# 5. Skeleton stored on case chronology
d, cp = fixture()
r = run_turn("rae", os.path.join(d, "cases"), cp, "mean", FRI_EVEN)
case = json.load(open(os.path.join(d, "cases", "rae.json")))
t("last_send has source", case["last_send"].get("source") in ("retrieved", "edited"))
t("chronology logged", len(case["chronology"]) >= 4)

print("ALL SKELETON TESTS PASS")
