"""Acceptance suite — every named incident from the dc ledger must be caught."""
import sys, os, shutil, tempfile, json
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from datetime import datetime
from dc_extended.switchboard import run_turn

SRC = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "cases")
FRI_EVEN = datetime(2026,9,18,18,11)

def fixture():
    d = tempfile.mkdtemp()
    shutil.copytree(SRC, os.path.join(d,"cases"))
    corpus = os.path.join(d,"sends.jsonl")
    shutil.copy(os.path.join(os.path.dirname(SRC),"corpus","sends.jsonl"), corpus)
    return d, corpus

def t(name, cond, detail=""):
    print(("PASS " if cond else "FAIL ") + name + (" | " + str(detail) if detail and not cond else ""))
    if not cond: sys.exit(1)

# 1. RAE FRIDAY — the line that sat unsent into Friday evening
d, cp = fixture()
r = run_turn("rae", os.path.join(d,"cases"), cp, "wyd on friday, i'll pick the place", FRI_EVEN)
t("rae friday -> dead_clock stop", r["action"]=="stop" and r.get("gate") in ("dead_clock","restack"), r)

# 2. JISU GOAL SMUGGLE — booking Wednesday over a stated get-to-know want
d, cp = fixture()
r = run_turn("jisu", os.path.join(d,"cases"), cp, "book dinner, my pick", datetime(2026,9,19,12,0))
t("jisu booking -> goal_smuggle stop", r["action"]=="stop" and r.get("gate")=="goal_smuggle", r)

# 3. MARY PLATFORM — Hinge opener while thread lives on IG
d, cp = fixture()
mp = os.path.join(d,"cases","mary.json")
mc = json.load(open(mp)); mc["updated"]="2026-09-19T10:00:00"
json.dump(mc, open(mp,"w"))
r = run_turn("mary", os.path.join(d,"cases"), cp, "fresh opener on hinge", datetime(2026,9,19,12,0))
t("mary hinge -> wrong_platform stop", r["action"]=="stop" and r.get("gate")=="wrong_platform", r)

# 4. GHOST FACE — draft carrying another woman's fact
d, cp = fixture()
r = run_turn("rae", os.path.join(d,"cases"), cp, "tease her: jisu said the walk was nice", FRI_EVEN,
                 operator_draft="jisu said the walk was nice, tease her about it")
t("cross-thread -> stop", r["action"]=="stop" and r.get("gate") in ("cross_thread","wrong_person"), r)

# 5. CAMILLE — closed thread generates nothing
d, cp = fixture()
r = run_turn("camille", os.path.join(d,"cases"), cp, "wyd on friday", FRI_EVEN)
t("closed thread -> stale_goal stop", r["action"]=="stop" and r.get("gate")=="stale_goal", r)

# 6. UNCONFIRMED — building on a send that never left the phone
d, cp = fixture()
casep = os.path.join(d,"cases","rae.json")
c = json.load(open(casep)); c["last_send"]={"text":"the dm","ts":"2026-09-18T20:00:00","confirmed":False,"family":"dm"}
json.dump(c, open(casep,"w"))
r = run_turn("rae", os.path.join(d,"cases"), cp, "follow up on the dm", datetime(2026,9,19,9,0))
t("unconfirmed build -> stop", r["action"]=="stop" and r.get("gate")=="unconfirmed_send", r)

# 7. CEREMONY — an audit request is not an action
d, cp = fixture()
r = run_turn("rae", os.path.join(d,"cases"), cp, "audit all my cases and compare them", FRI_EVEN)
t("audit -> paste_without_send stop", r["action"]=="stop" and r.get("why")=="paste_without_send", r)

# 8. STALE STATE — >48h old case refuses to generate
d, cp = fixture()
r = run_turn("mary", os.path.join(d,"cases"), cp, "tease about the armored truck", datetime(2026,9,19,12,0))
t("stale case -> reverify stop", r["action"]=="stop" and r.get("why")=="reverify", r)

# 9. LEGAL TEASE — same case, live referent, in-window: ship, and write-forward happened
d, cp = fixture()
r = run_turn("rae", os.path.join(d,"cases"), cp, "mean", FRI_EVEN)
t("legal tease -> send", r["action"]=="send", r)
c = json.load(open(os.path.join(d,"cases","rae.json")))
t("write-forward atomic", c["last_send"]["text"]==r["text"] and c["last_send"]["confirmed"]==False and len(c["chronology"])>=4)

print("ALL INCIDENT REPLAYS PASS")
