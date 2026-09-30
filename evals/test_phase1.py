"""Phase 1 acceptance: split confidence, window reads, deep contract, grading."""
import sys, os, json, tempfile, shutil
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from dc_extended.readlayer import (read_bubble_v2, needs_deep, push_window,
                                   resolve_split, grade_pending_read, WINDOW_PATTERNS)
from dc_extended.switchboard import run_turn
from datetime import datetime

SRC = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "cases")
NOW = datetime(2026, 9, 18, 18, 11)
BASIL = {"id": "b", "text": "basil plant", "type": "bit",
         "introduced_by": "him", "touched_by_her": True, "status": "live"}


def t(name, cond, detail=""):
    print(("PASS " if cond else "FAIL ") + name + (" | " + str(detail) if detail and not cond else ""))
    if not cond:
        sys.exit(1)


def case_with_bit(**kw):
    c = {"id": "t", "name": "T", "stage": "chat", "channel": "hinge",
         "goal": {"type": "date"}, "updated": "2026-09-19T10:00:00",
         "banned_families": [], "boundaries": [], "unique_facts": ["basil plant"],
         "missed_windows": [], "referents": [dict(BASIL)], "her_recent_bubbles": []}
    c.update(kw)
    return c


# ---------- split detection (A4 competing reads) ----------
r = read_bubble_v2("busy this week but prove it", case_with_bit())
t("charge+withdraw -> split", r["confidence"] == "split", r)
t("alternates both present", set(r["alternates"]) == {"escalate", "withdraw"}, r)
t("challenge wins (read.md)", r["kind"] == "escalate", r)

r = read_bubble_v2("basil plant is fine I guess, been busy", case_with_bit())
t("tease+withdraw -> split", r["confidence"] == "split" and r["kind"] == "tease", r)

r = read_bubble_v2("what are you up to? not sure I'm free", case_with_bit())
t("question+withdraw -> split", r["confidence"] == "split" and "question" in r["alternates"], r)

r = read_bubble_v2("the basil plant is unhinged", case_with_bit())
t("clean tease stays high", r["confidence"] == "high" and r["kind"] == "tease", r)

r = read_bubble_v2("no, stop it", case_with_bit())
t("refusal still beats tease", r["kind"] == "refuse" and r["confidence"] == "high", r)

r = read_bubble_v2("come over tonight", case_with_bit())
t("clean charge stays high", r["kind"] == "escalate" and r["confidence"] == "high", r)

# ---------- window reads ----------
c = case_with_bit(her_recent_bubbles=["haha", "haha"])
push_window(c, "what are you doing?")
r = read_bubble_v2("what are you doing?", c)
t("dry,dry,question -> testing", r.get("window_signal") == "testing_investment", r)

c = case_with_bit(her_recent_bubbles=["the basil plant though", "haha"])
push_window(c, "been busy lately")
r = read_bubble_v2("been busy lately", c)
t("tease,dry,withdraw -> push_pull", r.get("window_signal") == "push_pull", r)

c = case_with_bit()
for b in ["a", "b", "c", "d"]:
    push_window(c, b)
t("window trims to 3", len(c["her_recent_bubbles"]) == 3, c["her_recent_bubbles"])
push_window(c, "d")
t("consecutive dedupe", c["her_recent_bubbles"] == ["b", "c", "d"], c["her_recent_bubbles"])

# ---------- needs_deep ----------
t("split needs deep", needs_deep({"confidence": "split"}))
t("high does not", not needs_deep({"confidence": "high"}))

# ---------- run_turn deep contract ----------
d = tempfile.mkdtemp()
shutil.copytree(SRC, os.path.join(d, "cases"))
cp = os.path.join(d, "sends.jsonl")
shutil.copy(os.path.join(os.path.dirname(SRC), "corpus", "sends.jsonl"), cp)

res = run_turn("rae", os.path.join(d, "cases"), cp, "mean", NOW,
               her_bubble="busy this week but prove it")
t("split -> deep action", res["action"] == "deep", res)
t("deep why", res["why"] == "split_read", res)
t("deep emits no line", res.get("text") is None, res)
t("skeleton present", "skeleton" in res and res["skeleton"]["read"]["confidence"] == "split")

rae = json.load(open(os.path.join(d, "cases", "rae.json")))
t("pending_read stored on case", rae.get("pending_read", {}).get("primary") == "escalate", rae.get("pending_read"))
t("window pushed", "busy this week but prove it" in (rae.get("her_recent_bubbles") or []), rae.get("her_recent_bubbles"))

# resolve then force_fast proceeds
from dc_extended.readlayer import resolve_split
resolve_split(rae, "escalate")
json.dump(rae, open(os.path.join(d, "cases", "rae.json"), "w"))
res2 = run_turn("rae", os.path.join(d, "cases"), cp, "mean", NOW,
                her_bubble="busy this week but prove it", force_fast=True)
t("resolved + force_fast sends", res2["action"] == "send" and res2["text"], res2)

# resolve_split guards
try:
    resolve_split({"pending_read": {"alternates": ["escalate", "withdraw"], "primary": "escalate"}}, "dry")
    t("resolve rejects unknown kind", False)
except ValueError:
    t("resolve rejects unknown kind", True)
try:
    resolve_split({}, "escalate")
    t("resolve requires pending", False)
except ValueError:
    t("resolve requires pending", True)

# ---------- grading (affordance-as-probe) ----------
stats = os.path.join(d, "stats.json")
c = {"pending_read": {"signal": "busy but prove it", "alternates": ["escalate", "withdraw"], "primary": "escalate"}}
w = grade_pending_read(c, "extends", stats_path=stats)
t("extends -> primary wins", w == "escalate", w)
s = json.load(open(stats))
t("stats counted", s["escalate"] == {"n": 1, "correct": 1}, s)

c = {"pending_read": {"signal": "busy but prove it", "alternates": ["escalate", "withdraw"], "primary": "escalate"}}
w = grade_pending_read(c, "dry", stats_path=stats)
t("dry -> secondary wins", w == "withdraw", w)
s = json.load(open(stats))
t("primary marked wrong", s["escalate"] == {"n": 2, "correct": 1}, s)

c = {"pending_read": {"signal": "busy but prove it", "alternates": ["escalate", "withdraw"], "primary": "escalate"}}
w = grade_pending_read(c, "answers", stats_path=stats)
t("answers -> no count", w is None and "pending_read" not in c, w)
s = json.load(open(stats))
t("stats untouched on no-count", s["escalate"]["n"] == 2, s)

t("no pending -> None", grade_pending_read({}, "extends", stats_path=stats) is None)

print("ALL PHASE 1 TESTS PASS")
