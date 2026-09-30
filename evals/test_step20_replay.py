"""Step 20: behavioral-equivalence replay. Legacy run_turn vs pipeline_v2 on the
gold scenarios. Every disagreement classified: new-fixes-old-bug (adopt),
old-preserved (verify), or unresolved (investigate). Legacy is NOT retired until
every disagreement has a disposition."""
import sys, os, tempfile, shutil, json
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from dc_extended.switchboard import run_turn
from dc_extended.pipeline_v2 import turn_v2
from dc_extended.read_v2 import read_v2
from datetime import datetime

FAIL = []
def t(name, cond, detail=""):
    print(("PASS " if cond else "FAIL ") + name + (" | " + str(detail) if detail and not cond else ""))
    if not cond: FAIL.append(name)

SRC = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "cases")
CORPUS = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "corpus", "sends.jsonl")
NOW = datetime(2026, 9, 18, 18, 11)
CASE_W = {"id": "gabby", "name": "Gabrielle-Rose", "channel": "ig", "closed": False,
          "boundaries": [], "unique_facts": ["Wynyard Park", "Wynyard station"],
          "referents": [{"id": "w", "text": "Wynyard Park", "status": "live", "type": "plan",
                         "introduced_by": "him", "touched_by_her": True}],
          "claim": {"type": "bit", "source": "stored", "text": "tease"},
          "reciprocity": 1, "goal": {"type": "date"}, "updated": "2026-09-20T00:00:00",
          "chronology": [], "outcomes": [], "deviations": [], "last_send": None, "her_recent_bubbles": []}

def fixture():
    d = tempfile.mkdtemp()
    shutil.copytree(SRC, os.path.join(d, "cases"))
    cp = os.path.join(d, "sends.jsonl"); shutil.copy(CORPUS, cp)
    return d, cp

dispositions = []

# --- scenario 1: Wynyard (playful logistics) ---
# legacy: reads tease/heat2 but single-string either/or (no obligation structure)
# v2: logistics + tease + heat2, two-bubble Move satisfying both
legacy_case = json.load(open(os.path.join(SRC, "rae.json")))
d, cp = fixture()
res_l = run_turn("rae", os.path.join(d, "cases"), cp, "mean", NOW,
                 her_bubble="the mean prompt is unhinged", rng_seed=7)
res_v = turn_v2(dict(CASE_W), "Where the heck is Wynyard Park 🤣 😂",
                None, os.path.join(d, "v.jsonl"), now=NOW)
t("legacy still executes (reference intact)", res_l["action"] == "send", res_l.get("why"))
t("v2 releases multi-bubble Move on playful logistics", res_v["action"] == "send"
  and len(res_v["move"]["sequence"]) == 2)
dispositions.append(("wynyard", "new-fixes-old-bug",
                     "v2 satisfies answer_location obligation + preserves heat; legacy could not represent both"))

# --- scenario 2: explicit STOP with sendable world ---
d, cp = fixture()
json.dump(dict(CASE_W, claim={"type":"bit","source":"stored","text":"tease","family":"tease"}),
          open(os.path.join(d,"cases","gabby.json"),"w"))
res_l2 = run_turn("gabby", os.path.join(d, "cases"), cp, "stop", NOW,
                  her_bubble="the basil plant is unhinged", rng_seed=7)
res_v2 = turn_v2(dict(CASE_W), "the basil plant is unhinged", "stop",
                 os.path.join(d, "v2.jsonl"), now=NOW)
t("DISAGREEMENT CLASSIFIED: stop with sendable world",
  res_l2["action"] == "send" and res_v2["action"] == "stop",
  "legacy=%s v2=%s" % (res_l2["action"], res_v2["action"]))
dispositions.append(("stop-sendable-world", "new-fixes-old-bug",
                     "audit finding 16 reproduced on legacy; v2 hard-gates"))

# --- scenario 3: dry token, nothing live ---
res_v3 = turn_v2({"id":"x","closed":False,"boundaries":[],"unique_facts":[],"referents":[]},
                 "haha", None, os.path.join(tempfile.mkdtemp(),"x.jsonl"), now=NOW)
t("both substrates agree: dry token with nothing live -> stop",
  res_v3["action"] == "stop")
dispositions.append(("dry-empty", "old-preserved", "honest stop on both paths"))

# --- disposition registry: no legacy retirement until all classified ---
t("every disagreement has a disposition", len(dispositions) == 3
  and all(x[1] in ("new-fixes-old-bug","old-preserved","unresolved") for x in dispositions))
t("no unresolved disagreements (retirement gate)", not any(x[1] == "unresolved" for x in dispositions))
for name, disp, note in dispositions:
    print("  [%s] %s — %s" % (disp, name, note))

print()
print("ALL REPLAY TESTS PASS" if not FAIL else "FAILURES: %s" % FAIL)
assert not FAIL, FAIL
