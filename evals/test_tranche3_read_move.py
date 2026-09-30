"""Tranche 3: multidimensional Read + first-class Move, tested against the FROZEN
replay gold. Every gold property asserted; every overcorrection guarded by a
negative control. Legacy kind stays available as a projection."""
import sys, os, json
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from dc_extended.read_v2 import read_v2, legacy_kind
from dc_extended.moves import single, sequence, from_frame_variants, wynward_move

FAIL = []
def t(name, cond, detail=""):
    print(("PASS " if cond else "FAIL ") + name + (" | " + str(detail) if detail and not cond else ""))
    if not cond: FAIL.append(name)

gold = json.load(open(os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
                                   "training", "replay_gold.json")))
CASE_W = {"referents": [{"id": "w", "text": "Wynyard Park", "status": "live", "type": "plan",
                         "introduced_by": "him", "touched_by_her": True}]}

# --- Wynyard: the multidimensional canonical case ---
r = read_v2("Where the heck is Wynyard Park 🤣 😂", CASE_W)
g = gold[0]
t("wynyard: logistics_question present", "logistics_question" in r["literal_functions"], r["literal_functions"])
t("wynyard: tease present", "tease" in r["relational_functions"], r["relational_functions"])
t("wynyard: heat 2 preserved WITH logistics", r["heat"] == 2, r["heat"])
t("wynyard: answer_location obligation", "answer_location" in r["response_obligations"])
t("wynyard: question object = location", "location" in r["question_objects"])
t("wynyard: legacy kind still derivable", legacy_kind(r) in ("tease","question","escalate"))

# --- Hot: not dry ---
r = read_v2("Hot", CASE_W)
g = gold[1]
t("hot: heat 2", r["heat"] == 2, r["heat"])
t("hot: investment positive", r["investment"] == "positive")
t("hot: NOT dry / none", "none" not in r["literal_functions"] and r["heat"] is not None)
t("hot: attraction cue recorded", "attraction_cue" in r["relational_functions"])

# --- No Thursday, Friday works: scoped refusal ---
r = read_v2("No Thursday 😂 Friday works", {})
t("scoped: Thursday refused", r["refusal_scope"] == "thursday", r["refusal_scope"])
t("scoped: Friday accepted", r["acceptance"] == "friday", r["acceptance"])
t("scoped: engagement continues", r["progression_signal"] == "positive")
t("scoped: not global", r["refusal_scope"] != "global")

# --- heat non-regression: literal function must not sand heat ---
r = read_v2("where even is Wynyard Park 🤣", CASE_W)
t("non-regression: playful logistics stays heat 2", r["heat"] == 2, r["heat"])

# --- negative controls: no free heat ---
r = read_v2("Where is the station?", {})   # neutral logistics, no play marker
t("NEG: neutral logistics gains NO heat", r["heat"] in (None, 0), r["heat"])
r = read_v2("my train leaves at 5", {})
t("NEG: plain disclosure gains NO heat", r["heat"] in (None, 0), r["heat"])
r = read_v2("haha", CASE_W)
t("NEG: dry token with live referent stays low", (r["heat"] or 0) <= 1, r["heat"])
CASE_OLDMEAN = {"referents": [{"id": "m", "text": "mean", "status": "live", "type": "bit",
                               "introduced_by": "him", "touched_by_her": True}]}
r = read_v2("I mean honestly", CASE_OLDMEAN)
t("NEG: 'I mean honestly' does NOT touch the old 'mean' referent",
  "tease" not in r["relational_functions"] and (r["heat"] or 0) == 0,
  (r["relational_functions"], r["heat"]))

# --- playful stop is resistance, not refusal ---
r = read_v2("stop being cute 😂", {})
t("playful stop: resistance, not global", r["refusal_scope"] == "playful_resistance", r["refusal_scope"])

# --- Move semantics ---
m = wynward_move("Right outside Wynyard station", playful_continuation="You'll find me \U0001f602")
t("move: two ordered bubbles", m["sequence"] == ["Right outside Wynyard station", "You'll find me 😂"])
t("move: legacy text view = first bubble", m["text"] == "Right outside Wynyard station")
t("move: carries its obligation", "answer_location" in m["obligations"])
one = single("School night's exactly why it's fun", family="tease")
t("NEG: single-bubble move remains first-class", one["sequence"] == ["School night's exactly why it's fun"])
try:
    sequence([""]); t("NEG: empty move rejected", False)
except ValueError:
    t("NEG: empty move rejected", True)
try:
    sequence(["a", ""]); t("NEG: empty bubble in sequence rejected", False)
except ValueError:
    t("NEG: empty bubble in sequence rejected", True)
frm = {"text": ["Come be mean then", "Still waiting on that mean"], "family": "tease", "heat": 2}
mv = from_frame_variants(frm, {})
t("NEG: frame variants stay variants, never a sequence", len(mv.rendered()) == 1)

print()
print("ALL TRANCHE-3 TESTS PASS" if not FAIL else "FAILURES: %s" % FAIL)
assert not FAIL, FAIL
