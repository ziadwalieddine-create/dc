"""Phase 0 acceptance: intake parser, unified deviation recovery, decision ports."""
import sys, os, json
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from dc_extended.intake import open_turn, bind_case, lock_claim_intake, trigger_has_her
from dc_extended.recovery import (record_deviation, forced_resolution, update_reciprocity,
                                  INTEGRITY_STOPS, EVENT_MUTATIONS)
from dc_extended.variety import is_synonym, derive_seed, weighted_sample

SRC = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "cases")


def t(name, cond, detail=""):
    print(("PASS " if cond else "FAIL ") + name + (" | " + str(detail) if detail and not cond else ""))
    if not cond:
        sys.exit(1)


# ---------- intake: open_turn ----------
c = open_turn("she: the basil plant is unhinged\nplay the bit")
t("her bubble parsed", c["her_last_bubble"] == "the basil plant is unhinged", c)
t("claim locked from instruction", c["claim"]["type"] == "bit", c["claim"])
t("newer bubble flagged", c.get("newer_bubble") is True)

c = open_turn("she: whatever, don't book thursday")
t("don't-book beats weekday", c["claim"]["type"] == "bit", c["claim"])

c = open_turn("play the bit")
t("no invented bubble", c.get("her_last_bubble") is None, c)

c = open_turn("she: come over\nme: thursday drinks")
t("his mark sets claim", c["claim"]["type"] == "plan", c["claim"])

c = open_turn("her last: haha ok")
t("her-last mark", c["her_last_bubble"] == "haha ok", c)

t("trigger_has_her", trigger_has_her("she said: hi") and not trigger_has_her("book it"))

# stored claim survives a null instruction turn
c = open_turn("", {"her_last_bubble": None, "claim": {"type": "bit", "source": "stored", "text": "x"},
                   "newer_bubble": False})
t("stored claim survives", c["claim"]["source"] == "stored")

# lock_claim_intake precedence table
t("stop -> close", lock_claim_intake("stop", None)["type"] == "close")
t("to ig -> channel", lock_claim_intake("move to ig", None)["type"] == "channel")
t("lock thursday -> plan", lock_claim_intake("lock thursday", None)["type"] == "plan")
t("null op + stored", lock_claim_intake(None, {"type": "plan", "text": "p"})["type"] == "plan")
t("null null -> null", lock_claim_intake(None, None)["type"] is None)

# ---------- intake: bind_case ----------
case = {"id": "rae"}
bind_case(case, ["rae", "jisu"])
t("bind ok", True)
try:
    bind_case(case, ["rae", "rae"])
    t("duplicate raises", False)
except ValueError:
    t("duplicate raises", True)
try:
    bind_case({"id": None}, ["rae"])
    t("missing id raises", False)
except ValueError:
    t("missing id raises", True)

# ---------- recovery: record_deviation ----------
case = {"referents": [{"id": "b", "text": "basil plant", "status": "live"}],
        "last_send": {"text": "the basil one", "confirmed": True},
        "amplitude": 1, "reciprocity": 0, "claim": {"type": "bit"}}

row = record_deviation(case, "joke_misunderstood", referent_id="b")
t("joke misunderstood closes referent", case["referents"][0]["status"] == "closed")
t("explain stays false", EVENT_MUTATIONS["joke_misunderstood"]["explain"] is False)

case2 = {"referents": [], "amplitude": 2, "reciprocity": 0, "claim": {"type": "bit"}}
record_deviation(case2, "too_forward")
t("too forward caps -1", case2["amp_cap"] == 1, case2["amp_cap"])

case3 = {"referents": [], "reciprocity": 0, "claim": {"type": "bit"},
         "last_send": {"text": "wyd on friday", "confirmed": True}, "banned": []}
row3 = record_deviation(case3, "operator_rejected", mechanism=None, family="wyd_plan", gate=None)
t("mechanism UNKNOWN stays null", row3["mechanism"] is None)
t("family banned on rejection", case3["banned"][-1]["family"] == "wyd_plan")

row4 = record_deviation(case3, "operator_rejected", mechanism="too_passive", family="plan", gate="too_passive")
t("named mechanism recorded", row4["mechanism"] == "too_passive")
t("named gate on ban", case3["banned"][-1]["gate"] == "too_passive")

try:
    record_deviation(case3, "not_a_real_event")
    t("unknown event raises", False)
except ValueError:
    t("unknown event raises", True)

# ---------- recovery: forced_resolution ----------
fr_case = {"referents": [{"id": "b", "status": "live"}], "claim": {"type": "bit"},
           "unique_facts": ["basil plant"]}
t("integrity stops never degrade", all(forced_resolution(fr_case, w) is None for w in INTEGRITY_STOPS))
r = forced_resolution(fr_case, "unconfirmed_send")
t("degraded hold -> amp0 Play", r["policy"] == "Play" and r["amplitude"] == 0, r)
fr_case["referents"] = []
r = forced_resolution(fr_case, "dead_clock")
t("no referent -> Trade on fact", r["policy"] == "Trade", r)
fr_case["unique_facts"] = []
t("nothing legal -> None", forced_resolution(fr_case, "dead_clock") is None)
fr_case2 = {"referents": [{"id": "b", "status": "live"}], "refused_object": True,
            "claim": {"type": "bit"}, "unique_facts": []}
t("refused object blocks forced Play", forced_resolution(fr_case2, "dead_clock") is None)

# ---------- recovery: reciprocity ----------
rc = {"last_send": {"text": "line", "confirmed": True}, "reciprocity": 0, "outcomes": []}
t("extends +1", update_reciprocity(rc, "extends") == 1)
update_reciprocity(rc, "refuses")
t("refuses -2", rc["reciprocity"] == -1, rc["reciprocity"])
rc["reciprocity"] = 3
t("clamp at 3", update_reciprocity(rc, "extends") == 3)
rc2 = {"last_send": {"text": "l", "confirmed": True}, "reciprocity": 0, "outcomes": []}
update_reciprocity(rc2, "dry")
update_reciprocity(rc2, "dry")
t("two dries cap amp", rc2["amp_cap"] == 0)
try:
    update_reciprocity({"last_send": {"text": "l", "confirmed": False}, "reciprocity": 0, "outcomes": []}, "extends")
    t("receipt required", False)
except ValueError:
    t("receipt required", True)
try:
    update_reciprocity(rc2, "lol")
    t("unknown class raises", False)
except ValueError:
    t("unknown class raises", True)

# ---------- decision: is_synonym ----------
t("trailing adjective", is_synonym("the basil plant is a lot", "the basil plant is a lot tall"))
t("not synonyms", not is_synonym("basil plant died", "thursday at eight"))
t("same length no", not is_synonym("come over", "come over here now"))

# ---------- schema: template fields ----------
tmpl = json.load(open(os.path.join(SRC, "_template.json")))
for field in ("endpoint", "reciprocity", "amp_cap", "deviations"):
    t("template has " + field, field in tmpl)

print("ALL PHASE 0 TESTS PASS")
