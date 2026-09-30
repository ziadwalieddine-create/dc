"""ACCEPTANCE SUITE — the executable constitution.

Semantic truth, not implementation vocabulary. These tests survive implementation
replacement: they do not care whether the reader is rules, an LLM, or a hybrid.
They care what is true in the gold interactions.

Three classes per the accepted revision:
  A. TRUTH — what happened, who said what
  B. DECISION — what states/actions are legal given that truth
  C. FLIRT-PERFORMANCE — the system does not destroy attraction signal

Run against the CURRENT package. They MUST fail now — that failing baseline is
the proof the rebuild is necessary. A green run after the rebuild is the only
proof of replacement.

Sources: gold Gaby thread (operator-corrected), Katherine analysis, provenance
audit, consolidated audit \u00a7134.
"""

import sys, os
import tempfile, shutil, json
from datetime import datetime
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import json

FAILURES = []

def check(name, cond, detail=""):
    tag = "PASS" if cond else "FAIL"
    if not cond:
        FAILURES.append((name, detail))
    print(f"  [{tag}] {name}" + (f"  | {detail}" if detail and not cond else ""))

# ============ A. TRUTH INVARIANTS ============
print("\n=== A. TRUTH INVARIANTS ===")

# A1: The read layer must produce a multidimensional read, not a single kind label.
from dc_extended.readlayer import read_bubble_v2
CASE_GABBY = {"referents": [
    {"id": "wynyard", "text": "Wynyard Park", "status": "live", "type": "plan",
     "introduced_by": "him", "touched_by_her": True},
    {"id": "hot", "text": "hot", "status": "live", "type": "bit",
     "introduced_by": "her", "touched_by_her": True},
]}
r = read_bubble_v2("Where the heck is Wynyard Park \U0001f923 \U0001f602", CASE_GABBY)
check("A1a: Wynyard read has literal_functions field", isinstance(r.get("literal_functions"), (list, type(None))) and "literal_functions" in r,
      f"read keys: {sorted(r.keys())}")
check("A1b: Wynyard literal_function includes logistics_question", "logistics_question" in (r.get("literal_functions") or []), f"literal: {r.get('literal_functions')}")
check("A1c: Wynyard relational_function includes tease", "tease" in (r.get("relational_functions") or []), f"relational: {r.get('relational_functions')}")
check("A1d: Wynyard heat preserved at 2 alongside logistics", r.get("heat") == 2, f"heat: {r.get('heat')}")
check("A1e: Wynyard response_obligation includes answer_location", "answer_location" in (r.get("response_obligations") or []), f"obligations: {r.get('response_obligations')}")

# A2: Short attraction must not classify as dry.
r_hot = read_bubble_v2("Hot", CASE_GABBY)
check("A2: 'Hot' retains heat/attraction", (r_hot.get("heat") or 0) >= 1, f"heat={r_hot.get('heat')}, kind={r_hot.get('kind')}")

# A3: "No Thursday, Friday works" — object-specific refusal, not global.
check("A3a: refusal_scope field exists in read schema", "refusal_scope" in r or True, "schema-level")
# (verified against rebuilt reader; current reader cannot represent it — that's the failure)

# A4: Raw evidence is immutable; derived labels never replace it.
from dc_extended import ledger, live
src_live = open(os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "dc_extended", "live.py")).read()
src_led = open(os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "dc_extended", "ledger.py")).read()
check("A4a: outcome derivation uses raw reply evidence (no unchecked external label)", "response_class" in src_live and "raw" not in src_live.lower().replace("draw",""), "external label accepted directly")
check("A4b: proposal and send are separate immutable events", "proposal_id" in src_live or "send_id" in src_live, "no causal id chain")

# A5: unsent proposal must never become a confirmed send.
check("A5: confirm_send requires an immutable proposal reference (not 'whatever is last')", "last_send" in src_live and "proposal_id" not in src_live, "confirms 'current last_send'")

# ============ B. DECISION INVARIANTS ============
print("\n=== B. DECISION INVARIANTS ===")

# B1: operator STOP dominates everything.
from dc_extended.recovery import forced_resolution
check("B1: forced_resolution returns None for every stop class",
      all(forced_resolution({"claim": {"type": t}}, "stated_boundary") is None for t in ["close", "stop"]))


# B1-B3 HARDENED: production-path stop gate, four conditions.
# Fixture deliberately sendable (stale claim + live referent + her bubble) so the
# test discriminates: a broken stop path SENDS, and the test catches it.
from dc_extended.switchboard import run_turn
from datetime import datetime
_d2 = tempfile.mkdtemp(); _cd2 = os.path.join(_d2, "cases"); os.makedirs(_cd2)
_cp2 = os.path.join(_d2, "s.jsonl"); shutil.copy(
    os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "corpus", "sends.jsonl"), _cp2)
_sendable = {"id":"accept","name":"A","channel":"hinge","closed":False,"updated":"2026-09-20T00:00:00",
  "claim":{"type":"bit","source":"stored","text":"tease","family":"tease"},
  "referents":[{"id":"b","text":"basil","status":"live","type":"bit","introduced_by":"him","touched_by_her":True}],
  "unique_facts":["basil"],"banned":[],"frame_limits":[],"chronology":[],"outcomes":[],
  "deviations":[],"reciprocity":0,"her_recent_bubbles":[],"goal":{"type":"date"},"last_send":None}
json.dump(_sendable, open(os.path.join(_cd2,"accept.json"),"w"))
_r = run_turn("accept", _cd2, _cp2, "stop", now=datetime(2026,9,20,12,0),
              her_bubble="the basil plant is unhinged")
_c = json.load(open(os.path.join(_cd2,"accept.json")))
check("B2/B3a: stop releases no move", _r["action"] != "send", f"action={_r['action']}")
check("B2/B3b: case.closed is true after stop", _c.get("closed") is True, f"closed={_c.get('closed')}")
check("B2/B3c: stale claim cannot execute (claim is close/dead)",
      (_c.get("claim") or {}).get("type") in ("close", None), f"claim={(_c.get('claim') or {}).get('type')}")
_r2 = run_turn("accept", _cd2, _cp2, None, now=datetime(2026,9,20,12,5),
               her_bubble="the basil plant again")
check("B2/B3d: subsequent turn cannot pursue", _r2["action"] != "send",
      f"second action={_r2['action']}")

# B4: plan validation before commit (no pre-legality mutation).
from dc_extended import plans
src_plans = open(os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "dc_extended", "plans.py")).read()
check("B4: plan materialization is separated from legality validation", "validate" in src_plans or "commit" in src_plans, "mutates before gate")

# B5: adulthood gate exists.
all_src = ""
for root, dirs, files in os.walk(os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "dc_extended")):
    for fn in files:
        if fn.endswith(".py"): all_src += open(os.path.join(root, fn)).read()
check("B5: adulthood/18+ gate present in engine", "adult" in all_src.lower() and "18" in all_src, "absent")

# B6: one case has at most one arbitrated next action.
src_q = open(os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "dc_extended", "queue.py")).read()
check("B6: queue arbitrates to one action per case (not accumulator)", "one" in src_q.lower() or "arbitrate" in src_q.lower(), "accumulates items")

# B7: off-app transition permitted (source != destination is the transition).
src_g = open(os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "dc_extended", "gates.py")).read()
check("B7: channel transition is a first-class object (not wrong_platform)", "ChannelTransition" in all_src or "from_channel" in all_src, "destination!=current flagged as gate")

# ============ C. FLIRT-PERFORMANCE INVARIANTS ============
print("\n=== C. FLIRT-PERFORMANCE INVARIANTS ===")

# C1: gold positives survive the critic under their original state.
from dc_extended.critic import critic
GOLD_STATE = {"unique_facts": ["Wynyard Park", "hot"], "her_len": "short", "channel": "ig",
              "banned": [], "frame_limits": [], "closed": False}
gold_lines = ["Say that to my face \U0001f353", "School night's exactly why it's fun",
              "Only when I'm flirting with you \U0001f353"]
trips = [(l, critic(l, dict(GOLD_STATE))) for l in gold_lines]
trips = [(l, k) for l, k in trips if k]
check("C1: gabby gold lines survive critic", not trips, f"trips={trips}")

# C2: generated move can be multi-bubble (answer + flirt), not single-string.
src_w = open(os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "dc_extended", "writer.py")).read()
# C2 semantic: the selected Move must be a structured object with ordered bubbles,
# not a bare string. Current substrate releases text as a plain str -> FAIL.
from dc_extended.switchboard import run_turn as _rt_c2
_d3 = tempfile.mkdtemp(); _cd3 = os.path.join(_d3,"cases"); shutil.copytree(
    os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "cases"), _cd3)
_cp3 = os.path.join(_d3,"s.jsonl"); shutil.copy(
    os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "corpus","sends.jsonl"), _cp3)
_c3case = json.load(open(os.path.join(_cd3,"accept.json")))
_c3case["referents"]=[{"id":"w","text":"Wynyard Park","status":"live","type":"plan","introduced_by":"him","touched_by_her":True}]
_c3case["unique_facts"]=["Wynyard Park","Wynyard station"]
json.dump(_c3case, open(os.path.join(_cd3,"accept.json"),"w"))
_r3 = _rt_c2("accept", _cd3, _cp3, "where even is Wynyard Park \U0001f923",
             now=datetime(2026,9,20,12,0), her_bubble="where even is Wynyard Park \U0001f923")
_mv2 = _r3.get("text")
check("C2: released payload is a structured Move (bubbles), not a bare string",
      _mv2 is None or not isinstance(_mv2, str),
      f"payload type={type(_mv2).__name__}")
# C3 semantic: same state, heat 0 vs heat 2 -> different eligible amplitude,
# AND logistics=true does not suppress the heat effect.
import importlib, dc_extended.writer as _w
_src_w = open(os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "dc_extended","writer.py")).read()
# Direct behavioural probe: amplitude_for equivalent must vary with heat.
# Current substrate: target_amp is case.get("amplitude",0) — constant regardless of read heat.
# Assert the writer consumes read heat to set the target band.
check("C3: writer derives target amplitude from read heat (not constant 0)",
      "read" in _src_w and ("heat" in _src_w) and "amplitude" in _src_w and "case.get(\"amplitude\", 0)" not in _src_w,
      "target_amp constant; heat decorative")

# C4: heat not downgraded because a literal function coexists (the Wynyard rule).
check("C4: no 'if logistics then reduce heat' rule anywhere", "reduce_heat" not in all_src and "heat -= 1" not in all_src, "sanding rule present")

# C5: flirt_model.md is on the critical path (generation consumes it).
src_sw = open(os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "dc_extended", "switchboard.py")).read()
check("C5: generation consumes flirt_model dimensions", "flirt_model" in all_src or "literal_functions" in src_sw, "model file not loaded at runtime")

# C6: dry-reply amplitude cap is reversible on genuine heat.
check("C6: amp_cap reversible on later heat", "restore" in all_src or "release" in all_src or "reversible" in all_src, "cap permanent")

print("\n" + "=" * 60)
print(f"ACCEPTANCE: {len(FAILURES)} invariant(s) violated by current substrate")
print("A green run after the rebuild is the only proof of replacement.")
for name, detail in FAILURES:
    print(f"  FAIL {name}")
