"""Tranche 1 acceptance: immutable truth stream + authority boundary (shadow)."""
import sys, os, tempfile, json
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from dc_extended.events import EventLog, project

def _try(log, meth, *a):
    try:
        getattr(log, meth)(*a)
        return True
    except Exception:
        return False

FAIL = []
def t(name, cond):
    print(("PASS " if cond else "FAIL ") + name)
    if not cond: FAIL.append(name)

d = tempfile.mkdtemp()
p = os.path.join(d, "case1.jsonl")
log = EventLog(p)

e1 = log.append("incoming_message", "her", "c1", {"text": "Hot"})
e2 = log.append("operator_directive", "operator", "c1", {"text": "stop", "directive": "stop"})
t("append: hash chain verifies", log.verify())

pid = log.propose("c1", "say that to my face")
t("NEG: confirm without a proposal is impossible",
  _try(log, "confirm_send", "c1", "nonexistent", "say that to my face") is False)
t("NEG: payload mismatch cannot confirm",
  _try(log, "confirm_send", "c1", pid, "a different line") is False)
sid = log.confirm_send("c1", pid, "say that to my face")
t("confirm with matching payload succeeds", isinstance(sid, str))
st = project(log)
t("projection: send recorded, proposal consumed", len(st["confirmed_sends"]) == 1
  and len(st["open_proposals"]) == 0)

rid = log.append("incoming_reply", "her", "c1", {"text": "leave me alone"})
log.interpret_outcome("c1", rid, "extends")
t("interpretation: chain still verifies", log.verify())
raw = [e for e in log.read() if e["type"] == "incoming_reply"][0]
t("raw reply text unchanged after interpretation", raw["payload"]["text"] == "leave me alone")
log.interpret_outcome("c1", rid, "refusal", source="operator")
t("correction appends; history intact", log.verify()
  and len([e for e in log.read() if e["type"] == "outcome_interpretation"]) == 2)

t("explicit operator stop closes the projection", st["closed"] is True)

p2 = os.path.join(d, "case2.jsonl"); log2 = EventLog(p2)
log2.append("incoming_message", "her", "c2", {"text": "stop being cute 😂"})
t("NEG CONTROL: her playful 'stop being cute' is NOT a directive",
  project(log2)["closed"] is False)
t("NEG: unknown actor rejected", _try(log2, "append", "incoming_message", "martian", "c2", {}) is False)
t("NEG: unknown event type rejected", _try(log2, "append", "bogus", "her", "c2", {}) is False)

print()
print("ALL TRANCHE-1 TESTS PASS" if not FAIL else "FAILURES: %s" % FAIL)
assert not FAIL, FAIL
