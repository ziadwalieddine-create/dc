"""Queue tests."""
import sys, os, shutil, tempfile
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from datetime import datetime
from dc_extended.queue import build_queue

NOW = datetime(2026,9,19,12,0)

def t(name, cond):
    print(("PASS " if cond else "FAIL ") + name)
    if not cond: sys.exit(1)

mary = {"id":"mary","name":"Mary","stage":"off_app","channel":"ig",
        "ig_handoff_ts":"2026-09-13T20:00:00","updated":"2026-09-13T20:05:00",
        "last_send":None,"goal":{"type":"date"}}
rae = {"id":"rae","name":"Rae","stage":"off_app","channel":"ig",
       "ig_handoff_ts":"2026-09-17T21:40:00","updated":"2026-09-19T10:00:00",
       "last_send":{"text":"Found you. Still owe me the mean.","confirmed":True,"family":"callback"},
       "open_plan":{"what":"dinner","when_day":"saturday","proposed_ts":"2026-09-18T20:00:00","confirmed":False},
       "goal":{"type":"fwb"}}
camille = {"id":"camille","name":"Camille","closed":True,"updated":"2026-09-20T09:00:00"}

q = build_queue([mary, rae, camille], NOW)
acts = {(i["id"], i["action"]) for i in q}

t("mary off_app stale at 48h -> reverify", ("mary","reverify") in acts)
# same staleness, different stage: a confirmed-send ig thread also gets reverify
rae2 = dict(rae); rae2["updated"] = "2026-09-13T20:05:00"
q2 = build_queue([rae2], NOW)
t("rae off_app stale -> reverify too", ("rae","reverify") in
  {(i["id"], i["action"]) for i in q2})
t("camille excluded", not any(i["id"]=="camille" for i in q))
t("rae plan confirm item", ("rae","confirm_or_repair_plan") in acts)
t("ordered by deadline", q[0]["deadline"] <= q[-1]["deadline"])
print("ALL QUEUE TESTS PASS")
