"""Gate unit tests."""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from datetime import datetime
from dc_extended.gates import run_gates, day_dead

def base_case(**kw):
    c = {"id":"t","name":"Test","stage":"off_app","channel":"ig",
         "goal":{"type":"fwb"},"updated":"2026-09-19T10:00:00",
         "banned_families":[],"boundaries":[],"unique_facts":[],
         "missed_windows":[],"referents":[]}
    c.update(kw); return c

NOW = datetime(2026,9,18,18,11)  # a Friday evening
assert NOW.weekday() == 4

def t(name, cond):
    print(("PASS " if cond else "FAIL ") + name)
    if not cond: sys.exit(1)

# dead_clock: friday token on friday evening
h = run_gates(base_case(), "wyd on friday", NOW, claim={"type":"plan","family":"wyd_plan"})
t("dead_clock friday-evening", any(g=="dead_clock" for g,_ in h))

# dead_clock silent earlier same day
h = run_gates(base_case(), "wyd on friday", datetime(2026,9,18,10,0), claim={"type":"plan","family":"plan"})
t("dead_clock quiet at 10am", not any(g=="dead_clock" for g,_ in h))

# rae rule: sunday synonym for missed saturday
h = run_gates(base_case(missed_windows=[{"day":"saturday"}]), "sunday then", NOW, claim={"type":"plan","family":"plan"})
t("dead_clock synonym", any(g=="dead_clock" for g,_ in h))

# restack
h = run_gates(base_case(banned_families=["wyd_plan"]), "wyd?", NOW, claim={"type":"plan","family":"wyd_plan"})
t("restack", any(g=="restack" for g,_ in h))

# dual_plan
h = run_gates(base_case(open_plan={"what":"dinner","when_day":"saturday","proposed_ts":"2026-09-17T20:00:00"}),
              "drinks instead", NOW, claim={"type":"plan","family":"plan"})
t("dual_plan", any(g=="dual_plan" for g,_ in h))

# goal_smuggle
h = run_gates(base_case(goal={"type":"get_to_know"}), "dinner wednesday", NOW, claim={"type":"plan","family":"plan"})
t("goal_smuggle", any(g=="goal_smuggle" for g,_ in h))

# unconfirmed_send
h = run_gates(base_case(last_send={"text":"the dm","confirmed":False}),
              "follow up", NOW, claim={"type":"answer","family":"follow_up","builds_on_last":True})
t("unconfirmed_send", any(g=="unconfirmed_send" for g,_ in h))

# wrong_platform
h = run_gates(base_case(channel="ig"), "new opener", NOW, claim={"type":"channel","family":"new_opener","channel":"hinge"})
t("wrong_platform", any(g=="wrong_platform" for g,_ in h))

# stated_boundary
h = run_gates(base_case(boundaries=["photos"]), "sent you photos", NOW, claim={"type":"bit","family":"tease"})
t("stated_boundary", any(g=="stated_boundary" for g,_ in h))

# cross_thread + wrong_person
other = {"id":"j","name":"Jisu","unique_facts":["skydiving"]}
h = run_gates(base_case(), "jisu told me skydiving", NOW, claim={"type":"bit","family":"tease"}, all_cases=[base_case(), other])
t("cross_thread", any(g=="cross_thread" for g,_ in h))
t("wrong_person", any(g=="wrong_person" for g,_ in h))

# stale_goal
h = run_gates(base_case(closed=True, close_cause="archived"), "hey", NOW, claim={"type":"bit","family":"tease"})
t("stale_goal", any(g=="stale_goal" for g,_ in h))

print("ALL GATE TESTS PASS")
