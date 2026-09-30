"""Read layer unit tests — every kind from engine/read.md, heat, reversal."""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from dc_extended.switchboard import read_bubble


def base_case(**kw):
    c = {"id": "t", "name": "Test", "stage": "chat", "channel": "hinge",
         "goal": {"type": "date"}, "updated": "2026-09-19T10:00:00",
         "banned_families": [], "boundaries": [], "unique_facts": [],
         "missed_windows": [], "referents": []}
    c.update(kw)
    return c


def t(name, cond):
    print(("PASS " if cond else "FAIL ") + name)
    if not cond:
        sys.exit(1)


BASIL = {"id": "b", "text": "basil plant", "type": "bit",
         "introduced_by": "him", "touched_by_her": True, "status": "live"}

# none — no bubble
r = read_bubble(None, base_case())
t("none", r["kind"] == "none" and r["heat"] is None)

# refuse — hard no
r = read_bubble("no, I'm not interested", base_case())
t("refuse", r["kind"] == "refuse" and r["heat"] == 0)

# refuse beats tease in same bubble
r = read_bubble("no stop it haha", base_case(referents=[BASIL]))
t("refuse beats tease", r["kind"] == "refuse")

# dry — token
r = read_bubble("haha", base_case(referents=[BASIL]))
t("dry", r["kind"] == "dry" and r["heat"] == 0 and not r["she_extended"])

# dry — bare hey
r = read_bubble("hey", base_case())
t("dry hey", r["kind"] == "dry")

# escalate — charge
r = read_bubble("come over tonight", base_case())
t("escalate", r["kind"] == "escalate" and r["heat"] == 2 and r["she_extended"])

# tease — touched live bit
r = read_bubble("the basil plant is unhinged", base_case(referents=[BASIL]))
t("tease", r["kind"] == "tease" and r["heat"] == 2 and r["she_extended"])

# tease — no live bit, no tease
r = read_bubble("the basil plant is unhinged", base_case())
t("no tease without live bit", r["kind"] != "tease")

# offer — she named the move
r = read_bubble("we should get drinks sometime", base_case())
t("offer", r["kind"] == "offer" and r["heat"] == 1)

# logistics — asked about a day
r = read_bubble("what day works for you?", base_case())
t("logistics", r["kind"] == "logistics" and r["heat"] == 1)

# question
r = read_bubble("what do you mean by that?", base_case(referents=[BASIL]))
t("question", r["kind"] == "question")

# withdraw — busy, no counter
r = read_bubble("I'm busy this week", base_case())
t("withdraw", r["kind"] == "withdraw" and r["heat"] == 0)

# withdraw — not sure
r = read_bubble("not sure yet", base_case())
t("withdraw not sure", r["kind"] == "withdraw")

# withdraw — busy with counter is NOT withdraw
r = read_bubble("busy today but how about tomorrow?", base_case())
t("counter not withdraw", r["kind"] != "withdraw")

# topic — live bit exists, bubble is not it
r = read_bubble("I just got back from the coast", base_case(referents=[BASIL]))
t("topic", r["kind"] == "topic" and r["heat"] == 0)

# invest — new specific
r = read_bubble("I just got back from the coast, it was freezing", base_case())
t("invest", r["kind"] == "invest" and r["heat"] == 1)

# every read has a reversal
r = read_bubble("haha", base_case())
t("reversal present", bool(r.get("reversal")))

# every read has confidence
t("confidence present", r.get("confidence") in ("high", "split"))

print("ALL READ TESTS PASS")
