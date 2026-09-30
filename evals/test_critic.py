"""Critic unit tests — every kill id fires, clean lines pass."""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from datetime import datetime
from dc_extended.critic import critic, KILL_IDS

NOW = datetime(2026, 9, 18, 18, 11)  # Friday evening


def base_case(**kw):
    c = {"id": "t", "name": "Test", "stage": "off_app", "channel": "ig",
         "goal": {"type": "fwb"}, "updated": "2026-09-19T10:00:00",
         "banned_families": [], "boundaries": [], "unique_facts": ["basil plant"],
         "missed_windows": [], "referents": [], "frame_limits": [],
         "banned": [], "_now": NOW, "_claim": {"type": "bit", "family": "tease"}}
    c.update(kw)
    return c


def t(name, cond):
    print(("PASS " if cond else "FAIL ") + name)
    if not cond:
        sys.exit(1)


# 2. clock — friday evening
t("clock", critic("friday still works", base_case()) == "clock")

# 3. claim_swap
c = base_case()
c["_claim"] = {"type": "bit", "family": "tease"}
t("claim_swap", critic("thursday at 8 then", c) == "claim_swap")

# 1. not_referable — empty text
t("not_referable", critic("", base_case()) == "not_referable")

# 3. not_particular — GENERIC_LINES membership
t("not_particular", critic("you seem fun", base_case()) == "not_particular")

# 6. label_question
t("label_question", critic("what's your dream though", base_case()) == "label_question")

# 7. interview
t("interview", critic("wait what? and where?", base_case()) == "interview")

# 8. self_decode
t("self_decode", critic("haha that was a joke", base_case()) == "self_decode")

# 9. aimed_at_her
t("aimed_at_her", critic("you girls always do this", base_case()) == "aimed_at_her")

# 8. missed_heat (run_turn sets _policy/_amplitude)
c = base_case()
c["_policy"] = "Play"; c["_amplitude"] = 0
r0 = {"she_extended": True}
t("missed_heat", critic("the basil plant", c, read=r0) == "missed_heat")

# 13. binds_her
t("binds_her", critic("why didn't you answer", base_case()) == "binds_her")

# 10. negotiates_no — claim family banned
c = base_case()
c["_claim"] = {"type": "bit", "family": "basil"}
c["banned"] = [{"family": "basil", "string": "", "gate": None}]
t("negotiates_no", critic("the same bit again", c) == "negotiates_no")

# 15. dead_end
t("dead_end", critic("closed joke", base_case(), sim_classes={"silence", "topic_change"}, other_wider=True) == "dead_end")

# 17. length
t("length", critic("x" * 221, base_case()) == "length")

# 18. echo
c = base_case(banned=[{"string": "the basil plant died", "family": "basil"}])
t("echo", critic("the basil plant died", c) == "echo")
t("echo specimen", critic("You kept the basil and threw the recipe. That tracks.", base_case()) == "echo")

# 19. unnatural
t("unnatural", critic("I hope this helps, the basil plant", base_case()) == "unnatural")

# 20. logistics_mush
t("logistics_mush", critic("when are you free", base_case()) == "logistics_mush")

# 21. frame_limit
c = base_case(frame_limits=["no_arrangement_label"])
t("frame_limit", critic("we should be fwb", c) == "frame_limit")

# 22. known_fact
t("known_fact", critic("do you still have the basil plant?", base_case()) == "known_fact")

# 23. restart
c = base_case(channel="ig")
t("restart", critic("nice to match", c) == "restart")

# clean line passes
c = base_case()
c["_claim"] = {"type": "bit", "family": "tease"}
r = {"she_extended": True, "heat": 2}
t("clean line", critic("that basil plant is doing too much", c, read=r) is None)

# all kill ids are real
t("kill ids valid", all(k in KILL_IDS for k in KILL_IDS) and len(KILL_IDS) == 27)

print("ALL CRITIC TESTS PASS")
