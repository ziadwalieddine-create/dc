"""dc_extended.critic — first-hit string kills. Code only; no model judgment."""
import re
from datetime import datetime

KILL_IDS = [
    "not_referable", "claim_swap", "clock", "not_particular",
    "label_question", "interview", "self_decode", "aimed_at_her",
    "missed_heat", "binds_her", "negotiates_no", "dead_end",
    "length", "echo", "unnatural", "logistics_mush",
    "frame_limit", "known_fact", "restart",
    "braided_prompts", "self_congratulatory", "commanding",
    "presumptive_access", "energy_monitoring", "ambiguity_reopen",
    "manufactured_pull", "homework_assignment",
]

from .lexicon import LABELS, MUSH, SELF_DECODE, BINDS, AIMED, SPECIMEN, HINGE_OPEN, GENERIC_LINES


def _norm(text):
    return " " + (text or "").lower() + " "


def critic(text, case, read=None, sim_classes=None, other_wider=False):
    """First kill or None. Does not return a replacement string."""
    if not text:
        return "not_referable"
    low = _norm(text)
    facts = [f.lower() for f in case.get("unique_facts", [])]

    # 1. claim_swap — booking under wrong claim type
    claim_type = (case.get("_claim") or {}).get("type")
    if claim_type != "plan" and any(w in low for w in ("thursday at", "let's do ", "book us")):
        return "claim_swap"
    if claim_type != "channel" and any(w in low for w in ("ig:", "instagram:", "dm me at")):
        return "claim_swap"

    # 2. clock — dead weekday token
    from .gates import day_dead, WEEKDAYS
    now = case.get("_now")
    if isinstance(now, str):
        now = datetime.fromisoformat(now)
    if now:
        for w in WEEKDAYS:
            if re.search(r"\b" + w + r"\b", low) and day_dead(w, now):
                return "clock"

    # 3. not_particular — GENERIC_LINES membership
    if any(g in text.lower() for g in GENERIC_LINES):
        return "not_particular"

    # 4. label_question
    if any(q in text.lower() for q in LABELS):
        return "label_question"

    # 5. interview
    if text.count("?") > 1:
        return "interview"

    # 6. self_decode
    if any(s in low for s in SELF_DECODE):
        return "self_decode"

    # 7. aimed_at_her
    if any(s in text.lower() for s in AIMED):
        return "aimed_at_her"

    # 8. missed_heat — flat while she extended (run_turn sets _policy/_amplitude)
    if (read and read.get("she_extended") and case.get("_amplitude", 0) == 0
            and case.get("_policy") in ("Play", "Trade")):
        return "missed_heat"

    # 9. binds_her
    if any(s in text.lower() for s in BINDS):
        return "binds_her"

    # 10. negotiates_no — claim family is banned (real state, no flag)
    fam = (case.get("_claim") or {}).get("family")
    if fam and any((b.get("family") == fam) for b in case.get("banned", [])):
        return "negotiates_no"

    # 11. dead_end
    if sim_classes == {"silence", "topic_change"} and other_wider:
        return "dead_end"

    # 12. length
    if case.get("her_len") in (None, "short") and len(text) > 220:
        return "length"

    # 13. echo
    if any(text.strip().lower() == (b.get("string") or "").lower() for b in case.get("banned", [])):
        return "echo"
    if SPECIMEN in text.lower():
        return "echo"

    # 14. unnatural
    if "not x but y" in text.lower() or text.lower().startswith("i hope this helps"):
        return "unnatural"

    # 15. logistics_mush
    if any(m in text.lower() for m in MUSH):
        return "logistics_mush"

    # 16. frame_limit
    if "no_arrangement_label" in case.get("frame_limits", []) and any(
            w in text.lower() for w in ("fwb", "situationship")):
        return "frame_limit"
    if "no_relationship_speech" in case.get("frame_limits", []) and "our relationship" in text.lower():
        return "frame_limit"

    # 17. known_fact
    if "?" in text:
        for fact in facts:
            if len(fact) > 8 and fact in text.lower():
                return "known_fact"

    # 18. restart
    if case.get("channel") in ("ig", "sms") and any(p in text.lower() for p in HINGE_OPEN):
        return "restart"

    # 19. braided_prompts — more than one profile object braided into one line
    profile_objects = sum(1 for obj in ("formula 1", "batman", "garlic bread", "pool",
                                         "camping", "shopping", "karaoke", "motorbike")
                          if obj in low)
    if profile_objects >= 2:
        return "braided_prompts"

    # 20. self_congratulatory — compliments himself instead of her
    if any(p in low for p in ("good taste in owners", "i'm not boring", "not my fault")):
        return "self_congratulatory"

    # 21. manufactured_pull — creates a same-day meet she didn't set up
    if any(p in low for p in ("come ruin my", "come over tonight", "come through tonight",
                               "right now", "tonight then")):
        return "manufactured_pull"

    # 22. commanding — tells her what to do
    # NOTE: "prove it" deliberately absent — verified state-dependent (katherine msg 43,
    # strong positive answering her explicit challenge). Governed by the licence ladder,
    # not a string match.
    if any(p in low for p in ("show me", "send me a pic", "send a pic",
                               "save me a spot", "come have")):
        return "commanding"

    # 23. presumptive_access — assumes she's already interested/invested
    if any(p in low for p in ("you're texting me", "you in my dms", "you're already",
                               "you said anytime", "you said i'd be welcome")):
        return "presumptive_access"

    # 24. energy_monitoring — watches her behavior without a meetup context
    if any(p in low for p in ("same energy over here", "don't go quiet", "don't disappear")):
        return "energy_monitoring"
    if "keep that same energy" in low and "sitting across" not in low and "in person" not in low:
        return "energy_monitoring"

    # 25. ambiguity_reopen — asks what she meant instead of answering
    if "what did you think i meant" in low or "what do you think i meant" in low:
        return "ambiguity_reopen"

    # 26. homework_assignment — asks for a photo/proof instead of flirting
    if any(p in low for p in ("send a pic", "send me a pic", "show me a pic",
                               "pic of", "photo of")):
        return "homework_assignment"

    return None
