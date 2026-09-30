"""dc_extended.simulate — trajectory prediction per candidate. Code only."""
from .gates import day_dead, WEEKDAYS


def _has_affordance(text):
    """A question, a callback hook, or an open loop she can pick up."""
    if "?" in text:
        return True
    hooks = ("your turn", "top that", "beat that", "try me", "go on",
             "the mean", "the recipe")
    return any(h in text.lower() for h in hooks)


def _matches_live_bit(text, case):
    live = [r for r in case.get("referents", []) if r.get("status") == "live"]
    for r in live:
        tokens = [w for w in r["text"].lower().split() if len(w) > 3]
        if tokens and any(t in text.lower() for t in tokens):
            return True
    return False


def _names_dead_day(text, case, now):
    low = text.lower()
    for w in WEEKDAYS:
        if w in low and day_dead(w, now):
            return True
    return False


def _is_plan(text, claim):
    return (claim or {}).get("type") == "plan" or any(
        w in text.lower() for w in ("dinner", "drinks", "book", "pick the place"))


def _is_generic(text):
    generic = ("you're gorgeous", "you seem fun", "you're cute", "you're so hot",
               "hey", "hi", "what's up", "wyd")
    return any(g in text.lower() for g in generic)


def simulate(candidate, case, read, claim, now):
    """Predict her response class set. Deterministic, testable."""
    classes = set()

    # no affordance → silence likely
    if not _has_affordance(candidate):
        classes.add("silence")

    # plan on a dead day → topic_change
    if _is_plan(candidate, claim) and _names_dead_day(candidate, case, now):
        classes.add("topic_change")

    # generic opener → topic_change or silence
    if _is_generic(candidate):
        classes.add("topic_change")
        classes.add("silence")

    # matches live bit + heat >= 1 → extends
    if _matches_live_bit(candidate, case) and ((read or {}).get("heat") or 0) >= 1:
        classes.add("extends")

    # question + heat >= 1 → answers
    if "?" in candidate and ((read or {}).get("heat") or 0) >= 1:
        classes.add("answers")

    # refusal read → silence
    if (read or {}).get("kind") == "refuse":
        classes.add("silence")

    if not classes:
        classes.add("answers")

    return classes
