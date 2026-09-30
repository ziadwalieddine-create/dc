"""dc_extended.read_v2 — multidimensional Read, ADDITIVE.

Produces the rich representation BESIDE the legacy kind, which is derived as a
compatibility projection. Orthogonal dimensions coexist: logistics_question +
tease + heat 2 is legal. Competing hypotheses are only for genuinely incompatible
interpretations, not orthogonal features.

Heat preservation is a non-regression constraint: representing more literal
function must never reduce evidence-supported heat. Heat is bound to the evidence
that supports it — never inherited from unrelated old referents.
"""
import re
from .readlayer import read_bubble as legacy_read

CHARGE = ("come over", "tonight", "prove it", "show me")
CHARGE_WORDS = ("now", "bet")
TEASE_WORDS = ("cute", "hot", "stop", "naughty", "trouble")
DAY_WORDS = ("monday","tuesday","wednesday","thursday","friday","saturday","sunday")
ATTRACTION_CUES = {"hot", "cute", "handsome", "sexy", "gorgeous"}


def _words(b):
    return set(re.findall(r"[a-z']+", b))


def read_v2(bubble, case=None):
    """Multidimensional read. Never reduces evidence-supported heat."""
    case = case or {}
    legacy = legacy_read(bubble, case)
    b = (bubble or "").strip().lower()
    w = _words(b)
    r = {
        "incoming_ref": bubble,
        "literal_functions": [],
        "relational_functions": [],
        "question_objects": [],
        "response_obligations": [],
        "heat": None,
        "investment": None,
        "refusal_scope": None,
        "acceptance": None,
        "progression_signal": None,
        "flirt_license": "none",
        "confidence": "high",
        "kind": legacy["kind"],          # compatibility projection
        "reversal_evidence": None,
    }
    if not b:
        r["literal_functions"] = ["none"]
        return r

    has_q = "?" in b or b.startswith(("where","what","when","how","why","is ","are ","do ","does "))
    playful = any(e in b for e in ("😂","🤣","😏","😉","haha","ahaha","lol","jk"))

    # --- literal functions (independent of relational) ---
    if has_q:
        r["literal_functions"].append("question")
        if any(d in b for d in DAY_WORDS) or "when" in b or "what day" in b:
            r["literal_functions"].append("logistics_question")
        if any(tok in b for tok in ("where", "location", "what area", "which bar", "what venue")):
            r["literal_functions"].append("logistics_question")
            r["question_objects"].append("location")
            r["response_obligations"].append("answer_location")
        elif has_q:
            r["response_obligations"].append("answer_question")
    if any(d in b for d in DAY_WORDS) and not has_q:
        r["literal_functions"].append("logistics_statement")
    if re.search(r"(work|gym|dinner|walk|family|busy|tired)", b) and len(w) >= 3:
        r["literal_functions"].append("disclosure")

    # --- relational functions (independent of literal) ---
    live = [x for x in case.get("referents", []) if x.get("status") == "live"]
    AMBIGUOUS = {"mean", "just", "like", "right", "back", "down"}  # audit finding 49:
    # loose substring overlap must not create false heat ("I mean honestly" != the bit)
    def touches(t):
        tt = t["text"].lower()
        if tt in b and tt not in AMBIGUOUS:
            return True
        return any(len(x) > 3 and x not in AMBIGUOUS and re.search(r"\b%s\b" % re.escape(x), b)
                   for x in tt.split())
    touched = any(touches(t) for t in live)
    if touched and (playful or len(w) >= 2):
        r["relational_functions"].append("tease")
    if any(wd in w for wd in ("naughty","trouble","dangerous","dare","prove","show","make me")):
        r["relational_functions"].append("challenge")

    # --- compliments: positive evidence about HIM (misread licence binds here) ---
    if re.search(r"\b(only good one|the best|love (that|it|this)|you're (funny|cute|great)|good one|love your)\b", b):
        r["relational_functions"].append("compliment")
        if r["heat"] is None: r["heat"] = 1
        r["investment"] = "positive"

    # --- heat: bound to evidence, never to token count alone ---
    if b in ATTRACTION_CUES or (b in ("hot","cute") ):
        r["heat"] = 2
        r["relational_functions"].append("attraction_cue")
        r["investment"] = "positive"
    elif any(c in b for c in CHARGE) or (any(re.search(r"\b%s\b" % x, b) for x in CHARGE_WORDS) and playful):
        neg = re.search(r"\b(don't|do not|won't)\s+come", b)
        if neg:
            r["refusal_scope"] = "object_specific"
            r["heat"] = 0
        else:
            r["heat"] = 2
    elif playful and touched:
        r["heat"] = 2
    elif playful and has_q and ("where" in b or "what" in b or "when" in b) and touched:
        r["heat"] = 2   # playful logistics WITH a live referent: the Wynyard case
    elif touched:
        r["heat"] = 1
    elif len(w) >= 4 and not has_q:
        ENGAGE = ("curious", "love", "excited", "amazing", "obsessed", "adore", "fun")
        LOGISTICAL = ("train", "station", "leaves", "arrives", "address", "pm", "am", "uber", "bus")
        if any(x in w for x in ENGAGE):
            r["heat"] = 1
            r["investment"] = "positive"
        elif not any(x in w for x in LOGISTICAL):
            r["heat"] = 1   # personal disclosure, not logistics

    # short answers are NOT low investment by themselves ("I'm in", "Deal", "Hot")
    if len(w) <= 3 and r["heat"] is None and not has_q:
        if any(x in w for x in ("in","deal","yes","yeah","down","hot","cute","definitely")):
            r["heat"] = 1
            r["investment"] = "positive"

    # --- refusal scope: object-specific vs global (the audit's paired invariant) ---
    if "not sure" in b or "not free" in b:
        m3 = re.search(r"not free\s+(\w+)", b)
        if m3:
            r["refusal_scope"] = m3.group(1)
            if "but" in b:
                r["progression_signal"] = "positive"
                r["heat"] = r["heat"] or 1
    elif re.search(r"\bno\b", b):
        m = re.search(r"no\s+(\w+).*?(but|however|,)?\s*(\w+)\s+works", b)
        if m:
            r["refusal_scope"] = m.group(1)
            r["acceptance"] = m.group(3)
            r["progression_signal"] = "positive"
            r["heat"] = r["heat"] or 1
        elif playful:
            r["refusal_scope"] = "playful_resistance"   # "no way 😂", "nooo"
        elif any(c in b for c in CHARGE) or "come over" in b:
            r["refusal_scope"] = "object_specific"
        else:
            r["refusal_scope"] = "global"
    elif any(s in b for s in ("stop being","stop it","stoppp","stop being cute")) and playful:
        r["relational_functions"].append("playful_resistance")
        r["refusal_scope"] = "playful_resistance"
        if r["heat"] is None: r["heat"] = 1

    # --- flirt licence from evidence ---
    if r["heat"] == 2 and (touched or "attraction_cue" in r["relational_functions"]):
        r["flirt_license"] = "high"
    elif r["heat"] and playful:
        r["flirt_license"] = "moderate"

    # heat must never be reduced by coexisting literal function (non-regression):
    # literal_functions are recorded above WITHOUT touching r["heat"].
    return r


def legacy_kind(r):
    """Compatibility projection for downstream code during migration."""
    if "attraction_cue" in r["relational_functions"]:
        return "escalate"
    if "tease" in r["relational_functions"] and r["heat"] == 2:
        return "tease"
    if "question" in r["literal_functions"]:
        return "question"
    if r["refusal_scope"] == "global":
        return "refuse"
    return "invest" if r["heat"] else "dry"
