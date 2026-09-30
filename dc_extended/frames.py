"""dc_extended.frames — frame certification. Certification IS the test suite.

A frame is a word-level blueprint: slots declare what each span does, forecast
carries build-time-certified response classes, protected lists the tokens that
must survive surface edits. The mechanical failure shapes are oracles here;
full word-level judgment stays with the model at authoring time.
"""
import re
from .writer import _fill
from .critic import critic
from .lexicon import FILLER_WORDS, FORECAST_CLASSES

FRAME_FIELDS = ("frame", "id", "family", "text", "outcome", "slots", "forecast", "protected")


def load_frames(path):
    import json
    rows = []
    with open(path) as f:
        for line in f:
            line = line.strip()
            if line:
                rows.append(json.loads(line))
    return rows


def fixture_case(referent="basil plant", name="Rae"):
    return {"id": "cert", "name": name, "stage": "chat", "channel": "ig",
            "goal": {"type": "date"}, "banned_families": [], "boundaries": [],
            "unique_facts": [referent], "missed_windows": [],
            "referents": [{"id": "b", "text": referent, "type": "bit",
                           "introduced_by": "him", "touched_by_her": True, "status": "live"}],
            "frame_limits": [], "banned": [], "her_len": "short",
            "_now": None, "_claim": {"type": "bit", "family": "tease"},
            "ig_handle": "r.tmzon", "amplitude": 1}


def filler_hits(text):
    low = " " + text.lower() + " "
    return [w for w in FILLER_WORDS if re.search(r"\b" + re.escape(w) + r"\b", low)]


def certify_frame(frame, case_a=None, case_b=None):
    """Return a list of defects. Empty list = certified.

    Oracle checks: required fields, variant critic-pass under fixtures,
    particularity (referent-swapped fills must differ when slots carry one),
    filler-lexicon lint, forecast validity.
    """
    defects = []
    for field in FRAME_FIELDS:
        if field not in frame:
            defects.append("missing_field:%s" % field)
    if defects:
        return defects
    if not isinstance(frame["text"], (str, list)) or not frame["text"]:
        return ["no_variants"]

    variants = frame["text"] if isinstance(frame["text"], list) else [frame["text"]]
    case_a = case_a or fixture_case()
    case_b = case_b or fixture_case(referent="rooftop bar", name="Jisu")

    has_referent_slot = any("{referent}" in v for v in variants)
    filled_a, filled_b = [], []
    for v in variants:
        fa = _fill(v, case_a)
        fb = _fill(v, case_b)
        filled_a.append(fa)
        filled_b.append(fb)
        kill = critic(fa, dict(case_a), read={"kind": "tease", "heat": 2, "she_extended": True})
        if kill:
            defects.append("critic:%s:%s" % (kill, fa[:40]))
        hits = filler_hits(fa)
        if hits:
            defects.append("filler:%s:%s" % (",".join(hits), fa[:40]))

    if has_referent_slot and filled_a == filled_b:
        defects.append("not_particular:swapped fills identical")

    for cls in frame["forecast"]:
        if cls not in FORECAST_CLASSES:
            defects.append("bad_forecast_class:%s" % cls)
    if not frame["forecast"]:
        defects.append("empty_forecast")
    return defects


def certify_corpus(path):
    """Run the oracle over every authored frame. Returns {frame_id: defects}."""
    frames = [r for r in load_frames(path) if r.get("frame") is True]
    seen = set()
    report = {}
    for fr in frames:
        if fr["id"] in seen:
            report[fr["id"]] = ["duplicate_id"]
        seen.add(fr["id"])
        report.setdefault(fr["id"], certify_frame(fr))
    return report


SLOT_REASONS = {
    "referent":  "anchors to the live specific she can point at",
    "challenge": "answers the charge she extended instead of changing topic",
    "tease":     "keeps the bit warm without explaining it",
    "callback":  "reuses an object she already touched",
    "channel-fact": "grounds the move in a FACT she gave, not an invention",
    "question":  "gives her a probe she can actually answer",
    "affordance": "leaves her something to add",
    "plan":      "names one executable step, he picks",
    "decisive":  "he decides; she may counter, not choose",
}


def justification_card(frame):
    """Every declared slot gets its one-line reason. Audit turns are lookups."""
    return {"frame": frame.get("id"),
            "lines": ["%s — %s" % (s, SLOT_REASONS.get(s, "declared slot, reason at authoring time"))
                      for s in frame.get("slots", [])],
            "forecast": frame.get("forecast"),
            "protected": frame.get("protected", [])}


def card_for(frame_id, corpus_path):
    for row in load_frames(corpus_path):
        if row.get("id") == frame_id:
            return justification_card(row)
    return None
