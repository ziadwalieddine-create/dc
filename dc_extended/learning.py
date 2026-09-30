"""dc_extended.learning — learning on the trustworthy identity chain.

HistoricalMove != ReusableFrame. A line that HAPPENED stays historical even if
poor, rejected retrospectively, or non-portable — deleting "bad" sends would lose
negative evidence. Promotion to ReusableFrame is a SEPARATE certification with
its own gate (portability), never automatic from success.

Outcomes are DIMENSIONAL: reply, investment, attraction, heat, channel
progression, date progression, execution success, operator approval. Never one
collapsed "success" score. One reply event creates one outcome regardless of
bubble count. Caller-supplied labels are inputs, never evidence.
"""
import re

OUTCOME_DIMENSIONS = ("reply", "investment", "attraction", "heat",
                      "channel_progression", "date_progression",
                      "execution_success", "operator_approval")

POS_REPLY = ("extends", "counters", "answers", "mirrors_length", "accepted")
NEG_REPLY = ("dry", "silence", "topic_change")


def derive_outcome(reply_text, reply_class=None):
    """Dimensional outcome derived from RAW evidence. reply_class (caller-supplied)
    is recorded as an input, never silently replaces the raw text."""
    low = (reply_text or "").lower()
    d = {k: None for k in OUTCOME_DIMENSIONS}
    # raw evidence dominates: a refusal in the text overrides any caller label
    d["reply"] = "unclassified"
    d["label_supplied"] = reply_class   # recorded as input, never as classification
    if any(x in low for x in ("haha","lol","ok","k","👍","oki")) and len(low.split()) <= 3:
        d["investment"] = "low"
    elif low.strip():
        d["investment"] = "positive" if len(low.split()) >= 4 else "neutral"
    if re.search(r"\b(hot|cute|gorgeous|sexy|love|adore|miss you)\b", low):
        d["attraction"] = "positive"
    if re.search(r"\b(hot|horny|come over|want you|tonight)\b", low):
        d["heat"] = "positive"
    if re.search(r"\b(ig|instagram|whatsapp|number|text me)\b", low):
        d["channel_progression"] = "positive"
    if re.search(r"\b(wednesday|thursday|friday|saturday|sunday|yes|i'm in|works|deal)\b", low):
        d["date_progression"] = "positive"
    if re.search(r"\b(leave me alone|stop|block|don't contact)\b", low):
        d["attraction"] = "negative"
        d["investment"] = "negative"
    # a refusal can NEVER record as successful extension (finding: A4)
    if d["attraction"] == "negative" and d["reply"] in POS_REPLY:
        d["reply"] = "contradicts_raw_evidence"
    return d


class HistoricalMove:
    """Immutable record of what happened. Independent properties — a line can be
    actually_sent=True, operator_later_rejected=True, received_positive_reply=True
    simultaneously. No forced positive/negative label."""
    def __init__(self, move_id, text, case_id, **props):
        self.move_id = move_id
        self.text = text
        self.case_id = case_id
        for k in ("actually_sent", "operator_authored", "operator_approved",
                  "observed_response", "style_grade", "verification_status",
                  "operator_later_rejected"):
            setattr(self, k, props.get(k))

    def as_dict(self):
        return dict(self.__dict__)


IDENTITY_PATTERNS = ("wynyard", "hirst", "katherine", "gabby", "gabrielle",
                     "cullen", "basil", "leo", "parramatta", "prospect")
GENERIC_OK = ("drinks", "dinner", "gym", "walk", "sunday", "weather", "city")


def certify_portability(hist, other_case_names=()):
    """Promotion gate. Returns (ok, reasons). A successful historical sentence is
    NOT automatically a portable frame (finding 117)."""
    reasons = []
    low = hist.text.lower()
    for pat in IDENTITY_PATTERNS:
        if pat in low:
            reasons.append("case-specific token: %s" % pat)
    for name in other_case_names:
        if name and name.lower() in low:
            reasons.append("identity: %s" % name)
    if hist.verification_status not in ("verified", "operator_corrected"):
        reasons.append("not verified: %s" % hist.verification_status)
    if not hist.actually_sent:
        reasons.append("never confirmed sent")
    return (not reasons, reasons)


def promote_to_frame(hist, other_case_names=()):
    """Separate certification step. Historical success alone never promotes."""
    ok, reasons = certify_portability(hist, other_case_names)
    if not ok:
        return None, reasons
    return {"frame_id": "frame_" + hist.move_id, "portable_structure": hist.text,
            "source_case": hist.case_id, "certified": True}, []
