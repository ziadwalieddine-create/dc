"""dc_extended.ingress — structured ingress. HER evidence and OPERATOR directives
occupy separate namespaces BEFORE any interpretation.

Design rules (replacing substring heuristics):
- Command words match WHOLE TOKENS, never substrings. "preview" does not contain the
  token "review"; "overboard" does not contain the token "board"; "tonight"/"big"/
  "handle" are never channel markers.
- Negation scopes to the phrase it precedes: "don't set the date" negates date-lock;
  "not off app" negates off-app. Negated intents are NOT intents.
- Channels are distinct: sms/number never maps to ig.
- Explicit operator STOP is locked first and dominates everything downstream.
- HER text is never parsed as control. "she: stop being cute 😂" is evidence.
"""
import re

STOP_TOKENS = {"stop", "close", "block"}
LEAVE_PHRASES = ("leave it", "do not contact", "don't contact")
NEGATORS = ("don't", "do not", "not", "never")
AUDIT_TOKENS = {"audit", "review", "analyse", "analyze", "board", "report", "ceremony"}
CHANNEL_TOKENS = {"ig": {"ig", "instagram", "insta"},
                  "sms": {"sms", "text", "number", "phone"}}
ENDPOINT_PHRASES = [("off app", "off_app"), ("off-app", "off_app"),
                    ("date lock", "date_lock"), ("lock the date", "date_lock"),
                    ("set the date", "date_lock"), ("book", "date_lock"),
                    ("long term", "date_lock"), ("banter", "banter"),
                    ("get to know", "banter"), ("play the bit", "banter")]
PLAN_TOKENS = {"wednesday", "thursday", "friday", "saturday", "sunday",
               "monday", "tuesday", "drinks", "dinner", "book", "lock"}


class Directive:
    """Typed operator intent. None means unspecified — never guessed."""
    def __init__(self):
        self.stop = False
        self.claim_type = None      # bit | plan | answer | close | channel
        self.endpoint = None        # off_app | date_lock | banter | close
        self.channel = None         # ig | sms (target; distinct)
        self.audit = False
        self.negated = []
        self.raw = None

    def __repr__(self):
        return ("Directive(stop=%s claim=%s endpoint=%s channel=%s audit=%s negated=%s)"
                % (self.stop, self.claim_type, self.endpoint, self.channel,
                   self.audit, self.negated))


def _tokens(text):
    return re.findall(r"[a-z0-9']+", (text or "").lower())


def parse_operator(text):
    """Parse operator text into a typed Directive. Whole-token matching; negation
    scopes to the following phrase; no substring collisions are possible."""
    d = Directive()
    d.raw = text
    if not text:
        return d
    low = text.lower()
    toks = _tokens(text)

    # --- explicit STOP: whole-token or fixed phrases, never negated ---
    neg_tok = any(n in toks for n in ("don't", "do", "not"))
    if any(p in low for p in LEAVE_PHRASES) or (
            any(tok in STOP_TOKENS for tok in toks) and not neg_tok):
        d.stop = True
        d.claim_type = "close"
        d.endpoint = "close"
        return d   # STOP locks; nothing else is parsed

    # --- negation scope: collect negated phrases ("don't set the date") ---
    for i, tok in enumerate(toks):
        if tok in ("don't", "do", "not", "never"):
            phrase = " ".join(toks[i + 1:i + 4])
            d.negated.append(phrase)

    def negated(*words):
        return any(all(w in n for w in words) for n in d.negated)

    # --- audit: whole tokens only ---
    if any(tok in AUDIT_TOKENS for tok in toks) and not negated("audit"):
        d.audit = True

    # --- channel: distinct targets, whole tokens ---
    for chan, words in CHANNEL_TOKENS.items():
        if any(tok in words for tok in toks) and not negated(chan):
            d.channel = chan
            d.claim_type = "channel"

    # --- endpoint phrases ---
    for phrase, ep in ENDPOINT_PHRASES:
        if phrase in low and not negated(*phrase.split()):
            d.endpoint = ep
            if ep == "date_lock":
                d.claim_type = "plan"

    # --- plan tokens ---
    if d.claim_type is None and any(tok in PLAN_TOKENS for tok in toks):
        d.claim_type = "plan"

    # --- bit ---
    if d.claim_type is None and any(tok in ("tease", "bit", "banter", "mean") for tok in toks):
        d.claim_type = "bit"

    return d


HER_MARKS = ("she:", "she said:", "her:", "her last:")


def split_trigger(text):
    """Structurally separate HER evidence from OPERATOR text before interpretation.
    Returns (her_text, operator_text). Her text is never parsed as control."""
    if not text:
        return None, None
    lines = text.splitlines()
    her, op = [], []
    for line in lines:
        l = line.strip().lower()
        if any(l.startswith(m) for m in HER_MARKS):
            her.append(line.split(":", 1)[-1].strip())
        else:
            op.append(line.strip())
    return (her[-1] if her else None), (" ".join(o for o in op if o) or None)


def parse_trigger(text):
    """Full ingress: split actors, then parse ONLY the operator side."""
    her, op = split_trigger(text)
    return her, parse_operator(op)
