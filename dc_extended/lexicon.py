"""dc_extended.lexicon — one wordlist, two consumers.

The read layer and the critic must never drift apart: a word added to one
list and not the other is a latent bug no single-side test can catch.
Identity (is) is the guarantee; test_phase3 asserts it.
"""

# --- read layer ---
DRY_TOKENS = {"haha", "ok", "okay", "lol", "k", "kk", "👍", "😂", "😊"}
REFUSAL_TOKENS = ("no", "stop", "block", "not interested", "leave me alone")
WITHDRAW = ("busy", "later", "not sure", "maybe", "we'll see", "another time")
CHARGE = ("come over", "tonight", "now", "prove it", "show me", "bet")

# --- claim map ---
AUDIT_WORDS = ("audit", "review", "analyze", "compare", "board", "report", "ceremony")

# --- critic ---
LABELS = (
    "what's your dream", "what is your dream", "biggest fear", "childhood",
    "what do you do", "where are you from", "how do you know you like",
    "most attractive", "perfect day", "crystal ball",
    "what do your weeknights", "what are you enjoying", "what does your",
    "what's your usual", "what are your weekends",
)
MUSH = ("when are you free", "your call", "pick a day", "you pick")
SELF_DECODE = ("that was a joke", "just kidding", " jk", "haha that", "lol that")
BINDS = ("you have to", "reply now", "why didn't you", "you ignored", "answer me")
AIMED = ("women always", "you girls", "so stupid", "your body", "you're dumb")
GENERIC_LINES = ("you're gorgeous", "you seem fun", "you're cute", "you're so hot")
HINGE_OPEN = ("nice to match", "found you on hinge", "we matched", "hey from hinge")
SPECIMEN = "you kept the basil and threw the recipe. that tracks."

# --- frame certification ---
FILLER_WORDS = ("just", "really", "literally", "tbh", "imo", "sort of", "kind of")

FORECAST_CLASSES = ("extends", "answers", "counters", "mirrors_length",
                    "literal", "dry", "topic_change", "silence", "refuses")
