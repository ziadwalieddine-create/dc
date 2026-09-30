"""dc_extended.moves — moves as plain dicts: text / variants / sequence.

No Move class (agreed settlement). A move is:
  {"text": str, "variants": [str], "sequence": [str]}
- text: the single-string view
- variants: alternative phrasings (pick one)
- sequence: ordered bubbles (a multi-bubble move)
A single-bubble move is {"text": s, "variants": [...], "sequence": [s]}.
Bubbles and variants are distinct: a sequence is never a variants list.
"""


def single(text, variants=None, family=None, heat=None, obligations=None):
    return {"text": text, "variants": variants or [text], "sequence": [text],
            "family": family, "heat": heat, "obligations": obligations or []}


def sequence(bubbles, family=None, heat=None, obligations=None):
    if not bubbles or not all(isinstance(b, str) and b.strip() for b in bubbles):
        raise ValueError("a sequence needs at least one non-empty bubble")
    return {"text": bubbles[0], "variants": [bubbles[0]], "sequence": list(bubbles),
            "family": family, "heat": heat, "obligations": obligations or []}


def from_frame_variants(frame, case=None, rng=None):
    """One move from a frame: variants stay variants; pick one, never a sequence."""
    variants = frame["text"] if isinstance(frame["text"], list) else [frame["text"]]
    idx = rng.randrange(len(variants)) if (rng and len(variants) > 1) else 0
    return single(variants[idx], variants=variants, family=frame.get("family"),
                  heat=frame.get("heat"))


def wynward_move(location_fact, playful_continuation="You'll find me"):
    """The invariant: answer the location + preserve tease/heat. One bubble or two —
    the move is correct when the obligation is satisfied and heat is carried,
    NOT because it is always two bubbles."""
    return sequence([location_fact, playful_continuation],
                    family="callback", heat=2, obligations=["answer_location"])
