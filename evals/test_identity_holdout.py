"""D5 hold-out: the derived critic must not reject the operator's own sent lines
from cases OUTSIDE the annotation set. If a hold-out line hard-fails, a kill is
overfit to Gaby-shaped data, not learning him."""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from dc_extended.critic import critic

CASE = {"unique_facts": [], "her_len": "short", "channel": "hinge",
        "banned": [], "frame_limits": [], "closed": False}

# Operator-sent lines, audit-verified, cases NOT in the annotation corpus.
HOLDOUT = [
    ("anthea", "Good, because I'm very chalant about taking you for a drink"),
    ("anthea", "Does the cat in your photo need to approve first?"),
    ("hanan", "How bad do these side quests get? I need to know if I\u2019m bringing snacks or bail money"),
    ("lisa", "I ride a Yamaha, so I guess I'll have to make up for the lack of German engineering somehow"),
    ("lana", "it is, now i need to know if you\u2019re into the bike or just the guy riding it"),
    ("sam", "You definitely count getting dessert at 11pm as an adventure, don't you?"),
    ("s", "Careful, Leo doesn\u2019t just accept anyone \U0001f609"),
    ("s", "What makes you worthy? \U0001f352"),
    ("marzia", "Probably a bad influence with a nice face"),
    ("marzia", "I'm not sure what you mean. You asked me what I\u2019m looking for"),
    ("elina", "This profile is suspiciously wholesome"),
    ("elina", "What's the catch?"),
    ("elina", "Should I be worried or excited?"),
    ("aisha", "Oh I kind of wanted to date a bit and maybe along the way hopefully find something long term Wbu?"),
    ("amely-rose", "You look far too innocent for the amount of chaos in this profile"),
    ("angel", "What's your insta? So I can get this match organised"),
    ("andie", "Hey how have you been?"),
    ("alrissa", "Slow as in drinks this week? \U0001f440"),
    ("abby", "Gran\u2019s officially off duty \U0001f609"),
    ("abby", "Let me take you for a drink Thursday night. Does that work for you?"),
    ("arielle", "We should go for that coffee and walk sometime"),
    ("arielle", "I\u2019m taking you for that coffee and walk \U0001f352"),
    ("foley", "I have a feeling that shy side of you won\u2019t last very long over here \U0001f609"),
    ("foley", "That laugh gave you away before the question even did \U0001f352"),
    ("rae", "Come be mean then"),
    ("kathe", "Then prove it"),
    ("kathe", "So you only judge artists who can argue back \U0001f602"),
    ("kathe", "You tell me. I think you already know the answer \U0001f609"),
    ("kathe", "Careful, that almost sounds like a compliment \U0001f352"),
    ("kathe", "I\u2019ll help with the followers and handle pigeon security \U0001f609"),
    ("kathe", "That works"),
]

fails = []
for case, line in HOLDOUT:
    kill = critic(line, dict(CASE))
    if kill:
        fails.append((case, line[:50], kill))

print(f"hold-out: {len(HOLDOUT)} operator lines, {len(fails)} hard-failed")
for c, l, k in fails:
    print(f"  OVERFIT: [{c}] '{l}' killed by {k}")
if fails:
    sys.exit(1)
print("IDENTITY HOLDOUT PASS — no kill rejects the operator's own sent lines")
