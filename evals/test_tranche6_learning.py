"""Tranche 6: learning on the identity chain."""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from dc_extended.learning import (HistoricalMove, derive_outcome, promote_to_frame,
                                  certify_portability, OUTCOME_DIMENSIONS)

FAIL = []
def t(name, cond, detail=""):
    print(("PASS " if cond else "FAIL ") + name + (" | " + str(detail) if detail and not cond else ""))
    if not cond: FAIL.append(name)

# --- raw evidence dominates labels ---
d = derive_outcome("leave me alone", reply_class="extends")
t("raw beats label: refusal evidence overrides supplied label",
  d["attraction"] == "negative" and d["investment"] == "negative", d["reply"])
t("supplied label recorded as input, never classification",
  d.get("label_supplied") == "extends" and d["reply"] in ("unclassified", None), d["reply"])
d = derive_outcome("Yes! Wednesday works, can't wait 😂")
t("dimensional: date progression positive", d["date_progression"] == "positive")
t("dimensional: investment positive", d["investment"] == "positive")
t("eight dimensions, not one score", len(OUTCOME_DIMENSIONS) == 8
  and all(k in d for k in OUTCOME_DIMENSIONS))

# --- one reply, one outcome regardless of bubble count ---
d = derive_outcome("omg yes!! \U0001f602 and also where are we going? and what time?")
t("multi-bubble reply = ONE outcome", d["date_progression"] == "positive")

# --- HistoricalMove: independent properties coexist ---
h = HistoricalMove("m1", "Where are you enjoying the sun?", "gabby",
                   actually_sent=True, operator_later_rejected=True,
                   verification_status="operator_corrected")
t("sent + later-rejected + corrected coexist, no forced label",
  h.as_dict()["actually_sent"] is True and h.as_dict()["operator_later_rejected"] is True)

# --- promotion: certification gate, not automatic ---
good = HistoricalMove("m2", "Say that to my face", "gabby", actually_sent=True,
                      verification_status="operator_corrected")
fr, reasons = promote_to_frame(good, other_case_names=["katherine"])
t("portable move promotes with certification", fr is not None and fr["certified"] is True)
leaky = HistoricalMove("m3", "Meet me at Wynyard Park at 5:50", "gabby",
                       actually_sent=True, verification_status="operator_corrected")
fr2, reasons2 = promote_to_frame(leaky, other_case_names=["katherine"])
t("NEG: case-specific token blocks promotion", fr2 is None and any("wynyard" in r for r in reasons2), reasons2)
unverified = HistoricalMove("m4", "some line", "x", actually_sent=True,
                            verification_status="pending")
ok3, r3 = certify_portability(unverified)
t("NEG: unverified material cannot promote", not ok3)
notsent = HistoricalMove("m5", "some line", "x", actually_sent=False,
                         verification_status="verified")
ok4, r4 = certify_portability(notsent)
t("NEG: unsent proposal never promotes", not ok4)

# --- historical record survives rejection (negative evidence kept) ---
poor = HistoricalMove("m6", "What changed?", "aditi", actually_sent=True,
                      verification_status="verified", style_grade="poor")
t("poor sent line stays historical (not deleted)", poor.as_dict()["actually_sent"] is True
  and poor.as_dict()["style_grade"] == "poor")

print()
print("ALL TRANCHE-6 TESTS PASS" if not FAIL else "FAILURES: %s" % FAIL)
assert not FAIL, FAIL
