"""Tranche 4: semantic licences + heat->amplitude. Every licence paired:
licensed context permits, unlicensed context forbids — same wording both ways."""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from dc_extended.read_v2 import read_v2
from dc_extended.licences import Licences
from dc_extended.amplitude import permitted_band, restore_band

FAIL = []
def t(name, cond, detail=""):
    print(("PASS " if cond else "FAIL ") + name + (" | " + str(detail) if detail and not cond else ""))
    if not cond: FAIL.append(name)

CASE_W = {"referents": [{"id": "w", "text": "Wynyard Park", "status": "live", "type": "plan",
                         "introduced_by": "him", "touched_by_her": True}]}

# --- reciprocal challenge: licensed by HER challenge, never by the phrase ---
r = read_v2("that type of game is mine, you need to be more original", {})
lc = Licences(r, {"last_her_challenge": "that type of game is mine, you need to be more original"})
t("kathe: her explicit challenge licenses challenge", lc.permits("challenge") is True)
r = read_v2("what are you up to", {})
lc = Licences(r, {})
t("NEG: cold 'prove it' with NO challenge frame is unlicensed", lc.permits("challenge") is False)
r = read_v2("haha", {})
t("NEG: same wording, no licence — still forbidden", lc.permits("challenge") is False)

# --- qualification: licensed ONLY by her deserve-frame ---
r = read_v2("Are you sure you deserve to read it?", {})
lc = Licences(r, {"last_her_challenge": "Are you sure you deserve to read it?"})
t("kathe: her deserve-frame licenses qualification", lc.permits("qualification") is True)
r = read_v2("I replied at 2am", {})
lc = Licences(r, {})
t("NEG: grading her ordinary reply is unlicensed (the '2am points' failure)",
  lc.permits("qualification") is False)

# --- misread: evidence + reversibility ---
r = read_v2("your comment was the only good one", {})
lc = Licences(r, {})
t("kathe: her genuine compliment licenses misread", lc.permits("misread") is True)
r = read_v2("I work in the city", {})
lc = Licences(r, {})
t("NEG: neutral disclosure licenses no misread", lc.permits("misread") is False)

# --- escalation: attraction cue or heat 2 ---
r = read_v2("Hot", {})
lc = Licences(r, {})
t("gabby: 'Hot' licenses escalation", lc.permits("escalation") is True)
r = read_v2("the basil plant is nice", CASE_W)
lc = Licences(r, {})
t("NEG: polite touch does not license escalation", lc.permits("escalation") is False)

# --- heat -> amplitude: causal, not monotonic, logistics never penalizes ---
r = read_v2("Hot", {})
f0 = permitted_band({"heat": 0}, Licences({}))
f2 = permitted_band(r, Licences(r))
t("causal: heat 2 band wider than heat 0", f2[1] > f0[1], (f0, f2))
rW = read_v2("Where the heck is Wynyard Park 🤣 😂", CASE_W)
fW = permitted_band(rW, Licences(rW), literal_obligations=["answer_location"])
t("non-regression: logistics does NOT reduce the heat-2 ceiling", fW[1] == 2, fW)
t("not monotonic: obligation recommends inside the band, ceiling intact",
  fW[2] <= fW[1] and fW[1] == 2, fW)
rN = read_v2("Where is the station?", {})
fN = permitted_band(rN, Licences(rN))
t("NEG: neutral logistics stays at band 0", fN[1] == 0, fN)

# --- reversibility: later heat lifts the dry-streak cap ---
case = {"amp_cap": 0}
r = read_v2("Hot", {})
t("reversible: genuine heat restores the band", restore_band(case, r) is True and "amp_cap" not in case)
case2 = {"amp_cap": 0}
t("NEG: no heat, cap stays", restore_band(case2, read_v2("haha", {})) is False and case2["amp_cap"] == 0)

print()
print("ALL TRANCHE-4 TESTS PASS" if not FAIL else "FAILURES: %s" % FAIL)
assert not FAIL, FAIL
