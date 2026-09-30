"""Fix 6b: every case conforms to template + enumerated runtime-allowed keys."""
import sys, os, json
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

SRC = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "cases")
RUNTIME_ALLOWED = {"last_read", "last_send", "her_last_bubble", "chronology", "outcomes",
                   "deviations", "reciprocity", "amp_cap", "endpoint", "her_recent_bubbles",
                   "pending_read", "last_split", "last_read_override", "frame_limits",
                   "newer_bubble", "close_cause", "last_candidate", "last_reply_key", "red_flags",
                   "banned", "baseline", "claim", "handle"}


def t(name, cond, detail=""):
    print(("PASS " if cond else "FAIL ") + name + (" | " + str(detail) if detail and not cond else ""))
    if not cond:
        sys.exit(1)


tmpl = json.load(open(os.path.join(SRC, "_template.json")))
tmpl_keys = set(tmpl.keys())
allowed = tmpl_keys | RUNTIME_ALLOWED

cases = [f for f in os.listdir(SRC) if f.endswith(".json") and not f.startswith("_")]
t("fixtures exist", len(cases) >= 1, cases)
for fn in cases:
    c = json.load(open(os.path.join(SRC, fn)))
    extra = set(c.keys()) - allowed
    t(fn + " keys conform", not extra, sorted(extra))
    missing = [k for k in ("id", "closed", "referents", "chronology") if k not in c]
    t(fn + " required keys", not missing, missing)
    for k in ("reciprocity", "amp_cap", "deviations", "endpoint", "outcomes"):
        if k not in c:
            print("  note: %s lacks %s (runtime default handles it)" % (fn, k))

# no underscore keys persisted in any fixture
for fn in cases:
    c = json.load(open(os.path.join(SRC, fn)))
    us = [k for k in c if k.startswith("_")]
    t(fn + " no transient keys", not us, us)

print("ALL SCHEMA TESTS PASS")
