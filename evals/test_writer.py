"""Writer tests — retrieval ranking, critic screening, heat fit, fallback."""
import sys, os, tempfile, json
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from datetime import datetime
from dc_extended.writer import draft
from dc_extended.ledger import record_send, record_outcome, load_corpus

NOW = datetime(2026, 9, 18, 18, 11)


def base_case(**kw):
    c = {"id": "t", "name": "Rae", "stage": "off_app", "channel": "ig",
         "goal": {"type": "fwb"}, "updated": "2026-09-19T10:00:00",
         "banned_families": [], "boundaries": [], "unique_facts": ["basil plant"],
         "missed_windows": [], "referents": [
             {"id": "b", "text": "basil plant", "type": "bit",
              "introduced_by": "him", "touched_by_her": True, "status": "live"}],
         "frame_limits": [], "banned": [],
         "_now": NOW, "_claim": {"type": "bit", "family": "tease"},
         "ig_handle": "r.tmzon", "amplitude": 1}
    c.update(kw)
    return c


def t(name, cond):
    print(("PASS " if cond else "FAIL ") + name)
    if not cond:
        sys.exit(1)


def make_corpus(rows):
    fd, path = tempfile.mkstemp(suffix=".jsonl")
    os.close(fd)
    with open(path, "w") as f:
        for r in rows:
            f.write(json.dumps(r) + chr(10))
    return path


# --- retrieval + parameter fill ---
corpus = make_corpus([
    {"case": "rae", "family": "tease", "text": "the {referent} is doing too much",
     "outcome": "extends", "ts": "2026-09-17T21:00:00", "heat": 2, "amplitude": 1},
])
text, source, kill = draft(base_case(), "tease", corpus)
t("retrieved + filled", source == "retrieved" and "basil plant" in text and kill is None)

# --- heat fit: heat-2 row beats heat-0 row for a heat-2 read ---
corpus2 = make_corpus([
    {"case": "rae", "family": "tease", "text": "ok cool", "outcome": "extends",
     "ts": "2026-09-17T21:00:00", "heat": 0, "amplitude": 0},
    {"case": "rae", "family": "tease", "text": "the {referent} is doing too much",
     "outcome": "extends", "ts": "2026-09-16T21:00:00", "heat": 2, "amplitude": 1},
])
read = {"kind": "tease", "heat": 2, "she_extended": True}
text, source, kill = draft(base_case(), "tease", corpus2, read=read)
t("heat fit ranks", "basil plant" in text)

# --- critic screens: generic line killed, falls to next ---
corpus3 = make_corpus([
    {"case": "rae", "family": "tease", "text": "you're gorgeous", "outcome": "extends",
     "ts": "2026-09-17T21:00:00", "heat": 2, "amplitude": 1},
    {"case": "rae", "family": "tease", "text": "the {referent} is doing too much",
     "outcome": "extends", "ts": "2026-09-16T21:00:00", "heat": 2, "amplitude": 1},
])
text, source, kill = draft(base_case(), "tease", corpus3, read=read)
t("critic screens generic", "basil plant" in text and kill is None)

# --- all corpus killed, operator draft fallback ---
corpus4 = make_corpus([
    {"case": "rae", "family": "tease", "text": "you're gorgeous", "outcome": "extends",
     "ts": "2026-09-17T21:00:00", "heat": 2, "amplitude": 1},
])
text, source, kill = draft(base_case(), "tease", corpus4, read=read,
                           operator_draft="the basil plant is unhinged")
t("operator draft fallback", source == "edited" and "unhinged" in text)

# --- no corpus, no draft, no cold invention ---
text, source, kill = draft(base_case(), "tease", make_corpus([]), read=read)
t("no cold invention", text is None and source == "none")

# --- outcome recording ---
fd, cp = tempfile.mkstemp(suffix=".jsonl")
os.close(fd)
record_send(cp, "rae", "tease", "line one", "2026-09-17T21:00:00")
record_outcome(cp, "rae", "line one", "extends", "2026-09-17T22:00:00")
rows = load_corpus(cp)
t("outcome recorded", rows[0]["outcome"] == "extends")

print("ALL WRITER TESTS PASS")
