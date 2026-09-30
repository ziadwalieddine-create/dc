"""dc_extended.state — case schema, atomic load/save, chronology."""
import json, os, tempfile
from datetime import datetime, timedelta

REQUIRED = ["id","name","closed","channel","goal","updated"]

def parse_ts(s):
    return datetime.fromisoformat(s)

def load_case(path):
    with open(path) as f:
        case = json.load(f)
    for k in REQUIRED:
        if k not in case:
            raise ValueError("case missing field: %s" % k)
    return case

def save_case(case, path):
    """Atomic write. Emit without this completing = failed run.
    Underscore keys are transient (derived per turn) and never persist."""
    case = {k: v for k, v in case.items() if not k.startswith("_")}
    case["updated"] = datetime.now().isoformat(timespec="seconds")
    fd, tmp = tempfile.mkstemp(dir=os.path.dirname(path), suffix=".tmp")
    with os.fdopen(fd, "w") as f:
        json.dump(case, f, indent=2)
    os.replace(tmp, path)

def note(case, text, ts=None):
    ts = ts or datetime.now().isoformat(timespec="seconds")
    case.setdefault("chronology", []).append({"ts": ts, "note": text})

def is_stale(case, now, hours=48):
    return (now - parse_ts(case["updated"])) > timedelta(hours=hours)

def load_all(cases_dir):
    out = []
    for fn in sorted(os.listdir(cases_dir)):
        if fn.endswith(".json") and not fn.startswith("_"):
            out.append(load_case(os.path.join(cases_dir, fn)))
    return out
