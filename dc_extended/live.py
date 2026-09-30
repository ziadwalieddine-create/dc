"""dc_extended.live — running real threads. Outcomes are the product.

Protocol per turn: run_turn -> operator confirms -> confirm_send ->
her_reply(bubble, response_class). One data point then does three jobs:
reciprocity, read grading, and frame statistics with shrinkage from day one.
"""
import json
import os
from pathlib import Path
from .state import load_case, save_case, note
from .recovery import update_reciprocity
from .readlayer import grade_pending_read, push_window

GOOD = ("extends", "answers", "counters", "mirrors_length")
_STATS = Path(__file__).resolve().parents[1] / "frame_stats.json"


def load_stats(path=None):
    p = Path(path) if path else _STATS
    try:
        return json.loads(p.read_text())
    except Exception:
        return {}


def save_stats(stats, path=None):
    p = Path(path) if path else _STATS
    p.write_text(json.dumps(stats, indent=2) + "\n")


def record_frame_outcome(frame_id, response_class, stats_path=None):
    if not frame_id:
        return None
    stats = load_stats(stats_path)
    row = stats.setdefault(frame_id, {"n": 0})
    row["n"] += 1
    row[response_class] = row.get(response_class, 0) + 1
    save_stats(stats, stats_path)
    return row


def shrunk_rate(frame_id, frame_row, family_rate, k=5):
    """n < k: mostly the family prior. n large: mostly the frame's own rate.
    Five early rows can never crown a frame."""
    if not frame_row or frame_row.get("n", 0) == 0:
        return family_rate
    good = sum(frame_row.get(c, 0) for c in GOOD)
    rate = good / frame_row["n"]
    return (frame_row["n"] * rate + k * family_rate) / (frame_row["n"] + k)


def confirm_send(cases_dir, case_id, text=None):
    """The operator says it left the phone. Receipt discipline: the text must
    be the candidate, exactly (dc: actually_sent only from his confirmation)."""
    path = os.path.join(cases_dir, case_id + ".json")
    case = load_case(path)
    last = case.get("last_send")
    if not last or not last.get("text"):
        raise ValueError("no candidate to confirm")
    if text is not None and text.strip() != last["text"].strip():
        raise ValueError("confirmed text is not the candidate that left the pipeline")
    last["confirmed"] = True
    note(case, "confirmed: sent")
    save_case(case, path)
    return last


def her_reply(cases_dir, case_id, bubble, response_class,
              stats_path=None, read_stats_path=None):
    """One confirmed reply, three updates: reciprocity, read grading, frame stats.

    stats_path is REQUIRED: production passes the real path, tests pass a temp.
    A None default pointing at production is how stats get contaminated.
    """
    if stats_path is None:
        raise ValueError("stats_path required: pass a production path or a temp path for tests")
    if read_stats_path is None:
        raise ValueError("read_stats_path required: pass a production path or a temp path for tests")
    path = os.path.join(cases_dir, case_id + ".json")
    case = load_case(path)
    key = (bubble, response_class, (case.get("last_send") or {}).get("ts"))
    if list(case.get("last_reply_key") or []) == list(key):
        return {"duplicate": True, "reciprocity": case.get("reciprocity", 0),
                "graded": None, "frame": (case.get("last_send") or {}).get("frame_id")}
    case["last_reply_key"] = key
    case["her_last_bubble"] = bubble
    push_window(case, bubble)
    rec = update_reciprocity(case, response_class)
    graded = grade_pending_read(case, response_class, stats_path=read_stats_path)
    frame_id = (case.get("last_send") or {}).get("frame_id")
    if frame_id:
        record_frame_outcome(frame_id, response_class, stats_path=stats_path)
    note(case, "reply: %s rec=%s graded=%s frame=%s" % (response_class, rec, graded, frame_id))
    save_case(case, path)
    return {"reciprocity": rec, "graded": graded, "frame": frame_id}
