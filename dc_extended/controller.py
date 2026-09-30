"""dc_extended.controller — the single release authority.

01-date-controller owns final selection/release. 02-date-state (events.py) owns
factual truth. Specialists advise; the controller decides. During migration the
legacy run_turn remains available as the reference path — the controller wraps
it and records the shadow event stream. Nothing is released without the
controller, and the controller releases nothing without an auditable proposal.
"""
import os
from .events import EventLog
from .ingress import parse_operator, split_trigger


class Controller:
    """The only path to a released move. Shadow mode: every turn records events;
    legacy state still executes."""

    def __init__(self, events_path):
        self.log = EventLog(events_path)

    def run_shadow(self, case_id, run_turn_fn, *args, **kw):
        """Execute the legacy path; record the truth stream beside it."""
        her = kw.get("her_bubble")
        if her is None:
            for a in args:
                if isinstance(a, str) and a == case_id:
                    continue   # positional case_id is not her bubble
                if isinstance(a, str) and not os.path.isdir(str(a)) and not str(a).endswith((".jsonl",".json")) \
                   and not str(a).endswith("cases") and "cases" not in str(a) \
                   and len(str(a)) < 200 and " " in str(a):
                    her = a; break
        if her:
            self.log.append("incoming_message", "her", case_id, {"text": her})
        # operator text: explicit kwarg, else the first short free-text positional
        op = kw.get("operator_text")
        if op is None:
            import datetime as _dt
            for a in args:
                if not isinstance(a, str) or a == case_id:
                    continue
                if os.path.isdir(str(a)) or str(a).endswith((".jsonl",".json")) or "cases" in str(a):
                    continue
                try:
                    _dt.datetime.fromisoformat(str(a)); continue   # a 'now' timestamp
                except ValueError:
                    pass
                if len(str(a)) < 200:
                    op = a; break
        trig = kw.get("trigger_text")
        her2, op2 = split_trigger(trig) if trig else (None, None)
        if her2 and not her:
            self.log.append("incoming_message", "her", case_id, {"text": her2})
        d = parse_operator(op2 or op)
        if op2 or op:
            self.log.append("operator_directive", "operator", case_id,
                            {"text": op2 or op, "directive": "stop" if d.stop else None,
                             "channel": d.channel, "endpoint": d.endpoint})
        result = run_turn_fn(*args, **kw)
        if result.get("action") == "send" and result.get("text"):
            pid = self.log.propose(case_id, result["text"],
                                   source=result.get("source"))
            result["shadow_proposal_id"] = pid
        return result


def _directive(text):
    t = (text or "").lower()
    if t.strip() in ("stop", "close", "leave it"):
        return "stop"
    return None


def _run_primary(self, case, her_bubble=None, operator_text=None, now=None, other_cases=None):
    """PRIMARY PATH: the v2 substrate. Legacy run_shadow remains the reference path."""
    from .pipeline_v2 import turn_v2
    return turn_v2(case, her_bubble, operator_text, self.log.path, now=now,
                   other_cases=other_cases)


Controller.run_primary = _run_primary
