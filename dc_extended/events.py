"""dc_extended.events — the immutable truth stream.

Append-only event log + deterministic reducer. SHADOW MODE: events are recorded
alongside the legacy case JSON; the legacy files remain the runtime authority
until the migration proves this projection. There is no mutation API. There is
no way to rewrite history — no update, no delete, no in-place outcome write.

Event identity chain: event -> proposal -> confirmed send -> reply -> outcome.
A proposal is never a send. Raw evidence is never overwritten by interpretation.
"""
import json, os, hashlib, uuid
from datetime import datetime

EVENT_TYPES = ("incoming_message", "operator_directive", "candidate_proposed",
               "send_confirmed", "incoming_reply", "outcome_interpretation",
               "correction", "channel_change", "boundary_change", "date_object_change")
ACTORS = ("her", "operator", "engine", "audit")


def _eid(prefix):
    return "%s_%s" % (prefix, uuid.uuid4().hex[:12])


def _now():
    return datetime.now().isoformat(timespec="seconds")


def _hash(ev):
    d = {k: v for k, v in ev.items() if k != "hash"}
    return hashlib.sha256(json.dumps(d, sort_keys=True).encode()).hexdigest()[:16]


class EventLog:
    """Append-only, hash-chained. The ONLY write API is append()."""

    def __init__(self, path):
        self.path = path

    def _last_hash(self):
        try:
            lines = open(self.path).read().strip().splitlines()
            return json.loads(lines[-1])["hash"] if lines else None
        except FileNotFoundError:
            return None

    def append(self, type_, actor, case_id, payload, occurred_at=None):
        if type_ not in EVENT_TYPES:
            raise ValueError("unknown event type: %r" % type_)
        if actor not in ACTORS:
            raise ValueError("unknown actor: %r" % actor)
        ev = {"event_id": _eid("ev"), "type": type_, "actor": actor,
              "case_id": case_id, "payload": payload,
              "occurred_at": occurred_at or _now(), "recorded_at": _now(),
              "prev_hash": self._last_hash(), "hash": None}
        ev["hash"] = _hash(ev)
        with open(self.path, "a") as f:
            f.write(json.dumps(ev) + "\n")
        return ev["event_id"]

    def read(self):
        try:
            return [json.loads(l) for l in open(self.path).read().strip().splitlines()]
        except FileNotFoundError:
            return []

    def verify(self):
        """Hash chain integrity. If this returns False, history was tampered with."""
        prev = None
        for ev in self.read():
            if ev["prev_hash"] != prev or _hash(ev) != ev["hash"]:
                return False
            prev = ev["hash"]
        return True

    # ---- causal chain helpers: proposal is never a send ----
    def propose(self, case_id, rendered_payload, **kw):
        return self.append("candidate_proposed", "engine", case_id,
                           {"rendered_payload": rendered_payload, **kw})

    def _find(self, event_id, type_=None):
        for e in self.read():
            if e["event_id"] == event_id and (type_ is None or e["type"] == type_):
                return e
        return None

    def confirm_send(self, case_id, proposal_id, exact_payload, sent_at=None):
        """A ConfirmedSend REQUIRES an existing proposal IN THE SAME CASE, with an
        exact payload match. Cross-case confirmation is structurally impossible."""
        prop = self._find(proposal_id, "candidate_proposed")
        if prop is None:
            raise ValueError("no such proposal: %r" % proposal_id)
        if prop["case_id"] != case_id:
            raise ValueError("cross-case confirmation: proposal belongs to %r, not %r"
                             % (prop["case_id"], case_id))
        if prop["payload"]["rendered_payload"].strip() != exact_payload.strip():
            raise ValueError("confirmed payload does not match the proposed payload")
        return self.append("send_confirmed", "operator", case_id,
                           {"proposal_id": proposal_id, "exact_payload": exact_payload,
                            "actually_sent_at": sent_at or _now()})

    def interpret_outcome(self, case_id, reply_id, classification, source="engine"):
        """Interpretation is a DERIVED event. Raw reply stays immutable; a later
        correction appends another interpretation — it never edits this one."""
        rep = self._find(reply_id, "incoming_reply")
        if rep is None:
            raise ValueError("no such reply: %r" % reply_id)
        if rep["case_id"] != case_id:
            raise ValueError("cross-case outcome: reply belongs to %r" % rep["case_id"])
        return self.append("outcome_interpretation", source, case_id,
                           {"reply_id": reply_id, "classification": classification})


def project(log, case_id=None):
    """Deterministic reducer: event stream -> factual read model. Per-case:
    projections never merge different women. case_id=None projects all (admin)."""
    st = {"closed": False, "last_incoming": None, "channel": None,
          "confirmed_sends": [], "open_proposals": {}, "interpretations": []}
    for ev in log.read():
        if case_id is not None and ev["case_id"] != case_id:
            continue
        p, t = ev["payload"], ev["type"]
        if t == "incoming_message":
            st["last_incoming"] = p
        elif t == "operator_directive" and p.get("directive") in ("stop", "close"):
            st["closed"] = True
        elif t == "candidate_proposed":
            st["open_proposals"][ev["event_id"]] = p
        elif t == "send_confirmed":
            st["confirmed_sends"].append(p)
            st["open_proposals"].pop(p.get("proposal_id"), None)
        elif t == "channel_change":
            st["channel"] = p.get("to")
        elif t == "outcome_interpretation":
            st["interpretations"].append(p)
    return st
