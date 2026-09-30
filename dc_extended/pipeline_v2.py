"""dc_extended.pipeline_v2 — the composed new path (shadow).

One turn end-to-end: ingress -> read_v2 -> licences -> amplitude -> Move ->
state delta -> events -> final audit. Runs BESIDE the legacy run_turn; does not
replace it. The acceptance constitution's invariants are checked against THIS
path once composed.
"""
import os
from .ingress import parse_operator, split_trigger
from .read_v2 import read_v2
from .licences import Licences
from .amplitude import permitted_band
from .moves import single, sequence
from .dateobjects import StateDelta, final_audit, resolve_day, Plan
from .events import EventLog
from datetime import datetime


def turn_v2(case, her_bubble, operator_text, events_path, now=None, other_cases=None, corpus_path=None):
    """Serve one turn on the new substrate. Returns a full trace."""
    now = now or datetime.now()
    log = EventLog(events_path)

    # 1. ingress: actors separated before interpretation
    her, op = split_trigger(operator_text) if operator_text and "she:" in operator_text else (None, operator_text)
    directive = parse_operator(op)

    # 2. HARD GATE: operator STOP dominates everything
    if directive.stop:
        log.append("operator_directive", "operator", case["id"],
                   {"text": op, "directive": "stop"})
        return {"action": "stop", "why": "explicit_stop", "closed": True}

    # 3. record her evidence
    if her_bubble:
        log.append("incoming_message", "her", case["id"], {"text": her_bubble})

    # 4. multidimensional read (legacy kind preserved as projection)
    read = read_v2(her_bubble, case)

    # 5. licences + amplitude band
    lic = Licences(read, case)
    floor, ceiling, recommended = permitted_band(read, lic)

    # 6. O* governs legality: endpoint/channel intent constrains the move
    if directive.endpoint == "close" or directive.claim_type == "close":
        return {"action": "stop", "why": "endpoint_close", "closed": True}
    if directive.endpoint == "date_lock" and read.get("heat") == 2 and not lic.escalation_licensed:
        pass  # heat present but not licensed -> fall through to restrained generation
    if directive.endpoint == "off_app" or directive.claim_type == "channel":
        target = directive.channel or "ig"
        return {"action": "send", "move": sequence(
            ["you're more fun over here — what's your %s?" % target],
            family="channel", heat=0, obligations=[]),
            "read": read, "band": (floor, ceiling, recommended),
            "shadow_proposal_id": None, "delta_committed": False}
    if directive.endpoint == "date_lock" or directive.claim_type == "plan":
        plan_frame = _retrieve(dict(case, _force_family="plan"), read, directive, corpus_path)
        if plan_frame:
            pv = plan_frame["text"] if isinstance(plan_frame["text"], list) else [plan_frame["text"]]
            return {"action": "send", "move": single(pv[0], family="plan",
                    heat=0, obligations=["coordinate_plan"]),
                    "read": read, "band": (floor, ceiling, recommended),
                    "shadow_proposal_id": None, "delta_committed": False}
        return {"action": "send", "move": single(
            case.get("plan_line") or "drinks this week — I\u2019ll pick the spot, you pick the day",
            family="plan", heat=0, obligations=["coordinate_plan"]),
            "read": read, "band": (floor, ceiling, recommended),
            "shadow_proposal_id": None, "delta_committed": False}

    # 7. generation: prefer the existing certified frame corpus; fall back to a
    # minimal template that satisfies obligations. No canned third writer.
    obligations = list(read.get("response_obligations") or [])
    frame = _retrieve(case, read, directive, corpus_path=corpus_path)
    if frame:
        variants = frame["text"] if isinstance(frame["text"], list) else [frame["text"]]
        idx = (hash(str(her_bubble)) % len(variants)) if len(variants) > 1 else 0
        bubbles = [variants[idx]]
        if "answer_location" in obligations and len(bubbles) < 2:
            bubbles.append(case.get("location_fact") or "I'll send you the pin")
    elif "answer_location" in obligations:
        bubbles = [case.get("location_fact") or "Right outside the station"]
        if read.get("heat") and read.get("heat") >= 1:
            bubbles.append(case.get("playful_continuation") or "You'll find me")
    elif "answer_question" in obligations:
        bubbles = [case.get("answer_fact") or "good question — let me check"]
    elif read.get("heat") == 2 and lic.escalation_licensed:
        bubbles = [case.get("escalation_line") or "say that to my face"]
    elif (read.get("heat") or 0) >= 1 and read.get("relational_functions"):
        bubbles = [case.get("warm_line") or "that tracks"]
    else:
        return {"action": "stop", "why": "no_legal_move", "read": read}

    move = sequence(bubbles, family="composed", heat=min(read.get("heat") or 0, ceiling),
                        obligations=obligations) if len(bubbles) > 1 else single(
                        bubbles[0], family="composed", heat=min(read.get("heat") or 0, ceiling),
                        obligations=obligations)

    # 7. final audit on the EXACT rendered payload BEFORE any state mutation
    audit_ok, audit_fails = final_audit(move["sequence"], case, other_cases or [])
    if not audit_ok:
        return {"action": "stop", "why": "final_audit", "failures": audit_fails,
                "read": read, "case_unchanged": True}

    # 8. propose -> validate -> commit durable state (only after legality)
    delta = StateDelta({"last_read_heat": read.get("heat")})
    ok, reasons = delta.commit(case, gates=[])

    # 9. record proposal in the truth stream
    pid = log.propose(case["id"], " / ".join(move["sequence"]), source="pipeline_v2")
    return {"action": "send", "move": move, "read": read, "licences": {k: v for k, v in vars(lic).items()
            if isinstance(v, bool)}, "band": (floor, ceiling, recommended),
            "shadow_proposal_id": pid, "delta_committed": ok}


def _retrieve(case, read, directive, corpus_path=None):
    """Serve from the certified frame corpus by family/heat. Returns a frame dict
    or None. (Corpus wiring lands with the single-writer pass; template fallback
    above keeps the pipeline honest meanwhile.)"""
    if corpus_path and os.path.exists(corpus_path):
        import json as _json
        for line in open(corpus_path):
            fr = _json.loads(line)
            if not fr.get("frame"):
                continue
            if case.get("_force_family") == "plan":
                if fr.get("family") == "plan" and fr.get("outcome") in ("extends","answers","counters","mirrors_length"):
                    return fr
                continue   # force-family: skip everything else
            if read.get("heat") == 2 and fr.get("heat") == 2 and fr.get("family") in ("tease", "callback"):
                return fr
    return None
