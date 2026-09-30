"""dc_extended — one turn: load, read, classify, write, gate, critic, simulate, select, emit."""
import os
from datetime import datetime
from .state import load_case, save_case, load_all, note, is_stale, parse_ts
from .gates import run_gates, day_dead, WEEKDAYS
from .queue import build_queue
from .ledger import record_send
from .simulate import simulate
from .variety import derive_seed, weighted_sample, is_synonym
from .readlayer import read_bubble, read_bubble_v2, needs_deep, push_window
from .goals import lock_endpoint, sim_rank, triage
from .recovery import record_deviation
from .live import load_stats as _load_frame_stats, shrunk_rate
import random
import re

from .lexicon import AUDIT_WORDS

# --- read layer lives in dc_extended/readlayer.py ----------------------------

# --- claim classifier --------------------------------------------------------

def classify(operator_text):
    """Minimal claim map. Model craft lives above this; code only tags."""
    if not operator_text:
        return None
    t = operator_text.lower()
    if any(w in t for w in AUDIT_WORDS):
        return {"type": "audit"}
    claim = None
    if any(w in t for w in ("wyd", "dinner", "drinks", "book", "pick the place", "saturday", "friday", "wednesday", "thursday")):
        claim = {"type": "plan", "family": "plan"}
    elif any(w in t for w in ("ig", "instagram", "dm ", "handle")):
        claim = {"type": "channel", "family": "dm", "channel": "ig"}
    elif any(w in t for w in ("mean", "tease", "basil")) or re.search(r"\bbit\b", t):
        claim = {"type": "bit", "family": "tease"}
    elif "hinge" in t:
        claim = {"type": "channel", "family": "new_opener", "channel": "hinge"}
    elif "follow" in t or "still" in t:
        claim = {"type": "answer", "family": "follow_up", "builds_on_last": True}
    if claim is not None and "wyd" in t:
        claim["family"] = "wyd_plan"
    return claim


# --- skeleton builder --------------------------------------------------------

def _build_skeleton(case, read, claim, gate_log, candidates, winner_idx, now, seed=None,
                    triage_label=None, corpus_path=None):
    """Populate the full-output contract from actual pipeline state."""
    skel = {
        "seed": seed,
        "triage": triage_label,
        "endpoint": (case.get("endpoint") or {}).get("type"),
        "case": {
            "name": case.get("name"),
            "stage": case.get("stage"),
            "channel": case.get("channel"),
            "goal": (case.get("goal") or {}).get("type"),
        },
        "read": read,
        "claim": claim,
        "gates": gate_log,
        "candidates": candidates,
        "select": None,
        "surface": None,
        "release": None,
    }
    if winner_idx is not None and winner_idx < len(candidates):
        w = candidates[winner_idx]
        if winner_idx is not None and candidates:
            w = candidates[winner_idx]
            skel["frame"] = {"id": w.get("frame_id"), "variant": w.get("variant")}
            if w.get("frame_id") and corpus_path:
                from .frames import card_for
                card = card_for(w["frame_id"], corpus_path)
                if card:
                    skel["frame"]["card"] = card
        skel["select"] = {
            "winner": winner_idx,
            "because": "best simulation, then heat fit, then amplitude fit, then family rate",
        }
        skel["surface"] = {
            "protected_tokens": "PASS",
            "reverify": "PASS" if not is_stale(case, now) else "FAIL",
        }
        skel["release"] = w["text"]
    return skel


# --- main turn ---------------------------------------------------------------

def run_turn(case_id, cases_dir, corpus_path, operator_text=None, now=None,
             her_bubble=None, operator_draft=None, rng_seed=None, temp=1.0,
             force_fast=False, endpoint=None, frame_stats_path=None, trigger_text=None):
    now = now or datetime.now()
    seed = rng_seed if rng_seed is not None else derive_seed(case_id, now.isoformat())
    rng = random.Random(seed)
    path = os.path.join(cases_dir, case_id + ".json")
    case = load_case(path)
    all_cases = load_all(cases_dir)
    queue = build_queue(all_cases, now)

    # staleness overrides everything
    if is_stale(case, now):
        return {"action": "stop", "why": "reverify",
                "detail": "state stale >48h; confirm what actually left the phone before any line",
                "queue": queue}

    # raw trigger text parses at the boundary; explicit kwargs win
    claim = None
    if trigger_text is not None:
        from .intake import open_turn
        parsed = open_turn(trigger_text, {"her_last_bubble": case.get("her_last_bubble"),
                                          "claim": case.get("claim") or {"type": None},
                                          "newer_bubble": False})
        if parsed.get("newer_bubble"):
            her_bubble = parsed["her_last_bubble"]
        if (parsed.get("claim") or {}).get("type"):
            case["claim"] = parsed["claim"]
            claim = parsed["claim"]

    # --- read layer ---
    bubble = her_bubble if her_bubble is not None else case.get("her_last_bubble")
    read = read_bubble_v2(bubble, case)
    case["last_read"] = read
    push_window(case, bubble)

    # refusal stops before any claim work
    if read["kind"] == "refuse":
        case["closed"] = True
        case["close_cause"] = "her_refusal"
        note(case, "read: refuse — thread closed by her boundary")
        save_case(case, path)
        return {"action": "stop", "gate": "stated_boundary",
                "detail": "she refused; closure is a decision with a cause",
                "queue": queue,
                "skeleton": _build_skeleton(case, read, None, [], [], None, now, seed=seed)}

    # stored claim from case file is a legal fallback (dc: stored claims survive)
    if claim is None and case.get("claim", {}).get("type"):
        stored = case["claim"]
        if stored.get("source") in ("stored", "operator_this_turn"):
            claim = dict(stored)
            claim["source"] = "stored"
    if claim is None:
        claim = classify(operator_text)

    # plan materialization: store the day-word fact; derive the clock per turn
    from .plans import materialize_plan, derive_plan_window
    materialize_plan(case, claim)
    derive_plan_window(case, now)

    # paste_without_send — the ceremony is not an action
    if claim and claim.get("type") == "audit":
        return {"action": "stop", "why": "paste_without_send",
                "detail": "analysis is not a send. one line out, or one close.",
                "queue": queue}

    # compute contract: a split read means the deep pass, not a guessed line.
    # Deep mode writes state (pending_read in the skeleton), never essays.
    if needs_deep(read) and not force_fast:
        case["pending_read"] = {"signal": bubble, "alternates": read["alternates"],
                                "primary": read["kind"]}
        note(case, "deep: split read %s vs %s" % (read["alternates"][0], read["alternates"][1]))
        save_case(case, path)
        return {"action": "deep", "why": "split_read", "read": read, "queue": queue,
                "detail": "two action-different readings; resolve_split() then re-run force_fast",
                "skeleton": _build_skeleton(case, read, claim, [], [], None, now, seed=seed)}

    family_override = None

    # --- endpoint lock: his goal this turn wins (substituted_end guard) ---
    ep = lock_endpoint(endpoint or operator_text, case.get("endpoint"))
    case["endpoint"] = ep

    # --- triage: derived, never stored. cooling spends the one-move budget ---
    tri = triage(case)
    if tri == "cooling" and not read.get("she_extended") and not force_fast:
        note(case, "hold: cooling — two flat outcomes; wait for her investment")
        save_case(case, path)
        return {"action": "stop", "why": "cooling", "triage": tri, "queue": queue,
                "detail": "neutral thread, two flat outcomes; his one clean move is spent",
                "skeleton": _build_skeleton(case, read, claim, [], [], None, now,
                                            seed=seed, triage_label=tri)}

    # --- gates ---
    case["_now"] = now.isoformat()
    case["_claim"] = claim
    gate_log = []
    gate_hits = run_gates(case, operator_draft or "", now, claim=claim,
                          all_cases=all_cases, log=gate_log)
    if gate_hits:
        gate, detail = gate_hits[0]
        record_deviation(case, "gate_stop", mechanism=gate,
                         family=(claim or {}).get("family"),
                         context={"bubble": bubble, "claim_type": (claim or {}).get("type"),
                                  "channel": case.get("channel"),
                                  "operator_text": operator_text,
                                  "banned_families": case.get("banned_families", []),
                                  "open_plan": case.get("open_plan"),
                                  "missed_windows": case.get("missed_windows", [])})
        # forced resolution: only these two rows may degrade; integrity never does
        degraded = None
        if gate == "dead_clock":
            from .recovery import forced_resolution as _fr
            degraded = _fr(case, gate)
        if degraded is None:
            note(case, "gate %s: %s" % (gate, detail))
            save_case(case, path)
            return {"action": "stop", "gate": gate, "detail": detail, "queue": queue,
                    "skeleton": _build_skeleton(case, read, claim, gate_log, [], None, now, seed=seed)}
        # degraded: fall through to the writer at the degraded amplitude/policy
        case["amp_cap"] = degraded["amplitude"]
        if degraded.get("policy") == "Play":
            family_override = "tease"
        note(case, "gate %s degraded -> %s at amplitude %d"
             % (gate, degraded["policy"], degraded["amplitude"]))


    # --- write: gather candidates ---
    family = family_override if family_override else (claim or {}).get("family", "follow_up")
    from .writer import _fill, _heat_fit, _amplitude_fit
    from .critic import critic
    from .ledger import load_corpus, family_stats, RANK

    rows = load_corpus(corpus_path)
    stats = family_stats(rows)
    target_heat = read.get("heat")
    target_amp = case.get("amplitude", 0)

    pool = [r for r in rows if r.get("family") == family
            and r.get("outcome") in ("extends", "counters", "answers", "mirrors_length")]
    candidates = []
    if pool:
        def score(r):
            oc = RANK.get(r.get("outcome"), 0)
            rate = stats.get(family, {}).get("rate", 0)
            hf = _heat_fit(r.get("heat"), target_heat)
            af = _amplitude_fit(r.get("amplitude"), target_amp)
            return (oc + rate * 3 - hf * 2 - af, r.get("ts", ""))
        pool.sort(key=score, reverse=True)
        fstats_data = _load_frame_stats(frame_stats_path) if frame_stats_path else {}
        cap = case.get("amp_cap")
        for row in pool[:3]:  # top 3 candidates
            if cap is not None and (row.get("amplitude") or 0) > cap:
                continue
            variants = row["text"] if isinstance(row["text"], list) else [row["text"]]
            vi = rng.randrange(len(variants)) if len(variants) > 1 else 0
            text = _fill(variants[vi], case)
            case["_policy"] = {"tease": "Play", "callback": "Play", "follow_up": "Play",
                               "get_to_know": "Trade", "plan": "Bridge"}.get(family, "Play")
            case["_amplitude"] = row.get("amplitude") or 0
            sim_classes = simulate(text, case, read, claim, now)
            kill = critic(text, case, read=read, sim_classes=set(sim_classes),
                          other_wider=(len(pool) > 1))
            base_rate = stats.get(family, {}).get("rate", 0)
            fr = fstats_data.get(row.get("id")) if row.get("frame") else None
            candidates.append({
                "text": text,
                "source": "retrieved",
                "frame_id": row.get("id") if row.get("frame") else None,
                "variant": vi,
                "family_rate": shrunk_rate(row.get("id"), fr, base_rate),
                "heat_fit": _heat_fit(row.get("heat"), target_heat),
                "amp_fit": _amplitude_fit(row.get("amplitude"), target_amp),
                "critic": kill or "PASS",
                "simulate": sorted(sim_classes) if sim_classes else [],
            })
            if kill is not None:
                killed_text = candidates[-1]["text"]
                pool = [r for r in pool if not any(
                    is_synonym(killed_text, v)
                    for v in (r["text"] if isinstance(r["text"], list) else [r["text"]]))]

    if operator_draft:
        text = operator_draft.strip()
        kill = critic(text, case, read=read)
        sim_classes = simulate(text, case, read, claim, now) if kill is None else set()
        candidates.append({
            "text": text,
            "source": "edited",
            "family_rate": 0,
            "heat_fit": 0,
            "amp_fit": 0,
            "critic": kill or "PASS",
            "simulate": sorted(sim_classes) if sim_classes else [],
        })

    # filter to survivors
    survivors = [c for c in candidates if c["critic"] == "PASS"]
    if not survivors:
        best_kill = candidates[0]["critic"] if candidates else "no_candidates"
        note(case, "writer: no legal line (kill=%s)" % best_kill)
        save_case(case, path)
        return {"action": "stop", "why": best_kill or "no_legal_line",
                "detail": "no corpus match, no draft, or critic killed all candidates",
                "queue": queue,
                "skeleton": _build_skeleton(case, read, claim, gate_log, candidates, None, now, seed=seed)}

    # --- select: best simulation, then heat fit, then amplitude fit, then family rate ---
    SIM_RANK = sim_rank(ep["type"])
    def select_score(c):
        sim_score = max([SIM_RANK.get(s, 0) for s in c["simulate"]], default=0)
        return (sim_score, -c["heat_fit"], -c["amp_fit"], c["family_rate"])

    survivors.sort(key=select_score, reverse=True)

    def _numeric(c):
        s = select_score(c)
        return s[0] * 100 - s[1] * 10 - s[2] + s[3]

    idx = weighted_sample(survivors, [_numeric(c) for c in survivors], rng, temp)
    winner = survivors[idx]
    winner_idx = candidates.index(winner)

    # --- write-forward + emit as one transaction ---
    text = winner["text"]
    case["last_send"] = {"text": text, "ts": now.isoformat(timespec="seconds"),
                         "confirmed": False, "family": family,
                         "source": winner["source"], "frame_id": winner.get("frame_id")}
    case["her_last_bubble"] = bubble
    note(case, "send: %s [%s] read=%s heat=%s" % (text[:60], winner["source"], read["kind"], read["heat"]))
    save_case(case, path)
    record_send(corpus_path, case_id, family, text, now.isoformat(timespec="seconds"))

    skeleton = _build_skeleton(case, read, claim, gate_log, candidates, winner_idx, now, seed=seed, triage_label=tri, corpus_path=corpus_path)

    return {"action": "send", "text": text, "claim": (claim or {}).get("type"),
            "family": family, "source": winner["source"], "read": read["kind"],
            "heat": read["heat"], "queue": queue, "skeleton": skeleton,
            "frame": skeleton.get("frame")}
