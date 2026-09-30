"""dc_extended.writer — retrieval first, heat-aware, critic-screened.

Rank corpus by (family rate, heat fit, amplitude fit, recency).
Screen through critic. Fallback: operator draft edit. Never cold invention.
"""
from .ledger import load_corpus, RANK, family_stats
from .critic import critic
from .variety import weighted_sample

PARAMS = ("{name}", "{referent}", "{day}", "{handle}")


def _fill(template, case):
    refs = case.get("referents", [])
    live = next((r["text"] for r in refs if r.get("status") == "live"), None)
    out = template
    out = out.replace("{name}", case.get("name", ""))
    out = out.replace("{referent}", live or case.get("name", ""))
    out = out.replace("{handle}", case.get("ig_handle") or "")
    plan = case.get("open_plan") or {}
    out = out.replace("{day}", plan.get("when_day") or "")
    return out


def _heat_fit(row_heat, target_heat):
    """0 = exact match, 1 = adjacent, 2 = wrong tier."""
    if row_heat is None or target_heat is None:
        return 1
    return abs(row_heat - target_heat)


def _amplitude_fit(row_amp, target_amp):
    if row_amp is None or target_amp is None:
        return 0
    return abs(row_amp - target_amp)


def draft(case, family, corpus_path, operator_draft=None, read=None, now=None,
          rng=None, temp=1.0):
    """Returns (text, source, kill). source in {retrieved, edited, none}.
    kill is the critic id that stopped the best candidate, or None."""
    rows = load_corpus(corpus_path)
    stats = family_stats(rows)
    target_heat = (read or {}).get("heat")
    target_amp = case.get("amplitude", 0)

    pool = [r for r in rows if r.get("family") == family
            and r.get("outcome") in ("extends", "counters", "answers", "mirrors_length")]
    if pool:
        def score(r):
            oc = RANK.get(r.get("outcome"), 0)
            rate = stats.get(family, {}).get("rate", 0)
            hf = _heat_fit(r.get("heat"), target_heat)
            af = _amplitude_fit(r.get("amplitude"), target_amp)
            return (oc + rate * 3 - hf * 2 - af, r.get("ts", ""))

        pool.sort(key=score, reverse=True)
        passing = []
        for row in pool:
            _v = row["text"] if isinstance(row["text"], list) else [row["text"]]
            raw = _v[rng.randrange(len(_v))] if (rng and len(_v) > 1) else _v[0]
            candidate = _fill(raw, case)
            kill = critic(candidate, case, read=read)
            if kill is None:
                passing.append((score(row), candidate))
        if passing:
            idx = weighted_sample([p[1] for p in passing],
                                  [p[0] for p in passing], rng, temp)
            return passing[idx][1], "retrieved", None
        # all corpus candidates killed; fall through to operator draft
        _pv = pool[0]["text"] if isinstance(pool[0]["text"], list) else [pool[0]["text"]]
        best_kill = critic(_fill(_pv[rng.randrange(len(_pv))] if (rng and len(_pv) > 1) else _pv[0], case),
                           case, read=read)

    if operator_draft:
        candidate = operator_draft.strip()
        kill = critic(candidate, case, read=read)
        if kill is None:
            return candidate, "edited", None
        return None, "none", kill

    return None, "none", best_kill if pool else None
