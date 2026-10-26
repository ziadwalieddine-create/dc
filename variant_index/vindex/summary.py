"""Per-axis evidence summaries. A summary lists what exists on each axis; it never ranks or adds axes together."""
from . import schema


def summarise(rows):
    """Per-axis dict: n_rows, n_independent (distinct lineage groups), assertions (counts), sources, evidence_ids, conflict flag."""
    out = {}
    for ax in schema.AXES:
        sel = [r for r in rows if r['axis'] == ax]
        if not sel:
            continue
        a = {}
        for r in sel:
            a[r['effect_assertion']] = a.get(r['effect_assertion'], 0) + 1
        groups = {r['lineage_group'] or f"{r['source_id']}:{r['accession']}" for r in sel}
        out[ax] = dict(n_rows=len(sel), n_independent=len(groups), assertions=dict(sorted(a.items())),
                       sources=sorted({r['source_id'] for r in sel}), evidence_ids=sorted(r['evidence_id'] for r in sel),
                       conflict=('effect' in a and 'no_effect' in a) or 'conflicting' in a)
    return out


def contradictions(rows):
    """Pairs of rows about the same subject, axis and condition that assert effect and no_effect; both stay in the index."""
    out = []
    by = {}
    for r in rows:
        by.setdefault((r['subject_id'], r['axis'], r['claim_type']), []).append(r)
    for k, rs in by.items():
        eff = [r for r in rs if r['effect_assertion'] == 'effect']
        non = [r for r in rs if r['effect_assertion'] == 'no_effect']
        for e in eff:
            for n in non:
                if set(e['condition_ids'] or []) & set(n['condition_ids'] or []) or not (e['condition_ids'] and n['condition_ids']):
                    out.append((e['evidence_id'], n['evidence_id'], k))
    return out
