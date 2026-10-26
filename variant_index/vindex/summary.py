"""Per-axis evidence summaries. A summary lists what exists on each axis; it never ranks or adds axes together."""
from . import schema

AGGREGATE = 'aggregated_integration'     # summary rows: they restate other rows and never count as an independent source


def _lineage(r):
    return r['lineage_group'] or f"{r['source_id']}:{r['accession']}"


def summarise(rows):
    """Per-axis dict. Counts: n_rows; n_independent (distinct lineage groups, aggregates included); n_independent_by_assertion (distinct lineage groups per
    effect_assertion, aggregates excluded: replication of an effect is counted from the rows that assert it). Descriptors, each of one axis only:
    max_review_stars, min_p_value, max_pip, condition_ids_effect, condition_texts_effect. conflict is True when one axis holds both an effect and a no_effect assertion or a source-reported conflict."""
    out = {}
    for ax in schema.AXES:
        sel = [r for r in rows if r['axis'] == ax]
        if not sel:
            continue
        a, by = {}, {}
        for r in sel:
            a[r['effect_assertion']] = a.get(r['effect_assertion'], 0) + 1
            if r.get('evidence_kind') != AGGREGATE:
                by.setdefault(r['effect_assertion'], set()).add(_lineage(r))
        stars = [r['review_stars'] for r in sel if r.get('review_stars') is not None]
        pvals = [r['p_value'] for r in sel if r.get('p_value') is not None]
        pips = [r['pip'] for r in sel if r.get('pip') is not None]
        conds = sorted({c for r in sel if r['effect_assertion'] == 'effect' for c in (r.get('condition_ids') or ())})[:50]
        texts = sorted({r['condition_text'] for r in sel if r['effect_assertion'] == 'effect' and r.get('condition_text')})[:10]
        out[ax] = dict(n_rows=len(sel), n_independent=len({_lineage(r) for r in sel}), assertions=dict(sorted(a.items())),
                       n_independent_by_assertion={k: len(v) for k, v in sorted(by.items())},
                       sources=sorted({r['source_id'] for r in sel}), evidence_ids=sorted(r['evidence_id'] for r in sel),
                       conflict=('effect' in a and 'no_effect' in a) or 'conflicting' in a,
                       max_review_stars=max(stars) if stars else None, min_p_value=min(pvals) if pvals else None, max_pip=max(pips) if pips else None,
                       condition_ids_effect=conds, condition_texts_effect=texts)
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
