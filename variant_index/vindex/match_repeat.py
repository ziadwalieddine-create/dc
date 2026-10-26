"""Repeat-expansion loci: a repeat length needs a repeat genotyper.

SNV/indel calls cannot size a repeat longer than the reads. Each indexed locus is therefore UNRESOLVED unless the input holds a repeat
genotype (ExpansionHunter REPCN or TRGT AL) at that locus; then each allele is placed against the indexed thresholds (normal maximum,
pathologic minimum) and borderline calls (confidence interval crossing a threshold) are labelled as such. Thresholds come from the
index (STRchive) and are per locus; a locus without thresholds is genotyped but not classified.
"""
from dataclasses import dataclass, field

MISSING = ('a repeat-length genotype at this locus from a dedicated caller (ExpansionHunter, TRGT or equivalent) with its confidence interval; '
           'short-read SNV/indel calls cannot establish repeat length')


@dataclass
class RepeatResult:
    locus_id: str
    gene: str
    status: str                   # GENOTYPED | UNRESOLVED (indexed locus) | NO_INDEX_RECORD (a repeat genotype at a locus the index does not know)
    classes: tuple = ()
    counts: tuple = ()
    thresholds: dict = field(default_factory=dict)
    flags: tuple = ()
    missing: str = ''
    inheritance: str = ''
    disease: str = ''
    chrom: str = ''
    start: int = None
    end: int = None
    motif: str = ''
    caller: str = ''


def has_ci(call):
    """True when the caller stated a confidence interval for at least one allele (a placeholder '.' is not one)."""
    return any(c and c[0] is not None and c[1] is not None for c in call.ci or ())


def classify_allele(count, ci, normal_max, path_min):
    if count is None:
        return 'unknown'
    lo, hi = (ci if ci and ci[0] is not None else (count, count))
    if path_min is not None and lo >= path_min:
        return 'pathogenic_range'
    if normal_max is not None and hi <= normal_max:
        return 'normal_range'
    if path_min is not None and hi >= path_min > lo:
        return 'borderline_crossing_pathologic_minimum'
    if normal_max is not None and lo <= normal_max < hi:
        return 'borderline_crossing_normal_maximum'
    return 'intermediate'


def match_repeats(callset, loci):
    """loci: list of dicts(locus_id, gene, chrom, start, end, normal_max, pathologic_min, inheritance, disease)."""
    by_id = {r.locus_id: r for r in callset.repeats}
    out = []
    for L in loci:
        r = by_id.get(L['locus_id'])
        if r is None:
            r = next((x for x in callset.repeats if x.chrom == L['chrom'] and x.start <= L['end'] and x.end >= L['start']), None)
        th = dict(normal_max=L.get('normal_max'), pathologic_min=L.get('pathologic_min'))
        if r is None:
            out.append(RepeatResult(L['locus_id'], L.get('gene', ''), 'UNRESOLVED', thresholds=th, missing=MISSING, inheritance=L.get('inheritance', ''), disease=L.get('disease', ''),
                                    chrom=L['chrom'], start=L['start'], end=L['end']))
            continue
        flags = ()
        if not r.qc_pass:
            flags += ('repeat_call_fails_filter',)
        if not has_ci(r):
            flags += ('no_confidence_interval',)
        cls = tuple(classify_allele(c, r.ci[i] if i < len(r.ci) else None, th['normal_max'], th['pathologic_min']) for i, c in enumerate(r.counts))
        out.append(RepeatResult(L['locus_id'], L.get('gene', ''), 'GENOTYPED', cls, r.counts, th, flags, inheritance=L.get('inheritance', ''), disease=L.get('disease', ''),
                                chrom=r.chrom, start=r.start, end=r.end, motif=r.motif, caller=r.source))
    return out


def unindexed_repeats(callset, loci):
    """Repeat genotypes of the callset that match no indexed locus (by id or by overlap). They are kept, with the caller's counts and flags and no classification:
    the index holds no threshold for them, and a repeat call is never dropped for lacking one."""
    ids = {L['locus_id'] for L in loci}
    by_chrom = {}
    for L in loci:
        by_chrom.setdefault(L['chrom'], []).append(L)
    out = []
    for r in callset.repeats:
        if r.locus_id in ids or any(L['start'] <= r.end and L['end'] >= r.start for L in by_chrom.get(r.chrom, ())):
            continue
        flags = (() if r.qc_pass else ('repeat_call_fails_filter',)) + (() if has_ci(r) else ('no_confidence_interval',))
        out.append(RepeatResult(r.locus_id, '', 'NO_INDEX_RECORD', (), r.counts, {}, flags, chrom=r.chrom, start=r.start, end=r.end, motif=r.motif, caller=r.source))
    return out


def loci_from_index(index):
    """Locus dictionaries for match_repeats, rebuilt from the STRchive rows of an index: the region, the locus -> disease row and the two threshold rows.
    A threshold the source does not give stays None (such a locus is genotyped, not classified)."""
    regions = index.regions_by_kind(['strchive_locus'])
    ev = index.evidence_for('region', [r['region_id'] for r in regions])
    symbols = {g['gene_id']: g['symbol'] for g in index.table('entity_gene').to_pylist()} if index.table('entity_gene').num_rows else {}
    out = []
    for r in sorted(regions, key=lambda r: (r['chrom'], r['start'], r['region_id'])):
        by = {x['claim_type']: x for x in ev.get(r['region_id'], ())}
        gid = (r.get('gene_ids') or [None])[0]
        locus = by.get('strchive_disease_locus') or {}
        th = lambda k: int(by[k]['claim_num']) if k in by and by[k].get('claim_num') is not None else None
        out.append(dict(locus_id=r['name'], gene=symbols.get(gid) or (gid or '').replace('symbol:', ''), chrom=r['chrom'], start=r['start'], end=r['end'], normal_max=th('repeat_normal_max'),
                        pathologic_min=th('repeat_pathologic_min'), inheritance=locus.get('inheritance') or '', disease=locus.get('condition_text') or ''))
    return out
