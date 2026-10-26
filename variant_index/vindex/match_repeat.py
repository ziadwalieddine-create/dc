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
    status: str                   # GENOTYPED | UNRESOLVED
    classes: tuple = ()
    counts: tuple = ()
    thresholds: dict = field(default_factory=dict)
    flags: tuple = ()
    missing: str = ''
    inheritance: str = ''
    disease: str = ''


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
            out.append(RepeatResult(L['locus_id'], L.get('gene', ''), 'UNRESOLVED', thresholds=th, missing=MISSING, inheritance=L.get('inheritance', ''), disease=L.get('disease', '')))
            continue
        flags = ()
        if not r.qc_pass:
            flags += ('repeat_call_fails_filter',)
        if not r.ci:
            flags += ('no_confidence_interval',)
        cls = tuple(classify_allele(c, r.ci[i] if i < len(r.ci) else None, th['normal_max'], th['pathologic_min']) for i, c in enumerate(r.counts))
        out.append(RepeatResult(L['locus_id'], L.get('gene', ''), 'GENOTYPED', cls, r.counts, th, flags, inheritance=L.get('inheritance', ''), disease=L.get('disease', '')))
    return out
