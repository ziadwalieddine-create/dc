"""Known-locus checks: facts that are established in the literature and in public curation, tested against the built index.

The unit tests and the field-by-field comparisons prove that the index says what its sources say. These checks ask a different question: does the index, read the way a user
reads it, give the answer that the field already knows? An effect allele that was flipped, an rsID attached to the wrong site, a classification lost on the way, a star allele
without its defining position would all pass a comparison with the source file and fail here.

A spec is a list of loci; each locus has a public fact (`basis`) and checks:

  xref_key         the dbSNP identifier resolves to the stated normalised key, and to no other key at another position
  evidence         at least `min_rows` evidence rows of the subject (variant key, or the gene or region the check names) satisfy `where` (column == value or in a list), `regex`
                   (re.search on a column, case-insensitive) and the numeric bounds `bounds` {column: [low, high]}
  orientation      the effect-allele orientation of association rows: rows whose effect allele is `risk_allele` must raise the risk (odds ratio above 1, or beta above 0) and rows
                   whose effect allele is `other_allele` must lower it, in at least `min_fraction` of the rows that carry an effect size
  haplotype_allele the star allele (or named allele) exists for the gene, and its definition includes the stated variant key
  haplotype_unevaluable
                   the allele exists for the gene and the index keeps it as a definition that cannot be evaluated against a callset (`definition_scope` unresolved:<reason>, for
                   example a repeat-count notation), never silently dropped and never turned into a variant that does not define it
  region_evidence  a region of the stated kind overlaps the interval and carries evidence satisfying `where` / `regex`
  unmapped_evidence at least `min_rows` evidence rows that the index kept without a variant (mapping status and reason in `where`) under the author's identifier (`unmapped_key`) on
                   the chromosome's partition: records that cannot be placed are retained with their reason, never dropped
  assay_exists     at least `min_rows` assays (MaveDB score sets and other assays) target the gene (found by its symbol through the gene table: assays name their target by gene id)

A selection can also use `not_regex` (rows matching any of these patterns are excluded) and `list_any` ({column: values}: the list column holds at least one of the values).

A check is PASS, FAIL (the index contradicts the fact) or ABSENT (the index holds nothing about the subject: a gap, not a contradiction). Nothing here edits the index.
"""
import json
import os
import re

import pyarrow.compute as pc

DEFAULT_TEXT_COLUMNS = ('claim_value', 'condition_text', 'extra_json')


def _matches(row, where=None, regex=None, bounds=None, text_regex=None, not_regex=None, list_any=None, text_columns=DEFAULT_TEXT_COLUMNS):
    for col, want in (where or {}).items():
        have = row.get(col)
        if isinstance(want, list):
            if have not in want:
                return False
        elif have != want:
            return False
    for col, pat in (regex or {}).items():
        v = row.get(col)
        if v is None or not re.search(pat, str(v), re.I):
            return False
    for col, (lo, hi) in (bounds or {}).items():
        v = row.get(col)
        if v is None or (lo is not None and v < lo) or (hi is not None and v > hi):
            return False
    for col, pat in (not_regex or {}).items():
        v = row.get(col)
        if v is not None and re.search(pat, str(v), re.I):
            return False
    for col, values in (list_any or {}).items():
        if not set(row.get(col) or ()) & set(values):
            return False
    if text_regex is not None:
        joined = ' '.join(str(row.get(c)) for c in text_columns if row.get(c) is not None)
        if not re.search(text_regex, joined, re.I):
            return False
    return True


def select(rows, check):
    """The rows a check's selection keeps."""
    return [r for r in rows if _matches(r, check.get('where'), check.get('regex'), check.get('bounds'), check.get('text_regex'), check.get('not_regex'), check.get('list_any'))]


def _result(check, status, observed, **extra):
    return dict(id=check['id'], type=check['type'], status=status, observed=observed, **extra)


def _rsid_keys(index, rsid):
    cache = getattr(index, '_known_loci_xref', None)
    if cache is not None and rsid in cache:
        rows = cache[rsid]
    else:
        rows = index.filtered('xref', (pc.field('namespace') == 'dbsnp') & (pc.field('identifier') == rsid)).to_pylist()
    return rows, sorted({r['entity_id'] for r in rows if r['entity_level'] == 'variant'})


def check_xref_key(index, locus, check):
    rows, keys = _rsid_keys(index, locus['rsid'])
    want = locus['key']
    sources = sorted({r['relation'] for r in rows if r['entity_id'] == want})
    positions = sorted({tuple(k.split(':')[:2]) for k in keys})
    if not keys:
        return _result(check, 'ABSENT', f"{locus['rsid']} has no dbSNP cross-reference in the index")
    if want not in keys:
        return _result(check, 'FAIL', f"{locus['rsid']} resolves to {keys}, not {want}")
    if len(positions) > 1:
        return _result(check, 'FAIL', f"{locus['rsid']} resolves to {len(positions)} different positions: {keys}")
    return _result(check, 'PASS', f"{locus['rsid']} -> {want} (asserted by {', '.join(sources)})")


def _variant_rows(index, key):
    cache = getattr(index, '_known_loci_rows', None)
    if cache is not None and key in cache:
        return cache[key]
    return index.evidence_for_variants([key]).get(key, [])


def check_evidence(index, locus, check, rows=None):
    subject = check.get('subject', 'variant')
    if rows is None:
        rows = _variant_rows(index, locus['key']) if subject == 'variant' else _gene_rows(index, check['gene_symbol'])
    if not rows:
        return _result(check, 'ABSENT', 'the index holds no evidence rows for the subject')
    hits = select(rows, check)
    need = check.get('min_rows', 1)
    sources = sorted({r['source_id'] for r in hits})
    observed = f"{len(hits)} matching row(s) of {len(rows)} from {', '.join(sources) or 'no source'}"
    return _result(check, 'PASS' if len(hits) >= need else 'FAIL', observed, n_matching=len(hits), n_rows=len(rows))


def _gene_rows(index, symbol):
    genes = index.filtered('entity_gene', pc.field('symbol') == symbol).to_pylist()
    ids = [g['gene_id'] for g in genes]
    out = []
    for rows in index.evidence_for('gene', ids).values():
        out.extend(rows)
    return out


def raises_risk(row):
    """+1 when the row says the effect allele raises the risk, -1 when it lowers it, 0 when the row carries no usable effect size."""
    odds, beta = row.get('odds_ratio'), row.get('beta')
    if odds is not None and odds > 0 and odds != 1:
        return 1 if odds > 1 else -1
    if beta is not None and beta != 0:
        return 1 if beta > 0 else -1
    return 0


def check_orientation(index, locus, check):
    rows = select(_variant_rows(index, locus['key']), dict(check, bounds=None))
    risk, other = check['risk_allele'], check['other_allele']
    concordant = discordant = uninformative = 0
    examples = []
    for r in rows:
        sign = raises_risk(r)
        allele = r.get('effect_allele')
        if allele not in (risk, other) or sign == 0:
            uninformative += 1
            continue
        if (allele == risk) == (sign > 0):
            concordant += 1
        else:
            discordant += 1
            if len(examples) < 3:
                examples.append(dict(accession=r['accession'], effect_allele=allele, odds_ratio=r.get('odds_ratio'), beta=r.get('beta'), condition_text=r.get('condition_text')))
    informative = concordant + discordant
    if not rows:
        return _result(check, 'ABSENT', 'no association rows satisfy the selection')
    if informative < check.get('min_rows', 3):
        return _result(check, 'ABSENT', f'{len(rows)} selected row(s), only {informative} with an effect size and one of the two alleles', n_selected=len(rows))
    frac = concordant / informative
    return _result(check, 'PASS' if frac >= check.get('min_fraction', 0.9) else 'FAIL',
                   f'{concordant} of {informative} informative rows have the expected orientation for {risk} (raises) / {other} (lowers); {uninformative} uninformative',
                   concordant=concordant, discordant=discordant, uninformative=uninformative, discordant_examples=examples)


def check_haplotype_allele(index, locus, check):
    haps = [h for h in index.haplotypes(gene_symbol=check['gene_symbol']) if re.search(check['allele_regex'], h.get('allele_name') or '', re.I)]
    if not haps:
        return _result(check, 'ABSENT', f"no allele of {check['gene_symbol']} matches {check['allele_regex']!r}")
    defs = index.haplotype_definitions([h['haplotype_id'] for h in haps])
    with_key = sorted({h['allele_name'] for h in haps for d in defs.get(h['haplotype_id'], ()) if d['variant_key'] == locus['key']})
    names = sorted({h['allele_name'] for h in haps})
    if with_key:
        return _result(check, 'PASS', f"{locus['key']} defines {', '.join(with_key)} of {check['gene_symbol']}")
    return _result(check, 'FAIL', f"alleles {names[:6]} of {check['gene_symbol']} exist but none is defined by {locus['key']}")


def check_haplotype_unevaluable(index, locus, check):
    haps = [h for h in index.haplotypes(gene_symbol=check['gene_symbol']) if re.search(check['allele_regex'], h.get('allele_name') or '', re.I)]
    if not haps:
        return _result(check, 'ABSENT', f"no allele of {check['gene_symbol']} matches {check['allele_regex']!r}")
    defs = index.haplotype_definitions([h['haplotype_id'] for h in haps])
    want = check.get('reason')
    found = {}
    for h in haps:
        for d in defs.get(h['haplotype_id'], ()):
            scope = d.get('definition_scope') or ''
            if scope.startswith('unresolved:') and (not want or scope == 'unresolved:' + want):
                found.setdefault(h['allele_name'], set()).add(scope)
    if found:
        return _result(check, 'PASS', '; '.join(f"{name}: {', '.join(sorted(scopes))}" for name, scopes in sorted(found.items())))
    scopes = sorted({d.get('definition_scope') for h in haps for d in defs.get(h['haplotype_id'], ())})
    return _result(check, 'FAIL', f"alleles {sorted({h['allele_name'] for h in haps})[:6]} of {check['gene_symbol']} exist but no definition is kept as unresolved{' for ' + want if want else ''} (scopes: {scopes[:6]})")


def check_region_evidence(index, locus, check):
    regions = index.regions_overlapping(check['chrom'], check['start'], check['end'], check.get('kinds'))
    if not regions:
        return _result(check, 'ABSENT', f"no region of kind {check.get('kinds')} overlaps {check['chrom']}:{check['start']}-{check['end']}")
    ev = index.evidence_for('region', [r['region_id'] for r in regions])
    rows = [r for rs in ev.values() for r in rs]
    hits = select(rows, check)
    return _result(check, 'PASS' if len(hits) >= check.get('min_rows', 1) else 'FAIL', f'{len(regions)} overlapping region(s), {len(hits)} of {len(rows)} evidence rows match')


def check_unmapped_evidence(index, locus, check):
    from . import schema
    part = schema.partition_of(check['chrom'])
    rows = index.filtered('evidence', pc.field('unmapped_key') == check['unmapped_key'], sub=os.path.join('variant', part)).to_pylist()
    if not rows:
        return _result(check, 'ABSENT', f"no evidence row is kept under {check['unmapped_key']!r} on chromosome {check['chrom']}")
    hits = select(rows, check)
    states = describe_states(rows)
    return _result(check, 'PASS' if len(hits) >= check.get('min_rows', 1) else 'FAIL', f'{len(hits)} matching row(s) of {len(rows)} kept under the identifier ({states})')


def describe_states(rows):
    from collections import Counter
    c = Counter((r['mapping_status'], r.get('mapping_reason')) for r in rows)
    return '; '.join(f'{n} {st}{": " + why if why else ""}' for (st, why), n in sorted(c.items(), key=lambda kv: -kv[1]))


def check_assay_exists(index, locus, check):
    symbol = check['gene_symbol']
    ids = {g['gene_id'] for g in index.filtered('entity_gene', pc.field('symbol') == symbol).to_pylist()}               # assays name their target by gene id
    assays = [a for a in index.table('entity_assay').to_pylist() if a.get('target_gene') in ids or a.get('target_gene') == symbol]
    sources = sorted({a['source_id'] for a in assays})
    if not assays:
        return _result(check, 'ABSENT', f'no assay targets {symbol} ({len(ids)} gene id(s) carry that symbol)')
    return _result(check, 'PASS' if len(assays) >= check.get('min_rows', 1) else 'FAIL', f'{len(assays)} assay(s) from {", ".join(sources)}')


CHECKS = {'xref_key': check_xref_key, 'evidence': check_evidence, 'orientation': check_orientation, 'haplotype_allele': check_haplotype_allele,
          'haplotype_unevaluable': check_haplotype_unevaluable, 'region_evidence': check_region_evidence, 'assay_exists': check_assay_exists, 'unmapped_evidence': check_unmapped_evidence}


def validate_spec(spec):
    """Problems with a spec itself (a check that cannot run proves nothing)."""
    problems = []
    seen = set()
    for locus in spec.get('loci', ()):
        if not locus.get('id') or locus['id'] in seen:
            problems.append(f"locus id {locus.get('id')!r} missing or repeated")
        seen.add(locus.get('id'))
        for k in ('title', 'basis'):
            if not locus.get(k):
                problems.append(f"{locus.get('id')}: no {k}")
        for c in locus.get('checks', ()):
            if c.get('type') not in CHECKS:
                problems.append(f"{locus.get('id')}/{c.get('id')}: unknown check type {c.get('type')!r}")
            if c.get('type') in ('xref_key', 'haplotype_allele') and not locus.get('key'):
                problems.append(f"{locus.get('id')}/{c.get('id')}: needs a key")
            if c.get('type') == 'xref_key' and not locus.get('rsid'):
                problems.append(f"{locus.get('id')}/{c.get('id')}: needs an rsid")
            if c.get('type') == 'evidence' and c.get('subject', 'variant') == 'variant' and not locus.get('key'):
                problems.append(f"{locus.get('id')}/{c.get('id')}: a variant check needs a key")
            if c.get('type') == 'evidence' and c.get('subject') == 'gene' and not c.get('gene_symbol'):
                problems.append(f"{locus.get('id')}/{c.get('id')}: a gene check needs a gene_symbol")
            if c.get('type') == 'region_evidence' and not all(k in c for k in ('chrom', 'start', 'end')):
                problems.append(f"{locus.get('id')}/{c.get('id')}: a region check needs chrom, start and end")
            if c.get('type') in ('haplotype_allele', 'haplotype_unevaluable') and not (c.get('gene_symbol') and c.get('allele_regex')):
                problems.append(f"{locus.get('id')}/{c.get('id')}: a haplotype check needs gene_symbol and allele_regex")
            if c.get('type') == 'unmapped_evidence' and not (c.get('chrom') and c.get('unmapped_key')):
                problems.append(f"{locus.get('id')}/{c.get('id')}: an unmapped-evidence check needs chrom and unmapped_key")
            if c.get('type') == 'assay_exists' and not c.get('gene_symbol'):
                problems.append(f"{locus.get('id')}/{c.get('id')}: an assay check needs a gene_symbol")
            if not c.get('id'):
                problems.append(f"{locus.get('id')}: a check without an id")
            if c.get('type') == 'orientation' and not (c.get('risk_allele') and c.get('other_allele') and locus.get('key')):
                problems.append(f"{locus.get('id')}/{c.get('id')}: orientation needs a key, risk_allele and other_allele")
        if not locus.get('checks'):
            problems.append(f"{locus.get('id')}: no checks")
    return problems


def _prefetch(index, spec):
    """One read of the xref table and one of each chromosome partition for all the loci of the spec, instead of one scan per check (the index is tens of millions of rows)."""
    rsids = sorted({l['rsid'] for l in spec['loci'] if l.get('rsid')})
    keys = sorted({l['key'] for l in spec['loci'] if l.get('key')})
    index._known_loci_xref = {rs: [] for rs in rsids}
    if rsids:
        for r in index.filtered('xref', (pc.field('namespace') == 'dbsnp') & pc.field('identifier').isin(rsids)).to_pylist():
            index._known_loci_xref[r['identifier']].append(r)
    found = index.evidence_for_variants(keys) if keys else {}
    index._known_loci_rows = {k: found.get(k, []) for k in keys}


def run(index, spec):
    problems = validate_spec(spec)
    if problems:
        raise ValueError('invalid known-locus spec: ' + '; '.join(problems))
    loci, totals = [], dict(PASS=0, FAIL=0, ABSENT=0)
    _prefetch(index, spec)
    try:
        for locus in spec['loci']:
            results = []
            for c in locus['checks']:
                r = CHECKS[c['type']](index, locus, c)
                totals[r['status']] += 1
                results.append(r)
            loci.append(dict(id=locus['id'], title=locus['title'], basis=locus['basis'], key=locus.get('key'), rsid=locus.get('rsid'), checks=results))
    finally:
        index._known_loci_xref = index._known_loci_rows = None
    return dict(ok=totals['FAIL'] == 0 and totals['ABSENT'] == 0, totals=totals, n_loci=len(loci), loci=loci)


def load_spec(path):
    with open(path, encoding='utf-8') as f:
        return json.load(f)
