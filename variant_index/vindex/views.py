"""Rules and views over the complete catalogue of matches.

High-value rules say which indexed variants are worth asking "what does this sample say about this allele?" even when the sample's file holds no record
of it (the unresolved-loci report). They are filters over evidence rows, each stated in plain words, and a variant may satisfy several. They are not
scores and they do not rank.

A view is a filter and an explicit ordering over the matches of one sample. It never changes or replaces the per-axis evidence, no number here combines
axes, and a call that fails a filter stays in the full call table. The ordering of each view is part of its definition and is printed with it.
"""
import json
import os
from collections import Counter
from dataclasses import dataclass

import pyarrow as pa
import pyarrow.compute as pc
import pyarrow.parquet as pq

from . import schema

F = pc.field


@dataclass(frozen=True)
class Rule:
    name: str
    description: str
    expr: object          # Arrow expression over evidence rows (variant level, mapped)


HIGH_VALUE_RULES = (
    Rule('clinvar_pathogenic_reviewed', 'ClinVar aggregate germline classification asserting an effect (Pathogenic / Likely pathogenic) with at least two review stars '
         '(multiple submitters without conflict, expert panel or practice guideline)',
         (F('source_id') == 'clinvar') & (F('claim_type') == 'germline_classification') & (F('effect_assertion') == 'effect') & (F('review_stars') >= 2)),
    Rule('clingen_vcep_classification', 'ClinGen expert-panel (VCEP) classification of the variant asserting an effect',
         (F('source_id') == 'clingen') & (F('claim_type') == 'vcep_classification') & (F('effect_assertion') == 'effect')),
    Rule('clinpgx_level_1', 'ClinPGx clinical annotation at level of evidence 1A or 1B', (F('source_id') == 'pharmgkb') & (F('claim_type') == 'clinpgx_level_of_evidence') & F('claim_value').isin(['1A', '1B'])),
    Rule('functional_assay_abnormal', 'MaveDB calibrated class "abnormal" for the specific allele',
         (F('source_id') == 'mavedb') & (F('claim_type') == 'mavedb_calibrated_class') & (F('effect_assertion') == 'effect')),
    Rule('gnomad_lof_curated', 'gnomAD curators call the predicted loss-of-function variant LoF or likely LoF',
         (F('source_id') == 'gnomad') & (F('claim_type') == 'gnomad_lof_curation_verdict') & (F('effect_assertion') == 'effect')),
    Rule('fine_mapped_pip_ge_0_5', 'GWAS credible-set variant with a posterior inclusion probability of at least 0.5',
         (F('source_id') == 'opentargets') & (F('claim_type') == 'gwas_fine_mapping_pip') & (F('pip') >= 0.5)),
)


def high_value_keys(index, chrom, rules=HIGH_VALUE_RULES):
    """{variant_key: {'rules': [names], 'evidence_ids': [...]}} for the mapped variant-level evidence of one contig (or, for schema.OTHER_PARTITION, of all contigs outside the primary assembly) that satisfies a rule."""
    out = {}
    for rule in rules:
        expr = rule.expr & (F('level') == 'variant') & (F('mapping_status') == 'mapped')
        if chrom != schema.OTHER_PARTITION:
            expr = expr & (F('chrom') == chrom)                     # one contig; schema.OTHER_PARTITION asks for every alternate locus and scaffold at once
        sub = os.path.join('variant', schema.partition_of(chrom))
        t = index.filtered('evidence', expr, sub=sub, columns=['subject_id', 'evidence_id'])
        for k, e in zip(t['subject_id'].to_pylist(), t['evidence_id'].to_pylist()):
            d = out.setdefault(k, dict(rules=[], evidence_ids=[]))
            if rule.name not in d['rules']:
                d['rules'].append(rule.name)
            d['evidence_ids'].append(e)
    for d in out.values():
        d['rules'].sort()
        d['evidence_ids'].sort()
    return out


# ---- views over one sample's matches -------------------------------------------------------------------------------------------------------------
CLIN, FUNC, ASSOC, FINE, MAG, CTX, MECH = ('clinical_validity_actionability', 'empirical_functional_support', 'association_replication', 'fine_mapping_probability',
                                           'effect_magnitude_direction', 'applicability_context', 'mechanistic_confidence')


def axis(row, name):
    return (row['axes'] or {}).get(name) or {}


def asserts(row, name, what='effect'):
    return (axis(row, name).get('assertions') or {}).get(what, 0) > 0


def independent(row, name, what='effect'):
    return (axis(row, name).get('n_independent_by_assertion') or {}).get(what, 0)


@dataclass(frozen=True)
class View:
    name: str
    title: str
    where: object
    order: object
    ordering: str          # the ordering in words: printed with the view


VIEWS = (
    View('clinical_effect', 'Variant-level clinical classification asserts an effect', lambda r: asserts(r, CLIN),
         lambda r: (-(axis(r, CLIN).get('max_review_stars') or 0), -independent(r, CLIN), r['chrom'], r['pos'], r['key']),
         'review stars (highest first), then independent lineages asserting the effect, then position'),
    View('functional_effect', 'Assay evidence for the specific allele asserts an effect', lambda r: asserts(r, FUNC),
         lambda r: (-independent(r, FUNC), r['chrom'], r['pos'], r['key']), 'independent lineages asserting the effect, then position'),
    View('association', 'Statistical association asserts an effect', lambda r: asserts(r, ASSOC),
         lambda r: (-independent(r, ASSOC), axis(r, ASSOC).get('min_p_value') if axis(r, ASSOC).get('min_p_value') is not None else 1.0, r['chrom'], r['pos'], r['key']),
         'independent studies asserting the association (most first), then smallest p-value'),
    View('fine_mapped', 'Fine-mapping places posterior probability on the allele', lambda r: (axis(r, FINE).get('max_pip') or 0) >= 0.1,
         lambda r: (-(axis(r, FINE).get('max_pip') or 0), r['chrom'], r['pos'], r['key']), 'highest posterior inclusion probability first'),
    View('pharmacogenomic', 'Pharmacogenomic annotation of the allele', lambda r: bool({'pharmgkb', 'cpic'} & set(axis(r, CLIN).get('sources') or ())),
         lambda r: (-(axis(r, CLIN).get('max_review_stars') or 0), r['chrom'], r['pos'], r['key']), 'position'),
    View('conflicting_evidence', 'An axis holds both an effect and a no-effect assertion, or the source reports a conflict', lambda r: any(a.get('conflict') for a in (r['axes'] or {}).values()),
         lambda r: (r['chrom'], r['pos'], r['key']), 'position'),
    View('several_axes_assert_effect', 'Effect asserted on two or more separate axes (the axes are listed; their number is not a score)',
         lambda r: len([a for a, d in (r['axes'] or {}).items() if (d.get('assertions') or {}).get('effect')]) >= 2,
         lambda r: (r['chrom'], r['pos'], r['key']), 'position'),
)


def read_calls(results_dir, statuses=('EXACT_MATCH',)):
    """Rows of small_calls.parquet with the given match statuses, with the axis summary parsed."""
    t = pq.read_table(os.path.join(results_dir, 'small_calls.parquet'), filters=[('status', 'in', list(statuses))]) if statuses else pq.read_table(os.path.join(results_dir, 'small_calls.parquet'))
    for r in t.to_pylist():
        r['axes'] = json.loads(r['axes_json']) if r.get('axes_json') else {}
        yield r


def build_views(rows, views=VIEWS):
    """{view name: dict(title, ordering, rows)} over a list of call rows (as from read_calls). Rows are returned whole; nothing is dropped from `rows`."""
    out = {}
    for v in views:
        sel = sorted((r for r in rows if v.where(r)), key=v.order)
        out[v.name] = dict(title=v.title, ordering=v.ordering, rows=sel)
    return out


def trait_view(rows):
    """Group matches by the condition that effect-asserting evidence names: {condition id: dict(texts, axes, calls)}. Ids are the sources' own (EFO, MONDO, MedGen,
    OMIM, HP...): no ontology roll-up is applied, so two spellings of one disease stay two entries."""
    out = {}
    for r in rows:
        for ax, d in (r['axes'] or {}).items():
            for c in d.get('condition_ids_effect') or ():
                e = out.setdefault(c, dict(texts=set(), axes=set(), calls=set()))
                e['texts'].update(d.get('condition_texts_effect') or ())
                e['axes'].add(ax)
                e['calls'].add(r['key'])
    return {c: dict(texts=sorted(e['texts'])[:5], axes=sorted(e['axes']), n_calls=len(e['calls']), calls=sorted(e['calls'])) for c, e in sorted(out.items())}


# The prefix of a source's trait id (upper-case) -> the prefixes under which the ontology table lists the same concept as a cross-reference: a MedGen concept id is a UMLS CUI.
XREF_PREFIXES = {'MEDGEN': ('UMLS',), 'MESH': ('MESH',), 'OMIM': ('OMIM', 'OMIMPS'), 'ORPHANET': ('ORPHANET',), 'HP': ('HP',), 'NCIT': ('NCIT',), 'DOID': ('DOID',)}


def ontology_terms(path=None):
    """{term id: dict(name, ancestors, areas, xrefs, replaces)} from the Open Targets disease/phenotype ontology table (EFO, MONDO, HP, Orphanet, OBA, GO terms with their ancestors and
    therapeutic areas), keyed as the sources write the id (EFO:0004339). `xrefs` are the cross-references the table lists for the term (UMLS:C0011860, OMIM:125853, upper-cased) and
    `replaces` the retired terms it replaces (EFO:0001360). None when the table is not acquired."""
    import glob
    from .build.common import DATA
    files = sorted(glob.glob(os.path.join(path or os.path.join(DATA, 'opentargets', '26.09', 'disease'), '*.parquet')))
    if not files:
        return None
    out = {}
    for f in files:
        columns = [c for c in ('id', 'name', 'ancestors', 'therapeuticAreas', 'dbXRefs', 'obsoleteTerms') if c in pq.ParquetFile(f).schema_arrow.names]
        for r in pq.read_table(f, columns=columns).to_pylist():
            out[r['id'].replace('_', ':', 1)] = dict(name=r['name'], ancestors=[a.replace('_', ':', 1) for a in r['ancestors'] or ()], areas=[a.replace('_', ':', 1) for a in r['therapeuticAreas'] or ()],
                                                     xrefs=sorted({x.upper() for x in r.get('dbXRefs') or ()}), replaces=[x.replace('_', ':', 1) for x in r.get('obsoleteTerms') or ()])
    return out


def trait_crosswalk(ontology):
    """(by cross-reference, by retired term): {PREFIX:local id: {term ids}} and {retired term id: {term ids that replace it}}, built from the xrefs and replaces of ontology_terms."""
    by_xref, by_retired = {}, {}
    for tid, t in ontology.items():
        for x in t.get('xrefs') or ():
            by_xref.setdefault(x, set()).add(tid)
        for o in t.get('replaces') or ():
            by_retired.setdefault(o, set()).add(tid)
    return by_xref, by_retired


def cross_reference(cid, by_xref, by_retired):
    """(terms, how) for a trait id the ontology does not hold itself: the terms whose cross-references list it under the prefixes it can carry (XREF_PREFIXES; the prefix is part of the
    match, so OMIM:100100 never matches Orphanet:100100), and the terms that replace it when it is a retired term. Several terms is a many-to-one cross-reference, which the caller reports."""
    prefix, _, local = cid.partition(':')
    terms, how = set(), set()
    if cid in by_retired:
        terms |= by_retired[cid]
        how.add('replaces_a_retired_term')
    for p in XREF_PREFIXES.get(prefix.upper(), (prefix.upper(),)):
        hit = by_xref.get(f'{p}:{local}'.upper())
        if hit:
            terms |= hit
            how.add('dbXRef')
    return sorted(terms), sorted(how)


def trait_rollup(traits, ontology):
    """Roll the traits of trait_view up through an ontology: by therapeutic area, and the ancestors of each trait. A trait the ontology holds is rolled up through its own term. A trait it
    does not hold (MedGen, OMIM, MeSH, retired EFO / MONDO / HP / Orphanet terms) is rolled up through the terms that list it as a cross-reference or replace it, when there are any,
    and the mapping is kept under `cross_referenced` with how and whether it is many-to-one: a cross-reference groups, it is not an assertion that the concepts are the same. A trait that
    maps to nothing keeps its own entry and is listed as not rolled up; a trait in several areas is counted in each. The call lists are those of trait_view, unchanged."""
    by_xref, by_retired = trait_crosswalk(ontology)
    areas, unmapped, crossed = {}, [], {}
    for cid, e in traits.items():
        t = ontology.get(cid)
        if t is not None:
            terms = [cid]
        else:
            terms, how = cross_reference(cid, by_xref, by_retired)
            if not terms:
                unmapped.append(cid)
                continue
            crossed[cid] = dict(terms=terms, how=how, many_to_one=len(terms) > 1)
        for a in {a for x in terms for a in (ontology[x]['areas'] or [x])}:
            g = areas.setdefault(a, dict(name=(ontology.get(a) or {}).get('name') or a, traits=set(), crossed=set(), calls=set()))
            g['traits'].add(cid)
            if cid in crossed:
                g['crossed'].add(cid)
            g['calls'].update(e['calls'])
    return dict(by_therapeutic_area={a: dict(name=g['name'], n_traits=len(g['traits']), n_traits_cross_referenced=len(g['crossed']), n_calls=len(g['calls']), calls=sorted(g['calls']))
                                     for a, g in sorted(areas.items())},
                ancestors={cid: sorted({a for x in (crossed[cid]['terms'] if cid in crossed else [cid]) for a in ontology[x]['ancestors']}) for cid in sorted(traits) if cid in ontology or cid in crossed},
                cross_referenced=dict(sorted(crossed.items())), not_rolled_up=sorted(unmapped))


def pathway_view(rows, members):
    """Group matches by pathway through the genes they sit in. members: iterable of gene_set_member rows. A gene's pathway membership is an annotation of the
    gene; the matched variant's own evidence is unchanged and the pathway carries no evidence of its own."""
    by_gene = {}
    for m in members:
        by_gene.setdefault(m['gene_id'], []).append(m)
    out = {}
    for r in rows:
        for g in r.get('gene_ids') or ():
            for m in by_gene.get(g, ()):
                e = out.setdefault(m['set_id'], dict(name=m['name'], genes=set(), calls=set()))
                e['genes'].add(g)
                e['calls'].add(r['key'])
    return {p: dict(name=e['name'], n_genes=len(e['genes']), genes=sorted(e['genes']), n_calls=len(e['calls']), calls=sorted(e['calls'])) for p, e in sorted(out.items())}


VIEW_COLUMNS = ('key', 'chrom', 'pos', 'ref', 'alt', 'vtype', 'gt', 'zygosity', 'phased', 'status', 'flags', 'qc_flags', 'gene_symbols', 'gene_ids', 'region_kinds', 'n_evidence', 'axes_json', 'evidence_ids')


AXIS_SHORT = {'clinical_validity_actionability': 'clinical', 'empirical_functional_support': 'functional', 'mechanistic_confidence': 'mechanism', 'effect_magnitude_direction': 'effect size',
              'association_replication': 'association', 'fine_mapping_probability': 'fine-mapping', 'applicability_context': 'context'}


def describe(row):
    """One line per call: what each axis holds, in counts and descriptors of that axis alone (nothing is added up across axes)."""
    parts = []
    for ax, d in sorted((row['axes'] or {}).items()):
        bits = ', '.join(f'{n} {k}' for k, n in (d.get('assertions') or {}).items())
        extra = []
        if d.get('max_review_stars') is not None and (ax == CLIN or d['max_review_stars'] > 0):
            extra.append(f"review stars up to {d['max_review_stars']}")
        if d.get('min_p_value') is not None:
            extra.append(f"min p {d['min_p_value']:.1e}")
        if d.get('max_pip') is not None:
            extra.append(f"max PIP {d['max_pip']:.2f}")
        if d.get('conflict'):
            extra.append('conflicting assertions')
        parts.append(f"{AXIS_SHORT.get(ax, ax)}: {bits} ({d.get('n_independent', 0)} lineage{'s' if d.get('n_independent', 0) != 1 else ''})" + (f" [{'; '.join(extra)}]" if extra else ''))
    return ' | '.join(parts) or 'no variant-level evidence'


def structural_section(results_dir, load, top=25):
    """Report lines for the structural and copy-number calls: how many were called, how many have an equivalent event in the index, how many encompass a gene, and the calls whose
    equivalent events carry evidence (listed with the evidence's source, claim type and assertion). Frequencies of the events (gnomAD-SV, DGV, HGSVC3) are context, not evidence."""
    p = os.path.join(results_dir, 'sv_calls.parquet')
    if not os.path.exists(p):
        return []
    sv = pq.read_table(p, columns=['call_id', 'chrom', 'start', 'end', 'svtype', 'gt', 'dosage', 'n_equivalent', 'equivalent_event_ids', 'context_json', 'event_evidence_ids']).to_pylist()
    if not sv:
        return []
    by_type = Counter(r['svtype'] for r in sv)
    genes = sum(1 for r in sv if any(c.get('relation') == 'encompasses_gene' for c in json.loads(r['context_json'] or '[]')))
    with_event = sum(1 for r in sv if r['n_equivalent'])
    with_evidence = [r for r in sv if r['event_evidence_ids']]
    L = ['## Structural and copy-number calls', '', f"{len(sv):,} calls ({', '.join(f'{k} {v:,}' for k, v in sorted(by_type.items()))}); {with_event:,} have at least one equivalent event in the index (same type, breakpoints "
         f"within the stated tolerance, dosage consistent) and {genes:,} encompass at least one whole gene. Overlap alone is never equivalence: refused near-misses are in sv_calls.parquet with their metrics. "
         f"{len(with_evidence):,} call(s) have an equivalent event that carries evidence of its own (the frequency of a gnomAD-SV, DGV or HGSVC3 event is context, not evidence).", '']
    if with_evidence:
        links = {}
        for r in load('evidence_links.parquet'):
            if r.get('call_kind') == 'sv':
                links.setdefault(r['call_id'], []).append(r)
        L += [f'The first {min(top, len(with_evidence))} by position:', '', '| Call | Type | Genotype | Equivalent events | Evidence of those events |', '|---|---|---|---:|---|']
        for r in sorted(with_evidence, key=lambda r: (r['chrom'], r['start'], r['call_id']))[:top]:
            ev = Counter(f"{x['source_id']}: {x['claim_type']} ({x['effect_assertion']})" for x in links.get(r['call_id'], ()))
            L.append(f"| `{r['call_id']}` | {r['svtype']} | {r['gt']} | {r['n_equivalent']} | " + '; '.join(f'{k} x{n}' if n > 1 else k for k, n in sorted(ev.items())).replace('|', '/') + ' |')
        L.append('')
    return L


def reference_allele_section(results_dir, counts):
    """Report lines for the evidence whose effect allele is the reference allele: the indexed sites it concerns by what the sample says about them, and the links made at the
    sites whose genotype is a verified reference. A variant-only file without a callable-region mask or reference blocks verifies none: every such site is then absent from the callset."""
    sites = counts.get('reference_allele_effect_sites')
    if not sites:
        return []
    t = pq.read_table(os.path.join(results_dir, 'evidence_links.parquet'), columns=['call_kind', 'call_id', 'source_id', 'effect_assertion'], filters=[('call_kind', '=', 'verified_reference')])
    by = Counter(zip(t['source_id'].to_pylist(), t['effect_assertion'].to_pylist()))
    L = ['## Evidence about the reference allele', '',
         f"{sum(sites.values()):,} indexed single-nucleotide sites have evidence whose effect allele is the reference allele. A sample reaches such a row through an alternate call (its links state the "
         "copies of the effect allele) or through a verified reference genotype (a gVCF reference block, or a callable-region mask with its stated assumption). What the sample says about those sites: "
         + ', '.join(f'{k} {v:,}' for k, v in sorted(sites.items(), key=lambda kv: (-kv[1], kv[0]))) + '.', '',
         f"{len(set(t['call_id'].to_pylist())):,} verified-reference sites are linked to {t.num_rows:,} rows (call_kind `verified_reference` in evidence_links.parquet; two copies of the effect allele on an "
         "autosome, none stated elsewhere). A site the sample cannot answer for (absent from a variant-only file, no-call, filtered, inside a structural call, another allele called there) has no link; "
         "it is listed as unresolved only when a high-value rule selects it.", '']
    if by:
        L += ['| Source | Assertion | Links |', '|---|---|---:|'] + [f'| {s} | {a} | {n:,} |' for (s, a), n in sorted(by.items(), key=lambda kv: (-kv[1], kv[0]))]
        L.append('')
    return L


def read_evidence_section(un, counts):
    """Report lines for what aligned reads say about the unresolved loci: by the matcher's status and the state of the reads, and what they settle. Empty for a run without reads. The
    matcher's status is never changed by the reads (vindex/reads.py); the states mean fragments that carry the indexed allele, not a genotype."""
    rc = (counts or {}).get('reads')
    if not rc:
        return []
    by = Counter((r['status'], r['read_state']) for r in un if r.get('read_state'))
    states = sorted({st for _, st in by})
    L = ['### Read evidence for the unresolved loci', '',
         f"Aligned reads were asked about {rc['loci_evaluated']:,} of them: whether any fragment carries the indexed allele, not what the genotype is (fragment counts are not a dosage; no zygosity, "
         "copy number or phase is inferred). The status above is the callset's and stands; the reads answer the question it leaves open for "
         f"{rc['settles_the_open_question']['yes']:,} loci (absent from the callset, no-call, filtered, or a call that fails quality) and are context only for the other "
         f"{rc['settles_the_open_question']['no']:,} (inside a structural call, among others).", '']
    if by:
        L += ['| Callset status | ' + ' | '.join(st.replace('READS_', '') for st in states) + ' |', '|---|' + '---:|' * len(states)]
        for status in sorted({s for s, _ in by}):
            L.append(f'| {status} | ' + ' | '.join(f'{by.get((status, st), 0):,}' for st in states) + ' |')
        L.append('')
    meaning = Counter(r['read_interpretation'] for r in un if r.get('read_interpretation'))
    L += ['What they say: ' + ', '.join(f'{k.replace("_", " ")} {v:,}' for k, v in sorted(meaning.items(), key=lambda kv: (-kv[1], kv[0]))) + '.', '']
    if rc.get('not_evaluated_reasons'):
        L += ['Not asked of reads: ' + ', '.join(f'{k.replace("_", " ")} {v:,}' for k, v in rc['not_evaluated_reasons'].items()) + '.', '']
    return L


def render_report(results_dir, top=25):
    """report.md: what the run found, in the order the tables are produced. Views list the first `top` rows of their ordering; the full views are in views/*.parquet and
    every call is in small_calls.parquet. Written next to the results it summarises, so a sample's report stays where the sample's data are."""
    def load(name):
        p = os.path.join(results_dir, name)
        return pq.read_table(p).to_pylist() if os.path.exists(p) else []
    with open(os.path.join(results_dir, 'run_manifest.json')) as f:
        run = json.load(f)
    with open(os.path.join(results_dir, 'views', 'views_manifest.json')) as f:
        vm = json.load(f)
    c = run['counts']
    L = ['# Sample report', '',
         f"Run finished {run['finished_utc']} in {run['seconds']} s; sample `{run['input']['sample']}`, input `{os.path.basename(run['input']['path'])}` (sha256 {str(run['input'].get('sha256'))[:16] or 'not computed'}), "
         f"index schema {run['index'].get('schema_version')}, code {str((run.get('code') or {}).get('commit'))[:8]}{' plus uncommitted changes in variant_index/' if (run.get('code') or {}).get('uncommitted_changes_in_variant_index') else ''}{' (files under variant_index/ changed while the run was in progress)' if run.get('code_changed_during_run') else ''}. Assembly check: {json.dumps(run['input']['assembly'].get('ok'))}.", '',
         '## What was called and what the index said', '',
         f"{c['calls']:,} small-variant alleles with dosage above zero. Match status: " + ', '.join(f'{k} {v:,}' for k, v in sorted(c['by_status'].items(), key=lambda kv: -kv[1])) + '.',
         '`NO_INDEX_RECORD` means the index holds no record of the allele: that is absence from the index, not evidence of no effect. `FILTERED_CALL` is listed but its presence is not asserted. '
         f"{c['difficult_region_calls']:,} calls lie in GIAB difficult regions (flagged, not removed).", '',
         f"Structural calls: {sum(c['structural'].values()):,} ({', '.join(f'{k} {v}' for k, v in sorted(c['structural'].items())) or 'none'}); " +
         f"repeat genotypes: {', '.join(f'{k} {v}' for k, v in sorted(c['repeat_status'].items())) or 'none'}; star-allele genes: {', '.join(f'{k} {v}' for k, v in sorted(c['haplotype_status'].items())) or 'none'}.", '']
    rej = load('rejected_calls.parquet')
    if rej:
        L += ['Alleles the reader refused: ' + ', '.join(f'{k} {v}' for k, v in sorted(Counter(r['reason'] for r in rej).items())) + ' (listed in rejected_calls.parquet).', '']
    L += structural_section(results_dir, load, top)
    L += ['## High-value indexed alleles and what the sample says about them', '', '| Rule | ' + ' | '.join(sorted({k for v in c['high_value_rules'].values() for k in v})) + ' |', '|---|' + '---|' * len({k for v in c['high_value_rules'].values() for k in v})]
    cols = sorted({k for v in c['high_value_rules'].values() for k in v})
    for name, v in c['high_value_rules'].items():
        L.append(f"| `{name}` | " + ' | '.join(str(v.get(k, 0)) for k in cols) + ' |')
    L += ['', 'The rules are listed in run_manifest.json with their wording; they select rows, they do not rank them.', '']
    L += reference_allele_section(results_dir, c)
    un = load('unresolved_loci.parquet')
    L += [f'## Unresolved high-value loci: {len(un):,}', '']
    if un:
        by = Counter((r['status'], tuple(r['reasons'] or ())) for r in un)
        L += ['| Status | Reasons | Loci |', '|---|---|---|'] + [f"| {st} | {', '.join(rs)} | {n:,} |" for (st, rs), n in sorted(by.items(), key=lambda kv: -kv[1])]
        L.append('')
        L += read_evidence_section(un, c)
        reads_on = any(r.get('read_state') for r in un)
        L += [f'The first {min(top, len(un))} by position, with the exact evidence that would settle each:', '', '| Locus | Status | Rules | Missing evidence |' + (' Reads |' if reads_on else ''),
              '|---|---|---|---|' + ('---|' if reads_on else '')]
        for r in sorted(un, key=lambda r: (r['chrom'], r['pos'], r['key']))[:top]:
            reads = f" {(r.get('read_state') or 'not asked').replace('READS_', '')}" + (f" ({r['read_alt']}/{r['read_depth']} alternate)" if r.get('read_state') else '') + ' |' if reads_on else ''
            L.append(f"| `{r['key']}` | {r['status']} | {', '.join(r['rules'] or ())} | {(r['missing_evidence'] or ', '.join(r['reasons'] or ())).replace('|', '/')} |{reads}")
        L.append('')
    L += ['## Views', '', f"Calls considered: {vm['n_calls_considered']:,} (statuses {', '.join(vm['statuses'])}). A view selects and orders rows of small_calls.parquet; no view hides a call and none combines axes into one number. "
         "Evidence is attached to the site of the call: how many copies of a row's effect allele the sample carries (none, for a row about the reference allele and a homozygous alternate call) "
         "is in evidence_links.parquet (effect_allele_dosage), not in the views.", '']
    for name, v in vm['views'].items():
        L += [f"### {v['title']} ({v['n_calls']:,})", '', f"Ordering: {v['ordering']}.", '']
        rows = pq.read_table(os.path.join(results_dir, 'views', f'{name}.parquet')).to_pylist()[:top]
        for r in rows:
            r['axes'] = json.loads(r['axes_json']) if r.get('axes_json') else {}
            L.append(f"- `{r['key']}` {', '.join(r['gene_symbols'] or ()) or 'no gene'}, {r['zygosity'] if 'zygosity' in r and r['zygosity'] else r['gt']}"
                     f"{' (' + ', '.join(r['flags']) + ')' if r.get('flags') else ''}: {describe(r)}")
        if not rows:
            L.append('- none')
        L.append('')
    with open(os.path.join(results_dir, 'views', 'trait_view.json')) as f:
        traits = json.load(f)
    L += [f'### Traits named by effect-asserting evidence ({len(traits):,})', '', 'Identifiers are the sources\' own. The ten with the most calls:', '']
    for cid, e in sorted(traits.items(), key=lambda kv: (-kv[1]['n_calls'], kv[0]))[:10]:
        L.append(f"- `{cid}` {'; '.join(e['texts'][:2])}: {e['n_calls']} call(s) on {', '.join(e['axes'])}")
    rp = os.path.join(results_dir, 'views', 'trait_rollup.json')
    if os.path.exists(rp):
        with open(rp) as f:
            roll = json.load(f)
        crossed = roll.get('cross_referenced') or {}
        L += ['', f"Rolled up through the Open Targets ontology table into {len(roll['by_therapeutic_area'])} therapeutic areas: {len(roll['ancestors']) - len(crossed):,} trait id(s) through their own term, "
              f"{len(crossed):,} through a cross-reference or retired-term link of the table ({sum(1 for c in crossed.values() if c['many_to_one']):,} of them to more than one term, counted in the areas of each), "
              f"{len(roll['not_rolled_up']):,} not placed (no term holds, lists or replaces them) and kept as listed. A cross-reference groups; it does not assert that two concepts are the same. "
              'A trait in several areas is counted in each, so the rows do not add up to the number of traits:', '', '| Area | Traits | of them cross-referenced | Calls |', '|---|---:|---:|---:|']
        L += [f"| {g['name']} | {g['n_traits']:,} | {g.get('n_traits_cross_referenced', 0):,} | {g['n_calls']:,} |" for _a, g in sorted(roll['by_therapeutic_area'].items(), key=lambda kv: (-kv[1]['n_calls'], kv[0]))[:15]]
    with open(os.path.join(results_dir, 'views', 'pathway_view.json')) as f:
        paths = json.load(f)
    L += ['', f'### Pathways through the genes of matched calls ({len(paths):,})', '', 'Membership annotates genes; it is not evidence about any variant. The ten with the most calls:', '']
    for pid, e in sorted(paths.items(), key=lambda kv: (-kv[1]['n_calls'], kv[0]))[:10]:
        L.append(f"- `{pid}` {e['name']}: {e['n_calls']} call(s), {e['n_genes']} gene(s)")
    L.append('')
    text = '\n'.join(L)
    with open(os.path.join(results_dir, 'report.md'), 'w') as f:
        f.write(text)
    return text


CATALOGUE_VIEW = pa.schema([('key', pa.string()), ('chrom', pa.string()), ('pos', pa.int32()), ('ref', pa.string()), ('alt', pa.string()), ('rules', pa.list_(pa.string())),
                            ('n_evidence', pa.int32()), ('axes_json', pa.string()), ('evidence_ids', pa.list_(pa.string())), ('gene_ids', pa.list_(pa.string()))])


def write_catalogue_views(index, out_dir, rules=HIGH_VALUE_RULES):
    """The high-value variants of the catalogue itself, with no sample involved: every mapped variant that satisfies at least one rule, with the per-axis summary of all its
    variant-level evidence. One Parquet file per partition under <out_dir>/high_value_variants/, plus high_value_manifest.json (rule wording, counts). A variant may satisfy
    several rules and is listed once; the rules select, they do not rank, and the file is ordered by position."""
    from .normalise import parse_key
    from . import summary as _summary
    os.makedirs(os.path.join(out_dir, 'high_value_variants'), exist_ok=True)
    counts, total = Counter(), 0
    for part in index.variant_partitions():
        hv = high_value_keys(index, part, rules)
        if not hv:
            continue
        keys = sorted(hv, key=lambda k: (parse_key(k)[0], parse_key(k)[1], k))
        ev = index.evidence_for_variants(keys)
        genes = index.variant_genes(keys)
        rows = []
        for k in keys:
            c, pos, ref, alt = parse_key(k)
            rs = ev.get(k, [])
            rows.append(dict(key=k, chrom=c, pos=pos, ref=ref, alt=alt, rules=hv[k]['rules'], n_evidence=len(rs), axes_json=json.dumps(_summary.summarise(rs), sort_keys=True),
                             evidence_ids=sorted(r['evidence_id'] for r in rs), gene_ids=sorted({g['gene_id'] for g in genes.get(k, ())})))
            for r in hv[k]['rules']:
                counts[r] += 1
        pq.write_table(pa.Table.from_pylist(rows, schema=CATALOGUE_VIEW), os.path.join(out_dir, 'high_value_variants', f'{part}.parquet'), compression='zstd')
        total += len(rows)
    manifest = dict(variants=total, per_rule=dict(sorted(counts.items())), rules={r.name: r.description for r in rules},
                    note='A variant satisfying several rules is listed once. Gene-level and region-level evidence is not attached to these variants; nothing here is a score or a tier.')
    with open(os.path.join(out_dir, 'high_value_manifest.json'), 'w') as f:
        json.dump(manifest, f, indent=1, sort_keys=True)
    return manifest


def write_views(results_dir, index, statuses=('EXACT_MATCH',), ontology_dir=None):
    """Write views/<name>.parquet for every view, views/trait_view.json, views/pathway_view.json and views/views_manifest.json under the results directory.
    Views only select and order rows of small_calls; every call stays in small_calls.parquet whether or not any view lists it."""
    rows = list(read_calls(results_dir, statuses))
    schema = pq.read_schema(os.path.join(results_dir, 'small_calls.parquet'))
    vs = pa.schema([schema.field(c) for c in VIEW_COLUMNS])
    out = os.path.join(results_dir, 'views')
    os.makedirs(out, exist_ok=True)
    manifest = dict(statuses=list(statuses), n_calls_considered=len(rows), views={})
    for name, v in build_views(rows).items():
        pq.write_table(pa.Table.from_pylist([{c: r.get(c) for c in VIEW_COLUMNS} for r in v['rows']], schema=vs), os.path.join(out, f'{name}.parquet'), compression='zstd')
        manifest['views'][name] = dict(title=v['title'], ordering=v['ordering'], n_calls=len(v['rows']))
    traits = trait_view(rows)
    with open(os.path.join(out, 'trait_view.json'), 'w') as f:
        json.dump(traits, f, indent=1, sort_keys=True)
    genes = sorted({g for r in rows for g in (r.get('gene_ids') or ())})
    members = index.filtered('gene_set_member', F('gene_id').isin(genes)).to_pylist() if genes else []
    paths = pathway_view(rows, members)
    with open(os.path.join(out, 'pathway_view.json'), 'w') as f:
        json.dump(paths, f, indent=1, sort_keys=True)
    ontology = ontology_terms(ontology_dir)
    manifest.update(n_traits=len(traits), n_pathways=len(paths))
    if ontology:
        roll = trait_rollup(traits, ontology)
        with open(os.path.join(out, 'trait_rollup.json'), 'w') as f:
            json.dump(roll, f, indent=1, sort_keys=True)
        manifest.update(n_therapeutic_areas=len(roll['by_therapeutic_area']), n_traits_cross_referenced=len(roll['cross_referenced']), n_traits_not_rolled_up=len(roll['not_rolled_up']))
    manifest['note'] = ('trait ids are the sources\' own; trait_rollup.json groups them through the Open Targets ontology table (ancestors and therapeutic areas; ids the table does not hold are '
                        'grouped through its cross-references and retired terms where it has any, and the mapping is listed) and lists the ids it cannot place; '
                        'pathway membership annotates genes and is not evidence about any variant')
    with open(os.path.join(out, 'views_manifest.json'), 'w') as f:
        json.dump(manifest, f, indent=1, sort_keys=True)
    if os.path.exists(os.path.join(results_dir, 'run_manifest.json')):
        render_report(results_dir)
    return manifest
