"""Rules and views over the complete catalogue of matches.

High-value rules say which indexed variants are worth asking "what does this sample say about this allele?" even when the sample's file holds no record
of it (the unresolved-loci report). They are filters over evidence rows, each stated in plain words, and a variant may satisfy several. They are not
scores and they do not rank.

A view is a filter and an explicit ordering over the matches of one sample. It never changes or replaces the per-axis evidence, no number here combines
axes, and a call that fails a filter stays in the full call table. The ordering of each view is part of its definition and is printed with it.
"""
import json
import os
from dataclasses import dataclass

import pyarrow.compute as pc
import pyarrow.parquet as pq

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
    """{variant_key: {'rules': [names], 'evidence_ids': [...]}} for the mapped variant-level evidence of one chromosome that satisfies a rule."""
    out = {}
    for rule in rules:
        expr = rule.expr & (F('level') == 'variant') & (F('mapping_status') == 'mapped')
        sub = os.path.join('variant', chrom)
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
