"""Hand-made evidence rows for matcher tests. Every row is validated against the schema vocabulary."""
import hashlib

from vindex import schema


def ev(key, source='testsrc', axis='clinical_validity_actionability', assertion='effect', claim_type='germline_classification', claim='Pathogenic',
       kind='curated_expert_classification', lineage=None, accession=None, level='variant', **kw):
    acc = accession or f'{source}:{key}:{claim}'
    row = dict(evidence_id=hashlib.sha1(f'{source}|{acc}|{axis}|{claim_type}'.encode()).hexdigest()[:16], source_id=source, release_id=f'{source}@t:file',
               accession=acc, locator='line:1', raw_sha256='0' * 64, level=level, subject_id=key, axis=axis, evidence_kind=kind,
               claim_type=claim_type, claim_value=claim, effect_assertion=assertion, lineage_group=lineage, mapping_status='mapped')
    if level == 'variant':
        c, p, r, a = key.split(':')
        row.update(chrom=c, pos=int(p))
    row.update(kw)
    bad = schema.validate_evidence_row(row)
    assert not bad, bad
    return row


def index_from(rows, variants=(), genes=(), regions=()):
    from vindex.index_reader import Index
    ev_t = schema.rows_to_table('evidence', rows)
    kw = dict(evidence=ev_t)
    if variants:
        kw['entity_variant'] = schema.rows_to_table('entity_variant', [dict(variant_key=k, chrom=k.split(':')[0], pos=int(k.split(':')[1]), ref=k.split(':')[2], alt=k.split(':')[3],
                                                                         vtype='mnv' if len(k.split(':')[2]) == len(k.split(':')[3]) > 1 else 'snv') for k in variants])
    if genes:
        kw['entity_gene'] = schema.rows_to_table('entity_gene', list(genes))
    if regions:
        kw['entity_region'] = schema.rows_to_table('entity_region', list(regions))
    return Index.from_tables(**kw)
