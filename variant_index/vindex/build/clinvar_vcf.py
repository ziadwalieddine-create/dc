"""ClinVar GRCh38 VCF -> evidence rows (variant level) with cross-references and gene relations.

Rows per VCF record (the VCF is one record per variant allele; the ID column is the ClinVar VariationID):
  * germline aggregate classification (CLNSIG, CLNREVSTAT, CLNSIGCONF kept in extra_json), axis clinical_validity_actionability
  * oncogenicity classification (ONC) and somatic clinical impact (SCI) when present, origin somatic
  * population allele frequencies reported inside ClinVar (AF_ESP, AF_EXAC, AF_TGP), axis applicability_context
Records whose alleles do not validate against the reference stay in the index as unmapped/rejected rows with the reason.
Nothing is dropped: classes such as risk_factor, drug_response, Affects and conflicting stay as the source wrote them.
"""
import gzip
import json
import os

from ..normalise import normalise
from .common import DATA, evidence_row
from .classes import effect_of, stars, unknown_tokens, STARS   # noqa: F401  (effect_of/STARS stay importable from here)


def parse_info(s):
    d = {}
    for kv in s.split(';'):
        k, _, v = kv.partition('=')
        d[k] = v
    return d


def cond_ids(dis):
    ids = []
    for grp in dis.split('|'):
        for x in grp.split(','):
            if x and x != '.':
                ids.append(x)
    return sorted(set(ids))


def origin_of(code):
    try:
        v = int(code)
    except (TypeError, ValueError):
        return 'unknown'
    names = [n for bit, n in ((1, 'germline'), (2, 'somatic'), (4, 'inherited'), (8, 'paternal'), (16, 'maternal'), (32, 'de_novo'), (64, 'biparental'),
                              (128, 'uniparental'), (256, 'not_tested'), (512, 'tested_inconclusive')) if v & bit]
    return '+'.join(names) if names else 'unknown'


class ClinVarVcfAdapter:
    source_id = 'clinvar'
    targets = ['clinvar_vcf_grch38']

    def __init__(self, genes=None):
        self.genes = genes

    def run(self, ctx, path=None, limit=None):
        path = path or os.path.join(DATA, 'clinvar', 'clinvar.vcf.gz')
        ctx.adapter = self
        st = ctx.st()
        opener = gzip.open if path.endswith('.gz') else open
        version, genes = None, self.genes
        with opener(path, 'rt', encoding='utf-8') as f:
            hdr_lines = 0
            for ln, line in enumerate(f, start=1):
                if line.startswith('##'):
                    if line.startswith('##fileDate='):
                        version = line.strip().split('=', 1)[1]
                    continue
                if line.startswith('#'):
                    release = ctx.register_release('clinvar', 'clinvar_vcf_grch38', version or 'unknown')
                    continue
                if limit and st.c.get('records', 0) >= limit:
                    break
                st.inc('records')
                raw = line.rstrip('\n')
                f8 = raw.split('\t')
                chrom, pos, vid, ref, alt = f8[0], int(f8[1]), f8[2], f8[3], f8[4]
                info = parse_info(f8[7])
                canon = ctx.ref.contig(chrom)
                n = normalise(ctx.ref, canon, pos, ref, alt) if canon else None
                if n is not None and n.ok:
                    subject, mapping, reason = n.key, 'mapped', None
                else:
                    subject, mapping = None, ('rejected' if n is not None and n.status in ('ref_mismatch',) else 'unmapped')
                    reason = (n.status if n is not None else 'contig_unmapped')
                common = dict(chrom=(n.chrom if n is not None and n.ok else canon), pos=(n.pos if n is not None and n.ok else None),
                              unmapped_key=None if subject else f'{chrom}:{pos}:{ref}:{alt}', gene_ids=None, origin=origin_of(info.get('ORIGIN')))
                gene_ids = []
                for gi in (info.get('GENEINFO') or '').split('|'):
                    if ':' in gi:
                        sym, _, nid = gi.partition(':')
                        gene_ids.append(f'NCBIGene:{nid}')
                        if subject:
                            ens = genes.from_ncbi(nid) if genes else None
                            ctx.emit('variant_gene', dict(variant_key=subject, gene_id=ens or f'NCBIGene:{nid}', relation='clinvar_geneinfo',
                                                          consequence=info.get('MC', ''), release_id=release), name='clinvar')
                common['gene_ids'] = gene_ids
                acc = f'VariationID:{vid}|AlleleID:{info.get("ALLELEID", "")}'
                locator = f'L{ln}'
                lineage = f'clinvar:variation:{vid}'
                extra = {k: info[k] for k in ('CLNHGVS', 'CLNVC', 'CLNVCSO', 'MC', 'CLNVI', 'RS', 'CLNSIGCONF', 'CLNSIGSCV', 'DBVARID', 'ALLELEID') if k in info}
                for tag, ctype, origin_override in (('CLNSIG', 'germline_classification', None), ('ONC', 'oncogenicity_classification', 'somatic'),
                                                    ('SCI', 'somatic_clinical_impact', 'somatic')):
                    sig = info.get(tag)
                    if not sig:
                        continue
                    for t in unknown_tokens(sig):
                        st.inc(f'unrecognised_classification:{t}')
                    pre = {'CLNSIG': 'CLN', 'ONC': 'ONC', 'SCI': 'SCI'}[tag]
                    rev = info.get(pre + 'REVSTAT', '')
                    scv = (info.get(pre + 'SCV') or info.get('CLNSIGSCV') or '') if tag != 'CLNSIG' else info.get('CLNSIGSCV', '')
                    dn = info.get(pre + 'DN' if tag != 'CLNSIG' else 'CLNDN', '')
                    dis = info.get(pre + 'DISDB' if tag != 'CLNSIG' else 'CLNDISDB', '')
                    row = evidence_row(release, 'clinvar', acc, locator, raw, 'variant', subject, 'clinical_validity_actionability', 'curated_expert_classification',
                                       ctype, sig, effect_of(sig), mapping, reason, review_status=rev, review_stars=stars(rev), n_submitters=len([x for x in scv.split('|') if x]) or None,
                                       condition_ids=cond_ids(dis), condition_text=dn.replace('_', ' '), lineage_group=lineage,
                                       lineage_note='aggregate of ClinVar submissions (SCV); submitters may share cohorts', extra_json=json.dumps(extra, sort_keys=True), **{**common, **({'origin': origin_override} if origin_override else {})})
                    ctx.emit_evidence(row)
                if not any(info.get(t) for t in ('CLNSIG', 'ONC', 'SCI')):
                    # a ClinVar record with no classification of any kind still exists in ClinVar: keep that fact (it is not the same as being absent from ClinVar)
                    ctx.emit_evidence(evidence_row(release, 'clinvar', acc, locator, raw, 'variant', subject, 'clinical_validity_actionability', 'aggregated_integration', 'clinvar_record_without_classification',
                                                   'no_classification_in_vcf', 'not_applicable', mapping, reason, lineage_group=lineage, lineage_note='ClinVar record that carries no classification',
                                                   extra_json=json.dumps(extra, sort_keys=True), **common))
                    st.inc('records_without_classification')
                for tag, cohort in (('AF_ESP', 'GO-ESP'), ('AF_EXAC', 'ExAC'), ('AF_TGP', '1000G')):
                    if tag in info:
                        try:
                            af = float(info[tag])
                        except ValueError:
                            continue
                        row = evidence_row(release, 'clinvar', acc, locator, raw, 'variant', subject, 'applicability_context', 'population_observation',
                                           'population_allele_frequency', cohort, 'not_applicable', mapping, reason, discriminator=tag, score=af, score_unit='allele_frequency',
                                           eaf=af, lineage_group=f'cohort:{cohort}', lineage_note='frequency copied into ClinVar from the named cohort', **common)
                        ctx.emit_evidence(row)
                if subject:
                    ctx.offer('clinvar_variation', vid, subject)
                    if info.get('ALLELEID'):
                        ctx.offer('clinvar_allele', info['ALLELEID'], subject)
                    for rs in (info.get('RS') or '').split('|'):
                        if rs:
                            ctx.emit('xref', dict(namespace='dbsnp', identifier='rs' + rs, entity_level='variant', entity_id=subject, release_id=release, relation='clinvar_RS'), name='clinvar')
                            ctx.offer('dbsnp', 'rs' + rs, subject)
                    ctx.emit('xref', dict(namespace='clinvar_variation', identifier=vid, entity_level='variant', entity_id=subject, release_id=release, relation='vcf_id'), name='clinvar')
                    for vi in (info.get('CLNVI') or '').split('|'):
                        if vi.startswith('ClinGen:'):
                            ctx.offer('clingen_allele_registry', vi.split(':', 1)[1], subject)
                            ctx.emit('xref', dict(namespace='clingen_allele_registry', identifier=vi.split(':', 1)[1], entity_level='variant', entity_id=subject, release_id=release, relation='CLNVI'), name='clinvar')
                else:
                    ctx.offer('clinvar_variation_unmapped', vid, reason)        # this adapter kept the record as an unmapped/rejected row; later adapters must not repeat it
        st.inc('lines', ln if 'ln' in dir() else 0)
        return st.c
