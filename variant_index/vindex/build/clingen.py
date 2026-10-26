"""ClinGen -> evidence rows, each at its own level.

  ERepo variant classifications   level variant   VCEP assertion (effect / no_effect / uncertain) plus one row per applied ACMG/AMP criterion,
                                                  met or not met, placed on the axis the criterion speaks to. Variants are found through the
                                                  ClinVar Variation id, then the ClinGen Allele Registry id, both resolved against the ClinVar VCF
                                                  entities of the same build; HGVS alone is not resolved (unmapped, reason kept).
  Gene-disease validity           level gene      classification per gene and disease, moi, expert panel, date
  Dosage sensitivity (genes)      level gene      haploinsufficiency and triplosensitivity scores (0-3, 30 recessive, 40 unlikely)
  Dosage sensitivity (regions)    level region    the same scores for curated regions; recurrent-CNV regions are region entities
  Clinical actionability          level gene      consensus assertion and outcome/intervention scoring, adult and pediatric context

Gene-level rows stay on the gene. Nothing here is copied onto a variant.
"""
import csv
import json
import os
import re

from .common import DATA, evidence_row

CG = os.path.join(DATA, 'clingen')
RS_TAIL = 'A single record of a ClinGen expert panel; the underlying case reports and functional studies may also be cited by ClinVar submitters and in the literature'

PATHOGENIC_ASSERT = {'Pathogenic': 'effect', 'Likely Pathogenic': 'effect', 'Benign': 'no_effect', 'Likely Benign': 'no_effect', 'Uncertain Significance': 'uncertain'}
# ACMG/AMP criterion -> (axis, nature of the evidence the criterion speaks to)
CODE_AXIS = {
    'PVS1': ('mechanistic_confidence', 'null_variant_in_lof_gene'), 'PS1': ('mechanistic_confidence', 'same_amino_acid_change'), 'PM5': ('mechanistic_confidence', 'novel_missense_at_known_residue'),
    'PM1': ('mechanistic_confidence', 'hotspot_or_functional_domain'), 'PM4': ('mechanistic_confidence', 'protein_length_change'), 'PP2': ('mechanistic_confidence', 'missense_constraint'),
    'PP3': ('mechanistic_confidence', 'computational_evidence'), 'BP1': ('mechanistic_confidence', 'missense_in_truncating_gene'), 'BP3': ('mechanistic_confidence', 'inframe_in_repeat'),
    'BP4': ('mechanistic_confidence', 'computational_evidence'), 'BP7': ('mechanistic_confidence', 'synonymous_no_splice_impact'),
    'PS3': ('empirical_functional_support', 'functional_studies'), 'BS3': ('empirical_functional_support', 'functional_studies'),
    'PS2': ('association_replication', 'de_novo'), 'PM6': ('association_replication', 'assumed_de_novo'), 'PS4': ('association_replication', 'prevalence_in_affected'),
    'PP1': ('association_replication', 'segregation'), 'BS4': ('association_replication', 'non_segregation'), 'PM3': ('association_replication', 'in_trans_with_pathogenic'),
    'BP2': ('association_replication', 'observed_in_trans_or_cis'), 'PP4': ('association_replication', 'phenotype_specificity'), 'BP5': ('association_replication', 'alternate_molecular_basis'),
    'PM2': ('applicability_context', 'population_frequency'), 'BA1': ('applicability_context', 'population_frequency'), 'BS1': ('applicability_context', 'population_frequency'),
    'BS2': ('applicability_context', 'observed_in_healthy_adults'), 'PP5': ('clinical_validity_actionability', 'reputable_source'), 'BP6': ('clinical_validity_actionability', 'reputable_source'),
}


def code_family(code):
    return re.split(r'[_:\s(]', code.strip())[0]


def split_codes(s):
    return [c.strip() for c in re.split(r'[,;]', s or '') if c.strip()]


def resolve_variants(ctx, cv_id, ca_id):
    keys = ctx.found.get(('clinvar_variation', cv_id)) or ctx.found.get(('clingen_allele_registry', ca_id))
    return sorted(set(keys)) if keys else None


class ClinGenErepoAdapter:
    source_id = 'clingen'
    targets = ['clingen_erepo_variant_classifications']
    requires = ('clinvar',)      # variants are resolved through the ClinVar entities of the same build

    def __init__(self, path=None):
        self.path = path or os.path.join(CG, 'erepo_variant_classifications.tsv')

    def prepare(self, ctx):
        ids_cv, ids_ca = set(), set()
        with open(self.path, encoding='utf-8', newline='') as f:
            for r in csv.DictReader(f, delimiter='\t', quoting=csv.QUOTE_NONE):
                if (r['ClinVar Variation Id'] or '').isdigit():
                    ids_cv.add(r['ClinVar Variation Id'])
                if r['Allele Registry Id']:
                    ids_ca.add(r['Allele Registry Id'])
        ctx.want('clinvar_variation', ids_cv)
        ctx.want('clingen_allele_registry', ids_ca)

    def run(self, ctx, limit=None):
        ctx.adapter = self
        st = ctx.st()
        rel = None
        with open(self.path, encoding='utf-8', newline='') as f:
            lines = f.read().split('\n')
        hdr = lines[0].split('\t')
        rel = ctx.register_release('clingen', 'clingen_erepo_variant_classifications', 'live export 2026-10')
        for ln, line in enumerate(lines[1:], start=2):
            if not line.strip():
                continue
            if limit and st.c.get('records', 0) >= limit:
                break
            st.inc('records')
            r = dict(zip(hdr, line.split('\t')))
            keys = resolve_variants(ctx, r['ClinVar Variation Id'], r['Allele Registry Id'])
            if keys is None:
                keys, mapping, reason = [None], 'unmapped', ('no_clinvar_or_allele_registry_match' if (r['ClinVar Variation Id'] or r['Allele Registry Id']) else 'hgvs_only')
            elif len(keys) > 1:
                keys, mapping, reason = [None], 'ambiguous', 'identifier_resolves_to_several_alleles'
            else:
                mapping, reason = 'mapped', None
            key = keys[0]
            acc = f"ERepo:{r['Uuid']}"
            assertion = r['Assertion']
            mondo = [r['Mondo Id']] if r['Mondo Id'] else []
            chrom = pos = None
            if key:
                c, p, _, _ = key.split(':')
                chrom, pos = c, int(p)
            unmapped_key = None if key else (f"ClinVar:{r['ClinVar Variation Id']}" if r['ClinVar Variation Id'] else (f"CA:{r['Allele Registry Id']}" if r['Allele Registry Id'] else r['Variation']))
            common = dict(chrom=chrom, pos=pos, unmapped_key=unmapped_key, condition_ids=mondo, condition_text=r['Disease'], inheritance=r['Mode of Inheritance'] or None,
                          curator=r['Expert Panel'], evaluated=r['Approval Date'] or None, assertion_criteria=r['Guideline'] or None, gene_ids=[r['HGNC Gene Symbol']] if r['HGNC Gene Symbol'] else None,
                          lineage_group=f"clingen_erepo:{r['Uuid']}", lineage_note=RS_TAIL, origin='germline')
            met, notmet = split_codes(r['Applied Evidence Codes (Met)']), split_codes(r['Applied Evidence Codes (Not Met)'])
            extra = json.dumps(dict(met=met, not_met=notmet, pubmed=r['PubMed Articles'], link=r['Evidence Repo Link'], retracted=r['Retracted'], published=r['Published Date'],
                                    summary=(r['Summary of interpretation'] or '')[:1500]), sort_keys=True)
            ctx.emit_evidence(evidence_row(rel, 'clingen', acc, f'L{ln}', line, 'variant', key, 'clinical_validity_actionability', 'curated_expert_classification', 'vcep_classification', assertion,
                                           PATHOGENIC_ASSERT.get(assertion, 'uncertain'), mapping, reason, review_status='ClinGen expert panel', extra_json=extra, **common))
            for status, codes in (('met', met), ('not_met', notmet)):
                for code in codes:
                    fam = code_family(code)
                    axis, nature = CODE_AXIS.get(fam, ('clinical_validity_actionability', 'other_criterion'))
                    if status == 'met':
                        eff = 'effect' if fam.startswith('P') else ('no_effect' if fam.startswith('B') else 'uncertain')
                    else:
                        eff = 'not_applicable'
                    ctx.emit_evidence(evidence_row(rel, 'clingen', acc, f'L{ln}', line, 'variant', key, axis, 'curated_expert_classification', f'acmg_criterion:{code}', status, eff, mapping, reason,
                                                   discriminator=f'{code}|{status}', extra_json=json.dumps(dict(nature=nature, criterion_family=fam)), **common))
        return st.c


def _vrows(path):
    with open(path, encoding='utf-8', newline='') as f:
        hdr, out = None, []
        for ln, row in enumerate(csv.reader(f), start=1):
            if hdr is None:
                if row and row[0] == 'GENE SYMBOL':
                    hdr = row
                continue
            if not row or row[0].startswith('+'):
                continue
            out.append((ln, dict(zip(hdr, row)), json.dumps(row, ensure_ascii=False)))
    return out


VALIDITY_EFFECT = {'Definitive': 'effect', 'Strong': 'effect', 'Moderate': 'effect', 'Limited': 'uncertain', 'Disputed': 'conflicting', 'Refuted': 'no_effect', 'No Known Disease Relationship': 'no_effect'}
DOSAGE_EFFECT = {'3': 'effect', '2': 'uncertain', '1': 'uncertain', '0': 'uncertain', '30': 'no_effect', '40': 'no_effect'}
DOSAGE_TEXT = {'3': 'Sufficient evidence', '2': 'Emerging evidence', '1': 'Little evidence', '0': 'No evidence available', '30': 'Gene associated with autosomal recessive phenotype', '40': 'Dosage sensitivity unlikely'}
ACTION_EFFECT = {'Definitive Actionability': 'effect', 'Strong Actionability': 'effect', 'Moderate Actionability': 'effect', 'Limited Actionability': 'uncertain'}


class ClinGenGeneRegionAdapter:
    source_id = 'clingen'
    targets = ['clingen_gene_validity', 'clingen_gene_curation_list', 'clingen_region_curation_list', 'clingen_recurrent_cnv_v2_1',
               'clingen_actionability_adult_assertions', 'clingen_actionability_adult_scoring', 'clingen_actionability_pediatric_assertions', 'clingen_actionability_pediatric_scoring']

    def __init__(self, genes, root=None):
        self.genes, self.root = genes, root or CG

    def run(self, ctx, limit=None):
        ctx.adapter = self
        st = ctx.st()
        self.validity(ctx, st)
        self.dosage_genes(ctx, st)
        self.dosage_regions(ctx, st)
        self.recurrent(ctx, st)
        for context in ('adult', 'pediatric'):
            self.actionability(ctx, st, context)
        return st.c

    def gene_row(self, ctx, rel, acc, locator, raw, gene, axis, ctype, cval, eff, how, **kw):
        if gene is None:
            return evidence_row(rel, 'clingen', acc, locator, raw, 'gene', None, axis, 'curated_expert_classification', ctype, cval, eff, 'unmapped', f'gene_{how}', **kw)
        return evidence_row(rel, 'clingen', acc, locator, raw, 'gene', gene, axis, 'curated_expert_classification', ctype, cval, eff, 'mapped', None, **kw)

    def validity(self, ctx, st):
        rel = ctx.register_release('clingen', 'clingen_gene_validity', 'live export 2026-10-25')
        for ln, r, raw in _vrows(os.path.join(self.root, 'clingen_gene_disease_validity.csv')):
            st.inc('validity')
            g, how = self.genes.resolve(hgnc=r['GENE ID (HGNC)'], symbol=r['GENE SYMBOL'])
            acc = re.sub(r'.*/', '', r['ONLINE REPORT']) or f"{r['GENE SYMBOL']}|{r['DISEASE ID (MONDO)']}"
            ctx.emit_evidence(self.gene_row(ctx, rel, acc, f'L{ln}', raw, g, 'clinical_validity_actionability', 'gene_disease_validity', r['CLASSIFICATION'],
                                            VALIDITY_EFFECT.get(r['CLASSIFICATION'], 'uncertain'), how, condition_ids=[r['DISEASE ID (MONDO)']], condition_text=r['DISEASE LABEL'],
                                            inheritance=r['MOI'], curator=r['GCEP'], evaluated=r['CLASSIFICATION DATE'], assertion_criteria=r['SOP'],
                                            lineage_group=f'clingen_validity:{acc}', extra_json=json.dumps(dict(symbol=r['GENE SYMBOL'], report=r['ONLINE REPORT']))))

    def _dosage_rows(self, ctx, st, rel, level, subject, accession_base, locator, raw, row, name):
        for which, sc_key, desc_key, dis_key, pm_prefix in (('haploinsufficiency', 'Haploinsufficiency Score', 'Haploinsufficiency Description', 'Haploinsufficiency Disease ID', 'Haploinsufficiency PMID'),
                                                            ('triplosensitivity', 'Triplosensitivity Score', 'Triplosensitivity Description', 'Triplosensitivity Disease ID', 'Triplosensitivity PMID')):
            score = (row.get(sc_key) or '').strip()
            if not score or not score.lstrip('-').isdigit():
                continue
            pm = [row.get(f'{pm_prefix}{i}') for i in range(1, 7) if row.get(f'{pm_prefix}{i}')]
            dis = row.get(dis_key) or ''
            common = dict(claim_num=float(score), condition_ids=[dis] if dis else None, evaluated=row.get('Date Last Evaluated'), curator='ClinGen Dosage Sensitivity Working Group',
                          lineage_group=f'clingen_dosage:{accession_base}:{which}', extra_json=json.dumps(dict(description=row.get(desc_key) or DOSAGE_TEXT.get(score, ''), pmids=pm, name=name)))
            axis, kind = 'clinical_validity_actionability', 'curated_expert_classification'
            if subject is None:
                ctx.emit_evidence(evidence_row(rel, 'clingen', f'{accession_base}:{which}', locator, raw, level, None, axis, kind, f'{which}_score', DOSAGE_TEXT.get(score, score),
                                               DOSAGE_EFFECT.get(score, 'uncertain'), 'unmapped', 'gene_unresolved' if level == 'gene' else 'region_location_not_parsed',
                                               discriminator=which, unmapped_key=name, **common))
            else:
                ctx.emit_evidence(evidence_row(rel, 'clingen', f'{accession_base}:{which}', locator, raw, level, subject, axis, kind, f'{which}_score', DOSAGE_TEXT.get(score, score),
                                               DOSAGE_EFFECT.get(score, 'uncertain'), 'mapped', None, discriminator=which, **common))

    def _tsv_rows(self, path):
        hdr, out = None, []
        with open(path, encoding='utf-8', newline='') as f:
            for ln, line in enumerate(f, start=1):
                line = line.rstrip('\n')
                if line.startswith('#') and ('Gene Symbol' in line or 'ISCA ID' in line) and '\t' in line:
                    hdr = line.lstrip('#').split('\t')
                    continue
                if line.startswith('#') or not line.strip() or hdr is None:
                    continue
                out.append((ln, dict(zip(hdr, line.split('\t'))), line))
        return out

    def dosage_genes(self, ctx, st):
        rel = ctx.register_release('clingen', 'clingen_gene_curation_list', 'daily export 2026-10-25')
        for ln, r, raw in self._tsv_rows(os.path.join(self.root, 'ClinGen_gene_curation_list_GRCh38.tsv')):
            st.inc('dosage_genes')
            g, how = self.genes.resolve(symbol=r['Gene Symbol'], ncbi=(r.get('Gene ID') or '').strip())
            self._dosage_rows(ctx, st, rel, 'gene', g, f"gene:{r['Gene Symbol']}", f'L{ln}', raw, r, r['Gene Symbol'])

    def _region_bounds(self, loc):
        m = re.match(r'^(?:chr)?([0-9XYM]+):(\d+)-(\d+)$', (loc or '').strip())
        return (m.group(1).replace('M', 'MT') if m else None, int(m.group(2)) if m else None, int(m.group(3)) if m else None)

    def dosage_regions(self, ctx, st):
        rel = ctx.register_release('clingen', 'clingen_region_curation_list', 'daily export 2026-04-28')
        for ln, r, raw in self._tsv_rows(os.path.join(self.root, 'ClinGen_region_curation_list_GRCh38.tsv')):
            st.inc('dosage_regions')
            c, s, e = self._region_bounds(r['Genomic Location'])
            rid = r['ISCA ID']
            if c is None:
                self._dosage_rows(ctx, st, rel, 'region', None, f'region:{rid}', f'L{ln}', raw, r, r['ISCA Region Name'])
                continue
            ctx.emit('entity_region', dict(region_id=rid, kind='clingen_dosage_region', chrom=c, start=s, end=e, name=r['ISCA Region Name'], gene_ids=None, release_id=rel), name='clingen')
            self._dosage_rows(ctx, st, rel, 'region', rid, f'region:{rid}', f'L{ln}', raw, dict(r), r['ISCA Region Name'])

    def recurrent(self, ctx, st):
        rel = ctx.register_release('clingen', 'clingen_recurrent_cnv_v2_1', 'V2.1')
        with open(os.path.join(self.root, 'ClinGen_recurrent_CNV_V2.1-hg38.bed'), encoding='utf-8') as f:
            for ln, line in enumerate(f, start=1):
                if line.startswith(('track', '#')) or not line.strip():
                    continue
                c, s, e, name = line.rstrip('\n').split('\t')[:4]
                canon = c.replace('chr', '').replace('M', 'MT')
                ctx.emit('entity_region', dict(region_id=f'clingen_recurrent:{name}:{canon}:{int(s) + 1}-{e}', kind='clingen_recurrent_cnv', chrom=canon, start=int(s) + 1, end=int(e), name=name, gene_ids=None, release_id=rel), name='clingen')
                st.inc('recurrent_cnv_regions')

    def actionability(self, ctx, st, context):
        rel_a = ctx.register_release('clingen', f'clingen_actionability_{context}_assertions', 'live API 2026-10')
        rel_s = ctx.register_release('clingen', f'clingen_actionability_{context}_scoring', 'live API 2026-10')
        for fn, rel, kind in ((f'{context}_assertions_all.tsv', rel_a, 'assertion'), (f'{context}_scoring.tsv', rel_s, 'scoring')):
            with open(os.path.join(self.root, 'actionability', fn), encoding='utf-8', newline='') as f:
                hdr = None
                for ln, line in enumerate(f, start=1):
                    line = line.rstrip('\n')
                    if hdr is None:
                        hdr = line.lstrip('# ').split('\t')
                        continue
                    if not line.strip():
                        continue
                    r = dict(zip(hdr, line.split('\t')))
                    st.inc(f'actionability_{kind}')
                    subj = r['geneOrVariant']
                    variant_level = r['curationType'] != 'Gene-Condition'
                    g, how = (None, 'variant_curation') if variant_level else self.genes.resolve(symbol=subj)
                    common = dict(condition_ids=[x for x in (f"MONDO:{r['mondo'].split('_')[-1]}" if r.get('mondo', '').startswith('MONDO') else r.get('mondo'),) if x], condition_text=r['disease'],
                                  curator=r.get('scorer') or 'ClinGen Actionability Working Group', evaluated=r.get('releaseDate') or r.get('lastUpdated'),
                                  lineage_group=f"clingen_actionability:{r['docId']}", extra_json=json.dumps({k: v for k, v in r.items() if v and k not in ('disease',)}, sort_keys=True))
                    acc = f"{context}:{r['docId']}:{kind}:{ln}"
                    if kind == 'assertion':
                        val = r.get('consensusAssertion') or r.get('preliminaryAssertion') or r.get('suggestedAssertion') or ''
                        ctype, cval, eff = f'clinical_actionability_assertion_{context}', val or 'none recorded', ACTION_EFFECT.get(val, 'not_applicable')
                    else:
                        ctype, cval, eff = f'clinical_actionability_scoring_{context}', r.get('final-overall') or '', 'not_applicable'
                    if variant_level:
                        ctx.emit_evidence(evidence_row(rel, 'clingen', acc, f'L{ln}', line, 'variant', None, 'clinical_validity_actionability', 'curated_expert_classification', ctype, cval, eff,
                                                       'unmapped', 'variant_named_by_hgvs_only', unmapped_key=subj, **common))
                    else:
                        ctx.emit_evidence(self.gene_row(ctx, rel, acc, f'L{ln}', line, g, 'clinical_validity_actionability', ctype, cval, eff, how, **common))
