"""Open Targets 26.09 gene-level perturbation and burden evidence. Every row here is attached to the gene; none of it is ever copied onto a variant.

  DepMap essentiality (target_essentiality)   CRISPR knockout across cancer cell lines: one row per gene (essential in the pan-cancer call or not, share of
                                              screens with gene effect <= -1) and the tissues in which at least half of >= 3 screens depend on the gene.
  IMPC / MGI mouse models (evidence_impc)     one row per mouse model (gene x allelic composition x genetic background): zygosity read from the allelic
                                              composition, the mouse phenotypes the model shows and the human diseases the platform links it to.
                                              The source lists only phenotypes that were observed; an absent phenotype is not a negative result.
  CRISPR screens (evidence_crispr, _screen)   one row per gene x study contrast with library, cell type, contrast, log2 fold change and the resource's p-value.
  Gene burden (evidence_gene_burden)          rare-variant burden tests of qualifying variants in a gene: odds ratio or beta with interval, p-value, cohort,
                                              ancestry, allelic requirement.

A gene knockout or burden of qualifying variants says nothing about one particular allele; the perturbation field states what was done to the gene.
"""
import glob
import json
import os
import re
from collections import Counter

import pyarrow as pa
import pyarrow.parquet as pq

from .common import evidence_row
from .opentargets_genetics import OT

ALLELE = re.compile(r'^(?P<g1>[^<]+)<(?P<a1>[^>]*)>/(?P<g2>[^<]+)<(?P<a2>[^>]*)>$')


def zygosity(composition):
    """'Bmp4<tm1Blh>/Bmp4<+>' -> heterozygous (one wild-type '+' allele); '<a>/<a>' -> homozygous; 'X<em1> hom early' style -> as written."""
    c = (composition or '').strip()
    m = ALLELE.match(c)
    if m:
        a1, a2 = m.group('a1'), m.group('a2')
        if '+' in (a1, a2) and a1 != a2:
            return 'heterozygous'
        return 'homozygous' if a1 == a2 else 'compound heterozygous or hemizygous'
    low = c.lower()
    for w, z in (('hom', 'homozygous'), ('het', 'heterozygous'), ('hemi', 'hemizygous')):
        if re.search(rf'\b{w}', low):
            return z
    return 'unspecified'


def p_ok(v):
    return v if v is not None and 0 <= v <= 1 else None


def impc_locator(gene, composition, background):
    """The summary of one mouse model is built from every IMPC row with this target, allelic composition and genetic background, so the locator names the group: a verifier
    can recompute the summary from the rows of the files that the key selects."""
    return 'evidence_impc#' + json.dumps([gene, composition, background])


class OpenTargetsGeneAdapter:
    source_id = 'opentargets'
    targets = ['opentargets_26_09_target_essentiality', 'opentargets_26_09_evidence_impc', 'opentargets_26_09_evidence_crispr', 'opentargets_26_09_evidence_crispr_screen',
               'opentargets_26_09_evidence_gene_burden']

    def __init__(self, genes, root=None):
        self.genes, self.root = genes, root or OT

    def files(self, name):
        return sorted(glob.glob(os.path.join(self.root, name, '*.parquet')))

    def table(self, name, columns):
        ts = [pq.read_table(f, columns=columns) for f in self.files(name)]
        return pa.concat_tables(ts) if ts else None

    def subject(self, gene_id):
        known = getattr(self.genes, 'by_ensg', None)
        ok = bool(gene_id) and (known is None or gene_id in known)
        return (gene_id, 'mapped', None) if ok else (None, 'unmapped', 'gene_not_in_gene_map')

    def emit(self, ctx, rel, st, acc, locator, raw, gene, axis, kind, claim_type, claim_value, effect, **kw):
        gid, status, reason = self.subject(gene)
        extra = kw.pop('extra', {}) or {}
        ctx.emit_evidence(evidence_row(rel, 'opentargets', acc, locator, raw, 'gene', gid, axis, kind, claim_type, claim_value, effect, status, reason, discriminator=kw.pop('discriminator', ''),
                                       gene_ids=[gene] if gene else None, unmapped_key=None if gid else gene, extra_json=json.dumps(extra, sort_keys=True, default=str), **kw))
        st.inc(f'{claim_type}:{status}')

    def run(self, ctx, limit=None):
        ctx.adapter = self
        st = ctx.st()
        self.essentiality(ctx, st)
        self.impc(ctx, st)
        self.crispr(ctx, st)
        self.burden(ctx, st)
        return st.c

    # ---- DepMap -----------------------------------------------------------------------------------------------------------------
    def essentiality(self, ctx, st):
        rel = ctx.register_release('opentargets', 'opentargets_26_09_target_essentiality', '26.09')
        for f in self.files('target_essentiality'):
            t = pq.read_table(f)
            for i, r in enumerate(t.to_pylist()):
                screens_all, tissues = 0, []
                for ts in r['depMapEssentiality'] or ():
                    effects = [s['geneEffect'] for s in ts['screens'] or () if s.get('geneEffect') is not None]
                    screens_all += len(effects)
                    dep = sum(1 for e in effects if e <= -1.0)
                    if len(effects) >= 3 and dep / len(effects) >= 0.5:
                        tissues.append(dict(tissue=ts['tissueName'], screens=len(effects), dependent=dep))
                allv = [s['geneEffect'] for ts in (r['depMapEssentiality'] or ()) for s in (ts['screens'] or ()) if s.get('geneEffect') is not None]
                frac = (sum(1 for e in allv if e <= -1.0) / len(allv)) if allv else None
                raw = json.dumps(dict(targetId=r['targetId'], isEssential=r['isEssential'], n_screens=len(allv), frac=frac), sort_keys=True)
                self.emit(ctx, rel, st, f"depmap:{r['targetId']}", f"target_essentiality/{os.path.basename(f)}#{i}", raw, r['targetId'], 'empirical_functional_support', 'perturbation_experimental',
                          'crispr_knockout_essentiality', 'essential' if r['isEssential'] else 'not pan-essential', 'effect' if r['isEssential'] else 'no_effect', perturbation='CRISPR knockout (DepMap)',
                          model_system='human cancer cell lines', score=frac, score_unit='fraction_of_screens_with_gene_effect<=-1', sample_size=len(allv) or None, curator='DepMap via Open Targets',
                          lineage_group='depmap', lineage_note='one DepMap release; tissue rows are subsets of the same screens',
                          extra=dict(is_pan_essential=r['isEssential'], screens=len(allv), dependent_tissues=tissues, note='essentiality in cell lines is not evidence about a particular allele'))
                for tdep in tissues:
                    self.emit(ctx, rel, st, f"depmap:{r['targetId']}:{tdep['tissue']}", f"target_essentiality/{os.path.basename(f)}#{i}", raw, r['targetId'], 'empirical_functional_support', 'perturbation_experimental',
                              'crispr_knockout_tissue_dependency', f"{tdep['dependent']}/{tdep['screens']} screens", 'effect', perturbation='CRISPR knockout (DepMap)', tissue=tdep['tissue'],
                              model_system='human cancer cell lines', score=tdep['dependent'] / tdep['screens'], score_unit='fraction_of_screens_with_gene_effect<=-1', sample_size=tdep['screens'],
                              curator='DepMap via Open Targets', lineage_group='depmap', lineage_note='subset of the gene-level DepMap row', discriminator=tdep['tissue'])

    # ---- IMPC / MGI -------------------------------------------------------------------------------------------------------------
    def impc(self, ctx, st):
        rel = ctx.register_release('opentargets', 'opentargets_26_09_evidence_impc', '26.09')
        cols = ['targetId', 'biologicalModelAllelicComposition', 'biologicalModelGeneticBackground', 'biologicalModelId', 'diseaseId', 'diseaseFromSource', 'resourceScore',
                'diseaseModelAssociatedModelPhenotypes', 'literature', 'directionOnTrait']
        files = self.files('evidence_impc')
        if not files:
            return
        models = {}
        for f in files:                                  # 6.9M rows with 107M phenotype entries: stream in batches, keep only the per-model aggregates (28k models)
            for batch in pq.ParquetFile(f).iter_batches(batch_size=20_000, columns=cols):
                for r in batch.to_pylist():
                    k = (r['targetId'], r['biologicalModelAllelicComposition'], r['biologicalModelGeneticBackground'])
                    m = models.setdefault(k, dict(model_ids=set(), diseases={}, phenos=Counter(), lit=set(), rows=0))
                    m['rows'] += 1
                    if r['biologicalModelId']:
                        m['model_ids'].add(r['biologicalModelId'])
                    if r['diseaseId']:
                        m['diseases'][r['diseaseId']] = max(m['diseases'].get(r['diseaseId'], 0), r['resourceScore'] or 0)
                    for p in r['diseaseModelAssociatedModelPhenotypes'] or ():
                        m['phenos'][(p['id'], p['label'])] += 1
                    m['lit'].update(r['literature'] or ())
        for (gene, comp, bg), m in sorted(models.items(), key=lambda kv: tuple(str(x) for x in kv[0])):
            z = zygosity(comp)
            top = [f'{lab}' for (_i, lab), _n in m['phenos'].most_common(5)]
            raw = json.dumps(dict(targetId=gene, composition=comp, background=bg, rows=m['rows'], diseases=sorted(m['diseases'])), sort_keys=True)
            self.emit(ctx, rel, st, f"impc:{gene}:{comp}:{bg}", impc_locator(gene, comp, bg), raw, gene, 'empirical_functional_support', 'perturbation_experimental', 'mouse_model_phenotypes',
                      '; '.join(top) or 'no phenotype terms', 'effect', perturbation=f'mouse allele {comp} ({z})', model_system=f'mouse {bg}'.strip(), curator='IMPC / MGI via Open Targets',
                      condition_ids=sorted(m['diseases'])[:50] or None, inheritance=z, lineage_group=f"mgi_model:{sorted(m['model_ids'])[0] if m['model_ids'] else comp}",
                      lineage_note='a mouse model counts once however many human diseases it is matched to', discriminator=f'{comp}|{bg}',
                      extra=dict(zygosity=z, phenotype_counts=[dict(id=i, label=lab, n=n) for (i, lab), n in m['phenos'].most_common(40)], n_phenotypes=len(m['phenos']),
                                 disease_links={d: s for d, s in sorted(m['diseases'].items(), key=lambda kv: -kv[1])[:30]}, literature=sorted(m['lit'])[:20], source_rows=m['rows'],
                                 note='observed phenotypes only; an absent phenotype was not reported, which is not a negative result'))

    # ---- CRISPR screens ---------------------------------------------------------------------------------------------------------
    def crispr(self, ctx, st):
        rel_c = ctx.register_release('opentargets', 'opentargets_26_09_evidence_crispr', '26.09')
        rel_s = ctx.register_release('opentargets', 'opentargets_26_09_evidence_crispr_screen', '26.09')
        for f in self.files('evidence_crispr'):
            for i, r in enumerate(pq.read_table(f).to_pylist()):
                lines = [c.get('name') for c in r.get('diseaseCellLines') or ()]
                self.emit(ctx, rel_c, st, f"crispr:{r['id']}", f"evidence_crispr/{os.path.basename(f)}#{i}", json.dumps({k: r[k] for k in ('id', 'targetId', 'diseaseId', 'resourceScore')}, sort_keys=True, default=str), r['targetId'],
                          'empirical_functional_support', 'perturbation_experimental', 'crispr_screen_hit', r.get('diseaseFromSource') or r.get('diseaseId') or 'screen hit', 'effect', perturbation='CRISPR knockout screen',
                          model_system='; '.join(sorted(set(x for x in lines if x))[:5]) or None, p_value=p_ok(r.get('resourceScore')), condition_ids=[r['diseaseId']] if r.get('diseaseId') else None,
                          condition_text=r.get('diseaseFromSource'), curator=r.get('datasourceId'), lineage_group=f"pmid:{r['literature'][0]}" if r.get('literature') else f"crispr:{r['id']}",
                          discriminator=r['id'], extra=dict(cell_lines=lines[:20], literature=r.get('literature'), publication_date=r.get('publicationDate')))
        for f in self.files('evidence_crispr_screen'):
            for i, r in enumerate(pq.read_table(f).to_pylist()):
                self.emit(ctx, rel_s, st, f"crispr_screen:{r['id']}", f"evidence_crispr_screen/{os.path.basename(f)}#{i}", json.dumps({k: r[k] for k in ('id', 'targetId', 'studyId', 'log2FoldChangeValue', 'resourceScore')}, sort_keys=True, default=str),
                          r['targetId'], 'empirical_functional_support', 'perturbation_experimental', 'crispr_screen_hit', f"log2FC {r['log2FoldChangeValue']:.3g}" if r.get('log2FoldChangeValue') is not None else 'screen hit',
                          'effect', perturbation=f"CRISPR screen {r.get('crisprScreenLibrary') or ''}; {r.get('contrast') or ''}".strip('; '), model_system=r.get('geneticBackground'), cell_type=r.get('cellType'),
                          claim_num=r.get('log2FoldChangeValue'), p_value=p_ok(r.get('resourceScore')), condition_ids=[r['diseaseId']] if r.get('diseaseId') else None, condition_text=r.get('diseaseFromSource'),
                          curator=r.get('projectId'), lineage_group=f"pmid:{r['literature'][0]}" if r.get('literature') else f"crispr_screen:{r.get('studyId')}", discriminator=r['id'],
                          extra=dict(study_id=r.get('studyId'), overview=r.get('studyOverview'), tail=r.get('statisticalTestTail'), literature=r.get('literature'), datasource=r.get('datasourceId')))

    # ---- gene burden ------------------------------------------------------------------------------------------------------------
    def burden(self, ctx, st):
        rel = ctx.register_release('opentargets', 'opentargets_26_09_evidence_gene_burden', '26.09')
        for f in self.files('evidence_gene_burden'):
            for i, r in enumerate(pq.read_table(f).to_pylist()):
                pv = None
                if r.get('pValueMantissa') is not None and r.get('pValueExponent') is not None:
                    v = r['pValueMantissa'] * 10.0 ** r['pValueExponent']
                    pv = v if 1e-300 < v <= 1 else None
                eff = 'effect' if pv is not None and pv < 0.05 else 'uncertain'
                self.emit(ctx, rel, st, f"burden:{r['id']}", f"evidence_gene_burden/{os.path.basename(f)}#{i}", json.dumps({k: r[k] for k in ('id', 'targetId', 'diseaseId', 'oddsRatio', 'beta', 'pValueMantissa', 'pValueExponent')}, sort_keys=True, default=str),
                          r['targetId'], 'association_replication', 'statistical_association', 'gene_burden_association', r.get('diseaseFromSource') or r.get('diseaseId') or 'trait', eff,
                          odds_ratio=r.get('oddsRatio'), beta=r.get('beta'), ci_low=r.get('oddsRatioConfidenceIntervalLower') if r.get('oddsRatio') is not None else r.get('betaConfidenceIntervalLower'),
                          ci_high=r.get('oddsRatioConfidenceIntervalUpper') if r.get('oddsRatio') is not None else r.get('betaConfidenceIntervalUpper'), p_value=pv, ancestry=r.get('ancestry'),
                          sample_size=r.get('studySampleSize'), n_cases=r.get('studyCases'), inheritance=', '.join(r.get('allelicRequirements') or []) or None, condition_ids=[r['diseaseId']] if r.get('diseaseId') else None,
                          condition_text=r.get('diseaseFromSource'), curator=r.get('projectId'), lineage_group=f"burden_cohort:{r.get('projectId')}:{r.get('cohortId')}", discriminator=r['id'],
                          review_status=r.get('statisticalMethod'), extra=dict(cohort=r.get('cohortId'), release=r.get('releaseVersion'), cases_with_qualifying_variants=r.get('studyCasesWithQualifyingVariants'),
                                                                              statistical_method=r.get('statisticalMethod'), sex=r.get('sex'), direction_on_trait=r.get('directionOnTrait'),
                                                                              scope='burden of qualifying variants in the gene, not an effect of one allele'))
