"""Open Targets Genetics (platform 26.09): fine-mapped credible sets, studies and locus-to-gene predictions.

  credible-set member   level variant, axis fine_mapping_probability. One row per member that reaches the posterior inclusion probability (PIP)
                        threshold of its study type; beta / standard error / p-value are those of the member and are harmonised by the platform to the
                        alternate allele of the variant id (effect_allele = ALT). GWAS credible sets use PIP >= gwas_pip, molecular-QTL sets
                        (eQTL, sQTL, tuQTL, scEQTL, pQTL) use PIP >= qtl_pip, because the QTL sets alone hold ~70M members.
  L2G prediction        level gene, kind computational_prediction: the locus-to-gene model's score for a credible set (rows with score >= l2g_min);
                        the lead variant of the set is linked to the gene in variant_gene.
  study                 entity_study for every study that has an integrated row (GWAS Catalog studies keep the entity the GWAS Catalog adapter wrote).

Lineage: a published GWAS counts once per PubMed id (shared with the GWAS Catalog and ClinPGx rows for the same paper); molecular-QTL sets count once per
project (tissues of one project share donors). Listed but not integrated, and counted in the stats: members below the thresholds, sets with no member
at the threshold, L2G predictions below l2g_min, colocalisation tables (not acquired).
"""
import glob
import json
import os

import pyarrow as pa
import pyarrow.compute as pc
import pyarrow.parquet as pq

from ..normalise import normalise
from .common import DATA, evidence_row

OT = os.path.join(DATA, 'opentargets', '26.09')
SET_COLS = ['studyLocusId', 'studyId', 'variantId', 'chromosome', 'position', 'beta', 'standardError', 'pValueMantissa', 'pValueExponent', 'effectAlleleFrequencyFromSource', 'finemappingMethod',
            'credibleSetIndex', 'credibleSetlog10BF', 'purityMeanR2', 'sampleSize', 'confidence', 'studyType', 'isTransQtl', 'qualityControls', 'locus']
STUDY_COLS = ['studyId', 'geneId', 'projectId', 'studyType', 'traitFromSource', 'traitFromSourceMappedIds', 'pubmedId', 'nSamples', 'nCases', 'cohorts', 'discoverySamples',
              'ldPopulationStructure', 'hasSumstats', 'diseaseIds', 'biosampleId']


def p_value(m, e):
    if m is None or e is None:
        return None
    v = float(m) * (10.0 ** int(e))
    return v if 1e-300 < v <= 1 else None


def split_variant_id(vid):
    parts = (vid or '').split('_')
    return parts if len(parts) == 4 and parts[1].isdigit() else None


def ancestry_text(st):
    d = [f"{x['ancestry']}:{x['sampleSize']}" for x in (st.get('discoverySamples') or ()) if x.get('ancestry')]
    if d:
        return '; '.join(d)
    return '; '.join(f"{x['ldPopulation']}:{x['relativeSampleSize']:.2f}" for x in (st.get('ldPopulationStructure') or ()) if x.get('ldPopulation')) or None


class OpenTargetsGeneticsAdapter:
    source_id = 'opentargets'
    targets = ['opentargets_26_09_credible_set', 'opentargets_26_09_study', 'opentargets_26_09_l2g_prediction', 'opentargets_26_09_biosample']

    def __init__(self, genes, root=None, gwas_pip=0.1, qtl_pip=0.5, l2g_min=0.5):
        self.genes, self.root = genes, root or OT
        self.gwas_pip, self.qtl_pip, self.l2g_min = gwas_pip, qtl_pip, l2g_min

    def prepare(self, ctx):
        ctx.want('study')

    def files(self, name):
        return sorted(glob.glob(os.path.join(self.root, name, '*.parquet')))

    def lineage(self, st):
        if st['studyType'] == 'gwas':
            return f"pmid:{st['pubmedId']}" if st.get('pubmedId') else f"ot_project:{st['projectId']}:{st['studyId']}"
        return f"qtl_project:{st['projectId']}"

    def run(self, ctx, limit=None):
        ctx.adapter = self
        st = ctx.st()
        rel_cs = ctx.register_release('opentargets', 'opentargets_26_09_credible_set', '26.09')
        rel_st = ctx.register_release('opentargets', 'opentargets_26_09_study', '26.09')
        rel_l2g = ctx.register_release('opentargets', 'opentargets_26_09_l2g_prediction', '26.09')
        ctx.register_release('opentargets', 'opentargets_26_09_biosample', '26.09')
        biosample = {}
        for f in self.files('biosample'):
            t = pq.read_table(f, columns=['biosampleId', 'biosampleName'])
            biosample.update(zip(t['biosampleId'].to_pylist(), t['biosampleName'].to_pylist()))
        l2g = {}
        for f in self.files('l2g_prediction'):
            t = pq.read_table(f, columns=['studyLocusId', 'geneId', 'score'])
            t = t.filter(pc.greater_equal(t['score'], self.l2g_min))
            for sl, g, sc in zip(t['studyLocusId'].to_pylist(), t['geneId'].to_pylist(), t['score'].to_pylist()):
                l2g.setdefault(sl, []).append((g, sc))
            st.inc('l2g_rows_at_threshold', t.num_rows)
        study_tables = [pq.read_table(f, columns=STUDY_COLS) for f in self.files('study')]
        studies = pa.concat_tables(study_tables) if study_tables else None
        seen_studies, vg_seen = set(), set()
        for fi, path in enumerate(self.files('credible_set')):
            if limit and st.c.get('credible_sets', 0) >= limit:
                break
            t = pq.read_table(path, columns=SET_COLS)
            fname = os.path.basename(path)
            st.inc('credible_sets', t.num_rows)
            locus = t['locus']
            parent = pc.list_parent_indices(locus)
            flat = pc.list_flatten(locus)
            st.inc('members_total', len(flat))
            pip = pc.struct_field(flat, 'posteriorProbability')
            thr = pc.if_else(pc.equal(pc.take(t['studyType'], parent), 'gwas'), self.gwas_pip, self.qtl_pip)
            idx = pc.indices_nonzero(pc.fill_null(pc.greater_equal(pip, thr), False))
            kept_flat = pc.take(flat, idx).to_pylist()
            kept_parent = pc.take(parent, idx).to_pylist()
            st.inc('members_integrated', len(kept_flat))
            rows = t.drop(['locus']).to_pylist()
            wanted = {rows[p]['studyId'] for p in set(kept_parent)} | {r['studyId'] for r in rows if r['studyLocusId'] in l2g}
            smeta = {}
            if wanted and studies is not None:
                sub = studies.filter(pc.is_in(studies['studyId'], value_set=pa.array(sorted(wanted), pa.string())))
                smeta = {x['studyId']: x for x in sub.to_pylist()}
            norm_cache = {}

            def normalised(vid):
                if vid not in norm_cache:
                    parts = split_variant_id(vid)
                    canon = ctx.ref.contig(parts[0]) if parts else None
                    n = normalise(ctx.ref, canon, int(parts[1]), parts[2], parts[3]) if canon else None
                    norm_cache[vid] = n if n is not None else False
                return norm_cache[vid]

            for p, m in zip(kept_parent, kept_flat):
                row = rows[p]
                s = smeta.get(row['studyId'])
                if s is None:
                    st.inc('members_without_study_record')
                    continue
                self.study_entity(ctx, rel_st, s, biosample, seen_studies, st)
                n = normalised(m['variantId'])
                if n and n.ok:
                    key, mapping, reason = n.key, 'mapped', None
                else:
                    key, mapping, reason = None, ('rejected' if n and n.status == 'ref_mismatch' else 'unmapped'), (n.status if n else 'variant_id_not_parsable')
                is_lead = m['variantId'] == row['variantId']
                pv = p_value(m.get('pValueMantissa'), m.get('pValueExponent'))
                qtl = s['studyType'] != 'gwas'
                bio = biosample.get(s.get('biosampleId')) if qtl else None
                parts = split_variant_id(m['variantId']) or ['', '', '', '']
                condition = (s.get('traitFromSource') or '') if not qtl else f"{s['studyType']} of {s.get('geneId') or ''} in {bio or s.get('biosampleId') or 'unspecified biosample'}"
                raw = json.dumps(dict(studyLocusId=row['studyLocusId'], member=m), sort_keys=True, default=str)
                extra = dict(study_id=row['studyId'], project=s['projectId'], study_type=s['studyType'], credible_set_index=row['credibleSetIndex'], method=row['finemappingMethod'], confidence=row['confidence'],
                             lead_variant=row['variantId'], is_lead=is_lead, is95=m.get('is95CredibleSet'), is99=m.get('is99CredibleSet'), log_bf=m.get('logBF'), r2_to_lead=m.get('r2Overall'),
                             set_log10_bf=row['credibleSetlog10BF'], set_purity_mean_r2=row['purityMeanR2'], quality_controls=row['qualityControls'], trans_qtl=row['isTransQtl'],
                             p_text=None if pv is not None or m.get('pValueMantissa') is None else f"{m['pValueMantissa']}e{m['pValueExponent']}", beta_reference='alternate allele of the variant id')
                ctx.emit_evidence(evidence_row(rel_cs, 'opentargets', row['studyLocusId'], f"credible_set/{fname}#{row['studyLocusId']}:{m['variantId']}", raw, 'variant', key, 'fine_mapping_probability',
                                               'statistical_association', 'qtl_fine_mapping_pip' if qtl else 'gwas_fine_mapping_pip', f"{m['posteriorProbability']:.4g}", 'not_applicable', mapping, reason,
                                               discriminator=m['variantId'], claim_num=m['posteriorProbability'], pip=m['posteriorProbability'], beta=m.get('beta'), se=m.get('standardError'), p_value=pv,
                                               effect_allele=parts[3] or None, other_allele=parts[2] or None, eaf=row['effectAlleleFrequencyFromSource'] if is_lead else None,
                                               sample_size=row['sampleSize'] or s.get('nSamples'), n_cases=s.get('nCases'), ancestry=ancestry_text(s), tissue=bio, condition_ids=(s.get('traitFromSourceMappedIds') or []) + (s.get('diseaseIds') or []) or None,
                                               condition_text=condition or None, gene_ids=[s['geneId']] if qtl and s.get('geneId') else None, origin='germline', review_status=f"{row['finemappingMethod']}; {row['confidence']}",
                                               curator='Open Targets Genetics', lineage_group=self.lineage(s), lineage_note='credible sets of one study / project share samples', unmapped_key=None if key else m['variantId'],
                                               chrom=n.chrom if key else None, pos=n.pos if key else None, extra_json=json.dumps(extra, sort_keys=True, default=str)))
                st.inc(f'integrated_rows:{"qtl" if qtl else "gwas"}:{mapping}')
                if qtl and key and s.get('geneId') and (key, s['geneId']) not in vg_seen:
                    vg_seen.add((key, s['geneId']))
                    ctx.emit('variant_gene', dict(variant_key=key, gene_id=s['geneId'], relation=f"ot_{s['studyType']}_credible_set", consequence=None, release_id=rel_cs), name='opentargets')
            for r in rows:
                if r['studyLocusId'] not in l2g:
                    continue
                s = smeta.get(r['studyId'])
                if s is None:
                    st.inc('l2g_sets_without_study_record')
                    continue
                self.study_entity(ctx, rel_st, s, biosample, seen_studies, st)
                n = normalised(r['variantId'])
                lead_key = n.key if n and n.ok else None
                for gene, score in l2g[r['studyLocusId']]:
                    known = getattr(self.genes, 'by_ensg', None)
                    gid = gene if known is None or gene in known else None
                    condition = (s.get('traitFromSource') or '') if s['studyType'] == 'gwas' else f"{s['studyType']} of {s.get('geneId') or ''}"
                    raw = json.dumps(dict(studyLocusId=r['studyLocusId'], geneId=gene, score=score), sort_keys=True)
                    ctx.emit_evidence(evidence_row(rel_l2g, 'opentargets', r['studyLocusId'], f"l2g_prediction#{r['studyLocusId']}:{gene}", raw, 'gene', gid, 'fine_mapping_probability',
                                                   'computational_prediction', 'l2g_score', f'{score:.3f}', 'not_applicable', 'mapped' if gid else 'unmapped', None if gid else 'gene_not_in_gene_map',
                                                   discriminator=gene, claim_num=score, score=score, score_unit='l2g_score', condition_ids=(s.get('traitFromSourceMappedIds') or []) or None,
                                                   condition_text=condition or None, curator='Open Targets Genetics L2G', lineage_group=self.lineage(s), unmapped_key=None if gid else gene,
                                                   lineage_note='prediction for the credible set; not independent of the fine-mapping rows of the same set',
                                                   extra_json=json.dumps(dict(study_id=r['studyId'], lead_variant=r['variantId'], lead_key=lead_key, study_type=s['studyType']), sort_keys=True)))
                    st.inc('l2g_rows_integrated')
                    if lead_key and gid and (lead_key, gid, 'l2g') not in vg_seen:
                        vg_seen.add((lead_key, gid, 'l2g'))
                        ctx.emit('variant_gene', dict(variant_key=lead_key, gene_id=gid, relation='ot_l2g_lead_variant', consequence=f'{score:.3f}', release_id=rel_l2g), name='opentargets')
        return st.c

    def study_entity(self, ctx, rel_st, s, biosample, seen, st):
        sid = s['studyId']
        if sid in seen:
            return
        seen.add(sid)
        if ctx.found.get(('study', sid)):
            st.inc('studies_described_by_another_source')
            return
        qtl = s['studyType'] != 'gwas'
        ctx.emit('entity_study', dict(study_id=sid, source_id='opentargets', lineage_group=self.lineage(s), study_type=s['studyType'], pubmed=s.get('pubmedId'),
                                      trait_ids=(s.get('traitFromSourceMappedIds') or []) + (s.get('diseaseIds') or []) or None, trait_text=s.get('traitFromSource'), n_total=s.get('nSamples'),
                                      n_cases=s.get('nCases'), ancestry=ancestry_text(s), cohorts=', '.join(s.get('cohorts') or []) or None,
                                      biosample=(biosample.get(s.get('biosampleId')) or s.get('biosampleId')) if qtl else None, has_sumstats=s.get('hasSumstats'), release_id=rel_st), name='opentargets')
        st.inc('studies')
