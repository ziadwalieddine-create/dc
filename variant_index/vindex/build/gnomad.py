"""gnomAD v4.1 -> evidence rows.

  constraint (v4.1.2)      level gene, axis mechanistic_confidence, kind population_observation: one row for loss-of-function constraint (LOEUF = upper bound
                           of the observed/expected 90% interval, pLI, counts) and one for missense constraint, on the MANE Select transcript of each gene.
                           Constraint is a property of the gene's population variation; it says nothing about a particular allele.
  LoF curation             level variant, axis mechanistic_confidence, kind curated_expert_classification: the gnomAD curators' verdict (LoF, likely LoF,
                           uncertain, likely not LoF, not LoF) on predicted loss-of-function variants in genes with incomplete penetrance, with the
                           curation flags that were raised.
  joint frequency slice    level variant, axis applicability_context: allele count, number, frequency and per-ancestry frequencies of the joint exomes + genomes
                           release at the positions covered by the MaveDB score sets (the only gnomAD frequency slice acquired; whole-genome gnomAD is not mirrored).
                           A filtered site keeps its filter status; absence from the slice means 'outside the acquired regions', never 'not observed'.
  predictions (exomes)     level variant, kind computational_prediction: CADD, REVEL, SpliceAI, Pangolin, SIFT, PolyPhen scores gnomAD annotates onto the same
                           sites. They are predictions and are marked as such; none is used as an empirical result.
"""
import csv
import gzip
import json
import os

from ..normalise import normalise
from .common import DATA, evidence_row, retrieved

GN = os.path.join(DATA, 'gnomad')
RG = os.path.join(DATA, 'gnomad_regional')
VERDICT_EFFECT = {'LoF': 'effect', 'Likely LoF': 'effect', 'Uncertain': 'uncertain', 'Likely not LoF': 'no_effect', 'Not LoF': 'no_effect'}
ANC = ['afr', 'ami', 'amr', 'asj', 'eas', 'fin', 'mid', 'nfe', 'remaining', 'sas']


def num(s):
    try:
        v = float(s)
    except (TypeError, ValueError):
        return None
    return v if v == v else None


def i(s):
    v = num(s)
    return int(v) if v is not None else None


class GnomadConstraintAdapter:
    source_id = 'gnomad_constraint'
    targets = ['gnomad_constraint_v4_1_2']

    def __init__(self, genes, path=None):
        self.genes = genes
        self.path = path or os.path.join(GN, 'gnomad.v4.1.2.constraint_metrics.tsv.bgz')

    def run(self, ctx, limit=None):
        ctx.adapter = self
        st = ctx.st()
        rel = ctx.register_release('gnomad', 'gnomad_constraint_v4_1_2', 'v4.1.2')
        known = getattr(self.genes, 'by_ensg', None)
        rows = []
        with gzip.open(self.path, 'rt', encoding='utf-8', newline='') as f:
            hdr = f.readline().rstrip('\n').split('\t')
            ix = {n: k for k, n in enumerate(hdr)}
            for ln, line in enumerate(f, start=2):
                c = line.rstrip('\n').split('\t')
                st.inc('transcript_rows')
                if c[ix['mane_select']] != 'true':
                    continue
                if limit and len(rows) >= limit:
                    break
                rows.append((ln, line, c))
        # the file lists each MANE Select transcript twice: the Ensembl one (ENSG gene id) and its RefSeq partner (NCBI gene id); one gene is one measurement
        ensembl = {c[ix['gene_id']] for _ln, _l, c in rows if c[ix['gene_id']].startswith('ENSG')}
        for ln, line, c in rows:
            st.inc('mane_rows')
            g = lambda name: c[ix[name]] if name in ix and ix[name] < len(c) else ''
            gid = g('gene_id')
            if gid.isdigit():
                mapped = self.genes.from_ncbi(gid) if self.genes else None
                if mapped in ensembl:
                    st.inc('refseq_mane_rows_repeating_an_ensembl_row')
                    continue
                gid = mapped or f'NCBIGene:{gid}'
                st.inc('refseq_mane_rows_without_an_ensembl_partner')
            subject, status, reason = (gid, 'mapped', None) if gid and (known is None or gid in known) else (None, 'unmapped', 'gene_not_in_gene_map')
            common = dict(gene_ids=[gid] if gid else None, curator='gnomAD v4.1.2 constraint', lineage_group='gnomad_v4.1_constraint',
                          lineage_note='one model fitted to one release; the loss-of-function and missense rows share it', unmapped_key=None if subject else gid or g('gene'))
            flags = [x for x in (g('gene_flags') + g('constraint_flags')).replace('[', ' ').replace(']', ' ').replace('"', ' ').replace(',', ' ').split() if x]
            quality = {k.split('.', 1)[1]: num(g(k)) for k in ('gene_quality_metrics.exome_prop_bp_AN90', 'gene_quality_metrics.exome_mean_AS_MQ', 'gene_quality_metrics.exome_prop_segdup', 'gene_quality_metrics.exome_prop_LCR')}
            for kind, prefix, label in (('gnomad_lof_constraint', 'lof', 'LOEUF'), ('gnomad_missense_constraint', 'mis', 'upper oe')):
                upper, oe = num(g(f'{prefix}.oe_ci.upper')), num(g(f'{prefix}.oe'))
                if upper is None and oe is None:
                    st.inc(f'{kind}:no_metric')
                    continue
                extra = dict(transcript=g('transcript'), transcript_version=g('transcript_version'), symbol=g('gene'), canonical=g('canonical'), cds_length=i(g('cds_length')), coding_exons=i(g('num_coding_exons')),
                             observed=num(g(f'{prefix}.obs')), expected=num(g(f'{prefix}.exp')), oe=oe, oe_lower=num(g(f'{prefix}.oe_ci.lower')), oe_upper=upper, z_score=num(g(f'{prefix}.z_score')),
                             synonymous_z_score=num(g('syn.z_score')), flags=flags, gene_quality=quality)
                if prefix == 'lof':
                    extra.update(pLI=num(g('lof.pLI')), pNull=num(g('lof.pNull')), pRec=num(g('lof.pRec')), upper_decile=i(g('lof.oe_ci.upper_bin_decile')), upper_rank=i(g('lof.oe_ci.upper_rank')))
                ctx.emit_evidence(evidence_row(rel, 'gnomad_constraint', f"constraint:{gid}", f'gnomad.v4.1.2.constraint_metrics.tsv.bgz#L{ln}', line.rstrip('\n'), 'gene', subject, 'mechanistic_confidence', 'population_observation', kind,
                                               f'{label} {upper:.3f}' if upper is not None else f'oe {oe:.3f}', 'not_applicable', status, reason, discriminator=prefix, claim_num=upper, score=oe,
                                               score_unit='observed/expected', origin='germline', extra_json=json.dumps(extra, sort_keys=True), **common))
                st.inc(f'{kind}:{status}')
        return st.c


class GnomadLofCurationAdapter:
    source_id = 'gnomad'
    targets = ['gnomad_v4_1_lof_curation_incomplete_penetrance']

    def __init__(self, genes=None, path=None):
        self.genes = genes
        self.path = path or os.path.join(GN, 'lof_curation', 'incomplete_penetrance_curation_results.csv')

    def run(self, ctx, limit=None):
        ctx.adapter = self
        st = ctx.st()
        rel = ctx.register_release('gnomad', 'gnomad_v4_1_lof_curation_incomplete_penetrance', retrieved(ctx, 'gnomad_v4_1_lof_curation_incomplete_penetrance'))
        with open(self.path, encoding='utf-8', newline='') as f:
            text = f.read()
        lines = text.splitlines()
        for ln, row in enumerate(csv.DictReader(lines), start=2):
            st.inc('records')
            parts = row['Variant ID'].split('-')
            canon = ctx.ref.contig(parts[0]) if len(parts) == 4 else None
            n = normalise(ctx.ref, canon, int(parts[1]), parts[2], parts[3]) if canon and parts[1].isdigit() else None
            if n is not None and n.ok:
                key, mapping, reason = n.key, 'mapped', None
            else:
                key, mapping, reason = None, 'rejected' if n is not None and n.status == 'ref_mismatch' else 'unmapped', n.status if n is not None else 'variant_id_not_parsable'
            gid, _, sym = row['Gene'].partition(':')
            flags = sorted(k[5:] for k, v in row.items() if k.startswith('Flag ') and v == 'TRUE') + (['LoF in untranslated transcript'] if row.get('LoF in Untranslated Transcript') == 'TRUE' else [])
            verdict = row['Verdict']
            ctx.emit_evidence(evidence_row(rel, 'gnomad', f"lofcuration:{row['Variant ID']}:{gid}", f'incomplete_penetrance_curation_results.csv#L{ln}', lines[ln - 1], 'variant', key, 'mechanistic_confidence',
                                           'curated_expert_classification', 'gnomad_lof_curation_verdict', verdict, VERDICT_EFFECT.get(verdict, 'uncertain'), mapping, reason, discriminator=gid, gene_ids=[gid] if gid else None,
                                           curator='gnomAD LoF curation (incomplete penetrance genes)', lineage_group=f"gnomad_lof_curation:{row['Variant ID']}", chrom=n.chrom if key else None, pos=n.pos if key else None,
                                           unmapped_key=None if key else row['Variant ID'], origin='germline', extra_json=json.dumps(dict(gene=sym, flags=flags), sort_keys=True)))
            st.inc(f'curated_variants:{mapping}')
            if key and gid and (known := getattr(self.genes, 'by_ensg', None)) is not None and gid in known:
                ctx.emit('variant_gene', dict(variant_key=key, gene_id=gid, relation='gnomad_lof_curation_gene', consequence=None, release_id=rel), name='gnomad')
        return st.c


class GnomadFrequencyAdapter:
    source_id = 'gnomad'
    targets = ['gnomad_v4_1_joint_regional', 'gnomad_v4_1_exomes_regional']

    def __init__(self, joint=None, exomes=None):
        self.joint = joint or os.path.join(RG, 'gnomad_v4_1_joint_at_mavedb_positions.tsv.gz')
        self.exomes = exomes or os.path.join(RG, 'gnomad_v4_1_exomes_at_mavedb_positions.tsv.gz')

    def key(self, ctx, c):
        canon = ctx.ref.contig(c['chrom'])
        n = normalise(ctx.ref, canon, int(c['pos']), c['ref'], c['alt']) if canon and c['pos'].isdigit() else None
        if n is not None and n.ok:
            return n.key, 'mapped', None, n
        return None, ('rejected' if n is not None and n.status == 'ref_mismatch' else 'unmapped'), (n.status if n is not None else 'unusable_coordinates'), n

    def run(self, ctx, limit=None):
        ctx.adapter = self
        st = ctx.st()
        self.frequencies(ctx, st, limit)
        self.predictions(ctx, st, limit)
        return st.c

    def frequencies(self, ctx, st, limit):
        rel = ctx.register_release('gnomad', 'gnomad_v4_1_joint_regional', retrieved(ctx, 'gnomad_v4_1_joint_regional'))
        with gzip.open(self.joint, 'rt', encoding='utf-8', newline='') as f:
            hdr = f.readline().rstrip('\n').split('\t')
            for ln, line in enumerate(f, start=2):
                if limit and st.c.get('joint_rows', 0) >= limit:
                    break
                st.inc('joint_rows')
                c = dict(zip(hdr, line.rstrip('\n').split('\t')))
                key, mapping, reason, n = self.key(ctx, c)
                af = num(c.get('AF_joint'))
                by = {a: dict(ac=i(c.get(f'AC_joint_{a}')), an=i(c.get(f'AN_joint_{a}')), af=num(c.get(f'AF_joint_{a}'))) for a in ANC if c.get(f'AN_joint_{a}')}
                extra = dict(filter=c.get('filter'), ac=i(c.get('AC_joint')), an=i(c.get('AN_joint')), nhomalt=i(c.get('nhomalt_joint')), grpmax=c.get('grpmax_joint') or None, af_grpmax=num(c.get('AF_grpmax_joint')),
                             faf95_max=num(c.get('fafmax_faf95_max_joint')), faf95_max_group=c.get('fafmax_faf95_max_gen_anc_joint') or None, exomes=dict(ac=i(c.get('AC_exomes')), an=i(c.get('AN_exomes'))),
                             genomes=dict(ac=i(c.get('AC_genomes')), an=i(c.get('AN_genomes'))), exomes_filters=c.get('exomes_filters') or None, genomes_filters=c.get('genomes_filters') or None,
                             not_called_in_exomes=bool(c.get('not_called_in_exomes')), not_called_in_genomes=bool(c.get('not_called_in_genomes')), by_ancestry=by)
                ctx.emit_evidence(evidence_row(rel, 'gnomad', f"joint:{c['chrom']}-{c['pos']}-{c['ref']}-{c['alt']}", f'gnomad_v4_1_joint_at_mavedb_positions.tsv.gz#L{ln}', line.rstrip('\n'), 'variant', key, 'applicability_context',
                                               'population_observation', 'gnomad_allele_frequency', f'AF {af:.3g}' if af is not None else 'AF not stated', 'not_applicable', mapping, reason,
                                               claim_num=af, eaf=af, score=af, score_unit='allele_frequency', ancestry=c.get('grpmax_joint') or None, sample_size=(i(c.get('AN_joint')) // 2) if i(c.get('AN_joint')) else None,
                                               lineage_group='gnomad_v4.1_joint', lineage_note='one release; exome and genome components are not independent of the joint row', curator='gnomAD v4.1 joint',
                                               chrom=n.chrom if key else None, pos=n.pos if key else None, unmapped_key=None if key else f"{c['chrom']}:{c['pos']}:{c['ref']}:{c['alt']}", origin='germline',
                                               review_status=c.get('filter'), extra_json=json.dumps({k: v for k, v in extra.items() if v not in (None, {}, False)}, sort_keys=True)))
                st.inc(f'frequency_rows:{mapping}')

    PRED = [('cadd_phred', 'CADD', 'phred', None), ('revel_max', 'REVEL', 'score', None), ('spliceai_ds_max', 'SpliceAI', 'max delta score', None), ('pangolin_largest_ds', 'Pangolin', 'largest delta score', None),
            ('sift_max', 'SIFT', 'score', None), ('polyphen_max', 'PolyPhen-2', 'score', None)]

    def predictions(self, ctx, st, limit):
        rel = ctx.register_release('gnomad', 'gnomad_v4_1_exomes_regional', retrieved(ctx, 'gnomad_v4_1_exomes_regional'))
        with gzip.open(self.exomes, 'rt', encoding='utf-8', newline='') as f:
            hdr = f.readline().rstrip('\n').split('\t')
            for ln, line in enumerate(f, start=2):
                if limit and st.c.get('exome_rows', 0) >= limit:
                    break
                st.inc('exome_rows')
                c = dict(zip(hdr, line.rstrip('\n').split('\t')))
                key, mapping, reason, n = self.key(ctx, c)
                for col, name, unit, _ in self.PRED:
                    v = num(c.get(col))
                    if v is None:
                        continue
                    ctx.emit_evidence(evidence_row(rel, 'gnomad', f"exomes:{c['chrom']}-{c['pos']}-{c['ref']}-{c['alt']}", f'gnomad_v4_1_exomes_at_mavedb_positions.tsv.gz#L{ln}', line.rstrip('\n'), 'variant', key,
                                                   'mechanistic_confidence', 'computational_prediction', f'prediction:{name}', f'{name} {v:g}', 'not_applicable', mapping, reason, discriminator=name, claim_num=v,
                                                   score=v, score_unit=f'{name} {unit}', lineage_group=f'predictor:{name}', lineage_note='a model output annotated onto the site by gnomAD; a prediction, not an observation',
                                                   curator=name, chrom=n.chrom if key else None, pos=n.pos if key else None, unmapped_key=None if key else f"{c['chrom']}:{c['pos']}:{c['ref']}:{c['alt']}",
                                                   origin='germline'))
                    st.inc(f'prediction_rows:{name}:{mapping}')
