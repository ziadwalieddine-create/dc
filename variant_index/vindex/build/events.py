"""Structural-variant and copy-number events observed in population resources (GRCh38) -> entity_event rows.

An event is the unit a structural call is matched against (match_sv.compare): type, breakpoints with their confidence intervals, length, copy number,
precision and the frequency with which it was seen. Frequencies live on the event, not in evidence rows; per-ancestry detail and the source's own
annotations stay in extra_json. Intervals are 1-based, inclusive.

  gnomAD SV v4.1     short-read genome SVs (depth, paired-end, split-read evidence): DEL DUP INS INV CPX CTX BND and multiallelic CNVs. Calls from
                     depth alone are marked imprecise. Frequencies are for biallelic sites; multiallelic CNVs carry the non-reference copy-state
                     frequency.
  gnomAD CNV v4.1    exome-capture copy-number variants (DEL DUP) with the observed range of breakpoints (POSMIN/POSMAX, ENDMIN/ENDMAX -> confidence
                     intervals); always imprecise.
  DGV GRCh38         the merged variant records of the Database of Genomic Variants (array, sequencing and mixed studies): type from the source's
                     own subtype; gains and losses without a stated mechanism stay CNV_GAIN / CNV_LOSS; imprecise.
  HGSVC3 v1.0        haplotype-resolved long-read assemblies (insertions, deletions, inversions): precise breakpoints, allele frequency over the
                     phased haplotypes with a called allele.
Events from different resources are not merged: two resources describing the same event give two events, and the matcher reports whichever are
equivalent to a call. Counting both as independent replication would be wrong, so event frequencies are context, not evidence rows.
"""
import csv
import gzip
import json
import os
import re

from .common import DATA, retrieved

ANCESTRIES = ['afr', 'ami', 'amr', 'asj', 'eas', 'fin', 'mid', 'nfe', 'sas', 'remaining']


def info_dict(text):
    out = {}
    for kv in text.split(';'):
        if not kv:
            continue
        k, eq, v = kv.partition('=')
        out[k] = v if eq else True
    return out


def to_int(v):
    try:
        return int(float(v))
    except (TypeError, ValueError):
        return None


def to_float(v):
    try:
        f = float(v)
    except (TypeError, ValueError):
        return None
    return f if f == f else None


class _EventAdapter:
    source_id = 'events'
    targets = []

    def event(self, ctx, st, rel, **kw):
        row = dict(source_id=self.source_id, release_id=rel, origin='germline_population', assembly='GRCh38')
        row.update(kw)
        if not row.get('chrom'):
            st.inc('events_on_unplaced_or_unknown_contigs')
            return
        if row['end'] is None or row['start'] is None or row['end'] < row['start']:
            st.inc('events_with_unusable_interval')
            return
        ctx.emit('entity_event', row, name=self.source_id)
        st.inc(f"events:{row['event_type']}")


class GnomadSvAdapter(_EventAdapter):
    source_id = 'gnomad_sv'
    targets = ['gnomad_sv_v4_1_sites_vcf']

    def __init__(self, path=None):
        self.path = path or os.path.join(DATA, 'gnomad', 'gnomad.v4.1.sv.sites.vcf.gz')

    def run(self, ctx, limit=None):
        ctx.adapter = self
        st = ctx.st()
        rel = ctx.register_release('gnomad', 'gnomad_sv_v4_1_sites_vcf', 'v4.1')
        with gzip.open(self.path, 'rt', encoding='utf-8') as f:
            for line in f:
                if line.startswith('#'):
                    continue
                if limit and st.c.get('records', 0) >= limit:
                    break
                st.inc('records')
                c = line.rstrip('\n').split('\t')
                chrom = ctx.ref.contig(c[0])
                info = info_dict(c[7])
                svtype = info.get('SVTYPE') or c[4].strip('<>')
                pos = int(c[1])
                end = to_int(info.get('END')) or pos
                if svtype in ('BND', 'CTX') or info.get('CHR2') not in (None, c[0]):
                    end = pos                                    # a breakpoint on this contig only; the partner is in extra_json
                algos = set((info.get('ALGORITHMS') or '').split(',')) - {''}
                precise = bool(algos - {'depth'}) and svtype not in ('CNV',)
                multi = 'MULTIALLELIC' in info
                af = to_float(info.get('AF')) if not multi else to_float(info.get('CN_NONREF_FREQ'))
                svlen = to_int(info.get('SVLEN'))
                anc = {a: to_float(info.get(f'AF_{a}')) for a in ANCESTRIES if f'AF_{a}' in info}
                extra = dict(filter=c[6], algorithms=sorted(algos), evidence=info.get('EVIDENCE'), multiallelic=multi, af_by_ancestry=anc, ac=to_int(info.get('AC')), an=to_int(info.get('AN')),
                             n_homalt=to_int(info.get('N_HOMALT')), chr2=info.get('CHR2'), end2=to_int(info.get('END2')), cpx_type=info.get('CPX_TYPE'), unresolved_type=info.get('UNRESOLVED_TYPE'),
                             predicted_lof=info.get('PREDICTED_LOF'), predicted_copy_gain=info.get('PREDICTED_COPY_GAIN'), predicted_intragenic_exon_dup=info.get('PREDICTED_INTRAGENIC_EXON_DUP'),
                             predicted_partial_exon_dup=info.get('PREDICTED_PARTIAL_EXON_DUP'), predicted_promoter=info.get('PREDICTED_PROMOTER'), predicted_intergenic=bool(info.get('PREDICTED_INTERGENIC')),
                             bothsides_support=bool(info.get('BOTHSIDES_SUPPORT')), low_confidence_repetitive_large_dup=bool(info.get('LOW_CONFIDENCE_REPETITIVE_LARGE_DUP')),
                             outlier_sample_enriched=bool(info.get('OUTLIER_SAMPLE_ENRICHED_LENIENT')), cn_freq=info.get('CN_FREQ'), cn_status=info.get('CN_STATUS'))
                extra = {k: v for k, v in extra.items() if v not in (None, False, [], {}, '')}
                self.event(ctx, st, rel, event_id=f'gnomad_sv:{c[2]}', accession=c[2], event_type='CNV' if multi and svtype == 'CNV' else svtype, chrom=chrom, start=pos, end=end,
                           start_ci_lo=None, start_ci_hi=None, end_ci_lo=None, end_ci_hi=None, svlen=abs(svlen) if svlen and abs(svlen) > 1 else None,
                           copy_number=None, precise=precise, population_af=af, n_samples=to_int(info.get('N_BI_GENOS')) or to_int(info.get('CN_NUMBER')), raw_type=svtype,
                           extra_json=json.dumps(extra, sort_keys=True))
        return st.c


class GnomadCnvAdapter(_EventAdapter):
    source_id = 'gnomad_cnv'
    targets = ['gnomad_exome_cnv_v4_1']

    def __init__(self, path=None):
        self.path = path or os.path.join(DATA, 'gnomad', 'gnomad.v4.1.cnv.all.vcf.gz')

    def run(self, ctx, limit=None):
        ctx.adapter = self
        st = ctx.st()
        rel = ctx.register_release('gnomad', 'gnomad_exome_cnv_v4_1', 'v4.1')
        with gzip.open(self.path, 'rt', encoding='utf-8') as f:
            for line in f:
                if line.startswith('#'):
                    continue
                if limit and st.c.get('records', 0) >= limit:
                    break
                st.inc('records')
                c = line.rstrip('\n').split('\t')
                info = info_dict(c[7])
                pos, end = int(c[1]), to_int(info.get('END'))
                pmin, pmax, emin, emax = (to_int(info.get(k)) for k in ('POSMIN', 'POSMAX', 'ENDMIN', 'ENDMAX'))
                anc = {a: to_float(info.get(f'SF_{a}')) for a in ANCESTRIES if f'SF_{a}' in info}
                extra = dict(filter=c[6], genes=info.get('Genes'), median_exons_spanned=to_int(info.get('N_EXN_VAR')), sample_count=to_float(info.get('SC')), sample_total=to_float(info.get('SN')),
                             sf_by_ancestry=anc, note='exome-capture calls: breakpoints are the observed ranges, not resolved positions')
                self.event(ctx, st, rel, event_id=f'gnomad_cnv:{c[2]}', accession=c[2], event_type=info.get('SVTYPE') or c[4].strip('<>'), chrom=ctx.ref.contig(c[0]), start=pos, end=end,
                           start_ci_lo=(pmin - pos) if pmin is not None else None, start_ci_hi=(pmax - pos) if pmax is not None else None, end_ci_lo=(emin - end) if emin is not None and end else None,
                           end_ci_hi=(emax - end) if emax is not None and end else None, svlen=to_int(info.get('SVLEN')), copy_number=None, precise=False, population_af=to_float(info.get('SF')),
                           n_samples=to_int(info.get('SN')), raw_type=info.get('SVTYPE'), extra_json=json.dumps({k: v for k, v in extra.items() if v not in (None, {}, '')}, sort_keys=True))
        return st.c


DGV_TYPES = {'deletion': 'DEL', 'duplication': 'DUP', 'tandem duplication': 'DUP', 'gain': 'CNV_GAIN', 'loss': 'CNV_LOSS', 'insertion': 'INS', 'mobile element insertion': 'INS',
             'novel sequence insertion': 'INS', 'inversion': 'INV', 'gain+loss': 'CNV', 'complex': 'CPX', 'sequence alteration': 'CPX', 'insertion/deletion': 'CPX', 'deletion/duplication': 'CNV',
             'tandem duplication ': 'DUP', 'interchromosomal translocation': 'CTX', 'intrachromosomal translocation': 'CPX', 'translocation': 'CTX', 'herv insertion': 'INS', 'alu insertion': 'INS',
             'line1 insertion': 'INS', 'sva insertion': 'INS', 'alu deletion': 'DEL', 'line1 deletion': 'DEL', 'sva deletion': 'DEL', 'mobile element deletion': 'DEL', 'herv deletion': 'DEL'}


class DgvAdapter(_EventAdapter):
    source_id = 'dgv'
    targets = ['dgv_grch38_variants_2025_12_01']

    def __init__(self, path=None):
        self.path = path or os.path.join(DATA, 'dgv', 'GRCh38_hg38_variants_2025-12-01.txt')

    def run(self, ctx, limit=None):
        ctx.adapter = self
        st = ctx.st()
        rel = ctx.register_release('dgv', 'dgv_grch38_variants_2025_12_01', retrieved(ctx, 'dgv_grch38_variants_2025_12_01'))
        with open(self.path, encoding='utf-8', errors='replace', newline='') as f:
            hdr = f.readline().rstrip('\n').split('\t')
            for line in f:
                if limit and st.c.get('records', 0) >= limit:
                    break
                st.inc('records')
                r = dict(zip(hdr, line.rstrip('\n').split('\t')))
                if r.get('mergedorsample') != 'M':
                    st.inc('records_not_merged_variants')
                    continue
                sub = (r.get('variantsubtype') or '').strip().lower()
                etype = DGV_TYPES.get(sub)
                if etype is None:
                    st.inc('records_with_unmapped_subtype')
                    etype = 'CNV' if r.get('varianttype') == 'CNV' else (sub.upper() or 'UNKNOWN')
                gains, losses, n = to_int(r.get('observedgains')), to_int(r.get('observedlosses')), to_int(r.get('samplesize'))
                carriers = ((gains or 0) + (losses or 0)) / n if n and (gains is not None or losses is not None) else None
                extra = dict(study=r.get('reference'), pubmed=r.get('pubmedid'), method=r.get('method'), platform=r.get('platform'), observed_gains=gains, observed_losses=losses,
                             sample_size=n, carrier_fraction=carriers, genes=r.get('genes'), supporting_variants=(r.get('supportingvariants') or '')[:200] or None, subtype=sub,
                             note='DGV states no allele frequency for merged variants (the frequency column is empty); carrier_fraction is observed gains + losses over the study sample size')
                self.event(ctx, st, rel, event_id=f"dgv:{r['variantaccession']}", accession=r['variantaccession'], event_type=etype, chrom=ctx.ref.contig(r['chr']), start=to_int(r.get('start')),
                           end=to_int(r.get('end')), start_ci_lo=None, start_ci_hi=None, end_ci_lo=None, end_ci_hi=None, svlen=None, copy_number=None, precise=False, population_af=None, n_samples=n,
                           raw_type=r.get('variantsubtype'), extra_json=json.dumps({k: v for k, v in extra.items() if v not in (None, '')}, sort_keys=True))
        return st.c


class Hgsvc3Adapter(_EventAdapter):
    source_id = 'hgsvc3'
    targets = ['hgsvc3_grch38_sv_insdel_sym_vcf', 'hgsvc3_grch38_sv_inv_sym_vcf']

    def __init__(self, paths=None):
        d = os.path.join(DATA, 'hgsvc3', '1.0')
        self.paths = paths or [('hgsvc3_grch38_sv_insdel_sym_vcf', os.path.join(d, 'variants_GRCh38_sv_insdel_sym_HGSVC2024v1.0.vcf.gz')),
                               ('hgsvc3_grch38_sv_inv_sym_vcf', os.path.join(d, 'variants_GRCh38_sv_inv_sym_HGSVC2024v1.0.vcf.gz'))]

    def run(self, ctx, limit=None):
        ctx.adapter = self
        st = ctx.st()
        seen = set()
        for target, path in self.paths:
            rel = ctx.register_release('hgsvc3', target, '1.0')
            with gzip.open(path, 'rt', encoding='utf-8') as f:
                for line in f:
                    if line.startswith('##'):
                        continue
                    if line.startswith('#'):
                        n_samples = len(line.rstrip('\n').split('\t')) - 9
                        continue
                    if limit and st.c.get('records', 0) >= limit:
                        break
                    st.inc('records')
                    c = line.rstrip('\n').split('\t')
                    info = info_dict(c[7])
                    svtype = info.get('SVTYPE') or c[4].strip('<>')
                    pos = int(c[1])
                    svlen = to_int(info.get('SVLEN'))
                    end = to_int(info.get('END'))
                    if end is None:
                        end = pos + (abs(svlen) if svlen and svtype != 'INS' else 0)
                    ac = an = 0
                    for gt in c[9:]:
                        for a in re.split(r'[|/]', gt.split(':')[0]):
                            if a in ('0', '1'):
                                an += 1
                                ac += a == '1'
                    eid = f'hgsvc3:{c[2]}'
                    if eid in seen:
                        eid = f'hgsvc3:{c[2]}#{c[0]}:{pos}:{svtype}'
                    seen.add(eid)
                    extra = dict(ac=ac, an=an, haplotype_columns=n_samples, source_ids=info.get('ID'), note='allele frequency over the phased haplotypes with a called allele; assembly-based, breakpoints precise')
                    self.event(ctx, st, rel, event_id=eid, accession=c[2], event_type=svtype, chrom=ctx.ref.contig(c[0]), start=pos, end=end, start_ci_lo=None, start_ci_hi=None, end_ci_lo=None,
                               end_ci_hi=None, svlen=abs(svlen) if svlen else None, copy_number=None, precise=True, population_af=(ac / an) if an else None, n_samples=an // 2 if an else None,
                               raw_type=svtype, extra_json=json.dumps({k: v for k, v in extra.items() if v not in (None, '')}, sort_keys=True))
        return st.c
