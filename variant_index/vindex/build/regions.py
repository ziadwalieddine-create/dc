"""Genomic regions that change how a call or a negative result is read: GIAB difficult regions and STRchive repeat-disease loci.

GIAB genome stratifications (v3.6, GRCh38). Each acquired stratification becomes entity_region rows of one kind:
  giab_difficult_<name>   regions where short-read calling is unreliable or ambiguous (low mappability, segmental duplications, tandem repeats, satellites, falsely or
                          collapsed duplicated genes, MHC/KIR/VDJ, problematic promoters, assembly gaps). A call inside one is still reported; it carries the flag
                          'difficult_region' and, for an absent record, 'absent' means less than anywhere else.
  giab_xy_<name>          pseudo-autosomal, X-transposed and ampliconic regions: ploidy and mapping context for sex chromosomes, not a statement about call quality.
Intervals are BED (0-based, half-open) in the source and 1-based inclusive here. The union file is kept next to its components so a flag can say why.

STRchive (disease repeat loci, hg38). One entity_region (kind strchive_locus) per locus and region-level evidence rows, each its own claim:
  strchive_disease_locus   the curated locus -> disease association (inheritance, motifs, structure) as a context row
  repeat_normal_max        largest repeat count regarded as normal (repeat units)
  repeat_pathologic_min    smallest repeat count regarded as pathologic (repeat units)
A locus without a threshold keeps the locus row only: a measured length there is genotyped, not classified. Thresholds are per locus and per repeat unit as
published; alleles between the two are intermediate and are never called pathogenic from the thresholds alone. match_repeat.loci_from_index rebuilds the
locus dictionaries from these rows.
"""
import csv
import gzip
import json
import os

from .common import DATA, evidence_row, retrieved

GIAB = os.path.join(DATA, 'giab', 'v3.6', 'GRCh38')
STR = os.path.join(DATA, 'strchive')
XY = {'chrx_par', 'chry_par', 'chrx_xtr', 'chrx_ampliconic'}


def giab_kind(stem):
    return ('giab_xy_' if stem in XY else 'giab_difficult_') + stem


def num(v):
    try:
        return int(v)
    except (TypeError, ValueError):
        return None


class GiabRegionAdapter:
    source_id = 'giab_stratifications'
    targets = []

    def __init__(self, root=None):
        self.root = root or DATA            # files live at <root>/<dest of the acquisition target>

    def run(self, ctx, limit=None):
        ctx.adapter = self
        st = ctx.st()
        for tid, t in sorted(ctx.targets.items()):
            if t.get('source') != 'giab_stratifications' or t.get('kind') != 'bed.gz':
                continue
            path = os.path.join(self.root, t['dest'])
            if not os.path.exists(path):
                st.inc('stratification_files_not_acquired')
                continue
            stem = os.path.basename(path)[len('GRCh38_'):-len('.bed.gz')].lower()
            kind = giab_kind(stem)
            rel = ctx.register_release('giab_stratifications', tid, 'v3.6')
            with gzip.open(path, 'rt', encoding='utf-8') as f:
                for ln, line in enumerate(f, start=1):
                    if limit and st.c.get(f'intervals:{kind}', 0) >= limit:
                        break
                    c = line.rstrip('\n').split('\t')
                    if len(c) < 3:
                        continue
                    canon = ctx.ref.contig(c[0])
                    s, e = num(c[1]), num(c[2])
                    if canon is None or s is None or e is None or e <= s:
                        st.inc('intervals_on_other_contigs_or_unusable')
                        continue
                    ctx.emit('entity_region', dict(region_id=f'giab:{stem}:{canon}:{s + 1}-{e}', kind=kind, chrom=canon, start=s + 1, end=e, name=stem, gene_ids=None, release_id=rel),
                             name='giab_stratifications')
                    st.inc(f'intervals:{kind}')
        return st.c


class StrchiveAdapter:
    source_id = 'strchive'
    targets = ['strchive_hg38_general', 'strchive_hg38_trgt', 'strchive_hg38_stranger_json']

    def __init__(self, genes=None, root=None):
        self.genes, self.root = genes, root or STR

    def gene_id(self, symbol):
        if not symbol:
            return None
        return (self.genes.from_symbol(symbol) if self.genes else None) or f'symbol:{symbol}'

    def run(self, ctx, limit=None):
        ctx.adapter = self
        st = ctx.st()
        rel_g = ctx.register_release('strchive', 'strchive_hg38_general', retrieved(ctx, 'strchive_hg38_general'))
        rel_t = ctx.register_release('strchive', 'strchive_hg38_trgt', retrieved(ctx, 'strchive_hg38_trgt'))
        rel_j = ctx.register_release('strchive', 'strchive_hg38_stranger_json', retrieved(ctx, 'strchive_hg38_stranger_json'))
        trgt = {}
        with open(os.path.join(self.root, 'STRchive-disease-loci.hg38.TRGT.bed'), encoding='utf-8') as f:
            for ln, line in enumerate(f, start=1):
                c = line.rstrip('\n').split('\t')
                if len(c) >= 4:
                    kv = dict(p.partition('=')[::2] for p in c[3].split(';'))
                    trgt[kv.get('ID')] = dict(motifs=kv.get('MOTIFS'), struc=kv.get('STRUC'), line=ln, raw=line.rstrip('\n'))
        stranger = {}
        with open(os.path.join(self.root, 'STRchive-disease-loci.hg38.stranger.json'), encoding='utf-8') as f:
            for r in json.load(f):
                stranger[r['LocusId']] = r
        seen = set()
        with open(os.path.join(self.root, 'STRchive-disease-loci.hg38.general.bed'), encoding='utf-8', newline='') as f:
            hdr = None
            for ln, line in enumerate(f, start=1):
                line = line.rstrip('\n')
                if line.startswith('#'):
                    hdr = line.lstrip('#').split('\t')
                    continue
                if not line.strip():
                    continue
                if limit and st.c.get('loci', 0) >= limit:
                    break
                r = dict(zip(hdr, line.split('\t')))
                lid = r['id']
                seen.add(lid)
                st.inc('loci')
                canon = ctx.ref.contig(r['chrom'])
                gene = self.gene_id(r['gene'])
                j, tg = stranger.get(lid), trgt.get(lid)
                rid = f'strchive:{lid}'
                if canon is None:
                    st.inc('loci_on_other_contigs')
                    continue
                ctx.emit('entity_region', dict(region_id=rid, kind='strchive_locus', chrom=canon, start=int(r['start']) + 1, end=int(r['stop']), name=lid, gene_ids=[gene] if gene else None, release_id=rel_g),
                         name='strchive')
                lineage = dict(lineage_group=f'strchive:{lid}', lineage_note='one curated STRchive locus; its rows describe one locus, they are not independent studies')
                common = dict(chrom=canon, pos=int(r['start']) + 1, end=int(r['stop']), gene_ids=[gene] if gene else None, curator='STRchive', inheritance=r['inheritance'] or None, **lineage)
                extra = dict(reference_motif=r['reference_motif_reference_orientation'], pathogenic_motif=r['pathogenic_motif_reference_orientation'], bed_pathogenic_min=r['pathogenic_min'],
                             locus_structure=j and j['LocusStructure'], display_ru=j and j['DisplayRU'], variant_type=j and j['VariantType'], reference_region=j and j['ReferenceRegion'],
                             pathologic_region=j and j.get('PathologicRegion'), trgt_motifs=tg and tg['motifs'], trgt_structure=tg and tg['struc'])
                ctx.emit_evidence(evidence_row(rel_g, 'strchive', lid, f'STRchive-disease-loci.hg38.general.bed#L{ln}', line, 'region', rid, 'clinical_validity_actionability',
                                               'curated_expert_classification', 'strchive_disease_locus', r['disease'], 'effect', condition_text=r['disease'],
                                               extra_json=json.dumps({k: v for k, v in extra.items() if v not in (None, '')}, sort_keys=True), **common))
                if j is None:
                    st.inc('loci_without_stranger_record')
                    continue
                raw = json.dumps(j, sort_keys=True, ensure_ascii=False)
                loc = f"STRchive-disease-loci.hg38.stranger.json#{lid}"
                for ctype, key, assertion in (('repeat_normal_max', 'NormalMax', 'no_effect'), ('repeat_pathologic_min', 'PathologicMin', 'effect')):
                    v = j.get(key)
                    if v is None:
                        st.inc(f'{ctype}_missing')
                        continue
                    ctx.emit_evidence(evidence_row(rel_j, 'strchive', lid, loc, raw, 'region', rid, 'effect_magnitude_direction', 'curated_expert_classification', ctype, str(v), assertion,
                                                   score=float(v), claim_num=float(v), score_unit='repeat_units_of_' + (j['DisplayRU'] or 'unit'),
                                                   extra_json=json.dumps(dict(display_ru=j['DisplayRU'], locus_structure=j['LocusStructure'], pathologic_region=j.get('PathologicRegion'),
                                                                              note='threshold as published for this locus; between the two thresholds an allele is intermediate'), sort_keys=True),
                                                   **common))
                    st.inc(ctype)
        st.inc('stranger_records_without_general_row', len(set(stranger) - seen))
        return st.c
