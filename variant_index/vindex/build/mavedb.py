"""MaveDB multiplexed assays of variant effect -> assay entities and variant-level evidence.

Genomic coordinates come from Ensembl's VEP MaveDB plugin file (GRCh38; one row per candidate nucleotide allele of an assayed change, with the
assay's score); assay metadata, supersession, publications and score calibrations come from the MaveDB public dump (main.json). Every candidate
allele is checked against the reference sequence (REF mismatches are kept as rejected rows, never repaired).

Rows
  score row        axis empirical_functional_support, kind direct_experimental. claim_type 'mavedb_variant_score' when the assay resolved the
                   nucleotide allele (the author gave a nucleotide HGVS) and 'mavedb_protein_level_score' when it measured an amino-acid change:
                   then the subject is a nucleotide allele that encodes the change (derived from the genetic code by Ensembl), the evidence is about
                   the protein change, and effects outside the coding sequence (splicing, regulatory) were not tested.
  calibration row  one per calibration whose range contains the score: functional class (normal / abnormal / not_specified) and the ACMG/AMP
                   criterion and strength the calibration assigns. Same lineage group as the score row (an interpretation of it, not new evidence);
                   calibrations marked research-use-only stay marked.
Second source of genomic placement: MaveDB's own mapping (annotations.csv, mavedb.post_mapped_hgvs_g) joined to the score tables extracted from the dump
(acquisition/extract_mavedb_scores.py). These are single-variant, nucleotide-level assays on GRCh38 and are always 'mavedb_variant_score' rows; the Ensembl
plugin file holds only protein-level assays and multi-variant constructs.
Not integrated, counted: constructs with several changes in cis (the score belongs to the combination, not to one component), score sets that
MaveDB marks as superseded, rows whose score is not a number.
"""
import csv
import glob
import gzip
import hashlib
import json
import os

from ..normalise import normalise, parse_key
from .common import DATA, evidence_row, retrieved
from .hgvs import hgvs_to_key

MV = os.path.join(DATA, 'mavedb', 'dump_v4_extract', 'mavedb-dump.2026062418131')
VEP = os.path.join(DATA, 'ensembl_vep', 'MaveDB_variants.tsv.gz')
ANN = os.path.join(MV, 'csv')
SCORES = os.path.join(DATA, 'mavedb', 'scores_genomic_mapped')
STD = ['#chr', 'start', 'end', 'ref', 'alt', 'hgvs', 'urn', 'publish_date', 'refseq', 'pubmed', 'doi', 'url', 'accession', 'nt', 'splice', 'pro', 'score', 'sd', 'se', 'df']
KEYWORD_MODEL = 'Phenotypic Assay Model System'
CLASS_EFFECT = {'abnormal': 'effect', 'normal': 'no_effect', 'not_specified': 'uncertain'}


def in_range(x, lo, hi, inc_lo, inc_hi):
    if lo is not None and not (x > lo or (inc_lo and x == lo)):
        return False
    if hi is not None and not (x < hi or (inc_hi and x == hi)):
        return False
    return True


def csv_records(fh):
    """(first line number, row as dict, exact record text) for a header-first CSV; a record that spans several physical lines keeps all of them."""
    pulled = []

    def lines():
        for line in fh:
            pulled.append(line)
            yield line
    rd = csv.reader(lines())
    header = next(rd)
    pulled.clear()
    ln = 2
    for row in rd:
        text = ''.join(pulled)
        n = len(pulled)
        pulled.clear()
        yield ln, dict(zip(header, row)), text.rstrip('\r\n'), header
        ln += n


def num(s):
    try:
        v = float(s)
    except (TypeError, ValueError):
        return None
    return v if v == v and abs(v) != float('inf') else None


def std_error(x):
    """A standard error is positive. Some score sets write 0.0 for every variant, which says that no error estimate was given; the raw value stays in extra_json."""
    v = num(x)
    return v if v is not None and v > 0 else None


def short_hash(*parts):
    return int.from_bytes(hashlib.blake2b('|'.join(parts).encode(), digest_size=8).digest(), 'big')


class MaveDbAdapter:
    source_id = 'mavedb'
    targets = ['ensembl_vep_mavedb_variants', 'mavedb_zenodo_dump_v4', 'mavedb_scores_genomic_mapped']

    def __init__(self, genes, vep_path=None, main_json=None, ann_dir=None, scores_dir=None):
        self.genes = genes
        self.vep_path = vep_path or VEP
        self.main_json = main_json or os.path.join(MV, 'main.json')
        self.ann_dir = ann_dir if ann_dir is not None else ANN
        self.scores_dir = scores_dir if scores_dir is not None else SCORES
        self.sets = {}

    # ---- assay metadata ---------------------------------------------------------------------------------------------------------
    def load_metadata(self, ctx):
        with open(self.main_json, encoding='utf-8') as f:
            main = json.load(f)
        self.as_of = main.get('asOf', '')[:10]
        for es in main['experimentSets']:
            for ex in es['experiments']:
                kw = {}
                for k in ex.get('keywords') or ():
                    kk = k.get('keyword') or {}
                    kw.setdefault(kk.get('key'), []).append(kk.get('label'))
                pubs = [p['identifier'] for p in (ex.get('primaryPublicationIdentifiers') or ()) if p.get('dbName') == 'PubMed']
                for s in ex['scoreSets']:
                    spubs = [p['identifier'] for p in (s.get('primaryPublicationIdentifiers') or ()) if p.get('dbName') == 'PubMed'] or pubs
                    tg = (s.get('targetGenes') or [{}])[0]
                    ensg = None
                    for x in tg.get('externalIdentifiers') or ():
                        ident = x.get('identifier') or {}
                        if ident.get('dbName') == 'Ensembl' and str(ident.get('identifier', '')).startswith('ENSG'):
                            ensg = ident['identifier'].split('.')[0]
                    symbol = tg.get('mappedHgncName') or tg.get('name')
                    if not ensg and symbol:
                        ensg = self.genes.resolve(symbol=symbol)[0]
                    cals = []
                    for c in s.get('scoreCalibrations') or ():
                        for fc in c.get('functionalClassifications') or ():
                            rng = fc.get('range') or [None, None]
                            acmg = fc.get('acmgClassification') or {}
                            cals.append(dict(calibration=c.get('urn'), title=c.get('title'), research_use_only=c.get('researchUseOnly'), investigator_provided=c.get('investigatorProvided'),
                                             label=fc.get('label'), cls=fc.get('functionalClassification'), criterion=acmg.get('criterion'), strength=acmg.get('evidenceStrength'),
                                             points=acmg.get('points'), lo=rng[0], hi=rng[1], inc_lo=bool(fc.get('inclusiveLowerBound')), inc_hi=bool(fc.get('inclusiveUpperBound')),
                                             sources=[x.get('identifier') for x in c.get('thresholdSources') or ()]))
                    self.sets[s['urn']] = dict(urn=s['urn'], experiment=ex['urn'], title=s.get('title'), pubmed=spubs[0] if spubs else None, gene=symbol, ensg=ensg,
                                               model=', '.join(kw.get(KEYWORD_MODEL, [])) or None, mechanism=', '.join(kw.get('Phenotypic Assay Mechanism', [])) or None,
                                               method=', '.join(kw.get('Phenotypic Assay Method', [])) or None, library=', '.join(kw.get('Variant Library Creation Method', [])) or None,
                                               superseded_by=s.get('supersedingScoreSet'), mapping=s.get('mappingState'), n_variants=s.get('numVariants'), calibrations=cals,
                                               category=tg.get('category'), published=s.get('publishedDate'))

    def emit_assays(self, ctx, rel_dump, st):
        for urn, s in sorted(self.sets.items()):
            kind = 'endogenous_locus_variant_library' if (s['library'] or '').startswith('Endogenous') else 'construct_variant_library'
            ctx.emit('entity_assay', dict(assay_id=urn, source_id='mavedb', kind=kind, target_gene=s['ensg'], model_system=s['model'], tissue_or_cell=None, perturbation='variant library',
                                          dose_context=None, readout='; '.join(x for x in (s['method'], s['mechanism']) if x) or None, score_direction=None,
                                          lineage_group=self.lineage(s), release_id=rel_dump), name='mavedb')
            st.inc('assays')
            if s['superseded_by']:
                st.inc('assays_superseded')

    @staticmethod
    def lineage(s):
        return f"pmid:{s['pubmed']}" if s['pubmed'] else f"mavedb_experiment:{s['experiment']}"

    # ---- rows -------------------------------------------------------------------------------------------------------------------
    def run(self, ctx, limit=None):
        ctx.adapter = self
        st = ctx.st()
        self.load_metadata(ctx)
        rel_dump = ctx.register_release('mavedb', 'mavedb_zenodo_dump_v4', self.as_of or 'undated')
        rel = ctx.register_release('mavedb', 'ensembl_vep_mavedb_variants', retrieved(ctx, 'ensembl_vep_mavedb_variants'))
        self.emit_assays(ctx, rel_dump, st)
        opener = gzip.open if self.vep_path.endswith('.gz') else open
        seen, vg_seen = set(), set()
        with opener(self.vep_path, 'rt', encoding='utf-8', newline='') as f:
            hdr = f.readline().rstrip('\n').split('\t')
            ix = {n: i for i, n in enumerate(hdr)}
            missing = [n for n in STD if n not in ix]
            if missing:
                raise ValueError(f'unexpected MaveDB plugin file header, missing {missing}')
            extra_names = hdr[len(STD):]
            for ln, line in enumerate(f, start=2):
                if limit and st.c.get('rows', 0) >= limit:
                    break
                st.inc('rows')
                raw = line.rstrip('\n')
                c = raw.split('\t', len(STD))
                chrom, start, end, ref, alt, hgvs, urn, publish, _refseq, pubmed, _doi, _url, accession, nt, splice, pro, score_s, sd, se, df = c[:20] + [''] * (20 - len(c[:20]))
                tail = c[20] if len(c) > 20 else ''
                s = self.sets.get(urn)
                if s is None:
                    st.inc('rows_score_set_not_in_dump_metadata')
                    continue
                if s['superseded_by']:
                    st.inc('rows_superseded_score_set')
                    continue
                if ref.startswith('ga4gh:') or alt.startswith('ga4gh:') or ';' in nt or '[' in nt or '[' in pro or '[' in hgvs:
                    st.inc('rows_multi_variant_construct')
                    continue
                score = num(score_s)
                if score is None:
                    st.inc('rows_score_not_numeric')
                    continue
                canon = ctx.ref.contig(chrom)
                n = normalise(ctx.ref, canon, int(start), ref, alt) if canon and start.isdigit() else None
                if n is not None and n.ok:
                    key, mapping, reason = n.key, 'mapped', None
                else:
                    key, mapping, reason = None, ('rejected' if n is not None and n.status == 'ref_mismatch' else 'unmapped'), (n.status if n is not None else 'contig_or_position_unusable')
                dedupe = short_hash(accession, key or f'{chrom}:{start}:{ref}:{alt}')
                if dedupe in seen:
                    st.inc('rows_duplicate_allele_for_accession')
                    continue
                seen.add(dedupe)
                allele_level = bool(nt)
                claim_type = 'mavedb_variant_score' if allele_level else 'mavedb_protein_level_score'
                extras = {}
                if tail:
                    for name, v in zip(extra_names, tail.split('\t')):
                        if v.strip() and len(extras) < 12:
                            extras[name] = v
                extra = dict(protein_hgvs=hgvs if hgvs.startswith('NP_') else None, assay_hgvs_nt=nt or None, assay_hgvs_pro=pro or None, assay_hgvs_splice=splice or None, sd=num(sd), se=num(se), df=num(df),
                             score_set=urn, score_set_title=s['title'], published=publish, target_gene=s['gene'], mapping_basis='ensembl_vep_mavedb_plugin_file', other_score_columns=extras or None)
                if not allele_level:
                    extra['scope'] = 'amino-acid change measured; the nucleotide allele is one that encodes it (genetic code), effects outside the coding sequence were not tested'
                lineage = self.lineage(s)
                common = dict(chrom=n.chrom if key else canon, pos=n.pos if key else None, unmapped_key=None if key else f'{chrom}:{start}:{ref}:{alt}', gene_ids=[s['ensg']] if s['ensg'] else None,
                              condition_text=f"{s['gene'] or ''} {s['title'] or ''}".strip() or None, model_system=s['model'], perturbation='variant library', lineage_group=lineage,
                              lineage_note=f"MaveDB score set {urn}" + (f' (PMID {s["pubmed"]})' if s['pubmed'] else ''), evaluated=publish or None, origin='germline')
                self_reason = reason or ('protein_level_assay_allele_inferred_from_genetic_code' if not allele_level else None)
                ctx.emit_evidence(evidence_row(rel, 'mavedb', accession, f'MaveDB_variants.tsv.gz#L{ln}', raw, 'variant', key, 'empirical_functional_support', 'direct_experimental', claim_type,
                                               f'{score:.6g}', 'not_applicable', mapping, self_reason if mapping != 'mapped' or not allele_level else None, discriminator=key or f'{chrom}:{start}:{ref}:{alt}',
                                               claim_num=score, score=score, score_unit=f'mavedb_score:{urn}', se=std_error(se), extra_json=json.dumps(extra, sort_keys=True, ensure_ascii=False), **common))
                st.inc(f'score_rows:{"allele" if allele_level else "protein"}:{mapping}')
                if key and s['ensg'] and (key, s['ensg']) not in vg_seen:
                    vg_seen.add((key, s['ensg']))
                    ctx.emit('variant_gene', dict(variant_key=key, gene_id=s['ensg'], relation='mavedb_target_gene', consequence=None, release_id=rel), name='mavedb')
                self.calibrations(ctx, st, rel, s, score, accession, f'MaveDB_variants.tsv.gz#L{ln}', raw, key, mapping, self_reason if mapping != 'mapped' or not allele_level else None, extra, common)
        self.native_pass(ctx, st, seen, vg_seen)
        return st.c

    def calibrations(self, ctx, st, rel, s, score, accession, locator, raw, key, mapping, reason, extra, common):
        for cl in s['calibrations']:
            if cl['cls'] and in_range(score, cl['lo'], cl['hi'], cl['inc_lo'], cl['inc_hi']):
                value = '|'.join(str(x) for x in (cl['cls'], cl['criterion'] or '', cl['strength'] or ''))
                ctx.emit_evidence(evidence_row(rel, 'mavedb', accession, locator, raw, 'variant', key, 'empirical_functional_support', 'curated_expert_classification', 'mavedb_calibrated_class', value,
                                               CLASS_EFFECT.get(cl['cls'], 'uncertain'), mapping, reason, discriminator=f"{key or ''}|{cl['calibration']}|{cl['label']}", curator=cl['title'],
                                               review_status='research use only' if cl['research_use_only'] else None,
                                               extra_json=json.dumps(dict(extra, calibration=cl['calibration'], label=cl['label'], points=cl['points'], range=[cl['lo'], cl['hi']], threshold_sources=cl['sources'],
                                                                          investigator_provided=cl['investigator_provided']), sort_keys=True, ensure_ascii=False),
                                               **{**common, 'lineage_note': common['lineage_note'] + '; an interpretation of the score row, not an additional study'}))
                st.inc('calibrated_class_rows')

    # ---- MaveDB's own genomic mapping ---------------------------------------------------------------------------------------------
    def native_pass(self, ctx, st, seen, vg_seen):
        if not self.ann_dir or not os.path.isdir(self.ann_dir) or not os.path.isdir(self.scores_dir):
            st.inc('native_mapping_not_available')
            return
        rel = ctx.register_release('mavedb', 'mavedb_scores_genomic_mapped', self.as_of or 'undated')
        for sp in sorted(glob.glob(os.path.join(self.scores_dir, '*.scores.csv.gz'))):
            stem = os.path.basename(sp)[:-len('.scores.csv.gz')]
            urn = stem.replace('urn-mavedb-', 'urn:mavedb:', 1)
            s = self.sets.get(urn)
            ap = os.path.join(self.ann_dir, stem + '.annotations.csv')
            if s is None or not os.path.exists(ap):
                st.inc('native_score_files_without_metadata_or_annotations')
                continue
            if s['superseded_by']:
                st.inc('native_rows_superseded_score_set')
                continue
            scores = {}
            with gzip.open(sp, 'rt', encoding='utf-8', newline='') as f:
                for ln, row, _text, _hdr in csv_records(f):
                    scores[row['accession']] = (ln, row)
            with open(ap, encoding='utf-8', newline='') as f:
                for ln, row, raw, _hdr in csv_records(f):
                    g = row.get('mavedb.post_mapped_hgvs_g', 'NA')
                    if g in ('NA', ''):
                        continue
                    st.inc('native_rows')
                    if '[' in g or ';' in g or ',' in g:
                        st.inc('native_rows_multi_variant')
                        continue
                    acc = row['accession']
                    sc = scores.get(acc)
                    score = num(sc[1].get('scores.score', sc[1].get('score'))) if sc else None
                    if score is None:
                        st.inc('native_rows_without_numeric_score')
                        continue
                    contig, _, change = g.partition(':')
                    canon = ctx.ref.contig(contig)
                    status, keys, why = hgvs_to_key(ctx.ref, canon, change) if canon else ('unresolved', [], 'contig_unmapped')
                    if status == 'variant':
                        key, mapping, reason = keys[0], 'mapped', None
                    elif status == 'reference':
                        st.inc('native_rows_reference_sequence')
                        continue
                    else:
                        key, mapping, reason = None, 'rejected' if (why or '').startswith('reference_base_differs') else 'unmapped', why or 'unresolved'
                    dedupe = short_hash(acc, key or g)
                    if dedupe in seen:
                        st.inc('rows_duplicate_allele_for_accession')
                        continue
                    seen.add(dedupe)
                    lineage = self.lineage(s)
                    extra = dict(assay_hgvs_nt=row.get('hgvs_nt') if row.get('hgvs_nt') != 'NA' else None, assay_hgvs_pro=row.get('hgvs_pro') if row.get('hgvs_pro') != 'NA' else None,
                                 post_mapped_hgvs_g=g, post_mapped_hgvs_c=row.get('mavedb.post_mapped_hgvs_c') if row.get('mavedb.post_mapped_hgvs_c') != 'NA' else None,
                                 post_mapped_hgvs_p=row.get('mavedb.post_mapped_hgvs_p') if row.get('mavedb.post_mapped_hgvs_p') != 'NA' else None, vrs_digest=row.get('mavedb.post_mapped_vrs_digest'),
                                 vep_consequence=row.get('vep.vep_functional_consequence') if row.get('vep.vep_functional_consequence') != 'NA' else None, sd=num(sc[1].get('scores.sd', sc[1].get('sd'))), se=num(sc[1].get('scores.se', sc[1].get('se'))),
                                 score_set=urn, score_set_title=s['title'], target_gene=s['gene'], mapping_basis='mavedb_post_mapped_hgvs_g', score_locator=f'{stem}.scores.csv.gz#L{sc[0]}')
                    common = dict(chrom=None, pos=None, unmapped_key=None if key else g, gene_ids=[s['ensg']] if s['ensg'] else None, condition_text=f"{s['gene'] or ''} {s['title'] or ''}".strip() or None,
                                  model_system=s['model'], perturbation='variant library', lineage_group=lineage, lineage_note=f"MaveDB score set {urn}" + (f' (PMID {s["pubmed"]})' if s['pubmed'] else ''),
                                  evaluated=s['published'], origin='germline')
                    if key:
                        common['chrom'], common['pos'] = parse_key(key)[:2]
                    ctx.emit_evidence(evidence_row(rel, 'mavedb', acc, f'{stem}.annotations.csv#L{ln}', raw, 'variant', key, 'empirical_functional_support', 'direct_experimental', 'mavedb_variant_score',
                                                   f'{score:.6g}', 'not_applicable', mapping, reason, discriminator=key or g, claim_num=score, score=score, score_unit=f'mavedb_score:{urn}', se=std_error(sc[1].get('scores.se', sc[1].get('se'))),
                                                   extra_json=json.dumps(extra, sort_keys=True, ensure_ascii=False), **common))
                    st.inc(f'score_rows:native:{mapping}')
                    if key and s['ensg'] and (key, s['ensg']) not in vg_seen:
                        vg_seen.add((key, s['ensg']))
                        ctx.emit('variant_gene', dict(variant_key=key, gene_id=s['ensg'], relation='mavedb_target_gene', consequence=extra['vep_consequence'], release_id=rel), name='mavedb')
                    self.calibrations(ctx, st, rel, s, score, acc, f'{stem}.annotations.csv#L{ln}', raw, key, mapping, reason, extra, common)
