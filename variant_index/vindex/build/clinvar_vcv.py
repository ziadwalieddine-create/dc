"""ClinVar full VCV release (stream-derived Parquet tables) -> submission-level evidence, structural events and composite records.

The GRCh38 VCF adapter gives one aggregate row per simple variant. The VCV release adds what the aggregate hides or the VCF cannot hold:

  submissions (SCV)   one evidence row per submission: submitter, classification as submitted, review status, date last evaluated, condition, method,
                      citations, origin, affected status, number tested. A submission without a classification (ClinVar converts observations of
                      functional data and phenotype-only submissions to 'no classification provided') is kept as an observation row: functional
                      support when the method is in vitro / in vivo, otherwise context. Submissions and their aggregate share one lineage group, so a
                      view counting independent sources counts ClinVar once for a variant; the number of submitters and the disagreement between them
                      stay visible as rows. Rows attach to the variant key the VCF adapter assigned to the VariationID (that adapter must run first), to
                      the structural event below, or stay unmapped with the reason.
  structural events   GRCh38 locations of copy-number gains and losses (any size) and of deletions, duplications, insertions, inversions, indels and
                      complex rearrangements of at least 50 bp become entity_event rows. Uncertain extents: the inner extent is the nominal event and
                      the outer limits become confidence intervals; an unbounded limit ('?') leaves the interval open and the event imprecise. Copy
                      number is taken from the record name only when it states one value. Records that only have GRCh37 coordinates stay unmapped with
                      those coordinates in extra_json: no lift-over is performed.
  composite records   haplotypes, diplotypes, compound heterozygotes and phase-unknown sets: entity_haplotype (system 'clinvar', allele_level
                      'composite:<type>') whose definition rows are the member alleles that have a VCF key; unresolved members are listed, never
                      dropped. The classification belongs to the combination, never to one member. These are not star-allele systems; the haplotype
                      caller ignores them.
Every row's raw_sha256 is the SHA-256 of the canonical JSON (sorted keys) of the complete source row and the locator is '<table>/<part file>#<row index>'.
"""
import glob
import json
import os
import re

import pyarrow as pa
import pyarrow.compute as pc
import pyarrow.parquet as pq

from .classes import effect_of, stars, unknown_tokens
from .common import DATA, evidence_row
from .genes import entity_id

VCV = os.path.join(DATA, 'clinvar', 'vcv')
EVENT_TYPE = {'copy number gain': 'CNV_GAIN', 'copy number loss': 'CNV_LOSS', 'Deletion': 'DEL', 'Duplication': 'DUP', 'Tandem duplication': 'DUP', 'Insertion': 'INS',
              'Inversion': 'INV', 'Translocation': 'CTX', 'Indel': 'CPX', 'Complex': 'CPX'}
ANY_SIZE = {'copy number gain', 'copy number loss', 'Translocation'}
COMPOSITE = {'Haplotype': 'defining', 'Diplotype': 'member:genotype', 'CompoundHeterozygote': 'member:genotype', 'Phase unknown': 'member:phase_unknown',
             'Distinct chromosomes': 'member:distinct_chromosomes'}
MIN_EVENT_LENGTH = 50
FAMILIES = (('germline', 'germline_classification', 'germline'), ('somatic', 'somatic_clinical_impact', 'somatic'), ('oncogenicity', 'oncogenicity_classification', 'somatic'))
SCV_CLASS = {'GermlineClassification': ('scv_germline_classification', 'germline'), 'OncogenicityClassification': ('scv_oncogenicity_classification', 'somatic'),
             'SomaticClinicalImpact': ('scv_somatic_clinical_impact', 'somatic')}
AGG_NOTE = 'aggregate of ClinVar submissions (SCV); not independent of them, and submitters may share cohorts'
SCV_NOTE = 'one submission; counts with the aggregate and the other submissions of the same variant as one source; submitters may share cohorts'
OBS_NOTE = 'functional or phenotype observation without a classification; the cited study is the unit of independence'
CN_NAME = re.compile(r'\)x(\d+)(,\d+)*$')


def num(s):
    try:
        return int(float(s))
    except (TypeError, ValueError):
        return None


def canon(row):
    """Canonical text of one complete source row: what raw_sha256 digests, so a verifier can recompute it from the locator."""
    return json.dumps(row, sort_keys=True, ensure_ascii=False, default=str)


def condition_ids(text):
    """'HP:HP:0002282|MedGen:CN002072#Name' -> ['HP:0002282', 'MedGen:CN002072']; the source repeats the database prefix of some identifiers."""
    ids = (text or '').partition('#')[0]
    return [re.sub(r'^(\w+):\1:', r'\1:', x.strip()) for x in ids.split('|') if x.strip()]


def condition_name(text):
    return (text or '').partition('#')[2].strip() or None


def pubmed_ids(citations):
    return [c.split(':', 1)[1] for c in (citations or '').split('|') if c.lower().startswith(('pubmed:', 'pmid:'))]


def origin_label(s):
    t = [x.strip().replace(' ', '_') for x in (s or '').split('|') if x.strip() and x.strip() not in ('not applicable', 'not provided')]
    return '+'.join(t) or None


class ClinVarVcvAdapter:
    source_id = 'clinvar'
    targets = ['clinvar_vcv_xml_weekly_2026_1018']
    requires = ('clinvar',)

    def __init__(self, genes=None, root=None):
        self.genes, self.root = genes, root or VCV
        self.locs, self.legacy, self.events, self.composites = {}, {}, {}, {}

    def parts(self, table):
        return sorted(glob.glob(os.path.join(self.root, table, '*.parquet')))

    def prepare(self, ctx):
        ctx.want('clinvar_variation')
        ctx.want('clinvar_variation_unmapped')
        alleles = set()
        for f in self.parts('variation'):
            t = pq.read_table(f, columns=['variation_type', 'allele_ids'])
            t = t.filter(pc.is_in(t['variation_type'], value_set=pa.array(sorted(COMPOSITE))))
            for s in t['allele_ids'].to_pylist():
                alleles.update(x for x in (s or '').split('|') if x)
        ctx.want('clinvar_allele', alleles)

    def run(self, ctx, limit=None):
        ctx.adapter = self
        st = ctx.st()
        first = self.parts('variation')
        version = pq.read_table(first[0], columns=['release_date'])['release_date'][0].as_py() if first else None
        self.rel = ctx.register_release('clinvar', self.targets[0], version or 'undated')
        self.scan_locations(ctx, st)
        self.read_variations(ctx, st)
        self.read_scvs(ctx, st, limit)
        return st.c

    # ---- gene ids -----------------------------------------------------------------------------------------------------------------
    def gene_ids(self, text):
        out = []
        for tok in (text or '').split('|'):
            sym, _, nid = tok.partition(':')
            gid = entity_id(self.genes, symbol=sym or None, ncbi=nid or None) if nid else None
            if gid and gid not in out:
                out.append(gid)
        return out or None

    # ---- locations -> candidate events and legacy coordinates ---------------------------------------------------------------------
    def scan_locations(self, ctx, st):
        cols = ['variation_id', 'allele_type', 'assembly', 'chr', 'start', 'stop', 'inner_start', 'inner_stop', 'outer_start', 'outer_stop', 'variant_length']
        for f in self.parts('allele_location'):
            t = pq.read_table(f, columns=cols)
            t = t.filter(pc.is_in(t['allele_type'], value_set=pa.array(sorted(EVENT_TYPE))))
            for r in t.to_pylist():
                vid, typ = r['variation_id'], r['allele_type']
                s, e = num(r['start']), num(r['stop'])
                istart, istop, ostart, ostop = (num(r[k]) for k in ('inner_start', 'inner_stop', 'outer_start', 'outer_stop'))
                precise = s is not None and e is not None
                s0 = s if s is not None else (istart if istart is not None else ostart)
                e0 = e if e is not None else (istop if istop is not None else ostop)
                length = num(r['variant_length']) or (e0 - s0 + 1 if s0 is not None and e0 is not None and e0 >= s0 else None)
                if r['assembly'] != 'GRCh38':
                    if r['assembly'] == 'GRCh37' and vid not in self.legacy:
                        self.legacy[vid] = f"GRCh37:{r['chr']}:{s0}-{e0}"
                    continue
                st.inc('grch38_event_type_locations')
                if vid in self.locs:
                    st.inc('records_with_several_grch38_locations')
                    continue
                if typ not in ANY_SIZE and (length is None or length < MIN_EVENT_LENGTH):
                    st.inc('locations_below_event_length')
                    continue
                chrom = ctx.ref.contig(r['chr'])
                if chrom is None or s0 is None or e0 is None or e0 < s0:
                    st.inc('event_locations_unusable')
                    continue
                ci = None if precise else dict(start_ci_lo=(ostart - s0) if ostart is not None and istart is not None else None, start_ci_hi=0 if istart is not None else None,
                                               end_ci_lo=0 if istop is not None else None, end_ci_hi=(ostop - e0) if ostop is not None and istop is not None else None)
                self.locs[vid] = dict(type=typ, chrom=chrom, start=s0, end=e0, precise=precise, length=None if typ == 'Translocation' else length, ci=ci, inner=[istart, istop], outer=[ostart, ostop])
        for vid in list(self.legacy):
            if vid in self.locs:
                del self.legacy[vid]

    # ---- variation table: events, composite records, unmapped aggregates ------------------------------------------------------------
    def read_variations(self, ctx, st):
        for f in self.parts('variation'):
            fname = os.path.basename(f)
            t = pq.read_table(f)
            ids, types = t['variation_id'].to_pylist(), t['variation_type'].to_pylist()
            found = ctx.found
            need = [i for i, (v, ty) in enumerate(zip(ids, types)) if v in self.locs or ty in COMPOSITE
                    or (('clinvar_variation', v) not in found and ('clinvar_variation_unmapped', v) not in found)]
            st.inc('variation_rows', len(ids))
            for i, a in zip(need, t.take(need).to_pylist()):
                vid, typ, loc = a['variation_id'], a['variation_type'], f'variation/{fname}#{i}'
                raw = canon(a)
                if typ in COMPOSITE:
                    self.composite(ctx, st, a, loc, raw)
                elif vid in self.locs:
                    self.event(ctx, st, a, loc, raw, found.get(('clinvar_variation', vid)))
                elif a['record_type'] == 'included' and not any(a.get(f'{fam}_description') for fam, _c, _o in FAMILIES):
                    st.inc('included_records_without_classification_unkeyed')
                else:
                    self.unmapped_aggregate(ctx, st, a, loc, raw)

    def aggregates(self, ctx, st, a, loc, raw, level, subject, mapping, reason, **kw):
        n = 0
        for fam, ctype, origin in FAMILIES:
            desc = a.get(f'{fam}_description')
            if not desc:
                continue
            for tok in unknown_tokens(desc):
                st.inc(f'unrecognised_classification:{tok}')
            rev = a.get(f'{fam}_review_status') or ''
            ctx.emit_evidence(evidence_row(self.rel, 'clinvar', a['accession'], loc, raw, level, subject, 'clinical_validity_actionability', 'curated_expert_classification', ctype, desc,
                                           effect_of(desc), mapping, reason, discriminator=fam, review_status=rev or None, review_stars=stars(rev), n_submitters=num(a['number_of_submitters']),
                                           evaluated=a.get(f'{fam}_date_last_evaluated') or None, origin=origin, lineage_group=f"clinvar:variation:{a['variation_id']}", lineage_note=AGG_NOTE,
                                           gene_ids=self.gene_ids(a['genes']), **kw))
            n += 1
        st.inc(f'aggregate_rows:{level}:{mapping}', n)

    def event(self, ctx, st, a, loc, raw, keys):
        vid, l = a['variation_id'], self.locs[a['variation_id']]
        eid = f"clinvar:{a['accession']}"
        self.events[vid] = eid
        m = CN_NAME.search(a['name'] or '')
        copy_number = int(m.group(1)) if m and not m.group(2) else None
        one_key = sorted(set(keys))[0] if keys and len(set(keys)) == 1 else None
        extra = dict(variation_type=a['variation_type'], name=a['name'], genes=a['genes'], inner=l['inner'], outer=l['outer'], hgvs_grch38=a['hgvs_grch38'], variant_key=one_key,
                     copy_state_text=(m.group(0) if m and m.group(2) else None),
                     note='copy number appears only when the record name states a single value; breakpoints are uncertain unless precise')
        ci = l['ci'] or {}
        ctx.emit('entity_event', dict(event_id=eid, source_id='clinvar', release_id=self.rel, accession=a['accession'], event_type=EVENT_TYPE[l['type']], chrom=l['chrom'], start=l['start'],
                                      end=l['end'], start_ci_lo=ci.get('start_ci_lo'), start_ci_hi=ci.get('start_ci_hi'), end_ci_lo=ci.get('end_ci_lo'), end_ci_hi=ci.get('end_ci_hi'),
                                      svlen=l['length'], copy_number=copy_number, precise=l['precise'], population_af=None, n_samples=None, origin='clinical_submission',
                                      raw_type=a['variation_type'], assembly='GRCh38',
                                      extra_json=json.dumps({k: v for k, v in extra.items() if v not in (None, [None, None], '')}, sort_keys=True)), name='clinvar_vcv')
        ctx.emit('xref', dict(namespace='clinvar_variation', identifier=vid, entity_level='event', entity_id=eid, release_id=self.rel, relation='vcv_allele_location'), name='clinvar_vcv')
        st.inc(f"events:{EVENT_TYPE[l['type']]}")
        self.aggregates(ctx, st, a, loc, raw, 'event', eid, 'mapped', None, chrom=l['chrom'], pos=l['start'], end=l['end'])

    def composite(self, ctx, st, a, loc, raw):
        vid, typ = a['variation_id'], a['variation_type']
        hid = f"clinvar:{a['accession']}"
        self.composites[vid] = hid
        members = [m for m in (a['allele_ids'] or '').split('|') if m]
        keys, unresolved = [], []
        for m in members:
            ks = set(ctx.found.get(('clinvar_allele', m)) or ())
            if len(ks) == 1:
                k = next(iter(ks))
                if k not in keys:
                    keys.append(k)
            elif m not in unresolved:
                unresolved.append(m)
        symbols = {t.partition(':')[0] for t in (a['genes'] or '').split('|') if t}
        ctx.emit('entity_haplotype', dict(haplotype_id=hid, system='clinvar', gene_symbol=next(iter(symbols)) if len(symbols) == 1 else None, allele_name=a['name'], is_reference=False,
                                          structural=None, release_id=self.rel, accession=a['accession'], function=None, allele_level=f'composite:{typ}'), name='clinvar_vcv')
        for k in keys:
            ctx.emit('haplotype_definition', dict(haplotype_id=hid, variant_key=k, hgvs=None, release_id=self.rel, definition_scope=COMPOSITE[typ]), name='clinvar_vcv')
        for m in unresolved:
            ctx.emit('haplotype_definition', dict(haplotype_id=hid, variant_key=None, hgvs=f'AlleleID:{m}', release_id=self.rel, definition_scope='unresolved:member_has_no_vcf_key'), name='clinvar_vcv')
        ctx.emit('xref', dict(namespace='clinvar_variation', identifier=vid, entity_level='haplotype', entity_id=hid, release_id=self.rel, relation='vcv_composite_record'), name='clinvar_vcv')
        st.inc(f'composites:{typ}')
        st.inc('composite_members_resolved', len(keys))
        st.inc('composite_members_unresolved', len(unresolved))
        self.aggregates(ctx, st, a, loc, raw, 'haplotype', hid, 'mapped', None,
                        extra_json=json.dumps(dict(record_type=typ, member_alleles=members, member_keys=keys, unresolved_members=unresolved,
                                                   note='classification of the combination of alleles, not of any one member'), sort_keys=True))

    def unmapped_aggregate(self, ctx, st, a, loc, raw):
        vid, typ = a['variation_id'], a['variation_type']
        if not any(a.get(f'{fam}_description') for fam, _c, _o in FAMILIES):
            st.inc('unkeyed_records_without_classification')
            return
        level = 'event' if typ in EVENT_TYPE else 'variant'
        reason = ('only_grch37_location_no_liftover_performed' if vid in self.legacy else
                  'no_grch38_event_location_of_at_least_50bp' if level == 'event' else 'no_vcf_key_for_variation')
        extra = dict(variation_type=typ, name=a['name'], hgvs_grch38=a['hgvs_grch38'], canonical_spdi=a['canonical_spdi'], legacy_location=self.legacy.get(vid), genes=a['genes'])
        self.aggregates(ctx, st, a, loc, raw, level, None, 'unmapped', reason, unmapped_key=f'VariationID:{vid}', extra_json=json.dumps({k: v for k, v in extra.items() if v}, sort_keys=True))

    # ---- submissions ---------------------------------------------------------------------------------------------------------------
    def resolve(self, ctx, vid):
        """(level, subject, mapping, reason, candidate keys): composite records first, then the VCF key, then the structural event."""
        if vid in self.composites:
            return 'haplotype', self.composites[vid], 'mapped', None, None
        keys = set(ctx.found.get(('clinvar_variation', vid)) or ())
        if len(keys) == 1:
            return 'variant', next(iter(keys)), 'mapped', None, None
        if keys:
            return 'variant', None, 'ambiguous', 'variation_id_maps_to_several_keys', sorted(keys)
        if vid in self.events:
            return 'event', self.events[vid], 'mapped', None, None
        if vid in self.legacy:
            return 'event', None, 'unmapped', 'only_grch37_location_no_liftover_performed', None
        why = ctx.found.get(('clinvar_variation_unmapped', vid))
        if why:
            return 'variant', None, 'unmapped', f'vcf_record_unmapped:{why[0]}', None
        return 'variant', None, 'unmapped', 'no_vcf_key_or_structural_event_for_variation', None

    def scv_rows(self, f, batch=50_000):
        """(row index within the part file, row dict); batches bound memory (the free-text comments make this the widest table)."""
        i = 0
        for b in pq.ParquetFile(f).iter_batches(batch_size=batch):
            for r in b.to_pylist():
                yield i, r
                i += 1

    def read_scvs(self, ctx, st, limit):
        for f in self.parts('scv'):
            fname = os.path.basename(f)
            for i, r in self.scv_rows(f):
                if limit and st.c.get('scv_rows', 0) >= limit:
                    return
                st.inc('scv_rows')
                vid = r['variation_id']
                level, subject, mapping, reason, cands = self.resolve(ctx, vid)
                label = (r['classification'] or '').strip()
                pm = pubmed_ids(r['citations'])
                if r['classification_type']:
                    ctype, origin = SCV_CLASS[r['classification_type']]
                    for tok in unknown_tokens(label):
                        st.inc(f'unrecognised_classification:{tok}')
                    axis, kind, claim, effect, lineage = 'clinical_validity_actionability', 'curated_expert_classification', label or 'not stated', effect_of(label), f'clinvar:variation:{vid}'
                    origin, note = origin_label(r['origin']) or origin, SCV_NOTE
                else:
                    methods = set((r['method_types'] or '').split('|'))
                    axis, kind = ('empirical_functional_support', 'experimental_annotation') if methods & {'in vitro', 'in vivo'} else ('applicability_context', 'case_observation')
                    ctype, claim, effect, origin, note = 'scv_observation_without_classification', r['method_types'] or 'method not stated', 'not_applicable', origin_label(r['origin']), OBS_NOTE
                    lineage = f'pmid:{pm[0]}' if len(pm) == 1 else f"clinvar:scv:{r['scv']}"
                rev = (r['review_status'] or '').strip()
                chrom = pos = None
                if subject and level == 'variant':
                    c, p = subject.split(':')[:2]
                    chrom, pos = c, int(p)
                extra = dict(variation_id=vid, version=r['version'], org_id=r['org_id'], organization_category=r['organization_category'], method_types=r['method_types'],
                             affected_status=r['affected_status'], citations=r['citations'], submission_date=r['submission_date'], number_tested=r['number_tested'],
                             not_in_aggregate=True if r['contributes_to_aggregate'] == 'false' else None, record_status=r['record_status'] if r['record_status'] != 'current' else None,
                             candidate_keys=cands, comment_chars=num(r['comment_length']))
                nt = r['number_tested'] if (r['number_tested'] or '').isdigit() and len(r['number_tested']) < 10 else None
                ctx.emit_evidence(evidence_row(self.rel, 'clinvar', r['scv'], f'scv/{fname}#{i}', canon(r), level, subject, axis, kind, ctype, claim, effect, mapping, reason, discriminator=vid,
                                               curator=r['submitter'], review_status=rev or None, review_stars=stars(rev), assertion_criteria=r['assertion_method'] or None,
                                               evaluated=r['date_last_evaluated'] or None, condition_ids=condition_ids(r['conditions']) or None, condition_text=condition_name(r['conditions']),
                                               sample_size=int(nt) if nt else None, origin=origin, chrom=chrom, pos=pos, unmapped_key=None if subject else f'VariationID:{vid}',
                                               lineage_group=lineage, lineage_note=note, extra_json=json.dumps({k: v for k, v in extra.items() if v not in (None, '')}, sort_keys=True)))
                st.inc(f'scv_evidence:{kind}:{mapping}')
