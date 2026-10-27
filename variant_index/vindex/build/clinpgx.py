"""ClinPGx (formerly PharmGKB) -> evidence rows.

  summary annotations      one row per annotated variant / haplotype / gene descriptor: ClinPGx's own level of evidence (1A..4) and score, kind
                           'aggregated_integration' (it summarises the variant annotations below, so it never counts as an independent study)
  variant annotations      one association row per study finding (drug, phenotype or functional assay) and one effect-size row per reported ratio
                           (OR / HR / RR with interval, p-value, cases and controls, biogeographical group). Studies share the lineage group of their
                           PubMed id, so one paper cited by several annotations or databases counts once.
  guideline annotations    gene-level rows (and rows for the alleles a guideline names): CPIC / DPWG / CPNDS / RNPGx guideline exists, with its flags
  drug labels              gene-level rows: regulator, testing level (required / recommended / actionable / informative)

What an annotation is attached to:
  rsID with a GRCh38 location   the variant key(s) built from ClinPGx's own GRCh38 HGVS synonyms, REF-checked against the FASTA. A site with several
                                alternate alleles is resolved only when the alleles named in the annotation settle one of them (plus strand, or minus
                                strand with the orientation recorded); otherwise the row is ambiguous and keeps the candidates.
  star allele / named allele    the haplotype entity (CPIC definition first, then PharmVar); alleles that are not defined in the index stay unmapped
  'GENE poor metabolizer', ...  a gene-level descriptor: the row is attached to the gene, never to a variant
  anything else                 unmapped with the reason
Genotype and allele text is kept verbatim in extra_json; no effect allele or direction is derived from it.
"""
import json
import os
import re
import zipfile
from collections import Counter
from dataclasses import dataclass

from ..normalise import parse_key
from .common import DATA, evidence_row
from .hgvs import hgvs_to_key

PG = os.path.join(DATA, 'pharmgkb')
COMP = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
LEVEL_EFFECT = {'1A': 'effect', '1B': 'effect', '2A': 'effect', '2B': 'effect', '3': 'effect', '4': 'uncertain'}
GENE_DESCRIPTOR = re.compile(r'(metabolizer|acetylator|null|non-null|deficien|activity|haplotype|expression|function|status|genotype|carrier|variant)', re.I)
STAR = re.compile(r'^[A-Za-z0-9-]+\*[0-9A-Za-z:.\-]+$')
RS = re.compile(r'^rs\d+$')


def letters_of(*texts):
    """Base letters (A, C, G, T) in genotype text such as 'CT + TT' or 'G'; tokens with lowercase letters ('del') and star alleles are ignored."""
    out = set()
    for t in texts:
        for tok in re.findall(r'[A-Za-z]+', t or ''):
            if re.fullmatch(r'[ACGT]+', tok):
                out.update(tok)
    return out


def choose_alt(cands, letters):
    """Which candidate alternate allele of a site an annotation is about, from the bases the annotation names.

    ClinPGx documents that its alleles are mapped to the positive chromosomal strand, which is the orientation of the GRCh38 REF/ALT. The positive
    strand is therefore the reading used; a reading that fits only the opposite strand contradicts the convention and is not used to pick an allele.
    On the current export the convention holds (single-allele sites: positive-strand-only fits outnumber opposite-strand-only fits ~5000:1).
    Returns (status, key, reason, flags) with status mapped | ambiguous | rejected | unmapped."""
    if not cands:
        return 'unmapped', None, 'no_alternate_allele_known', ()
    parsed = {k: parse_key(k) for k in cands}
    snv = [k for k, (_c, _p, r, a) in parsed.items() if len(r) == 1 and len(a) == 1]
    refs = {parsed[k][2] for k in snv}
    if len(cands) == 1:
        flags = ()
        if snv and letters:
            r, a = parsed[snv[0]][2], parsed[snv[0]][3]
            plus = letters <= {r, a}
            minus = {COMP[x] for x in letters} <= {r, a}
            flags = ('annotation_alleles_plus_strand',) if plus and not minus else ('annotation_alleles_minus_strand',) if minus and not plus else \
                    ('annotation_alleles_strand_ambiguous',) if plus and minus else ('annotation_alleles_not_among_site_alleles',)
        return 'mapped', cands[0], None, flags
    if not letters:
        return 'ambiguous', None, 'several_alternate_alleles_annotation_names_no_bases', ()
    if len(refs) != 1:
        return 'ambiguous', None, 'several_alternate_alleles_no_single_base_reference', ()
    ref = next(iter(refs))
    alts = {parsed[k][3] for k in snv}
    plus = letters <= {ref} | alts
    minus = {COMP[x] for x in letters} <= {ref} | alts
    if not plus and not minus:
        return 'rejected', None, 'annotation_alleles_not_among_site_alleles', ()
    if not plus:
        return 'ambiguous', None, 'annotation_alleles_fit_only_the_opposite_strand', ()
    used = letters - {ref}
    if not used:
        return 'ambiguous', None, 'annotation_names_only_the_reference_allele', ()
    if len(used) > 1:
        return 'ambiguous', None, 'annotation_spans_several_alternate_alleles', ()
    alt = next(iter(used))
    flags = ('strand_by_clinpgx_positive_strand_convention',) + (('annotation_alleles_also_fit_the_opposite_strand',) if minus else ()) + \
            (('site_also_has_non_snv_alleles',) if len(snv) < len(cands) else ())
    return 'mapped', next(k for k in snv if k.endswith(f':{ref}:{alt}')), None, flags


@dataclass
class Subject:
    level: str                       # variant | haplotype | gene
    item: str                        # the name as written in the source
    entity: object = None            # entity id when mapped
    status: str = 'unmapped'         # mapped | ambiguous | rejected | unmapped
    reason: object = None
    chrom: object = None
    pos: object = None
    flags: tuple = ()
    candidates: tuple = ()
    also: tuple = ()                 # other entities for the same name (a star allele defined in both CPIC and PharmVar)


def zip_text(path, name):
    with zipfile.ZipFile(path) as z:
        return z.read(name).decode('utf-8', errors='replace')


def zip_version(path):
    with zipfile.ZipFile(path) as z:
        for n in z.namelist():
            m = re.match(r'^CREATED_(\d{4}-\d{2}-\d{2})\.txt$', n)
            if m:
                return m.group(1)
    return 'undated'


def tsv_rows(text):
    """(line number, row dict, raw line) for a header-first TSV whose values may contain double quotes."""
    lines = text.split('\n')
    hdr = lines[0].rstrip('\r').split('\t')
    for ln, line in enumerate(lines[1:], start=2):
        if not line.strip():
            continue
        cells = line.rstrip('\r').split('\t')
        if len(cells) < len(hdr):
            cells += [''] * (len(hdr) - len(cells))
        yield ln, dict(zip(hdr, cells)), line


def num(s):
    try:
        return float(str(s).strip())
    except (TypeError, ValueError):
        return None


def parse_p(s):
    """('= 0.231' | '< 0.001' | '= 1E-5' ...) -> (exact p or None, text)."""
    t = (s or '').strip()
    m = re.match(r'^=\s*([0-9.]+(?:[eE][-+]?\d+)?)$', t)
    if m:
        v = num(m.group(1))
        if v is not None and 0 <= v <= 1:
            return v, t
    return None, t or None


class ClinPgxAdapter:
    source_id = 'pharmgkb'
    targets = ['clinpgx_summaryAnnotations', 'clinpgx_clinicalAnnotations', 'clinpgx_variantAnnotations', 'clinpgx_variants', 'clinpgx_clinicalVariants',
               'clinpgx_guidelineAnnotations', 'clinpgx_drugLabels']
    requires = ('cpic',)

    def __init__(self, genes, root=None):
        self.genes = genes
        self.root = root or PG
        self.sites = {}

    def zp(self, name):
        return os.path.join(self.root, name + '.zip')

    # ---- preparation: which identifiers other adapters must offer ---------------------------------------------------------------
    def prepare(self, ctx):
        self.variants_text = zip_text(self.zp('variants'), 'variants.tsv')
        names, rsids = set(), set()
        self.rows_cache = {}
        for z, member in (('summaryAnnotations', 'summary_annotations.tsv'), ('variantAnnotations', 'var_drug_ann.tsv'), ('variantAnnotations', 'var_pheno_ann.tsv'),
                          ('variantAnnotations', 'var_fa_ann.tsv')):
            rows = list(tsv_rows(zip_text(self.zp(z), member)))
            self.rows_cache[member] = rows
            for _ln, r, _raw in rows:
                parts = [p.strip() for p in r['Variant/Haplotypes'].split(',') if p.strip()]
                for i in range(len(parts)):
                    for j in range(i + 1, min(len(parts), i + 8) + 1):
                        cand = ', '.join(parts[i:j])
                        if RS.match(cand):
                            rsids.add(cand)
                        else:
                            names.add(cand)
        ctx.want('pgx_allele_name', names)
        ctx.want('pgx_reference_allele')
        ctx.want('dbsnp', rsids)

    # ---- rsID sites -------------------------------------------------------------------------------------------------------------
    def load_sites(self, ctx):
        for _ln, r, _raw in tsv_rows(self.variants_text):
            name, loc = r['Variant Name'], r['Location']
            site = dict(variant_id=r['Variant ID'], location=loc, keys=[], unresolved={}, chrom=None, pos=None, reason=None, status='ambiguous', genes=[g for g in r['Gene Symbols'].split(',') if g.strip()],
                        synonyms=[s.strip() for s in r['Synonyms'].split(',') if s.strip()])
            m = re.match(r'^([^:]+):(\d+)(?:-(\d+))?$', loc or '')
            if not loc:
                site['reason'], site['status'] = 'no_grch38_location_in_clinpgx', 'unmapped'
            elif not m:
                site['reason'], site['status'] = 'location_not_parsed', 'unmapped'
            else:
                acc = m.group(1)
                canon = ctx.ref.contig(acc)
                if canon is None:
                    site['reason'], site['status'] = f'location_on_a_contig_of_another_assembly:{acc}', 'unmapped'
                else:
                    site['chrom'], site['pos'] = canon, int(m.group(2))
                    keys, bad = set(), Counter()
                    for syn in site['synonyms']:
                        if syn.startswith(acc + ':g.'):
                            status, ks, why = hgvs_to_key(ctx.ref, canon, syn[len(acc) + 1:])
                            if status == 'variant':
                                keys.update(ks)
                            elif status == 'unresolved':
                                bad[why.split(':')[0]] += 1
                    site['keys'], site['unresolved'] = sorted(keys), dict(bad)
                    if not keys and bad.get('repeat_length_notation'):
                        site['reason'], site['status'] = 'repeat_length_polymorphism', 'unmapped'
                    elif not keys:
                        site['reason'] = 'alternate_allele_not_given' if not bad else 'alternate_alleles_not_resolvable:' + ','.join(sorted(bad))
            self.sites[name] = site
        st = ctx.st()
        st.inc('rsid_sites', len(self.sites))
        st.inc('rsid_sites_with_keys', sum(1 for s in self.sites.values() if s['keys']))
        st.inc('rsid_sites_multiallelic', sum(1 for s in self.sites.values() if len(s['keys']) > 1))

    # ---- subject resolution -----------------------------------------------------------------------------------------------------
    def resolve_rs(self, ctx, rs, letters):
        site = self.sites.get(rs)
        if site is None:
            hits = sorted(set(ctx.found.get(('dbsnp', rs), ())))
            if len(hits) == 1:
                return Subject('variant', rs, hits[0], 'mapped', None, flags=('located_through_clinvar_rsid',))
            return Subject('variant', rs, status='unmapped', reason='rsid_not_in_clinpgx_variants' if not hits else 'rsid_resolves_to_several_alleles_in_clinvar', candidates=tuple(hits))
        if not site['keys']:
            return Subject('variant', rs, status=site['status'], reason=site['reason'], chrom=site['chrom'], pos=site['pos'])
        status, key, reason, flags = choose_alt(site['keys'], letters)
        if status == 'mapped':
            return Subject('variant', rs, key, 'mapped', None, site['chrom'], site['pos'], flags, tuple(site['keys']))
        return Subject('variant', rs, status=status, reason=reason, chrom=site['chrom'], pos=site['pos'], candidates=tuple(site['keys']))

    def items_of(self, ctx, text):
        """Split a Variant/Haplotypes cell on commas, re-joining fragments that together form an allele name known to the index."""
        parts = [p.strip() for p in (text or '').split(',') if p.strip()]
        out, i = [], 0
        while i < len(parts):
            for j in range(min(len(parts), i + 8), i, -1):
                cand = ', '.join(parts[i:j])
                if j - i > 1 and ctx.found.get(('pgx_allele_name', cand)):
                    out.append(cand)
                    i = j
                    break
            else:
                out.append(parts[i])
                i += 1
        return out

    def resolve_item(self, ctx, item, gene_symbols, letters):
        if RS.match(item):
            return self.resolve_rs(ctx, item, letters)
        hits = ctx.found.get(('pgx_allele_name', item))
        if hits:
            hits = sorted(set(hits), key=lambda h: (h.split(':', 1)[0] != 'cpic', h))
            return Subject('haplotype', item, hits[0], 'mapped', None, also=tuple(hits[1:]))
        if item.startswith('HLA-'):
            return Subject('haplotype', item, status='unmapped', reason='hla_allele_not_defined_in_index_needs_hla_typing')
        if STAR.match(item) or re.match(r'^[A-Za-z0-9-]+ \*', item):
            return Subject('haplotype', item, status='unmapped', reason='allele_not_defined_in_cpic_or_pharmvar')
        for g in gene_symbols:
            if item.startswith(g + ' ') and GENE_DESCRIPTOR.search(item[len(g):]):
                gid, _how = self.genes.resolve(symbol=g)
                if gid:
                    return Subject('gene', item, gid, 'mapped', 'gene_level_descriptor')
                return Subject('gene', item, status='unmapped', reason='gene_symbol_not_in_gene_map')
        if re.match(r'^[cgpm]\.', item):
            return Subject('variant', item, status='unmapped', reason='transcript_or_protein_hgvs_not_resolved')
        if any(item.startswith(g + ' ') for g in gene_symbols):
            return Subject('haplotype', item, status='unmapped', reason='named_polymorphism_without_definition_or_coordinates')
        return Subject('variant', item, status='unmapped', reason='unrecognised_variant_or_haplotype_name')

    def gene_ids(self, symbols):
        out = []
        for g in symbols:
            gid, _how = self.genes.resolve(symbol=g)
            out.append(gid or f'symbol:{g}')
        return out or None

    # ---- evidence helpers -------------------------------------------------------------------------------------------------------
    def subject_kw(self, s):
        if s.status == 'mapped':
            return dict(level=s.level, subject_id=s.entity, mapping_status='mapped', mapping_reason=s.reason, chrom=s.chrom, pos=s.pos)
        return dict(level=s.level, subject_id=None, mapping_status=s.status, mapping_reason=s.reason or 'unmapped', chrom=s.chrom, pos=s.pos, unmapped_key=s.item)

    def emit(self, ctx, rel, acc, locator, raw, s, axis, kind, claim_type, claim_value, effect, **kw):
        k = self.subject_kw(s)
        level, subject, status, reason = k.pop('level'), k.pop('subject_id'), k.pop('mapping_status'), k.pop('mapping_reason')
        extra = dict(kw.pop('extra', {}) or {})
        extra['subject_item'] = s.item
        if s.flags:
            extra['flags'] = list(s.flags)
        if s.candidates:
            extra['candidate_keys'] = list(s.candidates)
        if s.also:
            extra['also_defined_as'] = list(s.also)
        row = evidence_row(rel, 'pharmgkb', acc, locator, raw, level, subject, axis, kind, claim_type, claim_value, effect, status, reason, discriminator=s.item,
                           extra_json=json.dumps(extra, sort_keys=True, ensure_ascii=False, default=str), **k, **kw)
        ctx.emit_evidence(row)
        return row

    # ---- run --------------------------------------------------------------------------------------------------------------------
    def run(self, ctx, limit=None):
        ctx.adapter = self
        st = ctx.st()
        self.load_sites(ctx)
        self.emit_site_links(ctx)
        self.summary_annotations(ctx, st)
        self.variant_annotations(ctx, st)
        self.guidelines(ctx, st)
        self.labels(ctx, st)
        self.crosschecks(ctx, st)
        return st.c

    def emit_site_links(self, ctx):
        rel = ctx.register_release('pharmgkb', 'clinpgx_variants', zip_version(self.zp('variants')))
        for name, site in sorted(self.sites.items()):
            for key in site['keys']:
                ctx.emit('xref', dict(namespace='dbsnp', identifier=name, entity_level='variant', entity_id=key, release_id=rel, relation='clinpgx_variant_name'), name='pharmgkb')
                ctx.emit('xref', dict(namespace='clinpgx_variant', identifier=site['variant_id'], entity_level='variant', entity_id=key, release_id=rel, relation='clinpgx_variant_id'), name='pharmgkb')
                ctx.offer('dbsnp', name, key)
                for g in site['genes']:
                    gid, _how = self.genes.resolve(symbol=g.strip())
                    if gid:
                        ctx.emit('variant_gene', dict(variant_key=key, gene_id=gid, relation='clinpgx_listed_gene', consequence=None, release_id=rel), name='pharmgkb')

    def summary_annotations(self, ctx, st):
        path = self.zp('summaryAnnotations')
        rel = ctx.register_release('pharmgkb', 'clinpgx_summaryAnnotations', zip_version(path))
        alleles, evid = {}, {}
        for _ln, r, _raw in tsv_rows(zip_text(path, 'summary_ann_alleles.tsv')):
            alleles.setdefault(r['Summary Annotation ID'], []).append({k: r[k] for k in r if k != 'Summary Annotation ID'})
        for _ln, r, _raw in tsv_rows(zip_text(path, 'summary_ann_evidence.tsv')):
            evid.setdefault(r['Summary Annotation ID'], []).append({k: r[k] for k in r if k != 'Summary Annotation ID'})
        for ln, r, raw in self.rows_cache['summary_annotations.tsv']:
            st.inc('summary_annotations')
            genes = [g.strip() for g in re.split(r'[;,]', r['Gene']) if g.strip()]
            text_alleles = ' '.join(a.get('Genotype/Allele', '') for a in alleles.get(r['Summary Annotation ID'], ()))
            letters = letters_of(text_alleles)
            level = r['Level of Evidence']
            drugs = [d for d in r['Drug(s)'].split(';') if d]
            extra = dict(genes=genes, drugs=drugs, phenotypes=r['Phenotype(s)'], phenotype_category=r['Phenotype Category'], level_override=r['Level Override'], level_modifiers=r['Level Modifiers'],
                         pmid_count=r['PMID Count'], evidence_count=r['Evidence Count'], latest_history_date=r['Latest History Date (YYYY-MM-DD)'], url=r['URL'], specialty_population=r['Specialty Population'],
                         genotype_annotations=alleles.get(r['Summary Annotation ID'], []), evidence=evid.get(r['Summary Annotation ID'], []))
            for s in self.of_interest(ctx, self.subjects(ctx, r['Variant/Haplotypes'], genes, letters), text_alleles, genes):
                self.emit(ctx, rel, f"summary_annotation:{r['Summary Annotation ID']}", f'summary_annotations.tsv#L{ln}', raw, s, 'association_replication', 'aggregated_integration',
                          'clinpgx_level_of_evidence', level, LEVEL_EFFECT.get(level, 'uncertain'), claim_num=num(r['Score']), score=num(r['Score']), score_unit='clinpgx_annotation_score',
                          condition_text='; '.join(x for x in (', '.join(drugs), r['Phenotype(s)']) if x) or None, gene_ids=self.gene_ids(genes), curator='ClinPGx',
                          review_status=f'Level {level}' + (f' (override {r["Level Override"]})' if r['Level Override'] else ''), origin='germline',
                          lineage_group=f"clinpgx_summary:{r['Summary Annotation ID']}", lineage_note='ClinPGx summary of its variant annotations; not an independent study', extra=extra)
                st.inc(f'summary_rows:{s.status}')

    def subjects(self, ctx, text, genes, letters):
        out, seen = [], set()
        for item in self.items_of(ctx, text):
            if item in seen:
                continue
            seen.add(item)
            out.append(self.resolve_item(ctx, item, genes, letters))
        return out

    def variant_annotations(self, ctx, st):
        path = self.zp('variantAnnotations')
        rel = ctx.register_release('pharmgkb', 'clinpgx_variantAnnotations', zip_version(path))
        params = {}
        for ln, r, raw in tsv_rows(zip_text(path, 'study_parameters.tsv')):
            params.setdefault(r['Variant Annotation ID'], []).append((ln, r, raw))
        for member, kind_label in (('var_drug_ann.tsv', 'drug'), ('var_pheno_ann.tsv', 'phenotype'), ('var_fa_ann.tsv', 'functional')):
            for ln, r, raw in self.rows_cache[member]:
                self.one_variant_annotation(ctx, st, rel, member, kind_label, ln, r, raw, params.get(r['Variant Annotation ID'], []))

    def one_variant_annotation(self, ctx, st, rel, member, kind_label, ln, r, raw, params):
        st.inc(f'variant_annotations:{kind_label}')
        genes = [g.strip() for g in re.split(r'[;,]', r['Gene']) if g.strip()]
        letters = letters_of(r.get('Alleles'), r.get('Comparison Allele(s) or Genotype(s)'))
        vid, pmid = r['Variant Annotation ID'], r['PMID']
        functional = kind_label == 'functional'
        assoc = (r.get('Is/Is Not associated') or '').strip()
        sig = (r.get('Significance') or '').strip().lower()
        effect = 'no_effect' if assoc.startswith('Not') else {'yes': 'effect', 'no': 'no_effect'}.get(sig, 'uncertain')
        study_types = sorted({p[1]['Study Type'] for p in params if p[1]['Study Type']})
        ancestry = '; '.join(sorted({p[1]['Biogeographical Groups'] for p in params if p[1]['Biogeographical Groups'] and p[1]['Biogeographical Groups'] != 'Unknown'})) or None
        drugs = r.get('Drug(s)') or ''
        pheno = r.get('Phenotype') or r.get('Functional terms') or r.get('Side effect/efficacy/other') or ''
        n_cases = sum(int(num(p[1]['Study Cases'])) for p in params if num(p[1]['Study Cases']) is not None) if len(params) == 1 else None
        n_controls = sum(int(num(p[1]['Study Controls'])) for p in params if num(p[1]['Study Controls']) is not None) if len(params) == 1 else None
        ctx.emit('entity_study', dict(study_id=f'clinpgx_va:{vid}', source_id='pharmgkb', lineage_group=f'pmid:{pmid}', study_type=(study_types[0] if study_types else ('functional_assay' if functional else None)),
                                      pubmed=pmid, trait_ids=None, trait_text='; '.join(x for x in (drugs, pheno) if x) or None, n_total=(n_cases or 0) + (n_controls or 0) or None, n_cases=n_cases,
                                      ancestry=ancestry, cohorts=None, biosample=r.get('Cell type') or None, has_sumstats=False, release_id=rel), name='pharmgkb')
        group = {'alleles': r.get('Alleles'), 'comparison': r.get('Comparison Allele(s) or Genotype(s)'), 'metabolizer_types': r.get('Metabolizer types'),
                 'comparison_metabolizer_types': r.get('Comparison Metabolizer types'), 'direction_of_effect': r.get('Direction of effect'), 'notes': (r.get('Notes') or '')[:1000],
                 'phenotype_category': r.get('Phenotype Category'), 'population_types': r.get('Population types'), 'specialty_population': r.get('Specialty Population'),
                 'pd_pk_terms': r.get('PD/PK terms'), 'assay_type': r.get('Assay type'), 'significance': r.get('Significance'), 'is_is_not_associated': assoc,
                 'study_types': study_types}
        subjects = self.subjects(ctx, r['Variant/Haplotypes'], genes, letters)
        subjects = self.of_interest(ctx, subjects, r.get('Alleles'), genes)
        axis = 'empirical_functional_support' if functional else 'association_replication'
        kind = 'direct_experimental' if functional else ('case_observation' if study_types and all('case series' in t or 'case report' in t for t in study_types) else 'statistical_association')
        claim = (r.get('Sentence') or '')[:600]
        for s in subjects:
            row_kw = dict(condition_text='; '.join(x for x in (drugs, pheno) if x) or None, gene_ids=self.gene_ids(genes), curator='ClinPGx', ancestry=ancestry, sample_size=(n_cases or 0) + (n_controls or 0) or None,
                          n_cases=n_cases, origin='germline', lineage_group=f'pmid:{pmid}', lineage_note='a PubMed-indexed study; the same paper may be cited by other databases', model_system=r.get('Assay type') or None,
                          cell_type=r.get('Cell type') or None)
            self.emit(ctx, rel, f'variant_annotation:{vid}', f'{member}#L{ln}', raw, s, axis, kind, 'pharmacogenomic_functional_assay' if functional else 'pharmacogenomic_association', claim, effect,
                      extra=dict(group, studies=[{k: v for k, v in p[1].items() if k not in ('Study Parameters ID', 'Variant Annotation ID')} for p in params]), **row_kw)
            st.inc(f'variant_annotation_rows:{s.status}')
            for ln2, p, raw2 in params:
                ratio = num(p['Ratio Stat'])
                if ratio is None or p['Ratio Stat Type'] not in ('OR', 'HR', 'RR'):
                    continue
                pv, ptxt = parse_p(p['P Value'])
                cases, controls = num(p['Study Cases']), num(p['Study Controls'])
                ci_lo, ci_hi, ci_note = num(p['Confidence Interval Start']), num(p['Confidence Interval Stop']), {}
                if ci_lo is not None and ci_hi is not None and ci_lo > ci_hi:            # a start above the stop is a typing error in the source, not an interval: kept as text, not stored as bounds
                    ci_note = dict(interval_start_reported=p['Confidence Interval Start'], interval_stop_reported=p['Confidence Interval Stop'], interval_not_stored='start above stop')
                    ci_lo = ci_hi = None
                    st.inc('ratio_rows_with_inverted_interval')
                eff_kw = dict(odds_ratio=ratio if p['Ratio Stat Type'] == 'OR' else None, ci_low=ci_lo, ci_high=ci_hi, p_value=pv,
                              score=None if p['Ratio Stat Type'] == 'OR' else ratio, score_unit=None if p['Ratio Stat Type'] == 'OR' else p['Ratio Stat Type'],
                              ancestry=p['Biogeographical Groups'] if p['Biogeographical Groups'] not in ('', 'Unknown') else None, curator='ClinPGx', gene_ids=self.gene_ids(genes),
                              sample_size=int((cases or 0) + (controls or 0)) or None, n_cases=int(cases) if cases is not None else None, claim_num=ratio, origin='germline',
                              condition_text=row_kw['condition_text'], lineage_group=f'pmid:{pmid}', lineage_note='a PubMed-indexed study; the same paper may be cited by other databases')
                self.emit(ctx, rel, f"variant_annotation:{vid}:params:{p['Study Parameters ID']}", f'study_parameters.tsv#L{ln2}', raw2, s, 'effect_magnitude_direction', 'statistical_association',
                          'ratio_statistic:' + p['Ratio Stat Type'], f"{p['Ratio Stat Type']} {ratio:g}", effect,
                          extra=dict(p_value_text=ptxt, study_type=p['Study Type'], characteristics=p['Characteristics'], characteristics_type=p['Characteristics Type'],
                                     frequency_in_cases=p['Frequency In Cases'], frequency_in_controls=p['Frequency In Controls'], allele_of_frequency_cases=p['Allele Of Frequency In Cases'],
                                     allele_of_frequency_controls=p['Allele Of Frequency In Controls'], study_controls=p['Study Controls'], comparison=group['comparison'], alleles=group['alleles'],
                                     reference_unit='genotype_comparison_not_alt_allele', **ci_note), **eff_kw)
                st.inc('ratio_rows')

    def of_interest(self, ctx, subjects, alleles_text, genes):
        """The subjects a finding is about. For star alleles these are the alleles named in the 'Alleles' cell; reference alleles are kept only when
        nothing else remains (the finding is then about the reference allele)."""
        if not subjects or not any(s.level == 'haplotype' for s in subjects):
            return subjects
        named = set(re.findall(r'\*[0-9A-Za-z:.\-]+', alleles_text or ''))
        keep = [s for s in subjects if s.level != 'haplotype' or not named or any(s.item.endswith(n) for n in named)]
        non_ref = [s for s in keep if not (s.entity and s.status == 'mapped' and ctx.found.get(('pgx_reference_allele', s.entity)))]
        return non_ref or keep or subjects

    def guidelines(self, ctx, st):
        path = self.zp('guidelineAnnotations.json')
        rel = ctx.register_release('pharmgkb', 'clinpgx_guidelineAnnotations', zip_version(path))
        with zipfile.ZipFile(path) as z:
            for name in sorted(n for n in z.namelist() if n.endswith('.json')):
                raw = z.read(name).decode('utf-8', errors='replace')
                g = json.loads(raw).get('guideline') or {}
                gid, src = g.get('id'), g.get('source')
                if not gid:
                    st.inc('guideline_files_without_id')
                    continue
                st.inc('guideline_annotations')
                genes = [x['symbol'] for x in g.get('relatedGenes', ()) if x.get('symbol')]
                chems = [x['name'] for x in g.get('relatedChemicals', ()) if x.get('name')]
                flags = dict(recommendation=g.get('recommendation'), dosing_information=g.get('dosingInformation'), has_testing_info=g.get('hasTestingInfo'),
                             alternate_drug_available=g.get('alternateDrugAvailable'), other_prescribing_guidance=g.get('otherPrescribingGuidance'), cancer_genome=g.get('cancerGenome'),
                             pediatric=g.get('pediatric'), related_alleles=[x.get('symbol') for x in g.get('relatedAlleles', ())], chemicals=chems)
                value = f"{src}: " + ('recommendation provided' if g.get('recommendation') else 'no recommendation')
                common = dict(condition_text=', '.join(chems) or None, curator=src, review_status=src, lineage_group=f'clinpgx_guideline:{gid}', gene_ids=self.gene_ids(genes),
                              lineage_note='guideline text of the named body; the same guideline can also be read from the body\'s own publication (for CPIC, the CPIC API rows)')
                for sym in genes:
                    gene_id, _how = self.genes.resolve(symbol=sym)
                    s = Subject('gene', sym, gene_id, 'mapped' if gene_id else 'unmapped', None if gene_id else 'gene_symbol_not_in_gene_map')
                    self.emit(ctx, rel, f'guideline_annotation:{gid}', name, raw, s, 'clinical_validity_actionability', 'curated_expert_classification', 'pgx_guideline_annotation', value, 'not_applicable',
                              extra=dict(flags, name=g.get('name'), version=g.get('version')), **common)
                    st.inc('guideline_gene_rows')
                for a in g.get('relatedAlleles', ()):
                    s = self.resolve_item(ctx, a.get('symbol') or a.get('name') or '', genes, set())
                    if s.status == 'mapped' and s.level == 'gene':
                        continue
                    self.emit(ctx, rel, f'guideline_annotation:{gid}', name, raw, s, 'clinical_validity_actionability', 'curated_expert_classification', 'pgx_guideline_annotation', value, 'not_applicable',
                              extra=dict(flags, name=g.get('name'), allele_object=a.get('objCls')), **common)
                    st.inc(f'guideline_allele_rows:{s.status}')

    def labels(self, ctx, st):
        path = self.zp('drugLabels')
        rel = ctx.register_release('pharmgkb', 'clinpgx_drugLabels', zip_version(path))
        for ln, r, raw in tsv_rows(zip_text(path, 'drugLabels.tsv')):
            st.inc('drug_labels')
            genes = [g.strip() for g in r['Genes'].split(';') if g.strip()]
            value = r['Testing Level'] or 'no testing level stated'
            common = dict(condition_text=r['Chemicals'] or None, curator=r['Source'] or None, review_status=r['Source'] or None, lineage_group=f"drug_label:{r['PharmGKB ID']}", gene_ids=self.gene_ids(genes),
                          lineage_note='regulator label text annotated by ClinPGx')
            extra = {k: r[k] for k in ('Name', 'Source', 'Biomarker Flag', 'Has Prescribing Info', 'Has Dosing Info', 'Has Alternate Drug', 'Has Other Prescribing Guidance', 'Cancer Genome', 'Prescribing',
                                       'Variants/Haplotypes', 'Latest History Date (YYYY-MM-DD)')}
            for sym in genes:
                gene_id, _how = self.genes.resolve(symbol=sym)
                s = Subject('gene', sym, gene_id, 'mapped' if gene_id else 'unmapped', None if gene_id else 'gene_symbol_not_in_gene_map')
                self.emit(ctx, rel, f"drug_label:{r['PharmGKB ID']}", f'drugLabels.tsv#L{ln}', raw, s, 'clinical_validity_actionability', 'curated_expert_classification', 'drug_label_testing_level', value,
                          'not_applicable', extra=extra, **common)
                st.inc('drug_label_gene_rows')

    def crosschecks(self, ctx, st):
        """Row-level comparisons that integrate nothing: clinicalVariants must agree with the summary annotations, and the older clinical-annotation
        export is compared with the newer summary export so that release-to-release changes are counted (the change log reads these counters)."""
        summ = {r['Summary Annotation ID']: r for _ln, r, _raw in self.rows_cache['summary_annotations.tsv']}
        legacy_path = self.zp('clinicalAnnotations')
        if os.path.exists(legacy_path):
            legacy = {r['Clinical Annotation ID']: r for _ln, r, _raw in tsv_rows(zip_text(legacy_path, 'clinical_annotations.tsv'))}
            st.inc('legacy_clinical_annotations', len(legacy))
            st.inc('legacy_ids_not_in_summary', len(set(legacy) - set(summ)))
            st.inc('summary_ids_not_in_legacy', len(set(summ) - set(legacy)))
            for k in set(legacy) & set(summ):
                for f in ('Level of Evidence', 'Score', 'Drug(s)', 'Phenotype(s)', 'Variant/Haplotypes'):
                    if legacy[k][f] != summ[k][f]:
                        st.inc(f'legacy_vs_summary_changed:{f}')
            ctx.problems.append(f'clinpgx: clinicalAnnotations export is {zip_version(legacy_path)}, summaryAnnotations is {zip_version(self.zp("summaryAnnotations"))}; only the newer export is integrated')
        cv_path = self.zp('clinicalVariants')
        if os.path.exists(cv_path):
            cv = list(tsv_rows(zip_text(cv_path, 'clinicalVariants.tsv')))
            genes = lambda t: frozenset(g.strip() for g in re.split(r'[;,]', t) if g.strip())       # summary files separate genes with ';', clinicalVariants with ','
            have = {(r['Variant/Haplotypes'], genes(r['Gene']), r['Level of Evidence']) for r in summ.values()}
            st.inc('clinical_variants_rows', len(cv))
            st.inc('clinical_variants_without_matching_summary', sum(1 for _ln, r, _raw in cv if (r['variant'], genes(r['gene']), r['level of evidence']) not in have))
