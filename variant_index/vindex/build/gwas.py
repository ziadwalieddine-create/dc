"""GWAS Catalog associations (release 2026-10-22) -> evidence rows, with studies and ancestry as study entities.

Each association becomes up to two rows that share one lineage group (they are one measurement): the association itself (p-value, sample
metadata; axis association_replication) and its effect size (axis effect_magnitude_direction: odds ratio or beta, units, interval, effect
allele frequency). Rows are placed on a variant key only when the reported effect allele can be oriented against the GRCh38 reference
and the alleles known at the site (alleles.resolve_snv): wrong-strand and palindromic records without a frequency that settles them are
rejected, not guessed; records with an unknown effect allele map to the site's single alternate allele with no direction; haplotype,
interaction, HLA/amino-acid, indel-notation and multi-base records are kept as unmapped rows with the reason. The GWAS Catalog has no
association accession in this file, so the accession is study + SNP string + trait and the locator is the line number.
"""
import csv
import io
import json
import os
import re
import zipfile

from .alleles import resolve_snv
from .common import DATA, evidence_row

GW = os.path.join(DATA, 'gwas_catalog')
RS_ALLELE = re.compile(r'^(rs\d+)-([ACGT]+|\?)$')
CHR_ALLELE = re.compile(r'^chr([0-9XYM]+):(\d+)-([ACGT?]+)$')
GENOME_WIDE = 5e-8


def num(x):
    try:
        return float(x)
    except (TypeError, ValueError):
        return None


def efo_ids(uris):
    out = []
    for u in (uris or '').replace(';', ',').split(','):
        u = u.strip().rsplit('/', 1)[-1]
        if u:
            out.append(u.replace('_', ':', 1))
    return out


def parse_effect(or_beta, ci_text):
    """('beta' | 'odds_ratio' | None, signed value, ci_low, ci_high, unit text, direction). Beta is recognised by an 'increase'/'decrease' unit in
    the 95% CI text, as the Catalog writes it; the number for a decrease is stored negative; an odds ratio has no unit text."""
    v = num(or_beta)
    if v is None:
        return None, None, None, None, None, None
    t = (ci_text or '').strip()
    m = re.match(r'^\[\s*([-+]?[0-9.eE+-]+?)\s*-\s*([-+]?[0-9.eE+-]+)\s*\]\s*(.*)$', t)
    lo = hi = None
    unit = t
    if m:
        lo, hi, unit = num(m.group(1)), num(m.group(2)), m.group(3).strip()
    elif t.startswith('[NR]'):
        unit = t[4:].strip()
    low = unit.lower()
    if 'decrease' in low or 'increase' in low:
        sign = -1.0 if 'decrease' in low else 1.0
        return 'beta', sign * abs(v), lo, hi, unit, 'decreases' if sign < 0 else 'increases'
    if v > 0:
        return 'odds_ratio', v, lo, hi, unit or None, 'increases' if v > 1 else ('decreases' if v < 1 else None)
    return None, v, lo, hi, unit or None, None


class GwasCatalogAdapter:
    source_id = 'gwas_catalog'
    targets = ['gwas_catalog_associations', 'gwas_catalog_studies', 'gwas_catalog_ancestries']

    def __init__(self, genes, variants, root=None, assoc_zip=None):
        self.genes, self.variants = genes, variants
        self.root = root or GW
        self.assoc_zip = assoc_zip or os.path.join(self.root, 'gwas-catalog-associations_ontology-annotated-full.zip')

    # ---- study metadata ---------------------------------------------------------------------------------------------------------
    def studies(self, ctx, rel_studies, rel_anc, st):
        study = {}
        with open(os.path.join(self.root, 'gwas-catalog-download-studies-v1.0.3.1.tsv'), encoding='utf-8', errors='replace', newline='') as f:
            for r in csv.DictReader(f, delimiter='\t', quoting=csv.QUOTE_NONE):
                study[r['STUDY ACCESSION']] = dict(pubmed=r['PUBMED ID'], trait=r['DISEASE/TRAIT'], mapped_trait=r['MAPPED_TRAIT'], uris=r['MAPPED_TRAIT_URI'], cohort=r['COHORT'],
                                                   full=r['FULL SUMMARY STATISTICS'].lower() == 'yes', initial=r['INITIAL SAMPLE SIZE'], n_init=0, n_cases=0, anc=set(), anc_all=set())
        with open(os.path.join(self.root, 'gwas-catalog-download-ancestries-v1.0.3.1.tsv'), encoding='utf-8', errors='replace', newline='') as f:
            for r in csv.DictReader(f, delimiter='\t', quoting=csv.QUOTE_NONE):
                s = study.get(r['STUDY ACCESSION'])
                if s is None:
                    st.inc('ancestry_rows_without_study')
                    continue
                cat = r['BROAD ANCESTRAL CATEGORY'].strip()
                s['anc_all'].add(cat)
                if r['STAGE'] == 'initial':
                    s['anc'].add(cat)
                    n = num(r['NUMBER OF INDIVIDUALS'])
                    if n:
                        s['n_init'] += int(n)
                    c = num(r['NUMBER OF CASES'])
                    if c:
                        s['n_cases'] += int(c)
        for acc, s in sorted(study.items()):
            ctx.emit('entity_study', dict(study_id=acc, source_id='gwas_catalog', lineage_group=f'gwas_study:{acc}', study_type='gwas', pubmed=s['pubmed'], trait_ids=efo_ids(s['uris']),
                                          trait_text=s['trait'], n_total=s['n_init'] or None, n_cases=s['n_cases'] or None, ancestry='; '.join(sorted(s['anc'])), cohorts=s['cohort'],
                                          biosample=None, has_sumstats=s['full'], release_id=rel_studies), name='gwas_catalog')
            st.inc('studies')
        return study

    # ---- associations -----------------------------------------------------------------------------------------------------------
    def rows(self):
        with zipfile.ZipFile(self.assoc_zip) as z:
            name = z.namelist()[0]
            with z.open(name) as raw:
                f = io.TextIOWrapper(raw, encoding='utf-8', errors='replace', newline='')
                rd = csv.reader(f, delimiter='\t', quoting=csv.QUOTE_NONE)
                hdr = next(rd)
                for ln, row in enumerate(rd, start=2):
                    if len(row) < len(hdr):
                        row = row + [''] * (len(hdr) - len(row))
                    yield ln, dict(zip(hdr, row)), '\t'.join(row)

    def classify(self, r):
        """(kind, chrom, pos, rsid, risk) with kind in rs | chr | multi | interaction | other."""
        a = r['STRONGEST SNP-RISK ALLELE'].strip()
        if ' x ' in a:
            return 'interaction', None, None, None, None
        if ';' in a or ',' in a:
            return 'multi', None, None, None, None
        m = RS_ALLELE.match(a)
        if m:
            rs = r['SNP_ID_CURRENT'] and 'rs' + r['SNP_ID_CURRENT'] if r['MERGED'] == '1' and r['SNP_ID_CURRENT'] else m.group(1)
            chrom, pos = r['CHR_ID'].strip(), r['CHR_POS'].strip()
            if re.fullmatch(r'[0-9XYM]+', chrom) and pos.isdigit():
                return 'rs', chrom.replace('M', 'MT') if chrom == 'M' else chrom, int(pos), rs, m.group(2)
            return 'rs_nopos', None, None, rs, m.group(2)
        m = CHR_ALLELE.match(a)
        if m:
            return 'chr', m.group(1).replace('M', 'MT') if m.group(1) == 'M' else m.group(1), int(m.group(2)), None, m.group(3)
        if re.match(r'^chr[0-9XYM]+:\d+:[ID]-', a):
            return 'chr_indel', None, None, None, None
        return 'other', None, None, None, None

    def run(self, ctx, limit=None):
        ctx.adapter = self
        st = ctx.st()
        rel = ctx.register_release('gwas_catalog', 'gwas_catalog_associations', '2026-10-22')
        rel_s = ctx.register_release('gwas_catalog', 'gwas_catalog_studies', '2026-10-22')
        rel_a = ctx.register_release('gwas_catalog', 'gwas_catalog_ancestries', '2026-10-22')
        study = self.studies(ctx, rel_s, rel_a, st)
        # pass 1: positions the allele source must answer for
        need = set()
        pre = []
        for ln, r, raw in self.rows():
            if limit and ln - 1 > limit:
                break
            k = self.classify(r)
            pre.append(k)
            if k[0] in ('rs', 'chr') and ctx.ref.contig(k[1]):
                need.add((k[1], k[2]))
        self.variants.load(need)
        # pass 1b: frequencies for the sites whose strand needs them
        want = set()
        for k, (ln, r, raw) in zip(pre, self.rows()):
            if k[0] in ('rs', 'chr') and k[4] and k[4] not in ('?',) and len(k[4]) == 1:
                canon = ctx.ref.contig(k[1])
                if not canon:
                    continue
                cands = self.cands(k, canon, ctx)
                res = resolve_snv(self.refbase(ctx, canon, k[2]), k[4], num(r['RISK ALLELE FREQUENCY']), {c[0] for c in cands}, None)
                if res.status == 'rejected' and num(r['RISK ALLELE FREQUENCY']) is not None and len(cands) == 1:
                    want.add(cands[0][2])
            if limit and ln - 1 >= limit:
                break
        afs = self.variants.load_frequencies(want) if want else {}
        # pass 2: evidence
        for k, (ln, r, raw) in zip(pre, self.rows()):
            st.inc('records')
            self.emit(ctx, rel, study, ln, r, raw, k, afs, st)
        return st.c

    def refbase(self, ctx, canon, pos):
        try:
            return ctx.ref.fetch(canon, pos, pos)
        except ValueError:
            return 'N'

    def cands(self, k, canon, ctx):
        base = self.refbase(ctx, canon, k[2])
        c = self.variants.snv_alts(k[1], k[2], base)
        if k[3]:
            with_rs = [x for x in c if k[3] in x[1]]
            if with_rs:
                return with_rs
        return c

    def emit(self, ctx, rel, study, ln, r, raw, k, afs, st):
        kind, chrom, pos, rs, risk = k
        s = study.get(r['STUDY ACCESSION'], {})
        raf = num(r['RISK ALLELE FREQUENCY'])
        p = num(r['P-VALUE'])
        common = dict(lineage_group=f"gwas_study:{r['STUDY ACCESSION']}", lineage_note=f"GWAS Catalog study; PMID {r['PUBMEDID']}; cohort {s.get('cohort') or 'not stated'}; other databases that derive from this study are not independent replication",
                      condition_ids=efo_ids(r['MAPPED_TRAIT_URI']), condition_text=r['DISEASE/TRAIT'], gene_ids=[g.strip() for g in re.split(r'[;,]', r['SNP_GENE_IDS']) if g.strip().startswith('ENSG')] or None,
                      sample_size=s.get('n_init') or None, n_cases=s.get('n_cases') or None, ancestry='; '.join(sorted(s.get('anc', ()))) or None, origin='germline')
        subject, mapping, reason, flags, eff_allele, other_allele = None, 'unmapped', None, (), None, None
        unmapped_key = r['SNPS'] or r['STRONGEST SNP-RISK ALLELE']
        canon = ctx.ref.contig(chrom) if chrom else None
        if canon and pos is not None and not 1 <= pos <= ctx.ref.length(canon):
            st.inc('positions_beyond_contig_length')              # a typo in the catalogue (one real record has chr2:179179368916): kept as unmapped with its text, never coerced
            kind, pos = 'bad_position', None
        if kind in ('rs', 'chr'):
            if canon is None:
                reason = 'contig_unmapped'
            else:
                base = self.refbase(ctx, canon, pos)
                cands = self.cands(k, canon, ctx)
                if kind == 'chr' and not cands:
                    reason = 'author_position_assembly_unverified'
                elif len(risk) > 1 and risk != '?':
                    reason = 'multi_base_risk_allele'
                else:
                    alt_afs = {c[0]: afs.get(c[2]) for c in cands}
                    res = resolve_snv(base, risk, raf, {c[0] for c in cands}, alt_afs)
                    flags = res.flags
                    if res.status == 'mapped':
                        subject, mapping = f'{canon}:{pos}:{base}:{res.alt}', 'mapped'
                        eff_allele, other_allele = res.effect, res.other
                        if kind == 'rs' and cands and not any(rs in c[1] for c in cands):
                            flags = flags + ('rsid_not_listed_at_position',)
                    else:
                        mapping = {'ambiguous': 'ambiguous', 'rejected': 'rejected'}.get(res.status, 'unmapped')
                        reason = res.reason
        else:
            reason = {'interaction': 'interaction_between_variants', 'multi': 'multi_variant_record', 'rs_nopos': 'no_position_in_catalog', 'chr_indel': 'indel_notation_in_author_identifier',
                      'other': 'non_rsid_identifier_or_hla_or_amino_acid', 'bad_position': 'position_beyond_contig_length'}[kind]
        c0 = dict(chrom=canon if (subject or kind in ('rs', 'chr')) else None, pos=pos if (subject or kind in ('rs', 'chr')) else None, unmapped_key=None if subject else unmapped_key)
        acc = f"{r['STUDY ACCESSION']}:{unmapped_key}:{r['MAPPED_TRAIT'] or r['DISEASE/TRAIT']}"
        ev_assert = 'effect' if (p is not None and p <= GENOME_WIDE) else 'uncertain'
        extra = dict(snps=r['SNPS'], risk_allele_reported=r['STRONGEST SNP-RISK ALLELE'], context=r['CONTEXT'], flags=list(flags), p_text=r['P-VALUE (TEXT)'], pvalue_mlog=r['PVALUE_MLOG'],
                     replication=r['REPLICATION SAMPLE SIZE'], initial=r['INITIAL SAMPLE SIZE'], platform=r['PLATFORM [SNPS PASSING QC]'], pubmed=r['PUBMEDID'])
        if p == 0:
            extra['p_reported_as_zero'] = True
        ctx.emit_evidence(evidence_row(rel, 'gwas_catalog', acc, f'L{ln}', raw, 'variant', subject, 'association_replication', 'statistical_association', 'gwas_association',
                                       r['MAPPED_TRAIT'] or r['DISEASE/TRAIT'], ev_assert, mapping, reason, discriminator=f'L{ln}', p_value=p, eaf=raf, effect_allele=eff_allele, other_allele=other_allele,
                                       extra_json=json.dumps(extra, sort_keys=True), **c0, **common))
        kind_e, val, lo, hi, unit, direction = parse_effect(r['OR or BETA'], r['95% CI (TEXT)'])
        if kind_e:
            effect_known = eff_allele is not None
            ctx.emit_evidence(evidence_row(rel, 'gwas_catalog', acc, f'L{ln}', raw, 'variant', subject, 'effect_magnitude_direction', 'statistical_association', 'gwas_effect_size',
                                           kind_e, ev_assert, mapping, reason, discriminator=f'L{ln}:effect', p_value=p, eaf=raf, effect_allele=eff_allele, other_allele=other_allele,
                                           beta=val if kind_e == 'beta' else None, odds_ratio=val if kind_e == 'odds_ratio' else None, ci_low=lo, ci_high=hi, score_unit=unit,
                                           direction=direction if effect_known else None, extra_json=json.dumps(dict(raw_value=r['OR or BETA'], ci_text=r['95% CI (TEXT)'], direction_reported=direction,
                                           direction_attributed_to_allele=effect_known), sort_keys=True), **c0, **common))
        if subject:
            for g in common['gene_ids'] or ():
                ctx.emit('variant_gene', dict(variant_key=subject, gene_id=g, relation='gwas_catalog_snp_gene', consequence=r['CONTEXT'].split(';')[0].strip(), release_id=rel), name='gwas_catalog')
            if rs:
                ctx.emit('xref', dict(namespace='dbsnp', identifier=rs, entity_level='variant', entity_id=subject, release_id=rel, relation='gwas_catalog_snps'), name='gwas_catalog')
