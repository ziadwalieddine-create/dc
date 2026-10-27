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
import math
import os
import re
import zipfile

from .alleles import resolve_snv
from .common import DATA, evidence_row

GW = os.path.join(DATA, 'gwas_catalog')
RS_ALLELE = re.compile(r'^(rs\d+)-([ACGT]+|\?)$')
CHR_ALLELE = re.compile(r'^chr([0-9XYM]+):(\d+)-([ACGT?]+)$')
GENOME_WIDE = 5e-8
NOT_REPORTED = {'', 'NR', 'NA', 'N/A', '.', '-'}               # the Catalog writes 'NR' where a value was not reported: that is the absence of a value, not a value


def num(x):
    """A finite number, or None. float() also reads 'nan', 'inf' and 'infinity': the real file has such cells, and they are not measurements."""
    try:
        v = float(x)
    except (TypeError, ValueError):
        return None
    return v if math.isfinite(v) else None


def effect_type_basis(kind, ci_text):
    """Why a value is called a beta or an odds ratio: the Catalog writes the unit of a beta in the interval text and writes no unit for an odds ratio. A positive value with no interval
    text at all (12,039 real rows) is assumed to be an odds ratio, and the row says so."""
    if kind == 'beta':
        return 'unit_in_interval_text'
    return 'interval_text_without_unit' if (ci_text or '').strip() else 'assumed_odds_ratio_no_interval_text'


def lineage_of(pubmed, accession):
    """One lineage per publication. Open Targets, ClinPGx and the gene-evidence adapters key a paper as pmid:<id>; the Catalog's own study accession would make the same paper's
    associations and the fine-mapping built on them (98.7% of the Open Targets GWAS publications are in the Catalog) count as independent replication. Several Catalog studies of one
    paper count once, which can understate independent samples but never overstates them. A study without a PubMed id keeps its accession."""
    p = (pubmed or '').strip()
    return f'pmid:{p}' if p.isdigit() else f'gwas_study:{accession}'


def frequency(x):
    """An allele frequency: a finite number in [0, 1], else None. (The real file has 74 cells above 1.)"""
    v = num(x)
    return v if v is not None and 0 <= v <= 1 else None


def signed_interval(lo, hi, sign):
    """The interval of a beta that contains the signed beta. The Catalog writes the bounds of the magnitude and puts 'increase' or 'decrease' in the text, so a decrease of 0.032 reads
    [0.018-0.046]; some entries write the bounds already signed, in either order ([-0.05--0.087]). Bounds that straddle zero are taken as they are, sorted."""
    if lo is None or hi is None:
        return lo, hi
    a, b = min(lo, hi), max(lo, hi)
    return (-b, -a) if sign < 0 and a >= 0 else (a, b)


def efo_ids(uris):
    out = []
    for u in (uris or '').replace(';', ',').split(','):
        u = u.strip().rsplit('/', 1)[-1]
        if u:
            out.append(u.replace('_', ':', 1))
    return out


def text_or_none(x):
    """The Catalog's text, with its not-reported markers (blank, 'NR', 'NA') as None: a placeholder is not a trait, a cohort or a unit."""
    t = (x or '').strip()
    return None if t.upper() in NOT_REPORTED else t


NUMBER = r'[-+\u2212]?(?:\d+\.?\d*|\.\d+)(?:[eE][-+]?\d+)?'
# The interval as the Catalog writes it: normally '[lo-hi] unit', but also '[lo, hi]', '[lo\u2013hi]', 'lo-hi unit' and a closing brace; an upper bound of infinity is not a number
BRACKETED = re.compile(r'^\[\s*([-+]?[0-9.eE+-]+?)\s*-\s*([-+]?[0-9.eE+-]+)\s*\]\s*(.*)$')
LOOSE = re.compile(rf'^[\[(]?\s*({NUMBER})\s*(?:-|\u2013|\u2014|\u2212|,|;)\s*({NUMBER}|\u221e|inf(?:inity)?)\s*[\])}}]?\s*(.*)$', re.I)


def parse_effect(or_beta, ci_text):
    """('beta' | 'odds_ratio' | None, signed value, ci_low, ci_high, unit, direction). Beta is recognised by an 'increase'/'decrease' unit in the 95% CI text, as the Catalog
    writes it; the number for a decrease is stored negative and the direction word is taken out of the unit (it is the direction column); an odds ratio has no unit text."""
    v = num(or_beta)
    if v is None:
        return None, None, None, None, None, None
    t = (ci_text or '').strip()
    m = BRACKETED.match(t) or LOOSE.match(t)
    lo = hi = None
    unit = t
    if m:
        lo, hi, unit = num(m.group(1).replace('\u2212', '-')), num(m.group(2).replace('\u2212', '-')), m.group(3).strip()
    unit = re.sub(r'^\[?NR\]?(?=\s|$)\s*', '', unit).strip()          # 'NR' stands for the interval that was not reported; the unit that follows it is kept
    unit = None if unit.upper() in NOT_REPORTED else unit
    low = (unit or '').lower()
    if 'decrease' in low or 'increase' in low:
        sign = -1.0 if 'decrease' in low else 1.0
        lo, hi = signed_interval(lo, hi, sign)
        unit = re.sub(r'\s+', ' ', re.sub(r'\s*\b(?:increase|decrease)\b', '', unit, flags=re.I)).strip() or None
        return 'beta', sign * abs(v), lo, hi, unit, 'decreases' if sign < 0 else 'increases'
    if v > 0:
        lo, hi = signed_interval(lo, hi, 1.0)
        return 'odds_ratio', v, lo, hi, unit, 'increases' if v > 1 else ('decreases' if v < 1 else None)
    return None, v, lo, hi, unit, None


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
                study[r['STUDY ACCESSION']] = dict(pubmed=r['PUBMED ID'], trait=text_or_none(r['DISEASE/TRAIT']), mapped_trait=r['MAPPED_TRAIT'], uris=r['MAPPED_TRAIT_URI'], cohort=text_or_none(r['COHORT']),
                                                   full=r['FULL SUMMARY STATISTICS'].lower() == 'yes', initial=r['INITIAL SAMPLE SIZE'], n_init=0, n_cases=0, anc=set(), anc_all=set())
        with open(os.path.join(self.root, 'gwas-catalog-download-ancestries-v1.0.3.1.tsv'), encoding='utf-8', errors='replace', newline='') as f:
            for r in csv.DictReader(f, delimiter='\t', quoting=csv.QUOTE_NONE):
                s = study.get(r['STUDY ACCESSION'])
                if s is None:
                    st.inc('ancestry_rows_without_study')
                    continue
                cat = r['BROAD ANCESTRAL CATEGORY'].strip()
                if cat.upper() in NOT_REPORTED:
                    st.inc('ancestry_rows_with_category_not_reported')
                    cat = None
                if cat:
                    s['anc_all'].add(cat)
                if r['STAGE'] == 'initial':
                    if cat:
                        s['anc'].add(cat)
                    n = num(r['NUMBER OF INDIVIDUALS'])
                    if n:
                        s['n_init'] += int(n)
                    c = num(r['NUMBER OF CASES'])
                    if c:
                        s['n_cases'] += int(c)
        for acc, s in sorted(study.items()):
            ctx.emit('entity_study', dict(study_id=acc, source_id='gwas_catalog', lineage_group=lineage_of(s['pubmed'], acc), study_type='gwas', pubmed=s['pubmed'], trait_ids=efo_ids(s['uris']),
                                          trait_text=s['trait'], n_total=s['n_init'] or None, n_cases=s['n_cases'] or None, ancestry='; '.join(sorted(s['anc'])) or None, cohorts=s['cohort'],
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
                res = resolve_snv(self.refbase(ctx, canon, k[2]), k[4], frequency(r['RISK ALLELE FREQUENCY']), {c[0] for c in cands}, None)
                if res.status == 'rejected' and frequency(r['RISK ALLELE FREQUENCY']) is not None and len(cands) == 1:
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
        raf = frequency(r['RISK ALLELE FREQUENCY'])
        raf_reported = None
        if raf is None and num(r['RISK ALLELE FREQUENCY']) is not None:           # percentages and typing errors: not a frequency, and no evidence for a strand
            raf_reported = r['RISK ALLELE FREQUENCY']
            st.inc('risk_allele_frequency_outside_0_1')
        p = num(r['P-VALUE'])
        common = dict(lineage_group=lineage_of(r['PUBMEDID'], r['STUDY ACCESSION']), lineage_note=f"GWAS Catalog study {r['STUDY ACCESSION']}; PMID {r['PUBMEDID']}; cohort {s.get('cohort') or 'not stated'}; other databases that derive from this paper are not independent replication",
                      condition_ids=efo_ids(r['MAPPED_TRAIT_URI']), condition_text=text_or_none(r['DISEASE/TRAIT']), gene_ids=[g.strip() for g in re.split(r'[;,]', r['SNP_GENE_IDS']) if g.strip().startswith('ENSG')] or None,
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
        if raf_reported is not None:
            extra['risk_allele_frequency_reported'] = raf_reported
        ctx.emit_evidence(evidence_row(rel, 'gwas_catalog', acc, f'L{ln}', raw, 'variant', subject, 'association_replication', 'statistical_association', 'gwas_association',
                                       text_or_none(r['MAPPED_TRAIT']) or text_or_none(r['DISEASE/TRAIT']), ev_assert, mapping, reason, discriminator=f'L{ln}', p_value=p, eaf=raf, effect_allele=eff_allele, other_allele=other_allele,
                                       extra_json=json.dumps(extra, sort_keys=True), **c0, **common))
        kind_e, val, lo, hi, unit, direction = parse_effect(r['OR or BETA'], r['95% CI (TEXT)'])
        if kind_e:
            effect_known = eff_allele is not None
            ctx.emit_evidence(evidence_row(rel, 'gwas_catalog', acc, f'L{ln}', raw, 'variant', subject, 'effect_magnitude_direction', 'statistical_association', 'gwas_effect_size',
                                           kind_e, ev_assert, mapping, reason, discriminator=f'L{ln}:effect', p_value=p, eaf=raf, effect_allele=eff_allele, other_allele=other_allele,
                                           beta=val if kind_e == 'beta' else None, odds_ratio=val if kind_e == 'odds_ratio' else None, ci_low=lo, ci_high=hi, score_unit=unit,
                                           direction=direction if effect_known else None, extra_json=json.dumps(dict(raw_value=r['OR or BETA'], ci_text=r['95% CI (TEXT)'], direction_reported=direction,
                                           direction_attributed_to_allele=effect_known, effect_type_basis=effect_type_basis(kind_e, r['95% CI (TEXT)'])), sort_keys=True), **c0, **common))
        if subject:
            for g in common['gene_ids'] or ():
                ctx.emit('variant_gene', dict(variant_key=subject, gene_id=g, relation='gwas_catalog_snp_gene', consequence=r['CONTEXT'].split(';')[0].strip(), release_id=rel), name='gwas_catalog')
            if rs:
                ctx.emit('xref', dict(namespace='dbsnp', identifier=rs, entity_level='variant', entity_id=subject, release_id=rel, relation='gwas_catalog_snps'), name='gwas_catalog')
