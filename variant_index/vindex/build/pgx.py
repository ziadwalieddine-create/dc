"""Pharmacogenomic haplotype definitions and function: CPIC (API tables) and PharmVar (complete database 6.2.29).

Haplotypes are entities of their own (system:gene:allele). PharmVar and CPIC name the same star alleles, so the two systems are linked by a
shared star-allele name cross-reference, not merged. Allele function is haplotype-level evidence (axis effect_magnitude_direction); CPIC drug
recommendations are gene-level evidence (axis clinical_validity_actionability) that carries the phenotype it was written for. Nothing is
copied onto a variant: a variant belongs to a haplotype only through haplotype_definition rows.

CPIC definitions list, for every defining position, the base that each allele has there, including the reference base, so a definition row
is a variant only where the allele differs from the GRCh38 sequence. A row that cannot be turned into a normalised key is kept with
definition_scope 'unresolved:<reason>' and the matcher reports that allele as unevaluable (never silently absent): repeat-count notation,
IUPAC ambiguity codes, reference-base conflicts, alleles with no location rows.
"""
import json
import os
import re
import zipfile

from ..normalise import normalise
from .genes import entity_id
from .common import DATA, evidence_row
from .hgvs import change_to_key, split_changes

CPIC = os.path.join(DATA, 'cpic')
PV_ZIP = os.path.join(DATA, 'pharmvar', 'pharmvar-6.2.29.zip')

FUNC_EFFECT = {'normal function': 'no_effect', 'normal': 'no_effect', 'increased function': 'effect', 'decreased function': 'effect', 'no function': 'effect',
               'possible decreased function': 'effect', 'possible increased function': 'effect', 'uncertain function': 'uncertain', 'unknown function': 'uncertain',
               'uncertain': 'uncertain', 'variant of uncertain significance': 'uncertain', 'likely benign': 'no_effect', 'benign': 'no_effect',
               'likely pathogenic': 'effect', 'pathogenic': 'effect'}
FUNC_DIRECTION = {'increased function': 'gain_of_function', 'possible increased function': 'gain_of_function', 'decreased function': 'loss_of_function',
                  'possible decreased function': 'loss_of_function', 'no function': 'loss_of_function'}
_IUPAC = re.compile(r'^[ACGTRYMKSWBDHVN]+$')


_IUPAC_BASES = {'R': 'AG', 'Y': 'CT', 'M': 'AC', 'K': 'GT', 'S': 'CG', 'W': 'AT', 'B': 'CGT', 'D': 'AGT', 'H': 'ACT', 'V': 'ACG', 'N': 'ACGT', 'A': 'A', 'C': 'C', 'G': 'G', 'T': 'T'}


def _iupac(ref, chrom, parsed, allele):
    """An IUPAC code in a CPIC allele string marks a position where the allele tolerates more than one base. When the code is the reference base plus
    the listed alternative base(s) the position does not constrain the allele ('tolerated'); any other code is left unresolved."""
    ref_bases, alts = set(), set()
    for start, end, rest in parsed:
        sub = re.match(r'^([ACGT])>([ACGT])$', rest)
        if not sub or start != end:
            return 'unresolved', [], 'iupac_code_on_non_substitution'
        try:
            span = ref.fetch(chrom, start, end)
        except ValueError:
            return 'unresolved', [], 'outside_contig'
        if span != sub.group(1):
            return 'unresolved', [], f'reference_base_differs:{span}'
        ref_bases.add(sub.group(1))
        alts.add((start, sub.group(1), sub.group(2)))
    if len(ref_bases) != 1:
        return 'unresolved', [], 'iupac_code_on_mixed_reference_bases'
    rb = next(iter(ref_bases))
    allowed = set(''.join(_IUPAC_BASES[c] for c in allele))
    listed = {a for _, _, a in alts}
    if rb not in allowed:
        return 'unresolved', [], 'iupac_code_excludes_reference'
    if not allowed <= listed | {rb}:
        return 'unresolved', [], 'iupac_code_outside_listed_changes'
    keys = []
    for start, r, a in sorted(alts):
        if a in allowed:
            st, ks, why = change_to_key(ref, chrom, start, start, f'{r}>{a}')
            if st != 'variant':
                return st, ks, why
            keys += ks
    return 'tolerated', keys, None


def _allele_strings(rest, span):
    """(allele string of the variant, allele string of the reference sequence) for one change, in CPIC's allele-table conventions."""
    sub = re.match(r'^[ACGT]>([ACGT])$', rest)
    if sub:
        return sub.group(1), span
    if re.match(r'^del([ACGT]*|\d+)$', rest):
        return 'del' + span, span
    if re.match(r'^dup([ACGT]*|\d+)$', rest):
        return span * 2, span
    ins = re.match(r'^ins([ACGT]+)$', rest)
    if ins:
        return 'ins' + ins.group(1), 'del'
    di = re.match(r'^del[ACGT]*ins([ACGT]+)$', rest)
    if di:
        return 'ins' + di.group(1), 'del' + span
    return None, None


def cpic_variants(ref, chrom, loc, allele):
    """What a CPIC allele string means at one sequence location (GRCh38 HGVS g./m. notation).

    Returns (status, keys, reason): 'variant' with the normalised CHROM:POS:REF:ALT key, 'reference' (the allele has the GRCh38 sequence there),
    'tolerated' (an IUPAC code that allows the reference base and the listed alternative base(s); keys are those alternatives), or 'unresolved'
    with the reason. `chrom` is the contig name in `ref`."""
    allele = (allele or '').strip()
    loc = (loc or '').strip()
    if not [c for c in loc.split(';') if c.strip()]:
        return 'unresolved', [], 'empty_location'
    if '[' in loc or '(' in loc or '(' in allele or ';' in allele or ' or ' in allele:
        return 'unresolved', [], 'repeat_length_notation'
    parsed = split_changes(loc)
    if parsed is None:
        return 'unresolved', [], 'unparsed_notation'
    if _IUPAC.match(allele) and re.search(r'[RYMKSWBDHVN]', allele):
        return _iupac(ref, chrom, parsed, allele)
    is_reference = False
    for start, end, rest in parsed:
        status, keys, why = change_to_key(ref, chrom, start, end, rest)
        if status != 'variant':
            return 'unresolved', [], why or 'unsupported_notation'
        variant_string, reference_string = _allele_strings(rest, ref.fetch(chrom, start, end))
        if variant_string is None:
            return 'unresolved', [], 'unsupported_notation'
        if allele == variant_string:
            return 'variant', keys, None
        is_reference = is_reference or allele == reference_string
    if is_reference:
        return 'reference', [], None
    return 'unresolved', [], 'allele_string_matches_neither_reference_nor_a_change'


def snapshot(ctx, target):
    m = ctx.manifest.get(target, {})
    return 'api-' + (m.get('fetched_utc') or 'undated')[:10]


class CpicAdapter:
    source_id = 'cpic'
    targets = ['cpic_api_allele', 'cpic_api_allele_definition', 'cpic_api_allele_location_value', 'cpic_api_sequence_location', 'cpic_api_recommendation', 'cpic_genes', 'cpic_guidelines']

    def __init__(self, genes, root=None):
        self.genes, self.root = genes, root or CPIC

    def load(self, name, sub='api'):
        with open(os.path.join(self.root, sub, f'{name}.json') if sub else os.path.join(self.root, f'{name}.json'), encoding='utf-8') as f:
            return json.load(f)

    def run(self, ctx, limit=None):
        ctx.adapter = self
        st = ctx.st()
        rel_def = ctx.register_release('cpic', 'cpic_api_allele_definition', snapshot(ctx, 'cpic_api_allele_definition'))
        rel_fun = ctx.register_release('cpic', 'cpic_api_allele', snapshot(ctx, 'cpic_api_allele'))
        rel_rec = ctx.register_release('cpic', 'cpic_api_recommendation', snapshot(ctx, 'cpic_api_recommendation'))
        chrom_of = {g['symbol']: (g.get('chr') or '') for g in self.load('cpic_genes', None)}
        guideline = {g['id']: g for g in self.load('cpic_guidelines', None)}
        loc = {x['id']: x for x in self.load('sequence_location')}
        defs = {x['id']: x for x in self.load('allele_definition')}
        by_def, func_of = {}, {}
        for v in self.load('allele_location_value'):
            by_def.setdefault(v['alleledefinitionid'], []).append(v)
        for a in self.load('allele'):
            func_of.setdefault(a['definitionid'], []).append(a)
        for did, d in sorted(defs.items()):
            gene, name = d['genesymbol'], d['name']
            hid = f'cpic:{gene}:{name}'
            canon = ctx.ref.contig(chrom_of.get(gene) or '')
            keys, unresolved = [], []
            rows = by_def.get(did, ())
            for v in rows:
                l = loc.get(v['locationid'])
                if l is None:
                    unresolved.append(('?', 'location_record_missing'))
                    continue
                if canon is None:
                    unresolved.append((l['chromosomelocation'], 'gene_chromosome_unknown' if not chrom_of.get(gene) else 'contig_not_in_reference'))
                    continue
                status, ks, reason = cpic_variants(ctx.ref, canon, l['chromosomelocation'], v['variantallele'])
                if status == 'variant':
                    keys += [(k, f"{l['chromosomelocation']} allele={v['variantallele']}", 'defining') for k in ks]
                elif status == 'tolerated':
                    keys += [(k, f"{l['chromosomelocation']} allele={v['variantallele']}", 'tolerated') for k in ks]
                elif status == 'unresolved':
                    unresolved.append((f"{l['chromosomelocation']} allele={v['variantallele']}", reason))
            no_rows = not rows and not d['matchesreferencesequence']
            funcs = func_of.get(did, [])
            fn = next((f['functionalstatus'] or f['clinicalfunctionalstatus'] for f in funcs if f['functionalstatus'] or f['clinicalfunctionalstatus']), None)
            ctx.emit('entity_haplotype', dict(haplotype_id=hid, system='cpic', gene_symbol=gene, allele_name=name, is_reference=bool(d['matchesreferencesequence']),
                                              structural=bool(d['structuralvariation']), release_id=rel_def, accession=str(did), function=fn, allele_level='cpic'), name='cpic')
            for k, h, scope in keys:
                ctx.emit('haplotype_definition', dict(haplotype_id=hid, variant_key=k, hgvs=h, release_id=rel_def, definition_scope=scope), name='cpic')
            for h, reason in unresolved:
                ctx.emit('haplotype_definition', dict(haplotype_id=hid, variant_key=None, hgvs=h, release_id=rel_def, definition_scope=f'unresolved:{reason.split(":")[0]}'), name='cpic')
                st.inc(f'definition_rows_unresolved:{reason.split(":")[0]}')
            if no_rows:
                ctx.emit('haplotype_definition', dict(haplotype_id=hid, variant_key=None, hgvs=None, release_id=rel_def, definition_scope='unresolved:no_location_rows'), name='cpic')
                st.inc('alleles_without_location_rows')
            st.inc('alleles')
            st.inc('definition_rows_defining', sum(1 for k in keys if k[2] == 'defining'))
            st.inc('definition_rows_tolerated', sum(1 for k in keys if k[2] == 'tolerated'))
            if d.get('pharmvarid'):
                ctx.emit('xref', dict(namespace='pharmvar_id', identifier=str(d['pharmvarid']), entity_level='haplotype', entity_id=hid, release_id=rel_def, relation='cpic_definition_pharmvarid'), name='cpic')
            ident = f'{gene}{name}' if name.startswith('*') else f'{gene} {name}'
            ctx.emit('xref', dict(namespace='pgx_allele_name', identifier=ident, entity_level='haplotype', entity_id=hid, release_id=rel_def, relation='cpic_name'), name='cpic')
            ctx.offer('pgx_allele_name', ident, hid)
            if d['matchesreferencesequence']:
                ctx.offer('pgx_reference_allele', hid, hid)
            for f in funcs:
                fs = f['functionalstatus'] or f['clinicalfunctionalstatus']
                if not fs:
                    st.inc('allele_rows_without_function')
                    continue
                low = fs.lower()
                try:
                    avn = float(f.get('activityvalue'))
                except (TypeError, ValueError):
                    avn = None
                extra = dict(activity_value=f.get('activityvalue'), clinical_functional_status=f.get('clinicalfunctionalstatus'), functional_status=f.get('functionalstatus'),
                             clinical_functional_substrate=f.get('clinicalfunctionalsubstrate'), findings=f.get('findings'), comments=f.get('functioncomments'),
                             citations=f.get('citations'), frequency=f.get('frequency'), inferred_frequency=f.get('inferredfrequency'))
                ctx.emit_evidence(evidence_row(rel_fun, 'cpic', f"{gene}:{name}:{f['id']}", f"allele:{f['id']}", json.dumps(f, sort_keys=True, ensure_ascii=False), 'haplotype', hid,
                                               'effect_magnitude_direction', 'curated_expert_classification', 'allele_functional_status', fs, FUNC_EFFECT.get(low, 'not_applicable'),
                                               direction=FUNC_DIRECTION.get(low), claim_num=avn, review_status=f.get('strength'), curator='CPIC',
                                               lineage_group=f'cpic_allele:{gene}:{name}', gene_ids=[entity_id(self.genes, symbol=gene)], origin='germline', extra_json=json.dumps(extra, sort_keys=True, default=str)))
                st.inc('allele_function_rows')
        for r in self.load('recommendation'):
            g_info = guideline.get(r['guidelineid'], {})
            genes = sorted(set(r.get('phenotypes') or {}) | set(r.get('lookupkey') or {}) | set(r.get('implications') or {}) | set(r.get('allelestatus') or {}))
            raw = json.dumps(r, sort_keys=True, ensure_ascii=False)
            for g in genes or ['']:
                gid, _how = self.genes.resolve(symbol=g) if g else (None, 'unresolved')
                specific = {k: (r.get(k) or {}).get(g) for k in ('phenotypes', 'activityscore', 'lookupkey', 'implications', 'allelestatus')}
                extra = json.dumps(dict(gene_specific=specific, combination_genes=genes, drug_recommendation=r.get('drugrecommendation'), population=r.get('population'),
                                        comments=r.get('comments'), dosing_information=r.get('dosinginformation'), alternate_drug_available=r.get('alternatedrugavailable'),
                                        other_prescribing_guidance=r.get('otherprescribingguidance'), guideline=g_info.get('name')), sort_keys=True, ensure_ascii=False)
                kw = dict(curator='CPIC', lineage_group=f"cpic_guideline:{r['guidelineid']}", condition_ids=[r['drugid']], condition_text=g_info.get('name'), extra_json=extra, review_status=None if (r.get('classification') or '').strip().lower() in ('', 'n/a') else r.get('classification'))
                common = (rel_rec, 'cpic', f"rec:{r['id']}:{g}", f"rec:{r['id']}", raw)
                if gid:
                    ctx.emit_evidence(evidence_row(*common, 'gene', gid, 'clinical_validity_actionability', 'curated_expert_classification', 'cpic_drug_recommendation',
                                                   r.get('classification') or 'unclassified', 'not_applicable', 'mapped', None, discriminator=g, **kw))
                else:
                    ctx.emit_evidence(evidence_row(*common, 'gene', None, 'clinical_validity_actionability', 'curated_expert_classification', 'cpic_drug_recommendation',
                                                   r.get('classification') or 'unclassified', 'not_applicable', 'unmapped', 'gene_symbol_not_in_gene_map', discriminator=g, unmapped_key=g or None, **kw))
                st.inc('recommendation_rows')
        return st.c


class PharmVarAdapter:
    source_id = 'pharmvar'
    targets = ['pharmvar_6_2_29']

    def __init__(self, path=None):
        self.path = path or PV_ZIP

    def run(self, ctx, limit=None):
        ctx.adapter = self
        st = ctx.st()
        rel = ctx.register_release('pharmvar', 'pharmvar_6_2_29', '6.2.29')
        rs_seen = set()
        with zipfile.ZipFile(self.path) as z:
            for n in sorted(x for x in z.namelist() if x.endswith('.tsv') and '/GRCh38/' in x):
                lines = [l for l in z.read(n).decode('utf-8').split('\n') if l and not l.startswith('#')]
                hdr = lines[0].split('\t')
                defs = {}
                for l in lines[1:]:
                    r = dict(zip(hdr, l.split('\t')))
                    defs.setdefault(r['Haplotype Name'], []).append(r)
                for hap, rows in sorted(defs.items()):
                    gene = rows[0]['Gene']
                    short = hap[len(gene):] if hap.startswith(gene) else hap
                    hid = f'pharmvar:{gene}:{short}'
                    is_ref = all(r['ReferenceSequence'] == 'REFERENCE' for r in rows)
                    ctx.emit('entity_haplotype', dict(haplotype_id=hid, system='pharmvar', gene_symbol=gene, allele_name=short, is_reference=is_ref, structural=False, release_id=rel,
                                                      accession=hap, function=None, allele_level='core' if '.' not in short else 'suballele'), name='pharmvar')
                    st.inc('haplotypes')
                    for r in rows:
                        if r['ReferenceSequence'] == 'REFERENCE':
                            continue
                        key, why = self.variant(ctx, r)
                        h = f"{r['ReferenceSequence']}:{r['Variant Start']}-{r['Variant Stop']}:{r['Reference Allele']}>{r['Variant Allele']}:{r['rsID']}"
                        if key:
                            ctx.emit('haplotype_definition', dict(haplotype_id=hid, variant_key=key, hgvs=h, release_id=rel, definition_scope='defining'), name='pharmvar')
                            if r['rsID'] and (r['rsID'], key) not in rs_seen:
                                rs_seen.add((r['rsID'], key))
                                ctx.emit('xref', dict(namespace='dbsnp', identifier=r['rsID'], entity_level='variant', entity_id=key, release_id=rel, relation='pharmvar_rsID'), name='pharmvar')
                        else:
                            ctx.emit('haplotype_definition', dict(haplotype_id=hid, variant_key=None, hgvs=h, release_id=rel, definition_scope=f'unresolved:{why}'), name='pharmvar')
                            st.inc(f'definition_rows_unresolved:{why}')
                        st.inc('definition_rows')
                    if '.' not in short:
                        ident = f'{gene}{short}' if short.startswith('*') else f'{gene} {short}'
                        ctx.emit('xref', dict(namespace='pgx_allele_name', identifier=ident, entity_level='haplotype', entity_id=hid, release_id=rel, relation='pharmvar_name'), name='pharmvar')
                        ctx.offer('pgx_allele_name', ident, hid)
                    if is_ref:
                        ctx.offer('pgx_reference_allele', hid, hid)
        return st.c

    def variant(self, ctx, r):
        """Normalised key for one PharmVar definition row, or (None, reason). Deletions run Start..Stop with '-' as the allele; insertions give the two flanking positions."""
        chrom = ctx.ref.contig(r['ReferenceSequence'])
        if chrom is None:
            return None, 'reference_sequence_unmapped'
        try:
            s, e = int(r['Variant Start']), int(r['Variant Stop'])
        except ValueError:
            return None, 'position_not_numeric'
        ref, alt, typ = r['Reference Allele'], r['Variant Allele'], r['Type']
        try:
            if typ == 'substitution':
                n = normalise(ctx.ref, chrom, s, ref, alt)
            elif typ == 'deletion':
                if s < 2:
                    return None, 'deletion_at_contig_start'
                anchor = ctx.ref.fetch(chrom, s - 1, s - 1)
                span = ctx.ref.fetch(chrom, s, e)
                if ref not in ('-', '') and ref.upper() != span:
                    return None, 'deleted_sequence_differs_from_reference'
                n = normalise(ctx.ref, chrom, s - 1, anchor + span, anchor)
            elif typ == 'insertion':
                if e != s + 1:
                    return None, 'insertion_flanks_not_adjacent'
                b = ctx.ref.fetch(chrom, s, s)
                n = normalise(ctx.ref, chrom, s, b, b + alt)
            else:
                return None, f'type_{typ or "empty"}'
        except ValueError:
            return None, 'outside_contig'
        return (n.key, None) if n.ok else (None, n.status)
