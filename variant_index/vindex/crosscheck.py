"""Independent cross-check of a whole-genome run against established command-line tools.

The matcher normalises with its own code (vindex.normalise). This module recomputes two sets with `bcftools` only and compares them with the run's tables:

  calls      the sample's called small alleles, normalised with `bcftools norm -m -any -f FASTA`, genotype-filtered to the alleles the sample carries
             -> must equal the keys of small_calls.parquet
  evidence   the same keys intersected with a source VCF (ClinVar's) normalised the same way
             -> must equal the keys that have evidence_links from that source

A key present on one side only is reported with its explanation: a bcftools-only key that the matcher refused and listed in rejected_calls (REF disagrees with the
reference) is explained; anything else is a difference. Nothing here is sent anywhere: the sample VCF is read by a local bcftools process, and only counts and
example keys from the comparison are written to the report.
"""
import json
import os
import re
import shutil
import subprocess
import tempfile

import pyarrow.compute as pc
import pyarrow.parquet as pq

PLAIN_ALLELE = re.compile(r'^[ACGTN]+$')


def find_bcftools():
    d = os.environ.get('VINDEX_TOOLS_BIN')
    if d and os.path.exists(os.path.join(d, 'bcftools')):
        return os.path.join(d, 'bcftools')
    return shutil.which('bcftools')


def _chain(stages):
    """Run argv lists connected by pipes; yield the text lines of the last one. Raises with the failing stage's stderr (kept in files: a stage that warns for every
    record must not block on a full pipe)."""
    procs, errs, prev = [], [], None
    for i, argv in enumerate(stages):
        err = tempfile.TemporaryFile()
        p = subprocess.Popen(argv, stdin=prev.stdout if prev else None, stdout=subprocess.PIPE, stderr=err, text=(i == len(stages) - 1))
        if prev:
            prev.stdout.close()
        procs.append(p)
        errs.append(err)
        prev = p
    try:
        yield from procs[-1].stdout
    finally:
        procs[-1].stdout.close()
        failed = []
        for argv, p, err in zip(stages, procs, errs):
            rc = p.wait()
            if rc:
                err.seek(0)
                failed.append(f"{' '.join(argv[:3])}: exit {rc}: {err.read().decode(errors="replace")[-600:]}")
            err.close()
        if failed:
            raise RuntimeError('; '.join(failed))


class Normaliser:
    """bcftools-based normalised keys of one contig of one VCF."""

    def __init__(self, bcftools, reference, tmp):
        self.bcf, self.ref, self.tmp = bcftools, reference, tmp
        self.fai = reference.path + '.fai'
        self.fasta = reference.path
        self._idx = {}
        self._renames = 0

    def indexed(self, path):
        """`path` or `path##idx##index` so that region queries work on a file whose index is absent (the acquired directory is never written to)."""
        if path in self._idx:
            return self._idx[path]
        if any(os.path.exists(path + s) for s in ('.tbi', '.csi')):
            self._idx[path] = path
        else:
            idx = os.path.join(self.tmp, os.path.basename(path) + '.tbi')
            subprocess.run([self.bcf, 'index', '-t', '-f', '-o', idx, path], check=True, capture_output=True)
            self._idx[path] = f'{path}##idx##{idx}'
        return self._idx[path]

    def contig_names(self, path):
        """Contig names that hold records, as the (possibly built) index lists them."""
        r = subprocess.run([self.bcf, 'index', '-s', self.indexed(path)], capture_output=True, text=True, check=True)
        return [line.split('\t')[0] for line in r.stdout.splitlines() if line.strip()]

    def stream(self, path, name, canon, fmt, sample=None, carried_only=False, check_ref='w'):
        """Lines of `bcftools query -f fmt` over the records of contig `name` of the file, normalised against the reference (split, left-aligned and trimmed). check_ref
        'w' keeps a record whose REF disagrees with the reference; 'x' drops it."""
        stages = [[self.bcf, 'view', '-r', name] + (['-s', sample] if sample else []) + [self.indexed(path)]]
        if name != canon:
            self._renames += 1
            declare = os.path.join(self.tmp, f'declare_{self._renames}.fai')
            with open(declare, 'w') as f:                                       # annotate --rename-chrs refuses records on a contig its header does not define (the ClinVar header defines none)
                f.write(f'{name}\t{self.ref.length(canon)}\t0\t0\t0\n')
            rename = os.path.join(self.tmp, f'rename_{self._renames}.txt')
            with open(rename, 'w') as f:
                f.write(f'{name}\t{canon}\n')
            stages += [[self.bcf, 'reheader', '-f', declare], [self.bcf, 'annotate', '--rename-chrs', rename]]
        stages += [[self.bcf, 'reheader', '-f', self.fai], [self.bcf, 'norm', '-m', '-any', '-f', self.fasta, '-c', check_ref]]
        if carried_only:
            stages.append([self.bcf, 'view', '-i', 'GT="alt"'])
        stages.append([self.bcf, 'query', '-f', fmt])
        return _chain(stages)

    def keys(self, path, canon, sample=None, carried_only=False):
        """{normalised key} for the plain alleles of `canon` in the file; with carried_only, only those the sample carries."""
        out = set()
        for n in (n for n in self.contig_names(path) if self.ref.contig(n) == canon):
            for line in self.stream(path, n, canon, '%CHROM\t%POS\t%REF\t%ALT\n', sample=sample, carried_only=carried_only):
                c, pos, ref, alt = line.rstrip('\n').split('\t')
                if PLAIN_ALLELE.match(ref) and PLAIN_ALLELE.match(alt):
                    out.add(f'{c}:{pos}:{ref}:{alt}')
        return out


def read_column(results, table, column, where=None):
    """Set of the values of one column of a result table, optionally restricted to rows where one column equals a value."""
    cols = [column] + ([w for w, _ in where] if where else [])
    t = pq.read_table(os.path.join(results, f'{table}.parquet'), columns=cols)
    for w, v in (where or ()):
        t = t.filter(pc.equal(t[w], v))
    return set(t[column].to_pylist())


def crosscheck(vcf, results, reference, evidence_vcf=None, evidence_source='clinvar', contigs=None, sample=None, bcftools=None, examples=10):
    bcf = bcftools or find_bcftools()
    if not bcf:
        raise RuntimeError('bcftools not found (set VINDEX_TOOLS_BIN or put it on PATH)')
    from .wgs import chromosome_order
    order = chromosome_order(reference, contigs)
    tmp = tempfile.mkdtemp(prefix='crosscheck_')
    try:
        norm = Normaliser(bcf, reference, tmp)
        def by_contig(keys):
            out = {}
            for k in keys:
                out.setdefault(k.split(':', 1)[0], set()).add(k)
            return out
        calls_m = by_contig(read_column(results, 'small_calls', 'key'))
        src_m = by_contig(read_column(results, 'evidence_links', 'call_id', where=[('call_kind', 'small'), ('source_id', evidence_source)]))
        rejected = {(r['chrom'], r['pos']) for r in pq.read_table(os.path.join(results, 'rejected_calls.parquet')).to_pylist()}
        per, only_b, only_m, ev_b_all, ev_only_b, ev_only_m = {}, [], [], 0, [], []
        calls_b_all = 0
        for canon in order:
            calls_b = norm.keys(vcf, canon, sample=sample, carried_only=True)
            m_here, ev_here = calls_m.get(canon, set()), src_m.get(canon, set())
            ev_b = set()
            if evidence_vcf:
                ev_b = calls_b & norm.keys(evidence_vcf, canon)
            calls_b_all += len(calls_b)
            ev_b_all += len(ev_b)
            only_b += sorted(calls_b - m_here)
            only_m += sorted(m_here - calls_b)
            if evidence_vcf:
                ev_only_b += sorted(ev_b - ev_here)
                ev_only_m += sorted(ev_here - ev_b)
            per[canon] = dict(calls_bcftools=len(calls_b), calls_matcher=len(m_here), evidence_bcftools=len(ev_b), evidence_matcher=len(ev_here))
        refused = lambda k: (k.split(':')[0], int(k.split(':')[1])) in rejected
        explained = [k for k in only_b if refused(k)]
        unexplained_b = [k for k in only_b if not refused(k)]
        ev_explained = [k for k in ev_only_b if refused(k)]
        ev_unexplained = [k for k in ev_only_b if not refused(k)]
        report = dict(
            contigs=order, sample=sample, evidence_source=evidence_source if evidence_vcf else None, bcftools=subprocess.run([bcf, '--version'], capture_output=True, text=True).stdout.split('\n')[0],
            calls=dict(bcftools=calls_b_all, matcher=sum(v['calls_matcher'] for v in per.values()), matcher_on_contigs_not_compared=sum(len(v) for c, v in calls_m.items() if c not in set(order)),
                       only_bcftools=len(only_b), only_bcftools_explained_by_rejected_calls=len(explained),
                       only_bcftools_unexplained=len(unexplained_b), only_matcher=len(only_m)),
            evidence=dict(bcftools=ev_b_all, matcher=sum(v['evidence_matcher'] for v in per.values()), only_bcftools=len(ev_only_b), only_bcftools_explained_by_rejected_calls=len(ev_explained),
                          only_bcftools_unexplained=len(ev_unexplained), only_matcher=len(ev_only_m)) if evidence_vcf else None,
            examples=dict(only_bcftools_unexplained=unexplained_b[:examples], only_matcher=only_m[:examples], evidence_only_bcftools=ev_unexplained[:examples], evidence_only_matcher=ev_only_m[:examples]),
            per_contig=per)
        report['ok'] = not unexplained_b and not only_m and not ev_unexplained and not ev_only_m
        return report
    finally:
        shutil.rmtree(tmp, ignore_errors=True)


# ---- the ClinVar rows of the index against the ClinVar VCF, read by bcftools ---------------------------------------------------------------------------------------

# (claim type, classification tag, tag prefix of its review status / condition / submissions)
CLASSIFICATIONS = (('germline_classification', 'CLNSIG', 'CLN', 'CLNSIGSCV'), ('oncogenicity_classification', 'ONC', 'ONC', 'ONCSCV'),
                   ('somatic_clinical_impact', 'SCI', 'SCI', 'SCISCV'))
EXTRA_TAGS = ('CLNHGVS', 'CLNVC', 'CLNVCSO', 'MC', 'CLNVI', 'RS', 'CLNSIGCONF', 'CLNSIGSCV', 'DBVARID', 'ALLELEID')
FREQUENCIES = (('AF_ESP', 'GO-ESP'), ('AF_EXAC', 'ExAC'), ('AF_TGP', '1000G'))
# Restated from the ClinVar VCF header (ORIGIN) and from ClinVar's table of review statuses and stars, not imported from the adapter that is being checked.
ORIGIN_BITS = ((1, 'germline'), (2, 'somatic'), (4, 'inherited'), (8, 'paternal'), (16, 'maternal'), (32, 'de_novo'), (64, 'biparental'), (128, 'uniparental'),
               (256, 'not_tested'), (512, 'tested_inconclusive'), (1073741824, 'other'))
REVIEW_STARS = {'practice_guideline': 4, 'reviewed_by_expert_panel': 3, 'criteria_provided,_multiple_submitters,_no_conflicts': 2, 'criteria_provided,_multiple_submitters': 2,
                'criteria_provided,_single_submitter': 1, 'criteria_provided,_conflicting_classifications': 1, 'no_assertion_criteria_provided': 0,
                'no_classification_provided': 0, 'no_classification_for_the_single_variant': 0, 'no_classifications_from_unflagged_records': 0}
STRICT_ALLELE = re.compile(r'^[ACGT]+$')
RAW_TAGS = (['ALLELEID', 'ORIGIN', 'GENEINFO'] + [t for _, sig, pre, scv in CLASSIFICATIONS for t in (sig, pre + 'REVSTAT', pre + 'DN', pre + 'DISDB', scv)]
            + [t for t in EXTRA_TAGS if t not in ('ALLELEID', 'CLNSIGSCV')] + [t for t, _ in FREQUENCIES])


def _raw_format(declared):
    """The query format for the INFO tags the file declares; a tag its header lacks (older ClinVar releases have no ONC or SCI) reads as absent."""
    return '%ID\t' + '\t'.join(f'%INFO/{t}' if t in declared else '.' for t in RAW_TAGS) + '\n'
DIFFERENCES = ('missing_in_index', 'extra_in_index', 'claim_value_differs', 'review_status_differs', 'condition_text_differs', 'condition_ids_differ', 'n_submitters_differs',
               'review_stars_differ', 'origin_differs', 'gene_count_differs', 'extra_json_differs', 'frequency_differs', 'key_differs', 'mapped_but_bcftools_rejects',
               'unmapped_but_bcftools_accepts', 'gene_id_inconsistent')


def _origin(v):
    try:
        v = int(v)
    except (TypeError, ValueError):
        return 'unknown'
    names = [n for bit, n in ORIGIN_BITS if v & bit]
    rest = v
    for bit, _ in ORIGIN_BITS:
        rest &= ~bit
    names += [f'unrecognised_bit_{1 << i}' for i in range(rest.bit_length()) if rest >> i & 1]
    return '+'.join(names) if names else 'unknown'


def _ids(disdb):
    return sorted({x for grp in disdb.split('|') for x in grp.split(',') if x and x != '.'})


def clinvar_fields(index_dir, vcf, reference, contigs=None, bcftools=None, examples=10):
    """Compare the ClinVar-VCF-derived evidence rows of the index with the VCF itself, read by bcftools, record by record: the classification, review status, stars,
    condition text and identifiers, number of submissions, origin, number of genes, the carried-over INFO fields, the population frequencies and the normalised variant key
    (bcftools norm; a record bcftools drops for a REF that disagrees with the reference must be one the index refused). No record may be missing from the index and the
    index may hold no row the file lacks. The comparison is one contig at a time."""
    import collections
    import pyarrow.compute as pc
    import pyarrow.dataset as pds
    bcf = bcftools or find_bcftools()
    if not bcf:
        raise RuntimeError('bcftools not found (set VINDEX_TOOLS_BIN or put it on PATH)')
    from . import schema
    from .wgs import chromosome_order
    primary = chromosome_order(reference, contigs)
    passes = [(c, [c], c) for c in primary] + ([('other_contigs', None, schema.OTHER_PARTITION)] if contigs is None else [])
    tmp = tempfile.mkdtemp(prefix='clinvar_fields_')
    try:
        norm = Normaliser(bcf, reference, tmp)
        names_in_file = norm.contig_names(vcf)
        header = subprocess.run([bcf, 'view', '-h', norm.indexed(vcf)], capture_output=True, text=True, check=True).stdout
        raw_format = _raw_format(set(re.findall(r'^##INFO=<ID=([^,>]+)', header, re.M)))
        total = collections.Counter(vcf_records=0, vcf_classifications=0, index_rows=0, vcf_records_without_classification=0, mapped=0, unmapped_consistent=0)
        ex = {k: [] for k in DIFFERENCES}
        unmapped_reasons = collections.Counter()
        ncbi_gene = collections.defaultdict(set)
        per = {}
        for label, canon_list, part in passes:
            names = [n for n in names_in_file if (reference.contig(n) in canon_list if canon_list else not (reference.contig(n) and reference.is_primary(reference.contig(n))))]
            want, noclass, freq, keys, genes = {}, set(), {}, {}, {}
            n_records = 0
            for n in names:
                for line in _chain([[bcf, 'query', '-r', n, '-f', raw_format, norm.indexed(vcf)]]):
                    cols = line.rstrip('\n').split('\t')
                    f = dict(zip(RAW_TAGS, cols[1:]))
                    acc = f'VariationID:{cols[0]}|AlleleID:{f["ALLELEID"]}'
                    n_records += 1
                    found = False
                    for ctype, sig, pre, scv in CLASSIFICATIONS:
                        if f[sig] == '.':
                            continue
                        found = True
                        rev, dn, scvs = f[pre + 'REVSTAT'], f[pre + 'DN'], f[scv]
                        want[(ctype, acc)] = dict(
                            claim_value=f[sig], review_status='' if rev == '.' else rev.replace('_', ' '), condition_text='' if dn == '.' else dn.replace('_', ' '), condition_ids=_ids(f[pre + 'DISDB']),
                            n_submitters=(len([x for x in scvs.split('|') if x]) if scvs != '.' else 0) or None, review_stars=REVIEW_STARS.get(rev.strip().lower().replace(' ', '_')),
                            origin='somatic' if ctype != 'germline_classification' else _origin(f['ORIGIN'] if f['ORIGIN'] != '.' else None),
                            extra={t: f[t] for t in EXTRA_TAGS if f[t] != '.'})
                    if not found:
                        noclass.add(acc)
                    for tag, cohort in FREQUENCIES:
                        if f[tag] != '.':
                            freq[(acc, cohort)] = float(f[tag])
                    genes[acc] = [g for g in ('' if f['GENEINFO'] == '.' else f['GENEINFO']).split('|') if ':' in g]
                canon = reference.contig(n)
                if canon:
                    for line in norm.stream(vcf, n, canon, '%ID\t%INFO/ALLELEID\t%CHROM\t%POS\t%REF\t%ALT\n', check_ref='x'):
                        vid, aid, c, pos, ref, alt = line.rstrip('\n').split('\t')
                        if STRICT_ALLELE.match(ref) and STRICT_ALLELE.match(alt):
                            keys[f'VariationID:{vid}|AlleleID:{aid}'] = f'{c}:{pos}:{ref}:{alt}'
            counts = collections.Counter()

            def note(name, acc, ctype='germline_classification'):
                counts[name] += 1
                if len(ex[name]) < examples:
                    ex[name].append(acc if ctype == 'germline_classification' else f'{ctype}:{acc}')

            p = os.path.join(index_dir, 'evidence', 'variant', part)
            files = [os.path.join(p, f) for f in sorted(os.listdir(p)) if f.endswith('.parquet')] if os.path.isdir(p) else []
            seen, seen_freq, seen_noclass, key_checked = set(), set(), set(), set()
            if files:
                flt = (pc.field('source_id') == 'clinvar') & pc.match_substring_regex(pc.field('locator'), '^L[0-9]+$')
                if canon_list is not None:
                    flt = flt & (pc.field('chrom') == canon_list[0])
                cols = ['accession', 'claim_type', 'claim_value', 'review_status', 'condition_text', 'condition_ids', 'n_submitters', 'review_stars', 'origin', 'subject_id',
                        'mapping_status', 'mapping_reason', 'gene_ids', 'extra_json', 'score']
                for batch in pds.dataset(files, format='parquet').to_batches(columns=cols, filter=flt):
                    for r in batch.to_pylist():
                        acc, ctype = r['accession'], r['claim_type']
                        counts['index_rows'] += 1
                        if acc in genes and acc not in key_checked:                 # once per VCF record, whatever the row is
                            key_checked.add(acc)
                            gids = list(r['gene_ids'] or [])
                            if len(gids) != len(genes[acc]):
                                note('gene_count_differs', acc, ctype)
                            else:
                                for g, gid in zip(genes[acc], gids):
                                    ncbi_gene[g.partition(':')[2]].add(gid)
                            if r['mapping_status'] == 'mapped':
                                counts['mapped'] += 1
                                if acc not in keys:
                                    note('mapped_but_bcftools_rejects', acc, ctype)
                                elif keys[acc] != r['subject_id']:
                                    note('key_differs', acc, ctype)
                            else:
                                if acc in keys:
                                    note('unmapped_but_bcftools_accepts', acc, ctype)
                                else:
                                    counts['unmapped_consistent'] += 1
                                unmapped_reasons[r['mapping_reason'] or 'none'] += 1
                        if ctype == 'population_allele_frequency':
                            k = (acc, r['claim_value'])
                            seen_freq.add(k)
                            if k not in freq:
                                note('extra_in_index', f'{acc}|{r["claim_value"]}', ctype)
                            elif abs(freq[k] - (r['score'] if r['score'] is not None else float('nan'))) > 1e-9 * max(1.0, abs(freq[k])):
                                note('frequency_differs', f'{acc}|{r["claim_value"]}', ctype)
                            continue
                        if ctype == 'clinvar_record_without_classification':
                            seen_noclass.add(acc)
                            if acc not in noclass:
                                note('extra_in_index', acc, ctype)
                            continue
                        k = (ctype, acc)
                        seen.add(k)
                        w = want.get(k)
                        if w is None:
                            note('extra_in_index', acc, ctype)
                            continue
                        for field, name in (('claim_value', 'claim_value_differs'), ('review_status', 'review_status_differs'), ('condition_text', 'condition_text_differs'),
                                            ('n_submitters', 'n_submitters_differs'), ('review_stars', 'review_stars_differ'), ('origin', 'origin_differs')):
                            got = (r[field] or '') if isinstance(w[field], str) else r[field]
                            if got != w[field]:
                                note(name, acc, ctype)
                        if sorted(r['condition_ids'] or []) != w['condition_ids']:
                            note('condition_ids_differ', acc, ctype)
                        have = {t: v for t, v in json.loads(r['extra_json'] or '{}').items() if v not in ('', '.')}
                        if {t: have.get(t) for t in EXTRA_TAGS} != {t: w['extra'].get(t) for t in EXTRA_TAGS}:
                            note('extra_json_differs', acc, ctype)
            for k in sorted(set(want) - seen):
                note('missing_in_index', k[1], k[0])
            for k in sorted(set(freq) - seen_freq):
                note('missing_in_index', f'{k[0]}|{k[1]}', 'population_allele_frequency')
            for acc in sorted(noclass - seen_noclass):
                note('missing_in_index', acc, 'clinvar_record_without_classification')
            counts['vcf_records'] = n_records
            counts['vcf_classifications'] = len(want)
            counts['vcf_records_without_classification'] = len(noclass)
            per[label] = dict(counts)
            total.update(counts)
        inconsistent = sorted(n for n, g in ncbi_gene.items() if len(g) > 1)
        total['gene_id_inconsistent'] = len(inconsistent)
        ex['gene_id_inconsistent'] = [f'NCBI Gene {n}: {sorted(ncbi_gene[n])}' for n in inconsistent[:examples]]
        for k in DIFFERENCES:
            total.setdefault(k, 0)
        for c in per.values():
            for k in DIFFERENCES:
                c.setdefault(k, 0)
        report = dict(bcftools=subprocess.run([bcf, '--version'], capture_output=True, text=True).stdout.split('\n')[0], totals=dict(total), unmapped_by_reason=dict(unmapped_reasons),
                      examples=ex, per_contig=per)
        report['ok'] = all(total[k] == 0 for k in DIFFERENCES)
        return report
    finally:
        shutil.rmtree(tmp, ignore_errors=True)
