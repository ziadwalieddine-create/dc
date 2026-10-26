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

    def keys(self, path, canon, sample=None, carried_only=False):
        """{normalised key} for the plain alleles of `canon` in the file; with carried_only, only those the sample carries."""
        out = set()
        for n in (n for n in self.contig_names(path) if self.ref.contig(n) == canon):
            stages = [[self.bcf, 'view', '-r', n] + (['-s', sample] if sample else []) + [self.indexed(path)]]
            if n != canon:                                                      # a header without ##contig lines (ClinVar) cannot be renamed, and needs no renaming
                rename = os.path.join(self.tmp, f'rename_{len(self._idx)}.txt')
                with open(rename, 'w') as f:
                    f.write(f'{n}\t{canon}\n')
                stages.append([self.bcf, 'annotate', '--rename-chrs', rename])
            stages += [[self.bcf, 'reheader', '-f', self.fai], [self.bcf, 'norm', '-m', '-any', '-f', self.fasta, '-c', 'w']]
            if carried_only:
                stages.append([self.bcf, 'view', '-i', 'GT="alt"'])
            stages.append([self.bcf, 'query', '-f', '%CHROM\t%POS\t%REF\t%ALT\n'])
            for line in _chain(stages):
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
        calls_m = read_column(results, 'small_calls', 'key')
        src_m = read_column(results, 'evidence_links', 'call_id', where=[('call_kind', 'small'), ('source_id', evidence_source)])
        rejected = {(r['chrom'], r['pos']) for r in pq.read_table(os.path.join(results, 'rejected_calls.parquet')).to_pylist()}
        per, only_b, only_m, ev_b_all, ev_only_b, ev_only_m = {}, [], [], 0, [], []
        calls_b_all = 0
        for canon in order:
            calls_b = norm.keys(vcf, canon, sample=sample, carried_only=True)
            m_here = {k for k in calls_m if k.split(':', 1)[0] == canon}
            ev_here = {k for k in src_m if k.split(':', 1)[0] == canon}
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
            calls=dict(bcftools=calls_b_all, matcher=sum(v['calls_matcher'] for v in per.values()), only_bcftools=len(only_b), only_bcftools_explained_by_rejected_calls=len(explained),
                       only_bcftools_unexplained=len(unexplained_b), only_matcher=len(only_m)),
            evidence=dict(bcftools=ev_b_all, matcher=sum(v['evidence_matcher'] for v in per.values()), only_bcftools=len(ev_only_b), only_bcftools_explained_by_rejected_calls=len(ev_explained),
                          only_bcftools_unexplained=len(ev_unexplained), only_matcher=len(ev_only_m)) if evidence_vcf else None,
            examples=dict(only_bcftools_unexplained=unexplained_b[:examples], only_matcher=only_m[:examples], evidence_only_bcftools=ev_unexplained[:examples], evidence_only_matcher=ev_only_m[:examples]),
            per_contig=per)
        report['ok'] = not unexplained_b and not only_m and not ev_unexplained and not ev_only_m
        return report
    finally:
        shutil.rmtree(tmp, ignore_errors=True)


def clinvar_fields(index_dir, vcf, reference, contigs=None, bcftools=None, examples=10):
    """Compare every ClinVar-VCF-derived evidence row of the index with the VCF itself, read by bcftools: the classification string, review status and condition text of
    each (VariationID, AlleleID) must be the ones the file states, no record with a classification may be missing from the index, and the index may hold no such row
    that the file lacks. bcftools is the independent parser; the comparison is one contig at a time."""
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
        total = dict(vcf_records_with_classification=0, index_rows=0, missing_in_index=0, extra_in_index=0, claim_value_differs=0, review_status_differs=0, condition_text_differs=0)
        ex = {k: [] for k in ('missing_in_index', 'extra_in_index', 'claim_value_differs', 'review_status_differs', 'condition_text_differs')}
        per = {}
        for label, canon_list, part in passes:
            names = [n for n in names_in_file if (reference.contig(n) in canon_list if canon_list else not (reference.contig(n) and reference.is_primary(reference.contig(n))))]
            want = {}
            for n in names:
                stages = [[bcf, 'query', '-r', n, '-f', '%ID\t%INFO/ALLELEID\t%INFO/CLNSIG\t%INFO/CLNREVSTAT\t%INFO/CLNDN\n', norm.indexed(vcf)]]
                for line in _chain(stages):
                    vid, aid, sig, rev, dn = line.rstrip('\n').split('\t')
                    if sig not in ('.', ''):
                        want[f'VariationID:{vid}|AlleleID:{aid}'] = (sig, rev if rev != '.' else '', dn.replace('_', ' ') if dn != '.' else '')
            got = {}
            p = os.path.join(index_dir, 'evidence', 'variant', part)
            files = [os.path.join(p, f) for f in sorted(os.listdir(p)) if f.endswith('.parquet')] if os.path.isdir(p) else []
            if files:
                import pyarrow.compute as pc
                flt = (pc.field('source_id') == 'clinvar') & (pc.field('claim_type') == 'germline_classification') & pc.match_substring_regex(pc.field('locator'), '^L[0-9]+$')
                if canon_list is None:
                    pass
                else:
                    flt = flt & (pc.field('chrom') == canon_list[0])
                t = pds.dataset(files, format='parquet').to_table(columns=['accession', 'claim_value', 'review_status', 'condition_text', 'chrom'], filter=flt)
                for r in t.to_pylist():
                    got[r['accession']] = (r['claim_value'], r['review_status'] or '', r['condition_text'] or '')
            missing = sorted(set(want) - set(got))
            extra = sorted(set(got) - set(want))
            diff = {'claim_value_differs': [], 'review_status_differs': [], 'condition_text_differs': []}
            for k in set(want) & set(got):
                for i, name in enumerate(('claim_value_differs', 'review_status_differs', 'condition_text_differs')):
                    if want[k][i] != got[k][i]:
                        diff[name].append(k)
            per[label] = dict(vcf_records_with_classification=len(want), index_rows=len(got), missing_in_index=len(missing), extra_in_index=len(extra), **{k: len(v) for k, v in diff.items()})
            total['vcf_records_with_classification'] += len(want)
            total['index_rows'] += len(got)
            total['missing_in_index'] += len(missing)
            total['extra_in_index'] += len(extra)
            for k, v in diff.items():
                total[k] += len(v)
            for k, v in (('missing_in_index', missing), ('extra_in_index', extra), *[(k, sorted(v)) for k, v in diff.items()]):
                ex[k] += v[:examples]
        report = dict(bcftools=subprocess.run([bcf, '--version'], capture_output=True, text=True).stdout.split('\n')[0], totals=total, examples={k: v[:examples] for k, v in ex.items()}, per_contig=per)
        report['ok'] = all(total[k] == 0 for k in ('missing_in_index', 'extra_in_index', 'claim_value_differs', 'review_status_differs', 'condition_text_differs'))
        return report
    finally:
        shutil.rmtree(tmp, ignore_errors=True)
