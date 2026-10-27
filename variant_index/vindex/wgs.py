"""Whole-genome matching of one sample against the index, one chromosome at a time.

A 30x genome holds millions of calls. Each chromosome is read, matched and written before the next is opened, so memory is bounded by the largest
chromosome. Everything the sample's file says is accounted for in the output tables; nothing is filtered out of them:

  small_calls         one row per called SNV / indel allele (dosage > 0): normalised key, genotype, phase, quality, QC flags, match status, flags (difficult
                      region, multiallelic split), the per-axis summary of the variant-level evidence that matched, evidence ids, genes, region kinds.
                      A call with no index record is kept and says NO_INDEX_RECORD: absence from the index is not absence of effect.
  evidence_links      call -> evidence row, with the axis, source, assertion and lineage group of the row (the rows themselves stay in the index)
  gene_context        gene-level evidence for the genes that matched calls sit in, labelled gene_level_context_only: a gene's classification is never copied
                      onto the variants inside it
  sv_calls            every structural / copy-number call with equivalent events, refused near-misses (type, breakpoints, dosage) and context
  haplotype_calls     star-allele / haplotype calls per gene with the evidence they need and the alleles that could not be evaluated
  repeat_calls        repeat-expansion loci: genotyped only when a repeat genotyper reported them; repeat genotypes at loci the index does not know are kept as NO_INDEX_RECORD
  unresolved_loci     high-value indexed variants the sample cannot answer for (absent from a variant-only file, no-call, filtered, inside a structural call...),
                      each with the exact evidence that would settle it
  rejected_calls      alleles the reader refused (REF mismatch, unmapped contig...) with the reason
  run_manifest.json   input digest, assembly check, index digests, counts, versions, the commit and the SHA-256 of every module that produced the tables

Personal calls are read from a local path and written to a local directory; nothing is transmitted and nothing here is meant to be committed.
"""
import hashlib
import json
import os
import subprocess
import sys
import time
from collections import Counter

import pyarrow as pa
import pyarrow.parquet as pq

from . import SCHEMA_VERSION, schema, views
from .callability import CallableRegionCallability
from .index_reader import Index
from .match_hap import NON_STAR_SYSTEMS, call_gene
from .match_repeat import loci_from_index, match_repeats, unindexed_repeats
from .match_small import MatchContext, iter_interpretations, missing_evidence, status_of, zygosity
from .match_sv import match_svs
from .vcfio import OTHER_CONTIGS, CallSet, QCPolicy, read_vcf

S, LS, I32 = pa.string(), pa.list_(pa.string()), pa.int32()
SCHEMAS = {
    'small_calls': pa.schema([('key', S), ('chrom', S), ('pos', I32), ('ref', S), ('alt', S), ('vtype', S), ('gt', S), ('ploidy', I32), ('dosage', I32), ('zygosity', S), ('phased', pa.bool_()),
                              ('phase_set', S), ('qual', pa.float64()), ('filter', S), ('gq', I32), ('dp', I32), ('ad_ref', I32), ('ad_alt', I32), ('qc_pass', pa.bool_()), ('qc_flags', LS),
                              ('status', S), ('flags', LS), ('multiallelic', pa.bool_()), ('n_evidence', I32), ('axes_json', S), ('evidence_ids', LS), ('gene_ids', LS), ('gene_symbols', LS),
                              ('region_kinds', LS), ('mnv_json', S)]),
    'evidence_links': pa.schema([('call_id', S), ('call_kind', S), ('evidence_id', S), ('level', S), ('axis', S), ('source_id', S), ('claim_type', S), ('effect_assertion', S), ('lineage_group', S),
                                 ('relation', S)]),
    'gene_context': pa.schema([('gene_id', S), ('symbol', S), ('n_matched_calls', I32), ('evidence_id', S), ('axis', S), ('source_id', S), ('claim_type', S), ('claim_value', S),
                               ('effect_assertion', S), ('lineage_group', S), ('applies_as', S)]),
    'sv_calls': pa.schema([('call_id', S), ('chrom', S), ('start', I32), ('end', I32), ('svtype', S), ('svlen', I32), ('cn', I32), ('gt', S), ('dosage', I32), ('precise', pa.bool_()),
                           ('qc_pass', pa.bool_()), ('filter', S), ('n_equivalent', I32), ('equivalent_event_ids', LS), ('refused_json', S), ('context_json', S), ('unresolved', LS),
                           ('event_evidence_ids', LS)]),
    'haplotype_calls': pa.schema([('gene', S), ('system', S), ('status', S), ('diplotypes_json', S), ('basis', S), ('missing_json', S), ('flags', LS), ('unevaluable_json', S),
                                  ('position_states_json', S)]),
    'repeat_calls': pa.schema([('locus_id', S), ('gene', S), ('status', S), ('classes', LS), ('counts', pa.list_(I32)), ('thresholds_json', S), ('flags', LS), ('missing', S),
                               ('inheritance', S), ('disease', S), ('chrom', S), ('start', I32), ('end', I32), ('motif', S), ('caller', S)]),
    'unresolved_loci': pa.schema([('key', S), ('chrom', S), ('pos', I32), ('status', S), ('reasons', LS), ('missing_evidence', S), ('rules', LS), ('evidence_ids', LS), ('basis', S),
                                  ('detail_json', S)]),
    'rejected_calls': pa.schema([('chrom', S), ('pos', I32), ('ref', S), ('alt', S), ('reason', S), ('detail', S), ('rec_index', I32)]),
}
RESOLVED = ('EXACT_MATCH', 'VERIFIED_REFERENCE', 'OTHER_ALLELE_AT_SITE', 'DELETED_BY_UPSTREAM')


class AssemblyMismatch(ValueError):
    pass


class TableWriter:
    def __init__(self, path, schema, chunk=50_000):
        self.path, self.schema, self.chunk = path, schema, chunk
        self.rows, self.w, self.n = [], None, 0

    def add(self, row):
        self.rows.append(row)
        if len(self.rows) >= self.chunk:
            self.flush()

    def flush(self):
        if not self.rows:
            return
        t = pa.Table.from_pylist([{n: r.get(n) for n in self.schema.names} for r in self.rows], schema=self.schema)
        if self.w is None:
            self.w = pq.ParquetWriter(self.path, self.schema, compression='zstd')
        self.w.write_table(t)
        self.n += len(self.rows)
        self.rows = []

    def close(self):
        self.flush()
        if self.w is None:
            pq.write_table(self.schema.empty_table(), self.path, compression='zstd')
        else:
            self.w.close()


def sha256_file(path, chunk=1 << 22):
    h = hashlib.sha256()
    with open(path, 'rb') as f:
        for b in iter(lambda: f.read(chunk), b''):
            h.update(b)
    return h.hexdigest()


def _i(x):
    try:
        return None if x is None else int(x)
    except (TypeError, ValueError):
        return None


def load_bed_intervals(path, reference):
    """{canonical contig: [(start1, end1)]} from a BED (0-based half-open on disk, 1-based inclusive here); read once and filtered per chromosome."""
    out = {}
    with open(path) as f:
        for line in f:
            if not line.strip() or line.startswith(('#', 'track', 'browser')):
                continue
            c, a, b = line.split('\t')[:3]
            canon = reference.contig(c)
            if canon:
                out.setdefault(canon, []).append((int(a) + 1, int(b)))
    return out


def jdump(x):
    return json.dumps(x, sort_keys=True, default=str) if x not in (None, {}, [], ()) else None


def code_version(repo_root):
    try:
        head = subprocess.run(['git', '-C', repo_root, 'rev-parse', 'HEAD'], capture_output=True, text=True, timeout=10).stdout.strip() or None
        dirty = bool(subprocess.run(['git', '-C', repo_root, 'status', '--porcelain', '--', 'variant_index'], capture_output=True, text=True, timeout=30).stdout.strip())
        return dict(commit=head, uncommitted_changes_in_variant_index=dirty)
    except (OSError, subprocess.SubprocessError):
        return dict(commit=None, uncommitted_changes_in_variant_index=None)


def loaded_code():
    """{path relative to variant_index/: SHA-256} of the source of every vindex module (and the entry script) loaded in this process: the code that makes a run's tables, so a later
    audit can tell whether any of it has changed since."""
    root = os.path.realpath(os.path.join(os.path.dirname(__file__), '..'))
    out = {}
    for mod in list(sys.modules.values()):
        f = getattr(mod, '__file__', None)
        if f and f.endswith('.py') and os.path.realpath(f).startswith(root + os.sep):
            rel = os.path.relpath(os.path.realpath(f), root)
            if not rel.startswith('tests' + os.sep):                                                # the tests of a unit-test process are not the matcher
                out[rel] = sha256_file(f)
    return dict(sorted(out.items()))


def check_output_location(out_dir, repo_root):
    """Personal results are written under <repo>/data (git-ignored) or outside the repository, never inside the tracked tree."""
    out, repo = os.path.realpath(out_dir), os.path.realpath(repo_root)
    inside = os.path.commonpath([out, repo]) == repo
    if inside and not (out == os.path.join(repo, 'data') or out.startswith(os.path.join(repo, 'data') + os.sep)):
        raise ValueError(f'{out_dir} is inside the repository outside data/; personal results must go under data/ (git-ignored) or outside the repository')


def small_row(o):
    c = o.call
    ev = sorted(r['evidence_id'] for r in o.evidence)
    axes = {a: {k: v for k, v in d.items() if k != 'evidence_ids'} for a, d in o.axes.items()}
    return dict(key=c.key, chrom=c.chrom, pos=c.pos, ref=c.ref, alt=c.alt, vtype=c.vtype, gt=c.gt, ploidy=c.ploidy, dosage=c.dosage, zygosity=zygosity(c.dosage, c.ploidy), phased=c.phased,
                phase_set=None if c.phase_set is None else str(c.phase_set), qual=c.qual, filter=c.filter, gq=_i(c.gq), dp=_i(c.dp), ad_ref=_i(c.ad_ref), ad_alt=_i(c.ad_alt), qc_pass=c.qc_pass,
                qc_flags=list(c.qc_flags), status=o.status, flags=list(o.flags), multiallelic=c.multiallelic, n_evidence=len(ev), axes_json=jdump(axes), evidence_ids=ev or None,
                gene_ids=sorted({g['gene_id'] for g in o.genes}) or None, gene_symbols=sorted({g['symbol'] for g in o.genes if g.get('symbol')}) or None,
                region_kinds=sorted({r['kind'] for r in o.regions}) or None,
                mnv_json=jdump([dict(mnv=k, status=s.status, reasons=list(s.reasons)) for k, s in o.mnv_context]))


def link_rows(call_id, kind, rows, relation):
    return [dict(call_id=call_id, call_kind=kind, evidence_id=r['evidence_id'], level=r['level'], axis=r['axis'], source_id=r['source_id'], claim_type=r['claim_type'],
                 effect_assertion=r['effect_assertion'], lineage_group=r['lineage_group'], relation=relation) for r in rows]


def genes_by_chrom(index):
    """{chromosome: [gene symbols]} for the star-allele systems in the index, placed by the chromosome of their defining variants; '' for genes with none."""
    haps = [h for h in index.haplotypes() if h['system'] not in NON_STAR_SYSTEMS and h.get('gene_symbol')]
    defs = index.haplotype_definitions([h['haplotype_id'] for h in haps])
    where = {}
    for h in haps:
        for d in defs.get(h['haplotype_id'], ()):
            if d.get('variant_key'):
                where.setdefault(h['gene_symbol'], Counter())[d['variant_key'].split(':', 1)[0]] += 1
    out = {}
    for g in sorted({h['gene_symbol'] for h in haps}):
        out.setdefault(where[g].most_common(1)[0][0] if g in where else '', []).append(g)
    return out


OTHER_PASS = 'other_contigs'      # label of the pass over every record off the primary assembly (alternate loci, unplaced scaffolds, decoys, contigs the reference does not know)


def chromosome_order(reference, wanted=None):
    prim = [c for c in reference.lengths if reference.is_primary(c)]
    key = lambda c: (0, int(c)) if c.isdigit() else (1, {'X': 0, 'Y': 1, 'MT': 2}.get(c, 3))
    return [c for c in sorted(prim, key=key) if wanted is None or c in wanted]


def run_sample(vcf, index_dir, out_dir, reference, sample=None, qc=None, contigs=None, callable_bed=None, hash_input=True, external_haplotypes=None, rules=views.HIGH_VALUE_RULES,
               repo_root=None, allow_assembly_mismatch=False, progress=None):
    repo_root = repo_root or os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..'))
    check_output_location(out_dir, repo_root)
    os.makedirs(out_dir, exist_ok=True)
    t0 = time.time()
    code, files_at_start = code_version(repo_root), loaded_code()                                   # the code is loaded now: the tree at the end of a long run may differ
    qc = qc or QCPolicy()
    index = Index.open(index_dir)
    W = {n: TableWriter(os.path.join(out_dir, f'{n}.parquet'), s) for n, s in SCHEMAS.items()}
    counts = dict(by_status=Counter(), by_vtype=Counter(), matched_by_axis=Counter(), matched_by_source=Counter(), difficult_region_calls=0, calls=0, per_chromosome={}, vcf_stats=Counter(),
                  index_centric={r.name: Counter() for r in rules}, sv=Counter(), haplotypes=Counter(), repeats=Counter())
    gene_hits = Counter()
    assembly, sample_name = None, sample
    placed = genes_by_chrom(index)
    loci = loci_from_index(index) if index.regions_by_kind(['strchive_locus']) else []
    order = chromosome_order(reference, contigs)
    primary = set(chromosome_order(reference))
    passes = [(c, [c]) for c in order] + ([(OTHER_PASS, [OTHER_CONTIGS])] if contigs is None else [])        # a whole-file run also reads what lies off the primary assembly
    bed = load_bed_intervals(callable_bed, reference) if callable_bed else None
    for chrom, wanted in passes:
        cs = read_vcf(vcf, reference, sample, qc, contigs=wanted)
        if assembly is None:
            assembly, sample_name = cs.assembly, cs.sample
            if not assembly.get('ok', True) and not allow_assembly_mismatch:
                raise AssemblyMismatch(f'the VCF header does not match the reference ({assembly}); matching a different assembly would give wrong answers')
        counts['vcf_stats'].update({k: v for k, v in cs.stats.items() if isinstance(v, int)})
        providers = [CallableRegionCallability(((c, a, b) for c, ivs in bed.items() if (c not in primary if chrom == OTHER_PASS else c == chrom) for a, b in ivs), cs, 'callable_region')] if bed is not None else []
        ctx = MatchContext(cs, providers)
        n_chrom = 0
        for o in iter_interpretations(cs, index, ctx):
            row = small_row(o)
            W['small_calls'].add(row)
            for r in link_rows(o.call.key, 'small', o.evidence, 'exact_variant'):
                W['evidence_links'].add(r)
            n_chrom += 1
            counts['by_status'][o.status] += 1
            counts['by_vtype'][o.call.vtype] += 1
            counts['difficult_region_calls'] += 'difficult_region' in o.flags
            for a in o.axes:
                counts['matched_by_axis'][a] += 1
            for s in {r['source_id'] for r in o.evidence}:
                counts['matched_by_source'][s] += 1
            if o.evidence and o.status == 'EXACT_MATCH':          # presence asserted: a filtered call or a no-call does not make its gene a matched gene
                for g in {x['gene_id'] for x in o.genes}:
                    gene_hits[g] += 1
        counts['calls'] += n_chrom
        counts['per_chromosome'][chrom] = dict(records=cs.stats.get('records', 0), small_calls=n_chrom, svs=len(cs.svs), repeats=len(cs.repeats), refblocks=len(cs.refblocks))
        for rj in cs.rejected + cs.unplaced:
            W['rejected_calls'].add(dict(chrom=rj.contig, pos=rj.pos, ref=rj.ref, alt=rj.alt, reason=rj.reason, detail=rj.detail, rec_index=rj.rec_index))
        # structural and copy-number calls
        interpreted = match_svs(cs, index)
        eq_ids = sorted({m.event_id for si in interpreted for m in si.matches if m.equivalent})
        event_evidence = index.evidence_for('event', eq_ids) if eq_ids else {}                   # one read for the chromosome: a read per call scans the event evidence 25,000 times
        for si in interpreted:
            sv = si.sv
            cid = f'sv:{sv.chrom}:{sv.start}-{sv.end}:{sv.svtype}:{sv.rec_index}'
            eq = [m for m in si.matches if m.equivalent]
            refused = [dict(event_id=m.event_id, relation=m.relation, reasons=list(m.reasons), metrics=m.metrics) for m in si.matches if not m.equivalent]
            ev_rows = [r for m in eq for r in event_evidence.get(m.event_id, ())]
            for r in link_rows(cid, 'sv', ev_rows, 'equivalent_event'):
                W['evidence_links'].add(r)
            W['sv_calls'].add(dict(call_id=cid, chrom=sv.chrom, start=sv.start, end=sv.end, svtype=sv.svtype, svlen=_i(sv.svlen), cn=_i(sv.cn), gt=sv.gt, dosage=_i(sv.dosage), precise=bool(sv.precise),
                                   qc_pass=bool(sv.qc_pass), filter=sv.filter, n_equivalent=len(eq), equivalent_event_ids=[m.event_id for m in eq] or None, refused_json=jdump(refused),
                                   context_json=jdump(si.context), unresolved=list(si.unresolved) or None, event_evidence_ids=sorted(r['evidence_id'] for r in ev_rows) or None))
            counts['sv'][('equivalent' if eq else 'no_equivalent_event')] += 1
        # star-allele / haplotype calls for the genes on this chromosome
        for gene in placed.get(chrom, ()):
            counts['haplotypes'][_hap_row(W, call_gene(ctx, index, gene, external=external_haplotypes))] += 1
        # repeat loci on this chromosome: every indexed locus, and every repeat genotype in the file whether or not the index knows its locus
        here = [L for L in loci if L['chrom'] == chrom]
        for r in match_repeats(cs, here) + unindexed_repeats(cs, here):
            W['repeat_calls'].add(dict(locus_id=r.locus_id, gene=r.gene or None, status=r.status, classes=list(r.classes), counts=[c for c in r.counts if c is not None] or None, thresholds_json=jdump(r.thresholds),
                                       flags=list(r.flags), missing=r.missing or None, inheritance=r.inheritance or None, disease=r.disease or None, chrom=r.chrom or None, start=r.start, end=r.end,
                                       motif=r.motif or None, caller=r.caller or None))
            counts['repeats'][r.status] += 1
        # index-centric: what does the sample say about the high-value indexed alleles on this chromosome?
        for key, hv in sorted(views.high_value_keys(index, schema.OTHER_PARTITION if chrom == OTHER_PASS else chrom, rules).items(), key=lambda kv: (kv[0].split(':')[0], int(kv[0].split(':')[1]))):
            st = status_of(ctx, key)
            for rn in hv['rules']:
                counts['index_centric'][rn][st.status] += 1
            if st.status not in RESOLVED:
                W['unresolved_loci'].add(dict(key=key, chrom=key.split(':')[0], pos=int(key.split(':')[1]), status=st.status, reasons=list(st.reasons), rules=hv['rules'], evidence_ids=hv['evidence_ids'],
                                              missing_evidence=missing_evidence(st), basis=st.basis, detail_json=jdump({k: v for k, v in (st.detail or {}).items() if k != 'missing_evidence'})))
        if progress:
            progress(chrom, counts['per_chromosome'][chrom])
        del cs, ctx
    # sex-chromosome-free and unplaced genes: the haplotype caller states what is missing without calls
    empty = CallSet(sample=sample_name or '', path=vcf, assembly=assembly or {})
    for gene in placed.get('', ()):
        counts['haplotypes'][_hap_row(W, call_gene(MatchContext(empty), index, gene, external=external_haplotypes))] += 1
    # gene-level context for the genes that matched calls sit in
    if gene_hits:
        ev = index.evidence_for('gene', list(gene_hits))
        symbols = {g['gene_id']: g['symbol'] for g in index.filtered('entity_gene', pa_isin('gene_id', list(gene_hits))).to_pylist()}
        for gid in sorted(ev):
            for r in ev[gid]:
                W['gene_context'].add(dict(gene_id=gid, symbol=symbols.get(gid), n_matched_calls=gene_hits[gid], evidence_id=r['evidence_id'], axis=r['axis'], source_id=r['source_id'],
                                           claim_type=r['claim_type'], claim_value=r['claim_value'], effect_assertion=r['effect_assertion'], lineage_group=r['lineage_group'],
                                           applies_as='gene_level_context_only'))
    for w in W.values():
        w.close()
    files_at_end = loaded_code()
    changed = code_version(repo_root) != code or any(files_at_end.get(k) != v for k, v in files_at_start.items())
    manifest = dict(
        schema_version=SCHEMA_VERSION, finished_utc=time.strftime('%Y-%m-%dT%H:%M:%SZ', time.gmtime()), seconds=round(time.time() - t0, 1),
        input=dict(path=os.path.abspath(vcf), size_bytes=os.path.getsize(vcf), sha256=sha256_file(vcf) if hash_input else None, sample=sample_name, assembly=assembly,
                   callable_bed=os.path.abspath(callable_bed) if callable_bed else None, contigs=[c for c, _ in passes]),
        qc_policy=vars(qc), index=index_provenance(index_dir), code=code, code_files=files_at_end, code_changed_during_run=changed,
        rows={n: w.n for n, w in W.items()},
        counts=dict(calls=counts['calls'], by_status=dict(counts['by_status']), by_vtype=dict(counts['by_vtype']), matched_by_axis=dict(counts['matched_by_axis']),
                    matched_by_source=dict(counts['matched_by_source']), difficult_region_calls=counts['difficult_region_calls'], vcf=dict(counts['vcf_stats']),
                    per_chromosome=counts['per_chromosome'], structural=dict(counts['sv']), haplotype_status=dict(counts['haplotypes']), repeat_status=dict(counts['repeats']),
                    high_value_rules={n: dict(c) for n, c in counts['index_centric'].items()}),
        rules={r.name: r.description for r in rules})
    with open(os.path.join(out_dir, 'run_manifest.json'), 'w') as f:
        json.dump(manifest, f, indent=1, sort_keys=True, default=str)
    return manifest


def pa_isin(field, values):
    import pyarrow.compute as pc
    return pc.field(field).isin(values)


def _hap_row(W, r):
    W['haplotype_calls'].add(dict(gene=r.gene, system=r.system, status=r.status, diplotypes_json=jdump([list(d) for d in r.diplotypes]), basis=r.basis, missing_json=jdump(r.missing),
                                  flags=list(r.flags), unevaluable_json=jdump(r.unevaluable_alleles), position_states_json=jdump({k: list(v) for k, v in r.position_states.items()})))
    return r.status


def index_provenance(index_dir):
    out = dict(path=os.path.abspath(index_dir))
    for name in ('index_manifest.json', 'table_digests.json'):
        p = os.path.join(index_dir, name)
        if os.path.exists(p):
            out[name + '_sha256'] = sha256_file(p)
    p = os.path.join(index_dir, 'index_manifest.json')
    if os.path.exists(p):
        with open(p) as f:
            m = json.load(f)
        out.update(schema_version=m.get('schema_version'), built_utc=m.get('built_utc'), releases=m.get('releases'), adapters=[a['adapter'] for a in m.get('adapters', [])],
                   problems=m.get('problems'))
    return out
