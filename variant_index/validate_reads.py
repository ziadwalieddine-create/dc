#!/usr/bin/env python3
"""Check the read-evidence resolver (vindex/reads.py) against a benchmark on real alignments.

    python3 variant_index/validate_reads.py --calls data/results/HG002/small_calls.parquet --truth-vcf BENCHMARK.vcf.gz --truth-bed BENCHMARK.bed \\
        --reads https://.../HG002.GRCh38.60x.1.bam --reads-index HG002.GRCh38.60x.1.bam.bai --out data/results/HG002_reads/full [--fraction 0.5] [--processes 4] [--seed 20261028]
    python3 variant_index/validate_reads.py merge full=data/results/HG002_reads/full half=... quarter=... --record variant_index/reads_validation.json

The truth is a benchmark sample's call set (public GIAB HG002 by default in the project's runs) and the alignments are that sample's reads. Two sets of loci are asked of the reads:

  positives   benchmark SNVs, insertions and deletions inside the benchmark's confident regions (the normalised keys of the matcher's small_calls table), by variant type and
              zygosity: the truth is that the sample carries the allele
  negatives   random positions inside the confident regions, at least NEIGHBOUR bases from any benchmark record: the truth is that the sample does not carry the allele asked for (an
              alternate base chosen at random, or a one-base deletion or insertion of a random base, left-aligned by the project's normaliser)

What it measures is the resolver's rates by truth class and state: how often the reads state the truth (READS_ALTERNATE_PRESENT at a positive, READS_REFERENCE_ONLY at a negative), how
often they abstain (ambiguous, insufficient, not evaluated), and how often they state the opposite (a claim that contradicts the benchmark). It does not make the resolver a genotyper
and it says nothing outside the benchmark's confident regions, where a benchmark can be wrong, reads are harder to place, and the resolver abstains more. A negative is 'absent' as far as
the benchmark is complete inside its regions (about 99.5% of the variants of the sample by the benchmark's own estimate): a contradiction is examined, not assumed to be the resolver's.

Reads are read through pysam, over HTTPS when --reads is a URL (the index is then given with --reads-index). --fraction keeps that share of the fragments, chosen by a hash of the read
name so that both mates stay together and a rerun keeps the same ones: 0.5 of a 60x alignment is a 30x genome, 0.25 a 15x one, which shows how the rates move with depth. The results are
public benchmark data and are written under data/ (git-ignored); nothing personal is involved. Re-running with the same arguments and the same alignment gives the same loci and the same
states.
"""
import argparse
import itertools
import json
import math
import os
import random
import sys
import time
import zlib
from collections import Counter, defaultdict
from multiprocessing import Pool

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

import pyarrow as pa                                                # noqa: E402
import pyarrow.compute as pc                                        # noqa: E402
import pyarrow.parquet as pq                                        # noqa: E402
import pysam                                                        # noqa: E402

from vindex import reads as R                                       # noqa: E402
from vindex.build.common import DATA                                # noqa: E402
from vindex.callability import IntervalSet                          # noqa: E402
from vindex.normalise import normalise                              # noqa: E402
from vindex.refgenome import Reference                              # noqa: E402
from vindex.wgs import load_bed_intervals, sha256_file              # noqa: E402

FASTA = os.path.join(DATA, 'reference_sequence', 'Homo_sapiens.GRCh38.dna.toplevel.fa.gz')
REPORT = os.path.join(DATA, 'refseq', 'GCF_000001405.40_GRCh38.p14_assembly_report.txt')
POSITIVES = {('snv', 'heterozygous'): 400, ('snv', 'homozygous'): 200, ('ins', 'heterozygous'): 200, ('ins', 'homozygous'): 100, ('del', 'heterozygous'): 200, ('del', 'homozygous'): 100}
NEGATIVES = {'snv': 600, 'ins': 200, 'del': 200}
NEIGHBOUR = 60
AUTOSOMES = {str(i) for i in range(1, 23)}
CLAIMS = (R.ALTERNATE_PRESENT, R.REFERENCE_ONLY, R.OTHER_ALLELE_PRESENT)
ABSTENTIONS = (R.AMBIGUOUS, R.INSUFFICIENT, R.NOT_EVALUATED)


def wilson(k, n, z=1.96):
    if n == 0:
        return None
    p = k / n
    d = 1 + z * z / n
    centre = (p + z * z / (2 * n)) / d
    half = z * math.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / d
    return [round(max(0.0, centre - half), 4), round(min(1.0, centre + half), 4)]


def sample_positives(calls, mask, rng):
    """Reservoir sample per (type, zygosity) of the benchmark calls that are biallelic, diploid, quality-passing, autosomal and inside the confident regions."""
    pf = pq.ParquetFile(calls)
    seen, kept = Counter(), defaultdict(list)
    cols = ['key', 'chrom', 'pos', 'ref', 'vtype', 'zygosity', 'ploidy', 'qc_pass', 'multiallelic']
    for batch in pf.iter_batches(columns=cols, batch_size=200_000):
        keep = pc.and_(pc.and_(pc.equal(batch['ploidy'], 2), pc.invert(batch['multiallelic'])), pc.and_(batch['qc_pass'], pc.is_in(batch['vtype'], pa.array(['snv', 'ins', 'del']))))
        for r in batch.filter(keep).to_pylist():
            stratum = (r['vtype'], r['zygosity'])
            if stratum not in POSITIVES or r['chrom'] not in AUTOSOMES or not mask.covers(r['chrom'], r['pos'], r['pos'] + len(r['ref']) - 1):
                continue
            seen[stratum] += 1
            if len(kept[stratum]) < POSITIVES[stratum]:
                kept[stratum].append(r['key'])
            else:
                j = rng.randrange(seen[stratum])
                if j < POSITIVES[stratum]:
                    kept[stratum][j] = r['key']
    return {s: sorted(v) for s, v in kept.items()}, dict(seen)


def sample_negatives(reference, intervals, mask, truth_vcf, rng):
    tbx = pysam.TabixFile(truth_vcf)
    vcf_name = {reference.contig(n): n for n in tbx.contigs if reference.contig(n)}
    flat = [(c, a, b) for c in sorted(intervals, key=lambda c: int(c)) if c in AUTOSOMES for a, b in sorted(intervals[c]) if b - a > 2 * NEIGHBOUR]
    cumulative = list(itertools.accumulate(b - a for _, a, b in flat))
    out, tries = {k: [] for k in NEGATIVES}, 0
    used = set()
    while any(len(out[k]) < n for k, n in NEGATIVES.items()):
        tries += 1
        if tries > 5_000_000:
            raise RuntimeError('could not draw the negatives')
        kind = rng.choice([k for k, n in NEGATIVES.items() if len(out[k]) < n])
        c, a, b = rng.choices(flat, cum_weights=cumulative)[0]
        pos = rng.randint(a + NEIGHBOUR, b - NEIGHBOUR)
        if not mask.covers(c, pos - NEIGHBOUR, pos + NEIGHBOUR) or any(True for _ in tbx.fetch(vcf_name[c], pos - NEIGHBOUR - 1, pos + NEIGHBOUR)):
            continue
        seq = reference.fetch(c, pos, pos + 1)
        if 'N' in seq or len(seq) < 2:
            continue
        if kind == 'snv':
            n = normalise(reference, c, pos, seq[0], rng.choice([x for x in 'ACGT' if x != seq[0]]))
        elif kind == 'del':
            n = normalise(reference, c, pos, seq, seq[0])
        else:
            n = normalise(reference, c, pos, seq[0], seq[0] + rng.choice('ACGT'))
        if n.ok and n.key not in used:
            used.add(n.key)
            out[kind].append(n.key)
    return {k: sorted(v) for k, v in out.items()}


_EVIDENCE = None


class Subsample:
    """An alignment file that yields only the fragments whose read-name hash falls in the first `fraction` of the hash range (both mates of a pair share a name)."""

    def __init__(self, af, fraction):
        self.af, self.cut = af, int(fraction * 10_000)

    def fetch(self, *args, **kw):
        return (a for a in self.af.fetch(*args, **kw) if zlib.crc32(a.query_name.encode()) % 10_000 < self.cut)

    def __getattr__(self, name):
        return getattr(self.af, name)


_ARGS = None


def _make():
    global _EVIDENCE
    path, index, fasta, report, policy, fraction = _ARGS
    _EVIDENCE = R.ReadEvidence(path, Reference(fasta, report), policy, index=index)
    if fraction < 1:
        _EVIDENCE.af = Subsample(_EVIDENCE.af, fraction)


def _open(*args):
    global _ARGS
    _ARGS = args
    _make()


def _evaluate(key):
    for attempt in range(4):
        try:
            s = _EVIDENCE.evaluate(key)
            return key, dict(state=s.state, reasons=list(s.reasons), flags=list(s.flags), depth=s.depth, ref=s.ref_reads, alt=s.alt_reads, other=s.other_reads, alt_fraction=s.alt_fraction,
                             window=(s.detail or {}).get('window'), detail=json.dumps(s.detail, sort_keys=True) if s.detail else None)
        except (OSError, ValueError) as e:                            # a dropped connection: retry on a fresh handle (a failed seek can leave the old one unusable); a persistent failure is recorded, never counted as evidence
            err = f'{type(e).__name__}: {e}'
            time.sleep(2 * (attempt + 1))
            try:
                _EVIDENCE.close()
                _make()
            except (OSError, ValueError) as e2:
                err = f'{type(e2).__name__}: {e2}'
    return key, dict(state='ERROR', reasons=[err[:200]], flags=[], depth=0, ref=0, alt=0, other=0, alt_fraction=None, window=None, detail=None)


def summarise(rows):
    """Counts by truth class and state, the rates with Wilson intervals, the claims that contradict the truth, and the indels by the length of the window they were compared over."""
    by = defaultdict(Counter)
    for r in rows:
        by[r['truth']][r['state']] += 1
    out = dict(by_truth={t: dict(sorted(c.items())) for t, c in sorted(by.items())}, rates={})
    for t, c in sorted(by.items()):
        n = sum(c.values())
        rate = lambda k: dict(n=c[k], rate=round(c[k] / n, 4), ci95=wilson(c[k], n)) if n else None
        stated = c[R.ALTERNATE_PRESENT] if t.startswith('truth_alt') else c[R.REFERENCE_ONLY]
        claims = sum(c[x] for x in CLAIMS)
        out['rates'][t] = dict(loci=n, states={k: rate(k) for k in sorted(c)}, states_the_truth=dict(n=stated, rate=round(stated / n, 4), ci95=wilson(stated, n)),
                               abstains=dict(n=sum(c[x] for x in ABSTENTIONS), rate=round(sum(c[x] for x in ABSTENTIONS) / n, 4), ci95=wilson(sum(c[x] for x in ABSTENTIONS), n)),
                               contradicts_the_truth_among_claims=dict(n=claims - stated, claims=claims, rate=round((claims - stated) / claims, 4) if claims else None, ci95=wilson(claims - stated, claims)))
    ctx = defaultdict(lambda: defaultdict(Counter))
    for r in rows:
        if r['vtype'] in ('ins', 'del') and r['window']:
            width = r['window'][2] - r['window'][1] + 1
            ctx[r['truth_group']]['window ' + ('<= 7' if width <= 7 else '8-12' if width <= 12 else '> 12') + ' bases'][r['state']] += 1
    out['indels_by_window_length'] = {g: {w: dict(sorted(c.items())) for w, c in sorted(v.items())} for g, v in sorted(ctx.items())}
    frac = defaultdict(list)
    for r in rows:
        if r['alt_fraction'] is not None and r['truth'].startswith('truth_alt'):
            frac[r['truth']].append(r['alt_fraction'])
    q = lambda v, p: round(sorted(v)[min(len(v) - 1, int(p * len(v)))], 3)
    out['alternate_fraction_at_positives'] = {t: dict(n=len(v), q10=q(v, .1), median=q(v, .5), q90=q(v, .9)) for t, v in sorted(frac.items())}
    depth = defaultdict(list)
    for r in rows:
        if r['vtype'] == 'snv' and r['state'] != 'ERROR':
            depth[r['truth']].append(r['depth'])
    out['informative_fragments_at_snvs'] = {t: dict(n=len(v), q10=q(v, .1), median=q(v, .5), q90=q(v, .9)) for t, v in sorted(depth.items())}
    return out


def merge(argv):
    """variant_index/reads_validation.json from the summaries of runs at several depths: the rates and counts, the policy, the public inputs and the digest of reads.py the runs used (the gap
    audit refuses a record whose digest is not that of the module as it is now). Local paths are left out; the record is derived from public benchmark data only."""
    ap = argparse.ArgumentParser(prog='validate_reads.py merge')
    ap.add_argument('runs', nargs='+', help='label=output directory of a run')
    ap.add_argument('--record', required=True)
    a = ap.parse_args(argv)
    runs, digests, policies = {}, set(), []
    for spec in a.runs:
        label, _, d = spec.partition('=')
        with open(os.path.join(d, 'summary.json')) as f:
            s = json.load(f)
        digests.add(s['reads_module_sha256'])
        policies.append(s['policy'])
        runs[label] = dict(fraction=s['arguments']['fraction'], seed=s['arguments']['seed'], reads=s['arguments']['reads'], seconds=s['seconds'], by_truth=s['by_truth'], rates=s['rates'],
                           indels_by_window_length=s['indels_by_window_length'], alternate_fraction_at_positives=s['alternate_fraction_at_positives'],
                           informative_fragments_at_snvs=s['informative_fragments_at_snvs'])
    if len(digests) != 1 or any(p != policies[0] for p in policies):
        raise SystemExit('the runs were made with different versions of reads.py or different policies; run them again')
    record = dict(note=('Rates of vindex/reads.py against the GIAB HG002 GRCh38 v4.2.1 benchmark on public Illumina reads (Novoalign, 60x), by depth (fraction of the fragments kept). Positives are benchmark '
                        'SNVs and indels in the confident regions; negatives are random alleles at random positions in them. Produced by validate_reads.py; see its docstring for what it does not show.'),
                   reads_module_sha256=digests.pop(), policy=policies[0], runs=dict(sorted(runs.items())))
    with open(a.record, 'w') as f:
        json.dump(record, f, indent=1, sort_keys=True)
        f.write('\n')
    print(f'wrote {a.record}: runs {", ".join(sorted(runs))}')
    return 0


def main(argv=None):
    argv = sys.argv[1:] if argv is None else argv
    if argv and argv[0] == 'merge':
        return merge(argv[1:])
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--calls', required=True, help='small_calls.parquet of the matcher run on the benchmark call set')
    ap.add_argument('--truth-vcf', required=True, help='the benchmark VCF (tabix indexed): neighbours of the negatives are looked up in it')
    ap.add_argument('--truth-bed', required=True, help="the benchmark's confident regions")
    ap.add_argument('--reads', required=True, help='indexed BAM/CRAM of the benchmark sample, local or an https URL')
    ap.add_argument('--reads-index')
    ap.add_argument('--out', required=True)
    ap.add_argument('--processes', type=int, default=4)
    ap.add_argument('--seed', type=int, default=20261028)
    ap.add_argument('--fraction', type=float, default=1.0, help='share of the fragments to keep (by read-name hash): 0.5 turns a 60x alignment into 30x')
    ap.add_argument('--fasta', default=FASTA)
    ap.add_argument('--assembly-report', default=REPORT)
    a = ap.parse_args(argv)
    os.makedirs(a.out, exist_ok=True)
    ref = Reference(a.fasta, a.assembly_report)
    rng = random.Random(a.seed)
    intervals = load_bed_intervals(a.truth_bed, ref)
    mask = IntervalSet()
    for c, ivs in intervals.items():
        for s, e in ivs:
            mask.add(c, s, e)
    t0 = time.time()
    pos, seen = sample_positives(a.calls, mask, rng)
    neg = sample_negatives(ref, intervals, mask, a.truth_vcf, rng)
    print(f'loci drawn in {time.time() - t0:.0f}s: positives {{{", ".join(f"{k[0]}/{k[1]}: {len(v)} of {seen[k]:,}" for k, v in sorted(pos.items()))}}}, negatives {{{", ".join(f"{k}: {len(v)}" for k, v in neg.items())}}}', flush=True)
    loci = [(key, f"truth_alt_{z[:3]}", z, v) for (v, z), keys in sorted(pos.items()) for key in keys] + [(key, 'truth_absent', 'absent', v) for v, keys in sorted(neg.items()) for key in keys]
    policy = R.ReadPolicy()
    results = {}
    with Pool(a.processes, initializer=_open, initargs=(a.reads, a.reads_index, a.fasta, a.assembly_report, policy, a.fraction)) as pool:
        for i, (key, r) in enumerate(pool.imap_unordered(_evaluate, [k for k, *_ in loci], chunksize=8), 1):
            results[key] = r
            if i % 200 == 0:
                print(f'  {i:,} of {len(loci):,} ({time.time() - t0:.0f}s)', flush=True)
    rows = []
    for key, truth, zyg, vtype in loci:
        r = results[key]
        rows.append(dict(key=key, truth=truth, truth_group=f'{truth} {vtype}', zygosity=zyg, vtype=vtype, **r))
    schema = pa.schema([('key', pa.string()), ('truth', pa.string()), ('truth_group', pa.string()), ('zygosity', pa.string()), ('vtype', pa.string()), ('state', pa.string()),
                        ('reasons', pa.list_(pa.string())), ('flags', pa.list_(pa.string())), ('depth', pa.int32()), ('ref', pa.int32()), ('alt', pa.int32()), ('other', pa.int32()),
                        ('alt_fraction', pa.float64()), ('window_text', pa.string()), ('detail', pa.string())])
    table = [dict(r, window_text=f"{r['window'][0]}:{r['window'][1]}-{r['window'][2]}" if r['window'] else None) for r in rows]
    pq.write_table(pa.Table.from_pylist([{n: r.get(n) for n in schema.names} for r in table], schema=schema), os.path.join(a.out, 'loci.parquet'), compression='zstd')
    summary = dict(arguments=dict(calls=a.calls, truth_vcf=a.truth_vcf, truth_bed=a.truth_bed, reads=a.reads, seed=a.seed, fraction=a.fraction), sampled_from=dict(positives={f'{k[0]}/{k[1]}': v for k, v in sorted(seen.items())}),
                   policy=vars(policy), pysam=pysam.__version__, reads_module_sha256=sha256_file(os.path.join(HERE, 'vindex', 'reads.py')), seconds=round(time.time() - t0, 1), **summarise(rows))
    with open(os.path.join(a.out, 'summary.json'), 'w') as f:
        json.dump(summary, f, indent=1, sort_keys=True)
    print(json.dumps(summary['rates'], indent=1))
    return 0


if __name__ == '__main__':
    sys.exit(main())
