#!/usr/bin/env python3
"""Cross-check the base counts of the read-evidence resolver (vindex/reads.py) against `bcftools mpileup`, an established pileup tool, at the single-nucleotide loci of a reads run.

    python3 variant_index/crosscheck_reads_bcftools.py --results data/results/HG002_readrun --out data/results/HG002_readrun/crosscheck_reads.json
                                                       [--reads ALIGNMENT] [--reads-index BAI] [--limit-per-state N] [--processes 4] [--chunk 80]
                                                       [--record variant_index/reads_crosscheck.json]

The alignment and the policy are the run's own (the run manifest names them); --reads and --reads-index replace the alignment. Every SNV key of the run's unresolved_loci table that the reads
were asked about is looked up again with `bcftools mpileup` under the filters of the resolver's policy: mapping quality and base quality at the policy's minimum, no BAQ, unmapped,
secondary, supplementary, QC-failed and duplicate alignments excluded, pairs not marked as proper pairs kept (-A), no reference needed (the counts are by base). Two pileups are made:

  alignments  overlap detection off (-x): one count per alignment passing the filters. At every locus the total must equal the resolver's informative_alignments; each base must be
              counted at least as often as the resolver counted fragments carrying it (so that where every fragment is one alignment the three counts are equal); and the resolver's own
              bookkeeping must be consistent (informative alignments, fragments and conflicting fragments). There is no tolerance here.
  fragments   overlap detection on (the default): the two mates of a pair that overlap the site are counted once, as the resolver counts a fragment. A locus is `exact` when the reference,
              indexed alternate and other counts are equal. Where they differ, bcftools' own overlap handling is the difference, and the alignments say how much it can be: for the pairs
              (two alignments of one name that both pass the filters and cover the site) htslib (sam.c, tweak_overlap_quality) keeps the base of the better mate where the mates disagree (the
              resolver sets the fragment aside), keeps one base where two that agree are each below the base-quality minimum but add up to it, counts both mates of a pair that is not marked
              as a proper pair, scales by 0.8 or zeroes, by a hash of the read name, the base opposite a deletion in the other mate, and drops a good base whose disagreeing mate is
              below the minimum when 0.8 of it falls below the minimum too. A locus is `explained` when every count lies inside the range its pairs allow, read from the alignments and not from
              the resolver's bookkeeping, and the fragments the resolver set aside as conflicting are exactly the pairs whose mates show different bases of the minimum quality; it is
              `unexplained`, and listed with its pairs, otherwise. The effects are shown on synthetic pairs in tests/test_crosscheck_reads.py.

What the fragment comparison cannot show: the range is an allowance, not a prediction, so a wrong resolver count that stays inside the allowance of a locus that has such pairs is not
detected there (the number of conflicting fragments is the exception: it is compared exactly at every locus whose counts differ). The merging of mates is also tested on synthetic pairs (tests/test_reads.py) and against the truth set (reads_validation.json); the alignment comparison has no such freedom.

Indels are not compared: mpileup counts the reads that carry an indel at a position, not the sequence over a window, so the two do not answer the same question. --record is for public
alignments only (it names the alignment and the run's code commit); for a person's reads write --out under data/results, which Git ignores, as the per-locus table beside it holds
coordinates and counts. No sequence is written.
"""
import argparse
import json
import os
import random
import subprocess
import sys
import tempfile
from collections import Counter, defaultdict
from multiprocessing import Pool

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

import pyarrow.parquet as pq                                        # noqa: E402

from vindex import reads as R                                       # noqa: E402
from vindex.build.common import DATA                                # noqa: E402
from vindex.crosscheck import find_bcftools                         # noqa: E402
from vindex.refgenome import Reference                              # noqa: E402
from vindex.wgs import sha256_file                                  # noqa: E402

FASTA = os.path.join(DATA, 'reference_sequence', 'Homo_sapiens.GRCh38.dna.toplevel.fa.gz')
REPORT = os.path.join(DATA, 'refseq', 'GCF_000001405.40_GRCh38.p14_assembly_report.txt')
EXCLUDED = 'UNMAP,SECONDARY,QCFAIL,DUP,SUPPLEMENTARY'                # the alignments the resolver sets aside
FATAL = ('unrecognized option', 'Invalid option', 'Usage:', 'No such file', 'could not parse')
MISSING = object()
_CTX = {}


def _init(reads, index, bcftools, mapq, baseq):
    _CTX.update(reads=reads, index=index, bcftools=bcftools, mapq=mapq, baseq=baseq)


def parse_pileup(text, wanted):
    """{(contig, pos): ({base: n}, total, dp_equals_total)} from the VCF text of one `bcftools mpileup`; a position that is not listed has no read."""
    out = {}
    for line in text.splitlines():
        if line.startswith('#'):
            continue
        f = line.split('\t')
        site = (f[0], int(f[1]))
        if site not in wanted:
            continue
        fmt = dict(zip(f[8].split(':'), f[9].split(':')))
        ad = [int(x) for x in fmt['AD'].split(',')]
        counts = {a: n for a, n in zip(['N'] + f[4].split(','), ad) if n and a != '<*>'}
        out[site] = (counts, sum(ad), int(fmt['DP']) == sum(ad))
    return out


def pileup(positions, overlaps):
    """(parse_pileup of one call over the sorted (contig, pos) list, errors); a failing call is retried and then split, down to a single position."""
    c = _CTX
    target = c['reads'] + (f"##idx##{c['index']}" if c['index'] else '')
    with tempfile.NamedTemporaryFile('w', suffix='.tsv') as regions:
        regions.write(''.join(f'{contig}\t{pos}\n' for contig, pos in positions))
        regions.flush()
        cmd = [c['bcftools'], 'mpileup', '--no-reference', '-R', regions.name, '-q', str(c['mapq']), '-Q', str(c['baseq']), '-B', '-I', '-A', '-d', '1000000', '--ns', EXCLUDED,
               '-a', 'FORMAT/AD,FORMAT/DP', '--no-version'] + ([] if overlaps else ['-x']) + [target]
        for _ in range(3):
            p = subprocess.run(cmd, capture_output=True, text=True)
            if p.returncode == 0:
                return parse_pileup(p.stdout, set(positions)), []
            if any(m in p.stderr for m in FATAL):
                raise RuntimeError(f"bcftools mpileup cannot run: {p.stderr.strip()[-300:]}")
    if len(positions) == 1:
        return {}, [dict(site=list(positions[0]), error=p.stderr.strip()[-300:])]
    half = len(positions) // 2
    (a, ea), (b, eb) = pileup(positions[:half], overlaps), pileup(positions[half:], overlaps)
    return {**a, **b}, ea + eb


def work(task):
    overlaps, positions = task
    return (overlaps, positions, *pileup(positions, overlaps))


def alignments_verdict(row, detail, base, got):
    """The list of differences (empty: none) between the pileup of single alignments and the resolver's counts and bookkeeping at one locus."""
    counts, total, _ = got
    ref, alt = base
    mine = dict(ref=row['read_ref'], alt=row['read_alt'], other=row['read_other'])
    theirs = dict(ref=counts.get(ref, 0), alt=counts.get(alt, 0))
    theirs['other'] = total - theirs['ref'] - theirs['alt']
    informative, conflicts = detail.get('informative_alignments', 0), detail.get('conflicting_fragments', 0)
    problems = []
    if total != informative:
        problems.append(f'bcftools counts {total} alignments, the resolver {informative}')
    problems += [f'{k}: {theirs[k]} alignments but {mine[k]} fragments' for k in mine if theirs[k] < mine[k]]
    pairs_agreeing = informative - row['read_depth'] - 2 * conflicts                # each conflicting fragment has two informative alignments, each agreeing pair two for one fragment
    if pairs_agreeing < 0:
        problems.append(f'the resolver reports {informative} informative alignments for {row["read_depth"]} fragments and {conflicts} conflicting ones')
    return problems


def pair_patterns(af, contig, pos, min_mapq, min_baseq):
    """(up, down, up_effects, down_effects, kinds) at the 1-based site from the alignments alone: how many more fragments bcftools may count there than the resolver, how many fewer, the
    effects behind each, and the kinds of pair found (two alignments of a name that pass bcftools' read filters and cover the site). The resolver's counts are not used."""
    seen = defaultdict(list)
    for a in af.fetch(contig, pos - 1, pos):
        if a.is_unmapped or a.is_secondary or a.is_supplementary or a.is_qcfail or a.is_duplicate or a.mapping_quality < min_mapq:
            continue
        q = next((q for q, r in a.get_aligned_pairs() if r == pos - 1), MISSING)
        if q is MISSING:
            continue                                                                 # does not cover the site: clipped, or ends before it
        quals = a.query_qualities
        seen[a.query_name].append((None, 0, a.is_proper_pair) if q is None else (a.query_sequence[q].upper(), int(quals[q]) if quals is not None else 0, a.is_proper_pair))
    kinds, up, down = Counter(), Counter(), Counter()
    for mates in seen.values():
        if len(mates) < 2:
            continue
        if len(mates) > 2:                                                           # not a pair: allow what two alignments more than a fragment can make of it
            kinds['more_than_two_alignments_of_a_name'] += 1
            up['more_than_two_alignments_of_a_name'] += len(mates) - 1
            down['more_than_two_alignments_of_a_name'] += len(mates) - 1
            continue
        (b1, q1, p1), (b2, q2, p2) = mates
        proper = p1 and p2
        good = [q >= min_baseq for q in (q1, q2)]
        if b1 is None or b2 is None:
            kind = 'both_have_a_gap' if b1 is None and b2 is None else 'gap_mate'
            if kind == 'gap_mate' and proper and (q1 if b2 is None else q2) >= min_baseq:
                down['mate_has_a_gap_over_the_site'] += 1                            # the covering mate's base is zeroed or scaled by 0.8
        elif b1 == b2:
            kind = 'agree_both_low' if not any(good) else 'agree'
            if kind == 'agree_both_low' and proper and q1 + q2 >= min_baseq:
                up['low_quality_pair_agrees'] += 1                                   # the two qualities add up to the minimum
            if kind == 'agree' and not proper and all(good):
                up['agreeing_pair_not_marked_proper'] += 1                           # not de-overlapped: both mates counted
        elif all(good):
            kind = 'disagree'
            up['mates_disagree' if proper else 'disagreeing_pair_not_marked_proper'] += 1 if proper else 2
        elif any(good):
            kind = 'disagree_low_mate'
            if proper:
                down['low_quality_mate_disagrees'] += 1                              # 0.8 of the good base may fall below the minimum
        else:
            kind = 'disagree_both_low'
        kinds[kind if proper else kind + '/not_proper'] += 1
    return sum(up.values()), sum(down.values()), dict(sorted(up.items())), dict(sorted(down.items())), dict(sorted(kinds.items()))


def fragments_verdict(row, base, got, pairs, conflicts):
    """(verdict, resolver counts, bcftools counts, differences, effects that explain them, kinds of pair found); `pairs` is called only when the counts differ, and gives pair_patterns.
    Where they differ, the fragments the resolver set aside as conflicting (`conflicts`) must be the pairs of the alignments whose two mates show different bases, each of the minimum quality
    (a pair whose mates show two different third bases would count here and not for the resolver, which calls both 'another base'; there is none in the run this was made for)."""
    counts, total, _ = got
    ref, alt = base
    mine = dict(ref=row['read_ref'], alt=row['read_alt'], total=row['read_depth'])
    theirs = dict(ref=counts.get(ref, 0), alt=counts.get(alt, 0), total=total)
    diff = {k: theirs[k] - mine[k] for k in mine}
    diff['other'] = diff['total'] - diff['ref'] - diff['alt']
    if not any(diff.values()):
        return 'exact', mine, theirs, diff, [], {}
    up, down, up_effects, down_effects, kinds = pairs()
    if kinds.get('disagree', 0) + kinds.get('disagree/not_proper', 0) != conflicts:
        return 'unexplained', mine, theirs, diff, [], kinds
    if all(-down <= v <= up for v in diff.values()):
        used = (list(up_effects) if any(v > 0 for v in diff.values()) else []) + (list(down_effects) if any(v < 0 for v in diff.values()) else [])
        return 'explained', mine, theirs, diff, used, kinds
    return 'unexplained', mine, theirs, diff, [], kinds


def bcftools_version(bcftools):
    return subprocess.run([bcftools, '--version'], capture_output=True, text=True).stdout.splitlines()[0]


def snv_rows(table, limit_per_state=0, seed=0):
    """The rows of the unresolved-loci table that are SNV keys the reads were asked about, sorted by key; at most limit_per_state of each read state when asked."""
    rows = []
    for r in pq.read_table(table).to_pylist():
        if r['read_state'] and r['read_state'] != R.NOT_EVALUATED and all(len(x) == 1 for x in R.parse_key(r['key'])[2:]):
            rows.append(r)
    rows.sort(key=lambda r: r['key'])
    if not limit_per_state:
        return rows
    rng, by_state, out = random.Random(seed), defaultdict(list), []
    for r in rows:
        by_state[r['read_state']].append(r)
    for s in sorted(by_state):
        out += rng.sample(by_state[s], min(len(by_state[s]), limit_per_state))
    return sorted(out, key=lambda r: r['key'])


def compare(rows, reads, index, policy, reference, bcftools, processes=4, chunk=80, progress=None):
    """The summary dict and the per-locus lines of the comparison of `rows` (unresolved-loci rows of SNV keys with read evidence) with bcftools mpileup over `reads`."""
    ev = R.ReadEvidence(reads, reference, index=index)
    contig, order = dict(ev._contig), {n: i for i, n in enumerate(ev.af.references)}
    site = {r['key']: (contig[R.parse_key(r['key'])[0]], R.parse_key(r['key'])[1]) for r in rows if R.parse_key(r['key'])[0] in contig}
    positions = sorted(set(site.values()), key=lambda s: (order[s[0]], s[1]))
    tasks = [(overlaps, positions[i:i + chunk]) for i in range(0, len(positions), chunk) for overlaps in (False, True)]
    init, initargs = _init, (reads, index, bcftools, policy['min_mapq'], policy['min_baseq'])
    got, errors = {False: {}, True: {}}, []
    if processes > 1:
        pool = Pool(processes, initializer=init, initargs=initargs)
        results = pool.imap_unordered(work, tasks)
    else:
        init(*initargs)
        pool, results = None, map(work, tasks)
    try:
        for n, (overlaps, part, pile, errs) in enumerate(results, 1):
            for s in part:
                got[overlaps][s] = pile.get(s, ({}, 0, True))
            errors += [dict(e, overlaps_detected=overlaps) for e in errs]
            if progress:
                progress(f'{n} of {len(tasks)} calls')
    finally:
        if pool:
            pool.terminate()
    try:
        return _verdicts(rows, site, got, errors, ev, policy, bcftools)
    finally:
        ev.close()


def _verdicts(rows, site, got, errors, ev, policy, bcftools):
    failed = {tuple(e['site']) for e in errors}
    alignments, fragments, by_state, explainable_by = Counter(), Counter(), defaultdict(Counter), Counter()
    problems, unexplained, explained, not_dp, depth_resolver, depth_bcftools, lines = [], [], [], 0, 0, 0, []
    for r in rows:
        s = site.get(r['key'])
        if s is None or s in failed:
            continue
        detail = json.loads(r['read_detail_json']) if r['read_detail_json'] else {}
        base = R.parse_key(r['key'])[2:]
        pa = alignments_verdict(r, detail, base, got[False][s])
        verdict, mine, theirs, diff, used, kinds = fragments_verdict(r, base, got[True][s], lambda s=s: pair_patterns(ev.af, s[0], s[1], policy['min_mapq'], policy['min_baseq']),
                                                                  detail.get('conflicting_fragments', 0))
        alignments['exact' if not pa else 'unexplained'] += 1
        fragments[verdict] += 1
        by_state[r['read_state']][verdict] += 1
        explainable_by.update(used)
        not_dp += (not got[False][s][2]) + (not got[True][s][2])
        depth_resolver, depth_bcftools = depth_resolver + r['read_depth'], depth_bcftools + got[True][s][1]
        if pa:
            problems.append(dict(key=r['key'], state=r['read_state'], problems=pa))
        if verdict != 'exact':
            (unexplained if verdict == 'unexplained' else explained).append(dict(key=r['key'], state=r['read_state'], resolver=mine, bcftools=theirs, differences=diff, pairs=kinds,
                                                                               explainable_by=used))
        lines.append((r['key'], r['read_state'], r['read_ref'], r['read_alt'], r['read_depth'], detail.get('informative_alignments', 0), detail.get('conflicting_fragments', 0),
                      got[True][s][0].get(base[0], 0), got[True][s][0].get(base[1], 0), got[True][s][1], got[False][s][1], verdict, 'exact' if not pa else 'unexplained', ','.join(used)))
    compared = sum(alignments.values())
    summary = dict(loci=compared, snv_keys_of_the_run=len(rows), not_compared=len(rows) - compared, alignments=dict(sorted(alignments.items())), fragments=dict(sorted(fragments.items())),
                   explainable_by=dict(sorted(explainable_by.items())), by_state={k: dict(sorted(v.items())) for k, v in sorted(by_state.items())},
                   fragments_resolver=depth_resolver, fragments_bcftools=depth_bcftools, n_alignment_problems=len(problems), alignment_problems=problems[:50], n_unexplained=len(unexplained),
                   unexplained=unexplained[:50], explained=explained[:50], pileups_where_dp_differs_from_the_allele_counts=not_dp, errors=errors, bcftools=bcftools_version(bcftools),
                   filters=dict(min_mapq=policy['min_mapq'], min_baseq=policy['min_baseq'], baq=False, excluded_flags=EXCLUDED, proper_pair_required=False),
                   ok=bool(compared) and not problems and not unexplained and not errors and not not_dp)
    return summary, lines


HEADER = ('key', 'read_state', 'resolver_ref', 'resolver_alt', 'resolver_fragments', 'informative_alignments', 'conflicting_fragments', 'bcftools_ref', 'bcftools_alt', 'bcftools_fragments',
          'bcftools_alignments', 'fragments_verdict', 'alignments_verdict', 'explainable_by')


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--results', required=True, help='a run of match_wgs.py --reads')
    ap.add_argument('--out', required=True)
    ap.add_argument('--reads', help='the alignment (default: the one the run names)')
    ap.add_argument('--reads-index')
    ap.add_argument('--record', help='write the summary (counts only) here, for the repository')
    ap.add_argument('--limit-per-state', type=int, default=0, help='compare at most this many loci of each read state (default: all)')
    ap.add_argument('--seed', type=int, default=20261028)
    ap.add_argument('--processes', type=int, default=4)
    ap.add_argument('--chunk', type=int, default=80, help='positions per bcftools call')
    ap.add_argument('--fasta', default=FASTA)
    ap.add_argument('--assembly-report', default=REPORT)
    a = ap.parse_args(argv)
    bcftools = find_bcftools()
    if not bcftools:
        sys.exit('bcftools not found (set VINDEX_TOOLS_BIN)')
    with open(os.path.join(a.results, 'run_manifest.json')) as f:
        manifest = json.load(f)
    run = manifest['input']['reads']
    reads, index = a.reads or run['path'], a.reads_index or run.get('index')
    rows = snv_rows(os.path.join(a.results, 'unresolved_loci.parquet'), a.limit_per_state, a.seed)
    summary, lines = compare(rows, reads, index, run['policy'], Reference(a.fasta, a.assembly_report), bcftools, a.processes, a.chunk, lambda m: print(' ', m, file=sys.stderr, flush=True))
    summary.update(limit_per_state=a.limit_per_state, seed=a.seed if a.limit_per_state else None, reads=reads, run_commit=manifest['code']['commit'],
                   reads_module_sha256=sha256_file(os.path.join(HERE, 'vindex', 'reads.py')), tool_sha256=sha256_file(os.path.abspath(__file__)))
    with open(a.out, 'w') as f:
        json.dump(summary, f, indent=1, sort_keys=True)
    with open(os.path.splitext(a.out)[0] + '.loci.tsv', 'w') as f:
        f.write('\t'.join(HEADER) + '\n')
        f.writelines('\t'.join(map(str, line)) + '\n' for line in lines)
    if a.record:
        with open(a.record, 'w') as f:
            json.dump(summary, f, indent=1, sort_keys=True)
            f.write('\n')
    print(json.dumps({k: summary[k] for k in ('loci', 'snv_keys_of_the_run', 'alignments', 'fragments', 'explainable_by', 'by_state', 'n_alignment_problems', 'n_unexplained', 'ok')}, indent=1))
    return 0 if summary['ok'] else 1


if __name__ == '__main__':
    sys.exit(main())
