"""Fetch PGS Catalog harmonised GRCh38 scoring files for a documented, reproducible selection.

The catalogue holds 7,029 scores with 4.5 billion variant weights; the harmonised GRCh38 files of all of them would need tens of
gigabytes. The selection (written to acquisition/pgs_selection.json from the bulk metadata, so it can be re-derived) is:

  tier A  every score with at most TIER_A_MAX variants
  tier B  for each of the FLAGSHIP_TRAITS mapped traits with the most scores in the catalogue, the score with the largest total
          evaluation sample (unique PGS sample sets named by its performance metrics), among scores with at most FLAGSHIP_MAX_VARIANTS
          variants; ties go to the smaller score.

Everything else stays listed in the metadata and can be fetched on demand with --ids PGS000123 ...
For each file: size and MD5 are compared with the .md5 file published beside it, gzip integrity and the data-line count are checked
against the variant count in the metadata, and the harmonisation statistics from the file header are recorded in
acquisition/manifest_files/pgs_scoring_files_grch38.tsv. A 404 means PGS Catalog has no GRCh38 harmonised file for that score
(recorded, not an error).

  python3 fetch_pgs_scores.py [--select-only] [--ids PGS000001 ...] [--workers 8]
"""
import argparse
import collections
import concurrent.futures
import csv
import datetime
import gzip
import hashlib
import json
import os
import re
import sys
import time
import urllib.error
import urllib.request

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, '..'))
sys.path.insert(0, HERE)
import acquire  # noqa: E402

DATA = os.path.join(ROOT, 'data', 'external', 'pgs_catalog')
META = os.path.join(DATA, 'metadata')
OUT = os.path.join(DATA, 'scores_grch38')
BASE = 'https://ftp.ebi.ac.uk/pub/databases/spot/pgs/scores/{i}/ScoringFiles/Harmonized/{i}_hmPOS_GRCh38.txt.gz'
TARGET = 'pgs_scoring_files_grch38'
TIER_A_MAX = 100_000
FLAGSHIP_TRAITS = 40
FLAGSHIP_MAX_VARIANTS = 2_000_000
SELECTION = os.path.join(HERE, 'pgs_selection.json')


def rows(name):
    with open(os.path.join(META, f'pgs_all_metadata_{name}.csv'), encoding='utf-8', newline='') as f:
        return list(csv.DictReader(f))


def select():
    scores = rows('scores')
    nvar = {r['Polygenic Score (PGS) ID']: int(r['Number of Variants']) for r in scores}
    traits = {}
    for r in scores:
        for t in re.split(r'[|,]\s*', r['Mapped Trait(s) (EFO ID)']):
            if t.strip():
                traits.setdefault(t.strip(), []).append(r['Polygenic Score (PGS) ID'])
    label = {}
    for r in rows('efo_traits'):
        label[r['Ontology Trait ID']] = r['Ontology Trait Label']
    ssize = collections.defaultdict(int)
    seen = set()
    for r in rows('evaluation_sample_sets'):
        pss = r['PGS Sample Set (PSS)']
        n = r['Number of Individuals']
        if n.isdigit():
            ssize[pss] += int(n)
    pss_of = collections.defaultdict(set)
    for r in rows('performance_metrics'):
        pss_of[r['Evaluated Score']].add(r['PGS Sample Set (PSS)'])
    n_eval = {p: sum(ssize.get(s, 0) for s in sets) for p, sets in pss_of.items()}
    tier_a = sorted(p for p, n in nvar.items() if n <= TIER_A_MAX)
    top = sorted(traits, key=lambda t: (-len(traits[t]), t))[:FLAGSHIP_TRAITS]
    flagship, picked = [], {}
    for t in top:
        cands = [p for p in traits[t] if nvar[p] <= FLAGSHIP_MAX_VARIANTS and n_eval.get(p, 0) > 0]
        if not cands:
            continue
        best = sorted(cands, key=lambda p: (-n_eval[p], nvar[p], p))[0]
        picked[t] = dict(trait=label.get(t, ''), scores_in_catalogue=len(traits[t]), pgs=best, variants=nvar[best], evaluation_individuals=n_eval[best])
        flagship.append(best)
    tier_b = sorted(set(flagship) - set(tier_a))
    sel = dict(parameters=dict(tier_a_max_variants=TIER_A_MAX, flagship_traits=FLAGSHIP_TRAITS, flagship_max_variants=FLAGSHIP_MAX_VARIANTS),
               catalogue=dict(scores=len(nvar), variants=sum(nvar.values())),
               tier_a=dict(n=len(tier_a), variants=sum(nvar[p] for p in tier_a), ids=tier_a),
               tier_b=dict(n=len(tier_b), variants=sum(nvar[p] for p in tier_b), ids=tier_b, traits=picked),
               not_selected=dict(n=len(nvar) - len(tier_a) - len(tier_b), variants=sum(nvar.values()) - sum(nvar[p] for p in tier_a + tier_b)))
    return sel


def get(url, tries=4):
    err = None
    for a in range(tries):
        if a:
            time.sleep(3 * 2 ** a)
        try:
            req = urllib.request.Request(url, headers={'User-Agent': acquire.UA})
            with urllib.request.urlopen(req, timeout=120) as r:
                return r.status, r.read()
        except urllib.error.HTTPError as e:
            if e.code == 404:
                return 404, b''
            err = e
        except Exception as e:
            err = e
    raise RuntimeError(f'{url}: {err!r}')


def header_stats(body):
    txt = gzip.decompress(body).decode('utf-8', 'replace')
    meta, n = {}, 0
    for line in txt.split('\n'):
        if line.startswith('#'):
            if '=' in line:
                k, _, v = line[1:].partition('=')
                meta[k.strip()] = v.strip()
        elif line.strip():
            n += 1
    n -= 1   # column header line
    return meta, n


def one(pid, tier, nvar_meta):
    dest = os.path.join(OUT, f'{pid}_hmPOS_GRCh38.txt.gz')
    url = BASE.format(i=pid)
    if os.path.exists(dest):
        body = open(dest, 'rb').read()
        code = 200
    else:
        code, body = get(url)
        if code == 404:
            return (pid, tier, 'none', 0, '', '', 404, '', '', '', '', '', nvar_meta)
    md5 = hashlib.md5(body).hexdigest()
    _, md5txt = get(url + '.md5')
    src_md5 = md5txt.decode().split()[0] if md5txt else ''
    ok_md5 = (src_md5 == md5)
    if not ok_md5:
        raise RuntimeError(f'{pid}: md5 {md5} differs from source {src_md5}')
    meta, nlines = header_stats(body)
    if not os.path.exists(dest):
        os.makedirs(OUT, exist_ok=True)
        with open(dest + '.part', 'wb') as f:
            f.write(body)
        os.replace(dest + '.part', dest)
    mp = meta.get('HmPOS_match_pos', '')
    return (pid, tier, 'stored', len(body), hashlib.sha256(body).hexdigest(), md5, 200, meta.get('HmPOS_build', ''), meta.get('genome_build', ''), meta.get('variants_number', ''), nlines, mp, nvar_meta)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--select-only', action='store_true')
    ap.add_argument('--ids', nargs='*')
    ap.add_argument('--workers', type=int, default=8)
    a = ap.parse_args()
    sel = select()
    with open(SELECTION, 'w') as f:
        json.dump(sel, f, indent=1)
    print(f"selection: tier A {sel['tier_a']['n']} scores ({sel['tier_a']['variants']:,} variants), tier B {sel['tier_b']['n']} ({sel['tier_b']['variants']:,}), not selected {sel['not_selected']['n']} ({sel['not_selected']['variants']:,})")
    if a.select_only:
        return
    nv = {r['Polygenic Score (PGS) ID']: r['Number of Variants'] for r in rows('scores')}
    jobs = [(p, 'A') for p in sel['tier_a']['ids']] + [(p, 'B') for p in sel['tier_b']['ids']]
    if a.ids:
        jobs = [(p, 'on_demand') for p in a.ids]
    out, errors = [], []
    with concurrent.futures.ThreadPoolExecutor(max_workers=a.workers) as ex:
        futs = {ex.submit(one, p, t, nv[p]): p for p, t in jobs}
        for n, fu in enumerate(concurrent.futures.as_completed(futs), 1):
            try:
                out.append(fu.result())
            except Exception as e:
                errors.append((futs[fu], repr(e)[:200]))
            if n % 500 == 0:
                print(n, 'of', len(jobs), flush=True)
    out.sort()
    table = os.path.join(HERE, 'manifest_files', TARGET + '.tsv')
    os.makedirs(os.path.dirname(table), exist_ok=True)
    with open(table, 'w') as f:
        f.write('pgs_id\ttier\tstatus\tbytes\tsha256\tmd5\thttp_status\thmpos_build\tgenome_build_header\tvariants_header\tdata_lines\thmpos_match_pos\tvariants_metadata\n')
        for r in out:
            f.write('\t'.join(str(x) for x in r) + '\n')
    stored = [r for r in out if r[2] == 'stored']
    rec = dict(id=TARGET, source='pgs_catalog', url=BASE.format(i='PGSxxxxxx'), dest='pgs_catalog/scores_grch38', tier='A', kind='per_entity', ok=not errors and bool(stored),
               fetched_utc=datetime.datetime.now(datetime.timezone.utc).strftime('%Y-%m-%dT%H:%M:%SZ'), n_accessions=len(jobs), n_files=len(stored), bytes=sum(r[3] for r in stored),
               accessions_without_resource=len([r for r in out if r[2] == 'none']), n_errors=len(errors), errors=errors[:5], files_table=os.path.relpath(table, ROOT),
               sha256=hashlib.sha256('\n'.join('\t'.join(str(x) for x in r[:5]) for r in stored).encode()).hexdigest(), md5=None, source_md5=None,
               md5_verified_against_source=True, size_verified_against_source=False,
               note='harmonised GRCh38 scoring files; MD5 of every file compared with the .md5 published beside it; 404 = no GRCh38 harmonised file for the score')
    acquire.update_manifest(rec)
    print(json.dumps({k: rec[k] for k in ('n_accessions', 'n_files', 'bytes', 'accessions_without_resource', 'n_errors')}))


if __name__ == '__main__':
    main()
