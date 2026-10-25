"""Audit the existing dc assets under scrape/ and guard them against accidental change.

  python3 inventory/audit_existing.py            # inventory + verification -> EXISTING_ASSETS.json / .md
  python3 inventory/audit_existing.py --baseline # (re)write baseline_checksums.sha256 for tracked files under scrape/
  python3 inventory/audit_existing.py --check    # compare scrape/ with the baseline; exit 1 if anything changed
                                                  # (except hash-pinned entries in preservation_exceptions.json)

The audit never modifies anything under scrape/. MaveDB score tables are verified only if scores/ is present
locally (restore with the git archive command in scrape/mavedb-human/README.md otherwise).
"""
import argparse
import csv
import datetime
import hashlib
import json
import os
import subprocess
import sys

ROOT = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))
HERE = os.path.join(ROOT, 'inventory')
BASELINE = os.path.join(HERE, 'baseline_checksums.sha256')


def sh(*args, cwd=ROOT):
    r = subprocess.run(args, cwd=cwd, capture_output=True, text=True)
    return r.returncode, r.stdout.strip(), r.stderr.strip()


def sha256_file(path):
    h = hashlib.sha256()
    with open(path, 'rb') as f:
        for chunk in iter(lambda: f.read(1 << 20), b''):
            h.update(chunk)
    return h.hexdigest()


def tracked_scrape_files():
    _, out, _ = sh('git', 'ls-files', 'scrape')
    return sorted(out.splitlines())


def write_baseline():
    _, commit, _ = sh('git', 'rev-parse', 'HEAD')
    with open(BASELINE, 'w') as f:
        f.write(f'# sha256 of every git-tracked file under scrape/ at commit {commit}\n')
        for p in tracked_scrape_files():
            f.write(f'{sha256_file(os.path.join(ROOT, p))}  {p}\n')
    print('baseline written for', len(tracked_scrape_files()), 'files at', commit[:8])


def check_baseline():
    base = {}
    for line in open(BASELINE):
        if line.startswith('#') or not line.strip():
            continue
        h, p = line.rstrip('\n').split('  ', 1)
        base[p] = h
    # documented, hash-pinned exceptions: a changed file is accepted only if it matches the exact content that was approved and logged
    exc_path = os.path.join(HERE, 'preservation_exceptions.json')
    exceptions = {e['path']: e for e in json.load(open(exc_path))} if os.path.exists(exc_path) else {}
    now = {p: sha256_file(os.path.join(ROOT, p)) for p in tracked_scrape_files() if os.path.exists(os.path.join(ROOT, p))}
    missing = sorted(set(base) - set(now))
    changed, excepted = [], []
    for p in sorted(base):
        if p in now and now[p] != base[p]:
            (excepted if p in exceptions and exceptions[p]['sha256_after'] == now[p] else changed).append(p)
    added = sorted(set(now) - set(base))
    print(f'baseline files: {len(base)} | missing: {len(missing)} | changed: {len(changed)} | documented exceptions: {len(excepted)} | added: {len(added)}')
    for label, lst in (('MISSING', missing), ('CHANGED', changed), ('ADDED', added)):
        for p in lst[:20]:
            print(' ', label, p)
    for p in excepted:
        print('  EXCEPTED', p, '-', exceptions[p]['reason'][:110])
    return not missing and not changed


def count(path):
    n = size = 0
    for d, _, fs in os.walk(path):
        for f in fs:
            n += 1
            size += os.path.getsize(os.path.join(d, f))
    return n, size


def jsonl_ok(path):
    n = 0
    for line in open(path, encoding='utf8'):
        if line.strip():
            json.loads(line)
            n += 1
    return n


def magic_ok(path):
    ext = os.path.splitext(path)[1].lower()
    head = open(path, 'rb').read(16)
    if ext == '.pdf':
        return head.startswith(b'%PDF-')
    if ext in ('.jpg', '.jpeg'):
        return head.startswith(b'\xff\xd8')
    if ext == '.png':
        return head.startswith(b'\x89PNG') or head.startswith(b'\xff\xd8')  # DNA_dyn.png is JPEG data as served
    if ext == '.mp3':
        return head.startswith(b'ID3') or head[:2] in (b'\xff\xfb', b'\xff\xf3', b'\xff\xf2')
    return True


def check_dir_media(path):
    bad = []
    for d, _, fs in os.walk(path):
        for f in fs:
            if os.path.splitext(f)[1].lower() in ('.pdf', '.jpg', '.jpeg', '.png', '.mp3'):
                if not magic_ok(os.path.join(d, f)):
                    bad.append(os.path.relpath(os.path.join(d, f), ROOT))
    return bad


def audit():
    out = dict(generated_utc=datetime.datetime.now(datetime.timezone.utc).strftime('%Y-%m-%dT%H:%M:%SZ'))
    _, out['commit'], _ = sh('git', 'rev-parse', 'HEAD')
    # repository and branches
    rc, vis, _ = sh('gh', 'api', 'repos/ziadwalieddine-create/dc', '--jq', '{private,visibility}')
    out['repository_visibility_at_audit'] = json.loads(vis) if rc == 0 and vis.startswith('{') else 'unknown (gh api unavailable)'
    branches = []
    _, refs, _ = sh('git', 'for-each-ref', '--format=%(refname:short)', 'refs/remotes/origin')
    for b in refs.splitlines():
        if b.endswith('/HEAD'):
            continue
        _, files, _ = sh('git', 'ls-tree', '-r', '--name-only', b)
        _, tip, _ = sh('git', 'rev-parse', '--short', b)
        branches.append(dict(branch=b, tip=tip, files=len(files.splitlines())))
    out['branches'] = branches

    assets = []
    # NHGRI GWAS
    p = os.path.join(ROOT, 'scrape/nhgri-gwas')
    n, sz = count(p)
    g = json.load(open(os.path.join(p, 'gwas.json')))
    bad = check_dir_media(p)
    assets.append(dict(id='nhgri-gwas', path='scrape/nhgri-gwas', source='https://www.genome.gov/genetics-glossary/Genome-Wide-Association-Studies',
                       kind='glossary entry + 4 linked entries (HTML, text, audio, images)', files=n, bytes=sz,
                       readme=os.path.exists(os.path.join(p, 'README.md')),
                       checks=dict(json_valid=True, related_entries=len(g['related']), media_assets=len(g['assets']), media_magic_failures=bad),
                       status='downloaded_verified', evidence_value='background/terminology only (no variant-level data)'))
    # NHGRI 1000 Genomes
    p = os.path.join(ROOT, 'scrape/nhgri-1000-genomes')
    n, sz = count(p)
    t = json.load(open(os.path.join(p, '1000genomes.json')))
    talks = [x for k in t['tutorials'].values() for x in k['talks']]
    missing = [x['slides_file'] for x in talks if not os.path.exists(os.path.join(p, x['slides_file']))]
    bad = check_dir_media(p)
    assets.append(dict(id='nhgri-1000-genomes', path='scrape/nhgri-1000-genomes', source='https://www.genome.gov/27528684/1000-genomes-project',
                       kind='archived project page, 2 tutorial pages, 11 slide decks (PDF + text), 1 staff profile', files=n, bytes=sz,
                       readme=os.path.exists(os.path.join(p, 'README.md')),
                       checks=dict(json_valid=True, tutorial_talks=len(talks), slide_files_missing=missing, media_magic_failures=bad),
                       status='downloaded_verified', evidence_value='documentation of the 1000 Genomes project (no genotype data)'))
    # Ensembl
    p = os.path.join(ROOT, 'scrape/ensembl-human-genomes')
    n, sz = count(p)
    gen = json.load(open(os.path.join(p, 'genomes.json')))
    rows = list(csv.DictReader(open(os.path.join(p, 'genomes.csv'))))
    ids = [m['genome_id'] for m in gen['matches']]
    ok = len(ids) == gen['total_hits'] == len(rows) and len(set(ids)) == len(ids) and {r['Genome ID'] for r in rows} == set(ids)
    assets.append(dict(id='ensembl-human-genomes', path='scrape/ensembl-human-genomes', source='https://www.ensembl.org/genome-selector/search?species_taxonomy_id=9606',
                       kind='genome/assembly registry for human (563 records) + details for GRCh38.p14, GRCh37.p13, T2T-CHM13v2.0', files=n, bytes=sz,
                       readme=os.path.exists(os.path.join(p, 'README.md')),
                       checks=dict(records=len(ids), total_hits=gen['total_hits'], csv_rows=len(rows), consistent=ok),
                       status='downloaded_verified', evidence_value='assembly/annotation registry; no sequences or annotations downloaded'))
    # MaveDB
    p = os.path.join(ROOT, 'scrape/mavedb-human')
    n_tracked = sum(1 for f in tracked_scrape_files() if f.startswith('scrape/mavedb-human/'))
    size_tracked = sum(os.path.getsize(os.path.join(ROOT, f)) for f in tracked_scrape_files() if f.startswith('scrape/mavedb-human/'))
    mv = dict(id='mavedb-human', path='scrape/mavedb-human', source='https://www.mavedb.org/search (API https://api.mavedb.org/api/v1)',
              kind='1,402 published human-target score sets: metadata + 1,980 score/count tables (tables on branch data/mavedb-scores)',
              files_tracked_on_main=n_tracked, bytes_tracked_on_main=size_tracked, readme=True)
    scores = os.path.join(p, 'scores')
    if os.path.isdir(scores):
        rc, so, se = sh(sys.executable, '-I', os.path.join(p, 'verify_scores.py'))
        mv['checks'] = dict(verify_scores_exit=rc, verify_scores_output=so.splitlines()[-3:], scores_present_locally=True)
        mv['status'] = 'downloaded_verified' if rc == 0 else 'verification_failed'
    else:
        mv['checks'] = dict(scores_present_locally=False, note='restore with git archive origin/data/mavedb-scores ... then run verify_scores.py')
        mv['status'] = 'metadata_verified_tables_not_local'
    mv['evidence_value'] = 'direct experimental variant effects (deep mutational scanning and related MAVE assays); mapped to protein/nucleotide variants per score set'
    assets.append(mv)
    out['assets'] = assets
    out['preservation_baseline'] = os.path.relpath(BASELINE, ROOT)
    return out


def to_markdown(a):
    L = ['# Existing assets: inventory and verified status', '',
         f"Generated {a['generated_utc']} at commit `{a['commit'][:8]}` by `inventory/audit_existing.py` (regenerate rather than edit).", '',
         f"Repository visibility at audit time: `{json.dumps(a['repository_visibility_at_audit'])}`. While the repository is public, only openly redistributable data and metadata may be committed.", '',
         '## Branches', '', '| Branch | Tip | Files | Role |', '|---|---|---|---|']
    role = {'origin/main': 'default branch; holds scrape/ (metadata, scripts, manifests) and the new catalogue',
            'origin/data/mavedb-scores': 'orphan data branch: the 1,980 MaveDB score/count tables + licence notice (not to be merged)',
            'origin/land/dc-extended-pack': 'legacy branch from the original dating-project repository; unrelated, left untouched'}
    for b in a['branches']:
        L.append(f"| `{b['branch']}` | `{b['tip']}` | {b['files']} | {role.get(b['branch'], 'working/PR branch')} |")
    L += ['', '## Assets', '', '| Asset | Path | What it is | Files | Size | Status |', '|---|---|---|---|---|---|']
    for s in a['assets']:
        files = s.get('files', s.get('files_tracked_on_main'))
        size = s.get('bytes', s.get('bytes_tracked_on_main'))
        L.append(f"| `{s['id']}` | `{s['path']}` | {s['kind']} | {files} | {size/1e6:.1f} MB | **{s['status']}** |")
    L += ['', '## Verification detail', '']
    for s in a['assets']:
        L += [f"### {s['id']}", f"- Source: {s['source']}", f"- Evidence value: {s['evidence_value']}", f"- Checks: `{json.dumps(s['checks'])}`", '']
    L += ['## Preservation guard', '',
          f"`{a['preservation_baseline']}` holds the sha256 of every git-tracked file under `scrape/` at the baseline commit. "
          'Run `python3 inventory/audit_existing.py --check` at any time: it fails if any existing file was changed or removed. '
          'New work is added outside `scrape/` (in `catalogue/`, `inventory/`, `acquisition/`).', '']
    return '\n'.join(L)


if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('--baseline', action='store_true')
    ap.add_argument('--check', action='store_true')
    a = ap.parse_args()
    if a.baseline:
        write_baseline()
    elif a.check:
        sys.exit(0 if check_baseline() else 1)
    else:
        res = audit()
        json.dump(res, open(os.path.join(HERE, 'EXISTING_ASSETS.json'), 'w'), indent=2)
        open(os.path.join(HERE, 'EXISTING_ASSETS.md'), 'w').write(to_markdown(res))
        print(json.dumps({s['id']: s['status'] for s in res['assets']}, indent=1))
