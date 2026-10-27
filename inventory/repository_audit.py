"""Repository audit: what exists, what state it is in, and whether the documentation matches the implementation.

  python3 inventory/repository_audit.py             # writes REPOSITORY_AUDIT.json and REPOSITORY_AUDIT.md
  python3 inventory/repository_audit.py --check     # same, and exits 1 if a test or a documentation claim fails
  python3 inventory/repository_audit.py --fast      # reads inventory/EXISTING_ASSETS.json instead of re-running the MaveDB score verification

Dataset lifecycle (one row per acquisition target), each stage computed from files and never asserted:
  listed      the target is in acquisition/targets.json
  acquired    the manifest says ok AND the file or directory is present in data/external/ with the recorded size
  verified    the validator's checks passed (acquisition/validation_report.json)
  mapped      at least one passing compatibility test exercised it (acquisition/compat_report.json), or the built index holds an evidence row from it that is mapped to a normalised
              entity (the lifecycle report of the index, data/index/reports/lifecycle.json: a target that was integrated is mapped even if no join test names it)
  integrated  the built index (data/index/index_manifest.json) holds a release of the target, so its records are queryable as evidence or entities with provenance;
              without a built index in the checkout no target is reported integrated

Documentation claims are numbers and statements in the READMEs, the latest AUDIT_LOG entry and the curated gap texts. Each claim names the
document, a regular expression that must match the text, and a function that recomputes the value from the files. A claim whose text is no
longer in the document, or whose number differs from the recomputed value, fails. A claim that needs data absent from this checkout is
reported as skipped, not as passed.
"""
import argparse
import ast
import collections
import csv
import datetime
import functools
import glob
import json
import os
import re
import subprocess
import sys

ROOT = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))
HERE = os.path.join(ROOT, 'inventory')
DATA = os.path.join(ROOT, 'data', 'external')


def read_text(p, **kw):
    with open(p, **kw) as f:
        return f.read()


def read_lines(p):
    with open(p) as f:
        return f.readlines()


def read_csv(p):
    with open(p, newline='') as f:
        return list(csv.DictReader(f))


def sh(*args, cwd=ROOT):
    r = subprocess.run(args, cwd=cwd, capture_output=True, text=True)
    return r.returncode, r.stdout.strip(), r.stderr.strip()


def jload(rel, default=None):
    p = os.path.join(ROOT, rel)
    return json.loads(read_text(p)) if os.path.exists(p) else default


def jsonl(rel):
    p = os.path.join(ROOT, rel)
    return [json.loads(x) for x in read_lines(p) if x.strip()] if os.path.exists(p) else []


# ---- facts recomputed from the files ----------------------------------------------------------------------------------------------
@functools.lru_cache(None)
def mave_index():
    return read_csv(os.path.join(ROOT, 'scrape/mavedb-human/score_sets_index.csv'))


@functools.lru_cache(None)
def mave_details():
    out = []
    for fn in ('score_set_details.jsonl', 'score_set_details_accession_targets.jsonl'):
        out += jsonl('scrape/mavedb-human/' + fn)
    return out


@functools.lru_cache(None)
def mave_manifest():
    return jsonl('scrape/mavedb-human/scores_manifest.jsonl')


@functools.lru_cache(None)
def genomes():
    return jload('scrape/ensembl-human-genomes/genomes.json')['matches']


@functools.lru_cache(None)
def validation():
    return jload('acquisition/validation_report.json', {})


@functools.lru_cache(None)
def compat():
    return jload('acquisition/compat_report.json', {})


@functools.lru_cache(None)
def manifest():
    return {r['id']: r for r in jsonl('acquisition/manifest.jsonl')}


@functools.lru_cache(None)
def targets():
    return {t['id']: t for t in jload('acquisition/targets.json', [])}


def compat_counts():
    c = collections.Counter()
    for d in compat().get('sources', {}).values():
        for t in d['tests']:
            c['PASS' if t['passed'] is True else ('FAIL' if t['passed'] is False else 'INFO')] += 1
    return c


def vstat(tid, key):
    return validation().get(tid, {}).get('stats', {}).get(key)


def need_data(*rel):
    """Raise Skip when a needed file is absent from this checkout."""
    for r in rel:
        if not os.path.exists(os.path.join(ROOT, r)):
            raise Skip(r)


class Skip(Exception):
    pass


def count_files(*rels):
    n = 0
    for r in rels:
        need_data(r)
        for _, _, fs in os.walk(os.path.join(ROOT, r)):
            n += len(fs)
    return n


UNIT = {}


def unit_tests():
    if not UNIT:
        rc, out, err = sh(sys.executable, '-m', 'unittest', 'acquisition.test_helpers')
        m = re.search(r'Ran (\d+) tests?', err + out)
        UNIT.update(rc=rc, n=int(m.group(1)) if m else 0, ok=rc == 0 and 'OK' in (err + out))
    return UNIT


VIDX = {}


def vindex_tests():
    """Run the variant_index test suite when this interpreter has its dependencies (pysam, pyarrow); otherwise raise Skip."""
    if not VIDX:
        try:
            import pyarrow  # noqa: F401
            import pysam  # noqa: F401
        except ImportError as e:
            raise Skip(f'variant_index tests need {e.name}')
        rc, out, err = sh(sys.executable, '-W', 'error::ResourceWarning', '-m', 'unittest', 'discover', cwd=os.path.join(ROOT, 'variant_index', 'tests'))
        m = re.search(r'Ran (\d+) tests?', err + out)
        VIDX.update(rc=rc, n=int(m.group(1)) if m else 0, ok=rc == 0 and 'OK' in (err + out))
    return VIDX


def latest_audit_entry():
    t = read_text(os.path.join(ROOT, 'AUDIT_LOG.md'))
    parts = re.split(r'(?m)^## ', t)
    return parts[-1]


def mave_counts(field):
    return collections.Counter(r[field] for r in mave_index())


# (document, regular expression with one group per number, function returning the recomputed numbers in group order)
# A document of the form 'latest:AUDIT_LOG.md' means its last second-level section.
CLAIMS = [
    ('README.md', r'(\d+) genome/assembly records', lambda: [len(genomes())]),
    ('README.md', r'MaveDB: ([\d,]+) published score sets with a human target \(([\d,]+) by organism filter plus (\d+)',
     lambda: [len(mave_index()), mave_counts('selection')['organism_filter'], mave_counts('selection')['accession_target']]),
    ('README.md', r'the ([\d,]+) score and count tables', lambda: [len(mave_manifest())]),
    ('README.md', r'an additive set of (\d+) superseded MaveDB score sets', lambda: [len(read_csv(os.path.join(ROOT, 'acquisition/mavedb_superseded/superseded_index.csv')))]),
    ('README.md', r'(\d+) tutorial slide decks', lambda: [len([x for k in jload('scrape/nhgri-1000-genomes/1000genomes.json')['tutorials'].values() for x in k['talks']])]),
    ('scrape/mavedb-human/README.md', r'([\d,]+) in total', lambda: [len(mave_index())]),
    ('scrape/mavedb-human/README.md', r'\| `organism_filter` \|.*?\| ([\d,]+) \|', lambda: [mave_counts('selection')['organism_filter']]),
    ('scrape/mavedb-human/README.md', r'\| `accession_target` \|.*?\| (\d+) \|', lambda: [mave_counts('selection')['accession_target']]),
    ('scrape/mavedb-human/README.md', r'Together: ([\d,]+) experiments, ([\d,]+) distinct target names, ([\d,]+) variants',
     lambda: [len({r['experiment_urn'] for r in mave_index()}), len({t.get('name') for d in mave_details() for t in d.get('targetGenes', [])}), sum(int(r['num_variants']) for r in mave_index())]),
    ('scrape/mavedb-human/README.md', r'([\d,]+) files, about [\d]+ MB gzipped, ([\d,]+) score rows and ([\d,]+) count rows',
     lambda: [len(mave_manifest()), sum(m['rows'] for m in mave_manifest() if m['kind'] == 'scores'), sum(m['rows'] for m in mave_manifest() if m['kind'] == 'counts')]),
    ('scrape/mavedb-human/README.md', r'for the (\d+) sets with count columns', lambda: [sum(1 for r in mave_index() if r['count_columns'])]),
    ('scrape/mavedb-human/README.md', r'(\d+) sets have `processingState: incomplete`', lambda: [mave_counts('processing_state')['incomplete']]),
    ('scrape/mavedb-human/README.md', r'`complete` for only (\d+) of ([\d,]+) \((\d+) incomplete, (\d+) failed, (\d+) queued, (\d+) processing, (\d+) none\)',
     lambda: [mave_counts('mapping_state')[k] for k in ('complete',)] + [len(mave_index())] + [mave_counts('mapping_state')[k] for k in ('incomplete', 'failed', 'queued', 'processing', '')]),
    ('scrape/mavedb-human/README.md', r'([\d,]+) CC0, (\d+) CC BY 4\.0, (\d+) CC BY-SA 4\.0, (\d+) CC BY-NC-SA 4\.0', lambda: [mave_counts('license')[k] for k in ('CC0', 'CC BY 4.0', 'CC BY-SA 4.0', 'CC BY-NC-SA 4.0')]),
    ('scrape/mavedb-human/README.md', r'only `urn:mavedb:00000053-a-1`, ([\d,]+) variants', lambda: [max(int(r['num_variants']) for r in mave_index())]),
    ('scrape/ensembl-human-genomes/README.md', r'^(\d+) records, all in integrated release 2026-07', lambda: [len(genomes()) if {g['release']['name'] for g in genomes()} == {'2026-07'} else -1]),
    ('scrape/ensembl-human-genomes/README.md', r'\((\d+) have `hprc` in the assembly name', lambda: [sum(1 for g in genomes() if 'hprc' in json.dumps(g['assembly']).lower())]),
    ('scrape/ensembl-human-genomes/README.md', r'(\d+) records have variation data and (\d+) have regulation', lambda: [sum(1 for g in genomes() if g['has_variation']), sum(1 for g in genomes() if g['has_regulation'])]),
    ('acquisition/README.md', r'\((\d[\d,]*) files, no differences\)', lambda: [(need_data('data/external/mavedb/dump_v4_extract'), count_files('data/external/mavedb/dump_v4_extract', 'data/external/mavedb/dump_v4_sets'))[1]]),
    ('acquisition/README.md', r'\((\d[\d,]*) of ([\d,]+), none elsewhere\)', lambda: [vstat('alphamissense_hg38', 'class_boundary_rounding_rows'), vstat('alphamissense_hg38', 'rows')]),
    ('acquisition/README.md', r'# (\d+) tests, no data needed', lambda: [unit_tests()['n']]),
    ('catalogue/MASTER_INVENTORY.md', r'(\d+) sources catalogued', lambda: [len(glob.glob(os.path.join(ROOT, 'catalogue/sources/*.json')))]),
    ('latest:AUDIT_LOG.md', r'`validate\.py`: (\d+) of (\d+) targets pass', lambda: [sum(1 for v in validation().values() if v['passed']), len(validation())]),
    ('latest:AUDIT_LOG.md', r'(\d+) manifest records', lambda: [len(manifest())]),
    ('latest:AUDIT_LOG.md', r'`compat_tests\.py`: (\d+) tests, (\d+) PASS, (\d+) FAIL, (\d+) INFO', lambda: [sum(compat_counts().values()), compat_counts()['PASS'], compat_counts()['FAIL'], compat_counts()['INFO']]),
    ('latest:AUDIT_LOG.md', r'test_helpers`: (\d+) tests pass', lambda: [unit_tests()['n'] if unit_tests()['ok'] else -1]),
    ('catalogue/gaps.json', r'GTF: ([\d,]+) genes, ([\d,]+) transcripts, ([\d,]+) MANE Select', lambda: [vstat('ensembl_grch38_gtf_116', 'feature_counts')['gene'], vstat('ensembl_grch38_gtf_116', 'feature_counts')['transcript'], vstat('ensembl_grch38_gtf_116', 'tags')['MANE_Select']]),
    ('catalogue/gaps.json', r'toplevel FASTA: (\d+) sequences', lambda: [vstat('ensembl_grch38_dna_toplevel', 'records')]),
    ('catalogue/gaps.json', r'\(([\d.]+)M PDB chain segments, ([\d,]+) PDB entries, ([\d,]+) UniProt accessions', lambda: [round(vstat('sifts_pdb_chain_uniprot', 'rows') / 1e6, 2), vstat('sifts_pdb_chain_uniprot', 'pdb_entries'), vstat('sifts_pdb_chain_uniprot', 'uniprot_accessions')]),
    ('catalogue/gaps.json', r'\(([\d,]+) VCF files, ([\d,]+) variant records\)', lambda: [vstat('pharmvar_6_2_29', 'vcf_files'), vstat('pharmvar_6_2_29', 'vcf_records')]),
    ('catalogue/gaps.json', r'with ([\d,]+) haplotype names', lambda: [sum(vstat('pharmvar_6_2_29', 'haplotypes_refseqgene').values())]),
]


def to_num(x):
    return float(str(x).replace(',', ''))


def claim_text(doc):
    if doc.startswith('latest:'):
        return latest_audit_entry()
    p = os.path.join(ROOT, doc)
    return read_text(p) if os.path.exists(p) else None


def selftest_claims():
    """Adversarial test of the checker itself: change the first number of every claim in the text and require the claim to fail."""
    survivors = []
    for doc, rx, fn in CLAIMS:
        text = claim_text(doc)
        m = re.search(rx, text, re.S | re.M) if text is not None else None
        if not m:
            continue
        a, b = m.span(1)
        mutated = text[:a] + str(int(to_num(m.group(1))) + 7) + text[b:]
        m2 = re.search(rx, mutated, re.S | re.M)
        try:
            want = [round(float(v), 2) for v in fn()]
        except Skip:
            continue
        got = [round(to_num(g), 2) for g in m2.groups()] if m2 else None
        if got == want:
            survivors.append((doc, rx))
    return survivors


def check_claims():
    out = []
    for doc, rx, fn in CLAIMS:
        text = claim_text(doc)
        m = re.search(rx, text, re.S | re.M) if text is not None else None
        rec = dict(doc=doc, pattern=rx)
        if text is None:
            out.append(dict(rec, status='FAIL', detail='the document does not exist'))
            continue
        if not m:
            out.append(dict(rec, status='FAIL', detail='the text is no longer in the document'))
            continue
        try:
            want = [round(float(v), 2) for v in fn()]
        except Skip as e:
            out.append(dict(rec, status='SKIP', detail=f'needs {e} (not in this checkout)'))
            continue
        except Exception as e:  # noqa: BLE001
            out.append(dict(rec, status='FAIL', detail=f'recomputation failed: {e!r}'))
            continue
        got = [round(to_num(g), 2) for g in m.groups()]
        out.append(dict(rec, status='PASS' if got == want else 'FAIL', detail=f'document says {got}, recomputed {want}' if got != want else f'{got}'))
    return out


# ---- inventory sections -----------------------------------------------------------------------------------------------------------
def branches():
    rc, out, _ = sh('git', 'for-each-ref', '--format=%(refname:short)|%(objectname:short)|%(committerdate:short)|%(subject)', 'refs/heads', 'refs/remotes')
    res = []
    for line in out.splitlines():
        name, tip, date, subj = line.split('|', 3)
        if name.endswith('/HEAD') or name == 'origin':
            continue
        _, files, _ = sh('git', 'ls-tree', '-r', '--name-only', name)
        _, lr, _ = sh('git', 'rev-list', '--left-right', '--count', f'origin/main...{name}')
        behind, ahead = (lr.split() + ['?', '?'])[:2]
        res.append(dict(branch=name, tip=tip, date=date, subject=subj[:90], tracked_files=len(files.splitlines()), ahead_of_main=ahead, behind_main=behind))
    return res


def documents():
    _, out, _ = sh('git', 'ls-files', '*.md')
    res = []
    for p in out.splitlines():
        if p.startswith('catalogue/evidence/'):
            continue
        t = read_text(os.path.join(ROOT, p), encoding='utf-8', errors='replace')
        title = next((x[2:].strip() for x in t.splitlines() if x.startswith('# ')), '')
        _, d, _ = sh('git', 'log', '-1', '--format=%cs', '--', p)
        res.append(dict(path=p, title=title[:90], lines=t.count('\n'), bytes=len(t.encode()), last_commit=d))
    return res


def scripts():
    _, out, _ = sh('git', 'ls-files', '*.py', '*.js')
    res = []
    for p in out.splitlines():
        full = os.path.join(ROOT, p)
        src = read_text(full, encoding='utf-8', errors='replace')
        rec = dict(path=p, lines=src.count('\n'))
        if p.endswith('.py'):
            try:
                compile(src, full, 'exec')                   # syntax check without writing bytecode
                rec['compiles'] = True
                doc = ast.get_docstring(ast.parse(src)) or ''
                rec['purpose'] = doc.strip().splitlines()[0][:110] if doc else ''
            except Exception as e:  # noqa: BLE001
                rec['compiles'] = False
                rec['purpose'] = repr(e)[:110]
        else:
            rec['compiles'] = None
            rec['purpose'] = (next((x for x in src.splitlines() if x.strip().startswith(('//', '/*', '*'))), '') or '').strip()[:110]
        res.append(rec)
    return res


def manifests():
    res = []
    m = jsonl('acquisition/manifest.jsonl')
    res.append(dict(path='acquisition/manifest.jsonl', records=len(m), detail=f"{sum(1 for r in m if r.get('ok'))} ok"))
    for p in sorted(glob.glob(os.path.join(ROOT, 'acquisition/manifest_files/*.tsv'))):
        res.append(dict(path=os.path.relpath(p, ROOT), records=len(read_lines(p)) - 1, detail='per-file table'))
    sm = jsonl('scrape/mavedb-human/scores_manifest.jsonl')
    res.append(dict(path='scrape/mavedb-human/scores_manifest.jsonl', records=len(sm), detail='sha256 and row count of each score and count table'))
    res.append(dict(path='inventory/baseline_checksums.sha256', records=sum(1 for x in read_lines(os.path.join(HERE, 'baseline_checksums.sha256')) if not x.startswith('#') and x.strip()), detail='sha256 of every tracked file under scrape/ at the baseline'))
    res.append(dict(path='catalogue/evidence/evidence_index.jsonl', records=len(jsonl('catalogue/evidence/evidence_index.jsonl')), detail='append-only record of every probe (URL, time, status, sha256)'))
    return res


def index_targets():
    """Target ids that have a release in the built index, or None when no index is present in this checkout."""
    p = os.path.join(ROOT, 'data', 'index', 'index_manifest.json')
    if not os.path.exists(p):
        return None
    return {r.split('/', 1)[1].rsplit('@', 1)[0] for r in json.loads(read_text(p)).get('releases', []) if '/' in r}


def index_lifecycle():
    """target id -> row of the index's own lifecycle report (make_reports.py), or {} when this checkout holds no report."""
    p = os.path.join(ROOT, 'data', 'index', 'reports', 'lifecycle.json')
    if not os.path.exists(p):
        return {}
    return {r['target_id']: r for r in json.loads(read_text(p)).get('targets', [])}


def dataset_lifecycle():
    in_index = index_targets()
    by_index = index_lifecycle()
    comp_targets = collections.defaultdict(list)
    for sid, d in compat().get('sources', {}).items():
        for t in d['tests']:
            if t['passed'] is True:
                for x in t.get('targets', []):
                    comp_targets[x].append(t['name'][:60])
    res = []
    for tid, t in targets().items():
        m = manifest().get(tid, {})
        p = os.path.join(DATA, t['dest'])
        acquired, why = False, 'no manifest record'
        if m.get('ok'):
            if t.get('kind') == 'http_dir':
                size = sum(os.path.getsize(os.path.join(p, f)) for f in os.listdir(p) if f.endswith('.parquet')) if os.path.isdir(p) else None
                acquired, why = size == m.get('bytes'), f'{size} bytes on disk, manifest {m.get("bytes")}'
            elif t.get('kind') == 'stream_derive':          # the archive was streamed and parsed, not stored: the derived Parquet tables must be there with the recorded sizes
                tabs = m.get('tables') or {}
                bad = []
                for name, tm in tabs.items():
                    d = os.path.join(p, name)
                    files = [f for f in os.listdir(d) if f.endswith('.parquet')] if os.path.isdir(d) else []
                    size = sum(os.path.getsize(os.path.join(d, f)) for f in files)
                    if len(files) != tm.get('files') or size != tm.get('bytes'):
                        bad.append(f'{name}: {len(files)} files / {size} bytes, manifest {tm.get("files")} / {tm.get("bytes")}')
                acquired, why = bool(tabs) and not bad, '; '.join(bad) or f'{len(tabs)} derived tables present with the recorded sizes'
            elif t.get('kind') == 'per_entity':
                acquired, why = os.path.isdir(p) and any(True for _ in os.scandir(p)), 'directory present' if os.path.isdir(p) else 'directory missing'
            else:
                size = os.path.getsize(p) if os.path.isfile(p) else None
                acquired, why = size == m.get('bytes'), f'{size} bytes on disk, manifest {m.get("bytes")}'
        elif m:
            why = 'manifest record not ok'
        v = validation().get(tid)
        res.append(dict(id=tid, source=t['source'], tier=t.get('tier'), kind=t.get('kind'), supports=t.get('supports', []), listed=True, acquired=bool(acquired), acquired_detail=why if not acquired else '',
                        verified=bool(v and v['passed']), mapped=bool(comp_targets.get(tid)) or bool(by_index.get(tid, {}).get('mapped')), integrated=bool(in_index and tid in in_index), bytes=m.get('bytes')))
    return res


def run_checks(fast):
    res = []
    rc, out, err = sh(sys.executable, os.path.join(ROOT, 'inventory/audit_existing.py')) if not fast else (0, '', '')
    existing = jload('inventory/EXISTING_ASSETS.json', {})
    res.append(dict(check='existing assets re-audited (NHGRI, Ensembl registry, MaveDB score tables)', ok=rc == 0 and all(a['status'] == 'downloaded_verified' for a in existing.get('assets', [])),
                    detail=json.dumps({a['id']: a['status'] for a in existing.get('assets', [])})))
    rc, out, err = sh(sys.executable, os.path.join(ROOT, 'inventory/audit_existing.py'), '--check')
    res.append(dict(check='preservation guard: tracked files under scrape/ unchanged apart from documented exceptions', ok=rc == 0, detail=out.splitlines()[0] if out else err[:200]))
    u = unit_tests()
    res.append(dict(check='unit tests (acquisition.test_helpers)', ok=u['ok'], detail=f"{u['n']} tests"))
    rc, out, err = sh(sys.executable, os.path.join(ROOT, 'catalogue/build_catalogue.py'), '--check')
    try:
        v = vindex_tests()
        res.append(dict(check='unit tests (variant_index)', ok=v['ok'], detail=f"{v['n']} tests"))
    except Skip as e:
        res.append(dict(check='unit tests (variant_index)', ok=True, detail=f'skipped: {e}'))
    res.append(dict(check='catalogue records validate against the vocabulary, evidence files and acquisition targets', ok=rc == 0, detail=(out or err).splitlines()[0][:200] if (out or err) else ''))
    rc, out, err = sh(sys.executable, os.path.join(ROOT, 'inventory/check_repo_policy.py'))
    vis = next((x for x in out.splitlines() if x.startswith('visibility:')), 'visibility: unknown')
    res.append(dict(check='repository policy: size, personal-genome file types, credential patterns', ok=rc == 0 and 'policy check: PASS' in out, detail=vis))
    res.append(dict(check='all scripts compile', ok=all(s['compiles'] in (True, None) for s in scripts()), detail=f'{len(scripts())} scripts'))
    return res, vis


def markdown(a):
    L = ['# Repository audit', '', f"Generated {a['generated_utc']} at commit `{a['commit'][:8]}` by `inventory/repository_audit.py` (regenerate rather than edit). "
         'It records what exists, which lifecycle stage each dataset has reached, whether the checks pass, and whether the documentation says what the implementation does.', '']
    L += ['## Checks', '', '| Check | Result | Detail |', '|---|---|---|']
    for c in a['checks']:
        L.append(f"| {c['check']} | {'PASS' if c['ok'] else 'FAIL'} | {c['detail']} |")
    cl = collections.Counter(c['status'] for c in a['claims'])
    L += ['', f"## Documentation claims checked against the implementation: {cl['PASS']} pass, {cl['FAIL']} fail, {cl['SKIP']} skipped", '',
          '| Document | Claim (pattern) | Result | Detail |', '|---|---|---|---|']
    for c in a['claims']:
        L.append(f"| `{c['doc']}` | `{c['pattern'][:70]}` | {c['status']} | {c['detail'][:120]} |")
    L += ['', '## Branches', '', '| Branch | Tip | Date | Tracked files | Ahead of / behind origin/main | Last commit |', '|---|---|---|---|---|---|']
    for b in a['branches']:
        L.append(f"| `{b['branch']}` | `{b['tip']}` | {b['date']} | {b['tracked_files']} | +{b['ahead_of_main']} / -{b['behind_main']} | {b['subject']} |")
    s = a['dataset_summary']
    L += ['', '## Datasets (acquisition targets)', '',
          f"{s['listed']} listed · {s['acquired']} acquired (manifest ok and the file present with the recorded size) · {s['verified']} verified (validator passed) · {s['mapped']} mapped (a passing join or consistency test used it, or the built index holds an evidence row of it mapped to a normalised entity) · {s['integrated']} integrated "
          '(the built index holds a release of the target; `data/index` is git-ignored, so a checkout without it reports none).', '',
          '| Target | Source record | Tier | Listed | Acquired | Verified | Mapped | Integrated | Size |', '|---|---|---|---|---|---|---|---|---|']
    for d in a['datasets']:
        mk = lambda v: '✓' if v else '·'
        size = f"{d['bytes'] / 1e6:.1f} MB" if d['bytes'] else ''
        L.append(f"| `{d['id']}` | `{d['source']}` | {d['tier']} | {mk(d['listed'])} | {mk(d['acquired'])}{' (' + d['acquired_detail'][:40] + ')' if d['acquired_detail'] else ''} | {mk(d['verified'])} | {mk(d['mapped'])} | {mk(d['integrated'])} | {size} |")
    L += ['', '## Manifests', '', '| File | Records | What it records |', '|---|---|---|']
    for m in a['manifests']:
        L.append(f"| `{m['path']}` | {m['records']} | {m['detail']} |")
    L += ['', '## Documents', '', '| File | Title | Lines | Last commit |', '|---|---|---|---|']
    for d in a['documents']:
        L.append(f"| `{d['path']}` | {d['title']} | {d['lines']} | {d['last_commit']} |")
    L += ['', '## Scripts', '', '| File | Lines | Compiles | Purpose |', '|---|---|---|---|']
    for sc in a['scripts']:
        L.append(f"| `{sc['path']}` | {sc['lines']} | {sc['compiles'] if sc['compiles'] is not None else 'n/a'} | {sc['purpose']} |")
    L.append('')
    return '\n'.join(L)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--check', action='store_true')
    ap.add_argument('--fast', action='store_true')
    a = ap.parse_args()
    checks, vis = run_checks(a.fast)
    claims = check_claims()
    survivors = selftest_claims()
    checks.append(dict(check='the claim checker detects a changed number in every claim (mutation self-test)', ok=not survivors, detail=f'{len(CLAIMS)} claims' if not survivors else f'undetected: {survivors}'))
    ds = dataset_lifecycle()
    out = dict(generated_utc=datetime.datetime.now(datetime.timezone.utc).strftime('%Y-%m-%dT%H:%M:%SZ'), commit=sh('git', 'rev-parse', 'HEAD')[1], visibility=vis, checks=checks, claims=claims,
               branches=branches(), datasets=ds, dataset_summary=dict(listed=len(ds), **{k: sum(1 for d in ds if d[k]) for k in ('acquired', 'verified', 'mapped', 'integrated')}),
               manifests=manifests(), documents=documents(), scripts=scripts())
    with open(os.path.join(HERE, 'REPOSITORY_AUDIT.json'), 'w') as f:
        json.dump(out, f, indent=1)
    with open(os.path.join(HERE, 'REPOSITORY_AUDIT.md'), 'w') as f:
        f.write(markdown(out))
    cl = collections.Counter(c['status'] for c in claims)
    bad = [c for c in checks if not c['ok']] + [c for c in claims if c['status'] == 'FAIL']
    print(f"checks: {sum(1 for c in checks if c['ok'])}/{len(checks)} pass | claims: {dict(cl)} | datasets: {out['dataset_summary']} | {vis}")
    for c in bad:
        print('  FAIL', c.get('check') or (c['doc'], c['pattern'][:60]), '-', c['detail'])
    if a.check and bad:
        sys.exit(1)


if __name__ == '__main__':
    main()
