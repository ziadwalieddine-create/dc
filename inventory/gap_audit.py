#!/usr/bin/env python3
"""Adversarial gap audit of the phase requirements: every check either passes, fails, is skipped for a stated missing input, or is blocked on someone outside the session.

    python3 inventory/gap_audit.py [--index data/index] [--results data/results/HG002] [--no-tests]

Writes inventory/PHASE_AUDIT.md and PHASE_AUDIT.json and exits 1 if any check FAILS. A check that needs an artefact this checkout does not hold (the built index, a sample run,
the Python packages of the test suite) is reported SKIP with what it needs, never PASS. BLOCKED means only a person outside the session can settle it.

The checks attack the work from the outside rather than restate it: the requirements map must be followable, the documentation must not name files or commands that do not
exist, every number the audit log quotes is recomputed, no personal-genome file type is tracked now or anywhere in the history, the built index must have passed validation,
provenance verification, an idempotence comparison and the bcftools cross-checks, and every unresolved locus of the sample run must say what would settle it.
"""
import argparse
import datetime
import glob
import importlib
import json
import os
import re
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
VI = os.path.join(ROOT, 'variant_index')
sys.path.insert(0, HERE)

PERSONAL = re.compile(r'\.(vcf|vcf\.gz|gvcf|gvcf\.gz|bcf|bam|cram|sam|fastq|fastq\.gz|fq|fq\.gz|23andme\.txt)$', re.I)
CHECKS = []


def check(name):
    def deco(fn):
        CHECKS.append((name, fn))
        return fn
    return deco


class Skip(Exception):
    pass


class Blocked(Exception):
    pass


def sh(*args, cwd=ROOT, env=None):
    r = subprocess.run(args, cwd=cwd, capture_output=True, text=True, env=env)
    return r.returncode, r.stdout.strip(), r.stderr.strip()


def jload(path):
    with open(path) as f:
        return json.load(f)


def need(path):
    if not os.path.exists(path):
        raise Skip(f'needs {os.path.relpath(path, ROOT)}')
    return path


def can_import(*mods):
    for m in mods:
        try:
            importlib.import_module(m)
        except ImportError as e:
            raise Skip(f'needs the Python package {e.name}')


# ---- requirements and traceability ----------------------------------------------------------------------------------------------------------------
@check('requirements map is followable (paths and test modules exist, partial / missing statuses state their gap)')
def _requirements(a):
    rc, out, err = sh(sys.executable, os.path.join(VI, 'check_requirements.py'))
    if rc:
        raise AssertionError(out or err)
    return out.splitlines()[-1]


@check('REQUIREMENTS.md is the rendering of requirements.json (nobody edited it by hand)')
def _requirements_rendered(a):
    sys.path.insert(0, VI)
    import check_requirements as cr
    with open(os.path.join(VI, 'REQUIREMENTS.md')) as f:
        if f.read() != cr.render(cr.load()):
            raise AssertionError('REQUIREMENTS.md differs from the rendering of requirements.json; run check_requirements.py --write')
    return 'identical'


@check('no module or command-line tool of variant_index is outside the requirements map')
def _orphans(a):
    req = jload(os.path.join(VI, 'requirements.json'))['requirements']
    named = {c for r in req for c in r.get('code', [])}
    code = [os.path.relpath(p, ROOT) for p in glob.glob(os.path.join(VI, 'vindex', '*.py')) + glob.glob(os.path.join(VI, 'vindex', 'build', '*.py')) + glob.glob(os.path.join(VI, '*.py'))
            if not p.endswith('__init__.py')]
    text = ' '.join(named)
    orphans = sorted(c for c in code if c not in named and os.path.basename(c) not in text)
    if orphans:
        raise AssertionError('not named by any requirement: ' + ', '.join(orphans))
    return f'{len(code)} files, all named'


@check('every test module is used by at least one requirement')
def _test_orphans(a):
    req = jload(os.path.join(VI, 'requirements.json'))['requirements']
    used = {t for r in req for t in r.get('tests', [])}
    mods = {os.path.splitext(os.path.basename(p))[0] for p in glob.glob(os.path.join(VI, 'tests', 'test_*.py'))}
    extra = sorted(mods - used)
    if extra:
        raise AssertionError('test modules no requirement cites: ' + ', '.join(extra))
    return f'{len(mods)} test modules, all cited'


# ---- tests ----------------------------------------------------------------------------------------------------------------------------------------
@check('variant_index unit tests pass (including resource-warning errors)')
def _vindex_tests(a):
    if a.no_tests:
        raise Skip('--no-tests')
    can_import('pyarrow', 'pysam')
    rc, out, err = sh(sys.executable, '-W', 'error::ResourceWarning', '-m', 'unittest', 'discover', cwd=os.path.join(VI, 'tests'))
    m = re.search(r'Ran (\d+) tests?', err + out)
    if rc:
        raise AssertionError((err + out)[-1500:])
    return f"{m.group(1) if m else '?'} tests OK"


@check('acquisition helper tests pass')
def _acq_tests(a):
    if a.no_tests:
        raise Skip('--no-tests')
    rc, out, err = sh(sys.executable, '-m', 'unittest', 'acquisition.test_helpers')
    if rc:
        raise AssertionError((err + out)[-1500:])
    return re.search(r'Ran (\d+) tests?', err + out).group(0)


# ---- documentation against implementation -----------------------------------------------------------------------------------------------------------
@check('every number the READMEs and the latest audit-log entry quote is recomputed and agrees')
def _claims(a):
    ra = importlib.import_module('repository_audit')
    res = ra.check_claims()
    bad = [c for c in res if c['status'] == 'FAIL']
    if bad:
        raise AssertionError('; '.join(f"{c['doc']}: {c['detail']}" for c in bad[:6]))
    survivors = ra.selftest_claims()
    if survivors:
        raise AssertionError(f'the claim checker misses a changed number in {survivors[:3]}')
    skipped = sum(1 for c in res if c['status'] == 'SKIP')
    return f"{sum(1 for c in res if c['status'] == 'PASS')} pass, {skipped} skipped (data absent)"


@check('scripts and files that the documentation tells a reader to run or open exist')
def _doc_references(a):
    files = sh('git', 'ls-files', '*.md')[1].split()
    missing = []
    rx = re.compile(r'(?:python3?\s+(?:-m\s+)?|`)((?:[A-Za-z0-9_.-]+/)*[A-Za-z0-9_.-]+\.(?:py|sh|md))(?=[`\s]|$)')
    roots = ['', 'variant_index', 'acquisition', 'catalogue', 'catalogue/derive', 'inventory', 'scrape']
    for f in files:
        if f.startswith('scrape/') or os.path.basename(f) in ('ACQUIRED_FILES.md',) or f.startswith('catalogue/') and f.endswith(('MATRIX.md', 'INVENTORY.md', 'PROFILES.md', 'PLAN.md', 'RESOURCES.md', 'ATTRIBUTION.md')):
            continue                      # generated listings of acquired data files, not instructions
        with open(os.path.join(ROOT, f), errors='replace') as fh:
            text = fh.read()
        for m in rx.finditer(text):
            p = m.group(1)
            if p.startswith(('/', 'http', '.venv')) or '*' in p:
                continue
            here = os.path.dirname(f)
            if not any(os.path.exists(os.path.join(ROOT, r, here, p)) or os.path.exists(os.path.join(ROOT, r, p)) for r in roots):
                missing.append(f'{f} -> {p}')
    missing = sorted(set(missing))
    if missing:
        raise AssertionError('; '.join(missing[:8]) + (f' (+{len(missing) - 8} more)' if len(missing) > 8 else ''))
    return f'{len(files)} documents'


@check('catalogue records validate and every acquisition target is claimed by a source record')
def _catalogue(a):
    rc, out, err = sh(sys.executable, os.path.join(ROOT, 'catalogue', 'build_catalogue.py'), '--check')
    if rc:
        raise AssertionError((out or err)[:600])
    return out.splitlines()[0]


@check('statements in the catalogue are traceable to stored evidence')
def _evidence_support(a):
    rc, out, err = sh(sys.executable, os.path.join(ROOT, 'catalogue', 'derive', 'evidence_support.py'))
    last = [x for x in out.splitlines() if x.strip()][-1] if out else err[:200]
    if rc or not last.startswith('0 of'):
        raise AssertionError(last)
    return last


# ---- preservation and privacy ---------------------------------------------------------------------------------------------------------------------
@check('NHGRI, Ensembl, 1000 Genomes and MaveDB assets are unchanged apart from documented exceptions')
def _preservation(a):
    rc, out, err = sh(sys.executable, os.path.join(HERE, 'audit_existing.py'), '--check')
    if rc:
        raise AssertionError((out or err)[:500])
    return out.splitlines()[0] if out else 'ok'


@check('no personal-genome file type is tracked, and none was ever added to any branch')
def _no_personal_files(a):
    tracked = [f for f in sh('git', 'ls-files')[1].split('\n') if PERSONAL.search(f)]
    if tracked:
        raise AssertionError('tracked: ' + ', '.join(tracked[:5]))
    hist = sh('git', 'log', '--all', '--name-only', '--pretty=format:', '--diff-filter=A')[1].split('\n')
    ever = sorted({f for f in hist if f and PERSONAL.search(f)})
    allowed = {f for f in ever if f.startswith('variant_index/tests/fixtures/')}
    if set(ever) - allowed:
        raise AssertionError('added in history: ' + ', '.join(sorted(set(ever) - allowed)[:5]))
    return f'{len(ever)} synthetic fixture file(s) in history, no other'


@check('nothing under data/ is tracked and results directories are git-ignored')
def _data_untracked(a):
    tracked = [f for f in sh('git', 'ls-files')[1].split('\n') if f.startswith('data/')]
    if tracked:
        raise AssertionError(f'{len(tracked)} tracked files under data/, for example {tracked[:3]}')
    rc, out, _ = sh('git', 'check-ignore', 'data/index', 'data/results', 'data/pgs_registry')
    if len(out.split()) != 3:
        raise AssertionError('data/index, data/results and data/pgs_registry are not all ignored')
    return 'ignored'


@check('repository policy: file sizes, personal-genome types, credential patterns')
def _policy(a):
    rc, out, err = sh(sys.executable, os.path.join(HERE, 'check_repo_policy.py'))
    if rc or 'policy check: PASS' not in out:
        raise AssertionError(out[-600:] or err)
    vis = next((x for x in out.splitlines() if x.startswith('visibility:')), 'visibility: unknown')
    return vis


@check("the repository 'dc' is private")
def _private(a):
    rc, out, err = sh(sys.executable, os.path.join(HERE, 'check_repo_policy.py'))
    vis = next((x for x in out.splitlines() if x.startswith('visibility:')), 'visibility: unknown')
    if 'public' in vis:
        raise Blocked(f'{vis}: only the owner can change repository visibility (the session proxy refuses settings writes)')
    if 'private' not in vis:
        raise Skip(vis)
    return vis


# ---- the built index and the sample run -------------------------------------------------------------------------------------------------------------
@check('the built index passed validation (vocabulary, ranges, keys, alleles, partitions, reference alleles)')
def _validation(a):
    p = need(os.path.join(a.index, 'reports', 'validation.json'))
    v = jload(p)
    bad = {k: c['violations'] for k, c in v['checks'].items() if c['violations']}
    if not v['ok'] or bad:
        raise AssertionError(f'violations: {bad}')
    if os.path.getmtime(p) < os.path.getmtime(os.path.join(a.index, 'index_manifest.json')):
        raise AssertionError('validation.json is older than the index it describes')
    return f"{v['evidence_rows_checked']:,} evidence rows, {v.get('entity_variant_ref_checked', 0):,} REF alleles read back"


@check('provenance: release linkage, manifest agreement, recomputed record digests')
def _provenance(a):
    p = need(os.path.join(a.index, 'provenance_report.json'))
    v = jload(p)
    if not v['ok']:
        raise AssertionError(json.dumps({k: v[k] for k in ('orphan_release_ids', 'release_issues')})[:600])
    n = sum(d['verified'] for d in v['digests'].values())
    return f"{n:,} sampled records re-derived their digests across {len(v['digests'])} releases; {len(v['no_resolver'])} release(s) without a resolver are listed, not counted"


@check('every release that holds evidence rows has a resolver, and the sampled records of all of them re-derived their digests')
def _every_release_resolved(a):
    v = jload(need(os.path.join(a.index, 'provenance_report.json')))
    if v['no_resolver']:
        raise AssertionError(f"releases without a resolver: {[r['target_id'] for r in v['no_resolver']]}")
    bad = {d['target_id']: (d['n_mismatched'], d['n_unresolved']) for d in v['digests'].values() if d['n_mismatched'] or d['n_unresolved']}
    if bad:
        raise AssertionError(f'mismatched or unresolved records: {bad}')
    return f"{len(v['digests'])} releases, {sum(d['sampled'] for d in v['digests'].values()):,} records"


@check('lineage keys are shared between databases that describe the same study (the GWAS Catalog and Open Targets count a paper once)')
def _lineage(a):
    cov = jload(need(os.path.join(a.index, 'reports', 'coverage.json')))
    lin = cov.get('lineage')
    if not lin:
        raise Skip('coverage.json has no lineage section (duckdb was not available when it was made)')
    pair = lin['source_pairs'].get('gwas_catalog|opentargets', 0)
    if pair < 1000:
        raise AssertionError(f'only {pair} lineage groups are written by both the GWAS Catalog and Open Targets (4,328 papers were found in both before the keys were unified)')
    if not lin['source_pairs'].get('clingen|clinvar'):
        raise AssertionError('no lineage group is shared by ClinGen and ClinVar')
    return f"{lin['shared_between_sources']:,} of {lin['groups']:,} groups are written by more than one source; GWAS Catalog and Open Targets share {pair:,}"


BUILD_CODE_FIXED = ('__init__.py', 'normalise.py', 'refgenome.py', 'schema.py')


def build_code_at(commit):
    """(digest, number of files) of the build modules as committed at `commit`, computed the way vindex.build.runner.build_code computes it for the files on disk; None when a file is absent."""
    import hashlib
    pkg = 'variant_index/vindex'
    rc, out, _ = sh('git', 'ls-tree', '--name-only', commit, pkg + '/build/')
    if rc:
        return None
    rels = sorted(os.path.relpath(x, pkg) for x in out.splitlines() if x.endswith('.py')) + list(BUILD_CODE_FIXED)
    h = hashlib.sha256()
    for rel in rels:
        r = subprocess.run(['git', 'show', f'{commit}:{pkg}/{rel}'], cwd=ROOT, capture_output=True)
        if r.returncode:
            return None
        h.update(rel.encode() + b'\0')
        h.update(hashlib.sha256(r.stdout).digest())
    return h.hexdigest()[:16], len(rels)


def commit_with_build_code(want, max_commits=None):
    """The newest commit whose build modules have the digest `want` (HEAD first, then the commits that touched them), or None."""
    paths = ['variant_index/vindex/build'] + ['variant_index/vindex/' + f for f in BUILD_CODE_FIXED]
    rc, out, _ = sh('git', 'log', '--format=%H', '--', *paths)
    candidates = ['HEAD'] + (out.splitlines() if rc == 0 else [])
    for c in candidates[:max_commits]:
        got = build_code_at(c)
        if got and got[0] == want:
            return c, got[1]
    return None


@check('two builds from the same inputs with different hash seeds are byte-identical, and they ran the code that built the index')
def _idempotence(a):
    p = need(os.path.join(ROOT, 'data', 'idempotence', 'idempotence.json'))
    v = jload(p)
    if not v['identical']:
        raise AssertionError(json.dumps({k: v[k] for k in ('only_in_first', 'only_in_second', 'digest_differs')})[:600])
    examined = (v.get('build_code') or {}).get('sha256_16')
    if not examined:
        raise AssertionError('the idempotence report does not name the build code it examined (it predates that field); run check_idempotence.py again')
    built_by = (jload(need(os.path.join(a.index, 'index_manifest.json'))).get('build_code') or {}).get('sha256_16')
    if examined != built_by:
        raise AssertionError(f'the idempotence check examined build code {examined} but the index was built by {built_by}; run check_idempotence.py again')
    return f"{v['files']} files, limit {v['limit']} records per adapter, build code {built_by}"


@check('the index was built by committed code: its build-code digest is that of the build modules at a commit of this history')
def _index_code_committed(a):
    m = jload(need(os.path.join(a.index, 'index_manifest.json')))
    if sh('git', 'rev-parse', '--verify', 'HEAD')[0]:
        raise Skip('needs the git history of this repository')
    want = (m.get('build_code') or {}).get('sha256_16')
    if not want:
        raise AssertionError('the index manifest records no build_code (built by older code); rebuild the index')
    found = commit_with_build_code(want)
    if not found:
        raise AssertionError(f'no commit of this history has the build modules with digest {want}: the index was built from uncommitted or unknown code')
    commit, n = found
    head = build_code_at('HEAD')
    return f"commit {commit[:10] if commit != 'HEAD' else 'HEAD'} ({n} files)" + ('; HEAD has the same build modules' if head and head[0] == want else '; the build modules differ at HEAD')


@check('ClinVar-derived rows of the index equal the ClinVar VCF as read by bcftools')
def _clinvar_fields(a):
    v = jload(need(os.path.join(a.index, 'reports', 'clinvar_fields.json')))
    if not v['ok']:
        raise AssertionError(json.dumps(v['totals']) + ' ' + json.dumps(v['examples'])[:300])
    t = v['totals']
    return f"{t['vcf_records']:,} records, {t['vcf_classifications']:,} classifications, {t['mapped']:,} keys equal to bcftools norm"


@check('facts established outside the index (rsID positions, effect-allele orientation, classifications, star-allele definitions, curated assertions, retained unmapped rows) read back from it')
def _known_loci(a):
    v = jload(need(os.path.join(a.index, 'reports', 'known_loci.json')))
    if not v['ok']:
        bad = [f"{l['id']}/{c['id']}: {c['status']} {c['observed']}" for l in v['loci'] for c in l['checks'] if c['status'] != 'PASS']
        raise AssertionError('; '.join(bad)[:900])
    spec = jload(os.path.join(VI, 'known_loci.json'))
    wanted = {(l['id'], c['id']) for l in spec['loci'] for c in l['checks']}
    ran = {(l['id'], c['id']) for l in v['loci'] for c in l['checks']}
    if wanted != ran:
        raise AssertionError(f"the report was made from another version of known_loci.json ({len(wanted - ran)} checks of the spec are not in it, {len(ran - wanted)} checks of it are not in the spec); run verify_known_loci.py again")
    return f"{v['n_loci']} loci, {v['totals']['PASS']} checks, all pass"


@check('a sample run agrees with bcftools on its calls and on the ClinVar evidence of those calls')
def _crosscheck(a):
    v = jload(need(os.path.join(a.results, 'crosscheck.json')))
    if not v['ok']:
        raise AssertionError(json.dumps(v['calls']) + ' ' + json.dumps(v.get('evidence')) + ' ' + json.dumps(v['examples'])[:300])
    return f"{v['calls']['matcher']:,} calls; {v['evidence']['matcher']:,} with ClinVar evidence"


READS_BOUNDS = dict(false_alternate_among_absent=0.01, states_the_truth_at_positives=0.95)


@check('the read-evidence resolver was checked on real reads of a benchmark sample with the code as it is now, and its rates are inside the stated bounds')
def _reads_validation(a):
    import hashlib
    rec = jload(need(os.path.join(VI, 'reads_validation.json')))
    with open(os.path.join(VI, 'vindex', 'reads.py'), 'rb') as f:
        current = hashlib.sha256(f.read()).hexdigest()
    if rec.get('reads_module_sha256') != current:
        raise AssertionError('vindex/reads.py has changed since the validation record was made: run validate_reads.py again and merge the runs into reads_validation.json')
    if not rec.get('runs') or 'full' not in rec['runs']:
        raise AssertionError('the record holds no run at full depth')
    problems = []
    for label, run in sorted(rec['runs'].items()):
        for truth, c in run['by_truth'].items():
            if c.get('ERROR'):
                problems.append(f"{label}: {c['ERROR']} loci of {truth} could not be read")
        absent, n = run['rates']['truth_absent'], run['rates']['truth_absent']['loci']
        false_alt = run['by_truth']['truth_absent'].get('READS_ALTERNATE_PRESENT', 0)
        if n < 500 or false_alt / n > READS_BOUNDS['false_alternate_among_absent']:
            problems.append(f"{label}: {false_alt} of {n} loci absent from the benchmark were read as carrying the allele (bound {READS_BOUNDS['false_alternate_among_absent']:.0%}, and at least 500 loci)")
        for truth in ('truth_alt_het', 'truth_alt_hom'):
            r = run['rates'][truth]
            # the bound on how often the reads state the truth is for 30x and deeper (a 15x run abstains more and is reported); the bound on false claims is for every depth
            if r['loci'] < 300 or (run['fraction'] >= 0.5 and r['states_the_truth']['rate'] < READS_BOUNDS['states_the_truth_at_positives']):
                problems.append(f"{label}: {truth} states the truth at {r['states_the_truth']['rate']:.1%} of {r['loci']} loci (bound {READS_BOUNDS['states_the_truth_at_positives']:.0%} at half depth and above, and at least 300 loci)")
    if problems:
        raise AssertionError('; '.join(problems))
    full = rec['runs']['full']['rates']
    return (f"{', '.join(sorted(rec['runs']))}: at full depth {full['truth_alt_het']['states_the_truth']['rate']:.1%} of heterozygous and {full['truth_alt_hom']['states_the_truth']['rate']:.1%} of homozygous "
            f"benchmark alleles read as carried, {rec['runs']['full']['by_truth']['truth_absent'].get('READS_ALTERNATE_PRESENT', 0)} of {full['truth_absent']['loci']:,} absent ones read as carried")


@check('every unresolved high-value locus of the sample run states the evidence that would settle it')
def _missing_evidence(a):
    can_import('pyarrow')
    import pyarrow.parquet as pq
    rows = pq.read_table(need(os.path.join(a.results, 'unresolved_loci.parquet'))).to_pylist()
    empty = [r['key'] for r in rows if not (r['missing_evidence'] or '').strip()]
    if empty:
        raise AssertionError(f'{len(empty)} unresolved loci without missing_evidence, e.g. {empty[:3]}')
    return f'{len(rows):,} unresolved loci, all with text'


@check('the sample run accounts for every record it read (small calls, structural, repeat, rejected, no-call and reference counts add up)')
def _accounting(a):
    m = jload(need(os.path.join(a.results, 'run_manifest.json')))
    v = m['counts']['vcf']
    rows = m['rows']
    if rows['small_calls'] != m['counts']['calls']:
        raise AssertionError(f"small_calls rows {rows['small_calls']} != calls {m['counts']['calls']}")
    partial = v.get('small', 0) - m['counts']['calls']                    # alleles of a partly missing genotype that carry no copy of the allele are read but are not calls
    if partial < 0:
        raise AssertionError(f"reader produced {v.get('small')} small alleles but the table holds {m['counts']['calls']}")
    if rows['rejected_calls'] != v.get('rejected', 0) + v.get('unplaced_records', 0):
        raise AssertionError(f"rejected rows {rows['rejected_calls']} != rejected {v.get('rejected')} + unplaced {v.get('unplaced_records')}")
    if rows['repeat_calls'] < v.get('repeat_records', 0):
        raise AssertionError('fewer repeat rows than repeat records read')
    return f"{v.get('records', 0):,} records read; {m['counts']['calls']:,} small calls ({partial} read alleles carry no copy), {rows['sv_calls']:,} structural, {rows['rejected_calls']:,} rejected"


@check('the sample run was made by committed code that has not changed since, against the index that is on disk now')
def _run_provenance(a):
    import hashlib
    m = jload(need(os.path.join(a.results, 'run_manifest.json')))
    code = m.get('code') or {}
    if not code.get('commit'):
        raise AssertionError('the run manifest names no commit')
    if sh('git', 'cat-file', '-e', code['commit'] + '^{commit}')[0] != 0:
        raise AssertionError(f"commit {code['commit'][:10]} of the run is not in this repository")
    if code.get('uncommitted_changes_in_variant_index'):
        raise AssertionError('the run was made with uncommitted changes in variant_index/')
    if m.get('code_changed_during_run') is not False:
        raise AssertionError(f"code_changed_during_run is {m.get('code_changed_during_run')!r}: files under variant_index/ changed during the run, or the manifest predates the field")
    def sha(p):
        with open(p, 'rb') as f:
            return hashlib.sha256(f.read()).hexdigest()
    files = m.get('code_files') or {}
    if not files:
        raise AssertionError('the run manifest holds no per-module digests; rerun the sample')
    changed = sorted(k for k, d in files.items() if not os.path.exists(os.path.join(VI, k)) or sha(os.path.join(VI, k)) != d)
    if changed:
        raise AssertionError(f'{len(changed)} module(s) changed since the run (rerun the sample): {changed[:6]}')
    for name in ('index_manifest.json', 'table_digests.json'):
        recorded = (m.get('index') or {}).get(name + '_sha256')
        if not recorded or recorded != sha(os.path.join(a.index, name)):
            raise AssertionError(f'{name} of the index is not the one the run used')
    return f"commit {code['commit'][:10]}, {len(files)} modules unchanged, index manifests identical"


@check('the structural run on real calls (HG002 genotypes of the public HGSVC3 callset) accounts for every record and carries the summary of what it refused')
def _structural_run(a):
    can_import('pyarrow')
    import pyarrow.parquet as pq
    run_dir = a.results + '_sv'
    m = jload(need(os.path.join(run_dir, 'run_manifest.json')))
    v = m['counts']['vcf']
    n = m['rows']['sv_calls']
    if n + m['rows']['rejected_calls'] != v.get('records', 0) or n != sum(m['counts']['structural'].values()):
        raise AssertionError(f"records {v.get('records')}, structural rows {n}, rejected rows {m['rows']['rejected_calls']}, structural counts {m['counts']['structural']}")
    if m['code'].get('uncommitted_changes_in_variant_index') or m.get('code_changed_during_run') is not False:
        raise AssertionError('the structural run was made from a tree with uncommitted changes, or the tree changed during the run')
    t = pq.read_table(os.path.join(run_dir, 'sv_calls.parquet'), columns=['refused_summary_json', 'refused_json'])
    missing = sum(1 for x in t['refused_summary_json'].to_pylist() if x is None)
    biggest = max(len(x or '') for x in t['refused_json'].to_pylist())
    if missing or biggest > 100_000:
        raise AssertionError(f'{missing} calls without a refusal summary; the largest refused_json has {biggest:,} characters')
    return f"{n:,} structural calls in {m['seconds']:.0f} s; largest refusal list {biggest:,} characters"


@check('two matcher runs of one input with different hash seeds give identical result tables and manifests')
def _matcher_idempotence(a):
    v = jload(need(os.path.join(os.path.dirname(a.results), 'matcher_idempotence.json')))
    if not v['identical']:
        raise AssertionError(json.dumps({k: v[k] for k in ('tables_differ', 'only_in_a', 'only_in_b', 'manifest_fields_differ')})[:600])
    return f"{v['tables']} tables identical"


def run(a):
    out = []
    for name, fn in CHECKS:
        try:
            detail = fn(a)
            out.append(dict(check=name, status='PASS', detail=str(detail)))
        except Skip as e:
            out.append(dict(check=name, status='SKIP', detail=str(e)))
        except Blocked as e:
            out.append(dict(check=name, status='BLOCKED', detail=str(e)))
        except AssertionError as e:
            out.append(dict(check=name, status='FAIL', detail=str(e)[:900]))
        except Exception as e:  # noqa: BLE001
            out.append(dict(check=name, status='FAIL', detail=f'the check itself failed: {e!r}'[:900]))
    return out


def markdown(res, a):
    c = {s: sum(1 for r in res if r['status'] == s) for s in ('PASS', 'FAIL', 'SKIP', 'BLOCKED')}
    L = ['# Phase gap audit', '', f"Generated {datetime.datetime.now(datetime.timezone.utc).strftime('%Y-%m-%dT%H:%M:%SZ')} by `inventory/gap_audit.py` (regenerate rather than edit). "
         f"{c['PASS']} pass, {c['FAIL']} fail, {c['SKIP']} skipped, {c['BLOCKED']} blocked.", '',
         'PASS means the check ran and held. SKIP means an input this checkout does not hold was needed, so the check proves nothing here. BLOCKED means only a person outside the session can settle it.', '',
         '| Check | Result | Detail |', '|---|---|---|']
    L += [f"| {r['check']} | **{r['status']}** | {r['detail'].replace('|', '/').replace(chr(10), ' ')} |" for r in res]
    return '\n'.join(L) + '\n'


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--index', default=os.path.join(ROOT, 'data', 'index'))
    ap.add_argument('--results', default=os.path.join(ROOT, 'data', 'results', 'HG002'))
    ap.add_argument('--no-tests', action='store_true')
    a = ap.parse_args()
    res = run(a)
    with open(os.path.join(HERE, 'PHASE_AUDIT.json'), 'w') as f:
        json.dump(res, f, indent=1)
    with open(os.path.join(HERE, 'PHASE_AUDIT.md'), 'w') as f:
        f.write(markdown(res, a))
    for r in res:
        print(f"{r['status']:8} {r['check']}" + ('' if r['status'] == 'PASS' else f"\n         {r['detail'][:300]}"))
    return 1 if any(r['status'] == 'FAIL' for r in res) else 0


if __name__ == '__main__':
    sys.exit(main())
