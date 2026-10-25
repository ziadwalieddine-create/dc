"""Validate catalogue/sources/*.json and generate the inventory, coverage matrix, gap matrix and reuse plan.

Lifecycle stages are DERIVED, not asserted:
  catalogued  a record exists
  accessible  status.accessible is yes/partial (set from live probes recorded in evidence/)
  acquired    every acquisition target of the source is downloaded and verified in acquisition/manifest.jsonl
              (or the source is an existing dc asset whose audit status is downloaded_verified)
  validated   every acquired target passed its checks in acquisition/validation_report.json
  mapped      the source appears as 'mapped' in acquisition/compat_report.json (a demonstrated join to dc assets)
  usable      mapped AND validated AND status.usable is true in the record

  python3 build_catalogue.py [--check]     # --check validates only and exits non-zero on problems
"""
import argparse
import collections
import datetime
import glob
import json
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, '..'))
EXISTING_ASSET_IDS = ['nhgri-gwas', 'nhgri-1000-genomes', 'ensembl-human-genomes', 'mavedb-human']
EXISTING_SOURCE_MAP = {'mavedb': 'mavedb-human', 'ensembl_genomes': 'ensembl-human-genomes', 'nhgri_education': None}
REQUIRED = ['id', 'name', 'categories', 'maintainer', 'url', 'provenance', 'version', 'scope', 'variant_classes', 'assemblies', 'formats',
            'identifiers', 'mapping_methods', 'quantitative', 'licence', 'update_frequency', 'coverage_limitations', 'known_omissions',
            'overlap_with_existing', 'reuse_decision', 'status', 'acquisition_targets', 'evidence', 'unverified_fields']
STAGES = ['catalogued', 'accessible', 'acquired', 'validated', 'mapped', 'usable']


def jload(p, default=None):
    return json.load(open(p)) if os.path.exists(p) else default


def load_all():
    vocab = jload(os.path.join(HERE, 'vocab.json'))
    sources = [json.load(open(p)) for p in sorted(glob.glob(os.path.join(HERE, 'sources', '*.json')))]
    ev = {}
    idx = os.path.join(HERE, 'evidence', 'evidence_index.jsonl')
    for line in open(idx):
        r = json.loads(line)
        if r.get('excerpt'):
            ev[f"{r['source']}/{r['name']}"] = r
    targets = {t['id']: t for t in jload(os.path.join(ROOT, 'acquisition', 'targets.json'), [])}
    manifest = {}
    mp = os.path.join(ROOT, 'acquisition', 'manifest.jsonl')
    if os.path.exists(mp):
        for line in open(mp):
            r = json.loads(line)
            manifest[r['id']] = r
    validation = jload(os.path.join(ROOT, 'acquisition', 'validation_report.json'), {})
    compat = jload(os.path.join(ROOT, 'acquisition', 'compat_report.json'), {})
    existing = jload(os.path.join(ROOT, 'inventory', 'EXISTING_ASSETS.json'), {})
    gaps = jload(os.path.join(HERE, 'gaps.json'), {})
    return vocab, sources, ev, targets, manifest, validation, compat, existing, gaps


def validate(vocab, sources, ev, targets):
    errs = []
    dims = {d['key'] for d in vocab['dimensions']}
    ids = collections.Counter(s.get('id') for s in sources)
    for s in sources:
        i = s.get('id', '?')
        for k in REQUIRED:
            if k not in s:
                errs.append(f'{i}: missing field {k}')
        if errs and errs[-1].startswith(i + ': missing'):
            continue
        if ids[i] > 1:
            errs.append(f'{i}: duplicate id')
        for c in s['categories']:
            if c not in dims:
                errs.append(f'{i}: unknown category {c}')
        for v in s['variant_classes']:
            if v not in vocab['variant_classes']:
                errs.append(f'{i}: unknown variant class {v}')
        for a in s['assemblies']:
            if a not in vocab['assemblies']:
                errs.append(f'{i}: unknown assembly {a}')
        if s['licence'].get('access') not in vocab['access']:
            errs.append(f"{i}: bad licence.access {s['licence'].get('access')}")
        for blk, key in ((s['licence'], 'verified'), (s['version'], 'verified')):
            if blk.get(key) not in vocab['verification']:
                errs.append(f'{i}: bad verified flag {blk.get(key)}')
        if s['reuse_decision'].get('action') not in vocab['reuse_actions']:
            errs.append(f"{i}: bad reuse action {s['reuse_decision'].get('action')}")
        for o in s['overlap_with_existing']:
            if o not in EXISTING_ASSET_IDS:
                errs.append(f'{i}: unknown existing asset {o}')
        for e in s['evidence']:
            src, _, name = e.partition('/')
            if e not in ev and not os.path.exists(os.path.join(HERE, 'evidence', src, name + '.json')):
                errs.append(f'{i}: evidence not found: {e}')
        for t in s['acquisition_targets']:
            if t not in targets:
                errs.append(f'{i}: unknown acquisition target {t}')
        if s['licence'].get('verified') == 'unverified' and 'licence' not in ' '.join(s['unverified_fields']).lower():
            errs.append(f"{i}: licence is unverified but not listed in unverified_fields")
        if s['version'].get('verified') == 'unverified' and not any('version' in u or 'release' in u for u in s['unverified_fields']):
            errs.append(f"{i}: version is unverified but not listed in unverified_fields")
    claimed = {t for s in sources for t in s.get('acquisition_targets', [])}
    for t in targets:
        if t not in claimed:
            errs.append(f'target {t} is not claimed by any source record')
    return errs


def lifecycle(s, targets, manifest, validation, compat, existing):
    st = dict.fromkeys(STAGES, False)
    st['catalogued'] = True
    st['accessible'] = s['status']['accessible'] in ('yes', 'partial')
    tl = s['acquisition_targets']
    ok = [t for t in tl if manifest.get(t, {}).get('ok')]
    note = []
    if tl:
        st['acquired'] = 'partial' if 0 < len(ok) < len(tl) else (len(ok) == len(tl))
        if ok:
            vals = [validation.get(t, {}).get('passed') for t in ok]
            st['validated'] = all(vals) if all(v is not None for v in vals) else False
            if any(v is False for v in vals):
                note.append('validation failures: ' + ', '.join(t for t, v in zip(ok, vals) if v is False))
    ex = EXISTING_SOURCE_MAP.get(s['id'])
    if ex:
        status = next((a['status'] for a in existing.get('assets', []) if a['id'] == ex), None)
        if status == 'downloaded_verified':
            st['acquired'] = True if not tl or st['acquired'] is not False else st['acquired']
            st['validated'] = True if not tl or st['validated'] is not False else st['validated']
            note.append(f'existing asset {ex}: {status}')
    mapped = compat.get('sources', {}).get(s['id'], {})
    st['mapped'] = bool(mapped.get('mapped'))
    st['usable'] = bool(st['mapped'] and st['validated'] is True and s['status'].get('usable'))
    return st, note


def mark(v):
    return '✓' if v is True else ('~' if v == 'partial' else ('·' if v is False else '?'))


def build(check_only=False):
    vocab, sources, ev, targets, manifest, validation, compat, existing, gaps = load_all()
    errs = validate(vocab, sources, ev, targets)
    if errs:
        print('CATALOGUE VALIDATION ERRORS:')
        for e in errs:
            print(' -', e)
        return 1
    print(f'{len(sources)} source records valid; {len(targets)} acquisition targets; {len(ev)} evidence probes')
    if check_only:
        return 0
    dims = vocab['dimensions']
    now = datetime.datetime.now(datetime.timezone.utc).strftime('%Y-%m-%d %H:%M UTC')
    life = {}
    for s in sources:
        life[s['id']] = lifecycle(s, targets, manifest, validation, compat, existing)
    # consolidated JSON
    cat = dict(generated_utc=now, sources=[dict(s, lifecycle=life[s['id']][0], lifecycle_notes=life[s['id']][1]) for s in sources])
    json.dump(cat, open(os.path.join(HERE, 'catalogue.json'), 'w'), indent=1)
    gen = f'*Generated {now} by `catalogue/build_catalogue.py` from `catalogue/sources/*.json`, `acquisition/manifest.jsonl`, `acquisition/validation_report.json` and `acquisition/compat_report.json`. Edit the sources, not this file.*'
    legend = '`✓` done · `~` partial · `·` not done. Stages: **cat**alogued, **acc**essible from this environment, **acq**uired and checksum-verified, **val**idated against the source description, **map**ped to dc assets in a demonstrated join, **use**able.'

    # ---- MASTER INVENTORY -----------------------------------------------------------------------------------
    by_dim = collections.defaultdict(list)
    for s in sources:
        by_dim[s['categories'][0]].append(s)
    L = ['# Master inventory of reusable genomic resources', '', gen, '', legend, '',
         f'{len(sources)} sources catalogued. A source appears under its primary evidence dimension; the coverage matrix lists every dimension a source supports.', '']
    for d in dims:
        items = by_dim.get(d['key'], [])
        if not items:
            continue
        L += [f"## {d['id']} {d['label']}", '', '| Source | Version / date | Access · licence | Variants · assemblies | cat | acc | acq | val | map | use | Decision |', '|---|---|---|---|---|---|---|---|---|---|---|']
        for s in sorted(items, key=lambda x: x['id']):
            st = life[s['id']][0]
            lic = s['licence']
            lic_s = f"{lic['access']} · {lic['name'][:70]}{'' if lic['verified'] == 'live' else ' (unverified)'}"
            ver = s['version']
            ver_s = f"{ver['value'][:60]}{'' if ver['verified'] == 'live' else ' (unverified)'} · {ver['release_date']}"
            va = ', '.join(s['variant_classes']) or 'n/a'
            asm = ', '.join(s['assemblies'])
            L.append(f"| **{s['name']}** (`{s['id']}`) | {ver_s} | {lic_s} | {va} · {asm} | " + ' | '.join(mark(st[k]) for k in STAGES) + f" | {s['reuse_decision']['action']} |")
        L.append('')
    L += ['## Existing dc assets (preserved, unchanged)', '', 'See `inventory/EXISTING_ASSETS.md`. They are referenced here through `overlap_with_existing`.', '']
    open(os.path.join(HERE, 'MASTER_INVENTORY.md'), 'w').write('\n'.join(L))

    # ---- COVERAGE MATRIX ------------------------------------------------------------------------------------
    L = ['# Coverage matrix: evidence dimensions by lifecycle stage', '', gen, '', legend, '',
         'Counts are numbers of sources that support the dimension and have reached the stage. A dimension is **not** complete because a source is catalogued.', '',
         '| Dimension | Priority | Sources | cat | acc | acq | val | map | use |', '|---|---|---|---|---|---|---|---|---|']
    for d in dims:
        ss = [s for s in sources if d['key'] in s['categories']]
        cnt = {k: sum(1 for s in ss if life[s['id']][0][k] is True) for k in STAGES}
        L.append(f"| {d['id']} {d['label']} | {d['priority']} | {len(ss)} | " + ' | '.join(str(cnt[k]) for k in STAGES) + ' |')
    L.append('')
    for d in dims:
        ss = [s for s in sources if d['key'] in s['categories']]
        L += [f"### {d['id']} {d['label']}", '', d['description'], '']
        if not ss:
            L += ['*No source catalogued.*', '']
            continue
        L += ['| Source | cat | acc | acq | val | map | use | Notes |', '|---|---|---|---|---|---|---|---|']
        for s in sorted(ss, key=lambda x: x['id']):
            st, notes = life[s['id']]
            L.append(f"| `{s['id']}` | " + ' | '.join(mark(st[k]) for k in STAGES) + f" | {s['status']['notes'][:110]} {'; '.join(notes)[:90]} |")
        L.append('')
    open(os.path.join(HERE, 'COVERAGE_MATRIX.md'), 'w').write('\n'.join(L))

    # ---- GAP MATRIX -----------------------------------------------------------------------------------------
    L = ['# Gap matrix', '', gen, '', 'For each dimension: what exists, what is actually usable, what is missing and the next action. Curated text comes from `catalogue/gaps.json`; the status column is computed.', '',
         '| Dimension | Priority | Best stage reached | Usable evidence now | Main gap | Next action |', '|---|---|---|---|---|---|']
    for d in dims:
        ss = [s for s in sources if d['key'] in s['categories']]
        best = 'none'
        for k in reversed(STAGES):
            if any(life[s['id']][0][k] is True for s in ss):
                best = k
                break
        use = ', '.join(f"`{s['id']}`" for s in ss if life[s['id']][0]['usable']) or '—'
        g = gaps.get('dimensions', {}).get(d['key'], {})
        L.append(f"| {d['id']} {d['label']} | {d['priority']} | {best} | {use} | {g.get('gap', 'not curated yet')} | {g.get('next', '—')} |")
    L += ['', '## Resources that could not be accessed or reused as-is', '', '| Source | Reason | What is needed |', '|---|---|---|']
    for s in sorted(sources, key=lambda x: x['id']):
        if s['licence']['access'] in ('registration_required', 'restricted_agreement', 'blocked_for_automation', 'unverified') or s['status']['accessible'] in ('no', 'partial'):
            L.append(f"| `{s['id']}` | {s['status']['notes'][:160]} | {gaps.get('needs', {}).get(s['id'], s['reuse_decision']['rationale'][:140])} |")
    ns = gaps.get('not_yet_assessed', [])
    if ns:
        L += ['', '## Not yet assessed', '', 'Candidates for later evaluation (no verified facts recorded yet):', ''] + [f'- {x}' for x in ns]
    open(os.path.join(HERE, 'GAP_MATRIX.md'), 'w').write('\n'.join(L) + '\n')

    # ---- REUSE PLAN -----------------------------------------------------------------------------------------
    order = ['reuse_existing', 'adopt_curated', 'extend_integration', 'retrieve_source', 'build_original', 'defer', 'do_not_acquire']
    titles = {'reuse_existing': 'Reuse what is already verified', 'adopt_curated': 'Adopt an existing curated integration', 'extend_integration': 'Extend an existing integration',
              'retrieve_source': 'Retrieve the source data', 'build_original': 'Build original collection (only where justified)', 'defer': 'Deferred', 'do_not_acquire': 'Not acquired by design'}
    L = ['# Source-by-source reuse, acquisition and integration plan', '', gen, '',
         'Decision hierarchy: reuse existing verified material → adopt existing curated datasets → extend existing integrations → retrieve uncovered source data → build original methods only where justified.', '']
    for a in order:
        ss = [s for s in sources if s['reuse_decision']['action'] == a]
        if not ss:
            continue
        L += [f'## {titles[a]} ({len(ss)})', '']
        for s in sorted(ss, key=lambda x: x['id']):
            st, notes = life[s['id']]
            L += [f"### `{s['id']}`: {s['name']}", f"- **Why:** {s['reuse_decision']['rationale']}",
                  f"- **Status:** " + ' '.join(f'{k}{mark(st[k])}' for k in STAGES) + (' · ' + '; '.join(notes) if notes else ''),
                  f"- **Access:** {s['licence']['access']} ({s['licence']['name'][:120]})"]
            if s['acquisition_targets']:
                L.append('- **Targets:** ' + ', '.join(f"`{t}`" + ('' if manifest.get(t, {}).get('ok') else ' (not downloaded)') for t in s['acquisition_targets']))
            if s['unverified_fields']:
                L.append('- **Unverified:** ' + ', '.join(s['unverified_fields']))
            L.append('')
    open(os.path.join(HERE, 'REUSE_PLAN.md'), 'w').write('\n'.join(L))
    summary = {k: sum(1 for s in sources if life[s['id']][0][k] is True) for k in STAGES}
    print('stage counts:', summary)
    return 0


if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('--check', action='store_true')
    sys.exit(build(ap.parse_args().check))
