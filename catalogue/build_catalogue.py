"""Validate catalogue/sources/*.json and generate the inventory, coverage matrix, gap matrix and reuse plan.

Lifecycle stages are DERIVED from the acquisition record, never asserted, and are evaluated PER DIMENSION so that a
source is credited only for the evidence its downloaded files actually support (each acquisition target lists the
dimensions it supports; each compatibility test lists the targets it exercised):

  catalogued  a record exists
  accessible  status.accessible is yes/partial (from live probes recorded in evidence/)
  acquired    every target supporting the dimension is downloaded and verified in acquisition/manifest.jsonl
              ('partial' if only some are); existing dc assets count via inventory/EXISTING_ASSETS.json
  validated   all acquired targets passed their checks in acquisition/validation_report.json
  mapped      a compatibility test that exercised a supporting target passed (acquisition/compat_report.json)
  usable      acquired AND validated AND mapped, with a licence verified live that permits research use
              (computed; there is no hand-set flag, and it is not a claim that an analysis pipeline exists)

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
# catalogue source -> (existing asset id, dimensions that asset covers, is the base of the join tests)
EXISTING_SOURCES = {'mavedb': (['mavedb-human'], ['direct_experimental_variant_effects'], True),
                    'ensembl_genomes': (['ensembl-human-genomes'], ['reference_annotation'], False),
                    'nhgri_education': (['nhgri-gwas', 'nhgri-1000-genomes'], ['background_documentation'], False)}
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
    for line in open(os.path.join(HERE, 'evidence', 'evidence_index.jsonl')):
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
    return (vocab, sources, ev, targets, manifest, jload(os.path.join(ROOT, 'acquisition', 'validation_report.json'), {}),
            jload(os.path.join(ROOT, 'acquisition', 'compat_report.json'), {}), jload(os.path.join(ROOT, 'inventory', 'EXISTING_ASSETS.json'), {}),
            jload(os.path.join(HERE, 'gaps.json'), {}))


def validate(vocab, sources, ev, targets):
    errs = []
    dims = {d['key'] for d in vocab['dimensions']}
    ids = collections.Counter(s.get('id') for s in sources)
    for s in sources:
        i = s.get('id', '?')
        miss = [k for k in REQUIRED if k not in s]
        if miss:
            errs.append(f'{i}: missing fields {miss}')
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
        for blk in (s['licence'], s['version']):
            if blk.get('verified') not in vocab['verification']:
                errs.append(f'{i}: bad verified flag {blk.get("verified")}')
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
            else:
                miss = [d for d in targets[t].get('supports', []) if d not in s['categories']]
                if miss:
                    errs.append(f'{i}: target {t} supports {miss}, which are not among the record categories')
        if s['licence'].get('verified') == 'unverified' and 'licence' not in ' '.join(s['unverified_fields']).lower():
            errs.append(f'{i}: licence is unverified but not listed in unverified_fields')
        if s['version'].get('verified') == 'unverified' and not any('version' in u or 'release' in u for u in s['unverified_fields']):
            errs.append(f'{i}: version is unverified but not listed in unverified_fields')
        for k in ('usable', 'mapped_to_dc'):
            if k in s['status']:
                errs.append(f'{i}: status.{k} is hand-set; it is computed from the acquisition record')
    claimed = {t for s in sources for t in s.get('acquisition_targets', [])}
    for tid, t in targets.items():
        if tid not in claimed:
            errs.append(f'target {tid} is not claimed by any source record')
        if 'supports' not in t:
            errs.append(f'target {tid} has no "supports" list')
        for d in t.get('supports', []):
            if d not in dims:
                errs.append(f'target {tid} supports unknown dimension {d}')
    return errs


def lic_ok(s):
    return s['licence']['verified'] == 'live' and s['licence']['access'] in ('open', 'open_with_terms')


def lifecycle(s, dim, targets, manifest, validation, compat, existing):
    st = dict.fromkeys(STAGES, False)
    notes = []
    st['catalogued'] = True
    st['accessible'] = s['status']['accessible'] in ('yes', 'partial')
    tl = [t for t in s['acquisition_targets'] if dim is None or dim in targets[t].get('supports', [])]
    ok = [t for t in tl if manifest.get(t, {}).get('ok')]
    if tl:
        st['acquired'] = True if len(ok) == len(tl) else ('partial' if ok else False)
        if ok:
            vals = [validation.get(t, {}).get('passed') for t in ok]
            st['validated'] = bool(st['acquired'] is True and all(v is True for v in vals))
            bad = [t for t, v in zip(ok, vals) if v is False]
            if bad:
                notes.append('validation failures: ' + ', '.join(bad))
    base = False
    tests = compat.get('sources', {}).get(s['id'], {}).get('tests', [])
    ex = EXISTING_SOURCES.get(s['id'])
    if ex and (dim is None or dim in ex[1]):
        status = [next((a['status'] for a in existing.get('assets', []) if a['id'] == aid), None) for aid in ex[0]]
        if all(x == 'downloaded_verified' for x in status):
            st['acquired'] = True
            st['validated'] = True
            # the existing asset is the join base for ex[2]; otherwise any passed test of this source (they name no download target) maps it
            base = ex[2] or any(t.get('passed') is True for t in tests)
            notes.append('existing assets ' + ', '.join(ex[0]) + ': downloaded_verified')
    if dim is None:
        st['mapped'] = bool(base or any(t.get('passed') is True for t in tests))
    else:
        tset = set(tl)
        st['mapped'] = bool(base or any(t.get('passed') is True and tset & set(t.get('targets', [])) for t in tests))
    st['usable'] = bool(st['acquired'] is True and st['validated'] is True and st['mapped'] and lic_ok(s))
    return st, notes


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
    dim_label = {d['key']: d for d in dims}
    now = datetime.datetime.now(datetime.timezone.utc).strftime('%Y-%m-%d %H:%M UTC')
    life = {}      # (source, dim|None) -> (stages, notes)
    for s in sources:
        life[(s['id'], None)] = lifecycle(s, None, targets, manifest, validation, compat, existing)
        for d in s['categories']:
            life[(s['id'], d)] = lifecycle(s, d, targets, manifest, validation, compat, existing)
    usable_dims = {s['id']: [dim_label[d]['id'] for d in s['categories'] if life[(s['id'], d)][0]['usable']] for s in sources}

    cat = dict(generated_utc=now, sources=[dict(s, lifecycle_all_targets=life[(s['id'], None)][0], lifecycle_by_dimension={d: life[(s['id'], d)][0] for d in s['categories']},
                                                lifecycle_notes=life[(s['id'], None)][1]) for s in sources])
    json.dump(cat, open(os.path.join(HERE, 'catalogue.json'), 'w'), indent=1)
    gen = (f'*Generated {now} by `catalogue/build_catalogue.py` from `catalogue/sources/*.json`, `acquisition/targets.json`, `acquisition/manifest.jsonl`, '
           '`acquisition/validation_report.json` and `acquisition/compat_report.json`. Edit the sources, not this file.*')
    legend = ('`✓` done · `~` partial · `·` not done. Stages: **cat**alogued, **acc**essible from this environment, **acq**uired and checksum-verified, '
              '**val**idated against the source description, **map**ped to dc assets in a demonstrated join, **use**able (acquired + validated + mapped + licence verified). '
              'Stages beyond *acc* are evaluated per dimension, using only the downloaded files that support that dimension.')

    # ---- MASTER INVENTORY ---------------------------------------------------------------------------------
    by_dim = collections.defaultdict(list)
    for s in sources:
        by_dim[s['categories'][0]].append(s)
    L = ['# Master inventory of reusable genomic resources', '', gen, '', legend, '',
         f'{len(sources)} sources catalogued. A source appears under its primary evidence dimension. Columns *acq/val/map* here cover all of the source\'s downloaded files; '
         '"Usable for" lists the dimensions for which the source is actually usable. See `COVERAGE_MATRIX.md` for the per-dimension view.', '']
    for d in dims:
        items = by_dim.get(d['key'], [])
        if not items:
            continue
        L += [f"## {d['id']} {d['label']}", '', '| Source | Version / date | Access · licence | Variants · assemblies | cat | acc | acq | val | map | Usable for | Decision |', '|---|---|---|---|---|---|---|---|---|---|---|']
        for s in sorted(items, key=lambda x: x['id']):
            st = life[(s['id'], None)][0]
            lic = s['licence']
            lic_s = f"{lic['access']} · {lic['name'][:70]}{'' if lic['verified'] == 'live' else ' (unverified)'}"
            ver = s['version']
            ver_s = f"{ver['value'][:60]}{'' if ver['verified'] == 'live' else ' (unverified)'} · {ver['release_date']}"
            L.append(f"| **{s['name']}** (`{s['id']}`) | {ver_s} | {lic_s} | {', '.join(s['variant_classes']) or 'n/a'} · {', '.join(s['assemblies'])} | "
                     + ' | '.join(mark(st[k]) for k in STAGES[:5]) + f" | {', '.join(usable_dims[s['id']]) or '—'} | {s['reuse_decision']['action']} |")
        L.append('')
    L += ['## Existing dc assets (preserved, unchanged)', '', 'See `inventory/EXISTING_ASSETS.md`. They are referenced here through `overlap_with_existing`.', '']
    open(os.path.join(HERE, 'MASTER_INVENTORY.md'), 'w').write('\n'.join(L))

    # ---- COVERAGE MATRIX ----------------------------------------------------------------------------------
    L = ['# Coverage matrix: evidence dimensions by lifecycle stage', '', gen, '', legend, '',
         'Counts are numbers of sources reaching the stage **for that dimension**. A dimension is not complete because a source is catalogued, or even because some of its files are downloaded.', '',
         '| Dimension | Priority | Sources | cat | acc | acq | val | map | use |', '|---|---|---|---|---|---|---|---|---|']
    for d in dims:
        ss = [s for s in sources if d['key'] in s['categories']]
        cnt = {k: sum(1 for s in ss if life[(s['id'], d['key'])][0][k] is True) for k in STAGES}
        L.append(f"| {d['id']} {d['label']} | {d['priority']} | {len(ss)} | " + ' | '.join(str(cnt[k]) for k in STAGES) + ' |')
    L.append('')
    for d in dims:
        ss = [s for s in sources if d['key'] in s['categories']]
        L += [f"### {d['id']} {d['label']}", '', d['description'], '']
        if not ss:
            L += ['*No source catalogued.*', '']
            continue
        L += ['| Source | cat | acc | acq | val | map | use | Files supporting this dimension | Notes |', '|---|---|---|---|---|---|---|---|---|']
        for s in sorted(ss, key=lambda x: x['id']):
            st, notes = life[(s['id'], d['key'])]
            sup = [t for t in s['acquisition_targets'] if d['key'] in targets[t].get('supports', [])]
            got = [t for t in sup if manifest.get(t, {}).get('ok')]
            files = f"{len(got)}/{len(sup)} downloaded" if sup else 'none targeted'
            L.append(f"| `{s['id']}` | " + ' | '.join(mark(st[k]) for k in STAGES) + f" | {files} | {s['status']['notes'][:100]} {'; '.join(notes)[:80]} |")
        L.append('')
    open(os.path.join(HERE, 'COVERAGE_MATRIX.md'), 'w').write('\n'.join(L))

    # ---- GAP MATRIX ---------------------------------------------------------------------------------------
    L = ['# Gap matrix', '', gen, '', 'For each dimension: the best stage reached by any source, what is actually usable, what is missing and the next action. The first three columns are computed; the gap and next-action text is curated in `catalogue/gaps.json`.', '',
         '| Dimension | Priority | Best stage reached | Usable evidence now | Main gap | Next action |', '|---|---|---|---|---|---|']
    for d in dims:
        ss = [s for s in sources if d['key'] in s['categories']]
        best = 'none'
        for k in reversed(STAGES):
            if any(life[(s['id'], d['key'])][0][k] is True for s in ss):
                best = k
                break
        use = ', '.join(f"`{s['id']}`" for s in ss if life[(s['id'], d['key'])][0]['usable']) or '—'
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

    # ---- REUSE PLAN ---------------------------------------------------------------------------------------
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
            st, notes = life[(s['id'], None)]
            L += [f"### `{s['id']}`: {s['name']}", f"- **Why:** {s['reuse_decision']['rationale']}",
                  '- **Status (all files):** ' + ' '.join(f'{k}{mark(st[k])}' for k in STAGES[:5]) + f" · usable for: {', '.join(usable_dims[s['id']]) or 'nothing yet'}" + (' · ' + '; '.join(notes) if notes else ''),
                  f"- **Access:** {s['licence']['access']} ({s['licence']['name'][:120]})"]
            if s['acquisition_targets']:
                L.append('- **Targets:** ' + ', '.join(f"`{t}`" + ('' if manifest.get(t, {}).get('ok') else ' (not downloaded)') for t in s['acquisition_targets']))
            for d in s.get('deferred', []):
                L.append(f"- **Deferred ({d['approx_gib']} GiB):** {d['item']}: {d['reason']}")
            if s['unverified_fields']:
                L.append('- **Unverified:** ' + ', '.join(s['unverified_fields']))
            L.append('')
    open(os.path.join(HERE, 'REUSE_PLAN.md'), 'w').write('\n'.join(L))
    summary = {k: sum(1 for s in sources if life[(s['id'], None)][0][k] is True) for k in STAGES}
    print('sources reaching each stage (all files):', summary)
    print('sources usable for at least one dimension:', sum(1 for s in sources if usable_dims[s['id']]))
    return 0


if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('--check', action='store_true')
    sys.exit(build(ap.parse_args().check))
