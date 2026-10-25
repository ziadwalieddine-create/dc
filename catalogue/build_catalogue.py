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
            'overlap_with_existing', 'reuse_decision', 'status', 'acquisition_targets', 'evidence', 'unverified_fields', 'evidence_profile']
PROFILE_TEXT = ['biological_context', 'evidence_strength', 'mechanistic_confidence', 'phenotype_implication', 'absence_means']
GENE_LEVEL_CLASSES = {'gene_level', 'region_level', 'cnv', 'sv', 'str', 'complex_locus', 'mtdna'}
# the evidence kinds a mechanism basis may rest on
BASIS_KINDS = {'direct': {'direct_experimental', 'perturbation_experimental'}, 'curated': {'curated_expert_classification'}, 'association': {'statistical_association'},
               'prediction': {'computational_prediction'}, 'population': {'population_observation'}, 'case': {'case_observation'}}
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
    axes = {a['key'] for a in vocab['mechanism_axes']}
    all_ids = {s.get('id') for s in sources}
    declared = {}
    for s in sources:
        i = s.get('id', '?')
        pr = s.get('evidence_profile')
        if not isinstance(pr, dict):
            continue
        absent = [k for k in ('evidence_kind', 'allele_resolution', 'mechanism_evidence', 'overlaps_with_sources') + tuple(PROFILE_TEXT) if k not in pr]
        if absent:
            errs.append(f'{i}: evidence_profile lacks {absent}')
            continue
        kinds = pr['evidence_kind']
        if not kinds or any(k not in vocab['evidence_kinds'] for k in kinds):
            errs.append(f'{i}: evidence_profile.evidence_kind must be a non-empty list from the vocabulary: {kinds}')
            continue
        if pr['allele_resolution'] not in vocab['allele_resolutions']:
            errs.append(f"{i}: bad allele_resolution {pr['allele_resolution']}")
        for k in PROFILE_TEXT:
            if not isinstance(pr[k], str) or len(pr[k].strip()) < 8:
                errs.append(f'{i}: evidence_profile.{k} must state something')
        if pr['allele_resolution'] == 'none' and (pr['mechanism_evidence'] or not set(kinds) <= {'catalogue_metadata', 'documentation'}):
            errs.append(f'{i}: allele_resolution none is for catalogue and documentation records, which carry no mechanism evidence')
        for m in pr['mechanism_evidence']:
            if m.get('axis') not in axes:
                errs.append(f"{i}: unknown mechanism axis {m.get('axis')}")
            if m.get('basis') not in vocab['mechanism_bases']:
                errs.append(f"{i}: unknown mechanism basis {m.get('basis')}")
            vcs = m.get('variant_classes', [])
            if not vcs or any(v not in vocab['variant_classes'] for v in vcs):
                errs.append(f"{i}: mechanism entry {m.get('axis')} has bad variant_classes {vcs}")
            elif not set(vcs) <= set(s['variant_classes']):
                errs.append(f"{i}: mechanism entry {m.get('axis')} names variant classes {sorted(set(vcs) - set(s['variant_classes']))} that the record does not list")
            if m.get('dimension') is not None and m['dimension'] not in s['categories']:
                errs.append(f"{i}: mechanism entry {m.get('axis')} names dimension {m['dimension']}, which is not a category of the record")
            # a gene- or region-level resource must not be recorded as direct evidence about specific alleles or small variants
            if m.get('basis') == 'direct' and pr['allele_resolution'] in ('gene_level', 'region_level') and not set(vcs) <= GENE_LEVEL_CLASSES:
                errs.append(f"{i}: {pr['allele_resolution']} evidence cannot be direct evidence for {sorted(set(vcs) - GENE_LEVEL_CLASSES)} (a gene or region result does not describe a specific allele)")
            need = BASIS_KINDS.get(m.get('basis'))
            if need and not need & set(kinds):
                errs.append(f"{i}: mechanism basis {m.get('basis')} needs one of the evidence kinds {sorted(need)}")
        for c in pr.get('components', []):
            ck = c.get('evidence_kind')
            if not isinstance(c.get('name'), str) or len(c['name']) < 3 or not ck or any(k not in vocab['evidence_kinds'] for k in ck):
                errs.append(f"{i}: bad component {c.get('name')}")
                continue
            if not set(ck) <= set(kinds):
                errs.append(f"{i}: component {c['name']} has evidence kinds {sorted(set(ck) - set(kinds))} that the source profile does not list")
            if c.get('allele_resolution') not in vocab['allele_resolutions']:
                errs.append(f"{i}: component {c['name']} has bad allele_resolution {c.get('allele_resolution')}")
            o = c.get('origin', '')
            if not (o == 'native' or o in all_ids or o.startswith('external: ')):
                errs.append(f"{i}: component {c['name']} origin must be native, a catalogue source id or 'external: ...': {o!r}")
        if 'computational_prediction' in kinds and any(m['basis'] in ('direct', 'curated') for m in pr['mechanism_evidence']) and 'aggregated_integration' not in kinds:
            errs.append(f'{i}: a computational prediction resource cannot also carry direct or curated mechanism evidence unless it aggregates other sources')
        if set(kinds) <= {'reference_annotation', 'catalogue_metadata', 'documentation'} and any(m['basis'] != 'context' for m in pr['mechanism_evidence']):
            errs.append(f'{i}: reference, catalogue and documentation resources can only provide context')
        for o in pr['overlaps_with_sources']:
            if o.get('source') not in all_ids or o.get('source') == i:
                errs.append(f"{i}: overlap names unknown or own source {o.get('source')}")
                continue
            if o.get('relation') not in vocab['source_relations']:
                errs.append(f"{i}: overlap relation {o.get('relation')} not in the vocabulary")
            if o.get('relation') in ('identical_records', 'derived_from', 'may_conflict') and not o.get('note'):
                errs.append(f"{i}: overlap {o['relation']} with {o['source']} needs a note stating what is shared or contested")
            pair = frozenset((i, o['source']))
            if pair in declared:
                errs.append(f"{i}: overlap with {o['source']} is already declared by {declared[pair]}; declare each pair once")
            declared[pair] = i
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
            if not tl:                      # existing assets alone carry the dimension; when download targets also support it they must be ok as well
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
              'Stages beyond *acc* are evaluated per dimension, using only the downloaded files that support that dimension. '
              '`integrated` (a documented, tested query layer that returns records by variant, gene, protein or region) is defined in `vocab.json` and credited to no source, because no such layer exists yet.')

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
    held = [s for s in sorted(sources, key=lambda x: x['id']) if s['licence']['verified'] != 'live' and any(manifest.get(t, {}).get('ok') for t in s['acquisition_targets'])]
    L += ['', '## Held locally with a licence that is not verified', '',
          'These sources have files in `data/external/` (git-ignored) but no licence text was retrieved. Use them for local analysis with citation; do not redistribute until the terms are read.', '',
          '| Source | What the pages say | Files held |', '|---|---|---|']
    for s in held:
        got = [t for t in s['acquisition_targets'] if manifest.get(t, {}).get('ok')]
        L.append(f"| `{s['id']}` | {s['licence']['name'][:230]} | {len(got)} of {len(s['acquisition_targets'])} targets |")
    if not held:
        L.append('| none | | |')
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
    # ---- ACQUIRED FILES REGISTER --------------------------------------------------------------------------
    owner = {t: s['id'] for s in sources for t in s['acquisition_targets']}
    L = ['# Acquired files register', '', gen, '',
         'Every download target with a verified manifest record, as it stands in `data/external/` (git-ignored). *Source date* is the HTTP Last-Modified of the file where the server gave one, otherwise the time we fetched it. '
         '*SHA-256* is ours; the source publishes a checksum only where `acquisition/manifest.jsonl` says so. Row counts come from the validation report.', '',
         '| Target | Source record | Path under `data/external/` | Size | SHA-256 | Source date | Validated | Rows / records | Supports | Licence | Scope note |', '|---|---|---|---|---|---|---|---|---|---|---|']
    def human(n):
        return f'{n / 1e9:.2f} GB' if n >= 1e9 else (f'{n / 1e6:.1f} MB' if n >= 1e6 else f'{n / 1e3:.1f} KB')
    for tid in sorted(targets):
        m = manifest.get(tid)
        if not m or not m.get('ok'):
            continue
        t = targets[tid]
        v = validation.get(tid, {})
        st = v.get('stats', {})
        rows = next((st[k] for k in ('rows', 'records', 'n_accessions', 'accessions') if isinstance(st.get(k), int)), '')
        when = (m.get('last_modified') or m.get('fetched_utc') or '')[:31]
        sha = (m.get('sha256') or '')[:12]
        L.append(f"| `{tid}` | `{owner.get(tid, '?')}` | `{t['dest']}` | {human(m.get('bytes', 0))} | `{sha}` | {when} | {'✓' if v.get('passed') else ('✗' if v else '·')} | {rows:,}".replace(',', ',') if isinstance(rows, int) else
                 f"| `{tid}` | `{owner.get(tid, '?')}` | `{t['dest']}` | {human(m.get('bytes', 0))} | `{sha}` | {when} | {'✓' if v.get('passed') else ('✗' if v else '·')} | ")
        rec_s = next((x for x in sources if x['id'] == owner.get(tid)), None)
        lic_txt = (rec_s['licence']['name'][:90] + ('' if rec_s['licence']['verified'] == 'live' else ' (unverified)')) if rec_s else str(t.get('licence', ''))[:90]
        L[-1] += f" | {', '.join(t.get('supports', [])) or 'auxiliary'} | {lic_txt} | {t.get('scope_note', '')} |"
    open(os.path.join(HERE, 'ACQUIRED_FILES.md'), 'w').write('\n'.join(L) + '\n')
    # ---- ATTRIBUTION ---------------------------------------------------------------------------------------
    L = ['# Attribution and licence notes for the data we hold', '', gen, '',
         'One row per source whose files are in `data/external/` or whose tables are redistributed on the data branch. Licence and attribution wording is copied from the source records, which quote the pages they were read from; '
         'where a source asks for a particular citation, the source URL leads to it. A private repository would not remove any of these conditions.', '',
         '| Source | Licence | Attribution / redistribution note | Where the terms are stated |', '|---|---|---|---|']
    for s in sorted(sources, key=lambda x: x['id']):
        held = any(manifest.get(t, {}).get('ok') for t in s['acquisition_targets']) or s['id'] in EXISTING_SOURCES
        if not held:
            continue
        words = ('licen', 'polic', 'terms', 'disclaim', 'copyright', 'faq', 'download_dir', 'data_use', 'data_access', 'zenodo_record')
        urls = [ev[e]['url'] for e in s['evidence'] if e in ev and ev[e].get('http_status') == 200 and any(w in e.lower() for w in words)]
        where = ', '.join(f'<{u}>' for u in dict.fromkeys(urls))[:300] or f"<{s['url']}>"
        lic = s['licence']
        L.append(f"| `{s['id']}` | {lic['name'][:200]}{'' if lic['verified'] == 'live' else ' (unverified)'} | {lic['redistribution'][:260]} | {where} |")
    open(os.path.join(HERE, 'ATTRIBUTION.md'), 'w').write('\n'.join(L) + '\n')
    # ---- EVIDENCE PROFILES --------------------------------------------------------------------------------
    prof = {s['id']: s['evidence_profile'] for s in sources}
    L = ['# Evidence profiles', '', gen, '',
         'Every source record carries an `evidence_profile` that keeps apart what the resource is (a measurement, a perturbation, a population observation, an association, a curated classification, a prediction, a reference annotation), '
         'the resolution at which its records refer to alleles, the mechanisms it informs and on what basis, how strong the evidence is, how far a mechanistic reading can be trusted, what phenotype it implies and what absence means. '
         'These fields are independent: a curated classification can be strong and carry no mechanism; a strong association says nothing about the causal variant; a gene-level knockout or constraint metric is never allele-specific evidence; '
         'the absence of a study is never recorded as the absence of an effect. Hypotheses and competing explanations are listed per mechanism axis in `MECHANISM_MATRIX.md`.', '',
         'Rules enforced by `build_catalogue.py --check`: a gene- or region-level resource cannot be direct evidence for a specific allele or small variant; a prediction basis needs the `computational_prediction` kind; every mechanism basis needs a matching evidence kind; '
         'a prediction resource cannot also carry direct or curated mechanism evidence unless it aggregates other sources; reference, catalogue and documentation resources provide context only; mechanism entries cannot name variant classes the record does not list; '
         'every overlap between sources is declared once, and identical, derived or conflicting records must say what is shared or contested.', '',
         '## Evidence kinds', '', '| Kind | Meaning | Sources |', '|---|---|---|']
    for k, desc in vocab['evidence_kinds'].items():
        ss = [s['id'] for s in sources if k in prof[s['id']]['evidence_kind']]
        L.append(f"| `{k}` | {desc} | {', '.join(f'`{x}`' for x in ss) or '—'} |")
    L += ['', '## Allele resolution', '', 'The resolution at which a source\'s records refer to alleles. Gene-, region- and protein-level records describe a gene, interval or protein as a whole.', '', '| Resolution | Meaning | Sources |', '|---|---|---|']
    for k, desc in vocab['allele_resolutions'].items():
        ss = [s['id'] for s in sources if prof[s['id']]['allele_resolution'] == k]
        L.append(f"| `{k}` | {desc} | {', '.join(f'`{x}`' for x in ss) or '—'} |")
    L += ['', '## Per-source profiles', '']
    for s in sorted(sources, key=lambda x: x['id']):
        pr = prof[s['id']]
        me = '<br>'.join(f"`{m['axis']}` · {m['basis']} · {', '.join(m['variant_classes'])}" + (f" · {m['dimension']}" if m.get('dimension') else '') + (f" ({m['note']})" if m.get('note') else '') for m in pr['mechanism_evidence']) or 'none (no mechanism content)'
        comp = '<br>'.join(f"{c['name']} — {', '.join(c['evidence_kind'])} · {c['allele_resolution']} · from {c['origin']}" for c in pr.get('components', []))
        ov = '<br>'.join(f"{o['relation']} · `{o['source']}`" + (f" ({o['note']})" if o.get('note') else '') for o in pr['overlaps_with_sources'])
        L += [f"### `{s['id']}` · {s['name']}", '', '| | |', '|---|---|',
              f"| Evidence kinds | {', '.join(pr['evidence_kind'])} |", f"| Allele resolution | {pr['allele_resolution']} |", f'| Mechanism evidence | {me} |']
        if comp:
            L.append(f'| Components | {comp} |')
        L += [f"| Biological context | {pr['biological_context']} |", f"| Evidence strength | {pr['evidence_strength']} |", f"| Mechanistic confidence | {pr['mechanistic_confidence']} |",
              f"| Phenotype implication | {pr['phenotype_implication']} |", f"| Absence means | {pr['absence_means']} |"]
        if ov:
            L.append(f'| Overlaps (declared here) | {ov} |')
        L.append('')
    L += ['## Overlaps and deduplication', '',
          'Policy: identical records are combined once; independent experiments are never merged only because they concern the same gene or variant; assertions that can disagree are kept side by side, and disagreement is information. '
          'For an aggregate, reuse the aggregate and consult the original for fields the aggregate drops; count the underlying evidence once. A pair is declared once, by the source that depends on or restates the other.', '',
          '| Source | Relation | With | What is shared or contested |', '|---|---|---|---|']
    rels = [(s['id'], o) for s in sorted(sources, key=lambda x: x['id']) for o in prof[s['id']]['overlaps_with_sources']]
    for sid, o in sorted(rels, key=lambda x: (x[1]['relation'], x[0], x[1]['source'])):
        L.append(f"| `{sid}` | {o['relation']} | `{o['source']}` | {o.get('note', '')} |")
    if not rels:
        L.append('| none | | | |')
    L += ['', 'Relations: ' + '; '.join(f'**{k}** — {v}' for k, v in vocab['source_relations'].items()), '']
    aggs = [(s['id'], c) for s in sorted(sources, key=lambda x: x['id']) for c in prof[s['id']].get('components', [])]
    if aggs:
        L += ['## Components of mixed sources', '', 'Where one resource carries parts of different evidence kinds, each part keeps its own kind and resolution. The origin tells whether the part is native, taken from another catalogued source, or from an external method.', '',
              '| Source | Component | Evidence kind | Allele resolution | Origin |', '|---|---|---|---|---|']
        for sid, c in aggs:
            L.append(f"| `{sid}` | {c['name']} | {', '.join(c['evidence_kind'])} | {c['allele_resolution']} | {c['origin']} |")
        L.append('')
    open(os.path.join(HERE, 'EVIDENCE_PROFILES.md'), 'w').write('\n'.join(L))

    # ---- MECHANISM MATRIX ---------------------------------------------------------------------------------
    LET = {'direct': 'D', 'curated': 'C', 'association': 'A', 'prediction': 'P', 'population': 'F', 'context': 'X', 'case': 'K'}
    EMPIRICAL = ('direct', 'curated', 'association', 'case')
    TYPES = ['snv', 'indel', 'cnv', 'sv', 'str', 'complex_locus', 'mtdna']
    LEVELS = ['gene_level', 'region_level', 'protein_level']
    entries = [(s, m) for s in sources for m in prof[s['id']]['mechanism_evidence']]

    def entry_state(s, m):
        st = life[(s['id'], m['dimension'] if m.get('dimension') else None)][0]
        if st['usable']:
            return True, ''
        why = next((k for k in ('acquired', 'validated', 'mapped') if st[k] is not True), 'licence')
        return False, why

    def cell(axis, classes):
        es = [(s, m) for s, m in entries if m['axis'] == axis['key'] and set(classes) & set(m['variant_classes'])]
        use = collections.Counter(LET[m['basis']] for s, m in es if entry_state(s, m)[0])
        pend = sum(1 for s, m in es if not entry_state(s, m)[0])
        return (' '.join(f'{k}{v}' for k, v in sorted(use.items())) or '·') + (f' (+{pend})' if pend else '')
    legend_m = ('Cell = number of **usable** sources by basis: **D** direct measurement, **C** curated assertion, **A** statistical association, **K** case-level observation, **P** prediction, **F** population frequency or constraint, **X** context only; '
                '`(+n)` = n further sources declare the entry but are not yet usable (not acquired, validated or mapped, or without a verified licence); `·` = nothing; `n/a` = the axis does not apply to that class. '
                'Gene-, region- and protein-level entries are counted in the last column only and never fill an allele-level cell.')
    L = ['# Mechanism matrix: variant types by mechanism axis', '', gen, '',
         'Which mechanisms the catalogued evidence can inform for each type of called variant, on what basis, and what is usable now. The axes and their competing explanations are defined in `catalogue/vocab.json`; '
         'a hypothesis is never recorded as a finding, and a gene-level result is never recorded as an allele-specific one.', '', legend_m, '',
         '| Mechanism axis | ' + ' | '.join(TYPES) + ' | gene / region / protein level |', '|---|' + '---|' * (len(TYPES) + 1)]
    for a in vocab['mechanism_axes']:
        L.append(f"| {a['label']} (`{a['key']}`) | " + ' | '.join(cell(a, [t]) if t in a['variant_types'] else 'n/a' for t in TYPES) + ' | ' + (cell(a, [l for l in LEVELS if l in a['levels']]) if a['levels'] else 'n/a') + ' |')
    L += ['', '## Computed gaps', '',
          'An applicable (axis, class) cell is a gap when no usable source offers direct, curated, associational or case-level evidence for it. Context, population and prediction entries can inform but do not fill it.', '',
          '| Mechanism axis | Class | Usable but non-empirical | Declared, not yet usable (first missing stage) |', '|---|---|---|---|']
    ngap = 0
    for a in vocab['mechanism_axes']:
        for t in a['variant_types']:
            es = [(s, m) for s, m in entries if m['axis'] == a['key'] and t in m['variant_classes']]
            if any(entry_state(s, m)[0] and m['basis'] in EMPIRICAL for s, m in es):
                continue
            ngap += 1
            soft = sorted({f"`{s['id']}` ({LET[m['basis']]})" for s, m in es if entry_state(s, m)[0]})
            pend = sorted({f"`{s['id']}` ({LET[m['basis']]}; {entry_state(s, m)[1]})" for s, m in es if not entry_state(s, m)[0]})
            L.append(f"| {a['label']} | {t} | {', '.join(soft) or '—'} | {', '.join(pend) or '—'} |")
    if not ngap:
        L.append('| none | | | |')
    L += ['', f'{ngap} applicable cells have no usable empirical evidence.', '']
    for a in vocab['mechanism_axes']:
        L += [f"## {a['label']} (`{a['key']}`)", '', a['description'], '', f"- **Competing explanations:** {a['competing']}", f"- **What discriminates them:** {a['discriminating']}",
              f"- **Applies to:** {', '.join(a['variant_types'])}" + (f"; levels: {', '.join(a['levels'])}" if a['levels'] else ''), '']
        es = [(s, m) for s, m in entries if m['axis'] == a['key']]
        if not es:
            L += ['*No catalogued source declares evidence on this axis.*', '']
            continue
        L += ['| Source | Basis | Classes | Dimension | Usable now | Note |', '|---|---|---|---|---|---|']
        for s, m in sorted(es, key=lambda x: (x[1]['basis'], x[0]['id'])):
            ok, why = entry_state(s, m)
            L.append(f"| `{s['id']}` | {m['basis']} | {', '.join(m['variant_classes'])} | {m.get('dimension') or 'all'} | {'✓' if ok else 'no (' + why + ')'} | {m.get('note', '')} |")
        L.append('')
    none_ = [s['id'] for s in sources if not prof[s['id']]['mechanism_evidence']]
    L += ['## Sources with no mechanism content', '', ', '.join(f'`{x}`' for x in none_) or 'none', '']
    open(os.path.join(HERE, 'MECHANISM_MATRIX.md'), 'w').write('\n'.join(L))
    summary = {k: sum(1 for s in sources if life[(s['id'], None)][0][k] is True) for k in STAGES}
    print('sources reaching each stage (all files):', summary)
    print('sources usable for at least one dimension:', sum(1 for s in sources if usable_dims[s['id']]))
    return 0


if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('--check', action='store_true')
    sys.exit(build(ap.parse_args().check))
