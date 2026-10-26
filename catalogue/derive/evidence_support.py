"""Find statements in the source records that none of their stored evidence supports.

  python3 catalogue/derive/evidence_support.py [--source ID] [--limit N]

For every curated text item of a source record (scope, mapping_methods, quantitative, licence, version, coverage_limitations,
known_omissions, status.notes, and the free-text fields of evidence_profile) the script extracts distinctive tokens: numbers of three or more characters, and capitalised or
mixed-case names. A token counts as supported if it occurs in the record's own evidence excerpts, in the machine-generated
reports (validation, compatibility, manifest, target definitions) or in the record's own name, URL and maintainer. Tokens that occur
nowhere are listed so that a person can decide whether the statement came from the evidence under other words, or from memory
and must be reworded or listed in unverified_fields.

This is a heuristic filter for review, not a proof: a supported token does not make a sentence true, and an unsupported token
is a prompt to check, not a verdict.
"""
import argparse
import glob
import json
import os
import re

HERE = os.path.dirname(os.path.abspath(__file__))
CAT = os.path.abspath(os.path.join(HERE, '..'))
ROOT = os.path.abspath(os.path.join(CAT, '..'))

GENERIC = set('''The This That These Those With Without From Into Over Under Only Also Both Each Some Most Many Several Other Others Per Not No Yes
And But For Are Was Were Has Have Had Its Their There Than Then When Where Which While Will Would Can May Might Must Should Does Did
Open Data Database Databases Human Humans Gene Genes Genome Genomic Genomes Variant Variants Sequence Sequences Protein Proteins Tool Tools Table Tables File Files
Release Releases Version Versions Record Records Source Sources Study Studies Sample Samples Score Scores Level Levels Type Types Class Classes Method Methods
Access Terms Licence License Policy Use Used Using Users Page Pages Site Sites Server Download Downloads Listing Directory Index API REST JSON TSV CSV VCF BED GTF GFF FASTA
GRCh37 GRCh38 Mouse Cancer Cell Cells Tissue Tissues Disease Diseases Clinical Population Populations Frequency Frequencies Evidence Curated Predicted Prediction Predictions
Short Long Large Small High Low Full Whole Public Free Note Notes See Such Same Different New Old Recent Current Latest First Second Third Not None All Any'''.split())

TEXT_FIELDS = ('scope', 'mapping_methods')


_INDEX = None


def evidence_index():
    global _INDEX
    if _INDEX is None:
        _INDEX = {}
        for line in open(os.path.join(CAT, 'evidence', 'evidence_index.jsonl')):
            r = json.loads(line)
            _INDEX.setdefault(f"{r['source']}/{r['name']}", []).append(r)
    return _INDEX


def evidence_text(src_rec):
    out = []
    for e in src_rec.get('evidence', []):
        s, _, n = e.partition('/')
        for ext in ('.txt', '.json'):
            p = os.path.join(CAT, 'evidence', s, n + ext)
            if os.path.exists(p):
                out.append(open(p, errors='replace').read())
                break
        for r in evidence_index().get(e, []):
            out.append(' '.join(str(r.get(k) or '') for k in ('url', 'note')))
    return '\n'.join(out)


def manifest_bytes(src_rec):
    """Byte counts of the source's acquired files and their sum, so that totals quoted in the record can be matched."""
    sizes = []
    mp = os.path.join(ROOT, 'acquisition', 'manifest.jsonl')
    if os.path.exists(mp):
        for line in open(mp):
            r = json.loads(line)
            if r['id'] in src_rec.get('acquisition_targets', []) and r.get('ok') and isinstance(r.get('bytes'), int):
                sizes.append(r['bytes'])
    return sizes + ([sum(sizes)] if sizes else [])


def report_text():
    t = []
    for rel in ('acquisition/validation_report.json', 'acquisition/compat_report.json', 'acquisition/manifest.jsonl', 'acquisition/targets.json', 'inventory/EXISTING_ASSETS.json'):
        p = os.path.join(ROOT, rel)
        if os.path.exists(p):
            t.append(open(p, errors='replace').read())
    for p in sorted(glob.glob(os.path.join(ROOT, 'scrape', '*', 'README.md'))):
        t.append(open(p, errors='replace').read())
    return '\n'.join(t)


def tracked_names():
    import subprocess
    try:
        out = subprocess.run(['git', '-C', ROOT, 'ls-files'], capture_output=True, text=True, check=True).stdout.split('\n')
    except Exception:
        return set()
    return {os.path.basename(x).lower() for x in out if x}


def items(rec):
    for f in TEXT_FIELDS:
        if isinstance(rec.get(f), str):
            yield f, rec[f]
    for f in ('quantitative',):
        for k, v in rec.get(f, {}).items():
            if isinstance(v, str):
                yield f'{f}.{k}', v
    for f in ('licence', 'version', 'status'):
        for k, v in rec.get(f, {}).items():
            if isinstance(v, str) and k not in ('verified', 'accessible'):
                yield f'{f}.{k}', v
    for f in ('coverage_limitations', 'known_omissions'):
        for i, v in enumerate(rec.get(f, [])):
            yield f'{f}[{i}]', v
    pr = rec.get('evidence_profile', {})
    for f in ('biological_context', 'evidence_strength', 'mechanistic_confidence', 'phenotype_implication', 'absence_means'):
        if isinstance(pr.get(f), str):
            yield f'profile.{f}', pr[f]
    for i, m in enumerate(pr.get('mechanism_evidence', [])):
        if m.get('note'):
            yield f"profile.mechanism_evidence[{i}].note", m['note']
    for i, c in enumerate(pr.get('components', [])):
        yield f'profile.components[{i}].name', c['name']
    for i, o in enumerate(pr.get('overlaps_with_sources', [])):
        if o.get('note'):
            yield f"profile.overlaps[{i}].note", o['note']


UNITS = {'KB': 1e3, 'KIB': 1024, 'MB': 1e6, 'MIB': 1024 ** 2, 'GB': 1e9, 'GIB': 1024 ** 3, 'TB': 1e12, 'TIB': 1024 ** 4}
SIZE = re.compile(r'(\d[\d,]*\.?\d*)\s?(KiB|MiB|GiB|TiB|KB|MB|GB|TB)\b')


def tokens(text):
    sizes = [(float(m.group(1).replace(',', '')), m.group(2).upper()) for m in SIZE.finditer(text)]
    text_nosize = SIZE.sub(' ', text)
    ks = [float(m.group(1)) * 1000 for m in re.finditer(r'\b(\d+(?:\.\d+)?)k\b', text_nosize)]
    text_nosize = re.sub(r'\b\d+(?:\.\d+)?k\b', ' ', text_nosize)
    nums = set(re.findall(r'\d[\d,\.]*\d|\d{3,}', text_nosize))
    names = set()
    for m in re.finditer(r"\b([A-Za-z][A-Za-z0-9\+_\.]{2,})\b", text_nosize):
        w = m.group(1).strip('.')
        if w in GENERIC or len(w) < 3:
            continue
        mixed = any(c.isupper() for c in w[1:]) or any(c.isdigit() for c in w)
        capital = w[0].isupper() and m.start() > 0 and text_nosize[max(0, m.start() - 2):m.start()] not in ('. ', ': ', '; ')
        if mixed or capital:
            names.add(w)
    return nums, names, sizes, ks


def byte_values(corpus, extra=()):
    vals = {int(x) for x in re.findall(r'\b\d{5,}\b', corpus)}
    for num, u in re.findall(r'\b(\d+(?:\.\d+)?)\s?([KMGT])\b', corpus):
        vals.add(int(float(num) * 1024 ** ('KMGT'.index(u) + 1)))
    vals.update(int(x) for x in extra)
    return sorted(vals)


def thousand_supported(corpus_ints, n, tol=0.03):
    return any(abs(v - n) <= tol * n for v in corpus_ints)


def size_supported(vals, num, unit, tol=0.012):
    import bisect
    target = num * UNITS[unit]
    lo = bisect.bisect_left(vals, target * (1 - tol))
    return lo < len(vals) and vals[lo] <= target * (1 + tol)


def norm(s):
    return re.sub(r'[\s,]+', '', s.lower())


def declared_unverified(rec):
    """Field names that the record itself lists as unverified (the whole field is then exempt from this check and reported separately)."""
    words = ' '.join(rec.get('unverified_fields', [])).lower()
    return {f for f in ('scope', 'mapping_methods', 'quantitative', 'licence', 'version', 'coverage_limitations', 'known_omissions') if re.search(r'\b' + f.replace('_', '[ _]') + r'\b', words)}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--source')
    ap.add_argument('--limit', type=int, default=8)
    ap.add_argument('--json', help='write the flagged items to this file')
    a = ap.parse_args()
    reports = report_text()
    names_in_repo = tracked_names()
    total = flagged_items = exempt_items = 0
    out = []
    for p in sorted(glob.glob(os.path.join(CAT, 'sources', '*.json'))):
        rec = json.load(open(p))
        if a.source and rec['id'] != a.source:
            continue
        corpus = evidence_text(rec) + '\n' + reports + '\n' + ' '.join([rec['name'], rec['url'], rec['maintainer'], ' '.join(rec['provenance'])])
        low = corpus.lower()
        nlow = norm(corpus)
        vals = byte_values(corpus, manifest_bytes(rec))
        cints = sorted({int(x.replace(',', '')) for x in re.findall(r'\b\d[\d,]{2,}\b', corpus) if x.replace(',', '').isdigit()})
        exempt = declared_unverified(rec)
        lines = []
        for field, text in items(rec):
            total += 1
            nums, names, sizes, ks = tokens(text)
            miss = [n for n in sorted(names) if n.lower() not in low and n.lower().rstrip('s') not in low and n.lower() not in names_in_repo]
            miss += [n for n in sorted(nums) if norm(n) not in nlow]
            for v, u in sizes:
                lit = [f'{v:g}{u.lower()}', f'{v:g}{u[0].lower()}', f'{v:g}{u[0].lower()}ib', f'{v:g}{u[0].lower()}b']
                if not any(l in nlow for l in lit) and not size_supported(vals, v, u):
                    miss.append(f'{v:g} {u}')
            miss += [f'{k / 1000:g}k' for k in ks if not thousand_supported(cints, int(k))]
            if miss:
                if field.split('.')[0].split('[')[0] in exempt or re.search(r'unverified|not retrieved|not read|not examined', text):
                    exempt_items += 1
                    continue
                flagged_items += 1
                lines.append(f'  {field}: {text[:140]!r}\n      unsupported: {miss[:a.limit]}')
                out.append(dict(source=rec['id'], field=field, tokens=miss))
        if lines:
            print(f"== {rec['id']} ({len(lines)} flagged)")
            print('\n'.join(lines))
    print(f'\n{flagged_items} of {total} text items contain tokens found in no stored evidence; {exempt_items} further items sit in fields the record lists as unverified')
    if a.json:
        json.dump(out, open(a.json, 'w'), indent=1)


if __name__ == '__main__':
    main()
