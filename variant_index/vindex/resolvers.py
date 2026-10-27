"""Raw-record resolvers: the exact text an adapter hashed for an evidence row, recovered again from the acquired file by the row's locator.

A resolver is `fn(base, locators) -> {locator: text}`. `base` is the acquired file or directory of the release's target. A locator that cannot be resolved is left out,
and the verifier then reports the row as unresolved, never as verified. Each function restates the rule of one adapter (what the record is, how the locator counts), so a change
to the adapter that is not mirrored here makes the verification fail instead of passing silently.

Records that an adapter builds by aggregating several source rows (IMPC mouse-phenotype summaries) have no single source record to re-open and have no resolver.
"""
import csv
import json
import os
import re
import zipfile

import pyarrow.compute as pc
import pyarrow.dataset as pds
import pyarrow.parquet as pq

LINE = re.compile(r'^L(\d+)$')
MEMBER_LINE = re.compile(r'^(?P<member>.+)#L(?P<n>\d+)$')
OT_ROW = re.compile(r'^(?P<dataset>[a-z_0-9]+)/(?P<file>[^#]+)#(?P<i>\d+)$')


def clingen_csv_records(base, locators):
    """ClinGen gene-disease validity: the raw record is the JSON of the csv.reader row, and L<n> counts csv records from the top of the file (comment and separator rows included)."""
    wanted = {int(LINE.match(l).group(1)): l for l in locators if LINE.match(l)}
    out = {}
    if not wanted:
        return out
    with open(base, encoding='utf-8', newline='') as f:
        for n, row in enumerate(csv.reader(f), start=1):
            if n in wanted:
                out[wanted[n]] = json.dumps(row, ensure_ascii=False)
            if n >= max(wanted):
                break
    return out


def json_array(prefix, key='id', sort_keys=True):
    """A JSON array of objects (CPIC API tables); the locator is '<prefix>:<id>' and the raw record is the object as the adapter dumped it."""
    def resolve(base, locators):
        with open(base, encoding='utf-8') as f:
            index = {f'{prefix}:{r[key]}': r for r in json.load(f)}
        return {l: json.dumps(index[l], sort_keys=sort_keys, ensure_ascii=False) for l in locators if l in index}
    resolve.__name__ = f'json_array_{prefix}'
    return resolve


def strchive_stranger(base, locators):
    """STRchive-disease-loci.hg38.stranger.json#<LocusId>: the locus object."""
    with open(base, encoding='utf-8') as f:
        index = {r['LocusId']: r for r in json.load(f)}
    out = {}
    for l in locators:
        _, _, lid = l.partition('#')
        if lid in index:
            out[l] = json.dumps(index[lid], sort_keys=True, ensure_ascii=False)
    return out


def zip_members(base, locators):
    """ClinPGx zip exports. '<member>#L<n>' is line n of a TSV member (lines split on \\n, so a trailing \\r belongs to the record); a bare member name is a JSON member whose whole text
    is the record. Text is decoded as the adapter decoded it."""
    out, cache = {}, {}
    with zipfile.ZipFile(base) as z:
        names = set(z.namelist())
        for l in locators:
            m = MEMBER_LINE.match(l)
            if m and m.group('member') in names:
                member = m.group('member')
                if member not in cache:
                    cache[member] = z.read(member).decode('utf-8', errors='replace').split('\n')
                n = int(m.group('n'))
                if 0 < n <= len(cache[member]):
                    out[l] = cache[member][n - 1]
            elif l in names:
                out[l] = z.read(l).decode('utf-8', errors='replace')
    return out


def ot_files(base):
    return {os.path.basename(p): p for p in sorted(os.path.join(base, f) for f in os.listdir(base) if f.endswith('.parquet'))}


def ot_row_fields(fields):
    """Open Targets evidence tables: '<dataset>/<file>#<row index>' and the raw record is the JSON of selected fields of that row."""
    def resolve(base, locators):
        files, out, by_file = ot_files(base), {}, {}
        for l in locators:
            m = OT_ROW.match(l)
            if m and m.group('file') in files:
                by_file.setdefault(m.group('file'), {})[int(m.group('i'))] = l
        for fname, rows in by_file.items():
            t = pq.read_table(files[fname])
            for i, l in rows.items():
                if i < t.num_rows:
                    r = t.slice(i, 1).to_pylist()[0]
                    out[l] = json.dumps({k: r[k] for k in fields}, sort_keys=True, default=str)
        return out
    resolve.__name__ = 'ot_row_fields'
    return resolve


def ot_target_essentiality(base, locators):
    """DepMap essentiality per target: the raw record is derived from the row (its screens counted and the fraction with gene effect <= -1), so it is recomputed from the row."""
    files, out, by_file = ot_files(base), {}, {}
    for l in locators:
        m = OT_ROW.match(l)
        if m and m.group('file') in files:
            by_file.setdefault(m.group('file'), {})[int(m.group('i'))] = l
    for fname, rows in by_file.items():
        t = pq.read_table(files[fname])
        for i, l in rows.items():
            if i >= t.num_rows:
                continue
            r = t.slice(i, 1).to_pylist()[0]
            allv = [s['geneEffect'] for ts in (r['depMapEssentiality'] or ()) for s in (ts['screens'] or ()) if s.get('geneEffect') is not None]
            frac = (sum(1 for e in allv if e <= -1.0) / len(allv)) if allv else None
            out[l] = json.dumps(dict(targetId=r['targetId'], isEssential=r['isEssential'], n_screens=len(allv), frac=frac), sort_keys=True)
    return out


def ot_credible_set(base, locators):
    """credible_set/<file>#<studyLocusId>:<variantId>: the raw record is the study locus id and the member of its locus (all member fields)."""
    files, out, by_file = ot_files(base), {}, {}
    for l in locators:
        m = re.match(r'^credible_set/(?P<file>[^#]+)#(?P<sl>[^:]+):(?P<vid>.+)$', l)
        if m and m.group('file') in files:
            by_file.setdefault(m.group('file'), {})[(m.group('sl'), m.group('vid'))] = l
    for fname, wanted in by_file.items():
        ids = sorted({sl for sl, _ in wanted})
        t = pds.dataset(files[fname], format='parquet').to_table(columns=['studyLocusId', 'locus'], filter=pc.field('studyLocusId').isin(ids))
        for row in t.to_pylist():
            for member in row['locus'] or ():
                l = wanted.get((row['studyLocusId'], member['variantId']))
                if l is not None:
                    out[l] = json.dumps(dict(studyLocusId=row['studyLocusId'], member=member), sort_keys=True, default=str)
    return out


def ot_l2g(base, locators):
    """l2g_prediction#<studyLocusId>:<geneId>: the raw record is the prediction row (study locus, gene, score)."""
    files = ot_files(base)
    wanted = {}
    for l in locators:
        m = re.match(r'^l2g_prediction#(?P<sl>[^:]+):(?P<gene>.+)$', l)
        if m:
            wanted[(m.group('sl'), m.group('gene'))] = l
    out = {}
    if not wanted:
        return out
    ids = sorted({sl for sl, _ in wanted})
    for p in files.values():
        t = pds.dataset(p, format='parquet').to_table(columns=['studyLocusId', 'geneId', 'score'], filter=pc.field('studyLocusId').isin(ids))
        for sl, g, sc in zip(t['studyLocusId'].to_pylist(), t['geneId'].to_pylist(), t['score'].to_pylist()):
            l = wanted.get((sl, g))
            if l is not None:
                out[l] = json.dumps(dict(studyLocusId=sl, geneId=g, score=sc), sort_keys=True)
    return out


def mavedb_native(base, locators, annotations_dir=None):
    """<stem>.annotations.csv#L<n>: the exact text of the CSV record that starts on line n of the annotations file of a score set (a record may span lines)."""
    from .build.mavedb import ANN, csv_records
    root = annotations_dir or ANN
    out, by_file = {}, {}
    for l in locators:
        m = MEMBER_LINE.match(l)
        if m:
            by_file.setdefault(m.group('member'), {})[int(m.group('n'))] = l
    for name, rows in by_file.items():
        path = os.path.join(root, name)
        if not os.path.exists(path):
            continue
        with open(path, encoding='utf-8', newline='') as f:
            for ln, _row, raw, _hdr in csv_records(f):
                if ln in rows:
                    out[rows[ln]] = raw
    return out


def ot_impc_summary(base, locators):
    """evidence_impc#["target", "allelic composition", "genetic background"]: the raw record is the summary of all IMPC rows with that key (how many, and which diseases they link),
    recomputed from the rows of every file of the dataset."""
    wanted = {}
    for l in locators:
        if l.startswith('evidence_impc#'):
            try:
                gene, comp, bg = json.loads(l.partition('#')[2])
            except (ValueError, TypeError):
                continue
            wanted[(gene, comp, bg)] = l
    out = {}
    if not wanted:
        return out
    genes = sorted({k[0] for k in wanted})
    agg = {k: [0, set()] for k in wanted}
    for p in ot_files(base).values():
        t = pds.dataset(p, format='parquet').to_table(columns=['targetId', 'biologicalModelAllelicComposition', 'biologicalModelGeneticBackground', 'diseaseId'],
                                                      filter=pc.field('targetId').isin(genes))
        for g, c, b, d in zip(t['targetId'].to_pylist(), t['biologicalModelAllelicComposition'].to_pylist(), t['biologicalModelGeneticBackground'].to_pylist(), t['diseaseId'].to_pylist()):
            a = agg.get((g, c, b))
            if a is not None:
                a[0] += 1
                if d:
                    a[1].add(d)
    for k, l in wanted.items():
        rows, diseases = agg[k]
        if rows:
            out[l] = json.dumps(dict(targetId=k[0], composition=k[1], background=k[2], rows=rows, diseases=sorted(diseases)), sort_keys=True)
    return out
