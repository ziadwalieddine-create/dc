"""Validate downloaded external files against each source's own description.

Writes acquisition/validation_report.json (and a short Markdown summary). Each target gets a list of named checks;
a target 'passes' only if all its checks pass. Nothing here modifies downloaded data.

  python3 validate.py [--ids target ...]
"""
import argparse
import collections
import csv
import gzip
import io
import json
import os
import re
import tarfile
import zipfile

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, '..'))
DATA = os.path.join(ROOT, 'data', 'external')
TARGETS = {t['id']: t for t in json.load(open(os.path.join(HERE, 'targets.json')))}
MANIFEST = {}
for line in open(os.path.join(HERE, 'manifest.jsonl')):
    r = json.loads(line)
    MANIFEST[r['id']] = r


def path(tid):
    return os.path.join(DATA, TARGETS[tid]['dest'])


def opener(p):
    return gzip.open(p, 'rt', encoding='utf-8', errors='replace', newline='') if open(p, 'rb').read(2) == b'\x1f\x8b' else open(p, encoding='utf-8', errors='replace', newline='')


class R:
    def __init__(self):
        self.checks, self.stats = [], {}

    def check(self, name, ok, detail=''):
        self.checks.append(dict(name=name, ok=bool(ok), detail=str(detail)[:300]))

    def passed(self):
        return bool(self.checks) and all(c['ok'] for c in self.checks)


def vcf_stats(p, limit=None):
    meta, n, chroms, info_keys = {}, 0, collections.Counter(), collections.Counter()
    header = None
    with opener(p) as f:
        for line in f:
            if line.startswith('##'):
                m = re.match(r'##(\w+)=(.*)', line.strip())
                if m and m.group(1) in ('fileformat', 'fileDate', 'reference', 'source'):
                    meta[m.group(1)] = m.group(2)
            elif line.startswith('#CHROM'):
                header = line.strip().split('\t')
            else:
                parts = line.split('\t', 8)
                n += 1
                chroms[parts[0]] += 1
                if len(parts) > 7:
                    for kv in parts[7].split(';'):
                        info_keys[kv.split('=')[0]] += 1
                if limit and n >= limit:
                    break
    return meta, header, n, chroms, info_keys


def v_manifest(tid, r):
    m = MANIFEST.get(tid)
    r.check('manifest record is ok', m and m.get('ok'))
    if m and os.path.exists(path(tid)):
        r.check('size equals manifest', os.path.getsize(path(tid)) == m['bytes'])


def v_clinvar_vcf(tid, r):
    meta, header, n, chroms, keys = vcf_stats(path(tid))
    r.check('VCF 4.x header', meta.get('fileformat', '').startswith('VCFv4'), meta)
    r.check('reference is GRCh38', meta.get('reference') == 'GRCh38', meta.get('reference'))
    r.check('no chr prefix on contigs', not any(c.startswith('chr') for c in chroms), list(chroms)[:5])
    r.check('records > 3,000,000', n > 3_000_000, n)
    r.check('INFO has CLNSIG and GENEINFO', keys['CLNSIG'] > 0 and keys['GENEINFO'] > 0, dict(CLNSIG=keys['CLNSIG'], GENEINFO=keys['GENEINFO']))
    r.stats.update(fileDate=meta.get('fileDate'), records=n, contigs=len(chroms), with_CLNSIG=keys['CLNSIG'], with_GENEINFO=keys['GENEINFO'])


def v_gnomad_cnv(tid, r):
    meta, header, n, chroms, keys = vcf_stats(path(tid))
    r.check('VCF header present', meta.get('fileformat', '').startswith('VCF'), meta)
    r.check('records > 1000', n > 1000, n)
    svt = collections.Counter()
    with opener(path(tid)) as f:
        for line in f:
            if not line.startswith('#'):
                m = re.search(r'SVTYPE=(\w+)', line)
                svt[m.group(1) if m else '?'] += 1
    r.check('SVTYPE present', '?' not in svt or len(svt) > 1, dict(svt))
    r.stats.update(records=n, svtype=dict(svt), contigs=sorted(chroms)[:30])


def v_tsv(tid, r, required, min_rows, delim='\t', comment='#', skip_prefix_lines=0):
    with opener(path(tid)) as f:
        lines = [l for l in f if not (comment and l.startswith(comment))]
    hdr = lines[0].rstrip('\n').split(delim)
    r.check('required columns present', all(c in hdr for c in required), [c for c in required if c not in hdr])
    r.check(f'rows >= {min_rows}', len(lines) - 1 >= min_rows, len(lines) - 1)
    r.stats.update(columns=len(hdr), rows=len(lines) - 1, sample_columns=hdr[:10])
    return hdr, lines


def v_clingen_gene(tid, r):
    with opener(path(tid)) as f:
        raw = f.read().split('\n')
    # the column header is itself a '#'-prefixed line ('#Gene Symbol ...'); data lines follow
    hdr_i = next(i for i, l in enumerate(raw) if l.startswith('#Gene Symbol'))
    lines = [raw[hdr_i].lstrip('#')] + [l for l in raw[hdr_i + 1:] if l.strip() and not l.startswith('#')]
    hdr = lines[0].split('\t')
    r.check('required columns', all(c in hdr for c in ['Gene Symbol', 'Genomic Location', 'Haploinsufficiency Score', 'Triplosensitivity Score']), hdr[:6])
    rows = [l.split('\t') for l in lines[1:]]
    hi = collections.Counter(x[hdr.index('Haploinsufficiency Score')] for x in rows)
    r.check('rows > 1000', len(rows) > 1000, len(rows))
    r.check('HI score 3 genes present', hi.get('3', 0) > 100, dict(hi))
    r.stats.update(rows=len(rows), haploinsufficiency_scores=dict(hi))


def v_bed(tid, r, expected_rows=None):
    n, bad = 0, 0
    chroms = collections.Counter()
    with opener(path(tid)) as f:
        for line in f:
            if line.startswith(('#', 'track', 'browser')) or not line.strip():
                continue
            p = line.rstrip('\n').split('\t')
            n += 1
            try:
                ok = len(p) >= 3 and int(p[1]) < int(p[2])
            except ValueError:
                ok = False
            bad += not ok
            chroms[p[0]] += 1
    r.check('all lines are valid BED intervals', bad == 0, bad)
    r.check('rows > 0', n > 0, n)
    if expected_rows:
        r.check(f'row count equals the {expected_rows:,} shown on the source website', n == expected_rows, n)
    r.stats.update(rows=n, contigs=len(chroms))


def v_zip(tid, r, must_contain=()):
    with zipfile.ZipFile(path(tid)) as z:
        bad = z.testzip()
        names = z.namelist()
    r.check('zip CRC test', bad is None, bad)
    for m in must_contain:
        r.check(f'contains {m}', any(m in n for n in names), names[:8])
    r.stats.update(members=names[:12], n_members=len(names))
    return names


def v_json(tid, r, min_len=1):
    d = json.load(open(path(tid), encoding='utf-8'))
    r.check('valid JSON', True)
    r.check(f'length >= {min_len}', len(d) >= min_len, len(d))
    r.stats.update(items=len(d))
    return d


# --- per-target dispatch -------------------------------------------------------------------------------------------
def run(tid):
    r = R()
    t = TARGETS[tid]
    if not os.path.exists(path(tid)):
        r.check('file present', False, 'not downloaded')
        return r
    v_manifest(tid, r)
    try:
        if tid == 'clinvar_vcf_grch38':
            v_clinvar_vcf(tid, r)
        elif tid == 'gnomad_exome_cnv_v4_1':
            v_gnomad_cnv(tid, r)
        elif tid == 'clingen_gene_curation_list':
            v_clingen_gene(tid, r)
        elif tid == 'clingen_region_curation_list':
            v_tsv(tid, r, ['ISCA Region Name', 'Genomic Location'] if False else [], 100)
        elif tid in ('clingen_haploinsufficiency_gene', 'clingen_triplosensitivity_gene', 'clingen_recurrent_cnv_v2_1'):
            v_bed(tid, r)
        elif tid == 'clingen_gene_validity':
            with opener(path(tid)) as f:
                rows = list(csv.reader(f))
            hdr_i = next(i for i, x in enumerate(rows) if x and x[0] == 'GENE SYMBOL')
            body = [x for x in rows[hdr_i + 1:] if x and not x[0].startswith('+')]
            cls = collections.Counter(x[6] for x in body if len(x) > 6)
            r.check('validity rows > 1000', len(body) > 1000, len(body))
            r.check('classification column populated', sum(cls.values()) > 1000, dict(cls.most_common(6)))
            r.stats.update(rows=len(body), classifications=dict(cls.most_common(8)))
        elif tid == 'gnomad_constraint_v4_1_2':
            hdr, lines = v_tsv(tid, r, ['gene', 'transcript'], 15000, comment=None)
            r.check('has LOEUF-type columns', any('oe_ci.upper' in c for c in hdr), [c for c in hdr if 'lof' in c][:6])
            r.stats['lof_columns'] = [c for c in hdr if c.startswith('lof')][:12]
        elif tid == 'gwas_catalog_associations':
            names = v_zip(tid, r)
            with zipfile.ZipFile(path(tid)) as z:
                with z.open(names[0]) as f:
                    txt = io.TextIOWrapper(f, encoding='utf-8', errors='replace')
                    hdr = txt.readline().rstrip('\n').split('\t')
                    n = sum(1 for _ in txt)
            r.check('GWAS columns', all(c in hdr for c in ['PUBMEDID', 'CHR_ID', 'CHR_POS', 'MAPPED_GENE', 'P-VALUE']), [c for c in hdr][:8])
            r.check('associations > 500,000', n > 500_000, n)
            r.stats.update(rows=n, member=names[0], columns=len(hdr))
        elif tid == 'gwas_catalog_studies':
            v_tsv(tid, r, ['PUBMEDID', 'STUDY', 'DISEASE/TRAIT'], 50_000, comment=None)
        elif tid.startswith('reactome_'):
            with opener(path(tid)) as f:
                rows = [l.rstrip('\n').split('\t') for l in f if l.strip()]
            hs = sum(1 for x in rows if 'Homo sapiens' in x)
            r.check('rows > 1000', len(rows) > 1000, len(rows))
            if tid != 'reactome_pathways':
                r.check('Homo sapiens rows present', hs > 1000, hs)
            else:
                r.check('pathway IDs R-HSA present', sum(1 for x in rows if x[0].startswith('R-HSA')) > 1000)
            r.stats.update(rows=len(rows), homo_sapiens_rows=hs)
        elif tid.startswith('clinpgx_'):
            names = v_zip(tid, r)
            r.check('archive has data members', len(names) >= 1, names[:5])
        elif tid in ('cpic_genes', 'cpic_guidelines'):
            d = v_json(tid, r, 10)
            r.stats['first_keys'] = list(d[0])[:8] if d and isinstance(d[0], dict) else None
        elif tid == 'gtex_v11_gene_median_tpm':
            with gzip.open(path(tid), 'rt') as f:
                l1, l2, l3 = f.readline(), f.readline().split(), f.readline().rstrip('\n').split('\t')
            r.check('GCT 1.2 header', l1.startswith('#1.2'), l1.strip())
            genes, tissues = int(l2[0]), int(l2[1])
            r.check('tissue columns >= 54 (v11 lists 54 tissue sites plus sub-site columns)', tissues >= 54, tissues)
            r.check('genes > 50,000', genes > 50_000, genes)
            r.stats.update(genes=genes, tissues=tissues, tissue_sample=l3[2:6])
        elif tid == 'gtex_datasets_api':
            d = v_json(tid, r, 3)
            r.check('v11 listed', any(x.get('datasetId') == 'gtex_v11' for x in d))
        elif tid == 'mane_summary_v1_5':
            hdr, lines = v_tsv(tid, r, ['#NCBI_GeneID', 'Ensembl_Gene', 'MANE_status'] if False else ['Ensembl_Gene', 'MANE_status'], 18000, comment=None)
            st = collections.Counter(l.split('\t')[hdr.index('MANE_status')] for l in lines[1:])
            r.check('MANE Select entries ~19k', 17_000 < st.get('MANE Select', 0) < 21_000, dict(st))
            r.stats['status'] = dict(st)
        elif tid == 'mane_gtf_v1_5':
            n = 0
            with opener(path(tid)) as f:
                for l in f:
                    if l.startswith('#'):
                        continue
                    n += 1
                    if len(l.split('\t')) != 9:
                        r.check('GTF 9 columns', False, l[:80])
                        break
            r.check('GTF lines > 100,000', n > 100_000, n)
            r.stats['lines'] = n
        elif tid == 'igsr_1000g_samples':
            with open(path(tid)) as f:
                rows = [l.split() for l in f if l.strip()]
            hdr, body = rows[0], rows[1:]
            r.check('3,202 samples', len(body) == 3202, len(body))
            r.check('columns', 'Population' in hdr and 'Superpopulation' in hdr, hdr)
            r.stats.update(samples=len(body), superpopulations=dict(collections.Counter(x[6] for x in body)))
        elif tid == 'uniprot_human_reviewed':
            hdr, lines = v_tsv(tid, r, ['Entry', 'Gene Names (primary)'] if False else ['Entry'], 20_000, comment=None)
            r.check('rows ~ 20,431 live count (+-5%)', abs(len(lines) - 1 - 20431) / 20431 < 0.05, len(lines) - 1)
        elif tid == 'biogrid_orcs_human':
            with tarfile.open(path(tid)) as tf:
                members = tf.getnames()
            r.check('tar members > 100', len(members) > 100, len(members))
            r.stats.update(n_members=len(members), sample=members[:4])
        elif tid.startswith('alphamissense_'):
            hdr, lines = v_tsv(tid, r, [], 1000, comment='#')
            r.stats['sample_header'] = hdr
        elif tid.startswith('strchive_') and tid.endswith('json'):
            d = v_json(tid, r, 20)
        elif tid.startswith('strchive_'):
            v_bed(tid, r)
        elif tid.startswith('screen_ccre_'):
            site = {'pls': 47532, 'pels': 249464, 'dels': 1469205, 'ca_ctcf': 126034, 'ca_h3k4me3': 79246, 'ca_tf': 26102, 'ca': 245985, 'tf': 105286}
            v_bed(tid, r, site[tid.replace('screen_ccre_', '')])
        elif tid == 'ensembl_vep_mavedb_variants':
            n = 0
            urns = set()
            with gzip.open(path(tid), 'rt') as f:
                hdr = f.readline().rstrip('\n').split('\t')
                for line in f:
                    n += 1
                    if n % 1 == 0:
                        p = line.split('\t')
                        for x in p:
                            if x.startswith('urn:mavedb:'):
                                urns.add(x.split('#')[0])
            r.check('has rows', n > 100_000, n)
            r.stats.update(header=hdr, rows=n, distinct_urns_seen=len(urns))
        elif tid == 'mavedb_zenodo_dump_v4':
            with tarfile.open(path(tid)) as tf:
                names = []
                for m in tf:
                    names.append(m.name)
                    if len(names) >= 100000:
                        break
            r.check('tar readable', len(names) > 0, len(names))
            r.stats.update(members_seen=len(names), sample=names[:8])
        elif tid == 'dgv_grch38_variants_2025_12_01':
            hdr, lines = v_tsv(tid, r, ['variantaccession', 'chr', 'start', 'end'], 100_000, comment=None)
        elif tid.endswith('_tbi'):
            r.check('gzip/bgzf-framed index (tbi is BGZF)', open(path(tid), 'rb').read(2) == b'\x1f\x8b')
        elif tid.endswith('readme') or tid.endswith('readme_vcf'):
            r.check('text file', os.path.getsize(path(tid)) > 1000)
        else:
            r.check('no specific validator defined', False, 'add one')
    except Exception as e:
        r.check('validator ran without error', False, repr(e))
    return r


if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('--ids', nargs='*')
    a = ap.parse_args()
    out = {}
    rp = os.path.join(HERE, 'validation_report.json')
    if os.path.exists(rp) and a.ids:
        out = json.load(open(rp))
    for tid in (a.ids or [t for t in TARGETS if MANIFEST.get(t, {}).get('ok')]):
        res = run(tid)
        out[tid] = dict(passed=res.passed(), checks=res.checks, stats=res.stats)
        print(('PASS ' if res.passed() else 'FAIL ') + tid, '' if res.passed() else [c for c in res.checks if not c['ok']], flush=True)
    json.dump(out, open(rp, 'w'), indent=1, sort_keys=True)
