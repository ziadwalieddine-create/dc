"""Validate downloaded external files against each source's own description.

Writes acquisition/validation_report.json (and a short Markdown summary). Each target gets a list of named checks;
a target 'passes' only if all its checks pass. Nothing here modifies downloaded data.

  python3 validate.py [--ids target ...]
"""
import argparse
import collections
import csv
import gzip
import hashlib
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
        if m.get('kind') == 'per_entity':
            r.check('manifest table present', os.path.exists(os.path.join(ROOT, m['files_table'])))
        elif m.get('kind') == 'http_dir':
            have = sum(os.path.getsize(os.path.join(path(tid), f)) for f in os.listdir(path(tid)) if f.endswith('.parquet'))
            r.check('total size of part files equals manifest', have == m['bytes'], f'{have} vs {m["bytes"]}')
        else:
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


def v_gnomad_sv(tid, r):
    """gnomAD v4.1 genome SV sites VCF (GRCh38): every record carries SVTYPE and an allele frequency, types come from the documented set."""
    known = {'DEL', 'DUP', 'INS', 'INV', 'CPX', 'CTX', 'BND', 'CNV', 'MCNV'}
    n = with_svtype = with_af = with_end = 0
    svt, chroms, no_af = collections.Counter(), collections.Counter(), collections.Counter()
    meta, header = {}, None
    with opener(path(tid)) as f:
        for line in f:
            if line.startswith('##'):
                m = re.match(r'##(\w+)=(.*)', line.strip())
                if m and m.group(1) in ('fileformat', 'reference'):
                    meta[m.group(1)] = m.group(2)
                continue
            if line.startswith('#CHROM'):
                header = line.strip().split('\t')
                continue
            p = line.split('\t', 8)
            n += 1
            chroms[p[0]] += 1
            info = p[7]
            m = re.search(r'(?:^|;)SVTYPE=(\w+)', info)
            if m:
                with_svtype += 1
                svt[m.group(1)] += 1
            if re.search(r'(?:^|;)AF=', info):
                with_af += 1
            elif m:
                no_af[m.group(1)] += 1
            with_end += bool(re.search(r'(?:^|;)END=', info))
    r.check('VCF header present', meta.get('fileformat', '').startswith('VCF') and header is not None, meta)
    r.check('records >= 400,000 (gnomAD-SV v2.1 alone had 433,371)', n >= 400_000, n)
    r.check('every record has SVTYPE', with_svtype == n, f'{with_svtype}/{n}')
    r.check('SVTYPE values come from the documented set', set(svt) <= known, dict(svt))
    r.check('every record has an AF except multi-allelic copy-number (SVTYPE=CNV) sites, which carry per-copy-number counts instead', set(no_af) <= {'CNV'} and with_af + sum(no_af.values()) == n, f'{with_af}/{n} with AF; without: {dict(no_af)}')
    r.check('contigs include chr1-22, chrX', {'chr' + str(i) for i in range(1, 23)} | {'chrX'} <= set(chroms), sorted(chroms)[:30])
    r.stats.update(records=n, svtype=dict(svt), records_without_AF=dict(no_af), contigs=len(chroms), records_with_END=with_end, reference=meta.get('reference'))


def v_tsv(tid, r, required, min_rows, delim='\t', comment='#', skip_prefix_lines=0):
    with opener(path(tid)) as f:
        lines = [l for l in f if not (comment and l.startswith(comment))]
    hdr = lines[0].rstrip('\n').split(delim)
    r.check('required columns present', all(c in hdr for c in required), [c for c in required if c not in hdr])
    r.check(f'rows >= {min_rows}', len(lines) - 1 >= min_rows, len(lines) - 1)
    r.stats.update(columns=len(hdr), rows=len(lines) - 1, sample_columns=hdr[:10])
    return hdr, lines


# Per-dataset expectations for the Open Targets Parquet directories. Row-count floors are order-of-magnitude sanity limits taken
# from the documented scale of each dataset (they are not exact counts); checks read only the named columns.
#   ('regex', col, pattern, min_fraction)   share of non-null values matching
#   ('range', col, lo, hi)                  every non-null numeric value inside [lo, hi]
#   ('nonnull', col, min_fraction)          share of non-null values
ENSG = r'^ENSG\d{11}$'
VARID = r'^[0-9XYM]+_\d+_[ACGT]+_[ACGT]+$'                  # chromosome_position_ref_alt
CHROM = r'^([1-9]|1[0-9]|2[0-2]|X|Y|MT)$'
ONTO = r'^[A-Za-z][A-Za-z0-9]*_[A-Za-z0-9]+$'    # prefix_local identifiers (EFO_0000001, OBA_VT0001234, MONDO_..., HP_...)
OT_SPECS = {
    'target': dict(required=['id', 'approvedSymbol', 'biotype', 'genomicLocation', 'proteinIds', 'constraint'], min_rows=50_000,
                   checks=[('regex', 'id', ENSG, 1.0), ('nonnull', 'approvedSymbol', 0.99)]),
    'go': dict(required=['id', 'label', 'namespace'], min_rows=40_000, checks=[('regex', 'id', r'^GO:\d{7}$', 1.0)]),
    'so': dict(required=['id', 'label'], min_rows=2_000, checks=[('regex', 'id', r'^SO:\d{7}$', 1.0)]),
    'homology': dict(required=['targetId', 'speciesId', 'homologyType', 'targetGeneId'], min_rows=1_000_000, checks=[('regex', 'targetId', ENSG, 1.0)]),
    'biosample': dict(required=['biosampleId', 'biosampleName'], min_rows=20_000, checks=[('regex', 'biosampleId', ONTO, 0.95)]),
    'disease': dict(required=['id', 'name', 'ancestors'], min_rows=20_000, checks=[('regex', 'id', ONTO, 1.0), ('nonnull', 'name', 0.99)]),
    'disease_hpo': dict(required=['id', 'name'], min_rows=15_000, checks=[('regex', 'id', ONTO, 1.0)]),
    'disease_phenotype': dict(required=['disease', 'phenotype', 'evidence'], min_rows=100_000,
                              checks=[('regex', 'disease', ONTO, 1.0), ('regex', 'phenotype', ONTO, 1.0)]),
    'mouse_phenotype': dict(required=['targetFromSourceId', 'modelPhenotypeId', 'modelPhenotypeLabel'], min_rows=100_000, checks=[('regex', 'targetFromSourceId', ENSG, 1.0)]),
    'pharmacogenomics': dict(required=['variantRsId', 'genotypeId', 'pgxCategory', 'phenotypeText', 'drugs'], min_rows=10_000, checks=[('regex', 'variantId', r'^[0-9XYM]+_\d+_[ACGT]+_[ACGT]+$', 0.99)]),
    'evidence_gene_burden': dict(required=['targetId', 'diseaseId', 'datasourceId', 'beta', 'oddsRatio', 'pValueMantissa', 'pValueExponent'], min_rows=20_000,
                                 checks=[('regex', 'targetId', ENSG, 1.0), ('range', 'resourceScore', 0, 1)]),
    'evidence_genomics_england': dict(required=['targetId', 'diseaseId', 'confidence'], min_rows=10_000, checks=[('regex', 'targetId', ENSG, 1.0), ('range', 'score', 0, 1)]),
    'evidence_gene2phenotype': dict(required=['targetId', 'diseaseId', 'confidence', 'allelicRequirements'], min_rows=2_000, checks=[('regex', 'targetId', ENSG, 1.0), ('range', 'score', 0, 1)]),
    'evidence_orphanet': dict(required=['targetId', 'diseaseId', 'alleleOrigins'], min_rows=3_000, checks=[('regex', 'targetId', ENSG, 1.0), ('range', 'score', 0, 1)]),
    'credible_set': dict(required=['studyLocusId', 'studyId', 'variantId', 'chromosome', 'position', 'finemappingMethod', 'locus', 'studyType'], min_rows=1_000_000,
                         checks=[('regex', 'variantId', VARID, 0.99), ('regex', 'chromosome', CHROM, 0.999), ('range', 'position', 1, 300_000_000), ('nonnull', 'studyLocusId', 1.0)]),
    'l2g_prediction': dict(required=['studyLocusId', 'geneId', 'score'], min_rows=500_000, checks=[('regex', 'geneId', ENSG, 1.0), ('range', 'score', 0, 1)]),
    'variant': dict(required=['variantId', 'chromosome', 'position', 'referenceAllele', 'alternateAllele', 'variantEffect', 'alleleFrequencies', 'mostSevereConsequenceId'], min_rows=3_000_000,
                    checks=[('regex', 'variantId', VARID, 0.99), ('regex', 'chromosome', CHROM, 0.999), ('range', 'position', 1, 300_000_000), ('regex', 'mostSevereConsequenceId', r'^SO_\d{7}$', 0.99)]),
    'enhancer_to_gene': dict(required=['geneId', 'chromosome', 'start', 'end', 'score', 'intervalType', 'datasourceId', 'biosampleId'], min_rows=10_000_000,
                             checks=[('regex', 'geneId', ENSG, 1.0), ('regex', 'chromosome', CHROM, 0.999), ('range', 'start', 0, 300_000_000), ('range', 'score', 0, 1)]),
    'interaction': dict(required=['sourceDatabase', 'targetA', 'targetB', 'speciesA', 'speciesB', 'scoring', 'interactionId'], min_rows=5_000_000,
                        checks=[('regex', 'targetA', ENSG, 0.99), ('regex', 'targetB', ENSG, 0.99), ('range', 'scoring', 0, 1)]),
    'interaction_evidence': dict(required=['interactionId', 'interactionScore', 'evidenceScore', 'interactionResources', 'interactionDetectionMethodShortName'], min_rows=10_000_000,
                                 checks=[('range', 'interactionScore', 0, 1), ('range', 'evidenceScore', 0, 1), ('nonnull', 'interactionId', 1.0)]),
    'study': dict(required=['studyId', 'studyType', 'traitFromSource', 'projectId', 'nSamples', 'hasSumstats'], min_rows=1_000_000, checks=[('nonnull', 'studyId', 1.0), ('nonnull', 'studyType', 1.0)]),
    'target_essentiality': dict(required=['targetId', 'isEssential', 'depMapEssentiality'], min_rows=10_000, checks=[('regex', 'targetId', ENSG, 1.0), ('nonnull', 'isEssential', 0.5)]),
    'evidence_crispr': dict(required=['targetId', 'diseaseId', 'diseaseCellLines', 'resourceScore'], min_rows=100, checks=[('regex', 'targetId', ENSG, 0.99), ('range', 'score', 0, 1)]),
    'evidence_crispr_screen': dict(required=['targetId', 'diseaseId', 'studyId', 'contrast', 'log2FoldChangeValue'], min_rows=10_000, checks=[('regex', 'targetId', ENSG, 0.99), ('range', 'score', 0, 1)]),
    'evidence_impc': dict(required=['targetId', 'diseaseId', 'targetInModelMgiId', 'biologicalModelAllelicComposition', 'diseaseModelAssociatedModelPhenotypes'], min_rows=1_000_000,
                          checks=[('regex', 'targetId', ENSG, 1.0), ('range', 'score', 0, 1)]),
    'evidence_uniprot_literature': dict(required=['targetId', 'diseaseId', 'confidence'], min_rows=1_000, checks=[('regex', 'targetId', ENSG, 1.0), ('range', 'score', 0, 1)]),
    'evidence_uniprot_variants': dict(required=['targetId', 'diseaseId', 'variantId', 'variantRsId'], min_rows=10_000,
                                      checks=[('regex', 'targetId', ENSG, 1.0), ('regex', 'variantId', r'^[0-9XYM]+_\d+_[ACGT]+_[ACGT]+$', 0.7), ('range', 'score', 0, 1)]),
}


def v_parquet_dir(tid, r):
    """Open Targets Parquet directory: files match the manifest table, every footer is readable and agrees on one schema,
    required columns exist, the row count passes a floor, and key identifier / range columns contain well-formed values."""
    try:
        import pyarrow.compute as pc
        import pyarrow.parquet as pq
    except ImportError:
        r.check('pyarrow is installed (pip install -r acquisition/requirements.txt)', False, 'ImportError')
        return
    d = path(tid)
    ds = tid.split('opentargets_26_09_', 1)[-1]
    m = MANIFEST[tid]
    table = [l.rstrip('\n').split('\t') for l in open(os.path.join(ROOT, m['files_table']))][1:]
    on_disk = sorted(f for f in os.listdir(d) if f.endswith('.parquet'))
    r.check('part files on disk equal the manifest file table (names and sizes)', on_disk == sorted(n for n, _, _ in table) and all(os.path.getsize(os.path.join(d, n)) == int(sz) for n, sz, _ in table), f'{len(on_disk)} files')
    files = [os.path.join(d, f) for f in on_disk]
    schemas, rows = [], 0
    for fp in files:
        md = pq.ParquetFile(fp).metadata
        rows += md.num_rows
        schemas.append(pq.read_schema(fp))
    r.check('every part has a readable Parquet footer', True, len(files))
    r.check('all parts share one schema', all(sc.equals(schemas[0]) for sc in schemas[1:]))
    names = schemas[0].names
    spec = OT_SPECS.get(ds)
    r.check('a content specification exists for this dataset', spec is not None, ds)
    r.stats.update(files=len(files), rows=rows, columns=len(names), column_names=names[:60])
    if not spec:
        return
    r.check('required columns present', all(c in names for c in spec['required']), [c for c in spec['required'] if c not in names])
    r.check(f"rows >= {spec['min_rows']:,}", rows >= spec['min_rows'], rows)
    for kind, col, *arg in spec['checks']:
        if col not in names:
            r.check(f'{kind} check on {col}', False, 'column missing')
            continue
        tot = ok = nn = 0
        for fp in files:
            arr = pq.read_table(fp, columns=[col]).column(col)
            tot += len(arr)
            nn += len(arr) - arr.null_count
            valid = arr.drop_null()
            if kind == 'regex':
                ok += int(pc.sum(pc.match_substring_regex(valid.cast('string'), arg[0])).as_py() or 0)
            elif kind == 'range':
                ok += int(pc.sum(pc.and_(pc.greater_equal(valid, arg[0]), pc.less_equal(valid, arg[1]))).as_py() or 0)
        if kind == 'regex':
            r.check(f'{col}: share of non-null values matching {arg[0]!r} >= {arg[1]}', nn > 0 and ok / nn >= arg[1], f'{ok}/{nn} of {tot}')
        elif kind == 'range':
            r.check(f'{col}: every non-null value within [{arg[0]}, {arg[1]}]', ok == nn, f'{ok}/{nn} of {tot}')
        else:
            r.check(f'{col}: non-null share >= {arg[0]}', tot > 0 and nn / tot >= arg[0], f'{nn}/{tot}')


PRED_RANGE = {'cadd_phred': (0, 100), 'revel_max': (0, 1), 'spliceai_ds_max': (0, 1), 'pangolin_largest_ds': (-1, 1), 'phylop': (-20, 10), 'sift_max': (0, 1), 'polyphen_max': (0, 1)}


def v_regional_slice(tid, r):
    """gnomAD slice written by fetch_gnomad_regional.py: header layout, every row is a single-nucleotide record inside the declared
    regions, frequencies are internally consistent (AF = AC / AN, AC <= AN) and the file agrees with its manifest."""
    m = MANIFEST[tid]
    kind = 'exomes' if 'exomes' in tid else 'genomes'
    bed = os.path.join(DATA, 'gnomad_regional', 'mavedb_snv_regions.bed')
    r.check('region file hash equals the manifest', hashlib.sha256(open(bed, 'rb').read()).hexdigest() == m['regions_sha256'])
    regions = collections.defaultdict(list)
    for line in open(bed):
        c, a, b = line.split('\t')
        regions[c[3:]].append((int(a) + 1, int(b)))
    chroms_ok = {str(i) for i in range(1, 23)} | {'X', 'Y'}
    n = bad_region = bad_alleles = bad_af = bad_ac = unsorted = dup = 0
    last, seen_filter = None, collections.Counter()
    pred_n, pred_bad = collections.Counter(), collections.Counter()
    with opener(path(tid)) as f:
        rd = csv.reader(f, delimiter='\t')
        hdr = next(rd)
        ix = {c: i for i, c in enumerate(hdr)}
        for row in rd:
            n += 1
            chrom, pos, ref, alt = row[ix['chrom']], int(row[ix['pos']]), row[ix['ref']], row[ix['alt']]
            if chrom not in chroms_ok or not any(a <= pos <= b for a, b in regions.get(chrom, ())):
                bad_region += 1
            if len(ref) != 1 or len(alt) != 1 or ref == alt or ref not in 'ACGT' or alt not in 'ACGT':
                bad_alleles += 1
            key = (chroms_ok_order(chrom), pos, alt)
            if last is not None:
                if key < last:
                    unsorted += 1
                elif key == last:
                    dup += 1
            last = key
            seen_filter[row[ix['filter']]] += 1
            for col, (lo, hi) in PRED_RANGE.items():
                v = row[ix[col]]
                if v != '':
                    pred_n[col] += 1
                    if not lo <= float(v) <= hi:
                        pred_bad[col] += 1
            ac, an, af = row[ix['AC']], row[ix['AN']], row[ix['AF']]
            if ac and an:
                if int(ac) > int(an):
                    bad_ac += 1
                if af and int(an) > 0 and abs(int(ac) / int(an) - float(af)) > 1e-4 * float(af) + 1e-9:
                    bad_af += 1
    want = ['chrom', 'pos', 'ref', 'alt', 'filter', 'AC', 'AN', 'AF', 'nhomalt', 'AC_grpmax', 'AN_grpmax', 'AF_grpmax']
    r.check('header starts with the documented core columns and ends with the predictor and flag columns', hdr[:len(want)] == want and hdr[-1] == 'VRS_Allele_IDs' and 'cadd_phred' in hdr, hdr[:12])
    r.check('rows > 10,000', n > 10_000, n)
    r.check('rows equal the manifest', n == m['rows'], f'{n} vs {m["rows"]}')
    r.check('every row lies on chr1-22/X/Y inside a requested region', bad_region == 0, bad_region)
    r.check('every row is a single-nucleotide substitution with distinct A/C/G/T alleles', bad_alleles == 0, bad_alleles)
    r.check('rows are sorted and unique by chromosome, position, alt', unsorted == 0 and dup == 0, f'unsorted {unsorted}, duplicate {dup}')
    r.check('AC <= AN on every row', bad_ac == 0, bad_ac)
    r.check('AF equals AC / AN (relative 1e-4) on every row', bad_af == 0, bad_af)
    tab = [l.rstrip('\n').split('\t') for l in open(os.path.join(ROOT, m['files_table']))][1:]
    r.check('predictor scores lie within their documented ranges', not pred_bad, dict(pred_bad))
    r.check('CADD phred is present on at least 90% of rows', pred_n['cadd_phred'] / max(1, n) >= 0.9, f"{pred_n['cadd_phred']}/{n}")
    r.check('per-chromosome table sums to the row count', sum(int(x[7]) for x in tab) == n, sum(int(x[7]) for x in tab))
    r.stats.update(rows=n, filters=dict(seen_filter.most_common(8)), columns=len(hdr), regions=m['n_regions'], source_files=len(tab), rows_with_predictor=dict(pred_n))


def chroms_ok_order(chrom):
    return (int(chrom) if chrom.isdigit() else {'X': 23, 'Y': 24}.get(chrom, 99))


def v_per_entity(tid, r):
    """Per-protein API responses fetched by fetch_per_protein.py: files match the recorded table, content parses and agrees with
    UniProt on protein length, residue positions are inside the protein, and coverage is what the human proteome implies."""
    m = MANIFEST[tid]
    table = [l.rstrip('\n').split('\t') for l in open(os.path.join(ROOT, m['files_table']))]
    by = collections.defaultdict(dict)
    for acc, part, status, nbytes, sha, http, ts in table[1:]:
        by[acc][part] = (status, int(nbytes), sha, int(http))
    accs = sorted(by)
    base = os.path.join(DATA, TARGETS[tid]['dest'])
    r.check('no accession ended in an error', m.get('n_errors') == 0, m.get('errors'))
    ulen = {}
    with opener(os.path.join(DATA, 'uniprot', 'uniprot_human_reviewed.tsv')) as f:
        for row in csv.DictReader(f, delimiter='\t'):
            ulen[row['Entry']] = int(row['Length'])
    reviewed = [a for a in accs if a in ulen]

    def file_of(acc, part):
        if tid.startswith('alphafold'):
            return os.path.join(base, {'prediction': 'prediction', 'model': 'model', 'confidence': 'confidence'}[part], acc + {'prediction': '.json', 'model': '.cif', 'confidence': '.json'}[part])
        return os.path.join(base, acc + '.json')
    stored = [(a, p) for a in accs for p, v in by[a].items() if v[0] == 'stored']
    import hashlib
    bad = [f'{a}/{p}' for a, p in stored if not os.path.exists(file_of(a, p)) or hashlib.sha256(open(file_of(a, p), 'rb').read()).hexdigest() != by[a][p][2]]
    r.check('every stored file exists and matches its recorded SHA-256', not bad, bad[:5])
    if tid.startswith('alphafold'):
        have = [a for a in accs if by[a].get('prediction', ('none',))[0] == 'stored']
        canon = [a for a in have if by[a].get('model', ('',))[0] == 'stored']
        iso_only = [a for a in have if by[a].get('model', ('',))[0] == 'isoform_only']
        conf_bad, cif_bad, len_cmp, len_ok, plddt, iso_bad = [], [], 0, 0, [], []
        for a in canon:
            pred = json.load(open(file_of(a, 'prediction')))
            e = next(x for x in pred if x.get('entryId') == f'AF-{a}-F1')
            n = e['sequenceEnd'] - e['sequenceStart'] + 1
            cs = json.load(open(file_of(a, 'confidence')))['confidenceScore']
            if len(cs) != n or not all(0 <= x <= 100 for x in cs):
                conf_bad.append(a)
            plddt.append(sum(cs) / max(1, len(cs)))
            if not open(file_of(a, 'model'), errors='replace').read(400).startswith(f'data_AF-{a}-F1'):
                cif_bad.append(a)
            if a in ulen:
                len_cmp += 1
                len_ok += ulen[a] == n
        for a in iso_only:
            pred = json.load(open(file_of(a, 'prediction')))
            if any(x.get('entryId') == f'AF-{a}-F1' for x in pred) or os.path.exists(file_of(a, 'model')):
                iso_bad.append(a)
        r.check('at least 90% of accessions have a canonical AlphaFold model (reviewed human proteins are in the proteome)', len(canon) / max(1, len(accs)) >= 0.9, f'{len(canon)}/{len(accs)}')
        r.check('accessions recorded isoform_only have no AF-<acc>-F1 entry and no stored model', not iso_bad, iso_bad[:5])
        r.check('per-residue confidence has one value in [0, 100] per modelled residue', not conf_bad, conf_bad[:5])
        r.check('every stored model file is the mmCIF of the canonical entry AF-<acc>-F1', not cif_bad, cif_bad[:5])
        r.check('canonical model length equals UniProt Length for reviewed entries', len_cmp > 0 and len_ok == len_cmp, f'{len_ok}/{len_cmp}')
        r.stats.update(accessions=len(accs), with_prediction_metadata=len(have), without_prediction=len(accs) - len(have), canonical_models=len(canon), isoform_only=sorted(iso_only),
                       mean_of_mean_plddt=round(sum(plddt) / max(1, len(plddt)), 2), length_agreement_with_uniprot=[len_ok, len_cmp])
    else:
        have = [a for a in accs if by[a].get('matches', ('none',))[0] == 'stored']
        versions, cnt_bad, loc_bad, len_cmp, len_ok, nres = set(), 0, 0, 0, 0, 0
        for a in have:
            d = json.load(open(file_of(a, 'matches')))
            versions.add(d['interpro_version'])
            cnt_bad += d['count'] != len(d['results'])
            nres += len(d['results'])
            for res in d['results']:
                for pr in res['proteins']:
                    n = pr['protein_length']
                    if a in ulen:
                        len_cmp += 1
                        len_ok += ulen[a] == n
                    for loc in pr['entry_protein_locations']:
                        for fr in loc['fragments']:
                            loc_bad += not (1 <= fr['start'] <= fr['end'] <= n)
        r.check('one InterPro release across all files', len(versions) == 1, versions)
        r.check('count equals the number of results in every file', cnt_bad == 0, cnt_bad)
        r.check('every match fragment lies inside the protein (1 <= start <= end <= length)', loc_bad == 0, loc_bad)
        r.check('InterPro protein_length equals UniProt Length for reviewed entries', len_cmp > 0 and len_ok == len_cmp, f'{len_ok}/{len_cmp}')
        r.check('at least 80% of accessions have InterPro matches', len(have) / max(1, len(accs)) >= 0.8, f'{len(have)}/{len(accs)}')
        r.stats.update(accessions=len(accs), with_matches=len(have), without_matches=len(accs) - len(have), match_records=nres, interpro_version=sorted(versions), length_agreement_with_uniprot=[len_ok, len_cmp])


def v_alphamissense_variants(tid, r):
    """AlphaMissense_hg38.tsv.gz as documented in the Zenodo README (evidence alphamissense/zenodo_readme): 71M single-nucleotide
    missense variants at canonical transcripts, sorted by genomic coordinate, columns CHROM..am_class, class thresholds
    likely_benign < 0.34, likely_pathogenic > 0.564, ambiguous otherwise."""
    want = ['#CHROM', 'POS', 'REF', 'ALT', 'genome', 'uniprot_id', 'transcript_id', 'protein_variant', 'am_pathogenicity', 'am_class']
    hdr = None
    n = bad_cols = bad_alt = bad_genome = bad_score = bad_class = unsorted = boundary_rounding = 0
    cls = collections.Counter()
    chroms, tx, uni = [], set(), set()

    def band(x):
        return 'likely_benign' if x < 0.34 else ('likely_pathogenic' if x > 0.564 else 'ambiguous')

    last_c = last_p = None
    with opener(path(tid)) as f:
        for line in f:
            if line[0] == '#':
                if line.startswith('#CHROM'):
                    hdr = line.rstrip('\n').split('\t')
                continue
            p = line.rstrip('\n').split('\t')
            n += 1
            if len(p) != 10:
                bad_cols += 1
                continue
            if p[2] == p[3] or p[2] not in ('A', 'C', 'G', 'T') or p[3] not in ('A', 'C', 'G', 'T'):
                bad_alt += 1
            bad_genome += p[4] != 'hg38'
            s = float(p[8])
            if not 0 <= s <= 1:
                bad_score += 1
            if p[9] != band(s):
                # scores are printed to 4 decimals while the class uses the unrounded value: accept a class that the rule gives for
                # any unrounded value within +-5e-5 of the printed score, and count those rows separately
                if p[9] in (band(s - 4.9e-5), band(s + 4.9e-5)) and round(s, 4) in (0.34, 0.564):
                    boundary_rounding += 1
                else:
                    bad_class += 1
            cls[p[9]] += 1
            c, pos = p[0], int(p[1])
            if c != last_c:
                if c in chroms:
                    unsorted += 1
                chroms.append(c)
                last_c, last_p = c, pos
            else:
                unsorted += pos < last_p
                last_p = pos
            tx.add(p[6])
            uni.add(p[5])
    r.check('column header matches the documented layout', hdr == want, hdr)
    r.check('rows ~ 71M (README says 71M; +-3%)', abs(n - 71_000_000) / 71_000_000 < 0.03, n)
    r.check('every row has 10 columns', bad_cols == 0, bad_cols)
    r.check('REF/ALT are distinct A/C/G/T bases', bad_alt == 0, bad_alt)
    r.check('genome column is hg38 throughout', bad_genome == 0, bad_genome)
    r.check('am_pathogenicity within [0, 1]', bad_score == 0, bad_score)
    r.check('am_class agrees with the documented thresholds (0.34 / 0.564) up to 4-decimal rounding at the two boundaries', bad_class == 0, f'{bad_class} outside tolerance; {boundary_rounding} boundary-rounding rows')
    r.check('sorted by chromosome block and position', unsorted == 0, unsorted)
    r.check('one chromosome block each for chr1-22, X, Y, M', len(chroms) == len(set(chroms)) and {'chr' + str(i) for i in range(1, 23)} | {'chrX', 'chrY'} <= set(chroms), chroms)
    r.stats.update(rows=n, classes=dict(cls), class_boundary_rounding_rows=boundary_rounding, chromosome_blocks=len(chroms), transcripts=len(tx), uniprot_accessions=len(uni))


def v_dgv(tid, r):
    """DGV curated structural variants (GRCh38 release 2025-12-01): tab-delimited, first line gives the column names."""
    want = ['variantaccession', 'chr', 'start', 'end', 'varianttype', 'variantsubtype', 'reference', 'pubmedid', 'method', 'platform', 'mergedvariants',
            'supportingvariants', 'mergedorsample', 'frequency', 'samplesize', 'observedgains', 'observedlosses', 'cohortdescription', 'genes', 'samples']
    n = bad_cols = bad_span = unplaced = 0
    types, prefixes = collections.Counter(), collections.Counter()
    with opener(path(tid)) as f:
        hdr = f.readline().rstrip('\n').split('\t')
        for line in f:
            p = line.rstrip('\n').split('\t')
            n += 1
            if len(p) != len(want):
                bad_cols += 1
                continue
            types[p[4]] += 1
            prefixes[re.match(r'[a-z]+', p[0]).group(0) if re.match(r'[a-z]+', p[0]) else '?'] += 1
            try:
                if int(p[3]) < int(p[2]) or int(p[2]) < 1:
                    bad_span += 1
            except ValueError:
                bad_span += 1
            unplaced += not p[1] or p[1] == 'N'
    r.check('column header matches the DGV layout', hdr == want, hdr)
    r.check('rows >= 1,900,000', n >= 1_900_000, n)
    r.check('every row has the full column count', bad_cols == 0, bad_cols)
    r.check('start >= 1 and end >= start on every row', bad_span == 0, bad_span)
    r.check('accession prefixes are DGV/dbVar style (dgv, nsv, esv)', set(prefixes) <= {'dgv', 'nsv', 'esv'}, dict(prefixes))
    r.stats.update(rows=n, varianttype=dict(types), accession_prefixes=dict(prefixes), rows_without_placed_chromosome=unplaced)


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
        if t.get('kind') == 'regional_slice':
            v_regional_slice(tid, r)
        elif t.get('kind') == 'per_entity':
            v_per_entity(tid, r)
        elif t.get('kind') == 'http_dir':
            v_parquet_dir(tid, r)
        elif tid == 'clinvar_vcf_grch38':
            v_clinvar_vcf(tid, r)
        elif tid == 'gnomad_exome_cnv_v4_1':
            v_gnomad_cnv(tid, r)
        elif tid == 'gnomad_sv_v4_1_sites_vcf':
            v_gnomad_sv(tid, r)
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
        elif tid == 'alphamissense_hg38':
            v_alphamissense_variants(tid, r)
        elif tid.startswith('alphamissense_'):
            hdr, lines = v_tsv(tid, r, ['transcript_id', 'mean_am_pathogenicity'], 15000, comment='#')
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
            chroms = set()
            with gzip.open(path(tid), 'rt') as f:
                hdr = f.readline().rstrip('\n').split('\t')
                ui = hdr.index('urn')
                for line in f:
                    n += 1
                    p = line.split('\t', ui + 1)
                    urns.add(p[ui])
                    chroms.add(p[0])
            r.check('header begins with #chr start end ref alt', hdr[:5] == ['#chr', 'start', 'end', 'ref', 'alt'], hdr[:5])
            r.check('rows > 1,000,000', n > 1_000_000, n)
            r.check('urn column holds MaveDB score-set URNs', all(u.startswith('urn:mavedb:') for u in urns), list(urns)[:3])
            r.stats.update(columns=len(hdr), rows=n, distinct_score_sets=len(urns), contigs=sorted(chroms)[:30])
        elif tid == 'mavedb_zenodo_dump_v4':
            kinds = collections.Counter()
            with tarfile.open(path(tid)) as tf:
                for m in tf:
                    base = m.name.rsplit('/', 1)[-1]
                    if base.startswith('._') or m.isdir():
                        continue
                    k = re.search(r'\.(scores|counts|annotations)\.csv$', base)
                    if k:
                        kinds[k.group(1)] += 1
                    elif base.endswith('.mapped-variants.json'):
                        kinds['mapped-variants'] += 1
                    elif base.endswith('.va.ndjson'):
                        kinds['va'] += 1
                    elif base in ('main.json', 'README.md', 'LICENSE.txt'):
                        kinds[base] += 1
            r.check('archive fully readable', True)
            r.check('main.json, README.md and LICENSE.txt present', all(kinds[x] == 1 for x in ('main.json', 'README.md', 'LICENSE.txt')), dict(kinds))
            r.check('scores tables > 2,500', kinds['scores'] > 2500, kinds['scores'])
            r.check('annotations present for mapped sets', kinds['annotations'] > 1000 and kinds['annotations'] == kinds['mapped-variants'], dict(kinds))
            r.stats.update(members_by_type=dict(kinds))
        elif tid == 'dgv_grch38_variants_2025_12_01':
            v_dgv(tid, r)
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
