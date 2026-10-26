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
import struct
import sys
import tarfile
import zipfile

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import mavedb_snvs  # noqa: E402
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
    regions, frequencies are internally consistent (AF = AC / AN, AC <= AN) and the file agrees with its manifest. The joint callset
    must also add up: AC_joint = AC_exomes + AC_genomes and AN_joint = AN_exomes + AN_genomes wherever both callsets have the site."""
    m = MANIFEST[tid]
    kind = 'joint' if '_joint_' in tid else ('exomes' if 'exomes' in tid else 'genomes')
    sfx = '_joint' if kind == 'joint' else ''
    bed = os.path.join(DATA, 'gnomad_regional', 'mavedb_snv_regions.bed')
    r.check('region file hash equals the manifest', hashlib.sha256(open(bed, 'rb').read()).hexdigest() == m['regions_sha256'])
    regions = collections.defaultdict(list)
    for line in open(bed):
        c, a, b = line.split('\t')
        regions[c[3:]].append((int(a) + 1, int(b)))
    chroms_ok = {str(i) for i in range(1, 23)} | {'X', 'Y'}
    n = bad_region = bad_alleles = bad_af = bad_ac = unsorted = dup = 0
    bad_pool = pooled = 0
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
            if kind != 'joint':
                for col, (lo, hi) in PRED_RANGE.items():
                    v = row[ix[col]]
                    if v != '':
                        pred_n[col] += 1
                        if not lo <= float(v) <= hi:
                            pred_bad[col] += 1
            ac, an, af = row[ix['AC' + sfx]], row[ix['AN' + sfx]], row[ix['AF' + sfx]]
            if ac and an:
                if int(ac) > int(an):
                    bad_ac += 1
                if af and int(an) > 0 and abs(int(ac) / int(an) - float(af)) > 1e-4 * float(af) + 1e-9:
                    bad_af += 1
            if kind == 'joint':
                parts = [(row[ix['AC_' + k]], row[ix['AN_' + k]]) for k in ('exomes', 'genomes')]
                have = [(int(a_), int(b_)) for a_, b_ in parts if a_ != '' and b_ != '']
                if have and ac and an:
                    pooled += 1
                    if int(ac) != sum(a_ for a_, _ in have) or int(an) != sum(b_ for _, b_ in have):
                        bad_pool += 1
    if kind == 'joint':
        want = ['chrom', 'pos', 'ref', 'alt', 'filter', 'AC_joint', 'AN_joint', 'AF_joint', 'nhomalt_joint', 'AC_grpmax_joint', 'AN_grpmax_joint', 'AF_grpmax_joint', 'AC_exomes', 'AN_exomes', 'AF_exomes', 'AC_genomes', 'AN_genomes', 'AF_genomes']
        r.check('header starts with the documented joint columns and ends with the per-callset filter columns', hdr[:len(want)] == want and hdr[-1] == 'not_called_in_genomes', hdr[:18])
    else:
        want = ['chrom', 'pos', 'ref', 'alt', 'filter', 'AC', 'AN', 'AF', 'nhomalt', 'AC_grpmax', 'AN_grpmax', 'AF_grpmax']
        r.check('header starts with the documented core columns and ends with the predictor and flag columns', hdr[:len(want)] == want and hdr[-1] == 'VRS_Allele_IDs' and 'cadd_phred' in hdr, hdr[:12])
    r.check('rows > 10,000', n > 10_000, n)
    r.check('rows equal the manifest', n == m['rows'], f'{n} vs {m["rows"]}')
    r.check('every row lies on chr1-22/X/Y inside a requested region', bad_region == 0, bad_region)
    r.check('every row is a single-nucleotide substitution with distinct A/C/G/T alleles', bad_alleles == 0, bad_alleles)
    r.check('rows are sorted and unique by chromosome, position, alt', unsorted == 0 and dup == 0, f'unsorted {unsorted}, duplicate {dup}')
    r.check('AC <= AN on every row', bad_ac == 0, bad_ac)
    r.check('AF equals AC / AN (relative 1e-4) on every row', bad_af == 0, bad_af)
    if kind == 'joint':
        r.check('AC_joint and AN_joint equal the sums of the exome and genome counts present on the row', pooled > 0 and bad_pool == 0, f'{bad_pool} of {pooled}')
    else:
        r.check('predictor scores lie within their documented ranges', not pred_bad, dict(pred_bad))
        r.check('CADD phred is present on at least 90% of rows', pred_n['cadd_phred'] / max(1, n) >= 0.9, f"{pred_n['cadd_phred']}/{n}")
    tab = [l.rstrip('\n').split('\t') for l in open(os.path.join(ROOT, m['files_table']))][1:]
    r.check('per-chromosome table sums to the row count', sum(int(x[7]) for x in tab) == n, sum(int(x[7]) for x in tab))
    r.stats.update(rows=n, filters=dict(seen_filter.most_common(8)), columns=len(hdr), regions=m['n_regions'], source_files=len(tab), rows_with_predictor=dict(pred_n), rows_checked_for_pooling=pooled)


def chroms_ok_order(chrom):
    return (int(chrom) if chrom.isdigit() else {'X': 23, 'Y': 24}.get(chrom, 99))


HUMSAVAR_ROW = re.compile(r'^(\S+)\s+([OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9](?:[A-Z][A-Z0-9]{2}[0-9]){1,2})\s+(VAR_\d{6})\s+p\.([A-Z][a-z]{2})(\d+)([A-Z][a-z]{2})\s+(LP/P|LB/B|US)\s+(\S+)\s+(.*)$')


def v_humsavar(tid, r):
    """UniProt humsavar.txt: the header states the release and per-category counts; every data row must parse and the counts must agree with the header."""
    stated, cats, n, bad, dup = {}, collections.Counter(), 0, 0, 0
    release, footer = None, False
    seen, ftids = set(), set()
    in_table = False
    with opener(path(tid)) as f:
        for line in f:
            m = re.match(r'\s*(LP/P|LB/B|US|Total):\s+(\d+)', line)
            if m and not in_table:
                stated[m.group(1)] = int(m.group(2))
            m = re.match(r'Release:\s+(\S+ of \S+)', line)
            if m:
                release = m.group(1)
            if line.startswith('_________'):
                in_table = True
                continue
            if not in_table or not line.strip():
                continue
            if line.startswith('---'):
                in_table = False          # footer: licence statement
                footer = True
                continue
            g = HUMSAVAR_ROW.match(line.rstrip('\n'))
            if not g:
                bad += 1
                continue
            n += 1
            cats[g.group(7)] += 1
            key = (g.group(3), g.group(9))     # a variant with several diseases has one row per disease
            dup += key in seen
            seen.add(key)
            ftids.add(g.group(3))
    r.check('release line present', release is not None, release)
    r.check('every data row parses (gene, accession, VAR id, p. change, category, dbSNP, disease)', bad == 0 and n > 0, f'{bad} unparsed of {n + bad}')
    r.check('row count equals the Total stated in the file header', n == stated.get('Total'), f'{n} vs {stated.get("Total")}')
    r.check('per-category counts equal the header statistics', all(cats[k] == stated.get(k) for k in ('LP/P', 'LB/B', 'US')), f'{dict(cats)} vs {stated}')
    r.check('(FTId, disease) pairs are unique and the licence footer is present', dup == 0 and footer, f'{dup} repeated pairs; footer {footer}')
    r.stats.update(release=release, rows=n, distinct_ftids=len(ftids), categories=dict(cats), stated=stated)


def v_reference_slice(tid, r):
    """Reference sequence regions written by fetch_reference_regional.py: every row has the requested length and only A/C/G/T/N,
    rows are sorted and disjoint, the regions hash to the manifest, and every region lies inside a sequence listed in the saved .fai."""
    m = MANIFEST[tid]
    n = bases = bad_len = bad_chars = unsorted = overlap = beyond = unknown = 0
    last = None
    regions = []
    odd = collections.Counter()
    fai, fai_sha = {}, None
    fai_copy = os.path.join(DATA, 'reference_sequence', 'index_cache', os.path.basename(m['url']) + '.fai')
    if os.path.exists(fai_copy):
        raw = open(fai_copy, 'rb').read()
        fai_sha = hashlib.sha256(raw).hexdigest()
        for line in raw.decode().splitlines():
            p = line.split('\t')
            fai[p[0]] = int(p[1])
    with opener(path(tid)) as f:
        rd = csv.reader(f, delimiter='\t')
        hdr = next(rd)
        for c, s_, e_, seq in rd:
            n += 1
            s_, e_ = int(s_), int(e_)
            bases += len(seq)
            bad_len += len(seq) != e_ - s_ + 1 or s_ < 1
            odd.update(set(seq) - set('ACGTN'))
            key = (mavedb_snvs.contig_rank(c), s_)
            if last is not None:
                unsorted += key < last[:2]
                overlap += (key[0] == last[0] and s_ <= last[2])
            last = (key[0], s_, e_)
            unknown += c not in fai
            beyond += (c in fai and e_ > fai[c])
            regions.append(f'{c}\t{s_}\t{e_}')
    bad_chars = sum(odd.values())
    r.check('header is chrom, start, end, sequence', hdr == ['chrom', 'start', 'end', 'sequence'], hdr)
    r.check('rows and bases equal the manifest', n == m['rows'] and bases == m['bases'], f'{n} rows, {bases} bases vs {m["rows"]}, {m["bases"]}')
    r.check('the regions hash to the manifest value (the file is the set of regions that was requested and read)', hashlib.sha256('\n'.join(regions).encode()).hexdigest() == m.get('regions_sha256'), f'{n} regions')
    r.check('every sequence has the length of its region (1-based inclusive)', bad_len == 0, bad_len)
    r.check('sequences contain only A, C, G, T, N', bad_chars == 0, dict(odd))
    r.check('regions are sorted and disjoint', unsorted == 0 and overlap == 0, f'unsorted {unsorted}, overlapping {overlap}')
    r.check('the saved .fai equals the manifest hash and lists the 25 standard sequences with the GRCh38 lengths of chr1 and MT',
            fai_sha is not None and fai_sha == m.get('fai_sha256') and set(mavedb_snvs.STANDARD) <= set(fai) and fai.get('1') == 248956422 and fai.get('MT') == 16569,
            {'sha256_matches': fai_sha == m.get('fai_sha256'), '1': fai.get('1'), 'MT': fai.get('MT'), 'sequences': len(fai)})
    r.check('every region is on a sequence of the .fai and ends inside it', n > 0 and unknown == 0 and beyond == 0, f'on unknown sequence {unknown}, beyond the end {beyond}')
    r.stats.update(rows=n, bases=bases, sequences_read=m.get('sequences_read'), skipped_contigs=m.get('skipped_contigs'), fai_sequences=len(fai))


FASTA_ID = 'ensembl_grch38_dna_toplevel'


def read_fai(p):
    """samtools faidx index: name, length, byte offset of the first base, bases per line, bytes per line."""
    out = []
    for line in open(p):
        name, length, offset, bases, width = line.rstrip('\n').split('\t')[:5]
        out.append((name, int(length), int(offset), int(bases), int(width)))
    return out


def v_fai(tid, r):
    """The .fai published with the FASTA: unique sequence names, consistent line geometry, and byte-identical to the copy that
    fetch_reference_regional.py used for its remote reads (so the remote reads and this file describe the same index)."""
    fai = read_fai(path(tid))
    names = [x[0] for x in fai]
    r.check('one row per sequence with unique names and the 25 standard sequences', len(set(names)) == len(names) and set(mavedb_snvs.STANDARD) <= set(names), len(fai))
    r.check('line geometry is consistent (bytes per line = bases per line + 1)', all(x[4] == x[3] + 1 and x[3] > 0 for x in fai), {x[3] for x in fai})
    r.check('offsets are increasing', all(a[2] < b[2] for a, b in zip(fai, fai[1:])))
    r.check('identical to the .fai that the regional reads used', hashlib.sha256(open(path(tid), 'rb').read()).hexdigest() == MANIFEST.get('ensembl_grch38_reference_regional', {}).get('fai_sha256'))
    r.stats.update(sequences=len(fai), total_bases=sum(x[1] for x in fai), standard_bases=sum(x[1] for x in fai if x[0] in mavedb_snvs.STANDARD))


def v_gzi(tid, r):
    """bgzip block index: an 8-byte block count, then one (compressed offset, uncompressed offset) pair per block after the first.
    Every pair must describe the real block chain of the FASTA: each compressed offset is a BGZF block start, block sizes add up
    to the next offset, and the uncompressed sizes (ISIZE) add up to the next uncompressed offset."""
    raw = open(path(tid), 'rb').read()
    n = struct.unpack('<Q', raw[:8])[0]
    r.check('size equals 8 + 16 x block count', len(raw) == 8 + 16 * n, f'{len(raw)} bytes, {n} blocks')
    pairs = struct.unpack(f'<{2 * n}Q', raw[8:8 + 16 * n])
    comp, unc = [0] + list(pairs[0::2]), [0] + list(pairs[1::2])
    r.check('offsets are strictly increasing', all(a < b for a, b in zip(comp, comp[1:])) and all(a < b for a, b in zip(unc, unc[1:])))
    bad_magic = bad_chain = bad_isize = 0
    with open(path(FASTA_ID), 'rb') as f:
        for i, c in enumerate(comp):
            f.seek(c)
            h = f.read(18)
            if not (h[:4] == b'\x1f\x8b\x08\x04' and h[12:14] == b'BC'):
                bad_magic += 1
                continue
            bsize = struct.unpack('<H', h[16:18])[0] + 1
            f.seek(c + bsize - 4)
            isize = struct.unpack('<I', f.read(4))[0]
            if i + 1 < len(comp):
                bad_chain += comp[i + 1] != c + bsize
                bad_isize += unc[i + 1] != unc[i] + isize
    r.check('every compressed offset starts a BGZF block', bad_magic == 0, bad_magic)
    r.check('block sizes chain exactly to the next compressed offset', bad_chain == 0, bad_chain)
    r.check('uncompressed block sizes chain exactly to the next uncompressed offset', bad_isize == 0, bad_isize)
    r.stats.update(blocks=n + 1, last_uncompressed_offset=unc[-1])


def v_fasta_toplevel(tid, r):
    """Scan the whole bgzip FASTA and compare it with the published .fai: record order, names, lengths, byte offsets and line
    geometry must all agree; only upper-case A/C/G/T/N may occur (the file is unmasked); the two lengths that define the assembly
    (chr1 and MT) must match the GRCh38 values."""
    fai = read_fai(path(tid + '_fai'))
    r.check('companion .fai and .gzi are acquired and verified', all(MANIFEST.get(tid + x, {}).get('ok') and MANIFEST[tid + x].get('bsd_sum_verified_against_source') for x in ('_fai', '_gzi')))
    m = MANIFEST[tid]
    r.check('BSD checksum and size equal the values published in the Ensembl CHECKSUMS file', m.get('bsd_sum_verified_against_source') is True and m.get('size_verified_against_source') is True, m.get('bsd_sum'))
    records, odd, n_bases = [], collections.Counter(), collections.Counter()
    bad_line = bad_eol = 0
    cur = None
    pos = 0
    with gzip.open(path(tid), 'rb') as f:
        for line in f:
            n = len(line)
            if line[0] == 62:                                   # '>'
                if cur:
                    records.append(tuple(cur[:5]))
                cur = [line[1:].split(None, 1)[0].decode(), 0, pos + n, None, None, False]      # name, length, offset, bases per line, bytes per line, saw a short line
            else:
                bad_eol += line[-1] != 10
                w = n - 1
                if cur[3] is None:
                    cur[3], cur[4] = w, n
                if cur[5] or w > cur[3]:                           # a short line must be the last line of its record
                    bad_line += 1
                if w < cur[3]:
                    cur[5] = True
                cur[1] += w
                left = line.translate(None, b'ACGTN\n')
                if left:
                    odd.update(left.decode('ascii', 'replace'))
                n_bases['N'] += line.count(b'N')
            pos += n
    records.append(tuple(cur[:5]))
    n_bases['all'] = sum(x[1] for x in records)
    r.check('the sequence has the records, names, lengths, byte offsets and line geometry of the .fai', records == fai, f'{len(records)} records vs {len(fai)} in the .fai')
    r.check('every line has the index width except the last of its record, and every line ends with a newline', bad_line == 0 and bad_eol == 0, f'short lines inside a record {bad_line}, lines without newline {bad_eol}')
    r.check('only upper-case A, C, G, T, N occur (unmasked)', not odd, dict(odd))
    lens = {x[0]: x[1] for x in records}
    r.check('chr1 and MT have the GRCh38 lengths (248,956,422 and 16,569)', lens.get('1') == 248956422 and lens.get('MT') == 16569, {k: lens.get(k) for k in ('1', 'MT')})
    r.stats.update(records=len(records), bases=n_bases['all'], n_bases=n_bases['N'], standard_bases=sum(lens[k] for k in mavedb_snvs.STANDARD if k in lens), uncompressed_bytes=pos)


GTF_FEATURES = {'gene', 'transcript', 'exon', 'CDS', 'five_prime_utr', 'three_prime_utr', 'start_codon', 'stop_codon', 'Selenocysteine'}
GTF_GID, GTF_TID, GTF_TAG, GTF_BIO, GTF_TBIO = (re.compile(x) for x in (r'gene_id "([^"]+)"', r'transcript_id "([^"]+)"', r'tag "([^"]+)"', r'gene_biotype "([^"]+)"', r'transcript_biotype "([^"]+)"'))


def v_ensembl_gtf(tid, r):
    """Ensembl release-116 gene annotation: the header names GRCh38; every line has nine columns, a known feature type and a strand;
    every sequence name is in the FASTA's .fai and every feature lies inside its sequence; transcripts lie inside their genes, and
    exons, CDS, UTRs and codons inside their transcripts on the same sequence and strand; identifiers have the Ensembl form; and
    the CDS phase column follows from the cumulative CDS length in transcription order."""
    fai = {x[0]: x[1] for x in read_fai(path(FASTA_ID + '_fai'))}
    header, counts, tags_seen = [], collections.Counter(), collections.Counter()
    bad = collections.Counter()
    genes, tx = {}, {}
    cur = dict(id=None, cds=[], tags=set(), strand=None, bio=None, chrom=None)
    phase = collections.Counter()
    off_frame = collections.Counter()
    unexplained = []

    def close():
        if not cur['cds']:
            return
        segs = sorted(cur['cds'], reverse=cur['strand'] == '-')
        cum = 0
        ok = True
        for s_, e_, fr in segs:
            ok &= fr == (3 - cum % 3) % 3
            cum += e_ - s_ + 1
        if 'cds_start_NF' not in cur['tags']:
            phase['phase_checked'] += 1
            phase['phase_bad'] += not ok
        if 'cds_start_NF' in cur['tags'] or 'cds_end_NF' in cur['tags']:
            return
        if cur['bio'] != 'protein_coding':
            off_frame[cur['bio']] += cum % 3 != 0                  # for example protein_coding_LoF: a coding gene whose reference allele carries a loss-of-function change
        elif cur['tags'] & {'mRNA_start_NF', 'mRNA_end_NF'}:
            phase['incomplete_mrna'] += cum % 3 != 0                # the transcript end was not found, so the CDS can stop mid-codon
        elif cur['chrom'] == 'MT':
            phase['mitochondrial'] += cum % 3 != 0                  # mitochondrial stop codons are completed by polyadenylation
        else:
            phase['divisibility_tested'] += 1
            if cum % 3:
                phase['divisibility_off_frame'] += 1
                unexplained.append(cur['id'])

    with opener(path(tid)) as f:
        for line in f:
            if line.startswith('#'):
                header.append(line.strip())
                continue
            c = line.rstrip('\n').split('\t')
            if len(c) != 9:
                bad['columns'] += 1
                continue
            seq, src, feat, s_, e_, _, strand, fr, attr = c
            s_, e_ = int(s_), int(e_)
            counts[feat] += 1
            bad['feature_type'] += feat not in GTF_FEATURES
            bad['strand'] += strand not in ('+', '-')
            if seq not in fai:
                bad['sequence_not_in_fai'] += 1
            else:
                bad['span'] += not (1 <= s_ <= e_ <= fai[seq])
            gid = GTF_GID.search(attr)
            gid = gid.group(1) if gid else None
            bad['gene_id_form'] += not (gid and re.fullmatch(r'ENSGR?\d{11}', gid))
            if feat == 'gene':
                genes[gid] = (seq, s_, e_, strand)
                counts['gene:' + (GTF_BIO.search(attr).group(1) if GTF_BIO.search(attr) else 'none')] += 1
                continue
            tr = GTF_TID.search(attr)
            tr = tr.group(1) if tr else None
            bad['transcript_id_form'] += not (tr and re.fullmatch(r'ENSTR?\d{11}', tr))
            if feat == 'transcript':
                close()
                tags = set(GTF_TAG.findall(attr))
                tags_seen.update(tags)
                cur.update(id=tr, cds=[], tags=tags, strand=strand, chrom=seq, bio=GTF_TBIO.search(attr).group(1) if GTF_TBIO.search(attr) else None)
                tx[tr] = (seq, s_, e_, strand)
                g = genes.get(gid)
                bad['transcript_outside_gene'] += not (g and g[0] == seq and g[3] == strand and g[1] <= s_ and e_ <= g[2])
                continue
            t = tx.get(tr)
            bad['feature_outside_transcript'] += not (t and t[0] == seq and t[3] == strand and t[1] <= s_ and e_ <= t[2])
            if feat == 'CDS':
                bad['cds_frame_missing'] += fr not in ('0', '1', '2')
                if tr == cur['id'] and fr in ('0', '1', '2'):
                    cur['cds'].append((s_, e_, int(fr)))
    close()
    gb = {k.split(':', 1)[1]: v for k, v in counts.items() if k.startswith('gene:')}
    r.check('header names GRCh38', any(h.startswith('#!genome-version GRCh38') for h in header), header)
    r.check('every line has nine columns, a known feature type and a strand', not (bad['columns'] or bad['feature_type'] or bad['strand']), {k: bad[k] for k in ('columns', 'feature_type', 'strand')})
    r.check('every sequence name is in the FASTA index and every feature lies inside its sequence', not (bad['sequence_not_in_fai'] or bad['span']), {k: bad[k] for k in ('sequence_not_in_fai', 'span')})
    r.check('gene and transcript identifiers have the Ensembl form', not (bad['gene_id_form'] or bad['transcript_id_form']), {k: bad[k] for k in ('gene_id_form', 'transcript_id_form')})
    r.check('transcripts lie inside their genes and features inside their transcripts (same sequence and strand)', not (bad['transcript_outside_gene'] or bad['feature_outside_transcript']), {k: bad[k] for k in ('transcript_outside_gene', 'feature_outside_transcript')})
    r.check('every CDS line has a phase of 0, 1 or 2', bad['cds_frame_missing'] == 0, bad['cds_frame_missing'])
    r.check('CDS phase follows from the cumulative CDS length in transcription order (transcripts with a complete 5 prime end)', phase['phase_checked'] > 50000 and phase['phase_bad'] == 0, dict(phase))
    r.check('nuclear protein_coding sequences without an incompleteness tag have a length divisible by 3 (at most 0.01% exceptions)', phase['divisibility_tested'] > 50000 and phase['divisibility_off_frame'] <= 0.0001 * phase['divisibility_tested'],
            dict(tested=phase['divisibility_tested'], off_frame=phase['divisibility_off_frame'], ids=unexplained[:10]))
    r.check('genes, transcripts and MANE Select transcripts are present in plausible numbers', 50000 < counts['gene'] < 100000 and counts['transcript'] > 200000 and 17000 < tags_seen['MANE_Select'] < 21000, dict(genes=counts['gene'], transcripts=counts['transcript'], MANE_Select=tags_seen['MANE_Select']))
    r.stats.update(header=header, feature_counts={k: v for k, v in counts.items() if ':' not in k}, gene_biotypes=dict(collections.Counter(gb).most_common(8)), tags={k: tags_seen[k] for k in ('MANE_Select', 'MANE_Plus_Clinical', 'Ensembl_canonical', 'gencode_basic', 'gencode_primary') if k in tags_seen}, other_biotypes_with_length_not_divisible_by_3={k: v for k, v in off_frame.items() if v}, off_frame_without_incompleteness_tag=unexplained, **dict(phase))


SIFTS_COLS = ['PDB', 'CHAIN', 'SP_PRIMARY', 'RES_BEG', 'RES_END', 'PDB_BEG', 'PDB_END', 'SP_BEG', 'SP_END']
UNIPROT_ACC = re.compile(r'^(?:[OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9](?:[A-Z][A-Z0-9]{2}[0-9]){1,2})$')
PDB_ID = re.compile(r'^[0-9][a-z0-9]{3}$')


def v_sifts(tid, r):
    """SIFTS flat files (PDB chain segments mapped to UniProt): a dated comment line, the nine documented columns, PDB and UniProt
    identifiers of the documented forms, integer residue ranges, and the sanity floors of the PDB itself. Segments of engineered
    constructs can run backwards in UniProt numbering or differ in length between SEQRES and UniProt (internal deletions); those
    are counted, not rejected, and only a negligible number of reversed UniProt ranges is tolerated."""
    n = bad_cols = bad_pdb = bad_acc = bad_int = bad_res = bad_sp = len_diff = crlf = 0
    pdbs, accs = set(), set()
    reversed_examples = []
    with opener(path(tid)) as f:
        first = f.readline().strip()
        hdr = f.readline().rstrip('\r\n').split('\t')
        for line in f:
            n += 1
            crlf += line.endswith('\r\n')                       # pdb_chain_uniprot.tsv.gz uses CRLF line endings, uniprot_segments_observed.tsv.gz uses LF
            c = line.rstrip('\r\n').split('\t')
            if len(c) != 9:
                bad_cols += 1
                continue
            pdb, _ch, acc, rb, re_, _pb, _pe, sb, se = c
            bad_pdb += not PDB_ID.match(pdb)
            bad_acc += not UNIPROT_ACC.match(acc)
            try:
                rb, re_, sb, se = int(rb), int(re_), int(sb), int(se)
            except ValueError:
                bad_int += 1
                continue
            pdbs.add(pdb)
            accs.add(acc)
            bad_res += rb > re_ or min(rb, sb) < 1
            if sb > se:
                bad_sp += 1
                if len(reversed_examples) < 5:
                    reversed_examples.append(c)
            len_diff += (re_ - rb) != (se - sb)
    r.check('first line is a dated comment' + (' naming the PDB and UniProt releases' if tid == 'sifts_pdb_chain_uniprot' else ''),
            re.match(r'^# \d{4}/\d{2}/\d{2} - \d{2}:\d{2}', first) is not None and (tid != 'sifts_pdb_chain_uniprot' or re.search(r'PDB: [\d.]+ \| UniProt: [\d.]+', first) is not None), first)
    r.check('column header is the nine documented columns', hdr == SIFTS_COLS, hdr)
    r.check('every row has nine fields, a PDB identifier and a UniProt accession of the documented form', bad_cols == 0 and bad_pdb == 0 and bad_acc == 0 and bad_int == 0, dict(columns=bad_cols, pdb_id=bad_pdb, accession=bad_acc, integers=bad_int))
    r.check('RES ranges are ordered and positive (SEQRES numbering)', bad_res == 0, bad_res)
    r.check('UniProt ranges are ordered except for at most 0.01% of rows (reversed segments of engineered constructs)', bad_sp <= 0.0001 * n, dict(reversed_rows=bad_sp, rows=n, examples=reversed_examples))
    r.check('the file covers at least 200,000 PDB entries and 60,000 UniProt accessions', len(pdbs) > 200000 and len(accs) > 60000, dict(pdb_entries=len(pdbs), uniprot_accessions=len(accs)))
    r.stats.update(release_line=first, rows=n, pdb_entries=len(pdbs), uniprot_accessions=len(accs), reversed_uniprot_ranges=bad_sp, rows_with_unequal_seqres_and_uniprot_length=len_diff, rows_with_crlf_line_ending=crlf)


LOF_VERDICTS = {'LoF', 'Likely LoF', 'Uncertain', 'Likely not LoF', 'Not LoF'}


def v_lof_curation(tid, r):
    """gnomAD v4.1 LoF curation (incomplete penetrance): GRCh38 variant identifiers, Ensembl gene identifiers with symbols, one of
    the five documented verdicts, and TRUE/FALSE curation flags; the MD5 published for the object equals ours (checked at download)."""
    with open(path(tid), newline='', encoding='utf-8') as f:
        rows = list(csv.DictReader(f))
    cols = list(rows[0])
    flags = [c for c in cols if c.startswith('Flag ') or c.startswith('LoF in ')]
    r.check('columns begin Variant ID, Gene, Verdict and 22 TRUE/FALSE flags follow', cols[:3] == ['Variant ID', 'Gene', 'Verdict'] and len(flags) == len(cols) - 3 == 22, cols[:4])
    r.check('the MD5 equals the value published for the object', MANIFEST[tid].get('md5_verified_against_source') is True)
    r.check('variant identifiers are chrom-pos-ref-alt (GRCh38) and unique', all(re.match(r'^(?:[1-9]|1[0-9]|2[0-2]|X|Y)-\d+-[ACGT]+-[ACGT]+$', x['Variant ID']) for x in rows) and len({x['Variant ID'] for x in rows}) == len(rows), len(rows))
    r.check('genes are ENSG:SYMBOL', all(re.match(r'^ENSG\d{11}:\S+$', x['Gene']) for x in rows))
    r.check('verdicts are the five documented classes', {x['Verdict'] for x in rows} <= LOF_VERDICTS, collections.Counter(x['Verdict'] for x in rows))
    r.check('flags are TRUE or FALSE', all(x[c] in ('TRUE', 'FALSE') for x in rows for c in flags))
    r.check('more than 500 curated variants', len(rows) > 500, len(rows))
    r.stats.update(rows=len(rows), genes=len({x['Gene'] for x in rows}), verdicts=dict(collections.Counter(x['Verdict'] for x in rows)))


PV_COLS = ['Haplotype Name', 'Gene', 'rsID', 'ReferenceSequence', 'Variant Start', 'Variant Stop', 'Reference Allele', 'Variant Allele', 'Type']


def v_pharmvar(tid, r):
    """PharmVar complete-database archive: every member passes its CRC; one folder per gene with RefSeqGene, GRCh37 and GRCh38
    haplotype tables, one VCF per haplotype that has variants, and a haplotype FASTA; table, VCF and FASTA forms; and the exact
    relation between the three (VCF files = table haplotypes without variants, per gene and reference; FASTA records = haplotypes)."""
    z = zipfile.ZipFile(path(tid))
    r.check('every member passes its CRC', z.testzip() is None)
    names = z.namelist()
    ver = os.path.basename(path(tid))[:-len('.zip')]
    r.check('one top-level folder, named after the version', {n.split('/')[0] for n in names} == {ver})
    tsv, vcf, fasta = {}, collections.Counter(), {}
    for n in names:
        p = n.split('/')
        if n.endswith('.haplotypes.tsv') and len(p) == 4:
            tsv[(p[1], p[2])] = n
        elif n.endswith('.vcf') and len(p) == 4:
            vcf[(p[1], p[2])] += 1
        elif n.endswith('.haplotypes.fasta') and len(p) == 3:
            fasta[p[1]] = n
    genes = sorted({g for g, _ in tsv})
    r.check('every gene has a RefSeqGene, a GRCh37 and a GRCh38 table and a haplotype FASTA', len(genes) >= 15 and all(all((g, ref) in tsv for ref in ('RefSeqGene', 'GRCh37', 'GRCh38')) and g in fasta for g in genes), genes)
    bad_head = bad_cols = bad_type = bad_coord = bad_allele = 0
    haps, reference_haps = {}, {}
    n_rows = 0
    types = collections.Counter()
    for key, n in tsv.items():
        lines = z.read(n).decode().split('\n')
        bad_head += not (lines[0] == '#version=' + ver and lines[1].split('\t') == PV_COLS)
        hs, rs = set(), set()
        for l in lines[2:]:
            if not l.strip():
                continue
            c = l.split('\t')
            if len(c) != 9:
                bad_cols += 1
                continue
            n_rows += 1
            hs.add(c[0])
            types[c[8]] += 1
            if c[3] == 'REFERENCE':
                rs.add(c[0])
                continue
            bad_type += c[8] not in ('substitution', 'deletion', 'insertion')
            bad_coord += not (c[4].isdigit() and c[5].isdigit() and int(c[4]) <= int(c[5]))
            bad_allele += bool(set(c[6] + c[7]) - set('ACGTN-'))
        haps[key], reference_haps[key] = hs, rs
    r.check('every table has the version line and the nine documented columns', bad_head == 0, bad_head)
    r.check('every table row has nine fields, a documented type, ordered numeric coordinates and A/C/G/T/N/- alleles', bad_cols == 0 and bad_type == 0 and bad_coord == 0 and bad_allele == 0,
            dict(columns=bad_cols, type=bad_type, coordinates=bad_coord, alleles=bad_allele))
    off = [(k, len(haps[k]) - len(reference_haps[k]), vcf[k]) for k in tsv if len(haps[k]) - len(reference_haps[k]) != vcf[k]]
    r.check('VCF files per gene and reference equal the table haplotypes that have variants', not off, off[:5])
    off_fa = []
    for g in genes:
        recs = [l for l in z.read(fasta[g]).decode().split('\n') if l.startswith('>')]
        if len(recs) != len(haps[(g, 'RefSeqGene')]) or not all(('PharmVar Version:' + ver.split('-', 1)[1]) in l for l in recs):
            off_fa.append((g, len(recs), len(haps[(g, 'RefSeqGene')])))
    r.check('FASTA records equal the haplotypes of the RefSeqGene table and name this version', not off_fa, off_fa[:5])
    bad_vcf = no_var = 0
    n_var = 0
    for n in names:
        if n.endswith('.vcf'):
            lines = z.read(n).decode().split('\n')
            data = [l for l in lines if l and not l.startswith('#')]
            hdr = [l for l in lines if l.startswith('#CHROM')]
            ok = lines[0].startswith('##fileformat=VCFv4.') and len(hdr) == 1 and hdr[0].split('\t') == ['#CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO']
            for l in data:
                c = l.split('\t')
                ok &= len(c) == 8 and c[1].isdigit() and set(c[3]) <= set('ACGTN') and set(c[4]) <= set('ACGTN,')
            bad_vcf += not ok
            no_var += not data
            n_var += len(data)
    r.check('every VCF has a v4 header, the eight columns and well-formed records, and at least one variant', bad_vcf == 0 and no_var == 0, dict(malformed=bad_vcf, without_variants=no_var))
    r.stats.update(genes=len(genes), haplotypes_refseqgene={g: len(haps[(g, 'RefSeqGene')]) for g in genes}, table_rows=n_rows, row_types=dict(types), vcf_files=sum(vcf.values()), vcf_records=n_var, members=len(names))


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


# --- RefSeq GRCh38.p14, GIAB v3.6 stratifications, ProteinGym v1.3 ----------------------------------------------------
RS_REPORT = 'refseq_grch38p14_assembly_report'
RS_HEADER = ['Sequence-Name', 'Sequence-Role', 'Assigned-Molecule', 'Assigned-Molecule-Location/Type', 'GenBank-Accn', 'Relationship', 'RefSeq-Accn', 'Assembly-Unit', 'Sequence-Length', 'UCSC-style-name']
RS_ACC = re.compile(r'^(NC|NT|NW)_\d+\.\d+$')
PRIMARY_CHROMS = [str(i) for i in range(1, 23)] + ['X', 'Y']


def refseq_report(p):
    """(metadata, column names, rows) of an NCBI assembly report; the UCSC-style names and RefSeq accessions are the join keys used elsewhere."""
    meta, hdr, rows = {}, None, []
    with open(p, encoding='utf-8') as f:
        for l in f:
            l = l.rstrip('\n')
            if l.startswith('# Sequence-Name'):
                hdr = l[2:].split('\t')
            elif l.startswith('#'):
                if not l.startswith('##') and ':' in l:
                    k, _, v = l[1:].partition(':')
                    meta[k.strip()] = v.strip()
            elif l:
                rows.append(l.split('\t'))
    return meta, hdr, rows


def ucsc_lengths():
    """UCSC-style sequence name -> length, from the RefSeq assembly report (the bridge between GIAB's chr1 naming and the RefSeq and Ensembl naming)."""
    meta, hdr, rows = refseq_report(path(RS_REPORT))
    return {x[9]: int(x[8]) for x in rows if x[9] != 'na'}


def commented_table(p, header_start):
    """Column names from the commented header line starting with header_start, and the data rows of a '#'-commented TSV."""
    hdr, rows = None, []
    with open(p, encoding='utf-8') as f:
        for l in f:
            l = l.rstrip('\n')
            if l.startswith(header_start):
                hdr = l.lstrip('#').strip().split('\t')
            elif l and not l.startswith('#'):
                rows.append(l.split('\t'))
    return hdr, rows


def v_refseq(tid, r):
    p = path(tid)
    if tid == RS_REPORT:
        meta, hdr, rows = refseq_report(p)
        r.check('column header', hdr == RS_HEADER, hdr)
        r.check('GRCh38.p14, GCF_000001405.40', meta.get('Assembly name') == 'GRCh38.p14' and meta.get('RefSeq assembly accession') == 'GCF_000001405.40', {k: meta.get(k) for k in ('Assembly name', 'RefSeq assembly accession')})
        r.check('every row has 10 columns', all(len(x) == 10 for x in rows), len(rows))
        roles = collections.Counter(x[1] for x in rows)
        r.check('sequence roles are the known set', set(roles) <= {'assembled-molecule', 'alt-scaffold', 'fix-patch', 'novel-patch', 'unlocalized-scaffold', 'unplaced-scaffold'}, dict(roles))
        chrom = {x[0]: x for x in rows if x[1] == 'assembled-molecule' and x[3] == 'Chromosome'}
        r.check('chromosomes 1-22, X and Y are assembled molecules', set(PRIMARY_CHROMS) <= set(chrom), sorted(chrom))
        c1 = chrom.get('1', [''] * 10)
        r.check('chromosome 1 is NC_000001.11, 248,956,422 bases, named chr1', c1[6] == 'NC_000001.11' and c1[8] == '248956422' and c1[9] == 'chr1', c1)
        have = [x for x in rows if x[6] != 'na']
        r.check('RefSeq accessions have NC_, NT_ or NW_ form (a few GenBank-only sequences have none)', all(RS_ACC.match(x[6]) for x in have) and len(rows) - len(have) <= 10, [x[6] for x in have if not RS_ACC.match(x[6])][:5])
        r.check('RefSeq accessions are unique', len({x[6] for x in have}) == len(have), len(have))
        r.check('lengths are positive integers', all(x[8].isdigit() and int(x[8]) > 0 for x in rows))
        mt = [x for x in rows if x[3] == 'Mitochondrion']
        r.check('mitochondrial genome NC_012920.1 (16,569 bases) present', any(x[6] == 'NC_012920.1' and x[8] == '16569' for x in mt), [x[6] for x in mt])
        r.stats.update(rows=len(rows), roles=dict(roles), primary_chromosomes=len(chrom), mitochondrion=len(mt), without_refseq_accession=len(rows) - len(have))
    elif tid == 'refseq_grch38p14_assembly_regions':
        hdr, rows = commented_table(p, '# Region-Name')
        r.check('column header', hdr == ['Region-Name', 'Chromosome', 'Chromosome-Start', 'Chromosome-Stop', 'Scaffold-Role', 'Scaffold-GenBank-Accn', 'Scaffold-RefSeq-Accn', 'Assembly-Unit'], hdr)
        r.check('every row has 8 columns', all(len(x) == 8 for x in rows), len(rows))
        names = {x[0] for x in refseq_report(path(RS_REPORT))[2]}
        offchrom = [x[1] for x in rows if x[1] not in PRIMARY_CHROMS]
        r.check('chromosome is 1-22, X or Y, or (alternate loci placed on another alternate scaffold) a sequence of the assembly report', all(c in names for c in offchrom), offchrom[:5])
        r.check('start < stop', all(x[2].isdigit() and x[3].isdigit() and int(x[2]) < int(x[3]) for x in rows))
        roles = collections.Counter(x[4] for x in rows)
        r.check('roles are alt-scaffold, fix-patch, novel-patch, or the unplaced region types CEN, HET, PAR and na', set(roles) <= {'alt-scaffold', 'fix-patch', 'novel-patch', 'CEN', 'HET', 'PAR', 'na'}, dict(roles))
        sc = [x for x in rows if x[4] in ('alt-scaffold', 'fix-patch', 'novel-patch')]
        r.check('scaffold RefSeq accessions have NT_ or NW_ form for scaffold regions and are na for the others', all(re.match(r'^N[TW]_\d+\.\d+$', x[6]) for x in sc) and all(x[6] == 'na' for x in rows if x not in sc), len(sc))
        r.check('region rows > 400', len(rows) > 400, len(rows))
        r.stats.update(rows=len(rows), regions=len({x[0] for x in rows}), roles=dict(roles), scaffold_rows=len(sc))
    elif tid == 'refseq_grch38p14_alt_scaffold_placement':
        with open(p, encoding='utf-8') as f:
            rows = [l.rstrip('\n').split('\t') for l in f if l.strip()]
        hdr, rows = rows[0], rows[1:]
        want = ['#alt_asm_name', 'prim_asm_name', 'alt_scaf_name', 'alt_scaf_acc', 'parent_type', 'parent_name', 'parent_acc', 'region_name', 'ori', 'alt_scaf_start', 'alt_scaf_stop', 'parent_start', 'parent_stop', 'alt_start_tail', 'alt_stop_tail']
        r.check('column header', hdr == want, hdr)
        r.check('every row has 15 columns', all(len(x) == 15 for x in rows), len(rows))
        r.check('parent types are CHROMOSOME or SCAFFOLD', {x[4] for x in rows} <= {'CHROMOSOME', 'SCAFFOLD'}, collections.Counter(x[4] for x in rows))
        r.check('scaffold accessions have NT_ or NW_ form, parent accessions NC_, NT_ or NW_', all(re.match(r'^N[TW]_\d+\.\d+$', x[3]) and RS_ACC.match(x[6]) for x in rows))
        r.check('coordinates are integers and start <= stop', all(x[9].isdigit() and x[10].isdigit() and x[11].isdigit() and x[12].isdigit() and int(x[9]) <= int(x[10]) and int(x[11]) <= int(x[12]) for x in rows))
        r.check("orientation is '+', '-' or 'b' (mixed), as defined in the assembly_structure README", {x[8] for x in rows} <= {'+', '-', 'b'}, collections.Counter(x[8] for x in rows))
        r.check('scaffold accessions are unique', len({x[3] for x in rows}) == len(rows), len(rows))
        r.stats.update(rows=len(rows), assembly_units=dict(collections.Counter(x[0] for x in rows)), regions=len({x[7] for x in rows}), orientations=dict(collections.Counter(x[8] for x in rows)))
    elif tid == 'refseq_grch38p14_feature_count':
        hdr, rows = commented_table(p, '# Feature')
        r.check('column header', hdr == ['Feature', 'Class', 'Full Assembly', 'Assembly-unit accession', 'Assembly-unit name', 'Unique Ids', 'Placements'], hdr)
        r.check('every row has 7 columns', all(len(x) == 7 for x in rows), len(rows))
        pc = next((x for x in rows if x[0] == 'gene' and x[1] == 'protein_coding' and x[4] == 'all'), None)
        r.check('protein-coding genes (whole assembly) between 19,000 and 21,000', pc is not None and 19000 <= int(pc[5]) <= 21000, pc)
        r.stats.update(rows=len(rows), protein_coding_genes=int(pc[5]) if pc else None)
    elif tid == 'refseq_annotation_release_readme':
        kv = {}
        with open(p, encoding='utf-8') as f:
            for l in f:
                if ':\t' in l:
                    k, _, v = l.rstrip('\n').partition(':\t')
                    kv[k] = v
        r.check('annotation release GCF_000001405.40-RS_2025_08', kv.get('ANNOTATION RELEASE NAME') == 'GCF_000001405.40-RS_2025_08', kv.get('ANNOTATION RELEASE NAME'))
        r.check('assembly GRCh38.p14, GCF_000001405.40', kv.get('ASSEMBLY NAME') == 'GRCh38.p14' and kv.get('ASSEMBLY ACCESSION') == 'GCF_000001405.40')
        r.check('provider NCBI RefSeq, taxid 9606', kv.get('ANNOTATION PROVIDER') == 'NCBI RefSeq' and kv.get('TAXID') == '9606')
        r.stats.update({k.lower().replace(' ', '_'): kv.get(k) for k in ('ANNOTATION SOFTWARE VERSION', 'ANNOTATION EVIDENCE FREEZE DATE', 'ANNOTATION RELEASE DATE', 'NUMBER OF ASSEMBLY-UNITS')})
    elif tid == 'refseq_grch38p14_genomic_gff':
        meta, hdr, arows = refseq_report(path(RS_REPORT))
        accs = {x[6] for x in arows}
        head, types, unknown, bad, strands, nlines = [], collections.Counter(), collections.Counter(), 0, collections.Counter(), 0
        ids, parents, genes_hgnc, geneids = set(), set(), 0, set()
        with gzip.open(p, 'rt', encoding='utf-8') as f:
            for l in f:
                if l.startswith('#'):
                    head.append(l.rstrip('\n'))
                    continue
                c = l.rstrip('\n').split('\t')
                nlines += 1
                if len(c) != 9 or not c[3].isdigit() or not c[4].isdigit() or int(c[3]) > int(c[4]) or c[6] not in '+-.':
                    bad += 1
                    continue
                if c[0] not in accs:
                    unknown[c[0]] += 1
                types[c[2]] += 1
                strands[c[6]] += 1
                a = c[8]
                if c[2] in ('exon', 'CDS'):
                    m = re.search(r'Parent=([^;]+)', a)
                    if m:
                        parents.update(m.group(1).split(','))
                else:
                    m = re.search(r'(?:^|;)ID=([^;]+)', a)
                    if m:
                        ids.add(m.group(1))
                    if c[2] in ('mRNA', 'transcript', 'lnc_RNA'):
                        m = re.search(r'Parent=([^;]+)', a)
                        if m:
                            parents.update(m.group(1).split(','))
                    if c[2] == 'gene':
                        m = re.search(r'Dbxref=([^;]+)', a)
                        if m:
                            for x in m.group(1).split(','):
                                if x.startswith('GeneID:'):
                                    geneids.add(x)
                            genes_hgnc += 'HGNC:HGNC:' in m.group(1)
        r.check('GFF3 header names GRCh38.p14, GCF_000001405.40 and annotation release RS_2025_08', '##gff-version 3' in head and '#!genome-build GRCh38.p14' in head and '#!genome-build-accession NCBI_Assembly:GCF_000001405.40' in head and '#!annotation-source NCBI RefSeq GCF_000001405.40-RS_2025_08' in head, head[:6])
        r.check('every feature line has 9 columns, integer start <= end and a strand', bad == 0, bad)
        r.check('every sequence is in the assembly report', not unknown, dict(unknown.most_common(3)))
        r.check('every Parent refers to a defined feature ID', parents <= ids, sorted(parents - ids)[:5])
        fc = commented_table(path('refseq_grch38p14_feature_count'), '# Feature')[1]
        placed = sum(int(x[6]) for x in fc if x[0] == 'gene' and x[4] == 'all')
        r.check('gene and pseudogene features equal the gene placements of the feature-count file', types['gene'] + types['pseudogene'] == placed, (types['gene'] + types['pseudogene'], placed))
        r.check('mRNA > 100,000, exon > 1,000,000, CDS > 1,000,000', types['mRNA'] > 100000 and types['exon'] > 1000000 and types['CDS'] > 1000000, {k: types[k] for k in ('mRNA', 'exon', 'CDS')})
        r.stats.update(lines=nlines, feature_types=dict(types.most_common(12)), genes_with_geneid=len(geneids), genes_with_hgnc=genes_hgnc, strands=dict(strands))


def v_giab(tid, r):
    p = path(tid)
    lst = 'giab_v3_6_grch38_stratification_list'
    if tid == 'giab_v3_6_grch38_readme':
        t = open(p, encoding='utf-8').read()
        r.check('describes the GRCh38 stratifications', 'stratifications for the reference GRCh38' in t)
        r.check('has the Segmental Duplications and Other Difficult sections', '## Segmental Duplications' in t and '## Other Difficult' in t)
        r.stats['characters'] = len(t)
    elif tid == 'giab_v3_6_grch38_background':
        t = open(p, encoding='utf-8').read()
        r.check('states the reuse condition', 'freely available for reuse without embargo' in t)
        r.check('names the GIAB consortium', 'Genome in a Bottle' in t or 'GIAB' in t)
        r.stats['characters'] = len(t)
    elif tid == lst:
        rows = [l.rstrip('\n').split('\t') for l in open(p, encoding='utf-8') if l.strip()]
        r.check('two columns per row', all(len(x) == 2 for x in rows), len(rows))
        r.check('paths are category/GRCh38_*.bed.gz', all(re.match(r'^[A-Za-z]+/GRCh38_.+\.bed\.gz$', x[1]) for x in rows))
        r.check('names are unique', len({x[0] for x in rows}) == len(rows))
        have = {x[1] for x in rows}
        mine = [TARGETS[t]['url'].split('GRCh38@all/')[1] for t in TARGETS if t.startswith('giab_v3_6_grch38_') and TARGETS[t]['kind'] == 'bed.gz']
        r.check('every downloaded stratification is in the list', all(m in have for m in mine), [m for m in mine if m not in have][:5])
        r.check('list has > 150 stratifications', len(rows) > 150, len(rows))
        r.stats.update(rows=len(rows), categories=dict(collections.Counter(x[1].split('/')[0] for x in rows)), downloaded=len(mine))
    elif tid == 'giab_v3_6_grch38_md5s':
        rows = [l.rstrip('\n').split('  ', 1) for l in open(p, encoding='utf-8') if l.strip()]
        r.check('md5 and path per line', all(len(x) == 2 and re.match(r'^[0-9a-f]{32}$', x[0]) for x in rows), len(rows))
        have = {x[1]: x[0] for x in rows}
        mine = {TARGETS[t]['url'].split('GRCh38@all/')[1]: TARGETS[t].get('expected_md5') for t in TARGETS if t.startswith('giab_v3_6_grch38_') and TARGETS[t]['kind'] == 'bed.gz'}
        r.check('the md5 pinned for each downloaded stratification equals the published list', all(have.get(k) == v for k, v in mine.items()), [k for k, v in mine.items() if have.get(k) != v][:5])
        r.stats.update(rows=len(rows))
    else:
        lens = ucsc_lengths()
        n = bases = bad = past = 0
        chroms, overlap, unsorted, last = collections.Counter(), 0, 0, {}
        with gzip.open(p, 'rt', encoding='utf-8') as f:
            for l in f:
                if not l.strip() or l.startswith(('#', 'track', 'browser')):
                    continue
                c = l.rstrip('\n').split('\t')
                n += 1
                if len(c) < 3 or not c[1].isdigit() or not c[2].isdigit() or int(c[1]) >= int(c[2]):
                    bad += 1
                    continue
                s, e = int(c[1]), int(c[2])
                if c[0] not in lens or e > lens[c[0]]:
                    past += 1
                chroms[c[0]] += 1
                bases += e - s
                if c[0] in last:
                    unsorted += s < last[c[0]][0]
                    overlap += s < last[c[0]][1]
                last[c[0]] = (s, e)
        r.check('every line is a BED interval with start < end', bad == 0, bad)
        r.check('rows > 0', n > 0, n)
        r.check('chromosome names are UCSC-style names of the assembly report and every interval ends inside its sequence', past == 0, past)
        r.check('intervals are sorted by start within each chromosome', unsorted == 0, unsorted)
        r.stats.update(rows=n, bases=bases, chromosomes=len(chroms), overlapping_neighbours=overlap)


PG_COLS = ['mutant', 'mutated_sequence', 'DMS_score', 'DMS_score_bin']
PG_TOKEN = re.compile(r'^([A-Z])(\d+)([A-Z])$')


def v_proteingym(tid, r):
    if tid.endswith('_reference'):
        sub = 'indels' not in tid
        rows = list(csv.DictReader(open(path(tid), encoding='utf-8')))
        need = ['DMS_id', 'DMS_filename', 'UniProt_ID', 'taxon', 'target_seq', 'seq_len', 'DMS_total_number_mutants', 'DMS_binarization_cutoff', 'first_author', 'year', 'selection_type']
        if sub:
            need += ['includes_multiple_mutants', 'DMS_number_single_mutants', 'DMS_number_multiple_mutants']
        r.check('required columns present', all(c in rows[0] for c in need), [c for c in need if c not in rows[0]])
        r.check(f'assays equal the {217 if sub else 66} on the ProteinGym website', len(rows) == (217 if sub else 66), len(rows))
        r.check('DMS_id and DMS_filename are unique', len({x['DMS_id'] for x in rows}) == len(rows) and len({x['DMS_filename'] for x in rows}) == len(rows))
        r.check('target sequences are amino-acid strings whose length is seq_len', all(re.match(r'^[ACDEFGHIKLMNPQRSTVWYX]+$', x['target_seq']) and len(x['target_seq']) == int(x['seq_len']) for x in rows))
        if sub:
            r.check('single + multiple mutants = total, and includes_multiple_mutants agrees', all(int(x['DMS_number_single_mutants']) + int(x['DMS_number_multiple_mutants']) == int(x['DMS_total_number_mutants']) and (x['includes_multiple_mutants'].lower() == 'true') == (int(x['DMS_number_multiple_mutants']) > 0) for x in rows))
        total = sum(int(x['DMS_total_number_mutants']) for x in rows)
        r.stats.update(assays=len(rows), mutants=total, taxa=dict(collections.Counter(x['taxon'] for x in rows)), selection_types=dict(collections.Counter(x['selection_type'] for x in rows).most_common(8)))
        return
    sub = 'indels' not in tid
    ref = list(csv.DictReader(open(path(tid + '_reference'), encoding='utf-8')))
    prefix = 'DMS_ProteinGym_substitutions/' if sub else 'DMS_ProteinGym_indels/'
    with zipfile.ZipFile(path(tid)) as z:
        r.check('zip CRC test', z.testzip() is None)
        members = {n for n in z.namelist() if n.endswith('.csv')}
        r.check('one CSV per assay of the reference file', members == {prefix + x['DMS_filename'] for x in ref}, len(members))
        counts, hdrs, wt_bad, seq_bad, fmt_bad, nonnum, binbad, rows_total, extra = 0, collections.Counter(), 0, 0, 0, 0, 0, 0, []
        for x in ref:
            rd = csv.DictReader(io.StringIO(z.read(prefix + x['DMS_filename']).decode('utf-8')))
            hdrs[tuple(rd.fieldnames)] += 1
            n, ts = 0, x['target_seq']
            for row in rd:
                n += 1
                try:
                    float(row['DMS_score'])
                except (ValueError, KeyError):
                    nonnum += 1
                binbad += row.get('DMS_score_bin') not in ('0', '1')
                if sub:
                    s, ok = list(ts), True
                    for tk in row['mutant'].split(':'):
                        g = PG_TOKEN.match(tk)
                        if not g:
                            fmt_bad += 1
                            ok = False
                            break
                        pos = int(g.group(2))
                        if pos < 1 or pos > len(s) or s[pos - 1] != g.group(1):
                            wt_bad += 1
                            ok = False
                            break
                        s[pos - 1] = g.group(3)
                    if ok and ''.join(s) != row['mutated_sequence']:
                        seq_bad += 1
                elif not re.match(r'^[A-Z]+$', row['mutated_sequence']):
                    fmt_bad += 1
            rows_total += n
            if n != int(x['DMS_total_number_mutants']):
                counts += 1
        base_cols = set(PG_COLS if sub else PG_COLS[1:])
        r.check('every assay file has the benchmark columns (extra columns are listed in the stats)', all(base_cols <= set(h) for h in hdrs), [h for h in hdrs if not base_cols <= set(h)])
        r.check('rows per assay equal DMS_total_number_mutants of the reference file', counts == 0, counts)
        r.check('DMS_score is numeric and DMS_score_bin is 0 or 1 in every row', nonnum == 0 and binbad == 0, (nonnum, binbad))
        if sub:
            r.check('every mutation token is wild-type, position, mutant', fmt_bad == 0, fmt_bad)
            r.check('wild-type residue equals the target sequence at that position in every row', wt_bad == 0, wt_bad)
            r.check('applying the mutations to the target sequence gives mutated_sequence in every row', seq_bad == 0, seq_bad)
        else:
            r.check('mutated_sequence is an amino-acid string in every row', fmt_bad == 0, fmt_bad)
        r.stats.update(assays=len(ref), rows=rows_total, header_variants=[list(h) for h in hdrs], wild_type_mismatches=wt_bad, sequence_mismatches=seq_bad)


# --- HPA v25.1, BioGRID interactions, HGSVC3 SV calls, CPIC API tables --------------------------------------------------
ENSG_RE = re.compile(r'^ENSG\d{11}$')
BIOGRID_HEADER = ['#BioGRID Interaction ID', 'Entrez Gene Interactor A', 'Entrez Gene Interactor B', 'BioGRID ID Interactor A', 'BioGRID ID Interactor B', 'Systematic Name Interactor A', 'Systematic Name Interactor B',
                  'Official Symbol Interactor A', 'Official Symbol Interactor B', 'Synonyms Interactor A', 'Synonyms Interactor B', 'Experimental System', 'Experimental System Type', 'Author', 'Publication Source',
                  'Organism ID Interactor A', 'Organism ID Interactor B', 'Throughput', 'Score', 'Modification', 'Qualifications', 'Tags', 'Source Database', 'SWISS-PROT Accessions Interactor A',
                  'TREMBL Accessions Interactor A', 'REFSEQ Accessions Interactor A', 'SWISS-PROT Accessions Interactor B', 'TREMBL Accessions Interactor B', 'REFSEQ Accessions Interactor B',
                  'Ontology Term IDs', 'Ontology Term Names', 'Ontology Term Categories', 'Ontology Term Qualifier IDs', 'Ontology Term Qualifier Names', 'Ontology Term Types', 'Organism Name Interactor A', 'Organism Name Interactor B']


def v_hpa(tid, r):
    names = v_zip(tid, r)
    r.check('one TSV member', len(names) == 1 and names[0].endswith('.tsv'), names)
    with zipfile.ZipFile(path(tid)) as z, z.open(names[0]) as f:
        t = io.TextIOWrapper(f, encoding='utf-8')
        hdr = t.readline().rstrip('\n').split('\t')
        n, bad_gene, bad_num, genes, groups, pairs = 0, 0, 0, set(), set(), set()
        reliab = collections.Counter()
        for l in t:
            p = l.rstrip('\n').split('\t')
            n += 1
            genes.add(p[0])
            bad_gene += not ENSG_RE.match(p[0])
            if 'consensus' in tid or 'single_cell' in tid:
                groups.add(p[2])
                pairs.add((p[0], p[2]))
                try:
                    bad_num += float(p[3]) < 0
                except ValueError:
                    bad_num += 1
            else:
                reliab[p[2]] += 1
    if 'consensus' in tid:
        r.check('columns Gene, Gene name, Tissue, nTPM', hdr == ['Gene', 'Gene name', 'Tissue', 'nTPM'], hdr)
        r.check('20,162 genes and 51 consensus tissues, as stated on the HPA download page', len(genes) == 20162 and len(groups) == 51, (len(genes), len(groups)))
        r.check('no repeated gene-tissue pair, and at least 99% of all gene-tissue pairs present', len(pairs) == n and n >= 0.99 * len(genes) * len(groups), (n, len(pairs), len(genes) * len(groups)))
    elif 'single_cell' in tid:
        r.check('columns Gene, Gene name, Cell type, nCPM', hdr == ['Gene', 'Gene name', 'Cell type', 'nCPM'], hdr)
        r.check('20,162 genes, as stated on the HPA download page, and more than 100 cell types', len(genes) == 20162 and len(groups) > 100, (len(genes), len(groups)))
        r.check('no repeated gene-cell-type pair, and at least 99% of all gene-cell-type pairs present', len(pairs) == n and n >= 0.99 * len(genes) * len(groups), (n, len(pairs), len(genes) * len(groups)))
    else:
        r.check('columns start Gene, Gene name, Reliability, Main location', hdr[:4] == ['Gene', 'Gene name', 'Reliability', 'Main location'], hdr[:5])
        r.check('13,603 genes, as stated on the HPA download page', len(genes) == 13603 and n == 13603, (len(genes), n))
        r.check('reliability is Enhanced, Supported, Approved or Uncertain', set(reliab) <= {'Enhanced', 'Supported', 'Approved', 'Uncertain'}, dict(reliab))
    r.check('Gene is an Ensembl gene id in every row', bad_gene == 0, bad_gene)
    r.check('expression values are non-negative numbers', bad_num == 0, bad_num)
    r.stats.update(rows=n, genes=len(genes), groups=len(groups), reliability=dict(reliab) or None, grid_fraction=round(n / (len(genes) * len(groups)), 4) if groups else None)


def v_biogrid(tid, r):
    names = v_zip(tid, r)
    versions = {m.group(1) for m in (re.match(r'^BIOGRID-ORGANISM-.+-(\d+\.\d+\.\d+)\.tab3\.txt$', n) for n in names) if m}
    r.check('every member is BIOGRID-ORGANISM-<organism>-<version>.tab3.txt with one release number', len(versions) == 1 and all(n.endswith('.tab3.txt') for n in names), (sorted(versions), len(names)))
    hs = [n for n in names if 'Homo_sapiens' in n]
    r.check('a Homo sapiens member is present', len(hs) == 1, hs)
    with zipfile.ZipFile(path(tid)) as z, z.open(hs[0]) as f:
        t = io.TextIOWrapper(f, encoding='utf-8')
        hdr = t.readline().rstrip('\n').split('\t')
        r.check('37 documented columns', hdr == BIOGRID_HEADER, hdr[:5])
        i = {c: hdr.index(c) for c in ('#BioGRID Interaction ID', 'Experimental System Type', 'Organism ID Interactor A', 'Organism ID Interactor B', 'Throughput', 'Publication Source', 'Official Symbol Interactor A', 'Experimental System')}
        n = both = onehuman = 0
        partner = collections.Counter()
        types, system, thr, ids, pubs, short = collections.Counter(), collections.Counter(), collections.Counter(), set(), set(), 0
        for l in t:
            p = l.rstrip('\n').split('\t')
            if len(p) != len(hdr):
                short += 1
                continue
            n += 1
            ids.add(p[i['#BioGRID Interaction ID']])
            a, b = p[i['Organism ID Interactor A']], p[i['Organism ID Interactor B']]
            both += a == b == '9606'
            onehuman += '9606' in (a, b)
            if a != b:
                partner[b if a == '9606' else a] += 1
            types[p[i['Experimental System Type']]] += 1
            system[p[i['Experimental System']]] += 1
            thr[p[i['Throughput']]] += 1
            pubs.add(p[i['Publication Source']])
    r.check('every row has 37 fields', short == 0, short)
    r.check('interaction ids are unique', len(ids) == n, n)
    r.check('every interaction has a human interactor, and at least 90% are human-human (the rest are host-pathogen and cross-species interactions)', onehuman == n and both / max(1, n) >= 0.90, (onehuman, both, n))
    r.check('system type is physical or genetic', set(types) <= {'physical', 'genetic'}, dict(types))
    r.check('more than 1,000,000 interactions and 15,000 genetic interactions', n > 1_000_000 and types['genetic'] > 15000, (n, dict(types)))
    r.stats.update(release=sorted(versions)[0] if versions else None, organisms_in_archive=len(names), human_interactions=n, human_human=both, cross_species=n - both, cross_species_partner_taxa=dict(partner.most_common(6)), system_types=dict(types), throughput=dict(thr), publications=len(pubs), experimental_systems=len(system),
                   genetic_systems=[k for k, _ in system.most_common(40) if k in ('Synthetic Lethality', 'Synthetic Rescue', 'Dosage Rescue', 'Dosage Lethality', 'Negative Genetic', 'Positive Genetic', 'Synthetic Growth Defect', 'Phenotypic Enhancement', 'Phenotypic Suppression')])


def v_hgsvc3(tid, r):
    p = path(tid)
    if tid == 'hgsvc3_variant_calls_readme':
        t = open(p, encoding='utf-8').read()
        r.check('release 1.0, PAV, 65 phased genomes', '# Release 1.0 variant callset' in t and 'PAV' in t and '65 phased genomes' in t)
        r.stats['characters'] = len(t)
    elif tid == 'hgsvc3_variant_calls_manifest':
        rows = [l.split() for l in open(p, encoding='utf-8') if l.strip()]
        r.check('path, size and MD5 per line', all(len(x) == 3 and x[1].isdigit() and re.match(r'^[0-9a-f]{32}$', x[2]) for x in rows), len(rows))
        have = {x[0].lstrip('./'): (int(x[1]), x[2]) for x in rows}
        mine = {TARGETS[t]['url'].split('Variant_Calls/1.0/')[1]: (TARGETS[t]['expected_size'], TARGETS[t]['expected_md5']) for t in TARGETS if t.startswith('hgsvc3_') and TARGETS[t].get('expected_md5')}
        r.check('size and MD5 of the downloaded VCF files equal their manifest entries', all(have.get(k) == v for k, v in mine.items()), [k for k, v in mine.items() if have.get(k) != v])
        r.stats.update(rows=len(rows))
    else:
        want = 'INV' if '_inv_' in tid else None
        lens = ucsc_lengths()
        samples, n, bad_id, bad_len, bad_gt, past, bad_type, id_pos0 = 0, 0, 0, 0, 0, 0, 0, 0
        types, info_keys, contigs, ids = collections.Counter(), collections.Counter(), {}, set()
        hdr = None
        with gzip.open(p, 'rt', encoding='utf-8') as f:
            for l in f:
                if l.startswith('##contig'):
                    m = re.match(r'##contig=<ID=([^,]+),length=(\d+)', l)
                    contigs[m.group(1)] = int(m.group(2))
                elif l.startswith('#CHROM'):
                    hdr = l.rstrip('\n').split('\t')
                    samples = len(hdr) - 9
                elif l[0] != '#':
                    c = l.rstrip('\n').split('\t', 9)
                    n += 1
                    info = dict(kv.split('=', 1) if '=' in kv else (kv, True) for kv in c[7].split(';'))
                    info_keys.update(info.keys())
                    st = info.get('SVTYPE')
                    types[st] += 1
                    bad_type += st not in (('INV',) if want else ('INS', 'DEL'))
                    try:
                        sl = int(info['SVLEN'])
                    except (KeyError, ValueError):
                        sl = None
                    bad_len += sl is None or (st == 'DEL' and sl >= 0) or (st in ('INS', 'INV') and sl <= 0)
                    base = re.sub(r'\.\d+$', '', c[2])           # a repeated id carries a numeric suffix such as .1
                    pos1 = base == f"{c[0]}-{int(c[1]) + 1}-{st}-{abs(sl) if sl is not None else '?'}"
                    pos0 = base == f"{c[0]}-{int(c[1])}-{st}-{abs(sl) if sl is not None else '?'}"
                    id_pos0 += pos0 and not pos1
                    bad_id += not (pos1 or pos0) or c[2] != info.get('ID')
                    ids.add(c[2])
                    past += c[0] not in lens or int(c[1]) > lens[c[0]]
                    gts = c[9].split('\t', 1)[0] if len(c) > 9 else ''
                    bad_gt += not re.match(r'^[01.]\|[01.]$', gts)
        r.check('VCF header with 65 sample columns, as the release README states', hdr is not None and hdr[:9] == ['#CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT'] and samples == 65, samples)
        r.check('contig lengths in the header equal the RefSeq report lengths for chromosomes 1-22, X, Y', all(contigs.get(k) == v for k, v in lens.items() if k in contigs and re.match(r'^chr([0-9]+|X|Y)$', k)) and sum(1 for k in contigs if re.match(r'^chr([0-9]+|X|Y)$', k)) == 24, len(contigs))
        r.check('SVTYPE is INS or DEL' if not want else 'SVTYPE is INV', bad_type == 0, dict(types))
        r.check('SVLEN is present, negative for DEL and positive for INS' if not want else 'SVLEN is present and positive', bad_len == 0, bad_len)
        r.check('record id equals chromosome-position-SVTYPE-|SVLEN| (position POS+1, or POS for a few insertions; a repeated id has a numeric suffix) and the ID field of INFO', bad_id == 0, bad_id)
        r.check('ids are unique', len(ids) == n, n)
        r.check('positions lie inside the chromosome', past == 0, past)
        r.check('first genotype is a phased pair of 0, 1 or .', bad_gt == 0, bad_gt)
        r.check('records > 100,000' if not want else 'records > 100', n > (100_000 if not want else 100), n)
        r.stats.update(records=n, samples=samples, svtypes={str(k): v for k, v in types.items()}, info_keys=dict(info_keys.most_common(12)), id_position_equals_pos=id_pos0)


CPIC_KEYS = {'allele_definition': ['id', 'genesymbol', 'name', 'pharmvarid', 'matchesreferencesequence', 'structuralvariation'], 'allele_location_value': ['alleledefinitionid', 'locationid', 'variantallele'],
             'sequence_location': ['id', 'name', 'chromosomelocation', 'proteinlocation', 'genesymbol', 'dbsnpid', 'position'], 'allele': ['id', 'genesymbol', 'name', 'functionalstatus', 'clinicalfunctionalstatus', 'activityvalue', 'definitionid', 'citations', 'strength'],
             'gene_result': ['id', 'genesymbol', 'result', 'activityscore'], 'gene_result_lookup': ['id', 'phenotypeid', 'lookupkey', 'function1', 'function2'],
             'gene_result_diplotype': ['id', 'functionphenotypeid', 'diplotype', 'diplotypekey'], 'recommendation': ['id', 'guidelineid', 'drugid', 'drugrecommendation', 'classification']}


def v_cpic_api(tid, r):
    table = tid[len('cpic_api_'):]
    d = v_json(tid, r, 1)
    r.check('row count equals the count the API reported when the target was defined', len(d) == TARGETS[tid]['expected_rows'], (len(d), TARGETS[tid]['expected_rows']))
    r.check('documented keys present in every row', all(all(k in x for k in CPIC_KEYS[table]) for x in d), [k for k in CPIC_KEYS[table] if not all(k in x for x in d)])
    if 'id' in d[0]:
        r.check('ids are unique', len({x['id'] for x in d}) == len(d), len(d))
    if table == 'allele_definition':
        r.check('every allele has a gene symbol and a name', all(x['genesymbol'] and x['name'] for x in d))
        genes = collections.Counter(x['genesymbol'] for x in d)
        r.stats.update(genes=len(genes), with_pharmvar_id=sum(1 for x in d if x['pharmvarid']), structural=sum(1 for x in d if x['structuralvariation']))
    elif table == 'sequence_location':
        r.check('chromosomal notation starts with g. (nuclear) or m. (mitochondrial) and position is an integer', all((x['chromosomelocation'] or 'g.').startswith(('g.', 'm.')) and isinstance(x['position'], int) for x in d))
        r.stats.update(with_dbsnp=sum(1 for x in d if x['dbsnpid']), mitochondrial_locations=sum(1 for x in d if (x['chromosomelocation'] or '').startswith('m.')), mitochondrial_genes=sorted({x['genesymbol'] for x in d if (x['chromosomelocation'] or '').startswith('m.')}))
    elif table == 'allele_location_value':
        r.check('variant alleles are not empty', all(x['variantallele'] for x in d))
    elif table == 'gene_result_diplotype':
        single = collections.Counter(next(iter(x['diplotypekey'])) for x in d if '/' not in x['diplotype'])
        r.check('diplotypes are written allele/allele, except single alleles for the X-linked, mitochondrial and HLA genes', set(single) <= {'G6PD', 'MT-RNR1', 'HLA-A', 'HLA-B'}, dict(single))
        r.stats.update(single_allele_entries=dict(single))
    elif table == 'allele':
        r.stats.update(functional_status={str(k): v for k, v in collections.Counter(x['clinicalfunctionalstatus'] for x in d).most_common(8)})
    r.stats.update(rows=len(d))



# --- IMPC release 24.0 results -----------------------------------------------------------------------------------------
IMPC_GP_COLS = ['marker_accession_id', 'marker_symbol', 'colony_id', 'sex', 'zygosity', 'allele_accession_id', 'allele_symbol', 'mp_term_id', 'mp_term_name', 'p_value', 'effect_size', 'statistical_method']


def v_impc(tid, r):
    p = path(tid)
    if tid == 'impc_r24_readme':
        t = open(p, encoding='utf-8').read()
        r.check('describes the IMPC reports (viability, genotype-phenotype assertions)', t.startswith('# IMPC Reports') and 'Viability' in t and 'Genotype-Phenotype Assertions' in t)
        r.stats['characters'] = len(t)
        return
    with gzip.open(p, 'rt', encoding='utf-8', newline='') as f:
        rd = csv.DictReader(f)
        hdr = rd.fieldnames
        rows = list(rd)
    if tid == 'impc_r24_genotype_phenotype_assertions_all':
        r.check('documented columns present', all(c in hdr for c in IMPC_GP_COLS), [c for c in IMPC_GP_COLS if c not in hdr])
        r.check('rows > 50,000', len(rows) > 50000, len(rows))
        unrel = sum(1 for x in rows if x['allele_accession_id'].startswith('NOT-RELEASED-'))
        r.check('MGI accession for every marker; allele accession is MGI or a NOT-RELEASED placeholder', all(re.match(r'^MGI:\d+$', x['marker_accession_id']) and re.match(r'^(MGI:\d+|NOT-RELEASED-[0-9a-f]+)$', x['allele_accession_id']) for x in rows), unrel)
        onto = collections.Counter(x['mp_term_id'].split(':')[0] for x in rows)
        noterm = sum(1 for x in rows if not x['mp_term_id'])
        r.check('phenotype term ids are MP:nnnnnnn or (histopathology) MPATH:n, with at most 5 rows lacking a term', all(re.match(r'^(MP:\d{7}|MPATH:\d+)$', x['mp_term_id']) for x in rows if x['mp_term_id']) and noterm <= 5, dict(onto))
        zy = collections.Counter(x['zygosity'] for x in rows)
        r.check('zygosity is homozygote, heterozygote or hemizygote', set(zy) <= {'homozygote', 'heterozygote', 'hemizygote'}, dict(zy))
        bad_p = 0
        for x in rows:
            try:
                bad_p += not (0 <= float(x['p_value']) <= 1)
            except ValueError:
                bad_p += x['p_value'] != ''
        r.check('p-values lie in [0, 1]', bad_p == 0, bad_p)
        r.stats.update(rows=len(rows), genes=len({x['marker_accession_id'] for x in rows}), alleles=len({x['allele_accession_id'] for x in rows}), rows_with_unreleased_allele=unrel, rows_without_term=noterm, zygosity=dict(zy), ontologies=dict(onto), terms=len({x['mp_term_id'] for x in rows}), sexes=dict(collections.Counter(x['sex'] for x in rows)))
    elif tid == 'impc_r24_gene_mp_association_1':
        r.check('columns MP Term, Zygosity, # Genes, Gene Symbol (Gene Accession Id)', hdr == ['MP Term', 'Zygosity', '# Genes', 'Gene Symbol (Gene Accession Id)'], hdr)
        off = [(x['# Genes'], len(re.findall(r'MGI:\d+', x['Gene Symbol (Gene Accession Id)']))) for x in rows]
        r.check('# Genes is an integer and equals the number of listed genes, within one for at most 10 rows', all(a.isdigit() for a, _ in off) and all(abs(int(a) - b) <= 1 for a, b in off) and sum(int(a) != b for a, b in off) <= 10, [o for o in off if int(o[0]) != o[1]][:6])
        r.check('rows >= 50', len(rows) >= 50, len(rows))
        r.stats.update(rows=len(rows), rows_with_count_off_by_one=sum(int(a) != b for a, b in off))
    else:
        call = 'Viability Phenotype HOMs/HEMIs'
        r.check('documented columns present', all(c in hdr for c in ('Gene Symbol', 'Gene Accession Id', 'Allele Symbol', 'Zygosity', 'Total # Pups', call, 'Viability Call Method')), hdr[:6])
        r.check('rows > 5,000', len(rows) > 5000, len(rows))
        r.check('MGI accession for every gene', all(re.match(r'^MGI:\d+$', x['Gene Accession Id']) for x in rows))
        vc = collections.Counter(x[call] for x in rows)
        r.check('viability call is viable, subviable or lethal', set(vc) <= {'viable', 'subviable', 'lethal'}, dict(vc))
        r.check('pup counts are non-negative integers', all(x['Total # Pups'].isdigit() for x in rows))
        r.stats.update(rows=len(rows), genes=len({x['Gene Accession Id'] for x in rows}), calls=dict(vc))

# --- per-target dispatch -------------------------------------------------------------------------------------------
def run(tid):
    r = R()
    t = TARGETS[tid]
    if not os.path.exists(path(tid)):
        r.check('file present', False, 'not downloaded')
        return r
    v_manifest(tid, r)
    try:
        if tid == 'ensembl_grch38_reference_regional':
            v_reference_slice(tid, r)
        elif t.get('kind') == 'regional_slice':
            v_regional_slice(tid, r)
        elif t.get('kind') == 'per_entity':
            v_per_entity(tid, r)
        elif t.get('kind') == 'http_dir':
            v_parquet_dir(tid, r)
        elif tid.startswith('refseq_'):
            v_refseq(tid, r)
        elif tid.startswith('giab_v3_6_grch38'):
            v_giab(tid, r)
        elif tid.startswith('proteingym_v1_3_'):
            v_proteingym(tid, r)
        elif tid.startswith('hpa_v25_1_'):
            v_hpa(tid, r)
        elif tid.startswith('biogrid_interactions_'):
            v_biogrid(tid, r)
        elif tid.startswith('hgsvc3_'):
            v_hgsvc3(tid, r)
        elif tid.startswith('cpic_api_'):
            v_cpic_api(tid, r)
        elif tid.startswith('impc_r24_'):
            v_impc(tid, r)
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
            if tid != 'reactome_pathways':
                codes = collections.Counter(x[4] for x in rows if len(x) > 5 and x[5] == 'Homo sapiens')
                r.check('Homo sapiens evidence codes are IEA or TAS', bool(codes) and set(codes) <= {'IEA', 'TAS'}, dict(codes))
                r.stats['evidence_codes_homo_sapiens'] = dict(codes)
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
        elif tid == 'uniprot_humsavar':
            v_humsavar(tid, r)
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
        elif tid == FASTA_ID:
            v_fasta_toplevel(tid, r)
        elif tid == FASTA_ID + '_fai':
            v_fai(tid, r)
        elif tid == FASTA_ID + '_gzi':
            v_gzi(tid, r)
        elif tid == 'ensembl_grch38_gtf_116':
            v_ensembl_gtf(tid, r)
        elif tid == 'pharmvar_6_2_29':
            v_pharmvar(tid, r)
        elif tid.startswith('sifts_'):
            v_sifts(tid, r)
        elif tid == 'gnomad_v4_1_lof_curation_incomplete_penetrance':
            v_lof_curation(tid, r)
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
    text = json.dumps(out, indent=1, sort_keys=True)     # serialise first: a failure must not truncate the report
    open(rp, 'w').write(text)
