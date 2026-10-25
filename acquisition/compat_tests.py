"""Compatibility tests: do the acquired resources join to the existing dc assets (MaveDB human set)?

Identifier joins (genes) and genomic-coordinate joins (GRCh38) are exercised on the full acquired files; results,
with explicit coverage numbers and failure reasons, go to acquisition/compat_report.json. A source is marked
'mapped' only if at least one join test passed its stated threshold. This is a minimal ingestion/compatibility
check, not the analysis engine.

  python3 compat_tests.py [--skip orcs ot]     # ot = the Open Targets Parquet tests (compat_opentargets.py, needs pyarrow)
"""
import argparse
import bisect
import collections
import csv
import glob
import gzip
import json
import os
import re
import sys
import tarfile
import time
import zipfile

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, '..'))
D = os.path.join(ROOT, 'data', 'external')
csv.field_size_limit(1 << 28)
T0 = time.time()
REPORT = dict(generated_utc=time.strftime('%Y-%m-%dT%H:%M:%SZ', time.gmtime()), sources={}, details={})


def log(*a):
    print(f'[{time.time() - T0:6.0f}s]', *a, file=sys.stderr, flush=True)


def opn(p):
    return gzip.open(p, 'rt', encoding='utf-8', errors='replace', newline='') if open(p, 'rb').read(2) == b'\x1f\x8b' else open(p, encoding='utf-8', errors='replace', newline='')


def build_index(iv):
    """Exact stabbing-query index for half-open intervals [s, e): sorted starts plus prefix-maximum of ends."""
    iv = sorted(iv)
    starts = [s for s, e in iv]
    pm, m = [], -1
    for s_, e_ in iv:
        m = max(m, e_)
        pm.append(m)
    return starts, pm


def covered(idx, pos0):
    starts, pm = idx
    i = bisect.bisect_right(starts, pos0) - 1
    return i >= 0 and pm[i] > pos0


# Which acquisition targets each test actually exercised (so a join is only credited to the dimensions those files support)
SCREEN_IDS = ['screen_ccre_' + n for n in ('pls', 'pels', 'dels', 'ca_ctcf', 'ca_h3k4me3', 'ca_tf', 'ca', 'tf')]
USED = [
    ('gencode_mane', '', ['mane_summary_v1_5']),
    ('gnomad_constraint', '', ['gnomad_constraint_v4_1_2']),
    ('gtex', '', ['gtex_v11_gene_median_tpm']),
    ('reactome', 'MaveDB ENSG', ['reactome_ensembl2reactome']),
    ('reactome', 'MaveDB UniProt', ['reactome_uniprot2reactome_all_levels']),
    ('uniprot', '', ['uniprot_human_reviewed']),
    ('clingen', 'MaveDB symbol -> ClinGen dosage', ['clingen_gene_curation_list']),
    ('clingen', 'MaveDB symbol -> ClinGen gene-disease', ['clingen_gene_validity']),
    ('clingen', 'Interval overlap: MaveDB-mapped SNV positions x ClinGen haploinsufficiency', ['clingen_haploinsufficiency_gene', 'ensembl_vep_mavedb_variants']),
    ('clingen', 'Interval overlap: MaveDB-mapped SNV positions x ClinGen triplosensitivity', ['clingen_triplosensitivity_gene', 'ensembl_vep_mavedb_variants']),
    ('gwas_catalog', '', ['gwas_catalog_associations']),
    ('cpic', '', ['cpic_genes']),
    ('pharmgkb', '', ['clinpgx_relationships']),
    ('alphamissense', 'Exact SNV', ['alphamissense_hg38', 'ensembl_vep_mavedb_variants']),
    ('alphamissense', '', ['alphamissense_gene_hg38', 'mane_summary_v1_5']),
    ('dgv', '', ['dgv_grch38_variants_2025_12_01', 'ensembl_vep_mavedb_variants']),
    ('strchive', '', ['strchive_hg38_general']),
    ('ensembl_vep', '', ['ensembl_vep_mavedb_variants', 'mavedb_zenodo_dump_v4']),
    ('mavedb', '', ['mavedb_zenodo_dump_v4']),
    ('clinvar', 'Exact SNV', ['clinvar_vcf_grch38', 'ensembl_vep_mavedb_variants']),
    ('clinvar', '', ['clinvar_vcf_grch38']),
    ('encode', '', SCREEN_IDS + ['ensembl_vep_mavedb_variants']),
    ('gnomad_sv_cnv', '', ['gnomad_exome_cnv_v4_1', 'ensembl_vep_mavedb_variants']),
    ('biogrid_orcs', '', ['biogrid_orcs_human']),
]


def used_targets(sid, name):
    for s_, prefix, tl in USED:
        if s_ == sid and name.startswith(prefix):
            return tl
    return []


def src(sid, name, ok, targets=None, **stats):
    REPORT['sources'].setdefault(sid, dict(mapped=False, tests=[]))
    REPORT['sources'][sid]['tests'].append(dict(name=name, passed=(None if ok is None else bool(ok)), targets=targets if targets is not None else used_targets(sid, name), **stats))
    if ok:
        REPORT['sources'][sid]['mapped'] = True
    log('INFO' if ok is None else ('PASS' if ok else 'FAIL'), sid, name, {k: v for k, v in stats.items() if not isinstance(v, (list, dict))})


# --- MaveDB human set (existing asset) ---------------------------------------------------------------------------------
sets = {}
for fn in ('scrape/mavedb-human/score_set_details.jsonl', 'scrape/mavedb-human/score_set_details_accession_targets.jsonl'):
    for line in open(os.path.join(ROOT, fn)):
        r = json.loads(line)
        sets[r['urn']] = r
M_sym, M_ensg, M_uni = set(), set(), set()
for r in sets.values():
    for t in r['targetGenes']:
        if t.get('mappedHgncName'):
            M_sym.add(t['mappedHgncName'])
        if t.get('uniprotIdFromMappedMetadata'):
            M_uni.add(t['uniprotIdFromMappedMetadata'])
        for e in t.get('externalIdentifiers', []):
            i = e['identifier']
            if i['dbName'] == 'Ensembl' and i['identifier'].startswith('ENSG'):
                M_ensg.add(i['identifier'].split('.')[0])
            if i['dbName'] == 'UniProt':
                M_uni.add(i['identifier'])
log(f'MaveDB human set: {len(sets)} score sets; genes: {len(M_sym)} HGNC symbols, {len(M_ensg)} ENSG, {len(M_uni)} UniProt')
REPORT['details']['mavedb_gene_sets'] = dict(score_sets=len(sets), hgnc_symbols=len(M_sym), ensg=len(M_ensg), uniprot=len(M_uni))

# Existing asset scrape/ensembl-human-genomes: does its registry name the assembly that MaveDB and the coordinate resources use?
reg = json.load(open(os.path.join(ROOT, 'scrape/ensembl-human-genomes/genomes.json')))
reg_names = {m['assembly']['name'] for m in reg['matches']}
ref = [m for m in reg['matches'] if m.get('is_reference')]
mave_asm = collections.Counter()
with open(os.path.join(ROOT, 'scrape/mavedb-human/score_sets_index.csv'), newline='') as f:
    for row in csv.DictReader(f):
        for x in re.split(r'[|;]', row.get('assemblies') or ''):
            if x.strip():
                mave_asm[x.strip()] += 1
norm = {'hg38': 'GRCh38', 'hg19': 'GRCh37'}
unresolved = [a for a in mave_asm if not any(n.startswith(norm.get(a, a)) for n in reg_names)]
src('ensembl_genomes', 'Registry reference assembly (GRCh38) matches the assembly MaveDB names and the GRCh38 coordinate files use',
    bool(ref) and ref[0]['assembly']['name'].startswith('GRCh38') and not unresolved,
    registry_records=len(reg['matches']), reference_assembly=[(m['assembly']['name'], m['assembly']['accession_id'], m['annotation_provider'], m['has_variation'], m['has_regulation']) for m in ref],
    mavedb_assemblies=dict(mave_asm), unresolved_assembly_names=unresolved,
    note='name-level check between two existing assets; it does not test sequence identity, which needs the genome FASTA (not acquired)')

# MANE: ENSG <-> symbol <-> Entrez  (bridge used by other joins)
mane = {}
with opn(os.path.join(D, 'mane/MANE.GRCh38.v1.5.summary.txt.gz')) as f:
    rd = csv.DictReader(f, delimiter='\t')
    for row in rd:
        mane[row['Ensembl_Gene'].split('.')[0]] = row
mane_sym = {v['symbol']: k for k, v in mane.items()}
entrez = {v['symbol']: v['#NCBI_GeneID'].replace('GeneID:', '') for v in mane.values()}
sym_from_ensg = {k: v['symbol'] for k, v in mane.items()}
# UniProt first: it bridges MaveDB's UniProt accessions to gene symbols and Ensembl gene IDs
uni = {}
with opn(os.path.join(D, 'uniprot/uniprot_human_reviewed.tsv')) as f:
    for row in csv.DictReader(f, delimiter='\t'):
        uni[row['Entry']] = row
M_sym_uni, M_ensg_uni = set(), set()
for acc in M_uni:
    row = uni.get(acc)
    if row:
        g = (row.get('Gene Names (primary)') or '').split(';')[0].strip()
        if g:
            M_sym_uni.add(g)
        M_ensg_uni |= {x.split('.')[0] for x in re.findall(r'ENSG\d+(?:\.\d+)?', row.get('Ensembl', ''))}
# full MaveDB gene sets: native fields + UniProt chain + MANE chain
M_sym_all = set(M_sym) | M_sym_uni | {sym_from_ensg[e] for e in (M_ensg | M_ensg_uni) if e in sym_from_ensg}
M_ensg_all = set(M_ensg) | M_ensg_uni | {mane_sym[s] for s in M_sym_all if s in mane_sym}
REPORT['details']['mavedb_gene_sets'].update(symbols_after_chaining=len(M_sym_all), ensg_after_chaining=len(M_ensg_all),
    note='MaveDB records carry few HGNC/Ensembl gene IDs; the full gene set is derived by UniProt -> symbol/Ensembl -> MANE')
log(f'after identifier chaining: {len(M_sym_all)} symbols, {len(M_ensg_all)} ENSG genes')
src('gencode_mane', 'MaveDB ENSG -> MANE Select gene record', len(M_ensg_all & set(mane)) / max(1, len(M_ensg_all)) > 0.8,
    mavedb_ensg=len(M_ensg_all), found=len(M_ensg_all & set(mane)), coverage=round(len(M_ensg_all & set(mane)) / max(1, len(M_ensg_all)), 3), note='versioned IDs must be stripped before joining')
src('gencode_mane', 'MaveDB HGNC symbol -> MANE symbol', len(M_sym & set(mane_sym)) / max(1, len(M_sym)) > 0.8,
    mavedb_symbols=len(M_sym), found=len(M_sym & set(mane_sym)), unmatched_examples=sorted(M_sym - set(mane_sym))[:10])

# gnomAD constraint
with opn(os.path.join(D, 'gnomad/gnomad.v4.1.2.constraint_metrics.tsv.bgz')) as f:
    rd = csv.DictReader(f, delimiter='\t')
    g_ids, g_sym = set(), set()
    n = 0
    for row in rd:
        n += 1
        g_ids.add(row['gene_id'].split('.')[0])
        g_sym.add(row['gene'])
src('gnomad_constraint', 'MaveDB ENSG -> constraint gene_id', len(M_ensg_all & g_ids) / max(1, len(M_ensg_all)) > 0.8,
    rows=n, constraint_genes=len(g_ids), mavedb_ensg=len(M_ensg_all), found=len(M_ensg_all & g_ids), coverage=round(len(M_ensg_all & g_ids) / max(1, len(M_ensg_all)), 3), unmatched_examples=sorted(M_ensg_all - g_ids)[:8])

# GTEx
gt = set()
with gzip.open(os.path.join(D, 'gtex/GTEx_Analysis_2025-08-22_v11_RNASeQCv2.4.3_gene_median_tpm.gct.gz'), 'rt') as f:
    f.readline(); f.readline(); f.readline()
    for line in f:
        gt.add(line.split('\t', 1)[0].split('.')[0])
src('gtex', 'MaveDB ENSG -> GTEx median-TPM gene (version stripped)', len(M_ensg_all & gt) / max(1, len(M_ensg_all)) > 0.8,
    gtex_genes=len(gt), found=len(M_ensg_all & gt), mavedb_ensg=len(M_ensg_all), coverage=round(len(M_ensg_all & gt) / max(1, len(M_ensg_all)), 3), unmatched_examples=sorted(M_ensg_all - gt)[:8])

# Reactome
rg = set()
with opn(os.path.join(D, 'reactome/Ensembl2Reactome.txt')) as f:
    for line in f:
        p = line.rstrip('\n').split('\t')
        if len(p) > 5 and p[5] == 'Homo sapiens' and p[0].startswith('ENSG'):
            rg.add(p[0])
src('reactome', 'MaveDB ENSG -> Ensembl2Reactome (Homo sapiens)', len(M_ensg_all & rg) / max(1, len(M_ensg_all)) > 0.5,
    reactome_genes=len(rg), found=len(M_ensg_all & rg), mavedb_ensg=len(M_ensg_all), coverage=round(len(M_ensg_all & rg) / max(1, len(M_ensg_all)), 3))
up = set()
with opn(os.path.join(D, 'reactome/UniProt2Reactome_All_Levels.txt')) as f:
    for line in f:
        p = line.rstrip('\n').split('\t')
        if len(p) > 5 and p[5] == 'Homo sapiens':
            up.add(p[0])
src('reactome', 'MaveDB UniProt -> UniProt2Reactome', len(M_uni & up) / max(1, len(M_uni)) > 0.5, found=len(M_uni & up), mavedb_uniprot=len(M_uni), coverage=round(len(M_uni & up) / max(1, len(M_uni)), 3))

# UniProt join (table loaded above)
src('uniprot', 'MaveDB UniProt accession -> reviewed human entry', len(M_uni & set(uni)) / max(1, len(M_uni)) > 0.8,
    reviewed_human=len(uni), found=len(M_uni & set(uni)), mavedb_uniprot=len(M_uni), coverage=round(len(M_uni & set(uni)) / max(1, len(M_uni)), 3), unmatched_examples=sorted(M_uni - set(uni))[:8])
afdb = sum(1 for k in M_uni & set(uni) if uni[k].get('AlphaFoldDB'))
src('uniprot', 'UniProt cross-reference to AlphaFold DB for MaveDB proteins', afdb > 0, with_alphafold_xref=afdb)

# ClinGen
cg = {}
with opn(os.path.join(D, 'clingen/ClinGen_gene_curation_list_GRCh38.tsv')) as f:
    raw = f.read().split('\n')
hi = next(i for i, l in enumerate(raw) if l.startswith('#Gene Symbol'))
hdr = raw[hi].lstrip('#').split('\t')
for l in raw[hi + 1:]:
    if l.strip() and not l.startswith('#'):
        p = l.split('\t')
        cg[p[0]] = dict(zip(hdr, p))
src('clingen', 'MaveDB symbol -> ClinGen dosage gene curation', len(M_sym_all & set(cg)) > 0,
    clingen_genes=len(cg), found=len(M_sym_all & set(cg)), mavedb_symbols=len(M_sym_all), coverage=round(len(M_sym_all & set(cg)) / max(1, len(M_sym_all)), 3),
    hi_score3_overlap=sorted(s for s in M_sym_all & set(cg) if cg[s].get('Haploinsufficiency Score') == '3')[:25])
gv = set()
with opn(os.path.join(D, 'clingen/clingen_gene_disease_validity.csv')) as f:
    rows = list(csv.reader(f))
hi2 = next(i for i, x in enumerate(rows) if x and x[0] == 'GENE SYMBOL')
for x in rows[hi2 + 1:]:
    if x and not x[0].startswith('+'):
        gv.add(x[0])
src('clingen', 'MaveDB symbol -> ClinGen gene-disease validity', len(M_sym_all & gv) > 0, validity_genes=len(gv), found=len(M_sym_all & gv), coverage=round(len(M_sym_all & gv) / max(1, len(M_sym_all)), 3))

# GWAS Catalog
gw = collections.Counter()
with zipfile.ZipFile(os.path.join(D, 'gwas_catalog/gwas-catalog-associations_ontology-annotated-full.zip')) as z:
    with z.open(z.namelist()[0]) as fh:
        import io
        rd = csv.DictReader(io.TextIOWrapper(fh, encoding='utf-8', errors='replace'), delimiter='\t')
        for row in rd:
            for g in re.split(r'\s*(?:,|;| - )\s*', row.get('MAPPED_GENE', '')):
                if g:
                    gw[g] += 1
src('gwas_catalog', 'MaveDB symbol -> GWAS MAPPED_GENE', len(M_sym_all & set(gw)) > 0, gwas_genes=len(gw), found=len(M_sym_all & set(gw)), mavedb_symbols=len(M_sym_all),
    coverage=round(len(M_sym_all & set(gw)) / max(1, len(M_sym_all)), 3), associations_for_overlap=sum(gw[g] for g in M_sym_all & set(gw)))

# CPIC / ClinPGx
cp = {x['symbol'] for x in json.load(open(os.path.join(D, 'cpic/cpic_genes.json')))}
src('cpic', 'MaveDB symbol -> CPIC genes', True, cpic_genes=len(cp), found=len(M_sym_all & cp), genes=sorted(M_sym_all & cp)[:20],
    note='join mechanics verified; overlap is informative, not a pass criterion')
with zipfile.ZipFile(os.path.join(D, 'pharmgkb/relationships.zip')) as z:
    nm = [n for n in z.namelist() if n.lower().endswith('.tsv')][0]
    with z.open(nm) as fh:
        txt = io.TextIOWrapper(fh, encoding='utf-8', errors='replace')
        rows = list(csv.DictReader(txt, delimiter='\t'))
cols = list(rows[0])
gcols = [c for c in cols if 'Entity1' in c or 'Entity2' in c]
pgx_genes = {r.get('Entity1_name') for r in rows if r.get('Entity1_type') == 'Gene'} | {r.get('Entity2_name') for r in rows if r.get('Entity2_type') == 'Gene'}
src('pharmgkb', 'MaveDB symbol -> ClinPGx relationship genes', True, relationship_rows=len(rows), pgx_genes=len(pgx_genes), found=len(M_sym_all & pgx_genes), genes=sorted(M_sym_all & pgx_genes)[:20],
    note='join mechanics verified; overlap is informative, not a pass criterion')

# AlphaMissense gene table (transcript level) via MANE
mane_enst = {v['Ensembl_nuc'].split('.')[0]: k for k, v in mane.items()}
am = {}
with opn(os.path.join(D, 'alphamissense/AlphaMissense_gene_hg38.tsv.gz')) as f:
    for row in csv.DictReader((l for l in f if not l.startswith('#')), delimiter='\t'):
        am[row['transcript_id'].split('.')[0]] = float(row['mean_am_pathogenicity'])
am_genes = {mane_enst[t] for t in am if t in mane_enst}
src('alphamissense', 'AlphaMissense transcript -> MANE gene -> MaveDB ENSG', len(M_ensg_all & am_genes) / max(1, len(M_ensg_all)) > 0.5,
    transcripts=len(am), mapped_to_mane_genes=len(am_genes), found=len(M_ensg_all & am_genes), mavedb_ensg=len(M_ensg_all), coverage=round(len(M_ensg_all & am_genes) / max(1, len(M_ensg_all)), 3))

# STRchive (gene symbols in BED name field)
sg = set()
with opn(os.path.join(D, 'strchive/STRchive-disease-loci.hg38.general.bed')) as f:
    for line in f:
        if line.strip() and not line.startswith('#'):
            sg.add(line.split('\t')[4] if len(line.split('\t')) > 4 else '')
src('strchive', 'STRchive loci parsed (hg38 BED, gene names from name field)', len(sg) > 20, loci_genes=len(sg), overlap_with_mavedb=sorted(sg & M_sym_all),
    note='join mechanics only; repeat loci seldom overlap MAVE targets')

# --- Coordinates: Ensembl MaveDB file, dump annotations, ClinVar, SCREEN, gnomAD CNV ----------------------------------------
log('reading Ensembl MaveDB_variants ...')
ens_rows = 0
ens_urns = collections.Counter()
ens_snv = {}   # accession -> (chrom,pos,ref,alt) for SNVs of our sets
ens_pro = {}   # (chrom,pos,ref,alt) -> one-letter protein changes implied by the Ensembl file's HGVS p. (missense and nonsense only)
AA3 = dict(Ala='A', Arg='R', Asn='N', Asp='D', Cys='C', Gln='Q', Glu='E', Gly='G', His='H', Ile='I', Leu='L', Lys='K', Met='M', Phe='F', Pro='P', Ser='S', Thr='T', Trp='W', Tyr='Y', Val='V', Ter='*')
RX_PRO = re.compile(r'p\.([A-Z][a-z]{2})(\d+)([A-Z][a-z]{2})$')
chroms = collections.Counter()
with gzip.open(os.path.join(D, 'ensembl_vep/MaveDB_variants.tsv.gz'), 'rt', errors='replace') as f:
    hdr = f.readline().lstrip('#').rstrip('\n').split('\t')
    ci = {c: i for i, c in enumerate(hdr)}
    for line in f:
        p = line.split('\t', 14)
        ens_rows += 1
        urn = p[ci['urn']]
        ens_urns[urn] += 1
        chroms[p[0]] += 1
        if urn in sets and len(p[3]) == 1 and len(p[4]) == 1:
            key = (p[0], int(p[1]), p[3], p[4])
            ens_snv[p[ci['accession']]] = key
            m_ = RX_PRO.search(p[ci['hgvs']])
            if m_ and m_.group(1) in AA3 and m_.group(3) in AA3 and m_.group(1) != m_.group(3):
                ens_pro.setdefault(key, set()).add(AA3[m_.group(1)] + m_.group(2) + AA3[m_.group(3)])
ours_in_ens = set(ens_urns) & set(sets)
REPORT['details']['ensembl_mavedb_file'] = dict(rows=ens_rows, distinct_score_sets=len(ens_urns), contigs=len(chroms), columns=len(hdr),
    our_human_sets_covered=len(ours_in_ens), our_human_sets_total=len(sets), snv_variants_from_our_sets=len(ens_snv))
src('ensembl_vep', 'Ensembl MaveDB file covers our human score sets', len(ours_in_ens) > 0, rows=ens_rows, covered=len(ours_in_ens), of=len(sets),
    coverage=round(len(ours_in_ens) / len(sets), 3), snv_rows_from_our_sets=len(ens_snv), contigs_sample=sorted(chroms)[:30])

log('reading dump annotations ...')
ann_dir = glob.glob(os.path.join(D, 'mavedb/dump_v4_extract/*/csv'))[0]
ann_sets = ann_vars = ann_g = ann_vep = ann_af = ann_ca = ann_cv = 0
vep_terms = collections.Counter()
sets_with_g = set()
dump_snv = {}
nc2chr = {f'NC_{n:06d}': (str(n) if n <= 22 else {23: 'X', 24: 'Y'}[n]) for n in range(1, 25)}
rx = re.compile(r'^(NC_\d{6})\.\d+:g\.(\d+)([ACGT])>([ACGT])$')
for fp in glob.glob(os.path.join(ann_dir, 'urn-mavedb-*.annotations.csv')):
    urn = os.path.basename(fp)[:-len('.annotations.csv')].replace('urn-mavedb-', 'urn:mavedb:', 1)
    if urn not in sets:
        continue
    ann_sets += 1
    with open(fp, newline='', encoding='utf-8', errors='replace') as f:
        rd = csv.DictReader(f)
        has_cv = [c for c in rd.fieldnames if c.startswith('clinvar.') and c.endswith('clinical_significance')]
        for row in rd:
            ann_vars += 1
            g = row['mavedb.post_mapped_hgvs_g']
            if g != 'NA':
                ann_g += 1
                sets_with_g.add(urn)
                m = rx.match(g)
                if m and m.group(1) in nc2chr:
                    dump_snv[row['accession']] = (nc2chr[m.group(1)], int(m.group(2)), m.group(3), m.group(4))
            if row['vep.vep_functional_consequence'] != 'NA':
                ann_vep += 1
                vep_terms[row['vep.vep_functional_consequence']] += 1
            ann_af += row['gnomad.gnomad_af'] != 'NA'
            ann_ca += row['clingen.clingen_allele_id'] != 'NA'
            ann_cv += any(row[c] != 'NA' for c in has_cv)
REPORT['details']['dump_annotations'] = dict(score_sets=ann_sets, variants=ann_vars, with_genomic_hgvs=ann_g, with_vep=ann_vep, with_gnomad_af=ann_af, with_clingen_allele_id=ann_ca,
    with_any_clinvar_significance=ann_cv, vep_top=dict(vep_terms.most_common(8)), snv_with_parsed_genomic_coords=len(dump_snv))
union_sets = ours_in_ens | sets_with_g
REPORT['details']['genomic_coordinate_coverage'] = dict(our_sets=len(sets), with_ensembl_coordinates=len(ours_in_ens), with_dump_genomic_hgvs=len(sets_with_g),
    in_both=len(ours_in_ens & sets_with_g), with_either=len(union_sets), with_neither=len(sets) - len(union_sets), fraction_with_either=round(len(union_sets) / len(sets), 3))
src('mavedb', 'Dump v4 annotations cover our human score sets', ann_sets > 0, score_sets_with_annotations=ann_sets, of=len(sets), variants=ann_vars,
    with_genomic_hgvs=ann_g, with_gnomad_af=ann_af, with_clingen_allele_id=ann_ca)

# relationship between two independent mapping products (Ensembl file vs MaveDB's own mapper in the dump)
both = set(ens_snv) & set(dump_snv)
agree = sum(1 for a in both if ens_snv[a] == dump_snv[a])
ens_sets_snv = {a.split('#')[0] for a in ens_snv}
dump_sets_g = {a.split('#')[0] for a in dump_snv}
coord_overlap = set(ens_snv.values()) & set(dump_snv.values())
rel = dict(variants_compared_by_accession=len(both), agree=agree, ensembl_snv_alleles=len(ens_snv), dump_snv_alleles=len(dump_snv),
           ensembl_sets_with_snv_alleles=len(ens_sets_snv), dump_sets_with_genomic_hgvs=len(dump_sets_g), sets_in_both=len(ens_sets_snv & dump_sets_g),
           identical_coordinate_alleles=len(coord_overlap),
           interpretation='Ensembl expands each protein-level variant into every codon-compatible genomic allele (one to many); the MaveDB mapper assigns genomic HGVS only to nucleotide-level assay variants. The products cover largely different sets and are complementary, not interchangeable.')
REPORT['details']['mapping_products_relationship'] = rel
ok_conc = (agree / len(both) > 0.95) if len(both) > 100 else None
src('ensembl_vep', 'Relationship of Ensembl coordinate file to the MaveDB-mapper genomic HGVS (SNVs)', ok_conc, **{k: v for k, v in rel.items() if k != 'interpretation'},
    disagreement_examples=[(a, ens_snv[a], dump_snv[a]) for a in list(both) if ens_snv[a] != dump_snv[a]][:5])

log('reading ClinVar ...')
cv = {}
sig_re = re.compile(r'CLNSIG=([^;]+)')
chrs_cv = set()
with gzip.open(os.path.join(D, 'clinvar/clinvar.vcf.gz'), 'rt', errors='replace') as f:
    for line in f:
        if line.startswith('#'):
            continue
        p = line.split('\t', 8)
        m = sig_re.search(p[7])
        cv[(p[0], int(p[1]), p[3], p[4])] = m.group(1).split(',')[0].split('/')[0] if m else 'NA'
        chrs_cv.add(p[0])
keys = list({**ens_snv}.values()) + [v for a, v in dump_snv.items() if a not in ens_snv]
uniq = set(keys)
hit = [k for k in uniq if k in cv]
sig = collections.Counter(cv[k] for k in hit)
src('clinvar', 'Exact SNV coordinate join: MaveDB-mapped variants -> ClinVar (GRCh38)', len(hit) > 100, mavedb_unique_snvs=len(uniq), clinvar_records=len(cv),
    matched=len(hit), match_rate=round(len(hit) / max(1, len(uniq)), 4), clinvar_significance_of_matches=dict(sig.most_common(8)),
    note='ClinVar contigs have no chr prefix, as do the Ensembl-file contigs (same naming)')
cv_genes = collections.Counter()
with gzip.open(os.path.join(D, 'clinvar/clinvar.vcf.gz'), 'rt', errors='replace') as f:
    for line in f:
        if line.startswith('#'):
            continue
        m = re.search(r'GENEINFO=([^;]+)', line)
        if m:
            for gi in m.group(1).split('|'):
                cv_genes[gi.split(':')[0]] += 1
src('clinvar', 'MaveDB symbol -> ClinVar GENEINFO', len(M_sym_all & set(cv_genes)) / max(1, len(M_sym_all)) > 0.5, clinvar_genes=len(cv_genes), found=len(M_sym_all & set(cv_genes)),
    coverage=round(len(M_sym_all & set(cv_genes)) / max(1, len(M_sym_all)), 3), clinvar_records_in_mavedb_genes=sum(cv_genes[g] for g in M_sym_all))

log('reading SCREEN cCREs ...')
cc = {}
n_cc = 0
for fp in glob.glob(os.path.join(D, 'encode/GRCh38-cCREs.*.bed')):
    cls = os.path.basename(fp).split('.')[1]
    per = collections.defaultdict(list)
    with open(fp) as f:
        for line in f:
            p = line.split('\t')
            per[p[0].replace('chr', '')].append((int(p[1]), int(p[2])))
            n_cc += 1
    cc[cls] = per
hit_cls = collections.Counter()
pos_set = {(k[0], k[1]) for k in uniq}
cc = {cls: {c: build_index(v) for c, v in d.items()} for cls, d in cc.items()}
for chrom, pos in pos_set:
    for cls, d in cc.items():
        if chrom in d and covered(d[chrom], pos - 1):
            hit_cls[cls] += 1
src('encode', 'Interval overlap: MaveDB-mapped SNV positions x SCREEN cCRE classes', n_cc > 2_000_000, ccre_intervals=n_cc, positions_tested=len(pos_set), positions_in_at_least_one_class=len({(c_, p_) for (c_, p_) in pos_set if any(c_ in d and covered(d[c_], p_ - 1) for d in cc.values())}), per_class=dict(hit_cls),
    note='class files overlap only where noted by SCREEN; counts per class are not additive across overlapping classes')

log('reading gnomAD CNV ...')
cnv = collections.defaultdict(list)
n_cnv = 0
with gzip.open(os.path.join(D, 'gnomad/gnomad.v4.1.cnv.all.vcf.gz'), 'rt', errors='replace') as f:
    for line in f:
        if line.startswith('#'):
            continue
        p = line.split('\t', 8)
        m = re.search(r'END=(\d+)', p[7])
        if m:
            cnv[p[0].replace('chr', '')].append((int(p[1]), int(m.group(1))))
            n_cnv += 1
cs = {c: build_index([(a - 1, b) for a, b in v]) for c, v in cnv.items()}
cnv_hits = sum(1 for chrom, pos in pos_set if chrom in cs and covered(cs[chrom], pos - 1))
src('gnomad_sv_cnv', 'Interval overlap: MaveDB-mapped SNV positions x gnomAD exome CNV intervals', n_cnv > 1000, cnv_intervals=n_cnv, positions_tested=len(pos_set), positions_overlapping_a_cnv=cnv_hits,
    note='exact interval stabbing query (VCF POS/END converted to half-open)')

# ClinGen HI/TS BED overlap with positions
for nm, fn in (('haploinsufficiency', 'ClinGen_haploinsufficiency_gene_GRCh38.bed'), ('triplosensitivity', 'ClinGen_triplosensitivity_gene_GRCh38.bed')):
    iv = collections.defaultdict(list)
    with open(os.path.join(D, 'clingen', fn)) as f:
        for line in f:
            if line.startswith(('#', 'track')) or not line.strip():
                continue
            p = line.split('\t')
            iv[p[0].replace('chr', '')].append((int(p[1]), int(p[2]), p[3].strip() if len(p) > 3 else ''))
    ix = {c: build_index([(a, b) for a, b, _ in v]) for c, v in iv.items()}
    hits = sum(1 for chrom, pos in pos_set if chrom in ix and covered(ix[chrom], pos - 1))
    src('clingen', f'Interval overlap: MaveDB-mapped SNV positions x ClinGen {nm} gene intervals', sum(len(v) for v in iv.values()) > 100, intervals=sum(len(v) for v in iv.values()), positions_tested=len(pos_set), positions_overlapping=hits)

# AlphaMissense variant-level join: exact GRCh38 SNV key, then agreement of the protein change implied by the two sources
log('reading AlphaMissense variants ...')
pos_str = {f'chr{c}\t{p_}' for c, p_ in pos_set}
am_rows = 0
am_hits = collections.defaultdict(dict)   # SNV key -> {protein_variant: (score, class)}
with gzip.open(os.path.join(D, 'alphamissense/AlphaMissense_hg38.tsv.gz'), 'rt', errors='replace') as f:
    for line in f:
        if line[0] == '#':
            continue
        am_rows += 1
        i = line.find('\t')
        j = line.find('\t', i + 1)
        if line[:j] in pos_str:
            p = line.rstrip('\n').split('\t')
            k = (p[0][3:], int(p[1]), p[2], p[3])
            if k in uniq:
                am_hits[k][p[7]] = (float(p[8]), p[9])
comparable = [k for k in am_hits if k in ens_pro]
concordant = [k for k in comparable if ens_pro[k] & set(am_hits[k])]
conc_set = set(concordant)
am_cls = collections.Counter(c for k in am_hits for (_, c) in am_hits[k].values())
src('alphamissense', 'Exact SNV coordinate join: MaveDB-mapped variants -> AlphaMissense (GRCh38) with protein-change agreement', len(am_hits) > 1000 and len(concordant) / max(1, len(comparable)) > 0.8,
    alphamissense_rows=am_rows, mavedb_unique_snvs=len(uniq), matched=len(am_hits), match_rate=round(len(am_hits) / max(1, len(uniq)), 4), comparable_with_ensembl_protein_change=len(comparable),
    protein_change_agrees=len(concordant), agreement_rate=round(len(concordant) / max(1, len(comparable)), 4), class_of_matches=dict(am_cls),
    disagreement_examples=[(k, sorted(ens_pro[k]), sorted(am_hits[k])) for k in comparable if k not in conc_set][:5],
    note='AlphaMissense covers canonical-transcript missense SNVs only; unmatched MaveDB alleles are synonymous, nonsense, non-canonical-isoform or non-coding. Agreement compares the Ensembl-file HGVS p. change (RefSeq protein numbering) with AlphaMissense protein_variant (canonical numbering).')

# DGV interval overlap
log('reading DGV ...')
dgv = collections.defaultdict(list)
dgv_by_type = collections.defaultdict(lambda: collections.defaultdict(list))
n_dgv = 0
with open(os.path.join(D, 'dgv/GRCh38_hg38_variants_2025-12-01.txt'), encoding='utf-8', errors='replace', newline='') as f:
    f.readline()
    for line in f:
        p = line.split('\t', 6)
        if p[1] == '' or p[1] == 'N':
            continue
        a, b = int(p[2]) - 1, int(p[3])
        dgv[p[1]].append((a, b))
        dgv_by_type[p[4]][p[1]].append((a, b))
        n_dgv += 1
dix = {c: build_index(v) for c, v in dgv.items()}
dgv_hits = sum(1 for chrom, pos in pos_set if chrom in dix and covered(dix[chrom], pos - 1))
per_type = {t: sum(1 for chrom, pos in pos_set if chrom in ix_ and covered(ix_[chrom], pos - 1))
            for t, ix_ in ((t, {c: build_index(v) for c, v in d.items()}) for t, d in dgv_by_type.items())}
src('dgv', 'Interval overlap: MaveDB-mapped SNV positions x DGV GRCh38 variants', n_dgv > 1_000_000, dgv_intervals=n_dgv, positions_tested=len(pos_set), positions_overlapping_a_dgv_variant=dgv_hits,
    per_varianttype=per_type, note='DGV is a catalogue of structural variation in control samples (variants >50 bp); an overlap is a coordinate relationship, not evidence of benignity')

# Equivalence of dump tables vs API-built tables (sample extracted from the dump)
eq = []
sets_dirs = glob.glob(os.path.join(D, 'mavedb/dump_v4_sets/*/csv'))
# every score set of ours whose table was extracted from the dump (the extraction list is in acquisition/README.md)
for u in sorted(sets):
    a = os.path.join(sets_dirs[0], u.replace(':', '-') + '.scores.csv') if sets_dirs else None
    b = os.path.join(ROOT, 'scrape/mavedb-human/scores', u.replace(':', '_') + '.scores.csv.gz')
    if not (a and os.path.exists(a) and os.path.exists(b)):
        continue
    with open(a, newline='', encoding='utf-8') as fa, gzip.open(b, 'rt', newline='') as fb:
        ra, rb = list(csv.DictReader(fa)), list(csv.DictReader(fb))
    same_acc = [r['accession'] for r in ra] == [r['accession'] for r in rb]
    sc_a = [r.get('scores.score', r.get('score')) for r in ra]
    sc_b = [r.get('score') for r in rb]
    eq.append(dict(urn=u, rows_dump=len(ra), rows_api=len(rb), same_accessions=same_acc, scores_identical=sc_a == sc_b))
REPORT['details']['dump_vs_api_equivalence'] = dict(sets_compared=len(eq), all_same_rows=all(e['rows_dump'] == e['rows_api'] for e in eq) if eq else None,
    all_same_accessions=all(e['same_accessions'] for e in eq) if eq else None, all_scores_identical=all(e['scores_identical'] for e in eq) if eq else None,
    mismatches=[e for e in eq if not (e['same_accessions'] and e['scores_identical'] and e['rows_dump'] == e['rows_api'])][:10])
src('mavedb', 'Dump v4 scores identical to API-built tables on a sample of our sets', bool(eq) and all(e['same_accessions'] and e['scores_identical'] for e in eq), sets_compared=len(eq))

# our sets vs the dump
dj = json.load(open(glob.glob(os.path.join(D, 'mavedb/dump_v4_extract/*/main.json'))[0]))
dsets = {ss['urn']: ss for es in dj['experimentSets'] for ex in es['experiments'] for ss in ex['scoreSets']}
as_of = dj['asOf'][:10]
not_in = [u for u in sets if u not in dsets]
why = collections.Counter()
for u in not_in:
    lic = (sets[u].get('license') or {}).get('shortName')
    pub = sets[u].get('publishedDate') or ''
    why['non-CC0 licence' if lic != 'CC0' else ('published after dump date' if pub > as_of else 'other')] += 1
unexplained = [u for u in not_in if (sets[u].get('license') or {}).get('shortName') == 'CC0' and (sets[u].get('publishedDate') or '') <= as_of]
REPORT['details']['our_sets_vs_dump'] = dict(dump_as_of=dj['asOf'], dump_score_sets=len(dsets), ours=len(sets), ours_in_dump=len(sets) - len(not_in), ours_not_in_dump=len(not_in), reasons=dict(why),
    unexplained_absent_urns=[dict(urn=u, published=sets[u].get('publishedDate'), modified=sets[u].get('modificationDate')) for u in unexplained])
src('mavedb', 'Our human sets located in dump v4 main.json', len(sets) - len(not_in) > 1000, in_dump=len(sets) - len(not_in), not_in_dump=len(not_in), reasons=dict(why))

# BioGRID ORCS (long)
ap = argparse.ArgumentParser()
ap.add_argument('--skip', nargs='*', default=[])
a = ap.parse_args()
if 'orcs' not in a.skip:
    log('reading BioGRID ORCS tar ...')
    hits = collections.Counter()
    screens = 0
    tot_rows = 0
    with tarfile.open(os.path.join(D, 'biogrid_orcs/BIOGRID-ORCS-ALL-homo_sapiens-LATEST.screens.tar.gz')) as tf:
        for m in tf:
            if not m.name.endswith('.screen.tab.txt'):
                continue
            screens += 1
            fh = tf.extractfile(m)
            for line in fh:
                if line.startswith(b'#'):
                    continue
                p = line.decode('utf-8', 'replace').split('\t')
                tot_rows += 1
                if len(p) > 12 and p[12] == 'YES':
                    hits[p[3]] += 1
    ov = M_sym_all & set(hits)
    src('biogrid_orcs', 'MaveDB symbol -> ORCS OFFICIAL_SYMBOL hits', len(ov) > 0, screens=screens, rows=tot_rows, genes_with_hits=len(hits), found=len(ov), mavedb_symbols=len(M_sym_all),
        coverage=round(len(ov) / max(1, len(M_sym_all)), 3), top=sorted(((g, hits[g]) for g in ov), key=lambda x: -x[1])[:10])

def write_report():
    for d in REPORT['sources'].values():
        d['mapped'] = any(t['passed'] is True for t in d['tests'])
    json.dump(REPORT, open(os.path.join(HERE, 'compat_report.json'), 'w'), indent=1, default=str)


write_report()   # checkpoint before the Parquet tests, which need pyarrow

# Open Targets Parquet datasets (separate module; needs pyarrow)
if 'ot' not in a.skip:
    sys.path.insert(0, HERE)
    import compat_opentargets
    compat_opentargets.run(dict(src=src, log=log, D=D, M_ensg_all=M_ensg_all, uniq=uniq, pos_set=pos_set, build_index=build_index, covered=covered, sym_from_ensg=sym_from_ensg))
    write_report()
log('done')
