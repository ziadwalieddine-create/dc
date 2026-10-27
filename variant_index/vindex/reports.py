"""Reports over a built index: coverage statistics, validation, the lifecycle of every acquisition target, a change log between builds, and the gap report.

Nothing here is written by hand into the index: every number is read from the index tables, the acquisition manifest and the acquisition validation report, so a
claim in a report can be re-derived by running the report again. Statuses follow one vocabulary:
  listed      the target is named in acquisition/targets.json
  acquired    the manifest holds a record for it with ok = true
  verified    its acquisition validator passed (acquisition/validation_report.json)
  mapped      at least one index row from its release resolved to an index entity (mapping_status = mapped)
  integrated  the index holds rows from its release (mapped or retained as unmapped / ambiguous / rejected)
A target can be acquired and verified yet not integrated; that is a gap, and it is listed as one.
"""
import glob
import json
import os
import shutil
import tempfile
from collections import Counter

import pyarrow as pa
import pyarrow.compute as pc
import pyarrow.dataset as pds
import pyarrow.parquet as pq

from . import SCHEMA_VERSION, schema
from .build.common import REPO

# Modules of the system and what is true of each now. 'external' means the input must come from a dedicated tool; 'not_acquired' means the data are not in the environment.
MODULES = (
    dict(name='allele normalisation, identity keys, multiallelic splitting', status='implemented', tests='test_normalise, test_vcfio, test_match_small'),
    dict(name='small-variant matching with verified-reference and absence semantics', status='implemented', tests='test_match_small, test_build_clinvar, test_wgs'),
    dict(name='structural / copy-number matching (type, breakpoints, dosage; no overlap-only equivalence)', status='implemented', tests='test_match_sv, test_wgs'),
    dict(name='star-allele / haplotype calling from defining variants (CPIC, PharmVar)', status='implemented', tests='test_match_hap, test_build_pgx, test_wgs',
         note='structural genes (CYP2D6, HLA, KIR ...) require a dedicated caller; its result can be supplied with --external-haplotypes'),
    dict(name='repeat-expansion classification against STRchive thresholds', status='partial', tests='test_match_repeat, test_build_regions, test_wgs',
         note='classifies the output of ExpansionHunter / TRGT and keeps every repeat genotype of the file (NO_INDEX_RECORD where the index has no locus); sizing a repeat from short reads is not attempted'),
    dict(name='contigs outside the primary assembly (alternate loci, unplaced scaffolds, decoys)', status='implemented', tests='test_other_contigs',
         note='indexed under the shared NA partition and matched in one extra pass of a whole-file run; records on contigs the reference does not know are listed as rejected'),
    dict(name='independent cross-checks with bcftools (carried alleles, ClinVar evidence, ClinVar fields of the index)', status='implemented', tests='test_crosscheck'),
    dict(name='idempotence (two builds, separate processes, different hash seeds) and hash-seed independence of matching', status='implemented', tests='test_idempotence, test_wgs'),
    dict(name='polygenic scores with weight-coverage, ancestry and calibration gates', status='implemented', tests='test_pgs, test_pgs_run, test_pgs_registry',
         note='ancestry inference and reference score distributions are not part of the repository'),
    dict(name='ClinVar aggregate, submission, structural and composite-record evidence', status='implemented', tests='test_build_clinvar, test_build_clinvar_vcv'),
    dict(name='ClinGen, CPIC, PharmVar, ClinPGx, MaveDB, GWAS Catalog, Open Targets, gnomAD, DGV, HGSVC3, GIAB, STRchive, Reactome adapters', status='implemented', tests='test_build_*'),
    dict(name='provenance verification', status='implemented', tests='test_provenance, test_resolvers',
         note='release linkage and manifest agreement are checked for every release; the SHA-256 of the raw record is recomputed from the acquired file for a sample of 500 rows of every release that holds evidence (one rule per source in vindex/resolvers.py); entity rows carry the release of their source but no digest of their own'),
    dict(name='transcript HGVS (c., n.) to genomic coordinate mapping', status='unimplemented', note='records that carry only transcript HGVS stay unmapped with the reason'),
    dict(name='consequence annotation of called variants (VEP / bcftools csq)', status='unimplemented', note='the index stores source-stated consequences; none is computed for a sample'),
    dict(name='GRCh37 -> GRCh38 lift-over of legacy ClinVar CNV records', status='unimplemented', note='those records keep their GRCh37 coordinates and are unmapped with the reason'),
    dict(name='ancestry inference from genotypes (reference panel such as 1000 Genomes)', status='unimplemented', note='the 1000 Genomes metadata is acquired, the phased panels are not (size); polygenic-score ancestry gates use a declared group'),
    dict(name='annotation resources acquired for later use (AlphaMissense, ProteinGym, AlphaFold DB, InterPro, SIFTS, UniProt, GTEx, HPA, ENCODE, BioGRID, 1000 Genomes metadata)', status='unimplemented',
         note='acquired, validated and join-tested (acquisition/compat_report.json) but no adapter writes them into the index; the lifecycle report lists them as acquired and not integrated'),
    dict(name='ontology roll-up of traits (EFO / MONDO ancestors)', status='partial', tests='test_wgs',
         note='the trait ids of a sample\'s matches are grouped by ancestor and therapeutic area through the Open Targets ontology table (96% of the ids of the HG002 run); MedGen, OMIM and MeSH ids are not in it and stay as the sources wrote them; the index itself keeps the sources\' ids'),
    dict(name='GWAS Catalog full summary statistics', status='not_acquired', note='the association table and study metadata are integrated; 181,588 of 232,202 catalogue studies list full summary statistics (per-study files, mostly molecular-QTL sets), none is downloaded and no adapter for them exists; choosing the relevant studies is not done'),
    dict(name='Open Targets colocalisation', status='not_acquired', note='credible sets, fine-mapping and L2G are integrated; the colocalisation datasets (15.4 GiB) are not downloaded and have no adapter'),
    dict(name='personal 30x GRCh38 SNP/indel, CNV and SV calls', status='unavailable_input', note='no personal call file exists in this environment; matching is verified on synthetic fixtures and, when run, on a public truth set'),
)


def _evidence_dataset(index_dir):
    return pds.dataset(os.path.join(index_dir, 'evidence'), format='parquet', ignore_prefixes=[])


def coverage(index_dir):
    """Row counts of every table, and evidence counts by source, level, axis, kind, assertion and mapping status; unmapped reasons by source; distinct variants per axis."""
    out = dict(tables={}, evidence={})
    for name in schema.TABLES:
        p = os.path.join(index_dir, name)
        files = [os.path.join(dp, f) for dp, _, fs in os.walk(p) for f in fs if f.endswith('.parquet')] if os.path.isdir(p) else ([p + '.parquet'] if os.path.exists(p + '.parquet') else [])
        if files:
            out['tables'][name] = dict(files=len(files), rows=sum(pq.ParquetFile(f).metadata.num_rows for f in files))
    if not os.path.isdir(os.path.join(index_dir, 'evidence')):
        return out
    ds = _evidence_dataset(index_dir)
    cols = ['source_id', 'level', 'axis', 'evidence_kind', 'effect_assertion', 'mapping_status', 'mapping_reason']
    counts = {c: Counter() for c in cols}
    cross, reasons = Counter(), Counter()
    for b in ds.to_batches(columns=cols):
        t = pa.Table.from_batches([b])
        for c in cols[:6]:
            for v, n in zip(*[x.to_pylist() for x in pc.value_counts(t[c]).flatten()]):
                counts[c][v] += n
        for r in t.group_by(['source_id', 'axis']).aggregate([([], 'count_all')]).to_pylist():
            cross[(r['source_id'], r['axis'])] += r['count_all']
        for r in t.filter(pc.not_equal(t['mapping_status'], 'mapped')).group_by(['source_id', 'mapping_status', 'mapping_reason']).aggregate([([], 'count_all')]).to_pylist():
            reasons[(r['source_id'], r['mapping_status'], r['mapping_reason'])] += r['count_all']
    out['evidence'] = dict(total=sum(counts['source_id'].values()), **{c: dict(counts[c].most_common()) for c in cols[:6]},
                           source_by_axis={f'{s}|{a}': n for (s, a), n in sorted(cross.items())},
                           not_mapped={f'{s}|{m}|{r}': n for (s, m, r), n in sorted(reasons.items(), key=lambda kv: -kv[1])})
    axes_variants = {}
    vdir = os.path.join(index_dir, 'evidence', 'variant')
    if os.path.isdir(vdir):
        for chrom in sorted(os.listdir(vdir)):
            t = pds.dataset(os.path.join(vdir, chrom), format='parquet', ignore_prefixes=[]).to_table(columns=['subject_id', 'axis', 'mapping_status'], filter=pc.field('mapping_status') == 'mapped')
            for ax in set(t['axis'].to_pylist()):
                sub = t.filter(pc.equal(t['axis'], ax))
                axes_variants[ax] = axes_variants.get(ax, 0) + len(pc.unique(sub['subject_id']))
        out['evidence']['variants_with_evidence_by_axis'] = axes_variants
    out['lineage'] = lineage_sharing(index_dir)
    out['xref_coherence'] = xref_coherence(index_dir)
    return out


def lineage_sharing(index_dir):
    """Lineage groups, and how many of them are written by more than one source. A group two databases share is one study or record seen through both: it counts once in the
    independence tally. None when duckdb is not installed."""
    try:
        import duckdb
    except ImportError:
        return None
    files = sorted(glob.glob(os.path.join(index_dir, 'evidence', '*', '*', '*.parquet')))
    if not files:
        return dict(groups=0, shared_between_sources=0, source_pairs={})
    spill = tempfile.mkdtemp(prefix='vindex_duckdb_')
    try:
        con = duckdb.connect()
        con.execute(f"PRAGMA temp_directory='{spill}'")
        con.execute(f"SET memory_limit='{DUCKDB_MEMORY}'")
        con.execute(f'SET threads={DUCKDB_THREADS}')
        con.execute(f"create view pairs as select distinct lineage_group, source_id from read_parquet({files!r}) where lineage_group is not null")
        groups = con.execute("select count(distinct lineage_group) from pairs").fetchone()[0]
        shared = con.execute("select count(*) from (select lineage_group from pairs group by 1 having count(*) > 1)").fetchone()[0]
        rows = con.execute("select a.source_id, b.source_id, count(*) from pairs a join pairs b on a.lineage_group = b.lineage_group and a.source_id < b.source_id group by 1, 2 order by 3 desc").fetchall()
        con.close()
    finally:
        shutil.rmtree(spill, ignore_errors=True)
    return dict(groups=groups, shared_between_sources=shared, source_pairs={f'{a}|{b}': n for a, b, n in rows})


def xref_coherence(index_dir):
    """dbSNP identifiers that the index resolves to variants at more than one position, and which sources (xref relations) disagree. Sources reach an rsID by different routes
    (ClinVar's RS field, the GWAS Catalog's mapped SNP, PharmVar and ClinPGx variant names), so two sources placing one identifier at two sites is a mapping error, a merged
    identifier or a locus that exists twice (the pseudoautosomal regions): it is counted, not judged. None when duckdb is not installed."""
    try:
        import duckdb
    except ImportError:
        return None
    files = sorted(os.path.join(dp, f) for dp, _, fs in os.walk(os.path.join(index_dir, 'xref')) for f in fs if f.endswith('.parquet'))
    out = dict(identifiers=0, multi_position=0, multi_position_across_sources=0, multi_position_within_one_source=0, source_pairs={}, examples=[])
    if not files:
        return out
    spill = tempfile.mkdtemp(prefix='vindex_duckdb_')
    con = duckdb.connect()
    try:
        con.execute(f"PRAGMA temp_directory='{spill}'")
        con.execute(f"SET memory_limit='{DUCKDB_MEMORY}'")
        con.execute(f'SET threads={DUCKDB_THREADS}')
        con.execute(f"create view placed as select distinct identifier, relation, split_part(entity_id, ':', 1) || ':' || split_part(entity_id, ':', 2) as locus "
                    f"from read_parquet({files!r}) where namespace = 'dbsnp' and entity_level = 'variant'")
        out['identifiers'] = con.execute('select count(distinct identifier) from placed').fetchone()[0]
        con.execute('create table multi as select identifier, count(distinct locus) as loci, count(distinct relation) as relations from placed group by 1 having count(distinct locus) > 1')
        out['multi_position'] = con.execute('select count(*) from multi').fetchone()[0]
        con.execute('create table across as select distinct a.identifier from placed a join placed b on a.identifier = b.identifier and a.relation < b.relation and a.locus <> b.locus')
        out['multi_position_across_sources'] = con.execute('select count(*) from across').fetchone()[0]
        out['multi_position_within_one_source'] = out['multi_position'] - con.execute('select count(*) from multi join (select distinct identifier from placed '
                                                                                       'group by 1 having count(distinct relation) > 1) using (identifier)').fetchone()[0]
        out['source_pairs'] = {f'{a}|{b}': n for a, b, n in con.execute('select a.relation, b.relation, count(distinct a.identifier) as n from placed a join placed b on a.identifier = b.identifier '
                                                                         'and a.relation < b.relation and a.locus <> b.locus group by 1, 2 order by 3 desc').fetchall()}
        for ident, in con.execute('select identifier from across order by identifier limit 10').fetchall():
            out['examples'].append(dict(identifier=ident, placed={r: l for r, l in con.execute('select relation, string_agg(locus, \',\') from placed where identifier = ? group by 1', [ident]).fetchall()}))
    finally:
        con.close()
        shutil.rmtree(spill, ignore_errors=True)
    return out


# Text that a source uses for "no value". As a value in these columns it is a placeholder that was kept instead of dropped (found on the real index: '.' as a condition, 'NR' as
# an ancestry and as a unit, 'n/a' as a review status). 'unknown' and 'not provided' are not here: ClinVar and the Catalog use them as real categories.
MISSING_VALUE_MARKERS = ('', '.', '-', '--', '?', 'NR', 'NA', 'N/A', 'n/a', 'NaN', 'nan', 'None', 'none', 'null', 'NULL')
MARKER_COLUMNS = ('condition_text', 'ancestry', 'score_unit', 'review_status', 'tissue', 'cell_type', 'model_system', 'inheritance', 'penetrance', 'curator', 'assertion_criteria', 'dose', 'perturbation')
FLOAT_COLUMNS = ('claim_num', 'beta', 'se', 'ci_low', 'ci_high', 'odds_ratio', 'p_value', 'pip', 'score', 'eaf')
INTERVAL_TOLERANCE = 0.001           # share of rows with an interval whose estimate may lie outside it: rounding and typing errors in the sources themselves


# Reports must be able to run beside a whole-genome match (about 7 GB) on a 16 GB machine: DuckDB spills to disk beyond this instead of growing.
DUCKDB_MEMORY = os.environ.get('VINDEX_DUCKDB_MEMORY', '4GB')
DUCKDB_THREADS = int(os.environ.get('VINDEX_DUCKDB_THREADS', '2'))

# (table, column) pairs of entity tables whose string cells must not hold a marker text for "no value" (an empty string is allowed in an entity attribute)
ENTITY_MARKER_COLUMNS = (('entity_study', 'ancestry'), ('entity_study', 'cohorts'), ('entity_study', 'biosample'), ('entity_study', 'trait_text'), ('entity_assay', 'target_gene'),
                         ('entity_assay', 'model_system'), ('entity_assay', 'tissue_or_cell'), ('entity_assay', 'perturbation'), ('entity_assay', 'dose_context'), ('entity_assay', 'readout'),
                         ('entity_assay', 'score_direction'), ('entity_haplotype', 'function'), ('entity_haplotype', 'allele_level'), ('entity_gene', 'hgnc_id'), ('entity_gene', 'strand'),
                         ('entity_gene', 'mane_select'), ('entity_gene', 'biotype'), ('entity_region', 'name'), ('entity_event', 'origin'), ('entity_event', 'raw_type'),
                         ('variant_gene', 'consequence'), ('gene_set_member', 'name'))


def _table_files(index_dir, name):
    p = os.path.join(index_dir, name)
    files = sorted(os.path.join(dp, f) for dp, _, fs in os.walk(p) for f in fs if f.endswith('.parquet')) if os.path.isdir(p) else ([p + '.parquet'] if os.path.exists(p + '.parquet') else [])
    return files


def entity_checks(index_dir, check, note):
    """Referential integrity between the tables and invariants of the entity tables (DuckDB). `check(name, violations, detail)` records one result that can fail the validation; `note(name,
    count)` records a count that is reported but does not fail it (gene ids that a source names and the gene table does not hold: the row is kept, the gene just has no symbol). A cross-
    reference that points at nothing silently breaks every lookup that follows it, so each is counted: xref and variant_gene rows, haplotype definitions, pathway members, GWAS study accessions."""
    import duckdb
    spill = tempfile.mkdtemp(prefix='vindex_duckdb_')
    con = duckdb.connect()
    try:
        con.execute(f"PRAGMA temp_directory='{spill}'")
        con.execute(f"SET memory_limit='{DUCKDB_MEMORY}'")
        con.execute(f'SET threads={DUCKDB_THREADS}')
        have = {}
        for name in ('entity_variant', 'entity_event', 'entity_haplotype', 'entity_gene', 'entity_region', 'entity_study', 'entity_assay', 'xref', 'variant_gene', 'haplotype_definition', 'gene_set_member'):
            files = _table_files(index_dir, name)
            if files:
                con.execute(f'create view {name} as select * from read_parquet({files!r}, union_by_name=true)')
                have[name] = True
        one = lambda sql: int(con.execute(sql).fetchone()[0] or 0)
        if 'xref' in have:
            for level, table, col in (('variant', 'entity_variant', 'variant_key'), ('event', 'entity_event', 'event_id'), ('haplotype', 'entity_haplotype', 'haplotype_id'), ('gene', 'entity_gene', 'gene_id')):
                if table in have:
                    check(f'xref_{level}_not_an_entity', one(f"select count(*) from xref where entity_level = '{level}' and entity_id not in (select {col} from {table})"))
        for table, ref_table, ref_col, col, name, hard in (('variant_gene', 'entity_variant', 'variant_key', 'variant_key', 'variant_gene_variant_not_an_entity', True),
                                                            ('haplotype_definition', 'entity_haplotype', 'haplotype_id', 'haplotype_id', 'haplotype_definition_haplotype_not_an_entity', True),
                                                            ('variant_gene', 'entity_gene', 'gene_id', 'gene_id', 'variant_gene_gene_not_in_the_gene_table', False),
                                                            ('gene_set_member', 'entity_gene', 'gene_id', 'gene_id', 'gene_set_member_gene_not_in_the_gene_table', False)):
            if table in have and ref_table in have:
                n = one(f'select count(*) from {table} where {col} not in (select {ref_col} from {ref_table})')
                check(name, n) if hard else note(name, n)
        if 'entity_study' in have and os.path.isdir(os.path.join(index_dir, 'evidence')):
            ev = os.path.join(index_dir, 'evidence', '**', 'gwas_catalog.parquet')
            if glob.glob(ev, recursive=True):
                check('gwas_evidence_study_not_an_entity', one(f"select count(*) from (select distinct split_part(accession, ':', 1) as study from read_parquet('{ev}', union_by_name=true)) "
                                                               f"where study not in (select study_id from entity_study)"))
        if 'entity_event' in have:
            check('entity_event_start_after_end', one('select count(*) from entity_event where start > "end"'))
            check('entity_event_frequency_outside_0_1', one('select count(*) from entity_event where population_af < 0 or population_af > 1'))
            check('entity_event_negative_copy_number', one('select count(*) from entity_event where copy_number < 0'))
        if 'entity_study' in have:
            check('entity_study_cases_exceed_total', one('select count(*) from entity_study where n_cases > n_total'))
            check('entity_study_total_not_positive', one('select count(*) from entity_study where n_total <= 0'))
        for table in ('entity_gene', 'entity_region'):
            if table in have:
                check(f'{table}_start_after_end', one(f'select count(*) from {table} where start > "end"'))
        if 'entity_gene' in have:
            check('entity_gene_strand_not_plus_or_minus', one("select count(*) from entity_gene where strand is not null and strand not in ('+', '-')"))
        if 'entity_variant' in have:
            check('entity_variant_position_not_positive', one('select count(*) from entity_variant where pos < 1'))
            check('entity_variant_key_differs_from_its_parts', one("select count(*) from entity_variant where variant_key <> chrom || ':' || cast(pos as varchar) || ':' || ref || ':' || alt"))
        markers = ', '.join(f"'{m}'" for m in MISSING_VALUE_MARKERS if m)                 # an empty string in an entity attribute reads as "empty"; a marker text does not
        for table, col in ENTITY_MARKER_COLUMNS:
            if table in have and col in {r[0] for r in con.execute(f'describe {table}').fetchall()}:
                check(f'{table}_{col}_is_a_missing_value_marker', one(f'select count(*) from {table} where {col} in ({markers})'))
    finally:
        con.close()
        shutil.rmtree(spill, ignore_errors=True)


def validate(index_dir, reference=None, sample_keys=50000):
    """Checks that hold for every row: controlled vocabulary, mapped / unmapped consistency, key format, ranges. With a reference, the REF allele of a sample of keys
    is read back from the FASTA. With DuckDB available, evidence ids are checked for uniqueness and mapped subjects for existence."""
    res = dict(schema_version=SCHEMA_VERSION, checks={}, informational={}, ok=True)

    def check(name, bad, detail=None):
        res['checks'][name] = dict(violations=bad, detail=detail)
        if bad:
            res['ok'] = False

    def note(name, count):
        res['informational'][name] = count
    if not os.path.isdir(os.path.join(index_dir, 'evidence')):
        res['ok'] = False
        res['checks']['evidence_table_present'] = dict(violations=1)
        return res
    ds = _evidence_dataset(index_dir)
    vocab = dict(level=schema.LEVELS, axis=tuple(schema.AXES), evidence_kind=schema.EVIDENCE_KINDS, effect_assertion=tuple(schema.EFFECT_ASSERTIONS), mapping_status=tuple(schema.MAPPING_STATUS))
    bad = Counter()
    n = 0
    with_interval = outside_interval = 0
    marker_set = pa.array(MISSING_VALUE_MARKERS, pa.string())
    for b in ds.to_batches():
        t = pa.Table.from_batches([b])
        n += t.num_rows
        for col in MARKER_COLUMNS:
            bad[f'{col}_is_a_missing_value_marker'] += int(pc.sum(pc.fill_null(pc.is_in(t[col], value_set=marker_set), False)).as_py() or 0)
        for col in FLOAT_COLUMNS:
            bad[f'{col}_not_finite'] += int(pc.sum(pc.fill_null(pc.or_(pc.is_nan(t[col]), pc.is_inf(t[col])), False)).as_py() or 0)
        bad['eaf_outside_0_1'] += int(pc.sum(pc.fill_null(pc.or_(pc.less(t['eaf'], 0), pc.greater(t['eaf'], 1)), False)).as_py() or 0)
        bad['standard_error_not_positive'] += int(pc.sum(pc.fill_null(pc.less_equal(t['se'], 0), False)).as_py() or 0)
        bad['ci_low_above_ci_high'] += int(pc.sum(pc.fill_null(pc.greater(t['ci_low'], t['ci_high']), False)).as_py() or 0)
        est = pc.coalesce(t['beta'], t['odds_ratio'])
        both = pc.and_(pc.and_(pc.is_valid(est), pc.is_valid(t['ci_low'])), pc.is_valid(t['ci_high']))
        slack = pc.add(pc.multiply(pc.subtract(t['ci_high'], t['ci_low']), 0.01), 1e-9)
        away = pc.or_(pc.less(est, pc.subtract(t['ci_low'], slack)), pc.greater(est, pc.add(t['ci_high'], slack)))
        with_interval += int(pc.sum(both).as_py() or 0)
        outside_interval += int(pc.sum(pc.and_(both, pc.fill_null(away, False))).as_py() or 0)
        for col, allowed in vocab.items():
            bad[f'{col}_outside_vocabulary'] += int(pc.sum(pc.invert(pc.is_in(t[col], value_set=pa.array(allowed, pa.string())))).as_py() or 0)
        for col in ('evidence_id', 'source_id', 'release_id', 'accession', 'locator', 'raw_sha256', 'claim_type'):
            bad[f'{col}_missing'] += t[col].null_count + int(pc.sum(pc.equal(t[col], '')).as_py() or 0)
        mapped = pc.equal(t['mapping_status'], 'mapped')
        bad['mapped_without_subject'] += int(pc.sum(pc.and_(mapped, pc.is_null(t['subject_id']))).as_py() or 0)
        bad['not_mapped_with_subject'] += int(pc.sum(pc.and_(pc.invert(pc.equal(t['mapping_status'], 'ambiguous')), pc.and_(pc.invert(mapped), pc.is_valid(t['subject_id'])))).as_py() or 0)
        bad['not_mapped_without_reason'] += int(pc.sum(pc.and_(pc.invert(mapped), pc.is_null(t['mapping_reason']))).as_py() or 0)
        var_mapped = pc.and_(mapped, pc.equal(t['level'], 'variant'))
        colons = pc.count_substring(pc.fill_null(t['subject_id'], ''), ':')
        bad['variant_key_not_chrom_pos_ref_alt'] += int(pc.sum(pc.and_(var_mapped, pc.not_equal(colons, 3))).as_py() or 0)
        ex = pc.extract_regex(pc.fill_null(t['subject_id'], ''), r'^(?P<c>[^:]+):(?P<p>[0-9]+):(?P<r>[ACGTN]+):(?P<a>[ACGTN]+)$')
        kr, ka = pc.struct_field(ex, 'r'), pc.struct_field(ex, 'a')
        snv = pc.and_(var_mapped, pc.and_(pc.equal(pc.utf8_length(pc.fill_null(kr, '')), 1), pc.equal(pc.utf8_length(pc.fill_null(ka, '')), 1)))
        for col in ('effect_allele', 'other_allele'):
            given = pc.and_(snv, pc.and_(pc.is_valid(t[col]), pc.not_equal(t[col], '')))
            inside = pc.or_(pc.equal(t[col], kr), pc.equal(t[col], ka))
            bad[f'{col}_not_an_allele_of_the_snv_key'] += int(pc.sum(pc.and_(given, pc.invert(pc.fill_null(inside, False)))).as_py() or 0)
        pv, pip = t['p_value'], t['pip']
        bad['p_value_outside_0_1'] += int(pc.sum(pc.or_(pc.less(pv, 0), pc.greater(pv, 1))).as_py() or 0)
        bad['pip_outside_0_1'] += int(pc.sum(pc.or_(pc.less(pip, 0), pc.greater(pip, 1.0000001))).as_py() or 0)
        bad['raw_sha256_not_64_hex'] += int(pc.sum(pc.invert(pc.match_substring_regex(pc.fill_null(t['raw_sha256'], ''), '^[0-9a-f]{64}$'))).as_py() or 0)
    res['evidence_rows_checked'] = n
    for k, v in sorted(bad.items()):
        check(k, v)
    # an estimate beside an interval that does not contain it (a signed beta next to the bounds of its magnitude, for one) is a defect; a handful in the sources' own data is not
    check('estimate_outside_its_interval', max(0, outside_interval - int(INTERVAL_TOLERANCE * with_interval)),
          dict(rows_with_estimate_and_interval=with_interval, outside=outside_interval, tolerated=int(INTERVAL_TOLERANCE * with_interval)))
    spill = tempfile.mkdtemp(prefix='vindex_duckdb_')                  # DuckDB spills to its working directory by default; keep that out of the repository
    try:
        import duckdb
        con = duckdb.connect()
        con.execute(f"PRAGMA temp_directory='{spill}'")
        con.execute(f"SET memory_limit='{DUCKDB_MEMORY}'")
        con.execute(f'SET threads={DUCKDB_THREADS}')
        ev = os.path.join(index_dir, 'evidence', '**', '*.parquet')
        dup = con.execute(f"select count(*) - count(distinct evidence_id) from read_parquet('{ev}', union_by_name=true)").fetchone()[0]
        check('evidence_id_not_unique', int(dup))
        for level, table, col in (('variant', 'entity_variant', 'variant_key'), ('event', 'entity_event', 'event_id'), ('haplotype', 'entity_haplotype', 'haplotype_id')):
            tp = os.path.join(index_dir, table)
            has_table = os.path.isdir(tp) or os.path.exists(tp + '.parquet')
            if os.path.isdir(os.path.join(index_dir, 'evidence', level)) and has_table:
                tf = os.path.join(tp, '**', '*.parquet') if os.path.isdir(tp) else tp + '.parquet'
                miss = con.execute(f"select count(*) from read_parquet('{os.path.join(index_dir, 'evidence', level, '**', '*.parquet')}', union_by_name=true) e where e.mapping_status = 'mapped' "
                                   f"and e.subject_id not in (select {col} from read_parquet('{tf}', union_by_name=true))").fetchone()[0]
                check(f'mapped_{level}_subject_not_an_entity', int(miss))
        con.close()
        entity_checks(index_dir, check, note)
    except ImportError:
        res['checks']['duckdb_checks'] = dict(violations=0, detail='skipped: duckdb is not installed (evidence_id uniqueness and subject existence were not checked)')
    finally:
        shutil.rmtree(spill, ignore_errors=True)
    vdir = os.path.join(index_dir, 'evidence', 'variant')
    if os.path.isdir(vdir):
        misfiled = 0
        for c in sorted(os.listdir(vdir)):
            files = sorted(glob.glob(os.path.join(vdir, c, '*.parquet')))
            if not files:
                continue
            t = pds.dataset(files, format='parquet').to_table(columns=['subject_id', 'mapping_status'], filter=pc.field('mapping_status') == 'mapped')
            if t.num_rows:
                first = pc.list_element(pc.split_pattern(t['subject_id'], ':'), 0)
                if c == schema.OTHER_PARTITION:                      # alternate loci and scaffolds: the key's contig must be outside the primary set
                    misfiled += int(pc.sum(pc.is_in(first, value_set=pa.array(schema.PRIMARY_CONTIGS))).as_py() or 0)
                else:
                    misfiled += int(pc.sum(pc.not_equal(first, c)).as_py() or 0)
        check('mapped_variant_row_filed_under_another_chromosome', misfiled)
    if reference is not None:
        p = os.path.join(index_dir, 'entity_variant')
        files = sorted(glob.glob(os.path.join(p, '*.parquet'))) if os.path.isdir(p) else []
        total = sum(pq.ParquetFile(f).metadata.num_rows for f in files)
        if total:
            step = max(1, total // sample_keys)
            wrong, checked = [], 0
            for f in files:
                t = pq.read_table(f, columns=['variant_key', 'chrom', 'pos', 'ref'])
                for k in t.take(list(range(0, t.num_rows, step))).to_pylist():
                    checked += 1
                    canon = reference.contig(k['chrom'])
                    if canon and reference.check_ref(canon, k['pos'], k['ref'])[0] == 'mismatch':
                        wrong.append(k['variant_key'])
            check('entity_variant_ref_differs_from_reference', len(wrong), wrong[:5])
            res['entity_variant_ref_checked'] = checked
    return res


def lifecycle(index_dir, targets_path=None, manifest_path=None, validation_path=None, catalogue_dir=None):
    """One row per acquisition target with the five statuses and the evidence for each."""
    targets_path = targets_path or os.path.join(REPO, 'acquisition', 'targets.json')
    manifest_path = manifest_path or os.path.join(REPO, 'acquisition', 'manifest.jsonl')
    validation_path = validation_path or os.path.join(REPO, 'acquisition', 'validation_report.json')
    catalogue_dir = catalogue_dir or os.path.join(REPO, 'catalogue', 'sources')
    with open(targets_path) as f:
        targets = {t['id']: t for t in json.load(f)}
    manifest = {}
    if os.path.exists(manifest_path):
        with open(manifest_path) as f:
            manifest = {r['id']: r for r in (json.loads(l) for l in f if l.strip())}
    valid = {}
    if os.path.exists(validation_path):
        with open(validation_path) as f:
            valid = json.load(f)
    catalogued = {os.path.splitext(x)[0] for x in os.listdir(catalogue_dir)} if os.path.isdir(catalogue_dir) else set()
    rel_path = os.path.join(index_dir, 'source_release.parquet')
    releases = pq.read_table(rel_path).to_pylist() if os.path.exists(rel_path) else []
    rows_by_release = Counter()
    mapped_by_release = Counter()
    if os.path.isdir(os.path.join(index_dir, 'evidence')):
        for b in _evidence_dataset(index_dir).to_batches(columns=['release_id', 'mapping_status']):
            t = pa.Table.from_batches([b])
            for r in t.group_by(['release_id', 'mapping_status']).aggregate([([], 'count_all')]).to_pylist():
                rows_by_release[r['release_id']] += r['count_all']
                if r['mapping_status'] == 'mapped':
                    mapped_by_release[r['release_id']] += r['count_all']
    by_target = {}
    for r in releases:
        by_target.setdefault(r['target_id'], []).append(r['release_id'])
    rows = []
    for tid, t in sorted(targets.items()):
        m = manifest.get(tid)
        v = valid.get(tid)
        rids = by_target.get(tid, [])
        ev = sum(rows_by_release[r] for r in rids)
        mp = sum(mapped_by_release[r] for r in rids)
        rows.append(dict(target_id=tid, source=t['source'], listed=True, source_in_catalogue=t['source'] in catalogued, acquired=bool(m and m.get('ok')), verified=bool(v and v.get('passed')),
                         releases_in_index=rids, evidence_rows=ev, mapped_evidence_rows=mp, integrated=bool(rids), mapped=mp > 0, supports=t.get('supports')))
    summary = Counter()
    for r in rows:
        for k in ('listed', 'acquired', 'verified', 'integrated', 'mapped'):
            summary[k] += r[k]
    return dict(targets=rows, summary=dict(summary, total=len(rows)), acquired_not_integrated=[r['target_id'] for r in rows if r['acquired'] and not r['integrated']])


def change_log(prev_dir, cur_dir):
    """What differs between two builds: schema version, releases added / removed / re-versioned, adapter code digests, table row counts."""
    def load(d):
        with open(os.path.join(d, 'index_manifest.json')) as f:
            return json.load(f)
    a, b = load(prev_dir), load(cur_dir)
    ra, rb = set(a.get('releases', [])), set(b.get('releases', []))
    base = lambda r: r.rsplit('@', 1)[0]
    ver_a, ver_b = {base(r): r for r in ra}, {base(r): r for r in rb}
    code_a = {x['adapter']: x['code_sha256_16'] for x in a.get('adapters', [])}
    code_b = {x['adapter']: x['code_sha256_16'] for x in b.get('adapters', [])}
    ta, tb = a.get('tables', {}), b.get('tables', {})
    return dict(schema_version=dict(previous=a.get('schema_version'), current=b.get('schema_version')),
                releases=dict(added=sorted(set(ver_b) - set(ver_a)), removed=sorted(set(ver_a) - set(ver_b)),
                              reversioned={k: dict(previous=ver_a[k], current=ver_b[k]) for k in sorted(set(ver_a) & set(ver_b)) if ver_a[k] != ver_b[k]}),
                adapters=dict(added=sorted(set(code_b) - set(code_a)), removed=sorted(set(code_a) - set(code_b)), code_changed=sorted(k for k in set(code_a) & set(code_b) if code_a[k] != code_b[k])),
                build_code=dict(previous=a.get('build_code'), current=b.get('build_code'), changed=a.get('build_code') != b.get('build_code')),
                tables={k: dict(previous=ta.get(k, {}).get('rows', 0), current=tb.get(k, {}).get('rows', 0), delta=tb.get(k, {}).get('rows', 0) - ta.get(k, {}).get('rows', 0))
                        for k in sorted(set(ta) | set(tb)) if ta.get(k, {}).get('rows', 0) != tb.get(k, {}).get('rows', 0)},
                problems=dict(previous=a.get('problems'), current=b.get('problems')))


def gap_report(index_dir, life=None, cov=None):
    life = life or lifecycle(index_dir)
    cov = cov or coverage(index_dir)
    nm = (cov.get('evidence') or {}).get('not_mapped', {})
    by_source = {}
    for k, n in nm.items():
        s, status, reason = k.split('|', 2)
        by_source.setdefault(s, {}).setdefault(status, {})[reason] = n
    return dict(modules=[dict(m) for m in MODULES], acquired_not_integrated=life['acquired_not_integrated'],
                targets_not_acquired=[r['target_id'] for r in life['targets'] if not r['acquired']], targets_not_verified=[r['target_id'] for r in life['targets'] if r['acquired'] and not r['verified']],
                rows_kept_but_not_mapped={s: dict(total=sum(sum(v.values()) for v in d.values()), by_status={st: dict(sorted(v.items(), key=lambda kv: -kv[1])[:10]) for st, v in d.items()})
                                          for s, d in sorted(by_source.items())})
