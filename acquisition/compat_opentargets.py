"""Joins between the acquired Open Targets 26.09 Parquet datasets and the MaveDB human set (called from compat_tests.py).

Every test reads only the columns it needs, part file by part file, so memory stays small. A test is reported as PASS or FAIL
against a stated criterion, or as INFO where the result describes a relationship rather than a requirement. A test that
raises is recorded as FAIL with the error, never silently skipped. Requires pyarrow (acquisition/requirements.txt).
"""
import collections
import csv
import glob
import gzip
import os


def run(ctx):
    src, log, D = ctx['src'], ctx['log'], ctx['D']
    genes_py, uniq, pos_set = ctx['M_ensg_all'], ctx['uniq'], ctx['pos_set']
    build_index, covered = ctx['build_index'], ctx['covered']
    try:
        import pyarrow as pa
        import pyarrow.compute as pc
        import pyarrow.parquet as pq
    except ImportError:
        src('opentargets', 'Open Targets Parquet tests not run: pyarrow is not installed', None, note='pip install -r acquisition/requirements.txt')
        return
    OT = os.path.join(D, 'opentargets', '26.09')
    if not os.path.isdir(OT):
        src('opentargets', 'Open Targets Parquet tests not run: datasets not downloaded', None)
        return
    genes = pa.array(sorted(genes_py))
    snv_ids = pa.array(sorted(f'{c}_{p}_{r}_{a}' for c, p, r, a in uniq))

    def parts(ds):
        return sorted(glob.glob(os.path.join(OT, ds, '*.parquet')))

    def scan(ds, cols):
        for fp in parts(ds):
            yield pq.read_table(fp, columns=cols)

    def count(mask):
        return int(pc.sum(mask).as_py() or 0)

    def test(name, datasets, fn):
        log('OT test:', name)
        try:
            ok, stats = fn()
            src('opentargets', name, ok, targets=['opentargets_26_09_' + d for d in datasets], **stats)
        except Exception as e:           # a broken test must show up in the report
            src('opentargets', name, False, targets=['opentargets_26_09_' + d for d in datasets], error=repr(e)[:300])

    # ---- gene annotation and its consistency with sources we hold directly ---------------------------------------------
    def t_target():
        ids, sym = set(), {}
        for t in scan('target', ['id', 'approvedSymbol']):
            for i, s in zip(t.column('id').to_pylist(), t.column('approvedSymbol').to_pylist()):
                ids.add(i)
                sym[i] = s
        found = genes_py & ids
        mane_sym = ctx.get('sym_from_ensg', {})
        both = [g for g in found if g in mane_sym]
        same = sum(1 for g in both if sym.get(g) == mane_sym[g])
        return len(found) / max(1, len(genes_py)) > 0.9, dict(target_genes=len(ids), mavedb_ensg=len(genes_py), found=len(found), coverage=round(len(found) / max(1, len(genes_py)), 3),
                                                             symbol_equal_to_MANE_symbol=same, compared=len(both), unmatched_examples=sorted(genes_py - ids)[:8])
    test('OT gene join: MaveDB ENSG -> target annotation', ['target'], t_target)

    def t_constraint():
        gn = {}
        with gzip.open(os.path.join(D, 'gnomad/gnomad.v4.1.2.constraint_metrics.tsv.bgz'), 'rt', errors='replace', newline='') as f:
            for row in csv.DictReader(f, delimiter='\t'):
                if row['gene_id'].startswith('ENSG') and row['mane_select'] == 'true' and row['lof.oe'] not in ('NA', ''):
                    gn[row['gene_id']] = float(row['lof.oe'])
        ot = {}
        for t in scan('target', ['id', 'constraint']):
            for i, cons in zip(t.column('id').to_pylist(), t.column('constraint').to_pylist()):
                for c in cons or []:
                    if c['constraintType'] == 'lof' and c['oe'] is not None:
                        ot[i] = c['oe']
        both = sorted(set(gn) & set(ot))
        same = sum(1 for g in both if abs(gn[g] - ot[g]) < 1e-3)
        diffs = sorted(abs(gn[g] - ot[g]) for g in both)
        return None, dict(genes_in_both=len(both), lof_oe_equal_within_1e3=same, share_equal=round(same / max(1, len(both)), 4),
                          median_abs_difference=round(diffs[len(diffs) // 2], 5) if diffs else None,
                          interpretation='Open Targets repackages gnomAD constraint; the share of equal LoF o/e values shows whether its release matches the v4.1.2 file we hold.')
    test('OT constraint (target.constraint, LoF o/e) vs gnomAD v4.1.2 MANE-select constraint', ['target'], t_constraint)

    # ---- variants, credible sets, studies, locus-to-gene --------------------------------------------------------------
    def t_variant():
        tot = matched = with_af = 0
        methods, pops = collections.Counter(), collections.Counter()
        for t in scan('variant', ['variantId', 'alleleFrequencies', 'variantEffect']):
            tot += t.num_rows
            mask = pc.is_in(t.column('variantId'), value_set=snv_ids)
            matched += count(mask)
            sub = t.filter(mask)
            with_af += sum(1 for x in sub.column('alleleFrequencies').to_pylist() if x)
            eff = pc.list_flatten(t.column('variantEffect').combine_chunks())
            if len(eff):
                for v in pc.value_counts(pc.struct_field(eff, 'method')).to_pylist():
                    methods[v['values']] += v['counts']
            af = pc.list_flatten(t.column('alleleFrequencies').combine_chunks())
            if len(af):
                for v in pc.value_counts(pc.struct_field(af, 'populationName')).to_pylist():
                    pops[v['values']] += v['counts']
        return matched > 0, dict(variants_in_index=tot, mavedb_unique_snvs=len(uniq), matched=matched, match_rate=round(matched / max(1, len(uniq)), 4), matched_with_allele_frequencies=with_af,
                                 variant_effect_methods=dict(methods.most_common(12)), frequency_populations=dict(pops.most_common(12)),
                                 note='the index holds variants that occur in OT credible sets or evidence, so most MaveDB candidate alleles are expected to be absent')
    test('OT exact SNV join: MaveDB-mapped variants -> variant index (chrom_pos_ref_alt)', ['variant'], t_variant)

    cs_ids, cs_study = [], []

    def t_credible():
        n = lead = member = csets = 0
        types = collections.Counter()
        for t in scan('credible_set', ['studyLocusId', 'studyId', 'variantId', 'locus', 'studyType']):
            t = t.combine_chunks()
            n += t.num_rows
            cs_ids.append(t.column('studyLocusId').combine_chunks())
            cs_study.append(pc.unique(t.column('studyId').combine_chunks()))
            lead += count(pc.is_in(t.column('variantId'), value_set=snv_ids))
            for v in pc.value_counts(t.column('studyType')).to_pylist():
                types[v['values']] += v['counts']
            loc = t.column('locus').combine_chunks()
            flat = pc.list_flatten(loc)
            if len(flat):
                hit = pc.is_in(pc.struct_field(flat, 'variantId'), value_set=snv_ids)
                member += count(hit)
                csets += len(set(pc.filter(pc.list_parent_indices(loc), hit).to_pylist()))
        return n > 100_000, dict(credible_sets=n, by_study_type=dict(types.most_common(8)), lead_variants_matching_mavedb_snvs=lead, member_variant_hits=member,
                                 credible_sets_containing_a_mavedb_snv=csets)
    test('OT credible sets: lead and member variants matching MaveDB-mapped SNVs', ['credible_set'], t_credible)

    def t_study():
        studies, types, qtl_genes = [], collections.Counter(), set()
        for t in scan('study', ['studyId', 'studyType', 'geneId']):
            studies.append(t.column('studyId').combine_chunks())
            for v in pc.value_counts(t.column('studyType')).to_pylist():
                types[v['values']] += v['counts']
            qtl_genes |= set(pc.filter(t.column('geneId'), pc.is_in(t.column('geneId'), value_set=genes)).to_pylist())
        sid = pa.concat_arrays(studies)
        cs = pa.concat_arrays([pc.unique(c) for c in cs_study]) if cs_study else pa.array([], pa.string())
        cs = pc.unique(cs)
        resolved = count(pc.is_in(cs, value_set=sid)) if len(cs) else 0
        return (resolved / max(1, len(cs)) >= 0.99), dict(studies=len(sid), by_study_type=dict(types.most_common(10)), mavedb_genes_with_molecular_qtl_studies=len(qtl_genes),
                                                         credible_set_study_ids=len(cs), resolved_in_study_table=resolved, resolution_rate=round(resolved / max(1, len(cs)), 5))
    test('OT study table: molecular-QTL studies for MaveDB genes and credible-set studyId resolution', ['study', 'credible_set'], t_study)

    def t_l2g():
        n = hit = strong = 0
        genes_hit = set()
        loci_in_cs = 0
        cs_all = pa.concat_arrays(cs_ids) if cs_ids else pa.array([], pa.string())
        for t in scan('l2g_prediction', ['studyLocusId', 'geneId', 'score']):
            n += t.num_rows
            m = pc.is_in(t.column('geneId'), value_set=genes)
            hit += count(m)
            sub = t.filter(m)
            s = pc.greater_equal(sub.column('score'), 0.5)
            strong += count(s)
            genes_hit |= set(pc.filter(sub.column('geneId'), s).to_pylist())
            if len(cs_all):
                loci_in_cs += count(pc.is_in(t.column('studyLocusId'), value_set=cs_all))
        return (n > 100_000 and loci_in_cs / max(1, n) >= 0.99), dict(l2g_rows=n, rows_for_mavedb_genes=hit, rows_for_mavedb_genes_score_ge_0_5=strong, mavedb_genes_with_strong_l2g=len(genes_hit),
                                                                       rows_whose_locus_is_in_credible_set=loci_in_cs, locus_resolution_rate=round(loci_in_cs / max(1, n), 5))
    test('OT locus-to-gene: MaveDB genes and credible-set locus resolution', ['l2g_prediction', 'credible_set'], t_l2g)

    # ---- regulation ----------------------------------------------------------------------------------------------------
    def t_e2g():
        n = linked = 0
        iv = collections.defaultdict(list)
        kinds = collections.Counter()
        for t in scan('enhancer_to_gene', ['geneId', 'chromosome', 'start', 'end', 'intervalType']):
            n += t.num_rows
            sub = t.filter(pc.is_in(t.column('geneId'), value_set=genes))
            linked += sub.num_rows
            for c, s, e, k in zip(sub.column('chromosome').to_pylist(), sub.column('start').to_pylist(), sub.column('end').to_pylist(), sub.column('intervalType').to_pylist()):
                iv[c].append((s, e))
                kinds[k] += 1
        ix = {c: build_index(v) for c, v in iv.items()}
        hits = sum(1 for chrom, pos in pos_set if chrom in ix and covered(ix[chrom], pos - 1))
        return n > 1_000_000 and linked > 0, dict(intervals=n, intervals_linked_to_mavedb_genes=linked, by_interval_type=dict(kinds.most_common(8)), positions_tested=len(pos_set),
                                                  positions_inside_a_linked_interval=hits, note='intervals read as half-open [start, end); positions are MaveDB-mapped SNVs, overwhelmingly coding')
    test('OT enhancer-to-gene: intervals linked to MaveDB genes x MaveDB-mapped SNV positions', ['enhancer_to_gene'], t_e2g)

    # ---- interactions ----------------------------------------------------------------------------------------------------
    int_ids = []

    def t_interaction():
        n = rows = bad = 0
        found = set()
        src_db = collections.Counter()
        for t in scan('interaction', ['targetA', 'targetB', 'sourceDatabase', 'scoring', 'interactionId']):
            n += t.num_rows
            int_ids.append(t.column('interactionId').combine_chunks())
            ma, mb = pc.is_in(t.column('targetA'), value_set=genes), pc.is_in(t.column('targetB'), value_set=genes)
            found |= set(pc.filter(t.column('targetA'), ma).to_pylist()) | set(pc.filter(t.column('targetB'), mb).to_pylist())
            m = pc.or_(ma, mb)
            rows += count(m)
            for v in pc.value_counts(pc.filter(t.column('sourceDatabase'), m)).to_pylist():
                src_db[v['values']] += v['counts']
            sc = t.column('scoring').drop_null()
            bad += len(sc) - count(pc.and_(pc.greater_equal(sc, 0), pc.less_equal(sc, 1)))
        return (n > 1_000_000 and bad == 0 and len(found) > 0), dict(interactions=n, rows_with_a_mavedb_gene=rows, mavedb_genes_with_interactions=len(found), coverage=round(len(found) / max(1, len(genes_py)), 3),
                                                                     by_source_database=dict(src_db.most_common(8)), scores_outside_0_1=bad)
    test('OT interaction: MaveDB genes as interaction partners (scores within [0, 1])', ['interaction'], t_interaction)

    def t_interaction_evidence():
        n = ok = 0
        ref = pa.concat_arrays(int_ids) if int_ids else pa.array([], pa.int64())
        for t in scan('interaction_evidence', ['interactionId']):
            n += t.num_rows
            ok += count(pc.is_in(t.column('interactionId'), value_set=ref))
        return (n > 1_000_000 and ok / max(1, n) >= 0.99), dict(evidence_rows=n, interaction_ids_in_interaction_table=ok, resolution_rate=round(ok / max(1, n), 5))
    test('OT interaction_evidence: interactionId resolves in the interaction table', ['interaction_evidence', 'interaction'], t_interaction_evidence)

    # ---- gene-level evidence sets ------------------------------------------------------------------------------------------
    def gene_join(ds, col):
        def f():
            n = 0
            found = set()
            for t in scan(ds, [col]):
                n += t.num_rows
                found |= set(pc.filter(t.column(col), pc.is_in(t.column(col), value_set=genes)).to_pylist())
            return len(found) >= 10, dict(rows=n, mavedb_genes_with_records=len(found), coverage=round(len(found) / max(1, len(genes_py)), 3), join_column=col)
        return f
    for ds, col in [('target_essentiality', 'targetId'), ('evidence_impc', 'targetId'), ('mouse_phenotype', 'targetFromSourceId'), ('evidence_crispr', 'targetId'),
                    ('evidence_crispr_screen', 'targetId'), ('evidence_gene_burden', 'targetId'), ('evidence_genomics_england', 'targetId'), ('evidence_gene2phenotype', 'targetId'),
                    ('evidence_orphanet', 'targetId'), ('evidence_uniprot_variants', 'targetId'), ('evidence_uniprot_literature', 'targetId'), ('homology', 'targetId')]:
        if parts(ds):
            test(f'OT gene join: MaveDB ENSG -> {ds}', [ds], gene_join(ds, col))

    def t_disease_ref():
        dis = set()
        for t in scan('disease', ['id']):
            dis |= set(t.column('id').to_pylist())
        dis_arr = pa.array(sorted(dis))
        n = ok = 0
        for ds in ('evidence_gene_burden', 'evidence_genomics_england', 'evidence_gene2phenotype', 'evidence_orphanet', 'evidence_uniprot_variants', 'evidence_uniprot_literature', 'evidence_impc', 'evidence_crispr', 'evidence_crispr_screen'):
            for t in scan(ds, ['diseaseId']):
                col = t.column('diseaseId').drop_null()
                n += len(col)
                ok += count(pc.is_in(col, value_set=dis_arr))
        return (n > 100_000 and ok / max(1, n) >= 0.99), dict(evidence_rows_checked=n, disease_ids_resolving_in_ontology=ok, resolution_rate=round(ok / max(1, n), 5), ontology_terms=len(dis))
    test('OT evidence diseaseId resolves in the disease ontology', ['disease', 'evidence_gene_burden', 'evidence_genomics_england', 'evidence_gene2phenotype', 'evidence_orphanet', 'evidence_uniprot_variants',
                                                                  'evidence_uniprot_literature', 'evidence_impc', 'evidence_crispr', 'evidence_crispr_screen'], t_disease_ref)

    def t_pgx():
        n = hit = 0
        rs_pgx = set()
        for t in scan('pharmacogenomics', ['variantRsId', 'targetFromSourceId']):
            n += t.num_rows
            rs_pgx |= set(pc.filter(t.column('targetFromSourceId'), pc.is_in(t.column('targetFromSourceId'), value_set=genes)).to_pylist())
        return None, dict(rows=n, mavedb_genes_with_pgx_records=len(rs_pgx), note='join column targetFromSourceId; overlap is informative, not a pass criterion')
    if parts('pharmacogenomics'):
        test('OT gene join: MaveDB ENSG -> pharmacogenomics', ['pharmacogenomics'], t_pgx)
