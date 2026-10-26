"""Polygenic scores for one sample, one chromosome at a time, with the registry's gates beside every result.

A score is a model applied to the sample's dosages, not evidence attached to variants. Each weight is scored from a called genotype or from positive reference evidence;
a weight the file cannot answer for is missing, never zero, and the reasons are counted. Three gates are reported with every score and none changes the arithmetic:
weight coverage (share of absolute weight actually scored), ancestry (development / evaluation ancestry from the registry against the sample's inferred group, which
this repository cannot infer: it needs a reference panel) and calibration (a z-score only when a reference distribution for the same score and ancestry is supplied).
The numbers are research outputs, not risk predictions.
"""
import json
import os

import pyarrow as pa
import pyarrow.parquet as pq

from . import pgs
from .callability import CallableRegionCallability
from .match_small import MatchContext
from .pgs_registry import ancestry_meta
from .vcfio import CallSet, QCPolicy, read_vcf
from .wgs import chromosome_order, code_version, load_bed_intervals, sha256_file

S, LS, F, I32 = pa.string(), pa.list_(pa.string()), pa.float64(), pa.int32()
SCORES = pa.schema([('pgs_id', S), ('name', S), ('reported_trait', S), ('status', S), ('reason', S), ('scoring_file_sha256', S), ('n_weights', I32), ('n_scored', I32), ('n_missing', I32),
                    ('scored_weight_fraction', F), ('score_sum', F), ('missing_reasons_json', S), ('ancestry_status', S), ('ancestry_json', S), ('sample_ancestry', S), ('calibration_status', S),
                    ('z_score', F), ('flags', LS), ('registry_gates', LS)])
WEIGHTS = pa.schema([('row', I32), ('chrom', S), ('pos', I32), ('effect_allele', S), ('other_allele', S), ('weight', F), ('dosage', I32), ('basis_or_reason', S)])


def scoring_path(scoring_dir, pgs_id):
    return os.path.join(scoring_dir, f'{pgs_id}_hmPOS_GRCh38.txt.gz')


def weight_groups(models, ids, budget):
    """Split the selected scores into groups whose scoring files hold at most `budget` weights together (a larger score is a group of its own), so memory stays bounded
    however many scores are selected; the sample's calls are read once per group and chromosome."""
    groups, cur, n = [], [], 0
    for pid in ids:
        size = (models.get(pid) or {}).get('scoring_file_variants') or (models.get(pid) or {}).get('n_variants') or 1_000_000
        if cur and n + size > budget:
            groups.append(cur)
            cur, n = [], 0
        cur.append(pid)
        n += size
    if cur:
        groups.append(cur)
    return groups


def score_models(vcf, reference, scoring_dir, models, pgs_ids, out_dir, sample=None, qc=None, sample_ancestry=None, reference_distributions=None, callable_bed=None, build='GRCh38',
                 per_weight=False, repo_root=None, weights_budget=4_000_000):
    """models: {pgs_id: registry row} (see pgs_registry); reference_distributions: {pgs_id: {ancestry group: {'mean': m, 'sd': s}}} from a reference panel, optional."""
    from .wgs import check_output_location
    repo_root = repo_root or os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..'))
    check_output_location(out_dir, repo_root)
    os.makedirs(out_dir, exist_ok=True)
    qc = qc or QCPolicy()
    rows, scorable = [], []
    bases = {}
    for pid in sorted(set(pgs_ids)):
        m = models.get(pid, {})
        bases[pid] = dict(pgs_id=pid, name=m.get('name'), reported_trait=m.get('reported_trait'), scoring_file_sha256=m.get('scoring_file_sha256'), registry_gates=list(m.get('gates') or []),
                          sample_ancestry=sample_ancestry)
        if os.path.exists(scoring_path(scoring_dir, pid)):
            scorable.append(pid)
        else:
            rows.append(dict(bases[pid], status='not_scored', reason='scoring_file_not_acquired', flags=['not_scored']))
    bed = load_bed_intervals(callable_bed, reference) if callable_bed else None
    empty = CallSet(sample=sample or '', path=vcf, assembly={})
    sample_name = sample
    order = chromosome_order(reference)
    for group in weight_groups(models, scorable, weights_budget):
        loaded = {}
        for pid in group:
            base = bases[pid]
            sf = pgs.read_scoring_file(scoring_path(scoring_dir, pid))
            if sf.meta.get('pgs_id', pid) != pid:
                rows.append(dict(base, status='not_scored', reason='scoring_file_id_mismatch', flags=['not_scored']))
                continue
            if sf.meta.get('HmPOS_build') != build:
                rows.append(dict(base, status='not_scored', reason=f"scoring_file_harmonised_to_{sf.meta.get('HmPOS_build') or 'nothing'}_not_{build}", flags=['not_scored']))
                continue
            by = {}
            for w in sf.weights:
                canon = (reference.contig(w.chrom) or reference.contig('chr' + w.chrom)) if w.chrom else None
                by.setdefault(canon or '', []).append(w)
            sf.weights = None                                                                  # the per-chromosome lists hold them now
            loaded[pid] = (base, sf, by, [])
        needed = [c for c in order if any(c in by for _b, _s, by, _p in loaded.values())]
        for chrom in needed:
            autosome = chrom in pgs.AUTOSOMES
            cs = read_vcf(vcf, reference, sample, qc, contigs=[chrom]) if autosome else empty        # non-autosomal weights are excluded by rule and need no genotypes
            sample_name = sample_name or cs.sample
            providers = [CallableRegionCallability(((chrom, a, b) for a, b in bed.get(chrom, ())), cs, 'callable_region')] if bed is not None and autosome else []
            ctx = MatchContext(cs, providers)
            for pid, (_base, sf, by, per) in loaded.items():
                per.extend(pgs.dosages(by.pop(chrom, ()), ctx, reference))
            del cs, ctx
        for pid, (base, sf, by, per) in loaded.items():
            per.extend(pgs.dosages(by.get('', ()), MatchContext(empty), reference))              # weights with no usable position: counted as missing with the reason
            dist = ((reference_distributions or {}).get(pid) or {}).get(sample_ancestry) if sample_ancestry else None
            r = pgs.finish(pid, per, sample_ancestry, ancestry_meta(models.get(pid, {})), dist)
            rows.append(dict(base, status='scored', reason=None, n_weights=r.n_weights, n_scored=r.n_scored, n_missing=r.n_missing, scored_weight_fraction=r.abs_weight_scored_fraction,
                             score_sum=r.score_sum, missing_reasons_json=json.dumps(r.reasons, sort_keys=True), ancestry_status=r.ancestry['status'],
                             ancestry_json=json.dumps({k: v for k, v in r.ancestry.items() if k != 'status'}, sort_keys=True), calibration_status=r.calibration['status'],
                             z_score=r.calibration.get('z'), flags=list(r.flags)))
            if per_weight:
                os.makedirs(os.path.join(out_dir, 'pgs_weights'), exist_ok=True)
                pq.write_table(pa.Table.from_pylist([dict(row=p[0], chrom=p[1], pos=p[2], effect_allele=p[3], other_allele=p[4], weight=p[5], dosage=p[6], basis_or_reason=p[7]) for p in sorted(per, key=lambda x: x[0])],
                                                    schema=WEIGHTS), os.path.join(out_dir, 'pgs_weights', f'{pid}.parquet'), compression='zstd')
        loaded.clear()
    rows.sort(key=lambda r: r['pgs_id'])
    pq.write_table(pa.Table.from_pylist([{n: r.get(n) for n in SCORES.names} for r in rows], schema=SCORES), os.path.join(out_dir, 'pgs_scores.parquet'), compression='zstd')
    manifest = dict(input=dict(path=os.path.abspath(vcf), sha256=sha256_file(vcf), sample=sample_name), models=len(rows), scored=sum(1 for r in rows if r['status'] == 'scored'),
                    sample_ancestry=sample_ancestry, code=code_version(repo_root),
                    note='raw weighted sums are research outputs, not risk predictions; they are interpretable only when the coverage, ancestry and calibration gates are met')
    with open(os.path.join(out_dir, 'pgs_manifest.json'), 'w') as f:
        json.dump(manifest, f, indent=1, sort_keys=True)
    return manifest
