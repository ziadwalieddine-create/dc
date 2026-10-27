"""PGS Catalog model registry: the metadata that gates the interpretation of a polygenic score, kept apart from the variant index.

A polygenic score is a whole model, not a variant-level fact: it is applied to a sample with vindex.pgs (allele dosage, missingness, ancestry, calibration),
never by attaching its weights to indexed variants as evidence. The registry lists every score in the catalogue (acquired or not) so that coverage is stated:

  pgs_model        one row per score: names, traits (reported and EFO), development method, original genome build, number of variants and interaction terms, weight
                   type, publication, the three ancestry distributions (variant-association source, development, evaluation) as the catalogue gives them and parsed,
                   licence text, and the state of its GRCh38 harmonised scoring file: acquired (with digest and variant counts) or listed_not_acquired
  pgs_performance  one row per reported metric: odds ratio, hazard ratio, beta, AUROC, C-index and the free-form metrics, with the confidence interval when given,
                   the covariates and the evaluation sample set; the original text is kept
  pgs_sample       development and evaluation samples: size, cases, controls, sex, age, broad ancestry, country, cohorts

A score that is not additive (interaction terms), has no GRCh38 file, or was developed and evaluated without the sample's ancestry stays in the registry with the reason
that gates it. The catalogue publishes no reference distribution: without one a raw score is uncalibrated and the scorer says so.
"""
import csv
import hashlib
import json
import os
import re

import pyarrow as pa
import pyarrow.parquet as pq

S, LS, F, I32 = pa.string(), pa.list_(pa.string()), pa.float64(), pa.int32()
SCHEMAS = {
    'pgs_model': pa.schema([('pgs_id', S), ('name', S), ('reported_trait', S), ('efo_ids', LS), ('efo_labels', LS), ('development_method', S), ('development_details', S), ('original_genome_build', S),
                            ('n_variants', I32), ('n_interaction_terms', I32), ('weight_type', S), ('pgp_id', S), ('pmid', S), ('doi', S), ('matches_publication', S),
                            ('ancestry_gwas', S), ('ancestry_development', S), ('ancestry_evaluation', S), ('ancestry_gwas_json', S), ('ancestry_development_json', S),
                            ('ancestry_evaluation_json', S), ('ftp_link', S), ('release_date', S), ('licence', S), ('scoring_file', S), ('scoring_file_tier', S), ('scoring_file_sha256', S),
                            ('scoring_file_variants', I32), ('gates', LS)]),
    'pgs_performance': pa.schema([('ppm_id', S), ('pgs_id', S), ('pss_id', S), ('pgp_id', S), ('reported_trait', S), ('metric_family', S), ('metric_name', S), ('value', F), ('ci_low', F),
                                  ('ci_high', F), ('text', S), ('covariates', S), ('other_information', S), ('pmid', S)]),
    'pgs_sample': pa.schema([('stage', S), ('sample_set_id', S), ('pgs_id', S), ('n_individuals', I32), ('n_cases', I32), ('n_controls', I32), ('pct_male', F), ('age', S), ('ancestry_broad', S),
                             ('ancestry_detail', S), ('country', S), ('additional_ancestry', S), ('phenotype_definition', S), ('followup', S), ('gwas_study', S), ('pmid', S), ('cohorts', S)]),
}
METRICS = (('Hazard Ratio (HR)', 'HR'), ('Odds Ratio (OR)', 'OR'), ('Beta', 'Beta'), ('Area Under the Receiver-Operating Characteristic Curve (AUROC)', 'AUROC'), ('Concordance Statistic (C-index)', 'Cindex'))
NOT_LOG_SCALE = {'or', 'hr', 'or/hr', 'odds ratio', 'hazard ratio', 'odds ratio over expected risk'}
NUM = r'-?\d+(?:\.\d+)?(?:[eE][-+]?\d+)?'
VALUE = re.compile(rf'^\s*(?:(?P<name>[^:=\[\]]*?)\s*[:=]\s*)?(?P<val>{NUM})\s*(?:\[\s*(?P<lo>{NUM})\s*[,;]\s*(?P<hi>{NUM})\s*\])?')


def num(s):
    try:
        return float(s)
    except (TypeError, ValueError):
        return None


def integer(s):
    v = num(s)
    return int(v) if v is not None and v == v else None


def ancestry(text):
    """'European:72.7|Not Reported:18.2' -> {'European': 72.7, 'Not Reported': 18.2}; entries without a number are dropped, an empty text gives {}."""
    out = {}
    for part in (text or '').split('|'):
        k, _, v = part.rpartition(':')
        if k.strip() and num(v) is not None:
            out[k.strip()] = num(v)
    return out


def metric(text):
    """(name, value, ci_low, ci_high) of one reported metric string such as '1.55 [1.52,1.58]' or 'R2: 0.07 [0.05, 0.09]'; None when it holds no number."""
    m = VALUE.match(text or '')
    if not m:
        return None
    return (m.group('name') or None, num(m.group('val')), num(m.group('lo')), num(m.group('hi')))


def read_csv(path):
    with open(path, encoding='utf-8', errors='replace', newline='') as f:
        yield from csv.DictReader(f)


def gates_for(row, acquired, scoring):
    """Reasons a score must not be interpreted without a further step: stated, never silently applied."""
    g = []
    if not acquired:
        g.append('no_grch38_scoring_file_acquired')
    if (integer(row['Number of Interaction Terms']) or 0) > 0:
        g.append('interaction_terms_not_additive')
    wt = row['Type of Variant Weight'].strip()
    if wt.lower() in NOT_LOG_SCALE:
        g.append('weights_not_on_log_scale:' + wt)          # a sum of dosage x weight is the model only when the weights add
    elif wt in ('', 'NR'):
        g.append('weight_type_not_reported')
    if not ancestry(row['Ancestry Distribution (%) - PGS Evaluation']) and not ancestry(row['Ancestry Distribution (%) - Score Development/Training']):
        g.append('no_ancestry_information')
    if scoring and scoring.get('variants_header') and scoring.get('data_lines') and scoring['variants_header'] != scoring['data_lines']:
        g.append('scoring_file_variant_count_differs_from_header')
    g.append('uncalibrated_without_reference_distribution')
    return g


def build_registry(meta_dir, scoring_manifest, out_dir):
    """Write the three tables and registry_manifest.json under out_dir; returns the manifest. `scoring_manifest`: the acquisition TSV of stored scoring files."""
    os.makedirs(out_dir, exist_ok=True)
    scoring = {}
    if scoring_manifest and os.path.exists(scoring_manifest):
        with open(scoring_manifest, newline='') as f:
            for r in csv.DictReader(f, delimiter='\t'):
                scoring[r['pgs_id']] = r
    models = []
    for r in read_csv(os.path.join(meta_dir, 'pgs_all_metadata_scores.csv')):
        pid = r['Polygenic Score (PGS) ID']
        sc = scoring.get(pid)
        acquired = bool(sc and sc.get('status') == 'stored')
        a_g, a_d, a_e = (r[f'Ancestry Distribution (%) - {k}'] for k in ('Source of Variant Associations (GWAS)', 'Score Development/Training', 'PGS Evaluation'))
        models.append(dict(pgs_id=pid, name=r['PGS Name'], reported_trait=r['Reported Trait'], efo_ids=[x.strip() for x in r['Mapped Trait(s) (EFO ID)'].split('|') if x.strip()],
                           efo_labels=[x.strip() for x in r['Mapped Trait(s) (EFO label)'].split('|') if x.strip()], development_method=r['PGS Development Method'],
                           development_details=r['PGS Development Details/Relevant Parameters'], original_genome_build=r['Original Genome Build'], n_variants=integer(r['Number of Variants']),
                           n_interaction_terms=integer(r['Number of Interaction Terms']), weight_type=r['Type of Variant Weight'], pgp_id=r['PGS Publication (PGP) ID'], pmid=r['Publication (PMID)'],
                           doi=r['Publication (doi)'], matches_publication=r['Score and results match the original publication'], ancestry_gwas=a_g, ancestry_development=a_d,
                           ancestry_evaluation=a_e, ancestry_gwas_json=json.dumps(ancestry(a_g), sort_keys=True), ancestry_development_json=json.dumps(ancestry(a_d), sort_keys=True),
                           ancestry_evaluation_json=json.dumps(ancestry(a_e), sort_keys=True), ftp_link=r['FTP link'], release_date=r['Release Date'], licence=r['License/Terms of Use'],
                           scoring_file='acquired' if acquired else 'listed_not_acquired', scoring_file_tier=sc.get('tier') if sc else None, scoring_file_sha256=sc.get('sha256') if sc else None,
                           scoring_file_variants=integer(sc.get('data_lines')) if sc else None, gates=gates_for(r, acquired, sc)))
    perf = []
    for r in read_csv(os.path.join(meta_dir, 'pgs_all_metadata_performance_metrics.csv')):
        base = dict(ppm_id=r['PGS Performance Metric (PPM) ID'], pgs_id=r['Evaluated Score'], pss_id=r['PGS Sample Set (PSS)'], pgp_id=r['PGS Publication (PGP) ID'], reported_trait=r['Reported Trait'],
                    covariates=r['Covariates Included in the Model'] or None, other_information=r['PGS Performance: Other Relevant Information'] or None, pmid=r['Publication (PMID)'])
        for col, fam in METRICS:
            if r[col].strip():
                m = metric(r[col])
                perf.append(dict(base, metric_family=fam, metric_name=fam, value=m[1] if m else None, ci_low=m[2] if m else None, ci_high=m[3] if m else None, text=r[col]))
        for part in [p for p in re.split(r';\s*', r['Other Metric(s)']) if p.strip()]:
            m = metric(part)
            perf.append(dict(base, metric_family='other', metric_name=(m[0] if m else None) or part.split(':')[0][:60], value=m[1] if m else None, ci_low=m[2] if m else None,
                             ci_high=m[3] if m else None, text=part))
    samples = []
    for r in read_csv(os.path.join(meta_dir, 'pgs_all_metadata_evaluation_sample_sets.csv')):
        samples.append(dict(stage='evaluation', sample_set_id=r['PGS Sample Set (PSS)'], pgs_id=r['Polygenic Score (PGS) ID'], n_individuals=integer(r['Number of Individuals']),
                            n_cases=integer(r['Number of Cases']), n_controls=integer(r['Number of Controls']), pct_male=num(r['Percent of Participants Who are Male']), age=r['Sample Age'] or None,
                            ancestry_broad=r['Broad Ancestry Category'] or None, ancestry_detail=r['Ancestry (e.g. French, Chinese)'] or None, country=r['Country of Recruitment'] or None,
                            additional_ancestry=r['Additional Ancestry Description'] or None, phenotype_definition=r['Phenotype Definitions and Methods'] or None, followup=r['Followup Time'] or None,
                            gwas_study=r['GWAS Catalog Study ID (GCST...)'] or None, pmid=r['Source PubMed ID (PMID)'] or None, cohorts=r['Cohort(s)'] or None))
    for r in read_csv(os.path.join(meta_dir, 'pgs_all_metadata_score_development_samples.csv')):
        samples.append(dict(stage=r['Stage of PGS Development'] or 'development', sample_set_id=None, pgs_id=r['Polygenic Score (PGS) ID'], n_individuals=integer(r['Number of Individuals']),
                            n_cases=integer(r['Number of Cases']), n_controls=integer(r['Number of Controls']), pct_male=num(r['Percent of Participants Who are Male']), age=r['Sample Age'] or None,
                            ancestry_broad=r['Broad Ancestry Category'] or None, ancestry_detail=r['Ancestry (e.g. French, Chinese)'] or None, country=r['Country of Recruitment'] or None,
                            additional_ancestry=r['Additional Ancestry Description'] or None, phenotype_definition=r['Phenotype Definitions and Methods'] or None, followup=r['Followup Time'] or None,
                            gwas_study=r['GWAS Catalog Study ID (GCST...)'] or None, pmid=r['Source PubMed ID (PMID)'] or None, cohorts=r['Cohort(s)'] or None))
    rows = dict(pgs_model=models, pgs_performance=perf, pgs_sample=samples)
    digests = {}
    for name, rs in rows.items():
        path = os.path.join(out_dir, f'{name}.parquet')
        pq.write_table(pa.Table.from_pylist([{n: r.get(n) for n in SCHEMAS[name].names} for r in rs], schema=SCHEMAS[name]), path, compression='zstd')
        with open(path, 'rb') as f:
            digests[name] = hashlib.sha256(f.read()).hexdigest()
    gate_counts = {}
    for m in models:
        for g in m['gates']:
            gate_counts[g] = gate_counts.get(g, 0) + 1
    manifest = dict(models=len(models), scoring_files_acquired=sum(1 for m in models if m['scoring_file'] == 'acquired'),
                    acquired_not_in_catalogue_metadata=sorted(set(scoring) - {m['pgs_id'] for m in models}), performance_rows=len(perf), sample_rows=len(samples), gates=dict(sorted(gate_counts.items())),
                    ancestry_evaluation_has_non_european=sum(1 for m in models if any(k not in ('European', 'Not Reported') and v > 0 for k, v in json.loads(m['ancestry_evaluation_json']).items())),
                    table_sha256=digests)
    with open(os.path.join(out_dir, 'registry_manifest.json'), 'w') as f:
        json.dump(manifest, f, indent=1, sort_keys=True)
    return manifest


def models_for_trait(models, trait):
    """Ids of the acquired scores mapped to an ontology term. EFO_0004339 and EFO:0004339 are one term: the PGS Catalog writes the first, the index and its views the second."""
    key = (trait or '').strip().replace(':', '_', 1)
    return sorted(i for i, m in models.items() if m.get('scoring_file') == 'acquired' and any((x or '').strip().replace(':', '_', 1) == key for x in m.get('efo_ids') or ()))


def ancestry_meta(model_row):
    """The structure vindex.pgs.score expects for its ancestry gate: development and evaluation distributions as 'group:pct|group:pct' text."""
    return dict(development=model_row.get('ancestry_development'), evaluation=model_row.get('ancestry_evaluation'))
