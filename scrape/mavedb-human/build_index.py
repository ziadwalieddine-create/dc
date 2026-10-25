"""Build score_sets_index.csv from both detail files (organism-filter and accession-target selections)."""
import csv
import json
import os

HERE = os.path.dirname(os.path.abspath(__file__))
SOURCES = [('score_set_details.jsonl', 'organism_filter'),
           ('score_set_details_accession_targets.jsonl', 'accession_target')]


def ext(t, db):
    return '|'.join(sorted({e['identifier']['identifier'] for e in t.get('externalIdentifiers', []) if e['identifier']['dbName'] == db}))


rows = []
for fn, selection in SOURCES:
    for line in open(os.path.join(HERE, fn)):
        r = json.loads(line)
        r['_selection'] = selection
        rows.append(r)
rows.sort(key=lambda r: r['urn'])
assert len({r['urn'] for r in rows}) == len(rows), 'duplicate URN across selections'

with open(os.path.join(HERE, 'score_sets_index.csv'), 'w', newline='') as f:
    w = csv.writer(f)
    w.writerow(['urn', 'selection', 'title', 'published_date', 'num_variants', 'target_genes', 'target_categories',
                'target_accessions', 'assemblies', 'ensembl_ids', 'uniprot_ids', 'refseq_ids', 'experiment_urn',
                'publications', 'license', 'processing_state', 'mapping_state', 'score_columns', 'count_columns'])
    for r in rows:
        T = r['targetGenes']
        pubs = '|'.join(f"{p['dbName']}:{p['identifier']}" for p in (r.get('primaryPublicationIdentifiers') or []) + (r.get('secondaryPublicationIdentifiers') or []))
        dc = r.get('datasetColumns') or {}
        acc = [t['targetAccession'] for t in T if t.get('targetAccession')]
        w.writerow([r['urn'], r['_selection'], r['title'], r['publishedDate'], r['numVariants'],
                    '; '.join(t['name'] for t in T), '; '.join(sorted({t['category'] for t in T if t.get('category')})),
                    '|'.join(sorted({a['accession'] for a in acc})), '|'.join(sorted({a['assembly'] for a in acc if a.get('assembly')})),
                    ';'.join(filter(None, (ext(t, 'Ensembl') for t in T))), ';'.join(filter(None, (ext(t, 'UniProt') for t in T))),
                    ';'.join(filter(None, (ext(t, 'RefSeq') for t in T))), r['experiment']['urn'], pubs,
                    (r.get('license') or {}).get('shortName'), r.get('processingState'), r.get('mappingState'),
                    '|'.join(dc.get('scoreColumns') or []), '|'.join(dc.get('countColumns') or [])])
print(len(rows), 'score sets indexed')
