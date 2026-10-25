"""Fetch per-protein resources for the UniProt proteins that the MaveDB human score sets are about.

Entity set: the UniProt accessions named by our score sets (the mapped-metadata UniProt ID plus UniProt external
identifiers), the same set compat_tests.py joins on. For each accession:

  alphafold_db_mavedb_proteins   GET https://alphafold.ebi.ac.uk/api/prediction/<ACC>      -> data/external/alphafold_db/prediction/<ACC>.json
                                 the canonical entry's cifUrl and plddtDocUrl            -> .../model/<ACC>.cif, .../confidence/<ACC>.json
  interpro_mavedb_proteins       GET https://www.ebi.ac.uk/interpro/api/entry/all/protein/uniprot/<ACC>/ (all pages)
                                                                                          -> data/external/interpro/<ACC>.json

A 404 (AlphaFold) or 204 (InterPro) means the resource has nothing for that protein and is recorded as such, not as a failure.
The canonical AlphaFold model is AF-<ACC>-F1; accessions that only have isoform models are recorded as isoform_only.
Every stored file has its byte count and SHA-256 listed in acquisition/manifest_files/<target>.tsv; the target's record in
manifest.jsonl is ok only when no accession ended in an error. Re-running skips files already stored. Requests are polite:
identifying User-Agent, a few connections, a pause between calls.

  python3 fetch_per_protein.py [--limit N]
"""
import argparse
import concurrent.futures
import datetime
import hashlib
import json
import os
import sys
import time
import urllib.error
import urllib.request

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, '..'))
sys.path.insert(0, HERE)
import acquire  # noqa: E402  (shares the manifest lock and update function)

DATA = os.path.join(ROOT, 'data', 'external')
UA = acquire.UA


def accessions():
    out = set()
    for fn in ('scrape/mavedb-human/score_set_details.jsonl', 'scrape/mavedb-human/score_set_details_accession_targets.jsonl'):
        for line in open(os.path.join(ROOT, fn)):
            for t in json.loads(line)['targetGenes']:
                if t.get('uniprotIdFromMappedMetadata'):
                    out.add(t['uniprotIdFromMappedMetadata'])
                for e in t.get('externalIdentifiers', []):
                    if e['identifier']['dbName'] == 'UniProt':
                        out.add(e['identifier']['identifier'])
    return sorted(out)


def get(url, tries=4):
    """Returns (status, body, headers). 404 and 204 are results, anything else after retries raises."""
    err = None
    for attempt in range(tries):
        if attempt:
            time.sleep(4 * 2 ** attempt)
        try:
            req = urllib.request.Request(url, headers={'User-Agent': UA, 'Accept': 'application/json, */*'})
            with urllib.request.urlopen(req, timeout=90) as r:
                return r.status, r.read(), {k.lower(): v for k, v in r.headers.items()}
        except urllib.error.HTTPError as e:
            if e.code in (404, 204):
                return e.code, b'', {}
            err = e
        except Exception as e:
            err = e
    raise RuntimeError(f'{url}: {err!r}')


def store(path, body):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path + '.part', 'wb') as f:
        f.write(body)
    os.replace(path + '.part', path)
    return len(body), hashlib.sha256(body).hexdigest()


def now():
    return datetime.datetime.now(datetime.timezone.utc).strftime('%Y-%m-%dT%H:%M:%SZ')


def mtime(path):
    return datetime.datetime.fromtimestamp(os.path.getmtime(path), datetime.timezone.utc).strftime('%Y-%m-%dT%H:%M:%SZ')


def row(acc, part, path, code=200):
    b = open(path, 'rb').read()
    return (acc, part, 'stored', len(b), hashlib.sha256(b).hexdigest(), code, mtime(path))


def one_alphafold(acc):
    """Prediction metadata for any accession; the model and per-residue confidence only for the canonical entry AF-<acc>-F1.
    Very long proteins (about 2,700 residues and more) have no canonical model, only models of isoforms; those are recorded
    as isoform_only and their isoform models are not stored, so that no file passes for the canonical protein by mistake."""
    pred = os.path.join(DATA, 'alphafold_db', 'prediction', acc + '.json')
    if not os.path.exists(pred):
        code, body, _ = get(f'https://alphafold.ebi.ac.uk/api/prediction/{acc}')
        time.sleep(0.2)
        if code == 404:
            return [(acc, 'prediction', 'none', 0, '', 404, now())]
        store(pred, body)
    entries = json.load(open(pred))
    rows = [row(acc, 'prediction', pred)]
    canon = next((e for e in entries if e.get('entryId') == f'AF-{acc}-F1'), None)
    for kind, key, ext in (('model', 'cifUrl', 'cif'), ('confidence', 'plddtDocUrl', 'json')):
        dest = os.path.join(DATA, 'alphafold_db', kind, f'{acc}.{ext}')
        if canon is None:
            if os.path.exists(dest):
                os.remove(dest)      # left by an earlier version of this script that took an isoform model for the canonical one
            rows.append((acc, kind, 'isoform_only', 0, '', 200, now()))
            continue
        if not os.path.exists(dest):
            code, b, _ = get(canon[key])
            time.sleep(0.2)
            store(dest, b)
        rows.append(row(acc, kind, dest))
    return rows


def one_interpro(acc):
    dest = os.path.join(DATA, 'interpro', acc + '.json')
    if os.path.exists(dest):
        return [row(acc, 'matches', dest)]
    url = f'https://www.ebi.ac.uk/interpro/api/entry/all/protein/uniprot/{acc}/?page_size=200'
    pages, version = [], None
    while url:
        code, body, hdr = get(url)
        time.sleep(0.4)
        if code == 204:
            return [(acc, 'matches', 'none', 0, '', 204, now())]
        d = json.loads(body)
        version = version or hdr.get('interpro-version')
        pages.append(d)
        url = d.get('next')
    merged = dict(interpro_version=version, accession=acc, count=pages[0]['count'], results=[r for p in pages for r in p['results']])
    store(dest, json.dumps(merged, sort_keys=True).encode())
    return [row(acc, 'matches', dest)]


def run(target_id, source, fn, accs, workers, subdir):
    rows, errors = [], []
    with concurrent.futures.ThreadPoolExecutor(max_workers=workers) as ex:
        futs = {ex.submit(fn, a): a for a in accs}
        for fu in concurrent.futures.as_completed(futs):
            try:
                rows += fu.result()
            except Exception as e:
                errors.append((futs[fu], repr(e)[:200]))
    rows.sort()
    table = os.path.join(HERE, 'manifest_files', target_id + '.tsv')
    os.makedirs(os.path.dirname(table), exist_ok=True)
    with open(table, 'w') as f:
        f.write('accession\tpart\tstatus\tbytes\tsha256\thttp_status\tfetched_utc\n')
        for r in rows:
            f.write('\t'.join(str(x) for x in r) + '\n')
    stored = [r for r in rows if r[2] == 'stored']
    none = sorted({r[0] for r in rows if r[2] == 'none'})
    rec = dict(id=target_id, source=source, url=None, dest=subdir, tier='A', kind='per_entity', ok=not errors and bool(stored),
               fetched_utc=now(), n_accessions=len(accs), n_files=len(stored), bytes=sum(r[3] for r in stored), accessions_without_resource=len(none), isoform_only=sorted({r[0] for r in rows if r[2] == 'isoform_only'}), errors=errors[:10], n_errors=len(errors),
               files_table=os.path.relpath(table, ROOT), sha256=hashlib.sha256('\n'.join('\t'.join(str(x) for x in r[:5]) for r in stored).encode()).hexdigest(), md5=None, source_md5=None,
               md5_verified_against_source=False, size_verified_against_source=False,
               note='per-protein API responses; the services publish no checksums, sha256 values are ours; 404/204 = nothing for that protein')
    acquire.update_manifest(rec)
    return rec


if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('--limit', type=int)
    a = ap.parse_args()
    accs = accessions()[: a.limit] if a.limit else accessions()
    print(len(accs), 'UniProt accessions', flush=True)
    for tid, source, fn, w, sub in (('alphafold_db_mavedb_proteins', 'alphafold_db', one_alphafold, 4, 'alphafold_db/'), ('interpro_mavedb_proteins', 'interpro', one_interpro, 2, 'interpro/')):
        rec = run(tid, source, fn, accs, w, sub)
        print(tid, 'ok' if rec['ok'] else 'ERRORS', {k: rec[k] for k in ('n_accessions', 'n_files', 'bytes', 'accessions_without_resource', 'n_errors')}, 'isoform_only:', len(rec['isoform_only']), flush=True)
        if rec['errors']:
            print('  first errors:', rec['errors'][:3], flush=True)
