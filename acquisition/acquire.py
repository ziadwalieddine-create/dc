"""Download external resources listed in targets.json, verify them and record a manifest.

For each target: stream to data/external/<dest> (resumable via HTTP Range), compute md5 and sha256, compare with
the size/checksum the source publishes (expected_size, expected_md5, md5_url, expected_md5_base64, expected_bsd_sum), and append a
record to acquisition/manifest.jsonl. Files already present and recorded as verified are skipped.

  python3 acquire.py --list                       # show targets and status
  python3 acquire.py --tier A                     # download all tier-A targets
  python3 acquire.py --ids clinvar_vcf_grch38     # specific targets
  python3 acquire.py --tier A --max-mb 700        # skip targets whose expected size exceeds the cap

Data lives under data/external/ (git-ignored). Large files must not be committed; licences differ per source
(see catalogue/sources/*.json) and a private repository does not override redistribution terms.
"""
import argparse
import base64
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
DATA = os.path.join(ROOT, 'data', 'external')
MANIFEST = os.path.join(HERE, 'manifest.jsonl')
TARGETS = os.path.join(HERE, 'targets.json')
UA = 'dc-catalogue-acquire/1.0 (research; github ziadwalieddine-create/dc)'


def load_manifest():
    m = {}
    if os.path.exists(MANIFEST):
        for line in open(MANIFEST):
            if line.strip():
                r = json.loads(line)
                m[r['id']] = r
    return m


def save_manifest(m):
    with open(MANIFEST, 'w') as f:
        for k in sorted(m):
            f.write(json.dumps(m[k], sort_keys=True) + '\n')


def http_get_text(url):
    req = urllib.request.Request(url, headers={'User-Agent': UA})
    return urllib.request.urlopen(req, timeout=60).read().decode('utf-8', 'replace')


def expected_md5(t):
    if t.get('expected_md5'):
        return t['expected_md5'].lower()
    if t.get('expected_md5_base64'):
        return base64.b64decode(t['expected_md5_base64']).hex()
    if t.get('md5_url'):
        try:
            return http_get_text(t['md5_url']).split()[0].lower()
        except Exception:
            return None
    return None


def download(t, dest, tries=4):
    part = dest + '.part'
    os.makedirs(os.path.dirname(dest), exist_ok=True)
    err = None
    for attempt in range(tries):
        if attempt:
            time.sleep(5 * 2 ** attempt)
        try:
            start = os.path.getsize(part) if os.path.exists(part) else 0
            headers = {'User-Agent': UA}
            if start:
                headers['Range'] = f'bytes={start}-'
            req = urllib.request.Request(t['url'], headers=headers)
            with urllib.request.urlopen(req, timeout=120) as r:
                hdrs = {k.lower(): v for k, v in r.headers.items()}
                mode = 'ab' if (start and r.status == 206) else 'wb'
                if mode == 'wb':
                    start = 0
                with open(part, mode) as f:
                    while True:
                        chunk = r.read(1 << 20)
                        if not chunk:
                            break
                        f.write(chunk)
            os.replace(part, dest)
            return hdrs, None
        except Exception as e:
            err = repr(e)
    return {}, err


def digest(path):
    md5, sha = hashlib.md5(), hashlib.sha256()
    with open(path, 'rb') as f:
        for chunk in iter(lambda: f.read(1 << 20), b''):
            md5.update(chunk)
            sha.update(chunk)
    return md5.hexdigest(), sha.hexdigest()


def bsd_sum(path):
    """Unix `sum` (BSD algorithm, 1 KiB blocks), as used by Ensembl FTP CHECKSUMS files.

    Uses coreutils `sum` when present (fast); otherwise a pure-Python implementation of the same algorithm."""
    import shutil
    import subprocess
    if shutil.which('sum'):
        out = subprocess.run(['sum', path], capture_output=True, text=True, check=True).stdout.split()
        return int(out[0]), int(out[1])
    checksum, size = 0, 0
    with open(path, 'rb') as f:
        for chunk in iter(lambda: f.read(1 << 20), b''):
            size += len(chunk)
            for b in chunk:
                checksum = ((checksum >> 1) + ((checksum & 1) << 15) + b) & 0xFFFF
    return checksum, (size + 1023) // 1024


def signature_problem(path, kind):
    """Reject HTML landing/error pages saved under a binary or data file name, and wrong container types."""
    with open(path, 'rb') as f:
        head = f.read(512)
    low = head.lstrip().lower()
    is_html = low.startswith(b'<!doctype html') or low.startswith(b'<html')
    if kind != 'html' and is_html:
        return 'content is an HTML page, not ' + kind
    if kind.endswith(('gz', 'bgz')) and head[:2] != b'\x1f\x8b':
        return 'missing gzip signature for ' + kind
    if kind == 'zip' and head[:2] != b'PK':
        return 'missing zip signature'
    return None


def process(t, manifest, force=False):
    dest = os.path.join(DATA, t['dest'])
    prev = manifest.get(t['id'])
    if not force and prev and prev.get('ok') and os.path.exists(dest) and os.path.getsize(dest) == prev['bytes']:
        return prev, 'skipped (already verified)'
    hdrs, err = download(t, dest)
    rec = dict(id=t['id'], source=t['source'], url=t['url'], dest=t['dest'], tier=t.get('tier'), ok=False,
               fetched_utc=datetime.datetime.now(datetime.timezone.utc).strftime('%Y-%m-%dT%H:%M:%SZ'))
    if err:
        rec['error'] = err
        return rec, 'FAILED ' + err[:80]
    md5, sha = digest(dest)
    size = os.path.getsize(dest)
    exp_md5 = expected_md5(t)
    problems = []
    bad_sig = signature_problem(dest, t.get('kind', ''))
    if bad_sig:
        problems.append(bad_sig)
    if t.get('expected_size') and size != t['expected_size']:
        problems.append(f"size {size} != expected {t['expected_size']}")
    if exp_md5 and md5 != exp_md5:
        problems.append(f'md5 {md5} != source {exp_md5}')
    rec.update(bytes=size, md5=md5, sha256=sha, source_md5=exp_md5, md5_verified_against_source=bool(exp_md5) and md5 == exp_md5,
               size_verified_against_source=bool(t.get('expected_size')) and size == t['expected_size'],
               last_modified=hdrs.get('last-modified'), etag=hdrs.get('etag'), content_type=hdrs.get('content-type'))
    if t.get('expected_bsd_sum'):
        got = bsd_sum(dest)
        rec['bsd_sum'] = list(got)
        rec['bsd_sum_verified_against_source'] = list(got) == list(t['expected_bsd_sum'])
        if not rec['bsd_sum_verified_against_source']:
            problems.append(f"bsd sum {got} != source {t['expected_bsd_sum']}")
    if problems:
        rec['error'] = '; '.join(problems)
        os.replace(dest, dest + '.rejected')
        return rec, 'REJECTED ' + rec['error']
    rec['ok'] = True
    return rec, f'ok {size/1e6:.1f} MB' + (' md5 verified' if rec['md5_verified_against_source'] else (' bsd-sum verified' if rec.get('bsd_sum_verified_against_source') else ' (no source checksum)'))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--tier')
    ap.add_argument('--ids', nargs='*')
    ap.add_argument('--max-mb', type=float)
    ap.add_argument('--list', action='store_true')
    ap.add_argument('--force', action='store_true')
    a = ap.parse_args()
    targets = json.load(open(TARGETS))
    manifest = load_manifest()
    sel = [t for t in targets if (not a.tier or t.get('tier') == a.tier) and (not a.ids or t['id'] in a.ids)]
    if a.list:
        for t in sel:
            m = manifest.get(t['id'])
            st = 'verified' if m and m.get('ok') else ('failed' if m else 'not downloaded')
            print(f"{t['id']:34s} tier {t.get('tier')}  {st:14s} {t.get('expected_size', 0)/1e6:9.1f} MB  {t['url'][:70]}")
        return
    for t in sel:
        if a.max_mb and t.get('expected_size', 0) > a.max_mb * 1e6:
            print(f"{t['id']}: skipped (expected {t['expected_size']/1e6:.0f} MB > cap)")
            continue
        rec, msg = process(t, manifest, a.force)
        manifest[t['id']] = rec
        save_manifest(manifest)
        print(f"{t['id']}: {msg}", flush=True)


if __name__ == '__main__':
    main()
