"""Collect provenance evidence for catalogue entries.

For each probe in a probes JSON file, fetch the URL (whole body up to --max-bytes, or a byte range),
record status, selected headers, size, sha256 of the fetched bytes and a trimmed text excerpt, and
append a record to evidence/evidence_index.jsonl. Excerpts are stored under evidence/<source_id>/.

A probe is: {"source": id, "name": slug, "url": str, "range": [start, end] (optional),
             "ua": str (optional; override the User-Agent, only where a public API rejects identifying agents; say so in "note"),
             "render": true (optional; render a JavaScript-only page with Playwright and keep its visible text),
             "method": "GET"|"POST" (default GET), "json": body (POST), "gunzip": bool (optional),
             "note": str (optional)}

Usage: python3 probe.py probes.json [--only SOURCE_ID ...] [--probe SOURCE_ID/NAME ...] [--max-bytes N]
Re-running refreshes evidence; the index keeps every fetch (append-only) so changes between runs are auditable.
"""
import argparse
import datetime
import gzip
import hashlib
import html
import json
import os
import re
import sys
import urllib.error
import urllib.request

HERE = os.path.dirname(os.path.abspath(__file__))
EVID = os.path.join(HERE, 'evidence')
INDEX = os.path.join(EVID, 'evidence_index.jsonl')
UA = 'Mozilla/5.0 (dc-catalogue-probe)'
KEEP_HEADERS = ('x-total-results', 'content-type', 'content-length', 'last-modified', 'etag', 'content-range', 'x-amz-meta-md5', 'date')
EXCERPT_CAP = 60_000


def to_text(raw, ctype):
    if raw[:2] == b'\x1f\x8b':
        try:
            raw = gzip.decompress(raw)
        except Exception:
            try:
                import zlib
                raw = zlib.decompressobj(31).decompress(raw)
            except Exception:
                return '[gzip data, not decodable as text in this excerpt]'
    txt = raw.decode('utf-8', 'replace')
    if 'html' in (ctype or '') or txt.lstrip().lower().startswith('<!doctype html') or '<html' in txt[:500].lower():
        txt = re.sub(r'(?is)<(script|style)[^>]*>.*?</\1>', ' ', txt)
        txt = re.sub(r'(?s)<[^>]+>', ' ', txt)
        txt = html.unescape(txt)
        txt = re.sub(r'[ \t\r\f\v]+', ' ', txt)
        txt = re.sub(r'\n\s*\n+', '\n', txt)
    return txt


def render(p, max_bytes):
    """Render a JavaScript-only page with the pre-installed Playwright Chromium and return its visible text."""
    import subprocess
    try:
        r = subprocess.run(['node', os.path.join(HERE, 'render.js'), p['url'], str(max_bytes)],
                           capture_output=True, text=True, timeout=120)
        out = r.stdout.encode()
        err = None if out.strip() and not out.startswith(b'ERR') else (r.stderr.strip()[:200] or out.decode()[:200] or 'empty render')
        return (200 if err is None else 0), {'content-type': 'text/plain; rendered-by=playwright'}, out, err
    except Exception as e:
        return 0, {}, b'', repr(e)


def fetch(p, max_bytes):
    if p.get('render'):
        return render(p, max_bytes)
    headers = {'User-Agent': p.get('ua', UA), 'Accept': '*/*'}
    data = None
    if p.get('range'):
        headers['Range'] = 'bytes=%d-%d' % tuple(p['range'])
    if p.get('method') == 'POST':
        data = json.dumps(p.get('json')).encode()
        headers['Content-Type'] = 'application/json'
        headers['Accept'] = 'application/json'
    req = urllib.request.Request(p['url'], data=data, headers=headers, method=p.get('method', 'GET'))
    try:
        with urllib.request.urlopen(req, timeout=60) as r:
            body = r.read(max_bytes)
            return r.status, dict(r.headers), body, None
    except urllib.error.HTTPError as e:
        try:
            body = e.read(2000)
        except Exception:
            body = b''
        return e.code, dict(e.headers or {}), body, f'HTTPError {e.code}'
    except Exception as e:
        return 0, {}, b'', repr(e)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('probes')
    ap.add_argument('--only', nargs='*')
    ap.add_argument('--probe', nargs='*', help='run only these source/name probes')
    ap.add_argument('--max-bytes', type=int, default=400_000)
    a = ap.parse_args()
    probes = json.load(open(a.probes))
    os.makedirs(EVID, exist_ok=True)
    now = datetime.datetime.now(datetime.timezone.utc).strftime('%Y-%m-%dT%H:%M:%SZ')
    with open(INDEX, 'a') as idx:
        for p in probes:
            if a.only and p['source'] not in a.only:
                continue
            if a.probe and f"{p['source']}/{p['name']}" not in a.probe:
                continue
            status, hdrs, body, err = fetch(p, a.max_bytes)
            low = {k.lower(): v for k, v in hdrs.items()}
            text = to_text(body, low.get('content-type', '')) if body else ''
            d = os.path.join(EVID, p['source'])
            os.makedirs(d, exist_ok=True)
            fn = f"{p['name']}.txt"
            with open(os.path.join(d, fn), 'w', encoding='utf-8') as f:
                f.write(text[:EXCERPT_CAP])
            rec = dict(source=p['source'], name=p['name'], url=p['url'], fetched_utc=now, http_status=status,
                       error=err, bytes_fetched=len(body), sha256=hashlib.sha256(body).hexdigest(),
                       headers={k: low[k] for k in KEEP_HEADERS if k in low}, excerpt=f"evidence/{p['source']}/{fn}",
                       excerpt_truncated=len(text) > EXCERPT_CAP, range=p.get('range'), note=p.get('note'))
            idx.write(json.dumps(rec) + '\n')
            print(f"{status:>4} {p['source']}/{p['name']} ({len(body)} B){' ERR ' + err if err else ''}", file=sys.stderr)


if __name__ == '__main__':
    main()
