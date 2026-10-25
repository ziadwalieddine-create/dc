"""Store the text of a PDF as evidence (the probe collector decodes bodies as text and cannot read PDFs).

  python3 catalogue/derive/pdf_text.py SOURCE NAME URL

Downloads URL, extracts the text with pdftotext -layout, and writes catalogue/evidence/SOURCE/NAME.json with the URL, fetch time, HTTP headers of interest,
the SHA-256 of the PDF and the extracted text. A source record then cites it as "SOURCE/NAME" in its evidence list.
"""
import datetime
import hashlib
import json
import os
import subprocess
import sys
import tempfile
import urllib.request

HERE = os.path.dirname(os.path.abspath(__file__))


def main():
    src, name, url = sys.argv[1:4]
    req = urllib.request.Request(url, headers={'User-Agent': 'Mozilla/5.0 (dc-catalogue-probe)'})
    with urllib.request.urlopen(req, timeout=120) as r:
        raw = r.read()
        hdr = {k.lower(): v for k, v in r.headers.items() if k.lower() in ('content-type', 'content-length', 'last-modified', 'etag', 'x-goog-generation')}
        status = r.status
    with tempfile.TemporaryDirectory() as d:
        pdf = os.path.join(d, 'x.pdf')
        open(pdf, 'wb').write(raw)
        txt = subprocess.run(['pdftotext', '-layout', pdf, '-'], capture_output=True, text=True, check=True).stdout
    out = dict(source=src, name=name, url=url, http_status=status, headers=hdr, fetched_utc=datetime.datetime.now(datetime.timezone.utc).strftime('%Y-%m-%dT%H:%M:%SZ'),
               pdf_bytes=len(raw), pdf_sha256=hashlib.sha256(raw).hexdigest(), text=txt)
    os.makedirs(os.path.join(HERE, '..', 'evidence', src), exist_ok=True)
    json.dump(out, open(os.path.join(HERE, '..', 'evidence', src, name + '.json'), 'w'), indent=1)
    print(f'{src}/{name}: {len(raw)} bytes, {len(txt)} characters of text')


if __name__ == '__main__':
    main()
