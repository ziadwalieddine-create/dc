"""List the Ensembl VEP plugins from the documentation page and keep what the page says about their data and licences.

  python3 catalogue/derive/vep_plugins.py

Fetches https://www.ensembl.org/info/docs/tools/vep/script/vep_plugins.html and writes catalogue/evidence/ensembl_vep/plugin_list.json: one entry per plugin
(name, category, external libraries, developer, the plugin's description, and the sentences that mention a licence, academic or commercial use, or the
data files the plugin needs). The page is the only source; nothing is added from memory. The complete page is 270 KB and exceeds the 60 KB excerpt
limit of probe.py, which is why this derivation exists.
"""
import datetime
import hashlib
import html
import json
import os
import re
import urllib.request
from html.parser import HTMLParser

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, '..', 'evidence', 'ensembl_vep', 'plugin_list.json')
URL = 'https://www.ensembl.org/info/docs/tools/vep/script/vep_plugins.html'


def clean(s):
    s = re.sub(r'(?s)<(script|style)[^>]*>.*?</\1>', ' ', s)
    s = re.sub(r'(?s)<[^>]+>', ' ', s)
    return re.sub(r'\s+', ' ', html.unescape(s)).strip()


class Cells(HTMLParser):
    """Collect the text of the top-level cells of one plugin row, ignoring cells of tables nested inside a description."""

    def __init__(self):
        super().__init__(convert_charrefs=True)
        self.cells, self.depth, self.cur = [], 0, None
        self.skip = 0

    def handle_starttag(self, tag, attrs):
        if tag == 'table':
            self.depth += 1
        elif tag in ('script', 'style'):
            self.skip += 1
        elif tag == 'td' and self.depth == 0:
            self.cur = []
        elif tag in ('p', 'br', 'li') and self.cur is not None:
            self.cur.append(' ')

    def handle_endtag(self, tag):
        if tag == 'table':
            self.depth = max(0, self.depth - 1)
        elif tag in ('script', 'style'):
            self.skip = max(0, self.skip - 1)
        elif tag == 'td' and self.depth == 0 and self.cur is not None:
            self.cells.append(re.sub(r'\s+', ' ', ''.join(self.cur)).strip())
            self.cur = None

    def handle_data(self, data):
        if self.cur is not None and not self.skip:
            self.cur.append(data)


def main():
    req = urllib.request.Request(URL, headers={'User-Agent': 'Mozilla/5.0 (dc-catalogue-probe)'})
    raw = urllib.request.urlopen(req, timeout=120).read()
    page = raw.decode('utf-8', 'replace')
    starts = list(re.finditer(r'<tr id="([^"]+)" class="[^"]*plugin_row" data-category="([^"]*)">', page))
    plugins = []
    for i, m in enumerate(starts):
        end = starts[i + 1].start() if i + 1 < len(starts) else page.find('</tbody>', m.end())
        c = Cells()
        c.feed(page[m.end():end])
        cells = c.cells
        desc = re.sub(r'\s*\bmore\b\s*', ' ', cells[1], count=1).strip() if len(cells) > 1 else ''
        sents = re.split(r'(?<=[.!?])\s+', desc)
        key = [x for x in sents if re.search(r'licen[cs]e|academic|commercial|non-commercial|terms of use|freely|open access|CC BY|download|ftp\.|https?://', x, re.I)]
        plugins.append(dict(plugin=m.group(1), category=m.group(2), external_libraries=cells[3] if len(cells) > 3 else '', developer=cells[4] if len(cells) > 4 else '',
                            description=desc[:1800], licence_or_data_sentences=[x[:400] for x in key][:8]))
    out = dict(source_url=URL, fetched_utc=datetime.datetime.now(datetime.timezone.utc).strftime('%Y-%m-%dT%H:%M:%SZ'), page_bytes=len(raw), page_sha256=hashlib.sha256(raw).hexdigest(),
               n_plugins=len(plugins), categories=sorted({p['category'] for p in plugins}), plugins=plugins)
    json.dump(out, open(OUT, 'w'), indent=1)
    print(f"{len(plugins)} plugins in {len(out['categories'])} categories; {sum(1 for p in plugins if p['description'])} with a description; {sum(1 for p in plugins if p['licence_or_data_sentences'])} with licence or data sentences")


if __name__ == '__main__':
    main()
